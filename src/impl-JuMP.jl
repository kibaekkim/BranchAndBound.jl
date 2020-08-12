using JuMP

struct VariableBranch <: AbstractBranch
    lb::Dict{JuMP.VariableRef,Real}
    ub::Dict{JuMP.VariableRef,Real}
end

empty_variable_bound() = Dict{JuMP.VariableRef,Real}()

mutable struct JuMPNode{T<:AbstractBranch} <: AbstractNode
    @abstract_node_fields

    model::Union{Nothing,JuMP.Model}
    auxiliary_data::Dict

    function JuMPNode{T}(
        model = nothing,
        parent = nothing, 
        depth = 0, 
        bound = -Inf,
        id = -1, 
        branch = nothing,
        solution_status = MOI.OPTIMIZE_NOT_CALLED,
        solution = Dict{Any,Real}(),
        auxiliary_data = Dict()) where T<:AbstractBranch
        return new{T}(id, parent, depth, branch, solution_status, bound, solution, model, auxiliary_data)
    end
end

JuMPNode{T}(parent::JuMPNode{T}) where T<:AbstractBranch = JuMPNode{T}(nothing, parent, parent.depth + 1, parent.bound)
JuMPNode(parent::JuMPNode{T}) where T<:AbstractBranch = JuMPNode{T}(parent)

# Create child node with abstract branch object
function create_child_node(current_node::JuMPNode{T}, branch::AbstractBranch) where T<:AbstractBranch
    node = JuMPNode(current_node)
    node.branch = branch
    return node
end

# Create child node with variable bounds
create_child_node(current_node::JuMPNode{T}, variable::JuMP.VariableRef, lb::Real, ub::Real) where T<:AbstractBranch = create_child_node(current_node, VariableBranch(Dict(variable=>lb), Dict(variable=>ub)))
create_child_node_with_lb(current_node::JuMPNode{T}, variable::JuMP.VariableRef, lb::Real) where T<:AbstractBranch = create_child_node(current_node, VariableBranch(Dict(variable=>lb), empty_variable_bound()))
create_child_node_with_ub(current_node::JuMPNode{T}, variable::JuMP.VariableRef, ub::Real) where T<:AbstractBranch = create_child_node(current_node, VariableBranch(empty_variable_bound(), Dict(variable=>ub)))

# basic bounding function
function bound!(node::JuMPNode)
    JuMP.optimize!(node.model)
    node.solution_status = JuMP.termination_status(node.model)

    if node.solution_status in [MOI.OPTIMAL, MOI.LOCALLY_SOLVED]
        node.bound = JuMP.objective_value(node.model)
        vrefs = JuMP.all_variables(node.model)
        for v in vrefs
            node.solution[v] = JuMP.value(v)
        end
    elseif node.solution_status == MOI.INFEASIBLE
        node.bound = -Inf
    elseif node.solution_status == MOI.DUAL_INFEASIBLE
        node.bound = Inf
    else
        @warn "Unexpected node solution status: $(node.solution_status)"
        node.bound = -Inf
    end

    # Rollback the branch objects
    adjust_branch!(node.auxiliary_data["bounds_changed"])
end

function apply_changes!(node::JuMPNode)
    branch_objects = AbstractBranch[]
    bounds_changed = AbstractBranch[]

    # Bracktracking tree to the root and collecting branch objects
    pnode = node
    push!(branch_objects, pnode.branch)
    while !isnothing(pnode.parent)
        pnode = pnode.parent
        push!(branch_objects, pnode.branch)
    end
    @assert isnothing(pnode.parent)
    @assert !isnothing(pnode.model)

    # Shallow copy of the root model
    node.model = pnode.model

    # Mark bounds to be changed
    mark_bound_changes!(node, branch_objects)

    # Apply branch objects
    adjust_branch!(branch_objects)
end

function adjust_branch!(branch_objects::Array{AbstractBranch,1})
    for branch in branch_objects
        for (j,v) in branch.lb
            JuMP.set_lower_bound(j, v)
        end
        for (j,v) in branch.ub
            JuMP.set_upper_bound(j, v)
        end
    end
end

function mark_bound_changes!(node::JuMPNode, branch_objects::Array{AbstractBranch,1})
    node.auxiliary_data["bounds_changed"] = AbstractBranch[]
    for branch in branch_objects
        lb_changed = Dict{JuMP.VariableRef,Real}()
        ub_changed = Dict{JuMP.VariableRef,Real}()
        for (j,v) in branch.lb
            lb_changed[j] = JuMP.has_lower_bound(j) ? JuMP.lower_bound(j) : -Inf
        end
        for (j,v) in branch.ub
            ub_changed[j] = JuMP.has_upper_bound(j) ? JuMP.upper_bound(j) : Inf
        end
        Base.push!(node.auxiliary_data["bounds_changed"], VariableBranch(lb_changed, ub_changed))
    end
end
