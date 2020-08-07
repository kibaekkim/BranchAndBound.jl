using JuMP

struct VariableBranch <: AbstractBranch
    lb::Dict{JuMP.VariableRef,Real}
    ub::Dict{JuMP.VariableRef,Real}
end

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
    @show branch_objects

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
            lb_changed[j] = JuMP.lower_bound(j)
        end
        for (j,v) in branch.ub
            lb_changed[j] = JuMP.upper_bound(j)
        end
        Base.push!(node.auxiliary_data["bounds_changed"], VariableBranch(lb_changed, ub_changed))
    end
end
