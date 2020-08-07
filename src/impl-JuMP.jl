using JuMP

mutable struct JuMPNode <: AbstractNode
    id::Int
    parent::Union{Nothing,JuMPNode}
    depth::Int

    model::Union{Nothing,JuMP.Model}
    reference_map::JuMP.ReferenceMap

    variable_lower_bound_changes::Dict{JuMP.VariableRef,Real}
    variable_upper_bound_changes::Dict{JuMP.VariableRef,Real}

    solution_status::MOI.TerminationStatusCode
    bound::Real
    solution::Dict{JuMP.VariableRef,Real}

    function AbstractNode(parent = nothing, depth = 0, dual_bound = -Inf)
        return new(-1, parent, depth, nothing, nothing, Dict(), Dict(), Dict(), MOI.OPTIMIZE_NOT_CALLED, Inf, dual_bound, Dict())
    end
end

# basic bounding function
function bound!(node::JuMPNode)
    optimize!(node.model)
    node.solution_status = JuMP.termination_status(node.model)
    node.bound = JuMP.dual_objective_value(node.model)
    vrefs = JuMP.all_variables(node.model)
    for v in vrefs
        node.solution[v] = JuMP.value(v)
    end
end

# return the next search node
function next_node(tree::AbstractTree)::AbstractNode
    sort!(tree)
    node = Base.pop!(tree.nodes)
    node.model, node.reference_map = JuMP.copy_model(tree.model)
    JuMP.set_optimizer(node.model, tree.optimizer)
    apply_changes!(node)
    return node
end

# apply changes from the current and acestor nodes
function apply_changes!(node::JuMPNode, ancestor::JuMPNode)
    @assert !isnothing(node.model)

    for (v,bd) in ancestor.variable_lower_bound_changes
        JuMP.set_lower_bound(node.reference_map[v], bd)
    end

    for (v,bd) in ancestor.variable_upper_bound_changes
        JuMP.set_upper_bound(node.reference_map[v], bd)
    end

    # TODO: parent.constraint_changes

    if !isnothing(ancestor.parent)
        apply_changes!(node, ancestor.parent)
    end
end
