
abstract type AbstractBranch end

mutable struct AbstractNode
    id::Int
    parent::Union{Nothing,AbstractNode}
    depth::Int
    branch::AbstractBranch

    solution_status::MOI.TerminationStatusCode
    bound::Real
    solution::Dict{Int,Real}

    function AbstractNode(
            id = -1, 
            parent = nothing, 
            depth = 0, 
            solution_status = MOI.OPTIMIZE_NOT_CALLED,
            bound = -Inf,
            solution = Dict())
        N = new()
        N.id = id
        N.parent = parent
        N.depth = depth
        N.solution_status = solution_status
        N.bound = bound
        N.solution = solution
        return N
    end
end

AbstractNode(parent::AbstractNode) = AbstractNode(
    -1, 
    parent, 
    parent.depth + 1, 
    MOI.OPTIMIZE_NOT_CALLED,
    parent.bound,
    Dict())

# return node solution
node_solution(node::AbstractNode) = node.solution

# empty branching function
branch(node::AbstractNode)::Vector{AbstractNode} = []

# empty bounding function
function bound!(node::AbstractNode) end

# empty heuristics function
function heuristics!(node::AbstractNode) end

# apply changes (branched information) from ancestor to node
function apply_changes!(node::AbstractNode, ancestor::AbstractNode) end

apply_changes!(node::AbstractNode) = apply_changes!(node, node)
