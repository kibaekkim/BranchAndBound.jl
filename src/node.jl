abstract type AbstractNode end
abstract type AbstractBranch end

push!(branches::Array{AbstractBranch,1}, b::AbstractBranch) = Base.pushfirst!(branches, b)
push!(branches::Array{AbstractBranch,1}, b::Nothing) = branches

macro abstract_node_fields() 
    return esc(quote
        id::Int
        parent::Union{Nothing,AbstractNode}
        depth::Int
        branch::Union{Nothing,T}
        solution_status::MOI.TerminationStatusCode
        bound::Real
        solution::Dict{Any,Real}
    end)
end

# empty branching function
branch(node::AbstractNode)::Vector{AbstractNode} = []

# empty bounding function
function bound!(node::AbstractNode) end

# empty heuristics function
function heuristics!(node::AbstractNode) end

# apply changes (branched information) from ancestor to node
function apply_changes!(node::AbstractNode, ancestor::AbstractNode) end

apply_changes!(node::AbstractNode) = apply_changes!(node, node)
