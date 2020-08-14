
mutable struct AbstractTree
    node_counter::Int
    nodes::Vector{AbstractNode}
    processed::Vector{AbstractNode}
    best_bound::Real     # e.g., lower bound in minimizing MILP
    best_incumbent::Real # e.g., upper bound in minimizing MILP

    function AbstractTree(node_counter::Int = 0)
        return new(node_counter, [], [], -Inf, Inf)
    end
end

function run(tree::AbstractTree)
    while !termination(tree)
        node = next_node(tree)
        update_best_bound!(tree)

        bound!(node)
        heuristics!(node)
        processed!(tree, node)
        branch!(tree, node)
    end
end

# basic termination function
termination(tree::AbstractTree) = isempty(tree)
isempty(tree::AbstractTree) = Base.isempty(tree.nodes)

# add node to tree
function push!(tree::AbstractTree, node::AbstractNode)
    tree.node_counter += 1
    node.id = tree.node_counter
    Base.push!(tree.nodes, node)
    sort!(tree)
end

# add a set of nodes to tree
function push!(tree::AbstractTree, nodes::Vector{T}) where T<:AbstractNode
    for node in nodes
        push!(tree, node)
    end
end

# mark node as processed
processed!(tree::AbstractTree, node::AbstractNode) = Base.push!(tree.processed, node)

function update_best_bound!(tree::AbstractTree)
    if !isempty(tree)
        tree.best_bound = Base.minimum([node.bound for node in tree.nodes])
    end
end

function branch!(tree::AbstractTree, node::AbstractNode)
    children = branch(node)
    push!(tree, children)
end

# return the next search node
function next_node(tree::AbstractTree)
    sort!(tree)
    node = Base.pop!(tree.nodes)
    apply_changes!(node)
    return node
end

# best bound
sort!(tree::AbstractTree) = Base.sort!(tree.nodes, by=x->x.bound, rev=true)
