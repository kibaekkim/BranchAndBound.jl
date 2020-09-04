# This code is used to print out branching decisions for debugging purposes

# include("./spatial_bac.jl")

function analyze_tree_bounds(tree::BB.AbstractTree)
    all_nodes = union(tree.nodes, tree.processed)
    sort!(all_nodes, by = x->x.depth)
    curr_depth = -1
    for node in all_nodes
        if node.depth != curr_depth
            curr_depth = node.depth
            println("================== Depth $(curr_depth) =====================")
        end
        if node.depth > 0
            println("    Node id: $(node.id), parent id: $(node.parent.id), index: ($(node.branch.i), $(node.branch.j)), bounds: $(round.(node.branch.bounds, digits=3))")
        end
    end
end

function analyze_node_bounds(tree::BB.AbstractTree)
    nodes = tree.processed
    sort!(nodes, by = x->x.depth, rev = true)
    visited = []
    count = 0
    for node in nodes
        pnode = node
        while !isnothing(pnode.parent) && !(pnode.id in visited)
            push!(visited, pnode.id)
            if pnode.bound < pnode.parent.bound
                count = count + 1
                println("Parent node $(pnode.parent.id) has bound $(pnode.parent.bound), child node $(pnode.id) has bound $(pnode.bound)")
            end
            pnode = pnode.parent
        end
    end
    println("Percent of branches that have reverse bound order: $(float(count / (length(nodes) - 1)) * 100)%")
end

# analyze_tree_bounds(tree)
analyze_node_bounds(tree)