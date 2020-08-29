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

analyze_tree_bounds(tree)