using BranchAndBound
using JuMP
using GLPK

const BB = BranchAndBound

# continuous relaxation of MILP
m = Model(GLPK.Optimizer)

# The following will be supported in the later version of JuMP.
# @variable(m, x >= 0, Int)
# @variable(m, y >= 0, Int)
# undo_relax = JuMP.relax_integrality(m)

@variable(m, x >= 0) # assume integrality
@variable(m, y >= 0) # assume integrality

@objective(m, Min, -x - y)
@constraint(m, -2*x + 2*y >= 1)
@constraint(m, -8*x + 10*y <= 13)

# Initialize BNB tree
tree = BB.initialize_tree(m)

# User-defined branch function
function BB.branch!(tree::BB.AbstractTree, node::BB.AbstractNode)
    @info " Node id $(node.id), status $(node.solution_status), bound $(node.bound)"
    if node.bound >= tree.best_incumbent
        @info " Fathomed by bound"
    elseif node.solution_status in [MOI.OPTIMAL, MOI.LOCALLY_SOLVED]
        if !isinteger(node.solution[x])
            @info " Branch at x = $(node.solution[x]), x >= $(ceil(node.solution[x])), x <= $(floor(node.solution[x]))"
            child1 = BB.create_child_node_with_lb(node, x, ceil(node.solution[x]))
            child2 = BB.create_child_node_with_ub(node, x, floor(node.solution[x]))
            BB.push!(tree, child1)
            BB.push!(tree, child2)
        elseif !isinteger(node.solution[y])
            @info " Branch at y = $(node.solution[y]), y >= $(ceil(node.solution[y])), y <= $(floor(node.solution[y]))"
            child1 = BB.create_child_node_with_lb(node, y, ceil(node.solution[y]))
            child2 = BB.create_child_node_with_ub(node, y, floor(node.solution[y]))
            BB.push!(tree, child1)
            BB.push!(tree, child2)
        elseif node.bound < tree.best_incumbent
            @info " New incumbent bound: $(node.bound)"
            tree.best_incumbent = node.bound
        end
    else
        @info " Fathomed by solution status: $(node.solution_status)"
    end
end

function BB.termination(tree::BB.AbstractTree)
    @info "Tree nodes: processed $(length(tree.processed)), left $(length(tree.nodes)), total $(tree.node_counter), best bound $(tree.best_bound), best incumbent $(tree.best_incumbent)"
    if BB.isempty(tree)
        @info "Completed the tree search"
        return true
    end
    if length(tree.processed) >= 20
        @info "Reached node limit"
        return true
    end
    return false
end

BB.run(tree)
