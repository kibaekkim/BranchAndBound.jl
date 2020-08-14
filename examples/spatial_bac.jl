using BranchAndBound, JuMP, PowerModels, DeNet
using LinearAlgebra
const BB = BranchAndBound
const PM = PowerModels

mutable struct NodeWRPowerModel <: AbstractWRModel @pm_fields end

mutable struct SpatialBCBranch <: BB.AbstractBranch
    i::Int
    j::Int
    wii::JuMP.VariableRef
    wjj::JuMP.VariableRef
    wr::JuMP.VariableRef
    wi::JuMP.VariableRef
    mod_branch::Union{ComplexVariableBranch, VariableBranch} # this capture either a changed bound for Wij and Tij, or a changed bound for Wii/Wjj 
    bounds::NTuple{6, <:Real} # Lii, Uii, Ljj, Ujj, Lij, Uij
    valid_ineq_coeffs::NTuple{5, <:Real} # π coefficients
end

mutable struct ComplexVariableBranch <: BB.AbstractBranch
    lb::Dict{Tuple{JuMP.VariableRef, JuMP.VariableRef},Real}
    ub::Dict{Tuple{JuMP.VariableRef, JuMP.VariableRef},Real}
end

# Create child node with variable bounds for complex variables
# function BB.create_child_node(current_node::JuMPNode{ComplexVariableBranch},variable_pair::Tuple{JuMP.VariableRef, JuMP.VariableRef}, lb::Real, ub::Real)
#     create_child_node(current_node, ComplexVariableBranch(Dict(variable_pair=>lb), Dict(variable_pair=>ub)))
# end


# modified from PM.build_opf 
# get rid of voltage angle difference constraints
function PM.build_opf(pm::PM.NodeWRPowerModel)
    variable_bus_voltage(pm)
    variable_gen_power(pm)
    variable_branch_power(pm)
    variable_dcline_power(pm)

    objective_min_fuel_and_flow_cost_mod(pm)

    constraint_model_voltage(pm)

    for i in setdiff(ids(pm, :bus), ids(pm, :cut_bus))
        constraint_power_balance(pm, i)
    end

    for i in ids(pm, :branch)
        constraint_ohms_yt_from(pm, i)
        constraint_ohms_yt_to(pm, i)

        constraint_thermal_limit_from(pm, i)
        constraint_thermal_limit_to(pm, i)
    end

    for i in ids(pm, :dcline)
        constraint_dcline_power_losses(pm, i)
    end

    wr = var(pm, :wr)
    for x in wr
        if !JuMP.has_lower_bound(x) || JuMP.lower_bound(x) < 0
            @constraint(pm.model, x >= 0)
        end
    end
end

function sigmoid(x::Real)::Real
    x == 0 ? res = 0 : res = (sqrt(1 + x^2) - 1) / x
    return res
end

# read data, and formulate root model using PowerModels
file = "../data/case9.m"
data = prase_file(file)
pm = instantiate_model(data, NodeWRPowerModel, build_opf)

# collect data
Lii = Dict(i => bus["vmin"]^2 for (i, bus) in ref(pm, :bus))
Uii = Dict(i => bus["vmax"]^2 for (i, bus) in ref(pm, :bus))
Lij = Dict((i,j) => tan(branch["angmin"]) for ((i,j), branch) in ref(pm, :buspairs))
Uij = Dict((i,j) => tan(branch["angmax"]) for ((i,j), branch) in ref(pm, :buspairs))

# initialize branch-and-cut tree
node = BB.JuMPNode(pm.model)
node.auxiliary_data["PM"] = pm
node.auxiliary_data["Lii"] = Lii
node.auxiliary_data["Uii"] = Uii
node.auxiliary_data["Lij"] = Lij
node.auxiliary_data["Uij"] = Uij
tree = BB.initialize_tree()

function find_root(node::BB.AbstractNode)::BB.AbstractNode
    pnode = node
    while !isnothing(pnode.parent)
        pnode = pnode.parent
    end
    return pnode
end

function find_min_eigen(node::BB.AbstractNode)::Tuple{Int64, Int64}
    pm = find_root(node).auxiliary_data["PM"]
    w = var(pm, :w)
    wr = var(pm, :wr)
    wi = var(pm, :wi)
    min_lambda = Inf
    min_id = ()
    node_solution = node.solution
    for ((i,j),_) in ref(pm, :buspairs)
        lambda = 0.5 * (node_solution[w[i]] - node_solution[w[j]] - norm( [node_solution[w[i]] - node_solution[w[j]], 2 * wr[(i,j)], 2 * wi[(i,j)]] ) )
        if lambda < min_lambda
            min_id = (i,j)
            min_lambda = lambda
        end
    end
    return min_id
end

function get_LU_from_branches(node::BB.AbstractNode)::NTuple{2, Dict}
    Lii = Dict()
    Uii = Dict()
    Lij = Dict()
    Uij = Dict()
    pnode = node
    while !isnothing(pnode.parent)
        sbc_branch = pnode.branch
        mod_branch = sbc_branch.mod_branch
        if mod_branch isa ComplexVariableBranch
            Lij[(sbc_branch.i, sbc_branch.j)] = first(values(mod_branch.lb))
            Uij[(sbc_branch.i, sbc_branch.j)] = first(values(mod_branch.ub))
        elseif mod_branch isa BB.VariableBranch

        else
            error("Invalid branch type $(typeof(mod_branch))")
        end
        pnode = node.parent
    end

end

# User-defined branch function
function BB.branch!(tree::BB.AbstractTree, node::BB.AbstractNode)

    @info " Node id $(node.id), status $(node.solution_status), bound $(node.bound)"
    if node.bound >= tree.best_incumbent
        @info " Fathomed by bound"
    elseif node.solution_status in [MOI.OPTIMAL, MOI.LOCALLY_SOLVED]
        # Get current L, U matrices by back-tracing the branches
        # MVSB / MVWB / RBEB
        # MVSB:
            # Find the (i,j) that attains the greatest minimum eigenvalue of (wr + i wi)
            # pick the complex entry to branch on - need to solve up and down branch problems on three possible entries
            # Create the next node with branch
        
        (i,j) = find_min_eigen(node)
        
        # MVWB:
            # Solve WEV with up and down branch (2 times)
            # 


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

function BB.mark_bound_changes!(node::JuMPNode, branch_objects::Vector{SpatialBCBranch})
    node.auxiliary_data["new_cuts"] = JuMP.ConstraintRef[]
    for branch in branch_objects
        i = branch.i
        j = branch.j
        (Lii, Uii, Ljj, Ujj, Lij, Uij) = branch.bounds
        π = branch.valid_ineq_coeffs
        wii = branch.wii
        wjj = branch.wjj
        wr = branch.wr
        wi = branch.wi
        push!(node.auxiliary_data["new_cuts"], JuMP.@constraint(node.model, Lii <= wii <= Uii))
        push!(node.auxiliary_data["new_cuts"], JuMP.@constraint(node.model, Ljj <= wjj <= Ujj))
        push!(node.auxiliary_data["new_cuts"], JuMP.@constraint(node.model, Lij * wr <= wi))
        push!(node.auxiliary_data["new_cuts"], JuMP.@constraint(node.model, wi <= Uij * wi))
        push!(node.auxiliary_data["new_cuts"], JuMP.@constraint(node.model, π[1] + π[2] * wii + π[3] * wjj + π[4] * wr + π[5] * wi >= Ujj * wii + Uii * wjj - Uii * Ujj))
        push!(node.auxiliary_data["new_cuts"], JuMP.@constraint(node.model, π[1] + π[2] * wii + π[3] * wjj + π[4] * wr + π[5] * wi >= Ljj * wii + Lii * wjj - Lii * Ljj))
    end
end

BB.run(tree)