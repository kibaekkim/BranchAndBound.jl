using BranchAndBound, JuMP, PowerModels, DeNet
using Mosek, MosekTools, Ipopt
using LinearAlgebra
const BB = BranchAndBound
const PM = PowerModels


###########################################################################
#                  Data Structures with Helper Functions                  #
###########################################################################

mutable struct NodeWRMPowerModel <: AbstractWRMModel @pm_fields end

mutable struct ComplexVariableBranch <: BB.AbstractBranch
    lb::Dict{Tuple{JuMP.VariableRef, JuMP.VariableRef},Real} # Lij
    ub::Dict{Tuple{JuMP.VariableRef, JuMP.VariableRef},Real} # Uij
end

function compute_π(bounds::NTuple{6, <:Real})::NTuple{5, <:Real}
    function sigmoid(x::Real)::Real
        x == 0 ? res = 0 : res = (sqrt(1 + x^2) - 1) / x
        return res
    end    
    (Lii, Uii, Ljj, Ujj, Lij, Uij) = bounds
    π0 = -sqrt(Lii * Ljj * Uii * Ujj)
    π1 = -sqrt(Ljj * Ujj)
    π2 = -sqrt(Lii * Uii)
    π3 = (sqrt(Lii) + sqrt(Uii)) * (sqrt(Ljj) + sqrt(Ujj)) * (1 - sigmoid(Lij) * sigmoid(Uij)) / (1 + sigmoid(Lij) * sigmoid(Uij))
    π4 = (sqrt(Lii) + sqrt(Uii)) * (sqrt(Ljj) + sqrt(Ujj)) * (sigmoid(Lij) + sigmoid(Uij)) / (1 + sigmoid(Lij) * sigmoid(Uij))
    return (π0, π1, π2, π3, π4)
end

mutable struct SpatialBCBranch <: BB.AbstractBranch
    i::Int
    j::Int
    wii::JuMP.VariableRef
    wjj::JuMP.VariableRef
    wr::JuMP.VariableRef
    wi::JuMP.VariableRef
    mod_branch::Union{Nothing, ComplexVariableBranch, BB.VariableBranch} # this capture either a changed bound for Wij and Tij, or a changed bound for Wii/Wjj 
    bounds::NTuple{6, <:Real} # Lii, Uii, Ljj, Ujj, Lij, Uij
    valid_ineq_coeffs::NTuple{5, <:Real} # π coefficients # might not be needed
end

# create a SpatialBCBranch, in which mod_branch is Nothing (the exact matrix entry to branch is not decided yet)
# (maybe don't need this... instead decide the branching entry first...)
function create_sbc_branch(i::Int64, j::Int64, prev_node::BB.AbstractNode)
    root = find_root(prev_node)
    pm = root.auxiliary_data["PM"]
    wii = var(pm, :WR)[i,i]
    wjj = var(pm, :WR)[j,j]
    wr = var(pm, :WR)[i,j]
    wi = var(pm, :WI)[i,j]
    bounds = get_LU_from_branches(prev_node, i, j)
    πs = compute_π(bounds)
    return SpatialBCBranch(i, j, wii, wjj, wr, wi, nothing, bounds, πs)
end

# modified from PM.build_opf 
# get rid of voltage angle difference constraints
function PM.build_opf(pm::NodeWRMPowerModel)
    variable_bus_voltage(pm)
    variable_gen_power(pm)
    variable_branch_power(pm)
    variable_dcline_power(pm)

    DeNet.objective_min_fuel_and_flow_cost_mod(pm)

    constraint_model_voltage(pm)

    :cut_bus in keys(ref(pm)) ? cut_bus = ids(pm, :cut_bus) : cut_bus = []
    for i in setdiff(ids(pm, :bus), cut_bus)
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

    # wr = var(pm, :WR)
    # for x in wr
    #     if !JuMP.has_lower_bound(x) || JuMP.lower_bound(x) < 0
    #         @constraint(pm.model, x >= 0)
    #     end
    # end
end

function find_root(node::BB.AbstractNode)::BB.AbstractNode
    pnode = node
    while !isnothing(pnode.parent)
        pnode = pnode.parent
    end
    return pnode
end

function find_min_eigen(node::BB.AbstractNode)::Tuple{Int64, Int64}
    pm = find_root(node).auxiliary_data["PM"]
    wr = var(pm, :WR)
    wi = var(pm, :WI)
    min_lambda = Inf
    min_id = ()
    node_solution = node.solution
    for ((i,j),_) in ref(pm, :buspairs)
        lambda = 0.5 * (node_solution[wr[i,i]] - node_solution[wr[j,j]] - norm( [node_solution[wr[i,i]] - node_solution[wr[j,j]], 2 * node_solution[wr[i,j]], 2 * node_solution[wi[i,j]]] ) )
        if lambda < min_lambda
            min_id = (i,j)
            min_lambda = lambda
        end
    end
    return min_id
end

#=
function get_LU_from_branches(node::BB.AbstractNode)::NTuple{4, Dict}
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
            Lii[sbc_branch.i] = first(values(mod_branch.lb))
            Uii[sbc_branch.i] = first(values(mod_branch.ub))
        else
            error("Invalid branch type $(typeof(mod_branch))")
        end
        pnode = node.parent
    end
    Lii = merge(data.auxiliary_data["Lii"], Lii)
    Uii = merge(data.auxiliary_data["Uii"], Uii)
    Lij = merge(data.auxiliary_data["Lij"], Lij)
    Uij = merge(data.auxiliary_data["Uij"], Uij)
    return (Lii, Uii, Lij, Uij)
end
=#

# Output (Lii, Uii, Ljj, Ujj, Lij, Uij)
function get_LU_from_branches(node::BB.AbstractNode, i::Int64, j::Int64)::NTuple{6, <:Real}
    pnode = node
    (Lii, Uii, Ljj, Ujj, Lij, Uij) = Tuple([NaN for i in 1:6])
    while !isnothing(pnode.parent)
        sbc_branch = pnode.branch
        mod_branch = sbc_branch.mod_branch
        if mod_branch isa ComplexVariableBranch
            if (sbc_branch.i, sbc_branch.j) == (i, j) && isnan(Lij)
                Lij = first(values(mod_branch.lb))
                Uij = first(values(mod_branch.ub))
            end
        elseif mod_branch isa BB.VariableBranch
            if sbc_branch.i == i && isnan(Lii)
                Lii = first(values(mod_branch.lb))
                Uii = first(values(mod_branch.ub))
            elseif sbc_branch.i == j && isnan(Ljj)
                Ljj = first(values(mod_branch.lb))
                Ujj = first(values(mod_branch.ub))
            end
        else
            error("Invalid branch type $(typeof(mod_branch))")
        end
        pnode = node.parent
    end
    isnan(Lii) && (Lii = node.auxiliary_data["Lii"][i])
    isnan(Uii) && (Uii = node.auxiliary_data["Uii"][i])
    isnan(Ljj) && (Ljj = node.auxiliary_data["Lii"][j])
    isnan(Ujj) && (Ujj = node.auxiliary_data["Uii"][j])
    isnan(Lij) && (Lij = node.auxiliary_data["Lij"][(i,j)])
    isnan(Uij) && (Uij = node.auxiliary_data["Uij"][(i,j)])
    return (Lii, Uii, Ljj, Ujj, Lij, Uij)
end

# this function creates the 6 possible node candidates right after the new branch
# needed for MVSB, where we need to solve the problems of all 6 node candidates
# here the branch is the new branch next to the node
# This encodes one of MVSB/MVWB/RBEB
function find_branching_entry(branch::SpatialBCBranch, node::BB.AbstractNode)::Tuple{Int64, Int64}
    root = find_root(node)
    model = root.model
    opt_idx = solve_candidate_nodes(branch, node)
    return opt_idx
end

function solve_candidate_nodes(branch::SpatialBCBranch, node::BB.AbstractNode)
    function _make_fake_branch(branch::SpatialBCBranch, mode::String)
        (Lii, Uii, Ljj, Ujj, Lij, Uij) = branch.bounds
        if mode == "Wii_up" Lii = (Lii + Uii) / 2
        elseif mode == "Wii_down" Uii = (Lii + Uii) / 2
        elseif mode == "Wjj_up" Ljj = (Ljj + Ujj) / 2
        elseif mode == "Wjj_down" Ujj = (Ljj + Ujj) / 2
        elseif mode == "Wij_up" Lij = (Lij + Uij) / 2
        elseif mode == "Wij_down" Uij = (Lij + Uij) / 2
        end
        new_bounds = (Lii, Uii, Ljj, Ujj, Lij, Uij)
        new_π = compute_π(new_bounds)
        return SpatialBCBranch(branch.i, branch.j, branch.wii, branch.wjj, branch.wr, branch.wi, nothing, 
                                new_bounds, new_π )
    end

    function compute_score(up_solution::Dict{String, <:Real}, down_solution::Dict{String, <:Real})::Float64
        μ = 0.15
        λmin_up = 0.5 * (up_solution["Wii"] - up_solution["Wjj"] - norm([up_solution["Wii"] - up_solution["Wjj"], 2 * up_solution["Wij"], 2 * up_solution["Tij"]]) )
        λmin_down = 0.5 * (down_solution["Wii"] - down_solution["Wjj"] - norm([down_solution["Wii"] - down_solution["Wjj"], 2 * down_solution["Wij"], 2 * down_solution["Tij"]]) )
        return μ * max(-λmin_up, -λmin_down) + (1-μ) * min(-λmin_up, -λmin_down)
    end

    function delete_prev_branch_constr!(model::JuMP.Model, node::BB.AbstractNode)
        for cref in node.auxiliary_data["prev_branch_crefs"]
            delete(model, cref)
        end
    end

    i = branch.i
    j = branch.j
    root = find_root(node)
    model = root.model
    backtracking!(model, node)
    modes = ["Wii", "Wjj", "Wij"]
    best_score = -Inf
    best_mode = ""
    for mode in modes
        fake_branch_up = _make_fake_branch(branch, mode * "_up")
        local_constrs = add_constraints_from_branch!(model, fake_branch_up)
        optimize!(model)
        up_solution = Dict("Wii" => JuMP.value(fake_branch_up.wii), "Wjj" => JuMP.value(fake_branch_up.wjj), "Wij" => JuMP.value(fake_branch_up.wr), "Tij" => JuMP.value(fake_branch_up.wi))
        delete.(model, local_constrs)

        fake_branch_down = _make_fake_branch(branch, mode * "_down")
        local_constrs = add_constraints_from_branch!(model, fake_branch_down)
        optimize!(model)
        down_solution = Dict("Wii" => JuMP.value(fake_branch_down.wii), "Wjj" => JuMP.value(fake_branch_down.wjj), "Wij" => JuMP.value(fake_branch_down.wr), "Tij" => JuMP.value(fake_branch_down.wi))
        delete.(model, local_constrs)

        # compute score, update if better
        score = compute_score(up_solution, down_solution)
        if score > best_score
            best_score = score
            best_mode = mode
        end
    end

    delete_prev_branch_constr!(model, node)
    if best_mode == "Wii"
        return (i,i)
    elseif best_mode == "Wjj"
        return (j,j)
    else # "Wij"
        return (i,j)
    end
end

function backtracking!(model::JuMP.Model, node::BB.AbstractNode)
    node.auxiliary_data["prev_branch_crefs"] = ConstraintRef[]
    pnode = node
    while !isnothing(pnode.parent)
        crefs = add_constraints_from_branch!(model, pnode.branch)
        for cref in crefs
            push!(node.auxiliary_data["prev_branch_crefs"], cref)
        end
        pnode = pnode.parent
    end
end

function add_constraints_from_branch!(model::JuMP.Model, branch::SpatialBCBranch)::Vector{JuMP.ConstraintRef}
    i = branch.i
    j = branch.j
    (Lii, Uii, Ljj, Ujj, Lij, Uij) = branch.bounds
    πs = branch.valid_ineq_coeffs
    wii = branch.wii
    wjj = branch.wjj
    wr = branch.wr
    wi = branch.wi
    new_branches = []
    push!(new_branches, JuMP.@constraint(node.model, Lii <= wii <= Uii))
    push!(new_branches, JuMP.@constraint(node.model, Ljj <= wjj <= Ujj))
    push!(new_branches, JuMP.@constraint(node.model, Lij * wr <= wi))
    push!(new_branches, JuMP.@constraint(node.model, wi <= Uij * wi))
    push!(new_branches, JuMP.@constraint(node.model, πs[1] + πs[2] * wii + πs[3] * wjj + πs[4] * wr + πs[5] * wi >= Ujj * wii + Uii * wjj - Uii * Ujj))
    push!(new_branches, JuMP.@constraint(node.model, πs[1] + πs[2] * wii + πs[3] * wjj + πs[4] * wr + πs[5] * wi >= Ljj * wii + Lii * wjj - Lii * Ljj))
    return new_branches
end

function BB.branch!(tree::BB.AbstractTree, node::BB.AbstractNode)
    @info " Node id $(node.id), status $(node.solution_status), bound $(node.bound)"
    if node.bound >= tree.best_incumbent
        @info " Fathomed by bound"
    elseif node.solution_status in [MOI.OPTIMAL, MOI.LOCALLY_SOLVED, MOI.SLOW_PROGRESS]
        (i,j) = find_min_eigen(node)
        new_sbc_branch = create_sbc_branch(i, j, node)
        (new_i,new_j) = find_branching_entry(new_sbc_branch, node) # This function implements MVSB/MVWB/RBEB
        up_bounds_arr = [i for i in new_sbc_branch.bounds]
        down_bounds_arr = [i for i in new_sbc_branch.bounds]
        if new_i != new_j # branch on Wij
            up_bounds_arr[5] = (up_bounds_arr[5] + up_bounds_arr[6]) / 2
            down_bounds_arr[6] = (down_bounds_arr[5] + down_bounds_arr[6]) / 2
        elseif new_i == i # branch on Wii
            up_bounds_arr[1] = (up_bounds_arr[1] + up_bounds_arr[2]) / 2
            down_bounds_arr[2] = (down_bounds_arr[1] + down_bounds_arr[2]) / 2
        else # branch on Wjj
            up_bounds_arr[3] = (up_bounds_arr[3] + up_bounds_arr[4]) / 2
            down_bounds_arr[4] = (down_bounds_arr[3] + down_bounds_arr[4]) / 2
        end
        next_branch_up = deepcopy(new_sbc_branch)
        next_branch_down = new_sbc_branch
        next_branch_up.bounds = Tuple(up_bounds_arr)
        next_branch_down.bounds = Tuple(down_bounds_arr)
        child_up = BB.create_child_node(node, next_branch_up)
        child_down = BB.create_child_node(node, next_branch_down)
        BB.push!(tree, child_up)
        BB.push!(tree, child_down)
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

# This makes original adjust_branch! invalid
function BB.adjust_branch!(branch_objects::Array{SpatialBCBranch,1})
    return
end

# This makes original apply_changes! invalid
function BB.apply_changes!(node::BB.JuMPNode)
    return
end

function BB.bound!(node::BB.JuMPNode)
    node.model = find_root(node).model
    backtracking!(node.model, node)
    JuMP.optimize!(node.model)
    node.solution_status = JuMP.termination_status(node.model)

    if node.solution_status in [MOI.OPTIMAL, MOI.LOCALLY_SOLVED, MOI.SLOW_PROGRESS]
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
end

###########################################################################
#                                Main Code                                #
###########################################################################

# read data, and formulate root model using PowerModels
file = "/home/weiqizhang/.julia/dev/BranchAndBound/data/case9.m"
data = parse_file(file)
pm = instantiate_model(data, NodeWRMPowerModel, build_opf)
optimizer = optimizer_with_attributes(Mosek.Optimizer, "QUIET" => true)
# optimizer = optimizer_with_attributes(Ipopt.Optimizer, "print_level" => 0)
set_optimizer(pm.model, optimizer)

# collect data
Lii = Dict(i => bus["vmin"]^2 for (i, bus) in ref(pm, :bus))
Uii = Dict(i => bus["vmax"]^2 for (i, bus) in ref(pm, :bus))
Lij = Dict((i,j) => tan(branch["angmin"]) for ((i,j), branch) in ref(pm, :buspairs))
Uij = Dict((i,j) => tan(branch["angmax"]) for ((i,j), branch) in ref(pm, :buspairs))

# initialize branch-and-cut tree
node = BB.JuMPNode{SpatialBCBranch}(pm.model)
node.auxiliary_data["PM"] = pm
node.auxiliary_data["Lii"] = Lii
node.auxiliary_data["Uii"] = Uii
node.auxiliary_data["Lij"] = Lij
node.auxiliary_data["Uij"] = Uij
tree = BB.initialize_tree(node)

BB.run(tree)