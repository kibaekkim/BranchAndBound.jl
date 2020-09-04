using BranchAndBound, JuMP, PowerModels, DeNet
using Mosek, MosekTools, Ipopt, SCS
using LinearAlgebra
const BB = BranchAndBound
const PM = PowerModels

const MOSEK_OPTIMIZER = optimizer_with_attributes(Mosek.Optimizer, "QUIET" => true)
const SCS_OPTIMIZER = optimizer_with_attributes(SCS.Optimizer, "max_iters" => 100000, "eps" => 1e-4, "verbose" => 0)
const MAX_NODE_ALLOWED = 10000

###########################################################################
#                             Data Structures                             #
###########################################################################

# A new type specifically for spatial branch and bound, essentially SDPWRMPowerModel
mutable struct NodeWRMPowerModel <: AbstractWRMModel @pm_fields end

# A new branch type, mapping (Wij, Tij) to their lower and upper bounds
mutable struct ComplexVariableBranch <: BB.AbstractBranch
    lb::Dict{Tuple{JuMP.VariableRef, JuMP.VariableRef},Real} # Lij
    ub::Dict{Tuple{JuMP.VariableRef, JuMP.VariableRef},Real} # Uij
end

# This computes coefficients for cuts (3a), (3b) in the paper
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

# This is the main branch type used in our implementation
# (i,j) should always be one of the lines in the power system
mutable struct SpatialBCBranch <: BB.AbstractBranch
    i::Int
    j::Int
    wii::JuMP.VariableRef
    wjj::JuMP.VariableRef
    wr::JuMP.VariableRef
    wi::JuMP.VariableRef
    mod_branch::Union{Nothing, ComplexVariableBranch, BB.VariableBranch} # this captures either a changed bound for Wij and Tij, or a changed bound for Wii/Wjj 
    bounds::NTuple{6, <:Real} # Lii, Uii, Ljj, Ujj, Lij, Uij
    valid_ineq_coeffs::NTuple{5, <:Real} # π coefficients
end

# create a SpatialBCBranch, in which mod_branch is nothing (because the exact matrix entry to branch is not decided yet)
function create_sbc_branch(i::Int64, j::Int64, prev_node::BB.AbstractNode)
    root = find_root(prev_node)
    pm = root.auxiliary_data["PM"]
    bus_ids = ids(pm, :bus)
    lookup_w_index = Dict((bi,i) for (i,bi) in enumerate(bus_ids))
    widx_i = lookup_w_index[i]
    widx_j = lookup_w_index[j]
    wii = var(pm, :WR)[widx_i,widx_i]
    wjj = var(pm, :WR)[widx_j,widx_j]
    wr = var(pm, :WR)[widx_i,widx_j]
    wi = var(pm, :WI)[widx_i,widx_j]
    bounds = get_LU_from_branches(prev_node, i, j)
    πs = compute_π(bounds)
    return SpatialBCBranch(i, j, wii, wjj, wr, wi, nothing, bounds, πs)
end

# This function tracks the branches and update the bounds for selected entries (i,j)
function get_LU_from_branches(node::BB.AbstractNode, i::Int64, j::Int64)::NTuple{6, <:Real}
    pnode = node
    LU = [NaN for i in 1:6]
    while !isnothing(pnode.parent)
        sbc_branch = pnode.branch
        mod_branch = sbc_branch.mod_branch
        if sbc_branch.i == i && sbc_branch.j == j
            for i in eachindex(LU)
                if isnan(LU[i]) LU[i] = sbc_branch.bounds[i] end
            end
        elseif sbc_branch.i == i
            if isnan(LU[1]) LU[1] = sbc_branch.bounds[1] end
            if isnan(LU[2]) LU[2] = sbc_branch.bounds[2] end
        elseif sbc_branch.i == j
            if isnan(LU[3]) LU[3] = sbc_branch.bounds[1] end
            if isnan(LU[4]) LU[4] = sbc_branch.bounds[2] end
        elseif sbc_branch.j == i
            if isnan(LU[1]) LU[1] = sbc_branch.bounds[3] end
            if isnan(LU[2]) LU[2] = sbc_branch.bounds[4] end
        elseif sbc_branch.j == j
            if isnan(LU[3]) LU[3] = sbc_branch.bounds[3] end
            if isnan(LU[4]) LU[4] = sbc_branch.bounds[4] end
        end
        pnode = pnode.parent
    end
    if isnan(LU[1]) LU[1] = pnode.auxiliary_data["Lii"][i] end
    if isnan(LU[2]) LU[2] = pnode.auxiliary_data["Uii"][i] end
    if isnan(LU[3]) LU[3] = pnode.auxiliary_data["Lii"][j] end
    if isnan(LU[4]) LU[4] = pnode.auxiliary_data["Uii"][j] end
    if isnan(LU[5]) LU[5] = pnode.auxiliary_data["Lij"][(i,j)] end
    if isnan(LU[6]) LU[6] = pnode.auxiliary_data["Uij"][(i,j)] end

    return Tuple(LU)
end
###########################################################################
#                          PowerModels Extensions                         #
###########################################################################

# get rid of angle difference bounds, which will be included in other parts of code
# the valid inequality cuts for voltage products are still kept
function PM.constraint_voltage_angle_difference(pm::AbstractWModels, n::Int, f_idx, angmin, angmax)
    i, f_bus, t_bus = f_idx

    w_fr = var(pm, n, :w, f_bus)
    w_to = var(pm, n, :w, t_bus)
    wr   = var(pm, n, :wr, (f_bus, t_bus))
    wi   = var(pm, n, :wi, (f_bus, t_bus))

    cut_complex_product_and_angle_difference(pm.model, w_fr, w_to, wr, wi, angmin, angmax)
end

# This is essentially copied from original PM.build_opf 
function PM.build_opf(pm::NodeWRMPowerModel)
    variable_bus_voltage(pm)
    variable_gen_power(pm)
    variable_branch_power(pm)
    variable_dcline_power(pm)

    DeNet.objective_min_fuel_and_flow_cost_mod(pm)

    constraint_model_voltage(pm)

    :cut_bus in keys(ref(pm)) ? cut_bus = ids(pm, :cut_bus) : cut_bus = [] # This is needed in order to be compatible with network decomposition later
    for i in setdiff(ids(pm, :bus), cut_bus)
        constraint_power_balance(pm, i)
    end

    for i in ids(pm, :branch)
        constraint_ohms_yt_from(pm, i)
        constraint_ohms_yt_to(pm, i)
        constraint_voltage_angle_difference(pm, i)
        constraint_thermal_limit_from(pm, i)
        constraint_thermal_limit_to(pm, i)
    end

    for i in ids(pm, :dcline)
        constraint_dcline_power_losses(pm, i)
    end
end

###########################################################################
#                            Helper Functions                             #
###########################################################################

# This finds the root of the tree from any node in the tree
function find_root(node::BB.AbstractNode)::BB.AbstractNode
    pnode = node
    while !isnothing(pnode.parent)
        pnode = pnode.parent
    end
    return pnode
end

###########################################################################
#                        Key Steps Implementation                         #
###########################################################################

# This function finds the line index (i,j) that the corresponding 2-by-2 W matrix attains the maximum λmin
function find_min_eigen(node::BB.AbstractNode)::Tuple{Int64, Int64}
    pm = find_root(node).auxiliary_data["PM"]
    bus_ids = ids(pm, :bus)
    lookup_w_index = Dict((bi,i) for (i,bi) in enumerate(bus_ids))
    wr = var(pm, :WR)
    wi = var(pm, :WI)
    max_lambda = -Inf
    max_id = ()
    node_solution = node.solution
    eigenvalues = Dict()
    # for raw_k in bus_ids, raw_l in bus_ids
    #     k = lookup_w_index[raw_k]
    #     l = lookup_w_index[raw_l]
    #     if k < l
    #         lambda = 0.5 * (node_solution[wr[k,k]] + node_solution[wr[l,l]] - norm( [node_solution[wr[k,k]] - node_solution[wr[l,l]], 2 * node_solution[wr[k,l]], 2 * node_solution[wi[k,l]]] ) )
    #         println("Smallest eigenvalue for index pair ($(k), $(l)): ", lambda)
    #     end
    # end
    for (raw_i,raw_j) in ids(pm, :buspairs)
        i = lookup_w_index[raw_i]
        j = lookup_w_index[raw_j]
        lambda = 0.5 * (node_solution[wr[i,i]] + node_solution[wr[j,j]] - norm( [node_solution[wr[i,i]] - node_solution[wr[j,j]], 2 * node_solution[wr[i,j]], 2 * node_solution[wi[i,j]]] ) )
        eigenvalues[(raw_i, raw_j)] = lambda
        # println(lambda)
        if lambda > max_lambda
            max_id = (raw_i,raw_j)
            max_lambda = lambda
        end
    end
    node.auxiliary_data["eigenvalues"] = eigenvalues
    # println()
    # println(max_lambda)
    return max_id
end

# this is the interface for finding the complex entry to branch
function find_branching_entry(branch::SpatialBCBranch, node::BB.AbstractNode)::Tuple{Int64, Int64}
    root = find_root(node)
    model = root.model
    opt_idx = solve_candidate_nodes(branch, node)
    return opt_idx
end

# this implements weak branching to find the best, and returns the index of the complex entry to branch on
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
    
    function _add_constraints_from_fake_branch!(model::JuMP.Model, branch::SpatialBCBranch)
        i = branch.i
        j = branch.j
        (Lii, Uii, Ljj, Ujj, Lij, Uij) = branch.bounds
        πs = branch.valid_ineq_coeffs
        wii = branch.wii
        wjj = branch.wjj
        wr = branch.wr
        wi = branch.wi
        new_branches = []
        JuMP.set_lower_bound(wii, Lii)
        JuMP.set_upper_bound(wii, Uii)
        JuMP.set_lower_bound(wjj, Ljj)
        JuMP.set_upper_bound(wjj, Ujj)
        JuMP.@constraint(model, Lij * wr <= wi)
        JuMP.@constraint(model, Uij * wr >= wi)
        JuMP.@constraint(model, πs[1] + πs[2] * wii + πs[3] * wjj + πs[4] * wr + πs[5] * wi >= Ujj * wii + Uii * wjj - Uii * Ujj)
        JuMP.@constraint(model, πs[1] + πs[2] * wii + πs[3] * wjj + πs[4] * wr + πs[5] * wi >= Ljj * wii + Lii * wjj - Lii * Ljj)
    end    
    i = branch.i
    j = branch.j
    root = find_root(node)
    modes = ["Wii", "Wjj", "Wij"]
    best_score = -Inf
    best_mode = ""
    for mode in modes
        fake_branch = _make_fake_branch(branch, mode * "_up")
        model = JuMP.Model(MOSEK_OPTIMIZER)
        # model = JuMP.Model(SCS_OPTIMIZER)
        fake_branch.wii = JuMP.@variable(model, wii)
        fake_branch.wjj = JuMP.@variable(model, wjj)
        fake_branch.wr = JuMP.@variable(model, wr)
        fake_branch.wi = JuMP.@variable(model, wi)
        JuMP.@variable(model, λ)
        JuMP.@constraint(model, [wii + wjj - 2 * λ, wii - wjj, 2 * wr, 2 * wi] in SecondOrderCone())
        JuMP.@objective(model, Max, λ)
        _add_constraints_from_fake_branch!(model, fake_branch)
        optimize!(model)
        λ_up = JuMP.value(λ)

        fake_branch = _make_fake_branch(branch, mode * "_down")
        model = JuMP.Model(MOSEK_OPTIMIZER)
        # model = JuMP.Model(SCS_OPTIMIZER)
        fake_branch.wii = JuMP.@variable(model, wii)
        fake_branch.wjj = JuMP.@variable(model, wjj)
        fake_branch.wr = JuMP.@variable(model, wr)
        fake_branch.wi = JuMP.@variable(model, wi)
        JuMP.@variable(model, λ)
        JuMP.@constraint(model, [wii + wjj - 2 * λ, wii - wjj, 2 * wr, 2 * wi] in SecondOrderCone())
        JuMP.@objective(model, Max, λ)
        _add_constraints_from_fake_branch!(model, fake_branch)
        optimize!(model)
        λ_down = JuMP.value(λ)

        # compute score, update if better
        μ = 0.15
        score = μ * max(-λ_up, -λ_down) + (1-μ) * min(-λ_up, -λ_down)
        if score > best_score
            best_score = score
            best_mode = mode
        end
    end

    if best_mode == "Wii"
        return (i,i)
    elseif best_mode == "Wjj"
        return (j,j)
    else # "Wij"
        return (i,j)
    end
end

#=
# This implements strong branching, might be outdated

function solve_candidate_nodes(branch::SpatialBCBranch, node::BB.AbstractNode)
    function compute_score(up_solution::Dict{String, <:Real}, down_solution::Dict{String, <:Real})::Float64
        μ = 0.15
        λmin_up = 0.5 * (up_solution["Wii"] - up_solution["Wjj"] - norm([up_solution["Wii"] - up_solution["Wjj"], 2 * up_solution["Wij"], 2 * up_solution["Tij"]]) )
        λmin_down = 0.5 * (down_solution["Wii"] - down_solution["Wjj"] - norm([down_solution["Wii"] - down_solution["Wjj"], 2 * down_solution["Wij"], 2 * down_solution["Tij"]]) )
        return μ * max(-λmin_up, -λmin_down) + (1-μ) * min(-λmin_up, -λmin_down)
    end
    i = branch.i
    j = branch.j
    root = find_root(node)
    model = root.model
    # backtracking!(model, node)
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
        println(score)
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
=#

# this tracks all the branches above node and add constraints to model accordingly
# if there are multiple branches on the same complex entry, only add constraints for the one with largest depth (closest to node)
function backtracking!(model::JuMP.Model, node::BB.AbstractNode)
    node.auxiliary_data["prev_branch_crefs"] = ConstraintRef[]
    root = find_root(node)
    pnode = node
    already_modified = []
    while !isnothing(pnode.parent)
        (i,j) = (pnode.branch.i, pnode.branch.j)
        if !((i,j) in already_modified)
            crefs = add_constraints_from_branch!(model, pnode.branch, root)
            for cref in crefs
                push!(node.auxiliary_data["prev_branch_crefs"], cref)
            end
            push!(already_modified, (i,j))
        end
        pnode = pnode.parent
    end
end

# this modifies constraints in model based on branch
function add_constraints_from_branch!(model::JuMP.Model, branch::SpatialBCBranch, root::BB.AbstractNode)::Vector{JuMP.ConstraintRef}
    i = branch.i
    j = branch.j
    (Lii, Uii, Ljj, Ujj, Lij, Uij) = branch.bounds
    πs = branch.valid_ineq_coeffs
    wii = branch.wii
    wjj = branch.wjj
    wr = branch.wr
    wi = branch.wi
    new_branches = []
    JuMP.set_lower_bound(wii, Lii)
    JuMP.set_upper_bound(wii, Uii)
    JuMP.set_lower_bound(wjj, Ljj)
    JuMP.set_upper_bound(wjj, Ujj)
    JuMP.set_normalized_coefficient(root.auxiliary_data["Cuts"]["angle_lb"][(i,j)], wr, Lij)
    JuMP.set_normalized_coefficient(root.auxiliary_data["Cuts"]["angle_ub"][(i,j)], wr, Uij)
    push!(new_branches, JuMP.@constraint(model, πs[1] + πs[2] * wii + πs[3] * wjj + πs[4] * wr + πs[5] * wi >= Ujj * wii + Uii * wjj - Uii * Ujj))
    push!(new_branches, JuMP.@constraint(model, πs[1] + πs[2] * wii + πs[3] * wjj + πs[4] * wr + πs[5] * wi >= Ljj * wii + Lii * wjj - Lii * Ljj))
    return new_branches
end

# this resets constraints (2) and deletes constraints (3) in model
function delete_prev_branch_constr!(model::JuMP.Model, node::BB.AbstractNode)
    root = find_root(node)
    pm = root.auxiliary_data["PM"]
    w = var(pm, :w)
    wr = var(pm, :wr)
    for (i,_) in ref(pm, :bus)
        JuMP.set_lower_bound(w[i], root.auxiliary_data["Lii"][i])
        JuMP.set_upper_bound(w[i], root.auxiliary_data["Uii"][i])
    end
    for (pair,_) in ref(pm, :buspairs)
        JuMP.set_normalized_coefficient(root.auxiliary_data["Cuts"]["angle_lb"][pair], wr[pair], root.auxiliary_data["Lij"][pair])
        JuMP.set_normalized_coefficient(root.auxiliary_data["Cuts"]["angle_ub"][pair], wr[pair], root.auxiliary_data["Uij"][pair])
    end
    while !isempty(node.auxiliary_data["prev_branch_crefs"])
        cref = Base.pop!(node.auxiliary_data["prev_branch_crefs"])
        delete(model, cref)
    end
end

###########################################################################
#                 Extensions of BranchAndBound Functions                  #
#                   General Algorithm Flow Goes here                      #
###########################################################################

# return a "deepcopy" of SpatialBCBranch for everything except for variable references
# this is needed to ensure variable references in new branch copy are still valid
# the mod_branch will not be duplicated here if it is not nothing
function branch_copy(branch::SpatialBCBranch)
    return SpatialBCBranch(branch.i, branch.j, branch.wii, branch.wjj, branch.wr, branch.wi, 
                           branch.mod_branch, deepcopy(branch.bounds), deepcopy(branch.valid_ineq_coeffs))
end

function BB.branch!(tree::BB.AbstractTree, node::BB.AbstractNode)
    # @info " Node id $(node.id), status $(node.solution_status), bound $(node.bound)"
    root = find_root(node)
    model = root.model
    if node.bound >= tree.best_incumbent
        if isapprox(node.bound, tree.best_incumbent, rtol = 1e-5) && node.bound > root.auxiliary_data["best_bound"]
            root.auxiliary_data["best_bound"] = node.bound
            root.auxiliary_data["best_id"] = node.id
        end
        # @info " Fathomed by bound"
    elseif node.depth >= 100
        if node.bound > root.auxiliary_data["best_bound"]
            root.auxiliary_data["best_bound"] = node.bound
            root.auxiliary_data["best_id"] = node.id
        end
        # @info " Fathomed by maximum depth"
    elseif node.solution_status in [MOI.OPTIMAL, MOI.LOCALLY_SOLVED, MOI.SLOW_PROGRESS]
        # determine the complex entry to branch on based on solution
        if node.bound > root.auxiliary_data["best_bound"]
            root.auxiliary_data["best_bound"] = node.bound
            root.auxiliary_data["best_id"] = node.id
        end
        pm = root.auxiliary_data["PM"]
        (i,j) = find_min_eigen(node)
        WR = var(pm, :WR)
        WI = var(pm, :WI)
        # w = value.([WR -WI; WI WR])
        w = value.(WR) .+ im * value.(WI)
        println(eigvals(w))
        if node.auxiliary_data["eigenvalues"][(i,j)] <= 1e-5
            # update the node with best bound (including processed ones), recorded in root    
            @info " Fathomed by reaching rank-1 solution"
        else
            new_sbc_branch = create_sbc_branch(i, j, node)
            (new_i,new_j) = find_branching_entry(new_sbc_branch, node)
    
            # create branches and child nodes accordingly
            up_bounds_arr = [k for k in new_sbc_branch.bounds]
            down_bounds_arr = [k for k in new_sbc_branch.bounds]
            bus_ids = ids(pm, :bus)
            lookup_w_index = Dict((bi,i) for (i,bi) in enumerate(bus_ids))
            widx_new_i = lookup_w_index[new_i]
            widx_new_j = lookup_w_index[new_j]
            if new_i != new_j # branch on Wij
                up_bounds_arr[5] = (up_bounds_arr[5] + up_bounds_arr[6]) / 2
                down_bounds_arr[6] = (down_bounds_arr[5] + down_bounds_arr[6]) / 2
                vpair = (PM.var(pm, :WR)[widx_new_i,widx_new_j], PM.var(pm, :WI)[widx_new_i,widx_new_j])
                up_mod_branch = ComplexVariableBranch(Dict(vpair => up_bounds_arr[5]), Dict(vpair => up_bounds_arr[6]))
                down_mod_branch = ComplexVariableBranch(Dict(vpair => down_bounds_arr[5]), Dict(vpair => down_bounds_arr[6]))
                # @info " Branch at W$(new_i)$(new_j), [L, U] breaks into by [$(down_bounds_arr[5]),$(up_bounds_arr[5]),$(up_bounds_arr[6])]."
            elseif new_i == i # branch on Wii
                up_bounds_arr[1] = (up_bounds_arr[1] + up_bounds_arr[2]) / 2
                down_bounds_arr[2] = (down_bounds_arr[1] + down_bounds_arr[2]) / 2
                v = PM.var(pm, :WR)[widx_new_i, widx_new_j]
                up_mod_branch = BB.VariableBranch(Dict(v => up_bounds_arr[1]), Dict(v => up_bounds_arr[2]))
                down_mod_branch = BB.VariableBranch(Dict(v => down_bounds_arr[1]), Dict(v => down_bounds_arr[2]))
                # @info " Branch at W$(new_i)$(new_j), [L, U] breaks into by [$(down_bounds_arr[1]),$(up_bounds_arr[1]),$(up_bounds_arr[2])]."
            else # branch on Wjj
                up_bounds_arr[3] = (up_bounds_arr[3] + up_bounds_arr[4]) / 2
                down_bounds_arr[4] = (down_bounds_arr[3] + down_bounds_arr[4]) / 2
                v = PM.var(pm, :WR)[widx_new_i, widx_new_j]
                up_mod_branch = BB.VariableBranch(Dict(v => up_bounds_arr[3]), Dict(v => up_bounds_arr[4]))
                down_mod_branch = BB.VariableBranch(Dict(v => down_bounds_arr[3]), Dict(v => down_bounds_arr[4]))
                # @info " Branch at W$(new_i)$(new_j), [L, U] breaks into by [$(down_bounds_arr[3]),$(up_bounds_arr[3]),$(up_bounds_arr[4])]."
            end
            next_branch_up = branch_copy(new_sbc_branch)
            next_branch_down = new_sbc_branch
            next_branch_up.bounds = Tuple(up_bounds_arr)
            next_branch_down.bounds = Tuple(down_bounds_arr)
            next_branch_up.mod_branch = up_mod_branch
            next_branch_down.mod_branch = down_mod_branch
            child_up = BB.create_child_node(node, next_branch_up)
            child_down = BB.create_child_node(node, next_branch_down)
            BB.push!(tree, child_up)
            BB.push!(tree, child_down)    
        end
    else
        # @info " Fathomed by solution status: $(node.solution_status)"
    end
    push!(root.auxiliary_data["best_bounds"], root.auxiliary_data["best_bound"])
end

# implement depth first rule
function BB.next_node(tree::BB.AbstractTree)
    # sort! is redefined to find the node with maximum depth (consistent with the paper's implementation)
    # sort!(tree::BB.AbstractTree) = Base.sort!(tree.nodes, by=x->x.depth)
    BB.sort!(tree)
    node = Base.pop!(tree.nodes)
    return node
end

function BB.termination(tree::BB.AbstractTree)
    # @info "Tree nodes: processed $(length(tree.processed)), left $(length(tree.nodes)), total $(tree.node_counter), best bound $(tree.best_bound), best incumbent $(tree.best_incumbent)"
    if BB.isempty(tree)
        # @info "Completed the tree search"
        return true
    end
    if length(tree.processed) >= MAX_NODE_ALLOWED
        # @info "Reached node limit"
        return true
    end
    return false
end

# This makes original adjust_branch! invalid
BB.adjust_branch!(branch_objects::Array{SpatialBCBranch,1}) = nothing

# This makes original apply_changes! invalid
BB.apply_changes!(node::BB.JuMPNode) = nothing

function BB.bound!(node::BB.JuMPNode)
    model = find_root(node).model
    backtracking!(model, node)
    JuMP.optimize!(model)
    node.solution_status = JuMP.termination_status(model)
    if node.solution_status == MOI.INFEASIBLE || JuMP.dual_status(model) in [MOI.INFEASIBILITY_CERTIFICATE]
        node.bound = Inf
    elseif node.solution_status == MOI.DUAL_INFEASIBLE || JuMP.primal_status(model) in [MOI.INFEASIBILITY_CERTIFICATE]
        node.bound = -Inf
    elseif node.solution_status in [MOI.OPTIMAL, MOI.LOCALLY_SOLVED, MOI.SLOW_PROGRESS]
        node.bound = JuMP.objective_value(model)
        vrefs = JuMP.all_variables(model)
        for v in vrefs
            node.solution[v] = JuMP.value(v)
        end
    else
        @warn "Unexpected node solution status: $(node.solution_status)"
        node.bound = -Inf
    end
    delete_prev_branch_constr!(model, node)
end

#=
# This is supposed to find an incumbent for each node, by converting SDP problem into rank-1 form and solving with Ipopt
# This has not been tested

function BB.heuristics!(node::BB.JuMPNode)
    ipopt_optimizer = optimizer_with_attributes(Ipopt.Optimizer, "print_level" => 0)
    exact_model = copy(node.model)
    buses = collect(keys(ref(pm, :bus)))
    JuMP.@variable(exact_model, vr[buses])
    JuMP.@variable(exact_model, vi[buses])

    # get rid of PSDCone constraint (for Ipopt)
    psdcone = JuMP.all_constraints(exact_model, Array{GenericAffExpr{Float64, VariableRef}, 1}, MOI.PositiveSemidefiniteConeSquare)[1]
    JuMP.delete(exact_model, psdcone)

    # add rank 1 constraints
    for i in buses, j in buses
        if i == j
            Wi = JuMP.variable_by_name(exact_model, "0_WR[$(i),$(i)]")
            Wj = JuMP.variable_by_name(exact_model, "0_WR[$(j),$(j)]")
            JuMP.@constraint(exact_model, Wi == vr[i]^2 + vi[i]^2)
            JuMP.@constraint(exact_model, Wj == vr[j]^2 + vi[j]^2)
        else
            if !isnothing(JuMP.variable_by_name(exact_model, "0_WR[$(i),$(j)]"))
                WR = JuMP.variable_by_name(exact_model, "0_WR[$(i),$(j)]")
            else
                WR = JuMP.variable_by_name(exact_model, "0_WR[$(j),$(i)]")
            end
            WI = JuMP.variable_by_name(exact_model, "0_WI[$(i),$(j)]")
            JuMP.@constraint(exact_model, WR == vr[i] * vr[j] + vi[i] * vi[j])
            JuMP.@constraint(exact_model, WI == vr[j] * vi[i] - vi[j] * vr[i])
        end
    end

    # convert second order cone constraints into quadratic constraints
    soc_crefs = JuMP.all_constraints(exact_model, Array{GenericAffExpr{Float64, VariableRef}, 1}, MOI.SecondOrderCone)
    for cref in soc_crefs
        cobj = JuMP.constraint_object(cref)
        terms = cobj.func
        JuMP.@constraint(exact_model, sum(terms[i]^2 for i in 1:length(terms)) <= terms[1]^2)
        JuMP.delete(exact_model, cref)
    end

    JuMP.set_optimizer(exact_model, ipopt_optimizer)
    JuMP.optimize!(exact_model)
    node.best_incumbent = JuMP.objective_value(exact_model)
    JuMP.set_optimizer(node.model, optimizer)
end
=#

###########################################################################
#                     Helper Functions for Debugging                      #
###########################################################################
#=
function total_num_constraints(model::JuMP.Model)::Int64
    sum = 0
    for (i,j) in JuMP.list_of_constraint_types(model)
        sum += num_constraints(model, i, j)
    end
    return sum
end

function print_model_to_file(model::JuMP.Model, id::String)
    f = open(id, "w")
    print(f, model)
    close(f)
end
=#

###########################################################################
#                          Initialization Method                          #
###########################################################################

function initialize(pm::NodeWRMPowerModel)::Tuple{BB.AbstractTree, BB.AbstractNode}
    # collect data
    Lii = Dict(i => bus["vmin"]^2 for (i, bus) in ref(pm, :bus))
    Uii = Dict(i => bus["vmax"]^2 for (i, bus) in ref(pm, :bus))
    Lij = Dict((i,j) => tan(branch["angmin"]) for ((i,j), branch) in ref(pm, :buspairs))
    Uij = Dict((i,j) => tan(branch["angmax"]) for ((i,j), branch) in ref(pm, :buspairs))

    # initialize branch-and-cut tree
    node = BB.JuMPNode{SpatialBCBranch}(pm.model)
    node.auxiliary_data["best_id"] = 0
    node.auxiliary_data["best_bound"] = -Inf
    node.auxiliary_data["PM"] = pm
    node.auxiliary_data["Lii"] = Lii
    node.auxiliary_data["Uii"] = Uii
    node.auxiliary_data["Lij"] = Lij
    node.auxiliary_data["Uij"] = Uij
    node.auxiliary_data["Cuts"] = Dict( "angle_lb" => Dict() , "angle_ub" => Dict())
    wr = var(pm, :wr)
    wi = var(pm, :wi)
    for (i, branch) in ref(pm, :branch)
        f_bus = branch["f_bus"]
        t_bus = branch["t_bus"]
        node.auxiliary_data["Cuts"]["angle_lb"][(f_bus, t_bus)] = @constraint(pm.model, Lij[(f_bus, t_bus)] * wr[(f_bus, t_bus)] <= wi[(f_bus, t_bus)])
        node.auxiliary_data["Cuts"]["angle_ub"][(f_bus, t_bus)] = @constraint(pm.model, Uij[(f_bus, t_bus)] * wr[(f_bus, t_bus)] >= wi[(f_bus, t_bus)])
    end
    node.auxiliary_data["best_bounds"] = []

    tree = BB.initialize_tree(node)

    # set incumbent
    ipopt_optimizer = optimizer_with_attributes(Ipopt.Optimizer, "print_level" => 0)
    ipopt_solution = run_opf(file, ACRPowerModel, ipopt_optimizer)
    tree.best_incumbent = ipopt_solution["objective"]

    return tree, node
end
