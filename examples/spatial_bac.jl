using BranchAndBound, JuMP, PowerModels, DeNet
using Mosek, MosekTools, Ipopt
using LinearAlgebra
const BB = BranchAndBound
const PM = PowerModels

const MOSEK_OPTIMIZER = optimizer_with_attributes(Mosek.Optimizer, "QUIET" => true)
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
    valid_ineq_coeffs::NTuple{5, <:Real} # π coefficients
end

# create a SpatialBCBranch, in which mod_branch is Nothing (the exact matrix entry to branch is not decided yet)
# (maybe don't need this... instead decide the branching entry first...)
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

# get rid of angle difference bounds
# will include in other parts of code
function PM.constraint_voltage_angle_difference(pm::AbstractWModels, n::Int, f_idx, angmin, angmax)
    i, f_bus, t_bus = f_idx

    w_fr = var(pm, n, :w, f_bus)
    w_to = var(pm, n, :w, t_bus)
    wr   = var(pm, n, :wr, (f_bus, t_bus))
    wi   = var(pm, n, :wi, (f_bus, t_bus))

    cut_complex_product_and_angle_difference(pm.model, w_fr, w_to, wr, wi, angmin, angmax)
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
        constraint_voltage_angle_difference(pm, i)
        constraint_thermal_limit_from(pm, i)
        constraint_thermal_limit_to(pm, i)
    end

    for i in ids(pm, :dcline)
        constraint_dcline_power_losses(pm, i)
    end
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
    bus_ids = ids(pm, :bus)
    lookup_w_index = Dict((bi,i) for (i,bi) in enumerate(bus_ids))
    wr = var(pm, :WR)
    wi = var(pm, :WI)
    max_lambda = -Inf
    max_id = ()
    node_solution = node.solution
    for ((raw_i,raw_j),_) in ref(pm, :buspairs)
        i = lookup_w_index[raw_i]
        j = lookup_w_index[raw_j]
        lambda = 0.5 * (node_solution[wr[i,i]] - node_solution[wr[j,j]] - norm( [node_solution[wr[i,i]] - node_solution[wr[j,j]], 2 * node_solution[wr[i,j]], 2 * node_solution[wi[i,j]]] ) )
        if lambda > max_lambda
            max_id = (raw_i,raw_j)
            max_lambda = lambda
        end
    end
    return max_id
end

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

function delete_prev_branch_constr!(model::JuMP.Model, node::BB.AbstractNode)
    while !isempty(node.auxiliary_data["prev_branch_crefs"])
        cref = Base.pop!(node.auxiliary_data["prev_branch_crefs"])
        delete(model, cref)
    end
end

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

function solve_candidate_nodes(branch::SpatialBCBranch, node::BB.AbstractNode) # This implements weak branching
    i = branch.i
    j = branch.j
    root = find_root(node)
    model = root.model
    modes = ["Wii", "Wjj", "Wij"]
    best_score = -Inf
    best_mode = ""
    for mode in modes
        fake_branch = _make_fake_branch(branch, mode * "_up")
        model = JuMP.Model(MOSEK_OPTIMIZER)
        fake_branch.wii = JuMP.@variable(model, wii)
        fake_branch.wjj = JuMP.@variable(model, wjj)
        fake_branch.wr = JuMP.@variable(model, wr)
        fake_branch.wi = JuMP.@variable(model, wi)
        JuMP.@variable(model, λ)
        JuMP.@constraint(model, [wii + wjj - 2 * λ, wii - wjj, 2 * wr, 2 * wi] in SecondOrderCone())
        JuMP.@objective(model, Max, λ)
        add_constraints_from_branch!(model, fake_branch)
        optimize!(model)
        λ_up = JuMP.value(λ)

        fake_branch = _make_fake_branch(branch, mode * "_down")
        model = JuMP.Model(MOSEK_OPTIMIZER)
        fake_branch.wii = JuMP.@variable(model, wii)
        fake_branch.wjj = JuMP.@variable(model, wjj)
        fake_branch.wr = JuMP.@variable(model, wr)
        fake_branch.wi = JuMP.@variable(model, wi)
        JuMP.@variable(model, λ)
        JuMP.@constraint(model, [wii + wjj - 2 * λ, wii - wjj, 2 * wr, 2 * wi] in SecondOrderCone())
        JuMP.@objective(model, Max, λ)
        add_constraints_from_branch!(model, fake_branch)
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
    push!(new_branches, JuMP.@constraint(model, Lii <= wii <= Uii))
    push!(new_branches, JuMP.@constraint(model, Ljj <= wjj <= Ujj))
    push!(new_branches, JuMP.@constraint(model, Lij * wr <= wi))
    push!(new_branches, JuMP.@constraint(model, wi <= Uij * wr))
    push!(new_branches, JuMP.@constraint(model, πs[1] + πs[2] * wii + πs[3] * wjj + πs[4] * wr + πs[5] * wi >= Ujj * wii + Uii * wjj - Uii * Ujj))
    push!(new_branches, JuMP.@constraint(model, πs[1] + πs[2] * wii + πs[3] * wjj + πs[4] * wr + πs[5] * wi >= Ljj * wii + Lii * wjj - Lii * Ljj))
    return new_branches
end

# return a deepcopy of SpatialBCBranch for everything except for model
# this is needed to ensure variable references are not deepcopied
# the mod_branch will not be duplicated here if it is not nothing
function branch_copy(branch::SpatialBCBranch)
    return SpatialBCBranch(branch.i, branch.j, branch.wii, branch.wjj, branch.wr, branch.wi, 
                           branch.mod_branch, deepcopy(branch.bounds), deepcopy(branch.valid_ineq_coeffs))
end

function BB.branch!(tree::BB.AbstractTree, node::BB.AbstractNode)
    @info " Node id $(node.id), status $(node.solution_status), bound $(node.bound)"
    root = find_root(node)
    model = root.model
    if node.bound >= tree.best_incumbent
        if isapprox(node.bound, tree.best_incumbent, atol = 3)
            root.auxiliary_data["best_bound"] = node.bound
            root.auxiliary_data["best_id"] = node.id
        end
        delete_prev_branch_constr!(model, node)
        # @info " Fathomed by bound"
    elseif node.depth >= 10
        delete_prev_branch_constr!(model, node)
        # @info " Fathomed by maximum depth"
    elseif node.solution_status in [MOI.OPTIMAL, MOI.LOCALLY_SOLVED, MOI.SLOW_PROGRESS]
        if node.bound > root.auxiliary_data["best_bound"]
            root.auxiliary_data["best_bound"] = node.bound
            root.auxiliary_data["best_id"] = node.id
        end
        pm = root.auxiliary_data["PM"]
        (i,j) = find_min_eigen(node)
        new_sbc_branch = create_sbc_branch(i, j, node)
        (new_i,new_j) = find_branching_entry(new_sbc_branch, node) # This function implements MVSB/MVWB/RBEB
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
    else
        delete_prev_branch_constr!(model, node)
        # @info " Fathomed by solution status: $(node.solution_status)"
    end
end

# implement depth first rule
function BB.next_node(tree::BB.AbstractTree)
    # best bound
    sort!(tree::BB.AbstractTree) = Base.sort!(tree.nodes, by=x->x.depth)
    BB.sort!(tree)
    node = Base.pop!(tree.nodes)
    return node
end

function BB.termination(tree::BB.AbstractTree)
    @info "Tree nodes: processed $(length(tree.processed)), left $(length(tree.nodes)), total $(tree.node_counter), best bound $(tree.best_bound), best incumbent $(tree.best_incumbent)"
    if BB.isempty(tree)
        # @info "Completed the tree search"
        return true
    end
    if length(tree.processed) >= 1000
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

function BB.heuristics!(node::BB.JuMPNode)
    #=
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
    =#
end

###########################################################################
#                     Helper Functions for Debugging                      #
###########################################################################

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