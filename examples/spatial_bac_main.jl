include("./spatial_bac.jl")

# read data, and formulate root model using PowerModels
file = "/home/weiqizhang/.julia/dev/BranchAndBound/data/case5.m"
# file = "/home/weiqizhang/.julia/dev/BranchAndBound/data/case9.m"
# file = "/home/weiqizhang/anl/pglib-opf/sad/pglib_opf_case3_lmbd__sad.m"
# file = "/home/weiqizhang/anl/pglib-opf/sad/pglib_opf_case5_pjm__sad.m"
# file = "/home/weiqizhang/anl/pglib-opf/sad/pglib_opf_case14_ieee__sad.m" # this case is slow, and sees oscillating best bounds
# file = "/home/weiqizhang/anl/pglib-opf/api/pglib_opf_case14_ieee__api.m"
# file = "/home/weiqizhang/anl/pglib-opf/sad/pglib_opf_case24_ieee_rts__sad.m"
# file = "/home/weiqizhang/anl/pglib-opf/api/pglib_opf_case30_as__api.m"
data = parse_file(file)
pm = instantiate_model(data, NodeWRMPowerModel, build_opf)
optimizer = optimizer_with_attributes(Mosek.Optimizer, "QUIET" => true)
ipopt_optimizer = optimizer_with_attributes(Ipopt.Optimizer, "print_level" => 0)
set_optimizer(pm.model, optimizer)

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
tree = BB.initialize_tree(node)

# set incumbent
ipopt_solution = run_opf(file, ACRPowerModel, ipopt_optimizer)
tree.best_incumbent = ipopt_solution["objective"]

@time BB.run(tree)
