# Main file for the spatial branch-and-cut

include("./spatial_bac.jl")


# read data, and formulate root model using PowerModels
# file = "/home/weiqizhang/anl/pglib-opf/pglib_opf_case3_lmbd.m"
# file = "/home/weiqizhang/anl/pglib-opf/api/pglib_opf_case3_lmbd__api.m"
# file = "/home/weiqizhang/anl/pglib-opf/sad/pglib_opf_case3_lmbd__sad.m"
file = "/home/weiqizhang/anl/pglib-opf/pglib_opf_case5_pjm.m"
# file = "/home/weiqizhang/anl/pglib-opf/api/pglib_opf_case5_pjm__api.m"
# file = "/home/weiqizhang/anl/pglib-opf/api/pglib_opf_case14_ieee__api.m"
# file = "/home/weiqizhang/anl/pglib-opf/sad/pglib_opf_case14_ieee__sad.m" # this case is slow, and sees oscillating best bounds
# file = "/home/weiqizhang/anl/pglib-opf/api/pglib_opf_case24_ieee_rts__api.m"
# file = "/home/weiqizhang/anl/pglib-opf/sad/pglib_opf_case24_ieee_rts__sad.m"
data = parse_file(file)
pm = instantiate_model(data, NodeWRMPowerModel, build_opf)
optimizer = optimizer_with_attributes(Mosek.Optimizer, "QUIET" => true)
set_optimizer(pm.model, optimizer)

tree, node = initialize(pm)

@time BB.run(tree)

println("Best bound obtained at $(tree.processed[1].auxiliary_data["best_id"]), bound value $(tree.processed[1].auxiliary_data["best_bound"])")