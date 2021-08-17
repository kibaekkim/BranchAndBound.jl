# Main file for the spatial branch-and-cut

include("./spatial_bac.jl")
# include("./spatial_bac_sparse.jl")

function main(file::String)::Nothing
    data = parse_file(file)
    pm = instantiate_model(data, NodeWRMPowerModel, build_opf)
    optimizer = optimizer_with_attributes(Mosek.Optimizer, "QUIET" => true)
    # optimizer = optimizer_with_attributes(SCS.Optimizer, "max_iters" => 10000, "verbose" => 0)
    set_optimizer(pm.model, optimizer)

    tree, node = initialize(pm, file)

    @time BB.run(tree)

    if isempty(node.auxiliary_data["rank1"])
        println("No rank-1 solution found")
    else
        (bound, id) = findmin(node.auxiliary_data["rank1"])
        println("Best bound obtained at $(id), bound value $(bound)")
    end
end

# file = "/home/weiqizhang/anl/pglib-opf/pglib_opf_case3_lmbd.m"
# file = "/home/weiqizhang/anl/pglib-opf/api/pglib_opf_case3_lmbd__api.m"
# file = "/home/weiqizhang/anl/pglib-opf/sad/pglib_opf_case3_lmbd__sad.m"
# file = "/home/weiqizhang/anl/pglib-opf/pglib_opf_case5_pjm.m"
# file = "/home/weiqizhang/anl/pglib-opf/api/pglib_opf_case5_pjm__api.m"
# file = "/home/weiqizhang/anl/pglib-opf/api/pglib_opf_case14_ieee__api.m"
# file = "/home/weiqizhang/anl/pglib-opf/sad/pglib_opf_case14_ieee__sad.m" # this case is slow, and sees oscillating best bounds
# file = "/home/weiqizhang/anl/pglib-opf/api/pglib_opf_case24_ieee_rts__api.m"
# file = "/home/weiqizhang/anl/pglib-opf/sad/pglib_opf_case24_ieee_rts__sad.m"

# main(file)