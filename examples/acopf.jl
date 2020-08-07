using BranchAndBound
using PowerModels
using JuMP
using Ipopt

const BB = BranchAndBound
const PM = PowerModels

mutable struct NodeWRMModel <: PM.AbstractWRMModel PM.@pm_fields end

function PM.constraint_model_voltage(pm::NodeWRMModel, n::Int)

    WR = var(pm, n)[:WR]
    WI = var(pm, n)[:WI]

    nw = pm.cnw
    bus_ids = ids(pm, nw, :bus)
    lookup_w_index = Dict((bi,i) for (i,bi) in enumerate(bus_ids))
    
    # create voltage variables
    PM.variable_bus_voltage_real(pm, bounded=true)
    PM.variable_bus_voltage_imaginary(pm, bounded=true)

    vr = var(pm, n)[:vr]
    vi = var(pm, n)[:vi]

    # constraints on diagonal
    for (i, bus) in ref(pm, nw, :bus)
        w_idx = lookup_w_index[i]

        # vr
        JuMP.@constraint(pm.model, -bus["vmax"] * vr[i] - WR[i,i] <= 0)
        JuMP.@constraint(pm.model, bus["vmax"] * vr[i] - WR[i,i] >= 0)

        # vi
        JuMP.@constraint(pm.model, -bus["vmax"] * vi[i] - WI[i,i] <= 0)
        JuMP.@constraint(pm.model, bus["vmax"] * vi[i] - WI[i,i] >= 0)
    end

    # constraints on off-diagonal
    for (i,j) in ids(pm, nw, :buspairs)
        wi_idx = lookup_w_index[i]
        wj_idx = lookup_w_index[j]
        bus = ref(pm, nw, :bus)[i]

        # impose skew symmetric
        JuMP.@constraint(pm.model, WI[i,j] + WI[j,i] == 0)

        # vr
        JuMP.@constraint(pm.model, -bus["vmax"] * vr[j] - WR[i,j] <= 0)
        JuMP.@constraint(pm.model, bus["vmax"] * vr[j] - WR[i,j] >= 0)

        # vi
        JuMP.@constraint(pm.model, -bus["vmax"] * vi[j] - WI[i,j] <= 0)
        JuMP.@constraint(pm.model, bus["vmax"] * vi[j] - WI[i,j] >= 0)
    end
end

pm_data = PM.parse_file("../data/case9.m")
pm = PM.instantiate_model(pm_data, NodeWRMModel, PM.build_opf)

vr = var(pm, pm.cnw)[:vr]
vi = var(pm, pm.cnw)[:vi]
WR = var(pm, pm.cnw)[:WR]
WI = var(pm, pm.cnw)[:WI]

pm_data = PM.parse_file("../data/case9.m")
pm = PM.instantiate_model(pm_data, NodeWRMModel, PM.build_opf)

function BB.branch(node::BB.Node)
    vr_val = BB.node_solution.(node, vr)
    vi_val = BB.node_solution.(node, vi)
    WR_val = BB.node_solution.(node, WR)
    WI_val = BB.node_solution.(node, WI)

    for (i, bus) in ref(pm, pm.cnw, :bus)
        @show abs(vr_val[i]^2 - WR_val[i,i])
        @show abs(vi_val[i]^2 - WI_val[i,i])
    end
end
