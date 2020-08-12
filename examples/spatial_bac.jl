using BranchAndBound, JuMP, PowerModels, DeNet
const BB = BranchAndBound
const PM = PowerModels

mutable struct NodeWRPowerModel <: AbstractWRModel @pm_fields end

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
end

# read data, and formulate root model using PowerModels
file = "../data/case9.m"
data = prase_file(file)
pm = instantiate_model(data, NodeWRPowerModel, build_opf)
root_model = pm.model

# initialize branch-and-cut tree
tree = BB.AbstractTree()
root_node = BB.JuMPNode(root_model)

# I can extend heuristics! here...
