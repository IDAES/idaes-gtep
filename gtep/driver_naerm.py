#################################################################################
# The Institute for the Design of Advanced Energy Systems Integrated Platform
# Framework (IDAES IP) was produced under the DOE Institute for the
# Design of Advanced Energy Systems (IDAES).
#
# Copyright (c) 2018-2026 by the software owners: The Regents of the
# University of California, through Lawrence Berkeley National Laboratory,
# National Technology & Engineering Solutions of Sandia, LLC, Carnegie Mellon
# University, West Virginia University Research Corporation, et al.
# All rights reserved.  Please see the files COPYRIGHT.md and LICENSE.md
# for full copyright and license information.
#################################################################################

import pyomo.environ as pyo
from pyomo.contrib.appsi.solvers.highs import Highs
from pyomo.contrib.appsi.solvers.gurobi import Gurobi

from gtep.gtep_model import ExpansionPlanningModel
from gtep.gtep_data import ExpansionPlanningData
from gtep.gtep_solution import ExpansionPlanningSolution
from gtep.gtep_data_processing import DataProcessing

import xpress


# Add data
data_path = "./gtep/data/WECC_ADS_PNNL"
# data_path = "./data/5bus"
data_object = ExpansionPlanningData(
    stages=2,
    num_reps=2,
    num_commit=24,
    num_dispatch=1,
    duration_representative_period=24,
    save_period_structure_file=False,
    period_structure_json_file=None,
    # period_structure_json_file="period_structure_from_gtep.json",
)
data_object.load_prescient(data_path)

# print(data_object.md.data["elements"]["generator"]["AESO_solar"]["p_max"])
# quit()
# c = 0
# # print(data_object.md.data["elements"]["generator"].keys())
# for gen in data_object.md.data["elements"]["generator"].keys():
#     if data_object.md.data["elements"]["generator"][gen]["generator_type"] == "renewable":
#         c += 1
#         print(gen)
#         print(data_object.md.data["elements"]["generator"][gen]["p_max"])
#         if c == 2:
#             quit()

# [ESR WIP: Add bus and cost data files to be used on the
# DataProcessing class. This class processes data for the following
# type of generators: Natural Gas, Solar, Wind, and Coal. Note that
# during this data processing stage, the generator type is matched to
# existent generators in the data. The data contains the following
# types: (a) Natural Gas: Combustion Turbine (CT) and Fuel Efficiency
# (FE) and (b) Solar: Utility PV and Concentrated Solar Power (CSP)

bus_data_path = "./gtep/data/costs/Bus_data_gen_weights_mappings.csv"
cost_data_path = "./gtep/data/costs/2022_v3_Annual_Technology_Baseline_Workbook_Mid-year_update_2-15-2023_Clean.xlsx"
candidate_gens = [
    "Natural Gas_CT",
    "Natural Gas_FE",
    "Solar - Utility PV",
    "Land-Based Wind",
]

data_processing_object = DataProcessing()
data_processing_object.load_gen_data(
    bus_data_path=bus_data_path,
    cost_data_path=cost_data_path,
    candidate_gens=candidate_gens,
    save_csv=True,
)

# Populate and create GTEP model
mod_object = ExpansionPlanningModel(
    data=data_object,
    cost_data=data_processing_object,
)

mod_object.config["include_investment"] = True
mod_object.config["include_commitment"] = True
mod_object.config["include_redispatch"] = True
mod_object.config["scale_loads"] = False
mod_object.config["transmission"] = True
mod_object.config["storage"] = False
mod_object.config["flow_model"] = "transport"

mod_object.create_model()
print("model is created!")
#mod_object.model.investmentStage[1].representativePeriod[1].commitmentPeriod[1].dispatchPeriod[1].branchInUse["AESO_BCHA"].dc_power_flow.pprint()

pyo.TransformationFactory("gdp.bigm").apply_to(mod_object.model)
print("model is transformed!")

#mod_object.model.investmentStage[1].representativePeriod[1].commitmentPeriod[1].dispatchPeriod[1].branchInUse["AESO_BCHA"].dc_power_flow.pprint()

solver = "xpress"
opt = pyo.SolverFactory(solver)
xpress.init("/Users/jkskolf/naerm_xpauth.xpr")
if solver == "xpress":
    log_folder = "xpress_log_files"
    options_dict = {
        "logfile": log_folder + "/" + solver + ".log",
    }

    mod_object.results = opt.solve(
        mod_object.model,
        tee=True,
        #logfile=log_folder + "/" + solver + ".log",
    )

import pyomo.environ as pyo
import pyomo.gdp as gdp

import json
import os
from pyomo.core.base.expression import ScalarExpression, IndexedExpression

folder_name = "NAERM_initial_testing"

valid_names = ["Inst", "Oper", "Disa", "Ext", "Ret"]
# thermal_names = ["genInst", "genOper", "genDisa", "genExt", "genRet"]
renewable_investments = {}
dispatchable_investments = {}
load_shed = {}
power_flow = {}
for var in mod_object.model.component_objects(pyo.Var, descend_into=True):
    for index in var:
        if "Shed" in var.name:
            if pyo.value(var[index]) >= 0.001:
                load_shed[var.name + "." + str(index)] = pyo.value(var[index])
        elif "Flow" in var.name:
            if pyo.value(var[index]) >= 0.001:
                power_flow[var.name + "." + str(index)] = pyo.value(var[index])
        for name in valid_names:
            if name in var.name:
                if pyo.value(var[index]) >= 0.001:
                    renewable_investments[var.name + "." + str(index)] = pyo.value(
                        var[index]
                    )
for var in mod_object.model.component_objects(gdp.Disjunct, descend_into=True):
    for index in var:
        for name in valid_names:
            if name in var.name:
                if pyo.value(var[index].indicator_var) == True:
                    dispatchable_investments[var.name + "." + str(index)] = pyo.value(
                        var[index].indicator_var
                    )

costs = {}
for exp in mod_object.model.component_objects(pyo.Expression, descend_into=True):
    if "Cost" in exp.name or "cost" in exp.name:
        if type(exp) is ScalarExpression:
            costs[exp.name] = pyo.value(exp)
        if type(exp) is IndexedExpression:
            for e in exp:
                costs[exp[e].name] = pyo.value(exp[e])

renewable_investment_name = folder_name + "/renewable_investments.json"
dispatchable_investment_name = folder_name + "/dispatchable_investments.json"
load_shed_name = folder_name + "/load_shed.json"
costs_name = folder_name + "/costs.json"
flow_name = folder_name + "/flows.json"

if not os.path.exists(folder_name):
    os.makedirs(folder_name)

with open(renewable_investment_name, "w") as fil:
    json.dump(renewable_investments, fil)
with open(dispatchable_investment_name, "w") as fil:
    json.dump(dispatchable_investments, fil)
with open(load_shed_name, "w") as fil:
    json.dump(load_shed, fil)
with open(costs_name, "w") as fil:
    json.dump(costs, fil)
with open(flow_name, "w") as fil:
    json.dump(power_flow, fil)