from gtep.gtep_model import ExpansionPlanningModel
from gtep.gtep_data import ExpansionPlanningData
from gtep.gtep_solution import ExpansionPlanningSolution
from pyomo.core import TransformationFactory
from pyomo.environ import SolverFactory
from pyomo.contrib.appsi.solvers.highs import Highs
from pyomo.contrib.appsi.solvers.gurobi import Gurobi
from pyomo.contrib.solver.solvers.gurobi_direct import GurobiDirect
import gurobipy as gp
import gc

gc.disable()

data_path = "./gtep/data/Texas_2000"
data_object = ExpansionPlanningData()
data_object.load_prescient(data_path)


# for data in data_object.representative_data:
#     print(data.data["elements"]['load'])

load_scaling_path = data_path + "/ERCOT-Adjusted-Forecast.xlsb"
data_object.import_load_scaling(load_scaling_path)
outage_path = data_path + "/feb_10.csv"
data_object.import_outage_data(outage_path)

data_object.texas_case_study_updates(data_path)


# for data in data_object.representative_data:
#     for gen in data.data["elements"]["generator"].keys():
#         print(data.data["elements"]['generator'][gen]["lifetime"])
# import sys
# sys.exit()
# Initial goal:
# 3 investment periods (now, 5 years, 10 years)
# 4 representative days (1 per season)
# Hourly commitment and dispatch

## RMA:
## Change num_reps from 4 to 5 to include extreme days

mod_object = ExpansionPlanningModel(
    stages=3, data=data_object, num_reps=5, len_reps=24, num_commit=24, num_dispatch=1
)
# print(mod_object.data.data["elements"]["generator"]["1"])
# import sys
# sys.exit()
mod_object.config["include_investment"] = True
mod_object.config["scale_loads"] = False
mod_object.config["scale_texas_loads"] = True
mod_object.config["transmission"] = False

mod_object.config["flow_model"] = "DC"

# mod_object.config["thermal_investment"] = True
# mod_object.config["renewable_investment"] = True
mod_object.create_model()
mod_object.timer.toc("horrible")

# import sys
# sys.exit()
# TransformationFactory("gdp.bound_pretransformation").apply_to(mod_object.model)
mod_object.timer.toc("double horrible")
# import sys
# sys.exit()
TransformationFactory("gdp.bigm").apply_to(mod_object.model)
mod_object.timer.toc("triple horrible")

# import sys
# sys.exit()

opt = SolverFactory("gurobi_direct_v2")
# opt = Gurobi()
# opt = GurobiDirect()
mod_object.timer.toc("Actually, I think this is garbage collection")
# opt.gurobi_options['LogFile'] = "basic_logging.log"
# opt.gurobi_options['LogToConsole'] = 1
# opt = Highs()
mod_object.timer.toc(
    "let's start to solve -- this is really the start of the handoff to gurobi"
)
mod_object.results = opt.solve(
    mod_object.model,
    tee=True,
    solver_options={"LogFile": "t2k_logging.log", "MIPGap": 0.001},
)


# mod_object.model.write('bad_sol.sol')
# mod_object.results = opt.solve(mod_object.model)

# import sys
# sys.exit()

mod_object.timer.toc("we've solved, let's pull investment variables")
import pyomo.environ as pyo
import pyomo.gdp as gdp

valid_names = ["Inst", "Oper", "Disa", "Ext", "Ret"]
# thermal_names = ["genInst", "genOper", "genDisa", "genExt", "genRet"]
renewable_investments = {}
dispatchable_investments = {}
load_shed = {}
for var in mod_object.model.component_objects(pyo.Var, descend_into=True):
    for index in var:
        if "Shed" in var.name:
            if pyo.value(var[index]) >= 0.001:
                load_shed[var.name + "." + str(index)] = pyo.value(var[index])
        for name in valid_names:
            if name in var.name:
                # print(var, index, pyo.value(var[index]))
                if pyo.value(var[index]) >= 0.001:
                    renewable_investments[var.name + "." + str(index)] = pyo.value(
                        var[index]
                    )
for var in mod_object.model.component_objects(gdp.Disjunct, descend_into=True):
    for index in var:
        for name in valid_names:
            if name in var.name:
                # print(var.name)
                # print(var, index, pyo.value(var[index].indicator_var))
                if pyo.value(var[index].indicator_var) == True:
                    dispatchable_investments[var.name + "." + str(index)] = pyo.value(
                        var[index].indicator_var
                    )

## RMA:
## You may want to save a few more things from the Expressions

costs = {}
for exp in mod_object.model.component_objects(pyo.Expression, descend_into=True):
    if "Cost" in exp.name or "cost" in exp.name:
        if type(exp) is ScalarExpression:
            costs[exp.name] = pyo.value(exp)
        if type(exp) is IndexedExpression:
            for e in exp:
                costs[exp[e].name] = pyo.value(exp[e])

import json
import os

## RMA:
## You can change where results are saved down here

folder_name = "informs_optimization_results"
renewable_investment_name = folder_name + "/renewable_investments.json"
dispatchable_investment_name = folder_name + "/dispatchable_investments.json"
load_shed_name = folder_name + "/load_shed.json"
costs_name = folder_name + "/costs.json"

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

mod_object.timer.toc("we've dumped; get everybody and the stuff together")
