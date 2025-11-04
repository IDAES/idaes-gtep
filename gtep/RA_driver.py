from gtep.gtep_model import ExpansionPlanningModel
from gtep.gtep_data import ExpansionPlanningData
from gtep.gtep_solution import ExpansionPlanningSolution
from pyomo.core import TransformationFactory
from pyomo.contrib.appsi.solvers.gurobi import Gurobi
from IPython import embed
import pyomo.environ as pyo
import time
import matplotlib.pyplot as plt
import numpy as np

start_time = time.time()

data_path = "./data/5bus"
data_object = ExpansionPlanningData()
data_object.load_prescient(data_path)
# data_object.load_storage_csv(data_path)

mod_object = ExpansionPlanningModel(
    stages=2,
    data=data_object,
    num_reps=2,  # num rep days
    len_reps=24,  # len rep days
    num_commit=24,  # num commitment periods
    num_dispatch=2,  # num dispatch per commitment period
)

mod_object.config["include_commitment"] = False

mod_object.config["flow_model"] = "CP"

mod_object.config["transmission"] = True  # TRANSMISSION INVESTMENT FLAG
mod_object.config["thermal_generation"] = True  # THERMAL GENERATION INVESTMENT FLAG
mod_object.config["renewable_generation"] = True  # RENEWABLE GENERATION INVESTMENT FLAG
mod_object.config["scale_loads"] = False  # LEAVE AS FALSE
mod_object.config["scale_texas_loads"] = False  # LEAVE AS FALSE

mod_object.create_model()
TransformationFactory("gdp.bound_pretransformation").apply_to(mod_object.model)
TransformationFactory("gdp.bigm").apply_to(mod_object.model)
opt = Gurobi()
# opt = Highs()

mod_object.results = opt.solve(mod_object.model)
# sol_object = ExpansionPlanningSolution()
# sol_object.load_from_model(mod_object)

end_time = time.time()
elapsed_time = end_time - start_time
print(f"Elapsed time: {elapsed_time} seconds")
# sol_object.dump_json("./gtep_solution.json")

# sol_object.import_data_object(data_object)

# sol_object.plot_levels(save_dir="./plots/")

import pyomo.gdp as gdp

valid_names = ["Inst", "Oper", "Disa", "Ext", "Ret"]
thermal_names = ["genInst", "genOper", "genDisa", "genExt", "genRet"]
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

costs = {}
for exp in mod_object.model.component_objects(pyo.Expression, descend_into=True):
    if "operatingCost" in exp.name:
        costs[var.name] = pyo.value(exp)
    elif "investmentCost" in exp.name:
        costs[var.name] = pyo.value(exp)


import json

import os

if not os.path.exists("retirement_allowed_no_extreme_full_load"):
    os.makedirs("retirement_allowed_no_extreme_full_load")
folder_name = "retirement_allowed_no_extreme_full_load"
costs_name = folder_name + "costs.json"
with open(
    "retirement_allowed_no_extreme_full_load/renewable_investments.json", "w"
) as fil:
    json.dump(renewable_investments, fil)
with open(
    "retirement_allowed_no_extreme_full_load/dispatchable_investments.json", "w"
) as fil:
    json.dump(dispatchable_investments, fil)
with open("retirement_allowed_no_extreme_full_load/load_shed.json", "w") as fil:
    json.dump(load_shed, fil)
with open(costs_name, "w") as fil:
    json.dump(costs, fil)
