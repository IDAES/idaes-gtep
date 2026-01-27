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

# Add data
data_path = "./data/5bus"
data_object = ExpansionPlanningData()
data_object.load_prescient(data_path)

# [ESR WIP: Add bus and cost data files to be used on the
# DataProcessing class. This class processes data for the following
# type of generators: Natural Gas, Solar, Wind, and Coal. Note that
# during this data processing stage, the generator type is matched to
# existent generators in the data. The data contains the following
# types: (a) Natural Gas: Combustion Turbine (CT) and Fuel Efficiency
# (FE) and (b) Solar: Utility PV and Concentrated Solar Power (CSP)

bus_data_path = "data/costs/Bus_data_gen_weights_mappings.csv"
cost_data_path = "data/costs/2022_v3_Annual_Technology_Baseline_Workbook_Mid-year_update_2-15-2023_Clean.xlsx"
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
    stages=2,
    data=data_object,
    cost_data=data_processing_object,
    num_reps=2,
    len_reps=1,
    num_commit=6,
    num_dispatch=4,
    # [ESR: in min by default, for now]
    duration_dispatch=15,
)

mod_object.config["include_investment"] = True
mod_object.config["include_commitment"] = True
mod_object.config["include_redispatch"] = True
mod_object.config["scale_loads"] = True
mod_object.config["transmission"] = True
mod_object.config["storage"] = False
mod_object.config["flow_model"] = "DC"

mod_object.create_model()

# Apply transformations to logical terms
# pyo.TransformationFactory("gdp.bound_pretransformation").apply_to(mod_object.model)
pyo.TransformationFactory("gdp.bigm").apply_to(mod_object.model)

# Add solver
opt = pyo.SolverFactory("gurobi")
# opt = Gurobi()
# opt = Highs()
# opt = pyo.SolverFactory("highs")
mod_object.results = opt.solve(mod_object.model, tee=True)
# print(mod_object.results)
# mod_object.model.investmentStage.display()
# mod_object.report_model()

quit()

sol_object = ExpansionPlanningSolution()
sol_object.load_from_model(mod_object)
sol_object.dump_json("./gtep_solution.json")
sol_object.import_data_object(data_object)

# sol_object.read_json("./gtep_lots_of_buses_solution.json")  # "./gtep/data/WECC_USAEE"
# sol_object.read_json("./gtep_11bus_solution.json")  # "./gtep/data/WECC_Reduced_USAEE"
# sol_object.read_json("./gtep_solution.json")
# sol_object.read_json("./updated_gtep_solution_test.json")
# sol_object.read_json("./gtep_wiggles.json")
sol_object.plot_levels(save_dir="./plots/")

# save_numerical_results = False
# if save_numerical_results:

#     sol_object = ExpansionPlanningSolution()

#     sol_object.load_from_model(mod_object)
#     sol_object.dump_json()
# load_numerical_results = False

# if load_numerical_results:
#     # sol_object.read_json("./gtep_solution.json")
#     sol_object.read_json("./bigger_longer_wigglier_gtep_solution.json")
# plot_results = False

# if plot_results:
#     sol_object.plot_levels(save_dir="./plots/")


pass
