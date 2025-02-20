# -*- coding: utf-8 -*-
"""
Created on Tue Jun  4 14:56:46 2024

@author: jscelay

Test driver for IDAES-GTEP
"""

from pyomo.environ import ConcreteModel, Var, SolverFactory, value
from pyomo.environ import units as u
from gtep.gtep_model import ExpansionPlanningModel
from gtep.gtep_data import ExpansionPlanningData
from gtep.gtep_solution import ExpansionPlanningSolution
from egret.data.model_data import ModelData
from pyomo.core import TransformationFactory
from pyomo.contrib.appsi.solvers.highs import Highs
from pyomo.contrib.appsi.solvers.gurobi import Gurobi

# data_path = "./data/5bus"
data_path = "./gtep/data/5bus_jsc"
data_object = ExpansionPlanningData()
data_object.load_prescient(data_path)
# mod_object = ExpansionPlanningModel(
#    data=data_object.md, num_reps=2, len_reps=1, stages=2, num_commit=24, num_dispatch=12
# )
mod_object = ExpansionPlanningModel(
    data=data_object.md, num_reps=1, len_reps=1, stages=1, num_commit=24, num_dispatch=2
)
mod_object.create_model()
TransformationFactory("gdp.bound_pretransformation").apply_to(mod_object.model)
TransformationFactory("gdp.bigm").apply_to(mod_object.model)
opt = Highs()
# opt = SolverFactory("gurobi", solver_io="python")
# opt = Gurobi()
# mod_object.results = opt.solve(mod_object.model, tee=True)
mod_object.results = opt.solve(mod_object.model)


save_numerical_results = True
if save_numerical_results:

    sol_object = ExpansionPlanningSolution()

    sol_object.load_from_model(mod_object)
    sol_object.dump_json()
load_numerical_results = True
if load_numerical_results:
    sol_object.read_json("./gtep_solution_jscTest.json")
    # sol_object.read_json("./gtep_solution.json")
    # sol_object.read_json("./bigger_longer_wigglier_gtep_solution.json")
plot_results = True
if plot_results:
    sol_object.plot_levels(save_dir="./plots/")

pass
