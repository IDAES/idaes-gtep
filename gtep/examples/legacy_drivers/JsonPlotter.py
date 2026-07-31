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

from gtep.gtep_model import ExpansionPlanningModel
from gtep.gtep_data import ExpansionPlanningData
from gtep.gtep_solution import ExpansionPlanningSolution
from pyomo.core import TransformationFactory
from pyomo.environ import SolverFactory
from pyomo.contrib.appsi.solvers.highs import Highs
from pyomo.contrib.appsi.solvers.gurobi import Gurobi
from pyomo.contrib.appsi.solvers.ipopt import Ipopt

data_path = "./gtep/data/9bus"
data_object = ExpansionPlanningData()
data_object.load_prescient(data_path)

for key in data_object.md.data["elements"]["branch"]:
    data_object.md.data["elements"]["branch"][key]["capital_cost"] = 1000
    print(data_object.md.data["elements"]["branch"][key])

# mod_object = ExpansionPlanningModel(
#     stages=4,
#     data=data_object.md,
#     num_reps=1,
#     len_reps=1,
#     num_commit=24, # 24
#     num_dispatch=6, # 4
# )
# mod_object.create_model()


sol_object = ExpansionPlanningSolution()
# sol_object.load_from_model(mod_object)
sol_object.import_data_object(data_object)

load_numerical_results = True
if load_numerical_results:
    #     # sol_object.read_json("./gtep_solution.json")
    sol_object.read_json("./gtep_solution_DC_9busSanityTest.json")
plot_results = True
if plot_results:
    sol_object.plot_levels(save_dir="./gtep/DC9bus_sanityplots/")
