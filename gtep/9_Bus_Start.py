from gtep.gtep_model import ExpansionPlanningModel
from gtep.gtep_data import ExpansionPlanningData
from gtep.gtep_solution import ExpansionPlanningSolution
from pyomo.core import TransformationFactory
from pyomo.contrib.appsi.solvers.gurobi import Gurobi
from IPython import embed
import pyomo.environ as pyo
from forestlib.ph import stochastic_program
from forestlib.ph import ProgressiveHedgingSolver

data_path = "./data/9bus"
data_object = ExpansionPlanningData()
data_object.load_prescient(data_path)

mod_object = ExpansionPlanningModel(
    stages=2,
    data=data_object.md,
    num_reps=2,
    len_reps=1,
    num_commit=2,
    num_dispatch=2,
)
mod_object.create_model()
embed()
latitude = 34.0522  # Example latitude
longitude = -118.2437  # Example longitude

# Adding latitude and longitude to bus1
for bus in mod_object.data.data["elements"]["bus"].keys():
    mod_object.data.data["elements"]["bus"][bus]["latitude"] = latitude
    mod_object.data.data["elements"]["bus"][bus]["longitude"] = longitude


TransformationFactory("gdp.bound_pretransformation").apply_to(mod_object.model) 
TransformationFactory("gdp.bigm").apply_to(mod_object.model) 
#embed()
opt = Gurobi()
#opt = Highs()

mod_object.results = opt.solve(mod_object.model) 
#embed()
#sol_object = ExpansionPlanningSolution()
#sol_object.load_from_model(mod_object)
#sol_object.dump_json("./gtep_solution.json")

#sol_object.import_data_object(data_object)

#sol_object.plot_levels(save_dir="./plots/")

