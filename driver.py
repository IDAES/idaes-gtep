from pyomo.environ import ConcreteModel, Var, SolverFactory, value
from pyomo.environ import units as u
from gtep.gtep_model import ExpansionPlanningModel
from gtep.gtep_data import ExpansionPlanningData
from gtep.gtep_solution import ExpansionPlanningSolution
from egret.data.model_data import ModelData
from pyomo.core import TransformationFactory
from pyomo.contrib.appsi.solvers.highs import Highs


data_path = "./gtep/data/5bus"
data_object = ExpansionPlanningData()
data_object.load_prescient(data_path)
mod_object = ExpansionPlanningModel(
    data=data_object.md, num_reps=2, len_reps=3, num_commit=4, num_dispatch=5
)
mod_object.create_model()
TransformationFactory("gdp.bound_pretransformation").apply_to(mod_object.model)
TransformationFactory("gdp.bigm").apply_to(mod_object.model)
# opt = SolverFactory("gurobi")
opt = Highs()
# mod_object.results = opt.solve(mod_object.model, tee=True)
mod_object.results = opt.solve(mod_object.model)

sol_object = ExpansionPlanningSolution()

sol_object.load_from_model(mod_object)
sol_object.dump_json()

sol_object.read_json("./gtep_solution.json")
sol_object.plot_dispatch_level(save_dir="./plots/")
pass
