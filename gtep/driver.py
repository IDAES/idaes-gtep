from pyomo.environ import ConcreteModel, Var, SolverFactory, value
from pyomo.environ import units as u
from gtep.gtep_model import ExpansionPlanningModel
from gtep.gtep_data import ExpansionPlanningData
from egret.data.model_data import ModelData
from pyomo.core import TransformationFactory


data_path = "./gtep/data/5bus"
data_object = ExpansionPlanningData()
data_object.load_prescient(data_path)
mod_object = ExpansionPlanningModel(
    data=data_object.md, num_reps=1, len_reps=1, num_commit=1, num_dispatch=1
)
mod_object.create_model()
TransformationFactory("gdp.bound_pretransformation").apply_to(mod_object.model)
TransformationFactory("gdp.bigm").apply_to(mod_object.model)
opt = SolverFactory("gurobi")
mod_object.results = opt.solve(mod_object.model, tee=True)
