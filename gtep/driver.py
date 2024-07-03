from gtep.gtep_model import ExpansionPlanningModel
from gtep.gtep_data import ExpansionPlanningData
from gtep.gtep_solution import ExpansionPlanningSolution
from pyomo.core import TransformationFactory
from pyomo.contrib.appsi.solvers.highs import Highs
from pyomo.contrib.appsi.solvers.gurobi import Gurobi
from pyomo.environ import SolverFactory


data_path = "./gtep/data/WECC_Reduced_USAEE"
data_object = ExpansionPlanningData()
data_object.load_prescient(data_path)
mod_object = ExpansionPlanningModel(
    stages=2,
    data=data_object.md,
    num_reps=2,
    len_reps=1,
    num_commit=18,
    num_dispatch=4,
)
mod_object.create_model()
print("model build")
TransformationFactory("gdp.bound_pretransformation").apply_to(mod_object.model)
TransformationFactory("gdp.bigm").apply_to(mod_object.model)
print("model transformed")
opt = SolverFactory("gurobi")
mod_object.results = opt.solve(mod_object.model, tee=True)
# opt = Gurobi()
# opt = Highs()
# opt.config(visibility=True)
# with open("solver_log.txt", "w") as f:
#     opt.config.display(ostream=f)
#     mod_object.results = opt.solve(mod_object.model)


save_numerical_results = True
if save_numerical_results:

    sol_object = ExpansionPlanningSolution()

    sol_object.load_from_model(mod_object)
    sol_object.dump_json()
load_numerical_results = True
if load_numerical_results:
    # sol_object.read_json("./gtep_solution.json")
    sol_object.read_json("./oh_lawd_solution_wecc_reduced.json")
plot_results = True
if plot_results:
    sol_object.plot_levels(save_dir="./plots/")
