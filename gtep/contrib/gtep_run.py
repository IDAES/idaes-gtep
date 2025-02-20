from gtep.gtep_model import ExpansionPlanningModel
from gtep.gtep_data import ExpansionPlanningData
from gtep.gtep_solution import ExpansionPlanningSolution
from pyomo.core import TransformationFactory
from pyomo.contrib.appsi.solvers.highs import Highs
from pyomo.contrib.appsi.solvers.gurobi import Gurobi
import gurobipy as gp
from pyomo.environ import SolverFactory, Var, Expression
import csv, os


data_path = "./gtep/data/5bus"
data_object = ExpansionPlanningData()
data_object.load_prescient(data_path)

mod_object = ExpansionPlanningModel(
    stages=2,
    data=data_object.md,
    num_reps=2,
    len_reps=1,
    num_commit=6,
    num_dispatch=4,
)
mod_object.create_model()
TransformationFactory("gdp.bound_pretransformation").apply_to(mod_object.model)
TransformationFactory("gdp.bigm").apply_to(mod_object.model)
# opt = SolverFactory("gurobi")
# opt = Gurobi()
opt = Highs()
# mod_object.results = opt.solve(mod_object.model, tee=True)
mod_object.results = opt.solve(mod_object.model)

# Specify the output directory and filename
# output_dir = "output"  # Change to your desired directory
# os.makedirs(output_dir, exist_ok=True)  # Create directory if it doesn't exist
# output_file = os.path.join(output_dir, "optimal_variable_values_wo_reliability.csv")

# results = []
# for var in mod_object.model.component_objects(Var, active=True):
#     var_name = var.name
#     for index in var:
#         value = var[index].value
#         results.append((f"{var_name}[{index}]", value))

# for expr in mod_object.model.component_objects(Expression, active=True):
#     expr_name = expr.name
#     for index in expr:
#         try:
#             value = expr[index]()
#             results.append((f"{expr_name}[{index}]", value))
#         except ValueError:
#             results.append((f"{expr_name}[{index}]", None))

# with open(output_file, "w", newline="") as file:
#     writer = csv.writer(file)
#     writer.writerow(["Name", "Value"])  # Header row
#     for row in results:
#         writer.writerow(row)

sol_object = ExpansionPlanningSolution()
sol_object.load_from_model(mod_object)
sol_object.dump_json("./gtep_solution.json")
