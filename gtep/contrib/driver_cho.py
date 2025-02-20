from gtep.gtep_model_cho import ExpansionPlanningModelwithReliability
from gtep.gtep_data_cho import ExpansionPlanningDataforReliability
from gtep.gtep_solution_cho import ExpansionPlanningSolutionwithReliability
from pyomo.core import TransformationFactory
from pyomo.contrib.appsi.solvers.highs import Highs
from pyomo.contrib.appsi.solvers.gurobi import Gurobi
import gurobipy as gp
import pyomo.environ as pyo
from pyomo.environ import SolverFactory, Var, Expression
import csv, os


data_path = "./gtep/data/5bus"
data_object = ExpansionPlanningDataforReliability()
data_object.load_prescient(data_path)

mod_object = ExpansionPlanningModelwithReliability(
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
opt = SolverFactory("gurobi")
# opt = Gurobi()
# opt = Highs()
# mod_object.results = opt.solve(mod_object.model, tee=True)
mod_object.results = opt.solve(mod_object.model, tee=True)

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

# with open("optimal_variable_values_with_reliability.csv", "w", newline="") as file:
#     writer = csv.writer(file)
#     writer.writerow(["Name", "Value"])  # Header row
#     for row in results:
#         writer.writerow(row)

# sol_object = ExpansionPlanningSolutionwithReliability()
# sol_object.load_from_model(mod_object)
# sol_object.dump_json("./gtep_solution_reliability.json")


for stage in mod_object.model.investmentStage:
    for bus in mod_object.model.criticalBuses:
        for state in mod_object.model.states:
            print(
                "bus",
                bus,
                "state",
                state,
                "stage",
                stage,
                "production",
                mod_object.model.investmentStage[stage].prod_state[bus, state].value,
            )
