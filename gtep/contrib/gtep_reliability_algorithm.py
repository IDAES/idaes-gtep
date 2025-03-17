from gtep.gtep_model import ExpansionPlanningModel
from gtep.gtep_solution import ExpansionPlanningSolution
from gtep.gtep_data import ExpansionPlanningData
from gtep.contrib.gtep_reliability_model import ExpansionPlanningModelwithReliability
from gtep.contrib.gtep_model_result import solve_expansion_model
from pyomo.core import TransformationFactory
from pyomo.contrib.appsi.solvers.highs import Highs
from pyomo.contrib.appsi.solvers.gurobi import Gurobi
from pyomo.environ import Var, Expression, SolverFactory
import csv
import more_itertools


# Call dataset
data_path = "./gtep/data/SanDiego"
data_object = ExpansionPlanningData()
data_object.load_prescient(data_path)

num_planning_year = 2
num_rep_day = 2
length_rep_day = 1
num_commit_hour = 6
num_dispat_min = 4

# Call the expansion planning model without reliability
mod_object = ExpansionPlanningModel(
    stages=num_planning_year,
    data=data_object.md,
    num_reps=num_rep_day,
    len_reps=length_rep_day,
    num_commit=num_commit_hour,
    num_dispatch=num_dispat_min,
)


# Solve expansion planning model without reliability
# Export the results for the reliability-constrained model
expansion_model_results = solve_expansion_model(mod_object)


sol_object = ExpansionPlanningSolution()
sol_object.load_from_model(mod_object)
sol_object.dump_json("./gtep_pre_reliability_solution.json")


# Export the results of expansion planning without reliability
results_ref = []
for var in mod_object.model.component_objects(Var):
    var_name = var.name
    for index in var:
        values = var[index].value
        results_ref.append((f"{var_name}[{index}]", values))

for expr in mod_object.model.component_objects(Expression):
    expr_name = expr.name
    for index in expr:
        try:
            values = expr[index]()
            results_ref.append((f"{expr_name}[{index}]", values))
        except ValueError:
            results_ref.append((f"{expr_name}[{index}]", None))

with open("optimal_variable_values_without_reliability.csv", "w", newline="") as file:
    writer = csv.writer(file)
    writer.writerow(["Name", "Value"])  # Header row
    for row in results_ref:
        writer.writerow(row)

# Call expansion planning model with reliability
mod_object_rel = ExpansionPlanningModelwithReliability(
    stages=num_planning_year,
    data=data_object.md,
    num_reps=num_rep_day,
    len_reps=length_rep_day,
    num_commit=num_commit_hour,
    num_dispatch=num_dispat_min,
    rel_data=expansion_model_results,
)
mod_object_rel.create_model()


# Calculate the probability of capacity failure state based on the probability of failure
prob_state = {
    (bus, state): 1
    for bus in mod_object_rel.model.criticalBuses
    for state in mod_object_rel.model.states
}
failure = {
    gen: mod_object_rel.model.failureRate[gen]
    for gen in mod_object_rel.model.generators
}

for bus in mod_object_rel.model.criticalBuses:
    crit_generator_sets = list(
        more_itertools.powerset(mod_object_rel.model.criticalGenerators[bus])
    )
    state_idx = 1
    for cg_set in crit_generator_sets:
        fail_rate = 1
        for gen in mod_object_rel.model.criticalGenerators[bus]:
            if gen in cg_set:
                fail_rate *= failure[gen]
            else:
                fail_rate *= 1 - failure[gen]
        prob_state[bus, state_idx] = fail_rate
        state_idx += 1


# Update probability of failure of capacity failure state
for bus in mod_object_rel.model.criticalBuses:
    for state in mod_object_rel.model.states:
        mod_object_rel.model.prob[bus, state] = prob_state[bus, state]


# Transform GDP to MIP
TransformationFactory("gdp.bound_pretransformation").apply_to(mod_object_rel.model)
TransformationFactory("gdp.bigm").apply_to(mod_object_rel.model)

# Solve expansion planning model with reliability
opt = Gurobi()
mod_object_rel.results = opt.solve(mod_object_rel.model)


# Export results
results_rel = []
for var in mod_object_rel.model.component_objects(Var):
    var_name = var.name
    for index in var:
        values = var[index].value
        results_rel.append((f"{var_name}[{index}]", values))

for expr in mod_object_rel.model.component_objects(Expression):
    expr_name = expr.name
    for index in expr:
        try:
            values = expr[index]()
            results_rel.append((f"{expr_name}[{index}]", values))
        except ValueError:
            results_rel.append((f"{expr_name}[{index}]", None))

with open("optimal_variable_values_with_reliability.csv", "w", newline="") as file:
    writer = csv.writer(file)
    writer.writerow(["Name", "Value"])  # Header row
    for row in results_rel:
        writer.writerow(row)

sol_object = ExpansionPlanningSolution()
sol_object.load_from_model(mod_object_rel)
sol_object.dump_json("./gtep_reliability_solution.json")
