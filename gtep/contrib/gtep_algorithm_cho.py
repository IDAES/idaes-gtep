from gtep.gtep_model import ExpansionPlanningModel
from gtep.gtep_model_cho2 import ExpansionPlanningModelwithReliability
from gtep.gtep_data import ExpansionPlanningData
from gtep.gtep_data_cho import ExpansionPlanningDataforReliability
from pyomo.core import TransformationFactory
from pyomo.contrib.appsi.solvers.highs import Highs
from pyomo.contrib.appsi.solvers.gurobi import Gurobi
from collections import defaultdict
import gurobipy as gp
from pyomo.environ import *
from pyomo.environ import Set, Param, Var, Expression, Constraint, SolverFactory, value
import csv, os
import pandas as pd

# Call dataset
data_path = "./gtep/data/5bus"
data_object = ExpansionPlanningDataforReliability()
data_object.load_prescient(data_path)

num_planning_year = 2
num_rep_day = 2
length_rep_day = 1
num_commit_hour = 6
num_dispat_min = 4

# Call and solve expansion planning model without reliability
mod_object = ExpansionPlanningModel(
    stages=num_planning_year,
    data=data_object.md,
    num_reps=num_rep_day,
    len_reps=length_rep_day,
    num_commit=num_commit_hour,
    num_dispatch=num_dispat_min,
)
mod_object.create_model()
TransformationFactory("gdp.bound_pretransformation").apply_to(mod_object.model)
TransformationFactory("gdp.bigm").apply_to(mod_object.model)
mod_object_opt = SolverFactory("gurobi")
mod_object.results = mod_object_opt.solve(mod_object.model, tee=True)


# Call expansion planning model with reliability
mod_object_rel = ExpansionPlanningModelwithReliability(
    stages=num_planning_year,
    data=data_object.md,
    num_reps=num_rep_day,
    len_reps=length_rep_day,
    num_commit=num_commit_hour,
    num_dispatch=num_dispat_min,
)
mod_object_rel.create_model()

# Update critical buses and generators
# before running expansion planning model with reliability

# ----------------------- Step 1: Identify critical nodes -------------- #
# NOTE: The nodes with largest demand are selected, but it could be changed
loads = {
    bus: mod_object.model.loads[bus]
    for bus in mod_object.model.loads
    if mod_object.model.loads[bus] > 0
}
sorted_buses = sorted(loads.keys(), key=lambda bus: loads[bus], reverse=True)
critical_buses = sorted_buses[
    :2
]  # TODO: only select top 2 buses for now, but could be changed
noncritical_buses = list(set(mod_object_rel.model.buses) - set(critical_buses))
# print(critical_buses, noncritical_buses)

mod_object_rel.model.criticalBuses.update(critical_buses)
mod_object_rel.model.noncriticalBuses.update(noncritical_buses)
# print(list(mod_object_rel.model.criticalBuses))
# print(list(mod_object_rel.model.noncriticalBuses))


# ----------------------- Step 2: Identify critical generators --------- #
# Calculate the power output of each generators
total_thermal_generation = {gen: 0 for gen in mod_object.model.thermalGenerators}
for inv_stage in mod_object.model.investmentStage.values():
    for rep_period in inv_stage.representativePeriod.values():
        for commit_period in rep_period.commitmentPeriod.values():
            for dispatch in commit_period.dispatchPeriod.values():
                for gen in mod_object.model.thermalGenerators:
                    thermalgen_value = value(dispatch.thermalGeneration[gen])
                    total_thermal_generation[gen] += thermalgen_value


total_renewable_generation = {gen: 0 for gen in mod_object.model.renewableGenerators}
for inv_stage in mod_object.model.investmentStage.values():
    for rep_period in inv_stage.representativePeriod.values():
        for commit_period in rep_period.commitmentPeriod.values():
            for dispatch in commit_period.dispatchPeriod.values():
                for gen in mod_object.model.renewableGenerators:
                    renewablegen_value = value(dispatch.renewableGeneration[gen])
                    total_renewable_generation[gen] += renewablegen_value

total_generation = {}
for gen in mod_object.model.generators:
    if gen in mod_object.model.thermalGenerators:
        total_generation[gen] = total_thermal_generation[gen]
    elif gen in mod_object.model.renewableGenerators:
        total_generation[gen] = total_renewable_generation[gen]


# Calculate actual capacity factors of critical generators in critical nodes
average_capacity_factor = {}
for gen in mod_object.model.generators:
    if gen in mod_object.model.thermalGenerators:
        if total_thermal_generation[gen] != 0:
            average_capacity_factor[gen] = abs(
                total_thermal_generation[gen]
                / (
                    mod_object.model.thermalCapacity[gen]
                    * num_planning_year
                    * num_rep_day
                    * num_commit_hour
                    * num_dispat_min
                )
            )
        else:
            average_capacity_factor[gen] = 0

    elif gen in mod_object.model.renewableGenerators:
        if total_renewable_generation[gen] != 0:
            average_capacity_factor[gen] = abs(
                total_renewable_generation[gen]
                / (
                    mod_object.model.renewableCapacity[gen]
                    * num_planning_year
                    * num_rep_day
                    * num_commit_hour
                    * num_dispat_min
                )
            )
        else:
            average_capacity_factor[gen] = 0
# for gen in mod_object.model.thermalGenerators:
#     print(gen, "capacity", mod_object.model.thermalCapacity[gen])
# for gen in mod_object.model.renewableGenerators:
#     print(gen, "capacity", mod_object.model.renewableCapacity[gen])
# print(average_capacity_factor)


# Update average capacity factor
for gen in mod_object_rel.model.generators:
    mod_object_rel.model.averageCapacityFactor[gen] = average_capacity_factor[gen]


# Select top three generators from all generators as critical generators from critical buses
# Select non-critical generators from critical buses
# TODO: the number of critical generators selected should be flexible
filtered_generation = {}
for gen, output in total_generation.items():
    bus_id = gen.split("_")[0]
    bus_name = f"bus{bus_id}"

    if bus_name in critical_buses:
        filtered_generation[gen] = output

sorted_generators = sorted(
    filtered_generation.items(), key=lambda x: x[1], reverse=True
)
critical_generators = [
    gen[0] for gen in sorted_generators[:3]
]  # TODO: only select top 2 generators
non_critical_generators = [
    gen for gen in filtered_generation.keys() if gen not in critical_generators
]
# print(critical_generators, non_critical_generators)


# Update critical generators sets
critical_bus_gen = {bus: [] for bus in critical_buses}
critical_bus_thermalgen = {bus: [] for bus in critical_buses}
critical_bus_renewablegen = {bus: [] for bus in critical_buses}

for gen in critical_generators:
    bus_id = gen.split("_")[0]
    bus_name = f"bus{bus_id}"

    if bus_name in critical_buses:
        if gen in mod_object.model.generators:
            critical_bus_gen[bus_name].append(gen)

            if gen in mod_object.model.thermalGenerators:
                critical_bus_thermalgen[bus_name].append(gen)

            elif gen in mod_object.model.renewableGenerators:
                critical_bus_renewablegen[bus_name].append(gen)

for bus_name in mod_object_rel.model.criticalBuses:
    mod_object_rel.model.criticalGenerators[bus_name] = critical_bus_gen[bus_name]
    mod_object_rel.model.criticalthermalGenerators[bus_name] = critical_bus_thermalgen[
        bus_name
    ]
    mod_object_rel.model.criticalrenewableGenerators[bus_name] = (
        critical_bus_renewablegen[bus_name]
    )
    # print(bus_name, list(mod_object_rel.model.criticalGenerators[bus_name]),
    #       bus_name, list(mod_object_rel.model.criticalthermalGenerators[bus_name]),
    #       bus_name, list(mod_object_rel.model.criticalrenewableGenerators[bus_name])
    # )


# Update noncritical generators sets
noncritical_bus_gen = {bus: [] for bus in critical_buses}
noncritical_bus_thermalgen = {bus: [] for bus in critical_buses}
noncritical_bus_renewablegen = {bus: [] for bus in critical_buses}

for gen in non_critical_generators:
    bus_id = gen.split("_")[0]
    bus_name = f"bus{bus_id}"

    if bus_name in critical_buses:
        if gen in mod_object.model.generators:
            noncritical_bus_gen[bus_name].append(gen)

            if gen in mod_object.model.thermalGenerators:
                noncritical_bus_thermalgen[bus_name].append(gen)

            elif gen in mod_object.model.renewableGenerators:
                noncritical_bus_renewablegen[bus_name].append(gen)

for bus_name in mod_object_rel.model.criticalBuses:
    mod_object_rel.model.noncriticalthermalGenerators[bus_name] = (
        noncritical_bus_thermalgen[bus_name]
    )
    mod_object_rel.model.noncriticalrenewableGenerators[bus_name] = (
        noncritical_bus_renewablegen[bus_name]
    )
    # print(bus_name, list(mod_object_rel.model.noncriticalthermalGenerators[bus_name]),
    #       bus_name, list(mod_object_rel.model.noncriticalrenewableGenerators[bus_name])
    # )


# ----------------------- Step 3: Capacity failure states -------------- #
# Assign which critical generators are active in each failure state
failureStates = list(range(1, 9))
active_generator_state = {}
for bus in critical_buses:
    if len(critical_bus_gen[bus]) == 3:
        active_generator_state[bus, 1] = [
            critical_bus_gen[bus][0],
            critical_bus_gen[bus][1],
            critical_bus_gen[bus][2],
        ]
        active_generator_state[bus, 2] = [
            critical_bus_gen[bus][1],
            critical_bus_gen[bus][2],
        ]
        active_generator_state[bus, 3] = [
            critical_bus_gen[bus][0],
            critical_bus_gen[bus][2],
        ]
        active_generator_state[bus, 4] = [
            critical_bus_gen[bus][0],
            critical_bus_gen[bus][1],
        ]
        active_generator_state[bus, 5] = [critical_bus_gen[bus][2]]
        active_generator_state[bus, 6] = [critical_bus_gen[bus][1]]
        active_generator_state[bus, 7] = [critical_bus_gen[bus][0]]
        active_generator_state[bus, 8] = []

    elif len(critical_bus_gen[bus]) == 2:
        active_generator_state[bus, 1] = [
            critical_bus_gen[bus][0],
            critical_bus_gen[bus][1],
        ]
        active_generator_state[bus, 2] = [critical_bus_gen[bus][1]]
        active_generator_state[bus, 3] = [critical_bus_gen[bus][0]]
        active_generator_state[bus, 4] = []
        active_generator_state[bus, 5] = []
        active_generator_state[bus, 6] = []
        active_generator_state[bus, 7] = []
        active_generator_state[bus, 8] = []

    elif len(critical_bus_gen[bus]) == 1:
        active_generator_state[bus, 1] = [critical_bus_gen[bus][0]]
        active_generator_state[bus, 2] = []
        active_generator_state[bus, 3] = []
        active_generator_state[bus, 4] = []
        active_generator_state[bus, 5] = []
        active_generator_state[bus, 6] = []
        active_generator_state[bus, 7] = []
        active_generator_state[bus, 8] = []

    else:
        for state in failureStates:
            active_generator_state[bus, state] = []


active_critical_gen = {
    (bus, state): [] for bus in critical_buses for state in failureStates
}
active_critical_thermalgen = {
    (bus, state): [] for bus in critical_buses for state in failureStates
}
active_critical_renewablegen = {
    (bus, state): [] for bus in critical_buses for state in failureStates
}


for bus in critical_buses:
    for state in failureStates:
        if (bus, state) in active_generator_state:
            active_critical_gen[bus, state] = active_generator_state[bus, state]

for (bus, state), generators in active_critical_gen.items():
    thermal_gen = []
    renewable_gen = []

    for gen in generators:
        if gen in mod_object.model.thermalGenerators:
            thermal_gen.append(gen)
        elif gen in mod_object.model.renewableGenerators:
            renewable_gen.append(gen)

    active_critical_thermalgen[bus, state] = thermal_gen
    active_critical_renewablegen[bus, state] = renewable_gen
# print(active_critical_gen, active_critical_thermalgen, active_critical_renewablegen)

# Update sets for critical generators that are active at each failure state
for bus in mod_object_rel.model.criticalBuses:
    for state in mod_object_rel.model.states:
        mod_object_rel.model.activeCriticalthermalGenerators[bus, state] = (
            active_critical_thermalgen[bus, state]
        )
        mod_object_rel.model.activeCriticalrenewableGenerators[bus, state] = (
            active_critical_renewablegen[bus, state]
        )
        # print((bus, state), list(mod_object_rel.model.activeCriticalthermalGenerators[bus, state]),
        #       (bus, state), list(mod_object_rel.model.activeCriticalrenewableGenerators[bus, state]),
        # )


# Calculate the probability of capacity failure state based on the probability of failure
prob_state = {(bus, state): [] for bus in critical_buses for state in failureStates}
failure = {
    gen: mod_object_rel.model.failureRate[gen]
    for gen in mod_object_rel.model.generators
}

for bus in critical_buses:
    if len(critical_bus_gen[bus]) == 3:
        prob_state[bus, 1] = (
            (1 - failure[critical_bus_gen[bus][0]])
            * (1 - failure[critical_bus_gen[bus][1]])
            * (1 - failure[critical_bus_gen[bus][2]])
        )
        prob_state[bus, 2] = (
            failure[critical_bus_gen[bus][0]]
            * (1 - failure[critical_bus_gen[bus][1]])
            * (1 - failure[critical_bus_gen[bus][2]])
        )
        prob_state[bus, 3] = (
            (1 - failure[critical_bus_gen[bus][0]])
            * failure[critical_bus_gen[bus][1]]
            * (1 - failure[critical_bus_gen[bus][2]])
        )
        prob_state[bus, 4] = (
            (1 - failure[critical_bus_gen[bus][0]])
            * (1 - failure[critical_bus_gen[bus][1]])
            * failure[critical_bus_gen[bus][2]]
        )
        prob_state[bus, 5] = (
            failure[critical_bus_gen[bus][0]]
            * failure[critical_bus_gen[bus][1]]
            * (1 - failure[critical_bus_gen[bus][2]])
        )
        prob_state[bus, 6] = (
            failure[critical_bus_gen[bus][0]]
            * (1 - failure[critical_bus_gen[bus][1]])
            * failure[critical_bus_gen[bus][2]]
        )
        prob_state[bus, 7] = (
            (1 - failure[critical_bus_gen[bus][0]])
            * failure[critical_bus_gen[bus][1]]
            * failure[critical_bus_gen[bus][2]]
        )
        prob_state[bus, 8] = (
            failure[critical_bus_gen[bus][0]]
            * failure[critical_bus_gen[bus][1]]
            * failure[critical_bus_gen[bus][2]]
        )

    elif len(critical_bus_gen[bus]) == 2:
        prob_state[bus, 1] = (1 - failure[critical_bus_gen[bus][0]]) * (
            1 - failure[critical_bus_gen[bus][1]]
        )
        prob_state[bus, 2] = failure[critical_bus_gen[bus][0]] * (
            1 - failure[critical_bus_gen[bus][1]]
        )
        prob_state[bus, 3] = (1 - failure[critical_bus_gen[bus][0]]) * failure[
            critical_bus_gen[bus][1]
        ]
        prob_state[bus, 4] = (
            failure[critical_bus_gen[bus][0]] * failure[critical_bus_gen[bus][1]]
        )
        prob_state[bus, 5] = 1
        prob_state[bus, 6] = 1
        prob_state[bus, 7] = 1
        prob_state[bus, 8] = 1

    elif len(critical_bus_gen[bus]) == 1:
        prob_state[bus, 1] = 1 - failure[critical_bus_gen[bus][0]]
        prob_state[bus, 2] = failure[critical_bus_gen[bus][0]]
        prob_state[bus, 3] = 1
        prob_state[bus, 4] = 1
        prob_state[bus, 5] = 1
        prob_state[bus, 6] = 1
        prob_state[bus, 7] = 1
        prob_state[bus, 8] = 1

    else:
        for state in failureStates:
            prob_state[bus, state] = 1

# Update probability of failure of capacity failure state
for bus in mod_object_rel.model.criticalBuses:
    for state in mod_object_rel.model.states:
        mod_object_rel.model.prob[bus, state] = prob_state[bus, state]


# import pdb;pdb.set_trace()


# Solve reliability-constrained model
TransformationFactory("gdp.bound_pretransformation").apply_to(mod_object_rel.model)
TransformationFactory("gdp.bigm").apply_to(mod_object_rel.model)

mod_object_rel_opt = SolverFactory("gurobi")
mod_object_rel.results = mod_object_rel_opt.solve(mod_object_rel.model, tee=True)

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


# mod_object_rel.model.criticalBuses.pprint()
# mod_object_rel.model.noncriticalBuses.pprint()
# mod_object_rel.model.criticalGenerators.pprint()
# mod_object_rel.model.criticalthermalGenerators.pprint()
# mod_object_rel.model.criticalrenewableGenerators.pprint()
# mod_object_rel.model.noncriticalthermalGenerators.pprint()
# mod_object_rel.model.noncriticalrenewableGenerators.pprint()
# mod_object_rel.model.activeCriticalthermalGenerators.pprint()
# mod_object_rel.model.activeCriticalrenewableGenerators.pprint()

# mod_object_rel.model.prob.pprint()
# mod_object_rel.model.averageCapacityFactor.pprint()


for stage in mod_object_rel.model.investmentStage:
    for bus in mod_object_rel.model.criticalBuses:
        for state in mod_object_rel.model.states:
            print(
                "bus",
                bus,
                "state",
                state,
                "stage",
                stage,
                "production",
                mod_object_rel.model.investmentStage[stage]
                .prod_state[bus, state]
                .value,
            )
