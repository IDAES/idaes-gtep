from gtep.gtep_model import ExpansionPlanningModel
from gtep.gtep_data_cho import ExpansionPlanningDataforReliability
from pyomo.core import TransformationFactory
from pyomo.contrib.appsi.solvers.highs import Highs
from pyomo.contrib.appsi.solvers.gurobi import Gurobi
from pyomo.environ import SolverFactory, value
import gtep_result_save

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


def solve_expansion_model(mod_object):
    mod_object.create_model()
    TransformationFactory("gdp.bound_pretransformation").apply_to(mod_object.model)
    TransformationFactory("gdp.bigm").apply_to(mod_object.model)
    mod_object_opt = SolverFactory("gurobi")
    mod_object.results = mod_object_opt.solve(mod_object.model, tee=True)

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
    noncritical_buses = list(set(mod_object.model.buses) - set(critical_buses))
    # print(critical_buses, noncritical_buses)

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

    total_renewable_generation = {
        gen: 0 for gen in mod_object.model.renewableGenerators
    }
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
    # print(average_capacity_factor)

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
    # print(critical_bus_gen, critical_bus_thermalgen, critical_bus_renewablegen)

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
    # print(noncritical_bus_gen, noncritical_bus_thermalgen, noncritical_bus_renewablegen)

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

    results_sets = {
        "criticalBuses": critical_buses,
        "noncriticalBuses": noncritical_buses,
        "criticalGenerators": critical_bus_gen,
        "criticalthermalGenerators": critical_bus_thermalgen,
        "criticalrenewableGenerators": critical_bus_renewablegen,
        "noncriticalthermalGenerators": noncritical_bus_thermalgen,
        "noncriticalrenewableGenerators": noncritical_bus_renewablegen,
        "activeCriticalthermalGenerators": active_critical_thermalgen,
        "activeCriticalrenewableGenerators": active_critical_renewablegen,
        "averageCapacityFactor": average_capacity_factor,
    }

    gtep_result_save.results_sets = results_sets

    return


# expansion_planning_results = solve_expansion_model(mod_object)
