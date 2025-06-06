from gtep.gtep_model import ExpansionPlanningModel
from gtep.gtep_data import ExpansionPlanningData
from gtep.gtep_solution import ExpansionPlanningSolution
from gtep_model import ExpansionPlanningModel
from gtep_data import ExpansionPlanningData
from gtep_data_processing import DataProcessing
from gtep_solution import ExpansionPlanningSolution
from pyomo.core import TransformationFactory
from pyomo.contrib.appsi.solvers.highs import Highs
from pyomo.contrib.appsi.solvers.gurobi import Gurobi


data_path = "./data/5bus"
data_object = ExpansionPlanningData()
data_object.load_prescient(data_path)

# [ESR WIP: Add bus and cost data files to be used on the
# DataProcessing class. This class processes data for the following
# type of generators: Natural Gas, Solar, Wind, and Coal. Note that
# during this data processing stage, the generator type is matched to
# existent generators in the data. The data contains the following
# types: (a) Natural Gas: Combustion Turbine (CT) and Fuel Efficiency
# (FE) and (b) Solar: Utility PV and Concentrated Solar Power (CSP)

bus_data_path = "data/costs/Bus_data_gen_weights_mappings.csv"
cost_data_path = "data/costs/2022 v3 Annual Technology Baseline Workbook Mid-year update 2-15-2023 Clean.xlsx"
candidate_gens = ["Natural Gas_CT", "Natural Gas_FE", "Solar - Utility PV", "Land-Based Wind"]

data_processing_object = DataProcessing()
data_processing_object.load_gen_data(bus_data_path=bus_data_path,
                                     cost_data_path=cost_data_path,
                                     candidate_gens=candidate_gens,
                                     save_csv=True)

mod_object = ExpansionPlanningModel(
    stages=2,
    data=data_object.md,
    cost_data=data_processing_object,
    num_reps=2,
    len_reps=1,
    num_commit=6,
    num_dispatch=4,
)

mod_object.create_model()

TransformationFactory("gdp.bound_pretransformation").apply_to(mod_object.model)
TransformationFactory("gdp.bigm").apply_to(mod_object.model)
# opt = SolverFactory("gurobi")
opt = Gurobi()
# opt = Highs()
# mod_object.results = opt.solve(mod_object.model, tee=True)
mod_object.results = opt.solve(mod_object.model)
print(mod_object.results)
# mod_object.model.investmentStage.pprint()
# mod_object.report_model()

quit()

sol_object = ExpansionPlanningSolution()
sol_object.load_from_model(mod_object)
sol_object.dump_json("./gtep_solution.json")
sol_object.import_data_object(data_object)

# sol_object.read_json("./gtep_lots_of_buses_solution.json")  # "./gtep/data/WECC_USAEE"
# sol_object.read_json("./gtep_11bus_solution.json")  # "./gtep/data/WECC_Reduced_USAEE"
# sol_object.read_json("./gtep_solution.json")
# sol_object.read_json("./updated_gtep_solution_test.json")
# sol_object.read_json("./gtep_wiggles.json")
sol_object.plot_levels(save_dir="./plots/")

# save_numerical_results = False
# if save_numerical_results:

#     sol_object = ExpansionPlanningSolution()

#     sol_object.load_from_model(mod_object)
#     sol_object.dump_json()
# load_numerical_results = False

# if load_numerical_results:
#     # sol_object.read_json("./gtep_solution.json")
#     sol_object.read_json("./bigger_longer_wigglier_gtep_solution.json")
# plot_results = False

# if plot_results:
#     sol_object.plot_levels(save_dir="./plots/")


pass
