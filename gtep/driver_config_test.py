from gtep.gtep_model import ExpansionPlanningModel
from gtep.gtep_data import ExpansionPlanningData
from gtep.gtep_solution import ExpansionPlanningSolution
from pyomo.core import TransformationFactory
from pyomo.environ import SolverFactory
from pyomo.contrib.appsi.solvers.highs import Highs
from pyomo.contrib.appsi.solvers.gurobi import Gurobi
from pyomo.contrib.appsi.solvers.ipopt import Ipopt
from icecream import ic


data_path = "./gtep/data/5bus_jsc"
data_object = ExpansionPlanningData()
data_object.load_prescient(data_path)

for key in data_object.md.data["elements"]["branch"]:
    data_object.md.data["elements"]["branch"][key]["capital_cost"] = 1000
    print(data_object.md.data["elements"]["branch"][key])

mod_object = ExpansionPlanningModel(
    stages=3,
    data=data_object.md,
    num_reps=1,
    len_reps=1,
    num_commit=24, # 24
    num_dispatch=6, # 4
)
mod_object.config["flow_model"] = "DC"
for k,v in mod_object.config.items():
    ic(k,v)

# quit()

mod_object.create_model()

#ic(mod_object)


#quit()
TransformationFactory("gdp.bound_pretransformation").apply_to(mod_object.model)
TransformationFactory("gdp.bigm").apply_to(mod_object.model)
opt = SolverFactory("gurobi")
opt = Gurobi()
#opt = SolverFactory("ipopt")
#opt = Ipopt()
opt.config.logfile = "logfileDC_5bus_CompleteNetwork.txt"
# # mod_object.results = opt.solve(mod_object.model, tee=True)
mod_object.results = opt.solve(mod_object.model)

sol_object = ExpansionPlanningSolution()
sol_object.load_from_model(mod_object)
sol_object.dump_json("./gtep_solution_DC_5busCompleteNetwork.json")

sol_object.import_data_object(data_object)

# sol_object.read_json("./gtep_lots_of_buses_solution.json")  # "./gtep/data/WECC_USAEE"
# sol_object.read_json("./gtep_11bus_solution.json")  # "./gtep/data/WECC_Reduced_USAEE"
# sol_object.read_json("./gtep_solution.json")
# sol_object.read_json("./updated_gtep_solution_test.json")
# sol_object.read_json("./gtep_wiggles.json")
# sol_object.plot_levels(save_dir="./plots/")

# save_numerical_results = False
# if save_numerical_results:

#     sol_object = ExpansionPlanningSolution()

#     sol_object.load_from_model(mod_object)
#     sol_object.dump_json()
# load_numerical_results = False
# if load_numerical_results:
# #     # sol_object.read_json("./gtep_solution.json")
#     sol_object.read_json("./bigger_longer_wigglier_gtep_solution.json")
# plot_results = False
# if plot_results:
#     sol_object.plot_levels(save_dir="./gtep/ACreactive_plots/")



# pass

# # import sys
# # sys.exit()
# stages_range = range(1, 4)
# num_commit_range = range(24, 61, 6)
# flow_models = ["ACR", "DC"]
# for stage in stages_range:
#     for num_commit in num_commit_range:
#         for flow_model in flow_models:
#             mod_object = ExpansionPlanningModel(
#                 stages=stage, # 1
#                 data=data_object.md,
#                 num_reps=1, # 1
#                 len_reps=1,
#                 num_commit=num_commit, # 24
#                 num_dispatch=4, # 4
#             )
#             mod_object.config["flow_model"] = flow_model
#             for k,v in mod_object.config.items():
#                 ic(k,v)

#             # quit()

#             short_flow = "AC" if flow_model == "ACR" else "DC"
#             instance_suffix = f"_S{stage}_NC{num_commit}"
#             logfile_name = f"logfile12{short_flow}5BusJSC_{instance_suffix}.txt"
#             json_filename = f"./gtep_solution_{short_flow}_5BusJSC{instance_suffix}.json"

#             mod_object.create_model()

#             #quit()
#             TransformationFactory("gdp.bound_pretransformation").apply_to(mod_object.model)
#             TransformationFactory("gdp.bigm").apply_to(mod_object.model)
#             opt = SolverFactory("gurobi")
#             opt = Gurobi()
#             #opt = SolverFactory("ipopt")
#             #opt = Ipopt()
#             opt.config.logfile = logfile_name
#             # # mod_object.results = opt.solve(mod_object.model, tee=True)
#             mod_object.results = opt.solve(mod_object.model)

#             sol_object = ExpansionPlanningSolution()
#             sol_object.load_from_model(mod_object)
#             sol_object.dump_json(json_filename)

#             sol_object.import_data_object(data_object)

#             # sol_object.read_json("./gtep_lots_of_buses_solution.json")  # "./gtep/data/WECC_USAEE"
#             # sol_object.read_json("./gtep_11bus_solution.json")  # "./gtep/data/WECC_Reduced_USAEE"
#             # sol_object.read_json("./gtep_solution.json")
#             # sol_object.read_json("./updated_gtep_solution_test.json")
#             # sol_object.read_json("./gtep_wiggles.json")
#             # sol_object.plot_levels(save_dir="./plots/")

#             # save_numerical_results = False
#             # if save_numerical_results:

#             #     sol_object = ExpansionPlanningSolution()

#             #     sol_object.load_from_model(mod_object)
#             #     sol_object.dump_json()
#             # load_numerical_results = False
#             # if load_numerical_results:
#             #     # sol_object.read_json("./gtep_solution.json")
#             #     sol_object.read_json("./bigger_longer_wigglier_gtep_solution.json")
#             plot_results = False
#             if plot_results:
#                 sol_object.plot_levels(save_dir="./gtep/ACreactive_plots/")



#             pass
