import pyomo.environ as pyo
from pyomo.contrib.appsi.solvers.highs import Highs
from pyomo.contrib.appsi.solvers.gurobi import Gurobi

from gtep.gtep_model import ExpansionPlanningModel
from gtep.gtep_data import ExpansionPlanningData
from gtep.gtep_solution import ExpansionPlanningSolution
from gtep.gtep_data_processing import DataProcessing

# Add data
data_path = "./gtep/data/5bus"
data_object = ExpansionPlanningData(
    stages=2,
    num_reps=2,
    len_reps=2,
    num_commit=2,
    num_dispatch=2,
    duration_dispatch=30,
)
data_object.load_prescient(data_path)

# [ESR WIP: Add bus and cost data files to be used on the
# DataProcessing class. This class processes data for the following
# type of generators: Natural Gas, Solar, Wind, and Coal. Note that
# during this data processing stage, the generator type is matched to
# existent generators in the data. The data contains the following
# types: (a) Natural Gas: Combustion Turbine (CT) and Fuel Efficiency
# (FE) and (b) Solar: Utility PV and Concentrated Solar Power (CSP)
bus_data_path = "./gtep/data/costs/Bus_data_gen_weights_mappings.csv"
cost_data_path = "./gtep/data/costs/2022_v3_Annual_Technology_Baseline_Workbook_Mid-year_update_2-15-2023_Clean.xlsx"
ng_cost_path = "./gtep/data/costs/Total_Energy_Supply_Disposition_and_Price_Summary.csv"
candidate_gens = ["Natural Gas_CT", "Natural Gas_FE", "Solar - Utility PV"]

data_processing_object = DataProcessing()
data_processing_object.load_gen_data(
    bus_data_path,
    cost_data_path,
    ng_cost_path,
    candidate_gens,
)

# Populate and create GTEP model
mod_object = ExpansionPlanningModel(data=data_object, cost_data=data_processing_object)

mod_object.config["include_investment"] = True
mod_object.config["include_commitment"] = True
mod_object.config["include_redispatch"] = True
mod_object.config["scale_loads"] = True
mod_object.config["transmission"] = True
mod_object.config["storage"] = False
mod_object.config["flow_model"] = "DC"

mod_object.create_model()
