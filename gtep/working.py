from gtep.model_library.investment import (
    add_investment_params_and_variables,
    add_investment_disjuncts,
    add_investment_constraints,
)
import pyomo.environ as pyo
from pyomo.environ import units as u
from gtep.gtep_model import ExpansionPlanningModel
from gtep.gtep_data import ExpansionPlanningData
from gtep.gtep_data_processing import DataProcessing
from pathlib import Path

curr_dir = Path(__file__).resolve().parent
input_data_source = (
    curr_dir / "data" / "9_bus_GTEP_dir"
).resolve()

bus_data_path = (
    curr_dir
    / "data"
    / "costs"
    / "Bus_data_gen_weights_mappings.csv"
).resolve()
cost_data_path = (
    curr_dir
    / "data"
    / "costs"
    / "2022_v3_Annual_Technology_Baseline_Workbook_Mid-year_update_2-15-2023_Clean.xlsx"
).resolve()
ng_data_path = (
    curr_dir
    / "data"
    / "costs"
    / "Total_Energy_Supply_Disposition_and_Price_Summary.csv"
).resolve()

config={"storage": True, "transmission": True}

def create_object(config):
    data_object = ExpansionPlanningData(
        stages=1,
        num_reps=1,
        num_commit=1,
        num_dispatch=1,
    )
    data_object.load_prescient(input_data_source)
    data_object.load_storage_csv(str(input_data_source))

    candidate_gens = [
        "Natural Gas_FE",
        "Solar - Utility PV",
        "Land-Based Wind",
    ]

    data_processing_object = DataProcessing()
    data_processing_object.load_gen_data(
        bus_data_path=bus_data_path,
        cost_data_path=cost_data_path,
        ng_cost_path=ng_data_path,
        candidate_gens=candidate_gens,
        save_csv=False,
    )

    mod_object = ExpansionPlanningModel(
        data=data_object,
        cost_data=data_processing_object,
    )

    for config_option, config_val in config.items():
        mod_object.config[config_option] = config_val

    mod_object.create_model()
    return mod_object

m = create_object(config).model
# self.b = self.m.investmentStage[self.investment_stage]

pass