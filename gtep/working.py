from gtep.gtep_model import ExpansionPlanningModel
from gtep.gtep_data import ExpansionPlanningData
import pandas as pd
import os
import pyomo.environ as pyo

data_path = "./gtep/data/123_Bus_Coal"
# data_path = "./gtep/data/5bus"

# df = pd.read_csv(os.path.join(data_path, "branch.csv"))

# new_df = df.astype({"From Bus": "int32", "To Bus": "int32"})

# new_path = os.path.join(data_path, "branch.csv")
# new_df.to_csv(new_path, header=True, index=False)
# pass

data_object = ExpansionPlanningData(
    stages=2, num_reps=4, len_reps=24, num_commit=24, num_dispatch=1
)

data_object.load_prescient(data_path)

load_scaling_path = data_path + "/ERCOT-Adjusted-Forecast.xlsb"
data_object.import_load_scaling(load_scaling_path)

data_object.texas_case_study_updates(data_path)

mod_object = ExpansionPlanningModel(data=data_object)

mod_object.create_model()

# CT 99.6
# CC 80.3
"""
Oakes, Matt, et al. "Cost and Performance Baseline for Fossil Energy Plants, Volume 5: Natural Gas Electricity Generating Units for Flexible Operation." , May. 2023. https://doi.org/10.2172/1973266
"""

# nuclear 0
# coal 8 (kyle will send reference)
pass
