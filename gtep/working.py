from gtep.gtep_model import ExpansionPlanningModel
from gtep.gtep_data import ExpansionPlanningData
import pandas as pd
import os

data_path = "./gtep/data/123_Bus_Resil_Week"
data_object = ExpansionPlanningData(stages=3, num_reps=5, len_reps=24, num_commit=24, num_dispatch=1)


data_object.load_prescient(data_path)

load_scaling_path = data_path + "/ERCOT-Adjusted-Forecast.xlsb"
data_object.import_load_scaling(load_scaling_path)

data_object.texas_case_study_updates(data_path)

mod_object = ExpansionPlanningModel(data=data_object)

pass