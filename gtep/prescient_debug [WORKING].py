from egret.data.model_data import ModelData
from egret.parsers import rts_gmlc_parser as parser
import pandas as pd
import os

data_path = "./gtep/data/123_Bus_Resil_Week"
# data_path = "./gtep/data/123_Bus_Coal"
# data_path = "./gtep/data/5bus"

### from egret _create_rtsgmlc_skeleton ###
model_data = ModelData.empty_model_data_dict()

elements = model_data["elements"]
system = model_data["system"]

system["name"] = "RTS-GMLC"
system["baseMVA"] = 100.0

bus_id_to_name = parser._read_buses_and_areas(
    data_path, elements, system
)  # no error either 123
parser._read_branches(
    data_path, elements, bus_id_to_name
)  # no error 123_resil #error 123_coal
parser._read_generators(
    data_path, elements, bus_id_to_name
)  # error 123_resil (ramp_q key error)

pass
