from egret.data.model_data import ModelData
from egret.parsers import rts_gmlc_parser as parser
import pandas as pd
import os


def _read_branches(base_dir: str, elements: dict, bus_id_to_name: dict) -> None:

    # add the branches
    elements["branch"] = {}
    branch_df = pd.read_csv(os.path.join(base_dir, "branch.csv"))
    for idx, row in branch_df.iterrows():

        branch_dict = {
            "from_bus": bus_id_to_name[str(row["From Bus"])],
            "to_bus": bus_id_to_name[str(row["To Bus"])],
            "in_service": True,
            "resistance": float(row["R"]),
            "reactance": float(row["X"]),
            "charging_susceptance": float(row["B"]),
            "rating_long_term": float(row["Cont Rating"]) or None,
            "rating_short_term": float(row["LTE Rating"]) or None,
            "rating_emergency": float(row["STE Rating"]) or None,
            "angle_diff_min": -90,
            "angle_diff_max": 90,
            "pf": None,
            "qf": None,
            "pt": None,
            "qt": None,
        }

        TAP = float(row["Tr Ratio"])
        if TAP != 0.0:
            branch_dict["branch_type"] = "transformer"
            branch_dict["transformer_tap_ratio"] = TAP
            branch_dict["transformer_phase_shift"] = 0.0
        else:
            branch_dict["branch_type"] = "line"

        name = str(row["UID"])
        elements["branch"][name] = branch_dict
    branch_df = None


# data_path = "./gtep/data/123_Bus_Resil_Week"
data_path = "./gtep/data/123_Bus_Coal"

# from egret _create_rtsgmlc_skeleton
model_data = ModelData.empty_model_data_dict()

elements = model_data["elements"]
system = model_data["system"]

system["name"] = "RTS-GMLC"
system["baseMVA"] = 100.0

bus_id_to_name = parser._read_buses_and_areas(
    data_path, elements, system
)  # no error either 123
_read_branches(
    data_path, elements, bus_id_to_name
)  # no error 123_resil #error 123_coal
parser._read_generators(
    data_path, elements, bus_id_to_name
)  # error 123_resil (ramp_q key error)

pass
