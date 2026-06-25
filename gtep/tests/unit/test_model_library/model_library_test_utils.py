#################################################################################
# The Institute for the Design of Advanced Energy Systems Integrated Platform
# Framework (IDAES IP) was produced under the DOE Institute for the
# Design of Advanced Energy Systems (IDAES).
#
# Copyright (c) 2018-2026 by the software owners: The Regents of the
# University of California, through Lawrence Berkeley National Laboratory,
# National Technology & Engineering Solutions of Sandia, LLC, Carnegie Mellon
# University, West Virginia University Research Corporation, et al.
# All rights reserved.  Please see the files COPYRIGHT.md and LICENSE.md
# for full copyright and license information.
#################################################################################
from gtep.gtep_model import ExpansionPlanningModel
from gtep.gtep_data import ExpansionPlanningData
from gtep.gtep_data_processing import DataProcessing
from pathlib import Path

curr_dir = Path(__file__).resolve().parent
input_data_path = (curr_dir / ".." / ".." / ".." / "data" / "5bus").resolve()

bus_data_path = (
    curr_dir
    / ".."
    / ".."
    / ".."
    / "data"
    / "costs"
    / "Bus_data_gen_weights_mappings.csv"
).resolve()
cost_data_path = (
    curr_dir
    / ".."
    / ".."
    / ".."
    / "data"
    / "costs"
    / "2022_v3_Annual_Technology_Baseline_Workbook_Mid-year_update_2-15-2023_Clean.xlsx"
).resolve()
ng_data_path = (
    curr_dir
    / ".."
    / ".."
    / ".."
    / "data"
    / "costs"
    / "Total_Energy_Supply_Disposition_and_Price_Summary.csv"
).resolve()


def create_model(input_data_path=input_data_path, planning_data_args={}, config={}):
    default_data_planning_args = dict(
        stages=1,
        num_reps=1,
        num_commit=1,
        num_dispatch=1,
    )
    for arg, val in default_data_planning_args.items():
        if arg not in planning_data_args:
            planning_data_args[arg] = val

    data_object = ExpansionPlanningData(**planning_data_args)
    data_object.load_prescient(input_data_path)
    if "storage" in config and config["storage"]:
        data_object.load_storage_csv(str(input_data_path))

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

    return mod_object.model
