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
input_data_path = (curr_dir / ".." / ".." / "data" / "5bus").resolve()
path_9_bus = (curr_dir / ".." / ".." / "data" / "9_bus_GTEP_dir").resolve()

bus_data_path = (
    curr_dir / ".." / ".." / "data" / "costs" / "Bus_data_gen_weights_mappings.csv"
).resolve()
cost_data_path = (
    curr_dir
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
    / "data"
    / "costs"
    / "Total_Energy_Supply_Disposition_and_Price_Summary.csv"
).resolve()


def create_model(
    input_data_path: Path = input_data_path,
    planning_data_args: dict | None = None,
    prescient_data_args: dict | None = None,
    config: dict | None = None,
    candidate_gens: list[str] | None = None,
    include_cost_data: bool = True,
):
    """
    :param input_data_path:             Path to the input data. Defaults to 5bus model
    :param planning_data_args:          Keyword arguments to pass to the constructor for `ExpansionPlanningData`. Defaults to `{}`
    :param prescient_data_args:         Keyword arguments to pass to `load_prescient`. Defaults to `{}`
    :param config:                      Dictionary of model config options. Defaults to `{}`
    :param candidate_gens:              List of candidate generators passed to `DataProcessing.load_gen_data`. Defaults to
                                            `["Natural Gas_FE", "Solar - Utility PV", "Land-Based Wind"]`
    :param include_cost_data:           Whether to build a `DataProcessing` object and pass to the `ExpansionPlanningModel` constructor.
                                            If `False`, `candidate_gens` is ignored. Defaults to `True`
    :type input_data_path:              Path, optional
    :type planning_data_args:           dict, optional
    :type config:                       dict, optional
    :type candidate_gens:               list[str], optional
    :type include_cost_data:            bool, optional
    """

    if planning_data_args is None:
        planning_data_args = {}

    if prescient_data_args is None:
        prescient_data_args = {}

    if config is None:
        config = {}

    if candidate_gens is None:
        candidate_gens = [
            "Natural Gas_FE",
            "Solar - Utility PV",
            "Land-Based Wind",
        ]

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
    data_object.load_prescient(input_data_path, **prescient_data_args)
    if "storage" in config and config["storage"]:
        data_object.load_storage_csv(str(input_data_path))

    if include_cost_data:
        data_processing_object = DataProcessing()
        data_processing_object.load_gen_data(
            bus_data_path=bus_data_path,
            cost_data_path=cost_data_path,
            ng_cost_path=ng_data_path,
            candidate_gens=candidate_gens,
            save_csv=False,
        )
        mod_object = ExpansionPlanningModel(
            config=config,
            data=data_object,
            cost_data=data_processing_object,
        )
    else:
        mod_object = ExpansionPlanningModel(config=config, data=data_object)

    mod_object.create_model()

    return mod_object
