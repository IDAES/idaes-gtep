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

"""This script runs the GTEP model from a TOML configuration file.

This driver builds and solves a Generation and Transmission Expansion
Planning (GTEP) model using options specified in a TOML configuration
file located in the `examples/` directory. The config file defines the
logging, data inputs, cost data inputs, model options, model
transformations, solver settings, and results directory.

Example on how to run it:
    python driver_from_config.py --config examples/config_5bus.toml

"""

import argparse
import logging
from pathlib import Path

try:
    import tomllib
except ModuleNotFoundError:
    import tomli as tomllib

import pyomo.environ as pyo

from gtep.gtep_model import ExpansionPlanningModel
from gtep.gtep_data import ExpansionPlanningData
from gtep.gtep_solution import ExpansionPlanningSolution
from gtep.gtep_data_processing import DataProcessing

logger = logging.getLogger("gtep.driver_from_config")
logger.setLevel(logging.INFO)


def load_config(config_path):
    """This function loads a TOML configuration file."""
    config_path = Path(config_path).resolve()

    with open(config_path, "rb") as f:
        config = tomllib.load(f)

    return config


def main(config_path):
    """This function runs the GTEP model using options from a TOML
    config file.

    """

    # -------------------------------------------------------------------
    # Load configuration file.
    config = load_config(config_path)

    # Configure logging from the config file. Supported levels include
    # DEBUG, INFO, WARNING, and ERROR; defaults to INFO if not
    # provided.
    logging.basicConfig(
        level=getattr(logging, config.get("logging", {}).get("level", "INFO")),
        format="%(levelname)s:%(name)s:%(message)s",
    )

    # -------------------------------------------------------------------
    # Add data path from configuration file.
    data_config = config["data"]
    data_path = data_config["data_path"]

    # -------------------------------------------------------------------
    # Create directory to save results using the GTEP solution class.
    sol_object = ExpansionPlanningSolution(data_path)

    results_config = config.get("results", {})
    case_name = Path(data_path).name
    results_base_name = results_config.get("dir_name", "results")
    dir_name = sol_object.create_results_directory(f"{results_base_name}_{case_name}")

    # -------------------------------------------------------------------
    # Create the data object using values from the config file. If a
    # value is not provided, use the default values included below.
    data_object = ExpansionPlanningData(
        stages=data_config.get("stages", 2),
        num_reps=data_config.get("num_reps", 2),
        num_commit=data_config.get("num_commit", 6),
        num_dispatch=data_config.get("num_dispatch", 4),
        duration_representative_period=data_config.get(
            "duration_representative_period", 6
        ),
    )

    # Add optional representative-period inputs if they are provided
    # in the configuration file. If they are not provided,
    # load_prescient uses the default representative-period selection.
    load_prescient_kwargs = {}

    if "representative_dates" in data_config:
        load_prescient_kwargs["representative_dates"] = data_config[
            "representative_dates"
        ]

    if "representative_weights" in data_config:
        load_prescient_kwargs["representative_weights"] = data_config[
            "representative_weights"
        ]

    if "start_date" in data_config:
        load_prescient_kwargs["start_date"] = data_config["start_date"]

    data_object.load_prescient(data_path, **load_prescient_kwargs)

    # -------------------------------------------------------------------
    # Create cost data processing object.
    cost_config = config.get("cost_data", {})

    if cost_config.get("enabled", True):
        bus_data_path = cost_config["bus_data_path"]
        cost_data_path = cost_config["cost_data_path"]
        ng_cost_path = cost_config["ng_cost_path"]
        candidate_gens = cost_config["candidate_gens"]

        data_processing_object = DataProcessing()
        data_processing_object.load_gen_data(
            bus_data_path,
            cost_data_path,
            ng_cost_path,
            candidate_gens,
        )
    else:
        data_processing_object = None

    # -------------------------------------------------------------------
    # Populate and create GTEP model.
    mod_object = ExpansionPlanningModel(
        data=data_object,
        cost_data=data_processing_object,
    )

    # Apply model options from the config file. If an option is not
    # provided, use the default value shown below.
    model_config = config.get("model", {})

    model_defaults = {
        # Core model options
        "include_investment": True,
        "include_commitment": True,
        "include_redispatch": True,
        "flow_model": "DC",
        # Time-period options
        "time_period_subsets": [],
        "time_period_dict": {},
        "dispatch_randomization": False,
        # Common options
        "scale_loads": True,
        "scale_texas_loads": False,
        # Investment options
        "thermal_generation": False,
        "renewable_generation": False,
        "storage": False,
        "transmission": False,
        "transmission_switching": False,
        "advanced_hydro": False,
    }
    for key, default_value in model_defaults.items():
        mod_object.config[key] = model_config.get(key, default_value)

    # Save the model configuration settings to a CSV file.
    sol_object.save_model_config_to_csv(mod_object, dir_name)

    # -------------------------------------------------------------------
    # Create model and apply transformations to logical terms and
    # solve the model.

    mod_object.create_model()

    transformation_config = config.get("transformations", {})

    if transformation_config.get("bound_pretransformation", False):
        pyo.TransformationFactory("gdp.bound_pretransformation").apply_to(
            mod_object.model
        )

    gdp_transform = transformation_config.get("gdp_transform", "bigm")
    pyo.TransformationFactory(f"gdp.{gdp_transform}").apply_to(mod_object.model)

    # Solver options include "gurobi" and "highs". Default
    # is "highs", if none is provided.
    solver_config = config.get("solver", {})
    solver_name = solver_config.get("solver_name", "highs")
    tee = solver_config.get("tee", True)

    opt = pyo.SolverFactory(solver_name)
    mod_object.results = opt.solve(mod_object.model, tee=tee)

    # Print solver results stats
    print(mod_object.results)

    # -------------------------------------------------------------------
    # Save results in JSON files using the GTEP solution class.
    if results_config.get("save_json", True):
        value_threshold = results_config.get("value_threshold", 1e-3)

        sol_object.save_results_in_json_files(
            mod_object,
            dir_name,
            value_threshold=value_threshold,
        )

    # Create generation-mix plots and/or stackgraph if requested.
    plot_config = config.get("plots", {})

    if plot_config.get("enabled", False):
        plot_type = plot_config.get("plot_type", "all")
        case_json = plot_config.get("case_json", "combined")

        sol_object.create_plots(
            case_json,
            dir_name,
            data_path,
            plot_type,
        )

    if plot_config.get("stackgraph", False):
        rep_days = [
            pyo.value(mod_object.model.representativeDate[rep])
            for rep in mod_object.model.representativeDate
        ]

        sol_object.create_stackgraph(dir_name, rep_days)

    logger.info("GTEP run complete.")


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="Run GTEP from a TOML config file.")
    parser.add_argument(
        "--config",
        required=True,
        help="Path to the TOML configuration file.",
    )

    args = parser.parse_args()
    main(args.config)
