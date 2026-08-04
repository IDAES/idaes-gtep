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

import argparse
import logging
import tomllib
from pathlib import Path

import pyomo.environ as pyo
from pyomo.contrib.appsi.solvers.highs import Highs
from pyomo.contrib.appsi.solvers.gurobi import Gurobi

from gtep.gtep_model import ExpansionPlanningModel
from gtep.gtep_data import ExpansionPlanningData
from gtep.gtep_solution import ExpansionPlanningSolution
from gtep.gtep_data_processing import DataProcessing

logger = logging.getLogger("gtep.driver_config")
logger.setLevel(logging.INFO)

"""This script includes the main driver to run the GTEP model from a
TOML configuration file.

This driver builds and solves a Generation and Transmission Expansion
Planning (GTEP) model using options specified in a TOML configuration
file located in the `examples/` directory. The config file defines the
logging, data inputs, cost data inputs, model options, model
transformations, solver settings, and results directory.

Example on how to run it:
    python driver_from_config.py --config examples/config_5bus.toml

"""


def load_config(config_path):
    """Load a TOML configuration file."""
    config_path = Path(config_path).resolve()

    with open(config_path, "rb") as f:
        config = tomllib.load(f)

    return config


def main(config_path):

    # Load configuration file.
    config = load_config(config_path)

    # Configure logging from the config file. Supported levels include
    # DEBUG, INFO, WARNING, and ERROR; defaults to INFO if not
    # provided.
    logging.basicConfig(
        level=getattr(logging, config.get("logging", {}).get("level", "INFO")),
        format="%(levelname)s:%(name)s:%(message)s",
    )

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
        len_reps=data_config.get("len_reps", 1),
        num_commit=data_config.get("num_commit", 6),
        num_dispatch=data_config.get("num_dispatch", 4),
        duration_dispatch=data_config.get("duration_dispatch", 15),
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

    mod_object.config["include_investment"] = model_config.get(
        "include_investment", True
    )
    mod_object.config["include_commitment"] = model_config.get(
        "include_commitment", True
    )
    mod_object.config["include_redispatch"] = model_config.get(
        "include_redispatch", True
    )
    mod_object.config["scale_loads"] = model_config.get("scale_loads", True)
    mod_object.config["transmission"] = model_config.get("transmission", True)
    mod_object.config["storage"] = model_config.get("storage", False)
    mod_object.config["flow_model"] = model_config.get("flow_model", "DC")
    mod_object.config["advanced_hydro"] = model_config.get("advanced_hydro", False)

    # Save the model configuration settings to a CSV file.
    sol_object.save_model_config_to_csv(mod_object, dir_name)

    # Create model.
    mod_object.create_model()

    # -------------------------------------------------------------------
    # Apply transformations to logical terms and solve the model.
    transformation_config = config.get("transformations", {})

    if transformation_config.get("bound_pretransformation", False):
        pyo.TransformationFactory("gdp.bound_pretransformation").apply_to(
            mod_object.model
        )

    if transformation_config.get("bigm", True):
        pyo.TransformationFactory("gdp.bigm").apply_to(mod_object.model)

    # Add solver. Solver options include "gurobi" and "highs". Default
    # is "gurobi", if none is provided.
    solver_config = config.get("solver", {})
    solver_name = solver_config.get("name", "gurobi")
    tee = solver_config.get("tee", True)

    opt = pyo.SolverFactory(solver_name)
    mod_object.results = opt.solve(mod_object.model, tee=tee)

    print(mod_object.results)

    # -------------------------------------------------------------------
    # Save results in JSON files using the GTEP solution class.
    value_threshold = results_config.get("value_threshold", 1e-3)
    sol_object.save_results_in_json_files(
        mod_object,
        dir_name,
        value_threshold=value_threshold,
    )


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="Run GTEP from a TOML config file.")
    parser.add_argument(
        "--config",
        required=True,
        help="Path to the TOML configuration file.",
    )

    args = parser.parse_args()
    main(args.config)
