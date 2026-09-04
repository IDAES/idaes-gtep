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

"""This script gives a pure Python driver for running the GTEP model
without needing a separate configuration file.

This driver uses all the main model configuration options directly in
one place. It is useful for debugging, development, and small example
runs where editing a Python script is preferred over using a
configuration file.

"""

import logging
from pathlib import Path

import pyomo.environ as pyo

from gtep.gtep_model import ExpansionPlanningModel
from gtep.gtep_data import ExpansionPlanningData
from gtep.gtep_solution import ExpansionPlanningSolution
from gtep.gtep_data_processing import DataProcessing

logger = logging.getLogger("gtep.driver")
logger.setLevel(logging.INFO)

base_dir = Path(__file__).resolve().parent

# ---------------------------------------------------------------------
# Add data path

case_name = "5bus"
data_path = base_dir / "data" / case_name

# ---------------------------------------------------------------------
# Create directory to save results using the GTEP solution class

sol_object = ExpansionPlanningSolution(data_path)
dir_name = sol_object.create_results_directory(f"results_{case_name}")

# ---------------------------------------------------------------------
# Create data modeling object

data_object = ExpansionPlanningData(
    stages=2,
    num_reps=4,
    num_commit=6,
    num_dispatch=4,
    duration_representative_period=6,
)

data_object.load_prescient(
    data_path,
    representative_dates=[
        "2020-01-28 00:00",
        "2020-04-23 00:00",
        "2020-07-05 00:00",
        "2020-10-14 00:00",
    ],
    representative_weights={
        "2020-01-28 00:00": 115,
        "2020-04-23 00:00": 95,
        "2020-07-05 00:00": 50,
        "2020-10-14 00:00": 105,
    },
)

# ---------------------------------------------------------------------
# Create cost data processing object

include_cost_data = True

if include_cost_data:
    bus_data_path = base_dir / "data" / "costs" / "Bus_data_gen_weights_mappings.csv"
    cost_data_path = (
        base_dir
        / "data"
        / "costs"
        / "2022_v3_Annual_Technology_Baseline_Workbook_Mid-year_update_2-15-2023_Clean.xlsx"
    )
    ng_cost_path = (
        base_dir
        / "data"
        / "costs"
        / "Total_Energy_Supply_Disposition_and_Price_Summary.csv"
    )
    candidate_gens = [
        "Natural Gas_CT",
        "Natural Gas_FE",
        "Solar - Utility PV",
    ]

    data_processing_object = DataProcessing()
    data_processing_object.load_gen_data(
        bus_data_path,
        cost_data_path,
        ng_cost_path,
        candidate_gens,
    )
else:
    data_processing_object = None

# ---------------------------------------------------------------------
# Populate and configure GTEP model

mod_object = ExpansionPlanningModel(
    data=data_object,
    cost_data=data_processing_object,
)

# Core model options
mod_object.config["include_investment"] = True
mod_object.config["include_commitment"] = True
mod_object.config["include_redispatch"] = True
mod_object.config["flow_model"] = "DC"

# Time-period options
mod_object.config["time_period_subsets"] = []
mod_object.config["time_period_dict"] = {}
mod_object.config["dispatch_randomization"] = True

# Common options
mod_object.config["scale_loads"] = False
mod_object.config["scale_texas_loads"] = False

# Investment options
mod_object.config["thermal_generation"] = False
mod_object.config["renewable_generation"] = False
mod_object.config["storage"] = False
mod_object.config["transmission"] = True
mod_object.config["transmission_switching"] = False
mod_object.config["advanced_hydro"] = False

# ---------------------------------------------------------------------
# Save model configuration

sol_object.save_model_config_to_csv(mod_object, dir_name)

# ---------------------------------------------------------------------
# Create model

logger.info("Creating GTEP model.")
mod_object.create_model()

# ---------------------------------------------------------------------
# Apply GDP transformations

# Options include "bigm", "chull", etc.. "bigm" is set as the default.
apply_bound_pretransformation = True
gdp_transform = "bigm"

if apply_bound_pretransformation:
    pyo.TransformationFactory("gdp.bound_pretransformation").apply_to(mod_object.model)

pyo.TransformationFactory(f"gdp.{gdp_transform}").apply_to(mod_object.model)

# ---------------------------------------------------------------------
# Add solver

# Select Pyomo-compatible MILP solver for the transformed GDP
# model. highs is the default solver. Other options are gurobi,
# xpress, etc..
solver_name = "highs"
tee = True

opt = pyo.SolverFactory(solver_name)
mod_object.results = opt.solve(mod_object.model, tee=tee)

print(mod_object.results)

# ---------------------------------------------------------------------
# Save results in JSON files and create plots using the GTEP solution
# class

sol_object.save_results_in_json_files(
    mod_object,
    dir_name,
    value_threshold=1e-3,
)

# Create generation-mix plots from the saved investment JSON files.
# Set case_json to "dispatchables" or "renewables" to plot each group
# separately, or use "combined" to merge both groups in the same
# plots. Set plot_type to "piechart", "treemap", or "all"; "all"
# generates both pie chart and treemap plots.
plot_type = "all"
case_json = "combined"
sol_object.create_plots(case_json, dir_name, data_path, plot_type)


# Create stackgraph
rep_days = [
    pyo.value(mod_object.model.representativeDate[rep])
    for rep in mod_object.model.representativeDate
]
sol_object.create_stackgraph(dir_name, rep_days)

logger.info("GTEP run complete.")
