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

"""This script runs only the GTEP analysis tools

This driver does not build or solve the GTEP optimization model.  It
is built for cases where a full GTEP run has already been completed
and the JSON result files produced by
`ExpansionPlanningSolution.save_results_in_json_files` are
available. When any of the required JSON files are missing the script
aborts with a clear `FileNotFoundError`. The script produces:

- Generation‑mix plots (piechart, treemap, or both)
- Stackgraph for representative days
- Quantitative analysis metrics
- Results for a specific given day (e.g., a single hour on a chosen date)

How to run:

    python driver_analysis.py --data-dir <path/to/input_csv_folder> --results-dir <path/to/json_results_folder>

"""

import argparse
import json
import logging
from pathlib import Path

import pyomo.environ as pyo

from gtep.gtep_model import ExpansionPlanningModel
from gtep.gtep_data import ExpansionPlanningData
from gtep.gtep_solution import ExpansionPlanningSolution
from gtep.gtep_data_processing import DataProcessing

logger = logging.getLogger("gtep.driver_analysis")
logger.setLevel(logging.INFO)

base_dir = Path(__file__).resolve().parent


def load_all_results(results_dir: Path) -> dict:
    """This function reads all JSON files written by the
    `save_results_in_json_files` function.  It raises an error if any
    of the required files is missing.

    """

    results_dir = results_dir.expanduser().resolve()

    # List of all the required files
    required = [
        "renewable_investments",
        "dispatchable_investments",
        "load_shed",
        "costs",
        "flows",
        "generation",
        "curtailment",
        "loads",
        "reserves",
        "charging",
        "discharging",
        "rep_days",
        "list_load_buses",
        "list_buses",
        "list_generators",
        "list_thermal_generators",
        "list_renewable_generators",
        "list_lines",
        "list_storage",
    ]

    data = {}
    for name in required:
        file_path = results_dir / f"{name}.json"
        if not file_path.is_file():
            raise FileNotFoundError(f"Required result file not found: {file_path}")
        with file_path.open("r", encoding="utf-8") as f:
            data[name] = json.load(f)

    return data


def main():

    parser = argparse.ArgumentParser(
        description="Run only the analysis / post‑processing tools for GTEP."
    )
    parser.add_argument(
        "--data-dir",
        required=True,
        help="Directory that contains the data CSV input files.",
    )
    parser.add_argument(
        "--results-dir",
        required=True,
        help="Directory that contains the results JSON files.",
    )
    args = parser.parse_args()

    data_path = Path(args.data_dir).resolve()
    results_path = Path(args.results_dir).resolve()

    logger.info("Loading JSON result files from %s", results_path)

    # Read all the results files and save the representative days
    # separately since we use it during the calculation of the
    # metrics.
    results = load_all_results(results_path)
    rep_days = results["rep_days"]

    # Create the solution object.
    sol_object = ExpansionPlanningSolution(data_path)

    # Plot generation‑mix. The options for `plot_type` are:
    # "piechart", "treemap", or "all". The options for `case_json`
    # are: "renewables", "dispatchables", or "combined".
    plot_type = "all"
    case_json = "combined"
    logger.info("Creating generation‑mix plots (%s, %s)", case_json, plot_type)
    sol_object.create_plots(
        case_json,
        results_path,
        data_path,
        plot_type,
    )

    # Stackgraph of the representative day(s)
    logger.info("Creating stack‑graph for %d representative days", len(rep_days))
    sol_object.create_stackgraph(results_path, rep_days)

    # Analysis metrics (CSV output set to true by default)
    logger.info("Running full analysis metrics")
    outputs, metrics_df = sol_object.create_analysis(
        results_path=results_path,
        save_csv=True,
        csv_path=results_path / "analysis_metrics.csv",
        hours_per_period=1.0,
        cost_units="million",
        data_dir=data_path,
    )
    logger.info("Analysis metrics saved to %s", results_path / "analysis_metrics.csv")

    # Save results for a single selected day and hour
    logger.info("Reading results for a single day / hour")
    sol_object.read_day_results(
        rep_days=rep_days,
        target_day="2020-04-23 00:00",
        target_hour=12,
        save_csv=True,
        results_path=results_path,
    )

    logger.info("Analysis run complete.")


if __name__ == "__main__":
    main()
