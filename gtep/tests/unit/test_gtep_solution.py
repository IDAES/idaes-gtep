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

import json
from pathlib import Path

import pandas as pd
import pyomo.common.unittest as unittest
import pyomo.environ as pyo
from pyomo.environ import TransformationFactory
from pyomo.common.tempfiles import TempfileManager

from gtep.gtep_solution import ExpansionPlanningSolution
from gtep.gtep_model import ExpansionPlanningModel
from gtep.tests.unit.utils_for_testing import create_model, input_data_path


def get_solved_model():
    """This function creates and solves a small GTEP model to test
    GTEP solution class.

    """

    mod_object = create_model(
        planning_data_args={
            "stages": 1,
            "num_reps": 1,
            "num_commit": 1,
            "num_dispatch": 1,
            "duration_representative_period": 1,
        },
        prescient_data_args={
            "representative_dates": ["2020-01-28 00:00"],
        },
        include_cost_data=False,
    )

    TransformationFactory("gdp.bound_pretransformation").apply_to(mod_object.model)
    TransformationFactory("gdp.bigm").apply_to(mod_object.model)

    opt = pyo.SolverFactory("highs")

    if not opt.available():
        raise unittest.SkipTest("Solver not available")

    mod_object.results = opt.solve(mod_object.model)

    return mod_object


class TestExpansionPlanningSolution(unittest.TestCase):
    def test_solution_init_and_generation_types(self):
        # Test that the solution class initializes with the given
        # input data.
        sol_object = ExpansionPlanningSolution(input_data_path)

        self.assertIsInstance(sol_object.gen_df, pd.DataFrame)
        self.assertFalse(sol_object.gen_df.empty)

        self.assertIn("CC", sol_object.gen_types)
        self.assertIn("CT", sol_object.gen_types)
        self.assertIn("PV", sol_object.gen_types)
        self.assertIn("WIND", sol_object.gen_types)
        self.assertIn("BATTERY", sol_object.gen_types)
        self.assertIn("OTHER", sol_object.gen_types)

        self.assertTrue(hasattr(sol_object.gen_types["CC"], "label"))
        self.assertTrue(hasattr(sol_object.gen_types["CC"], "color"))

    def test_create_results_directory_and_save_model_config(self):
        # Test the creation of the results directory and the model
        # configuration CSV output file.
        mod_object = create_model(
            planning_data_args={
                "stages": 1,
                "num_reps": 1,
                "num_commit": 1,
                "num_dispatch": 1,
                "duration_representative_period": 1,
            },
            include_cost_data=False,
        )

        sol_object = ExpansionPlanningSolution(input_data_path)

        with TempfileManager.new_context() as tempfile:
            temp_dir = Path(tempfile.mkdtemp())
            results_dir = temp_dir / "results_test"

            dir_name = sol_object.create_results_directory(results_dir)

            self.assertTrue(Path(dir_name).exists())
            self.assertTrue(Path(dir_name).is_dir())

            sol_object.save_model_config_to_csv(mod_object, dir_name)

            config_file = Path(dir_name) / "model_config.csv"
            self.assertTrue(config_file.exists())

            config_df = pd.read_csv(config_file)
            self.assertIn("config_key", config_df.columns)
            self.assertIn("config_value", config_df.columns)
            self.assertIn("value_type", config_df.columns)

    def test_save_results_in_json_files(self):
        # Test that the results of the solved model are written to
        # JSON files.
        mod_object = get_solved_model()
        sol_object = ExpansionPlanningSolution(input_data_path)

        with TempfileManager.new_context() as tempfile:
            temp_dir = Path(tempfile.mkdtemp())
            results_dir = sol_object.create_results_directory(temp_dir / "results")

            sol_object.save_results_in_json_files(
                mod_object,
                results_dir,
                value_threshold=1e-3,
            )

            expected_files_with_dicts = [
                "renewable_investments.json",
                "dispatchable_investments.json",
                "load_shed.json",
                "costs.json",
                "flows.json",
                "generation.json",
                "curtailment.json",
                "loads.json",
                "reserves.json",
                "charging.json",
                "discharging.json",
            ]

            for filename in expected_files_with_dicts:
                file_path = Path(results_dir) / filename
                self.assertTrue(file_path.exists())

                with open(file_path, "r") as f:
                    data = json.load(f)

                self.assertIsInstance(data, dict)

            expected_files_with_lists = [
                "list_load_buses.json",
                "list_buses.json",
                "list_generators.json",
                "list_thermal_generators.json",
                "list_renewable_generators.json",
                "list_lines.json",
                "list_storage.json",
            ]

            for filename in expected_files_with_lists:
                file_path = Path(results_dir) / filename
                self.assertTrue(file_path.exists())

                with open(file_path, "r") as f:
                    data = json.load(f)

                self.assertIsInstance(data, list)

    def test_create_analysis_and_csv(self):
        # Test values saved using the create_analysis method.
        solved_mod = get_solved_model()
        sol_object = ExpansionPlanningSolution(input_data_path)

        # Create a temporary results directory and write all JSON
        # result files.
        with TempfileManager.new_context() as tempfile:
            temp_dir = Path(tempfile.mkdtemp())
            results_dir = sol_object.create_results_directory(temp_dir / "results")

            # Write the JSON files.
            sol_object.save_results_in_json_files(
                solved_mod,
                results_dir,
                value_threshold=1e-3,
            )

            # Run the create_analysis method and save .csv file.
            outputs, metrics_df = sol_object.create_analysis(
                results_path=results_dir,
                save_csv=True,
                csv_path=f"{results_dir}/analysis_metrics.csv",
                hours_per_period=1.0,
                cost_units="million",
                data_dir=input_data_path,
            )

            # Verify the CSV file exists
            csv_file = results_dir / "analysis_metrics.csv"
            self.assertTrue(
                csv_file.exists(), msg="analysis_metrics.csv was not created"
            )

            # Load the CSV and sanity‑check its structure
            df = pd.read_csv(csv_file)

            # Expected columns – exactly the ones the method builds
            expected_cols = {"metric", "category", "value", "units"}
            self.assertTrue(
                expected_cols.issubset(df.columns),
                msg=f"CSV missing columns: {expected_cols - set(df.columns)}",
            )

            # Check metrics that must be present
            key_metrics = {
                "buses",
                "load_buses",
                "branches",
                "existing_branches",
                "candidate_branches",
                "generators",
                "existing_generators",
                "candidate_generators",
                "thermal_generators",
                "existing_thermal_generators",
                "candidate_thermal_generators",
                "renewable_generators",
                "existing_renewable_generators",
                "candidate_renewable_generators",
                "storage",
                "existing_storage",
                "candidate_storage",
                "total_load",
                "total_load_shed",
                "total_generation",
                "total_generation",
                "total_curtailment",
                "total_storage_charging",
                "total_storage_discharging",
                "installed_generators",
                "installed_branches",
                "total_objective_cost",
                "total_operating_cost",
                "total_expansion_cost",
                "total_penalty_cost",
            }

            missing = key_metrics - set(df["metric"])
            self.assertFalse(
                missing,
                msg=f"The following expected metrics are absent from the CSV: {missing}",
            )

            # Expected values
            expected_values = [
                ("buses", "parameter_in_model", 5),
                ("load_buses", "parameter_in_model", 3),
                ("branches", "parameter_in_model", 7),
                ("existing_branches", "parameter_in_model", 7),
                ("candidate_branches", "parameter_in_model", 0),
                ("generators", "parameter_in_model", 12),
                ("existing_generators", "parameter_in_model", 8),
                ("candidate_generators", "parameter_in_model", 4),
                ("thermal_generators", "parameter_in_model", 4),
                ("existing_thermal_generators", "parameter_in_model", 4),
                ("candidate_thermal_generators", "parameter_in_model", 0),
                ("renewable_generators", "parameter_in_model", 8),
                ("existing_renewable_generators", "parameter_in_model", 4),
                ("candidate_renewable_generators", "parameter_in_model", 4),
                ("storage", "parameter_in_model", 0),
                ("existing_storage", "parameter_in_model", 0),
                ("candidate_storage", "parameter_in_model", 0),
            ]

            def _row_for(metric, category):
                matches = df[(df["metric"] == metric) & (df["category"] == category)]
                self.assertEqual(
                    len(matches),
                    1,
                    msg=f"Expected exactly one row for metric='{metric}', category='{category}'. "
                    f"Found {len(matches)}.",
                )
                return matches.iloc[0]

            # Perform the numeric assertions
            for metric, category, expected in expected_values:
                row = _row_for(metric, category)
                actual = row["value"]
                if isinstance(expected, float):
                    # Allow a tolerance for floating‑point rounding
                    self.assertAlmostEqual(
                        actual,
                        expected,
                        places=2,
                        msg=f"Metric '{metric}' (category={category}) value mismatch: "
                        f"expected {expected}, got {actual}",
                    )
                else:
                    # Integer values have an exact match
                    self.assertEqual(
                        actual,
                        expected,
                        msg=f"Metric '{metric}' (category={category}) value mismatch: "
                        f"expected {expected}, got {actual}",
                    )
