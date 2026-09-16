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

            expected_files = [
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

            for filename in expected_files:
                file_path = Path(results_dir) / filename
                self.assertTrue(file_path.exists())

                with open(file_path, "r") as f:
                    data = json.load(f)

                self.assertIsInstance(data, dict)
