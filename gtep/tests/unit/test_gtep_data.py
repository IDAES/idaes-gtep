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


import pyomo.common.unittest as unittest
from unittest.mock import create_autospec
import pytest
from gtep.gtep_data import ExpansionPlanningData
import prescient.simulator.config
import pandas as pd
from pathlib import Path

curr_dir = Path(__file__).resolve().parent
input_data_source = (curr_dir / ".." / ".." / "data" / "5bus").resolve()

load_scaling_file = (
    curr_dir / ".." / ".." / "data" / "Texas_2000" / "ERCOT-Adjusted-Forecast.xlsb"
).resolve()

storage_file = (curr_dir / ".." / ".." / "data" / "9_bus_GTEP_dir").resolve()

texas_data_path = (curr_dir / ".." / ".." / "data" / "Texas_2000").resolve

outage_data_path = (
    curr_dir / ".." / ".." / "data" / "123_Bus_Resil_Week" / "may_20.csv"
).resolve()


# Function Mocks (Load Prescient)
@pytest.fixture
def mock_prescient_config():
    mock_instance = create_autospec(
        prescient.simulator.config.PrescientConfig, instance=True
    )
    mock_instance.set_value.return_value = None
    mock_instance.num_days = 365
    mock_instance.sced_frequency_minutes = 60

    with unittest.mock.patch(
        "gtep.gtep_data.PrescientConfig", return_value=mock_instance
    ):
        yield mock_instance


class TestExpansionPlanningData(unittest.TestCase):

    def test_data_init(self):
        # Test that the ExpansionPlanningData object initializes properly with default values
        testObject = ExpansionPlanningData()
        self.assertIsInstance(testObject, ExpansionPlanningData)
        self.assertEqual(testObject.stages, 2)
        self.assertEqual(testObject.num_reps, 4)
        self.assertEqual(testObject.len_reps, 1)
        self.assertEqual(testObject.num_commit, 24)
        self.assertEqual(testObject.num_dispatch, 1)
        self.assertEqual(testObject.duration_dispatch, 60)

        # Test that the ExpansionPlanningData object initializes properly with input values
        testObject = ExpansionPlanningData(1, 2, 2, 2, 2, 15)
        self.assertEqual(testObject.stages, 1)
        self.assertEqual(testObject.num_reps, 2)
        self.assertEqual(testObject.len_reps, 2)
        self.assertEqual(testObject.num_commit, 2)
        self.assertEqual(testObject.num_dispatch, 2)
        self.assertEqual(testObject.duration_dispatch, 15)

        # Test that the ExpansionPlanningData object initializes properly with partial input values
        testObject = ExpansionPlanningData(duration_dispatch=15)
        self.assertEqual(testObject.stages, 2)
        self.assertEqual(testObject.num_reps, 4)
        self.assertEqual(testObject.len_reps, 1)
        self.assertEqual(testObject.num_commit, 24)
        self.assertEqual(testObject.num_dispatch, 1)
        self.assertEqual(testObject.duration_dispatch, 15)

    # -------------------------------------------------LOAD_PRESCIENT------------------------------------------------------------ #
    def test_options_dict(self, mock_prescient_config):
        # Test passing in an options dictionary
        testObject = ExpansionPlanningData()
        # set options that are not the defaults
        options = {
            "num_days": 100,
            "ruc_horizon": 12,
        }

        testObject.load_prescient(data_path=input_data_source, options_dict=options)

        mock_prescient_config.set_value.assert_called_once()
        passed_dict = mock_prescient_config.set_value.call_args[0][0]
        self.assertEqual(passed_dict["data_path"], str(input_data_source))
        self.assertEqual(passed_dict["num_days"], 100)
        self.assertEqual(passed_dict["ruc_horizon"], 12)

    def test_no_options_dict(self, mock_prescient_config):

        # Test not passing in an options dictionary
        testObject = ExpansionPlanningData()
        testObject.load_prescient(data_path=input_data_source)
        mock_prescient_config.set_value.assert_called_once()
        passed_dict = mock_prescient_config.set_value.call_args[0][0]
        # Default options should be set
        self.assertEqual(passed_dict["data_path"], str(input_data_source))
        self.assertEqual(passed_dict["num_days"], 365)
        self.assertEqual(passed_dict["ruc_horizon"], 36)

    def test_default_representative_dates(self):
        # Test no representative dates passed in, initializing with defaults
        testObject = ExpansionPlanningData()
        testObject.load_prescient(data_path=input_data_source)
        # default dates:
        expected_dates = [
            "2020-01-28 00:00",
            "2020-04-23 00:00",
            "2020-07-05 00:00",
            "2020-10-14 00:00",
        ]
        self.assertEqual(testObject.representative_dates, expected_dates)

    def test_passed_representative_dates(self):
        # Test new representative dates passed in, replacing defaults
        expected_dates = [
            "2020-01-28 00:00",
            "2020-04-23 00:00",
            "2020-07-05 00:00",
            "2020-12-14 00:00",
        ]

        testObject = ExpansionPlanningData()
        testObject.load_prescient(
            data_path=input_data_source,
            representative_dates=expected_dates,
        )

        self.assertEqual(testObject.representative_dates, expected_dates)

    def test_representative_date_not_in_time_keys(self):
        # Test passing invalid/not covered dates
        testObject = ExpansionPlanningData()
        bad_dates = [
            "2020-01-28 00:00",
            "2020-04-23 00:00",
            "2020-07-05 00:00",
            "2099-01-01 00:00",  # Not in time_keys
        ]
        with self.assertRaises(ValueError):
            testObject.load_prescient(
                data_path=input_data_source, representative_dates=bad_dates
            )

    def test_empty_representative_dates(self):
        # Test an empty list passed to overwrite defaults without setting new dates
        testObject = ExpansionPlanningData()
        with self.assertRaises(ZeroDivisionError):
            testObject.load_prescient(
                data_path=input_data_source, representative_dates=[]
            )

    def test_default_representative_weights(self):
        # Test no representative weights passed in, so the function calculates them
        testObject = ExpansionPlanningData()
        testObject.load_prescient(data_path=input_data_source)

        total_weight = 365 * testObject.stages
        expected_weight = int(total_weight / len(testObject.representative_dates))

        for w in testObject.representative_weights.values():
            self.assertEqual(w, expected_weight)

    def test_no_representative_weights_passed_5_dates(self):
        # Test no representative weights passed in and custom dates passed in, so the function calculates them
        dates = [
            "2020-01-28 00:00",
            "2020-04-23 00:00",
            "2020-07-05 00:00",
            "2020-10-14 00:00",
            "2020-12-14 00:00",
        ]

        testObject = ExpansionPlanningData()
        testObject.load_prescient(
            data_path=input_data_source, representative_dates=dates
        )

        total_weight = 365 * testObject.stages
        expected_weight = int(total_weight / len(testObject.representative_dates))

        for w in testObject.representative_weights.values():
            self.assertEqual(w, expected_weight)

    def test_passed_representative_weights(self):
        # Test representative weights passed in
        weights = {1: 91, 2: 91, 3: 91, 4: 91}
        testObject = ExpansionPlanningData()
        testObject.load_prescient(
            data_path=input_data_source, representative_weights=weights
        )
        self.assertEqual(testObject.representative_weights, weights)

    # FIXME
    def test_in_service_flags_set(self):
        # Test in service flags are being set properly
        testObject = ExpansionPlanningData()
        testObject.load_prescient(data_path=input_data_source)

        assert testObject.md.data["elements"]["generator"]["2-c"]["in_service"] == False
        assert (
            testObject.md.data["elements"]["branch"]["branch2-c"]["in_service"] == False
        )

    def test_missing_simulation_objects_csv(self, tmp_path):
        # test if a data path is missing the simulation objects, it should throw an error
        testObject = ExpansionPlanningData()
        # tmp_path is empty, no simulation_objects.csv
        with self.assertRaises(FileNotFoundError):
            testObject.load_prescient(data_path=tmp_path)

    def test_incorrect_simulation_objects_csv(self, tmp_path):
        # test if simulations_objects is missing key details, it should throw an error
        csv_file = tmp_path / "simulation_objects.csv"
        csv_file.write_text("index,DAY_AHEAD\nWrongKey,24\n")

        testObject = ExpansionPlanningData()
        with self.assertRaises(KeyError):
            testObject.load_prescient(data_path=tmp_path)

    # FIXME
    def test_clone_at_time_keys_called_correctly(self):
        testObject = ExpansionPlanningData()
        testObject.load_prescient(data_path=input_data_source)

        # Patch clone_at_time_keys on the existing model instance
        with unittest.mock.patch.object(
            testObject.md, "clone_at_time_keys", wraps=testObject.md.clone_at_time_keys
        ) as mock_clone:

            # Call load_prescient again to trigger cloning with patched method
            testObject.load_prescient(data_path=input_data_source)

            # Check that clone_at_time_keys was called once per representative date
            expected_calls = len(testObject.representative_dates)
            assert mock_clone.call_count == expected_calls

    # -------------------------------------------------IMPORT_LOAD_SCALING------------------------------------------------------------ #
    # FIXME
    def test_import_load_scaling_normal(self):
        # test successful passthrough of load scaling function
        testObject = ExpansionPlanningData()
        testObject.import_load_scaling(load_scaling_file)

        df = testObject.load_scaling
        assert isinstance(df, pd.DataFrame)
        expected_columns = ["year", "month", "day", "hour"] + [
            str(i) for i in range(1, 9)
        ]
        for col in expected_columns:
            assert col in df.columns
        assert not df.empty

    def test_import_load_scaling_incorrect_num_years(self):
        # Test value error raised if the length of forecast years is incorrect
        testObject = ExpansionPlanningData(stages=3)
        forecast_years = [2025, 2030]

        with self.assertRaises(ValueError):
            testObject.import_load_scaling(load_scaling_file, forecast_years)

    def test_import_load_scaling_incorrect_years_too_early(self):
        # Test value error raised if the forecast years are outside the supported ranges
        testObject = ExpansionPlanningData(stages=3)
        forecast_years = [2019, 2030, 2055]

        with self.assertRaises(ValueError):
            testObject.import_load_scaling(load_scaling_file, forecast_years)

    # -------------------------------------------------IMPORT_OUTAGE_DATA------------------------------------------------------------ #
    def test_import_outage_data(self):
        testObject = ExpansionPlanningData()

        testObject.import_outage_data(outage_data_path)

        df = testObject.bus_hours

        self.assertHasAttr(testObject, "bus_hours")
        self.assertIsInstance(df, pd.DataFrame)
        self.assertIn(["hour", "Bus Number"], df.columns)

    # -------------------------------------------------LOAD_DEFAULT_DATA_SETTINGS------------------------------------------------------------ #
    def test_load_default_data_settings(self):
        testObject = ExpansionPlanningData()
        testObject.load_prescient(data_path=input_data_source)

        testObject.load_default_data_settings()

        # Check generators
        for gen_name, gen in testObject.md.data["elements"]["generator"].items():
            if gen.get("fuel") == "C":
                if gen.get("in_service") is False:
                    self.assertEqual(gen["lifetime"], 1)
                else:
                    self.assertEqual(gen["lifetime"], 2)
            else:
                self.assertEqual(gen["lifetime"], 3)

            # Check other fixed attributes
            self.assertEqual(gen["spinning_reserve_frac"], 0.1)
            self.assertEqual(gen["quickstart_reserve_frac"], 0.1)
            self.assertEqual(gen["capital_multiplier"], 1)
            self.assertEqual(gen["extension_multiplier"], 0)
            self.assertEqual(gen["max_operating_reserve"], 1)
            self.assertEqual(gen["max_spinning_reserve"], 1)
            self.assertEqual(gen["max_quickstart_reserve"], 1)
            self.assertEqual(gen["ramp_up_rate"], 0.1)
            self.assertEqual(gen["ramp_down_rate"], 0.1)
            self.assertEqual(gen["emissions_factor"], 1)
            self.assertEqual(gen["start_fuel"], 1)
            self.assertEqual(gen["investment_cost"], 1)

        # Check branches
        for branch in testObject.md.data["elements"]["branch"].values():
            self.assertEqual(branch["loss_rate"], 0)
            self.assertEqual(branch["distance"], 1)
            self.assertEqual(branch["capital_cost"], 10000000)

        # Check system
        system = testObject.md.data["system"]
        self.assertEqual(system["min_operating_reserve"], 0.1)
        self.assertEqual(system["min_spinning_reserve"], 0.1)

    # -------------------------------------------------LOAD_STORAGE_CSV------------------------------------------------------------ #
    def test_load_storage_csv_success(self):
        testObject = ExpansionPlanningData()
        testObject.load_prescient(data_path=input_data_source)
        testObject.load_storage_csv(storage_file)

        # Check that storage data was loaded into md.data["elements"]["storage"]
        storage = testObject.md.data["elements"].get("storage", None)
        self.assertIsNotNone(storage)
        self.assertIsInstance(storage, dict)

        # Check some expected keys in storage data
        expected_keys = {
            "bus",
            "generator",
            "storage_type",
            "energy_capacity",
            "initial_state_of_charge",
            "investment_cost",
            "investment_cost_kwh",
        }
        for key in expected_keys:
            assert key in storage["100MW_400MWh_1"].keys()

    def test_load_storage_string_path(self):
        testObject = ExpansionPlanningData()
        testObject.load_prescient(data_path=input_data_source)
        testObject.load_storage_csv(str(storage_file))  # should not throw an error

    def test_load_storage_csv_file_not_found(self):
        testObject = ExpansionPlanningData()
        testObject.load_prescient(input_data_source)
        testObject.load_storage_csv(input_data_source)

        # Storage should be set to empty dict
        storage = testObject.md.data["elements"].get("storage", None)
        self.assertIsInstance(storage, dict)
        self.assertEqual(storage, {})

    # -------------------------------------------------TEXAS_CASE_STUDY_UPDATES----------------------------------------------------------- #
    def test_texas_case_study(self):
        testObject = ExpansionPlanningData()
        testObject.load_prescient(data_path=input_data_source)

        # Call the method under test
        testObject.texas_case_study_updates(texas_data_path)

        generator = testObject.md.data["elements"].get("generator", None)
        self.assertIsNotNone(generator)

        expected_columns = [
            "capex1",
            "capex2",
            "capex3",
            "fuel_cost1",
            "fuel_cost2",
            "fuel_cost3",
            "fixed_ops1",
            "fixed_ops2",
            "fixed_ops3",
            "var_ops1",
            "var_ops2",
            "var_ops3",
        ]

        # Check that each expected column is added to each generator
        for gen_name, gen_data in generator.items():
            for col in expected_columns:
                self.assertIn(
                    col, gen_data, f"Column {col} missing in generator {gen_name}"
                )

    def test_texas_case_study_invalid_data_path(self):
        # Test that an error is raised if not a Texas case Study
        testObject = ExpansionPlanningData()

        with self.assertRaises(ValueError, match="not a Texas case study"):
            testObject.texas_case_study_updates(input_data_source)
