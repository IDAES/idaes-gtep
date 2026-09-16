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

import tempfile
import pandas as pd
from pathlib import Path

import pyomo.common.unittest as unittest

from gtep.gtep_data import ExpansionPlanningData

curr_dir = Path(__file__).resolve().parent
input_data_source = (curr_dir / ".." / ".." / "data" / "5bus").resolve()
data_source_123 = (curr_dir / ".." / ".." / "data" / "123_Bus_Resil_Week").resolve()

load_scaling_file = (
    curr_dir / ".." / ".." / "data" / "Texas_2000" / "ERCOT-Adjusted-Forecast.xlsb"
).resolve()

storage_file = (curr_dir / ".." / ".." / "data" / "9_bus_GTEP_dir").resolve()

texas_data_path = (curr_dir / ".." / ".." / "data" / "123_Bus_Resil_Week").resolve()

outage_data_path = (
    curr_dir / ".." / ".." / "data" / "123_Bus_Resil_Week" / "may_20.csv"
).resolve()

requirements_file = (
    Path(__file__).resolve().parents[2] / "data" / "required_data_files.csv"
)


class TestExpansionPlanningData(unittest.TestCase):

    def test_data_init(self):
        # Test that the ExpansionPlanningData object initializes
        # properly with default values.
        testObject = ExpansionPlanningData()
        self.assertIsInstance(testObject, ExpansionPlanningData)
        self.assertEqual(testObject.stages, 2)
        self.assertEqual(testObject.num_reps, 4)
        self.assertEqual(testObject.num_commit, 24)
        self.assertEqual(testObject.num_dispatch, 1)
        self.assertEqual(testObject.duration_representative_period, 24)

        # Test that the ExpansionPlanningData object initializes
        # properly with input values.
        testObject = ExpansionPlanningData(1, 2, 2, 2, 2, 15)
        self.assertEqual(testObject.stages, 1)
        self.assertEqual(testObject.num_reps, 2)
        self.assertEqual(testObject.num_commit, 2)
        self.assertEqual(testObject.num_dispatch, 2)
        self.assertEqual(testObject.duration_representative_period, 2)

        # Test that the ExpansionPlanningData object initializes
        # properly with partial input values.
        testObject = ExpansionPlanningData()
        self.assertEqual(testObject.stages, 2)
        self.assertEqual(testObject.num_reps, 4)
        self.assertEqual(testObject.num_commit, 24)
        self.assertEqual(testObject.num_dispatch, 1)
        self.assertEqual(testObject.duration_representative_period, 24)

    def test_default_representative_dates(self):
        # Test no representative dates passed in, initializing with
        # defaults.
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
        # Test new representative dates passed in, replacing defaults.
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
        # Test passing invalid/not covered dates.
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

    def test_in_service_flags_set(self):
        # Test in service flags are being set properly.
        testObject = ExpansionPlanningData()
        testObject.load_prescient(data_path=input_data_source)

        self.assertEqual(
            testObject.md.data["elements"]["generator"]["3_CT"]["in_service"], True
        )
        self.assertEqual(
            testObject.md.data["elements"]["branch"]["branch_3_4_1"]["in_service"],
            True,
        )

    def test_missing_simulation_objects_csv(self):
        # Test if a data path is missing the simulation objects, it
        # should throw an error.
        with tempfile.TemporaryDirectory() as tmpdirname:
            testObject = ExpansionPlanningData()
            with self.assertRaises(FileNotFoundError):
                testObject.load_prescient(data_path=tmpdirname)

    def test_load_prescient_creates_expected_representative_data(self):
        representative_dates = [
            "2020-01-28 00:00",
            "2020-04-23 00:00",
        ]

        testObject = ExpansionPlanningData(
            num_reps=2,
            num_commit=24,
            num_dispatch=1,
            duration_representative_period=24,
        )

        testObject.load_prescient(
            data_path=input_data_source,
            representative_dates=representative_dates,
        )

        # Check that one representative ModelData object was created
        # for each requested representative date.
        self.assertEqual(
            len(testObject.representative_data),
            len(representative_dates),
        )

        # Check that the representative dates were stored.
        self.assertEqual(
            len(testObject.representative_dates),
            len(representative_dates),
        )

        for rep_md, rep_date in zip(
            testObject.representative_data,
            representative_dates,
        ):
            # Each representative data object should be a
            # cloned/sliced ModelData object, not the full base model
            # data object.
            self.assertIsNot(rep_md, testObject.md)

            time_keys = rep_md.data["system"]["time_keys"]

            # Each representative data object should contain one day
            # of hourly time keys.
            self.assertEqual(len(time_keys), testObject.num_commit)

            # The first time key should correspond to the requested
            # representative date.
            self.assertEqual(
                time_keys[0],
                pd.to_datetime(rep_date).strftime("%Y-%m-%d %H:%M"),
            )

            # Check that load time series were sliced to the
            # representative period length.
            for load_data in rep_md.data["elements"]["load"].values():
                p_load = load_data.get("p_load")

                if (
                    isinstance(p_load, dict)
                    and p_load.get("data_type") == "time_series"
                ):
                    self.assertEqual(
                        len(p_load["values"]), testObject.duration_representative_period
                    )

            # Check that renewable generator time series were sliced
            # to the representative period length.
            for gen_data in rep_md.data["elements"]["generator"].values():
                p_max = gen_data.get("p_max")

                if isinstance(p_max, dict) and p_max.get("data_type") == "time_series":
                    self.assertEqual(
                        len(p_max["values"]), testObject.duration_representative_period
                    )

    def test_import_load_scaling_normal(self):
        # Test successful passthrough of load scaling function.
        testObject = ExpansionPlanningData()
        testObject.import_load_scaling(load_scaling_file)

        df = testObject.load_scaling
        self.assertIsInstance(df, pd.DataFrame)
        expected_columns = ["year", "month", "day", "hour"] + [
            str(i) for i in range(1, 9)
        ]
        for col in expected_columns:
            self.assertIn(col, df.columns)
        self.assertFalse(df.empty)

    def test_import_load_scaling_incorrect_num_years(self):
        # Test value error raised if the length of forecast years is
        # incorrect.
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

    def test_import_outage_data(self):
        testObject = ExpansionPlanningData()

        testObject.import_outage_data(outage_data_path)

        df = testObject.bus_hours

        self.assertTrue(hasattr(testObject, "bus_hours"))
        self.assertIsInstance(df, pd.DataFrame)
        self.assertIn("hour", df.columns)
        self.assertIn("Bus Number", df.columns)

    def test_load_storage_csv_success(self):
        testObject = ExpansionPlanningData()
        testObject.load_prescient(data_path=input_data_source)
        testObject.load_storage_csv(storage_file)

        # Check that storage data was loaded into
        # md.data["elements"]["storage"]
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
            self.assertIn(key, storage["100MW_400MWh_1"].keys())

    def test_load_storage_string_path(self):
        testObject = ExpansionPlanningData()
        testObject.load_prescient(data_path=input_data_source)
        testObject.load_storage_csv(str(storage_file))  # should not throw an error
        self.assertIn("storage", testObject.md.data["elements"].keys())

        # Check that the storage data is not an empty dict
        self.assertTrue(testObject.md.data["elements"]["storage"])

    def test_load_storage_csv_file_not_found(self):
        testObject = ExpansionPlanningData()
        testObject.load_prescient(data_source_123)
        testObject.load_storage_csv(data_source_123)

        # Storage should be set to empty dict
        storage = testObject.md.data["elements"].get("storage", None)
        self.assertIsInstance(storage, dict)
        self.assertEqual(storage, {})

    def test_texas_case_study(self):
        testObject = ExpansionPlanningData()
        testObject.load_prescient(data_path=texas_data_path)

        # Call the method under test
        testObject.texas_case_study_updates(texas_data_path)

        generator = testObject.md.data["elements"].get("generator", None)
        self.assertIsNotNone(generator)

    def test_texas_case_study_invalid_data_path(self):
        # Test that an error is raised if not a Texas case Study
        testObject = ExpansionPlanningData()

        with self.assertRaises(ValueError):
            testObject.texas_case_study_updates(input_data_source)

    def test_get_start_date_from_simulation_objects(self):
        # Create the data object.
        testObject = ExpansionPlanningData()

        # Read the start date from the simulation_objects.csv in the
        # 5bus case.
        start_date = testObject.get_start_date_from_simulation_objects(
            input_data_source
        )

        # Check that only the date portion is returned, without the
        # hour.
        self.assertEqual(start_date, "01/01/2020")

    def test_assert_start_date_matches_day_ahead_files(self):
        # Create the data object.
        testObject = ExpansionPlanningData()

        # Read the start date from the simulation_objects.csv in the
        # 5bus case.
        start_date = testObject.get_start_date_from_simulation_objects(
            input_data_source
        )

        # Check that the start-date year matches the DAY_AHEAD files.
        testObject.assert_start_date_matches_day_ahead_files(
            input_data_source,
            start_date,
        )

        # If no AssertionError is raised, the year check passed.
        self.assertTrue(True)

    def test_get_required_columns_from_file(self):
        # Create the data object.
        testObject = ExpansionPlanningData()

        # Check that the required-data specification file exists.
        self.assertTrue(requirements_file.exists())

        # Read the required columns for gen.csv.
        required_columns = testObject.get_required_columns_from_file(
            requirements_file,
            "gen.csv",
        )

        # Check that selected required columns from gen.csv are
        # parsed.
        self.assertIn("GEN UID", required_columns)
        self.assertIn("Bus ID", required_columns)
        self.assertIn("Unit Type", required_columns)
        self.assertIn("Fuel", required_columns)
        self.assertIn("PMax MW", required_columns)

    def test_validate_required_columns(self):
        # Create the data object.
        testObject = ExpansionPlanningData()

        # Check that the test input files exist.
        self.assertTrue((input_data_source / "gen.csv").exists())

        # Validate that the existing 5bus gen.csv file has all
        # required columns and required row values.
        testObject.validate_required_columns(
            input_data_source,
            "gen.csv",
            requirements_file=requirements_file,
        )
