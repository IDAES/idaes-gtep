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


# Data Path Fixture (Load Prescient)
@pytest.fixture
def test_data_path(tmp_path):
    # Create a minimal simulation_objects.csv file with expected content
    csv_content = """index,DAY_AHEAD
Periods_per_Step,24
"""
    csv_file = tmp_path / "simulation_objects.csv"
    csv_file.write_text(csv_content)
    load_file = tmp_path / "storage.csv"
    load_file.write_text(
        """name,bus,generator,storage_type,energy_capacity,initial_state_of_charge,end_state_of_charge,minimum_state_of_charge,charge_efficiency,discharge_efficiency,max_discharge_rate,min_discharge_rate,max_charge_rate,min_charge_rate,initial_charge_rate,initial_discharge_rate,charge_cost,discharge_cost,retention_rate_60min,ramp_up_input_60min,ramp_down_input_60min,ramp_up_output_60min,ramp_down_output_60min,in_service,capital_multiplier,extension_multiplier,investment_cost,investment_cost_kwh
100MW_400MWh_1,3,None,'battery',400,320,320,80,1,1,100,10,400,1,0,0,175,0,1,400,400,100,100,FALSE,1,1,1490.89,372.7225"""
    )
    return tmp_path


# Load Data File Fixture (Import Load)
@pytest.fixture
def actual_load_path():
    test_dir = Path(__file__).parent
    # Navigate up to 'gtep' directory, then into 'data'
    data_dir = test_dir.parent.parent / "data"
    excel_path = data_dir / "Texas_2000" / "ERCOT-Adjusted-Forecast.xlsb"
    return excel_path


# Outage Data File Fixtures (Import Outages)
@pytest.fixture
def outage_data(tmp_path):
    outage_content = """fips_code,lim_timestamp,case_4b_prob,date
48001,5/20/25 0:00,0,20-May
48025,5/20/25 3:00,0,20-May
48025,5/20/25 4:00,0,20-May
48025,5/20/25 5:00,0,20-May
48025,5/20/25 6:00,0,20-May
48025,5/20/25 7:00,0,20-May
48025,5/20/25 8:00,0,20-May
48025,5/20/25 9:00,0,20-May
48025,5/20/25 10:00,0,20-May
48027,5/20/25 20:00,0,20-May
48163,5/20/25 20:00,0.006944444,20-May"""
    outage_csv = tmp_path / "outage.csv"
    outage_csv.write_text(outage_content)

    return outage_csv


@pytest.fixture
def county_fips_csv(tmp_path):
    dir_path = tmp_path / "gtep" / "data" / "123_Bus_Resil_Week"
    dir_path.mkdir(parents=True, exist_ok=True)
    content = """county_number,County,FIPS,,
1,Anderson,48001,,
13,Bee,48025,,
27,Burnet,48053,,
82,Frio,48163,,
"""
    file = dir_path / "county_fips_match.csv"
    file.write_text(content)
    return file


@pytest.fixture
def bus_data_csv(tmp_path):
    dir_path = tmp_path / "gtep" / "data" / "123_Bus_Resil_Week"
    dir_path.mkdir(parents=True, exist_ok=True)
    content = """Bus Number,County
1001,Anderson
1002,Bee
1003,Burnet
1004,Frio
"""
    file = dir_path / "Bus_data_gen_weights_mappings.csv"
    file.write_text(content)
    return file


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


@pytest.fixture
def mock_gmlc_provider():
    with unittest.mock.patch(
        "prescient.data.providers.gmlc_data_provider.GmlcDataProvider"
    ) as MockProvider:
        instance = MockProvider.return_value
        instance._start_time = "2020-01-01 00:00"
        # Mock get_initial_actuals_model to return a model with .data attribute
        model_mock = unittest.mock.Mock()
        model_mock.data = {
            "elements": {
                "generator": {
                    "gen1": {"fuel": "C", "in_service": True},
                    "gen2-c": {"fuel": "None", "in_service": True},
                },
                "branch": {"branch1": {}, "branch2-c": {}},
                "storage": {"stor1": {}},
            },
            "system": {
                "time_keys": [
                    "2020-01-28 00:00",
                    "2020-01-28 01:00",
                    "2020-01-28 02:00",
                    "2020-04-23 00:00",
                    "2020-07-05 00:00",
                    "2020-10-14 00:00",
                    "2020-12-14 00:00",
                ]
            },
        }
        model_mock.clone_at_time_keys.side_effect = lambda keys: f"clone_{keys}"
        instance.get_initial_actuals_model.return_value = model_mock
        instance.populate_with_actuals.return_value = None
        yield instance


@pytest.mark.usefixtures("mock_prescient_config", "mock_gmlc_provider")
class TestExpansionPlanningData:

    def test_data_init(self):
        # Test that the ExpansionPlanningData object initializes properly with default values
        testObject = ExpansionPlanningData()
        assert isinstance(testObject, ExpansionPlanningData)
        assert testObject.stages == 2
        assert testObject.num_reps == 4
        assert testObject.len_reps == 1
        assert testObject.num_commit == 24
        assert testObject.num_dispatch == 1
        assert testObject.duration_dispatch == 60

        # Test that the ExpansionPlanningData object initializes properly with input values
        testObject = ExpansionPlanningData(1, 2, 2, 2, 2, 15)
        assert isinstance(testObject, ExpansionPlanningData)
        assert testObject.stages == 1
        assert testObject.num_reps == 2
        assert testObject.len_reps == 2
        assert testObject.num_commit == 2
        assert testObject.num_dispatch == 2
        assert testObject.duration_dispatch == 15

        # Test that the ExpansionPlanningData object initializes properly with partial input values
        testObject = ExpansionPlanningData(duration_dispatch=15)
        assert isinstance(testObject, ExpansionPlanningData)
        assert testObject.stages == 2
        assert testObject.num_reps == 4
        assert testObject.len_reps == 1
        assert testObject.num_commit == 24
        assert testObject.num_dispatch == 1
        assert testObject.duration_dispatch == 15

    # -------------------------------------------------LOAD_PRESCIENT------------------------------------------------------------ #
    def test_options_dict(
        self, mock_prescient_config, mock_gmlc_provider, test_data_path
    ):
        # Test passing in an options dictionary
        testObject = ExpansionPlanningData()
        # set options that are not the defaults
        options = {
            "num_days": 100,
            "ruc_horizon": 12,
        }

        testObject.load_prescient(data_path=str(test_data_path), options_dict=options)
        mock_prescient_config.set_value.assert_called_once()
        passed_dict = mock_prescient_config.set_value.call_args[0][0]
        assert passed_dict["data_path"] == str(test_data_path)
        assert passed_dict["num_days"] == 100
        assert passed_dict["ruc_horizon"] == 12

    def test_no_options_dict(
        self, mock_prescient_config, mock_gmlc_provider, test_data_path
    ):

        # Test not passing in an options dictionary
        testObject = ExpansionPlanningData()
        testObject.load_prescient(data_path=str(test_data_path))
        mock_prescient_config.set_value.assert_called_once()
        passed_dict = mock_prescient_config.set_value.call_args[0][0]
        # Default options should be set
        assert passed_dict["data_path"] == str(test_data_path)
        assert passed_dict["num_days"] == 365
        assert passed_dict["ruc_horizon"] == 36

    def test_default_representative_dates(
        self, mock_prescient_config, mock_gmlc_provider, test_data_path
    ):
        # Test no representative dates passed in, initializing with defaults
        testObject = ExpansionPlanningData()
        testObject.load_prescient(data_path=str(test_data_path))
        # default dates:
        expected_dates = [
            "2020-01-28 00:00",
            "2020-04-23 00:00",
            "2020-07-05 00:00",
            "2020-10-14 00:00",
        ]
        assert testObject.representative_dates == expected_dates

    def test_passed_representative_dates(
        self, mock_prescient_config, mock_gmlc_provider, test_data_path
    ):
        # Test new representative dates passed in, replacing defaults
        testObject = ExpansionPlanningData()
        testObject.load_prescient(
            data_path=str(test_data_path),
            representative_dates=[
                "2020-01-28 00:00",
                "2020-04-23 00:00",
                "2020-07-05 00:00",
                "2020-12-14 00:00",
            ],
        )

        expected_dates = [
            "2020-01-28 00:00",
            "2020-04-23 00:00",
            "2020-07-05 00:00",
            "2020-12-14 00:00",
        ]
        assert testObject.representative_dates == expected_dates

    def test_representative_date_not_in_time_keys(
        self, mock_prescient_config, mock_gmlc_provider, test_data_path
    ):
        # Test passing invalid/not covered dates
        testObject = ExpansionPlanningData()
        bad_dates = [
            "2020-01-28 00:00",
            "2020-04-23 00:00",
            "2020-07-05 00:00",
            "2099-01-01 00:00",  # Not in time_keys
        ]
        with pytest.raises(ValueError):
            testObject.load_prescient(
                data_path=str(test_data_path), representative_dates=bad_dates
            )

    def test_empty_representative_dates(
        self, mock_prescient_config, mock_gmlc_provider, test_data_path
    ):
        # Test an empty list passed to overwrite defaults without setting new dates
        testObject = ExpansionPlanningData()
        with pytest.raises(ZeroDivisionError):
            testObject.load_prescient(
                data_path=str(test_data_path), representative_dates=[]
            )

    def test_invalid_representative_dates_type(
        self, mock_prescient_config, mock_gmlc_provider, test_data_path
    ):
        # Test setting None for the dates, overwriting defaults without adding new dates
        testObject = ExpansionPlanningData()
        with pytest.raises(TypeError):
            testObject.load_prescient(
                data_path=str(test_data_path), representative_dates=None
            )

    def test_default_representative_weights(
        self, mock_prescient_config, mock_gmlc_provider, test_data_path
    ):
        # Test no representative weights passed in, so the function calculates them
        testObject = ExpansionPlanningData()
        testObject.load_prescient(data_path=str(test_data_path))

        total_weight = mock_prescient_config.num_days * testObject.stages
        assert total_weight == (365 * 2)
        expected_weight = int(total_weight / len(testObject.representative_dates))
        assert expected_weight == int((365 * 2) / 4)

        for w in testObject.representative_weights.values():
            assert w == expected_weight

    def test_no_representative_weights_passed_5_dates(
        self, mock_prescient_config, mock_gmlc_provider, test_data_path
    ):
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
            data_path=str(test_data_path), representative_dates=dates
        )

        total_weight = mock_prescient_config.num_days * testObject.stages
        assert total_weight == (365 * 2)
        expected_weight = int(total_weight / len(testObject.representative_dates))
        assert expected_weight == int((365 * 2) / 5)

        for w in testObject.representative_weights.values():
            assert w == expected_weight

    def test_passed_representative_weights(
        self, mock_prescient_config, mock_gmlc_provider, test_data_path
    ):
        # Test representative weights passed in
        weights = {1: 91, 2: 91, 3: 91, 4: 91}
        testObject = ExpansionPlanningData()
        testObject.load_prescient(
            data_path=str(test_data_path), representative_weights=weights
        )
        assert testObject.representative_weights == weights

    def test_key_function_calls(
        self, mock_prescient_config, mock_gmlc_provider, test_data_path
    ):
        # Test that each function is called within the load_prescient function
        testObject = ExpansionPlanningData()

        # Patch internal methods to track calls
        with (
            unittest.mock.patch.object(
                testObject, "load_default_data_settings"
            ) as mock_load_defaults,
            unittest.mock.patch.object(
                testObject, "load_storage_csv"
            ) as mock_load_storage,
        ):

            testObject.load_prescient(data_path=str(test_data_path))

            # External calls
            mock_prescient_config.set_value.assert_called_once()
            mock_gmlc_provider.get_initial_actuals_model.assert_called_once()
            mock_gmlc_provider.populate_with_actuals.assert_called_once()

            # Internal calls
            mock_load_defaults.assert_called_once()
            mock_load_storage.assert_called_once_with(str(test_data_path))

    def test_total_num_steps_calculation(
        self, mock_prescient_config, mock_gmlc_provider, test_data_path
    ):
        # Check that the total number of time steps is calculated based in the number of days in the config
        mock_prescient_config.num_days = 365

        testObject = ExpansionPlanningData()
        testObject.load_prescient(data_path=str(test_data_path))

        # period_per_step is 24
        expected_total_steps = 365 * 24

        mock_gmlc_provider.get_initial_actuals_model.assert_called_once()
        call_args = mock_gmlc_provider.get_initial_actuals_model.call_args[1]
        assert call_args["num_time_steps"] == expected_total_steps

    def test_total_num_steps_calculation_2_years(
        self, mock_prescient_config, mock_gmlc_provider, test_data_path
    ):
        # Check that the total number of time steps is calculated based in the number of days in the config
        mock_prescient_config.num_days = 730

        testObject = ExpansionPlanningData()
        testObject.load_prescient(data_path=str(test_data_path))

        # period_per_step is 24
        expected_total_steps = 730 * 24

        mock_gmlc_provider.get_initial_actuals_model.assert_called_once()
        call_args = mock_gmlc_provider.get_initial_actuals_model.call_args[1]
        assert call_args["num_time_steps"] == expected_total_steps

    def test_in_service_flags_set(
        self, mock_prescient_config, mock_gmlc_provider, test_data_path
    ):
        # Test in service flags are being set properly
        testObject = ExpansionPlanningData()
        testObject.load_prescient(data_path=str(test_data_path))

        assert (
            testObject.md.data["elements"]["generator"]["gen2-c"]["in_service"] == False
        )
        assert (
            testObject.md.data["elements"]["branch"]["branch2-c"]["in_service"] == False
        )

    def test_missing_simulation_objects_csv(
        self, mock_prescient_config, mock_gmlc_provider, tmp_path
    ):
        # test if a data path is missing the simulation objects, it should throw an error
        testObject = ExpansionPlanningData()
        # tmp_path is empty, no simulation_objects.csv
        with pytest.raises(FileNotFoundError):
            testObject.load_prescient(data_path=str(tmp_path))

    def test_incorrect_simulation_objects_csv(
        self, mock_prescient_config, mock_gmlc_provider, tmp_path
    ):
        # test if simulations_objects is missing key details, it should throw an error
        csv_file = tmp_path / "simulation_objects.csv"
        csv_file.write_text("index,DAY_AHEAD\nWrongKey,24\n")

        testObject = ExpansionPlanningData()
        with pytest.raises(KeyError):
            testObject.load_prescient(data_path=str(tmp_path))

    def test_clone_at_time_keys_called_correctly(
        self, mock_prescient_config, mock_gmlc_provider, test_data_path
    ):
        testObject = ExpansionPlanningData()
        testObject.load_prescient(data_path=str(test_data_path))

        # Patch clone_at_time_keys on the existing model instance
        with unittest.mock.patch.object(
            testObject.md, "clone_at_time_keys", wraps=testObject.md.clone_at_time_keys
        ) as mock_clone:

            # Call load_prescient again to trigger cloning with patched method
            testObject.load_prescient(data_path=str(test_data_path))

            # Check that clone_at_time_keys was called once per representative date
            expected_calls = len(testObject.representative_dates)
            assert mock_clone.call_count == expected_calls

    # -------------------------------------------------IMPORT_LOAD_SCALING------------------------------------------------------------ #
    def test_import_load_scaling_normal_with_actual_file(self, actual_load_path):
        # test successful passthrough of load scaling function
        testObject = ExpansionPlanningData()
        testObject.import_load_scaling(actual_load_path)

        df = testObject.load_scaling
        assert isinstance(df, pd.DataFrame)
        expected_columns = ["year", "month", "day", "hour"] + [
            str(i) for i in range(1, 9)
        ]
        for col in expected_columns:
            assert col in df.columns
        assert not df.empty

    def test_import_load_scaling_incorrect_num_years(self, actual_load_path):
        # Test value error raised if the length of forecast years is incorrect
        testObject = ExpansionPlanningData(stages=3)
        forecast_years = [2025, 2030]

        with pytest.raises(ValueError):
            testObject.import_load_scaling(actual_load_path, forecast_years)

    def test_import_load_scaling_incorrect_years_too_early(self, actual_load_path):
        # Test value error raised if the forecast years are outside the supported ranges
        testObject = ExpansionPlanningData(stages=3)
        forecast_years = [2019, 2030, 2055]

        with pytest.raises(ValueError):
            testObject.import_load_scaling(actual_load_path, forecast_years)

    # -------------------------------------------------IMPORT_OUTAGE_DATA------------------------------------------------------------ #
    def test_import_outage_data(
        self, outage_data, county_fips_csv, bus_data_csv, monkeypatch
    ):
        testObject = ExpansionPlanningData()

        # Load DataFrames from temp CSV files for patching
        df_outage = pd.read_csv(outage_data)
        df_county_fips = pd.read_csv(county_fips_csv)
        df_bus_data = pd.read_csv(bus_data_csv)

        # Patch pd.read_csv to return appropriate DataFrame based on input path
        def mock_read_csv(filepath, *args, **kwargs):
            if filepath == str(outage_data):
                return df_outage
            elif filepath.endswith("county_fips_match.csv"):
                return df_county_fips
            elif filepath.endswith("Bus_data_gen_weights_mappings.csv"):
                return df_bus_data
            else:
                raise FileNotFoundError(f"Unexpected file path: {filepath}")

        monkeypatch.setattr(pd, "read_csv", mock_read_csv)

        # Patch to_csv to avoid actual file writes during test
        monkeypatch.setattr(pd.DataFrame, "to_csv", lambda self, *args, **kwargs: None)

        testObject.import_outage_data(str(outage_data))

        assert hasattr(testObject, "bus_hours")
        df = testObject.bus_hours
        assert isinstance(df, pd.DataFrame)
        assert set(df.columns) == {"hour", "Bus Number"}
        assert pd.api.types.is_integer_dtype(df["hour"])
        assert pd.api.types.is_integer_dtype(df["Bus Number"])
        assert (df["hour"] == 20).any()
        assert (df["Bus Number"] == 1004).any()

    # -------------------------------------------------LOAD_DEFAULT_DATA_SETTINGS------------------------------------------------------------ #
    def test_load_default_data_settings(
        self, mock_prescient_config, mock_gmlc_provider, test_data_path
    ):
        testObject = ExpansionPlanningData()
        testObject.load_prescient(data_path=str(test_data_path))

        testObject.load_default_data_settings()

        # Check generators
        for gen_name, gen in testObject.md.data["elements"]["generator"].items():
            if gen.get("fuel") == "C":
                if gen.get("in_service") is False:
                    assert gen["lifetime"] == 1
                else:
                    assert gen["lifetime"] == 2
            else:
                assert gen["lifetime"] == 3

            # Check other fixed attributes
            assert gen["spinning_reserve_frac"] == 0.1
            assert gen["quickstart_reserve_frac"] == 0.1
            assert gen["capital_multiplier"] == 1
            assert gen["extension_multiplier"] == 0
            assert gen["max_operating_reserve"] == 1
            assert gen["max_spinning_reserve"] == 1
            assert gen["max_quickstart_reserve"] == 1
            assert gen["ramp_up_rate"] == 0.1
            assert gen["ramp_down_rate"] == 0.1
            assert gen["emissions_factor"] == 1
            assert gen["start_fuel"] == 1
            assert gen["investment_cost"] == 1

        # Check branches
        for branch in testObject.md.data["elements"]["branch"].values():
            assert branch["loss_rate"] == 0
            assert branch["distance"] == 1
            assert branch["capital_cost"] == 10000000

        # Check system
        system = testObject.md.data["system"]
        assert system["min_operating_reserve"] == 0.1
        assert system["min_spinning_reserve"] == 0.1

    # -------------------------------------------------LOAD_STORAGE_CSV------------------------------------------------------------ #
    def test_load_storage_csv_success(
        self, test_data_path, mock_gmlc_provider, mock_prescient_config
    ):
        testObject = ExpansionPlanningData()
        # Setup md with empty data structure to avoid errors
        testObject.md = mock_gmlc_provider.get_initial_actuals_model.return_value
        testObject.load_storage_csv(str(test_data_path))

        # Check that storage data was loaded into md.data["elements"]["storage"]
        storage = testObject.md.data["elements"].get("storage", None)
        assert storage is not None
        assert isinstance(storage, dict)

        # Check that storage keys include the name from your CSV fixture
        assert "100MW_400MWh_1" in storage

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
            assert key in storage["100MW_400MWh_1"]

    def test_load_storage_csv_file_not_found(
        self, tmp_path, mock_gmlc_provider, mock_prescient_config
    ):
        testObject = ExpansionPlanningData()
        testObject.md = mock_gmlc_provider.get_initial_actuals_model.return_value

        # Use a directory without storage.csv
        missing_path = tmp_path

        with unittest.mock.patch("builtins.print") as mock_print:
            testObject.load_storage_csv(str(missing_path))

            # Should print warning about missing file
            mock_print.assert_called_once()
            args = mock_print.call_args[0][0]
            assert "does not exist" in args

        # Storage should be set to empty dict
        storage = testObject.md.data["elements"].get("storage", None)
        assert storage == {}
