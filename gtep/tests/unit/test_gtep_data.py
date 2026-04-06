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


@pytest.fixture
def simulation_objects_csv(tmp_path):
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
                "generator": {"gen1": {}, "gen2-c": {}},
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
        self, mock_prescient_config, mock_gmlc_provider, simulation_objects_csv
    ):
        # Test passing in an options dictionary
        testObject = ExpansionPlanningData()
        # set options that are not the defaults
        options = {
            "num_days": 100,
            "ruc_horizon": 12,
        }

        testObject.load_prescient(
            data_path=str(simulation_objects_csv), options_dict=options
        )
        mock_prescient_config.set_value.assert_called_once()
        passed_dict = mock_prescient_config.set_value.call_args[0][0]
        assert passed_dict["data_path"] == str(simulation_objects_csv)
        assert passed_dict["num_days"] == 100
        assert passed_dict["ruc_horizon"] == 12

    def test_no_options_dict(
        self, mock_prescient_config, mock_gmlc_provider, simulation_objects_csv
    ):

        # Test not passing in an options dictionary
        testObject = ExpansionPlanningData()
        testObject.load_prescient(data_path=str(simulation_objects_csv))
        mock_prescient_config.set_value.assert_called_once()
        passed_dict = mock_prescient_config.set_value.call_args[0][0]
        # Default options should be set
        assert passed_dict["data_path"] == str(simulation_objects_csv)
        assert passed_dict["num_days"] == 365
        assert passed_dict["ruc_horizon"] == 36

    def test_default_representative_dates(
        self, mock_prescient_config, mock_gmlc_provider, simulation_objects_csv
    ):
        # Test no representative dates passed in, initializing with defaults
        testObject = ExpansionPlanningData()
        testObject.load_prescient(data_path=str(simulation_objects_csv))
        # default dates:
        expected_dates = [
            "2020-01-28 00:00",
            "2020-04-23 00:00",
            "2020-07-05 00:00",
            "2020-10-14 00:00",
        ]
        assert testObject.representative_dates == expected_dates

    def test_passed_representative_dates(
        self, mock_prescient_config, mock_gmlc_provider, simulation_objects_csv
    ):
        # Test new representative dates passed in, replacing defaults
        testObject = ExpansionPlanningData()
        testObject.load_prescient(
            data_path=str(simulation_objects_csv),
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
        self, mock_prescient_config, mock_gmlc_provider, simulation_objects_csv
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
                data_path=str(simulation_objects_csv), representative_dates=bad_dates
            )

    def test_empty_representative_dates(
        self, mock_prescient_config, mock_gmlc_provider, simulation_objects_csv
    ):
        # Test an empty list passed to overwrite defaults without setting new dates
        testObject = ExpansionPlanningData()
        with pytest.raises(ZeroDivisionError):
            testObject.load_prescient(
                data_path=str(simulation_objects_csv), representative_dates=[]
            )

    def test_invalid_representative_dates_type(
        self, mock_prescient_config, mock_gmlc_provider, simulation_objects_csv
    ):
        # Test setting None for the dates, overwriting defaults without adding new dates
        testObject = ExpansionPlanningData()
        with pytest.raises(TypeError):
            testObject.load_prescient(
                data_path=str(simulation_objects_csv), representative_dates=None
            )

    def test_default_representative_weights(
        self, mock_prescient_config, mock_gmlc_provider, simulation_objects_csv
    ):
        # Test no representative weights passed in, so the function calculates them
        testObject = ExpansionPlanningData()
        testObject.load_prescient(data_path=str(simulation_objects_csv))

        total_weight = mock_prescient_config.num_days * testObject.stages
        assert total_weight == (365 * 2)
        expected_weight = int(total_weight / len(testObject.representative_dates))
        assert expected_weight == int((365 * 2) / 4)

        for w in testObject.representative_weights.values():
            assert w == expected_weight

    def test_no_representative_weights_passed_5_dates(
        self, mock_prescient_config, mock_gmlc_provider, simulation_objects_csv
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
            data_path=str(simulation_objects_csv), representative_dates=dates
        )

        total_weight = mock_prescient_config.num_days * testObject.stages
        assert total_weight == (365 * 2)
        expected_weight = int(total_weight / len(testObject.representative_dates))
        assert expected_weight == int((365 * 2) / 5)

        for w in testObject.representative_weights.values():
            assert w == expected_weight

    def test_passed_representative_weights(
        self, mock_prescient_config, mock_gmlc_provider, simulation_objects_csv
    ):
        # Test representative weights passed in
        weights = {1: 91, 2: 91, 3: 91, 4: 91}
        testObject = ExpansionPlanningData()
        testObject.load_prescient(
            data_path=str(simulation_objects_csv), representative_weights=weights
        )
        assert testObject.representative_weights == weights

    def test_key_function_calls(
        self, mock_prescient_config, mock_gmlc_provider, simulation_objects_csv
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

            testObject.load_prescient(data_path=str(simulation_objects_csv))

            # External calls
            mock_prescient_config.set_value.assert_called_once()
            mock_gmlc_provider.get_initial_actuals_model.assert_called_once()
            mock_gmlc_provider.populate_with_actuals.assert_called_once()

            # Internal calls
            mock_load_defaults.assert_called_once()
            mock_load_storage.assert_called_once_with(str(simulation_objects_csv))

    def test_total_num_steps_calculation(
        self, mock_prescient_config, mock_gmlc_provider, simulation_objects_csv
    ):
        # Check that the total number of time steps is calculated based in the number of days in the config
        mock_prescient_config.num_days = 365

        testObject = ExpansionPlanningData()
        testObject.load_prescient(data_path=str(simulation_objects_csv))

        # period_per_step is 24
        expected_total_steps = 365 * 24

        mock_gmlc_provider.get_initial_actuals_model.assert_called_once()
        call_args = mock_gmlc_provider.get_initial_actuals_model.call_args[1]
        assert call_args["num_time_steps"] == expected_total_steps

    def test_total_num_steps_calculation_2_years(
        self, mock_prescient_config, mock_gmlc_provider, simulation_objects_csv
    ):
        # Check that the total number of time steps is calculated based in the number of days in the config
        mock_prescient_config.num_days = 730

        testObject = ExpansionPlanningData()
        testObject.load_prescient(data_path=str(simulation_objects_csv))

        # period_per_step is 24
        expected_total_steps = 730 * 24

        mock_gmlc_provider.get_initial_actuals_model.assert_called_once()
        call_args = mock_gmlc_provider.get_initial_actuals_model.call_args[1]
        assert call_args["num_time_steps"] == expected_total_steps

    def test_in_service_flags_set(
        self, mock_prescient_config, mock_gmlc_provider, simulation_objects_csv
    ):
        # Test in service flags are being set properly
        testObject = ExpansionPlanningData()
        testObject.load_prescient(data_path=str(simulation_objects_csv))

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
        self, mock_prescient_config, mock_gmlc_provider, simulation_objects_csv
    ):
        testObject = ExpansionPlanningData()
        testObject.load_prescient(data_path=str(simulation_objects_csv))

        # Patch clone_at_time_keys on the existing model instance
        with unittest.mock.patch.object(
            testObject.md, "clone_at_time_keys", wraps=testObject.md.clone_at_time_keys
        ) as mock_clone:

            # Call load_prescient again to trigger cloning with patched method
            testObject.load_prescient(data_path=str(simulation_objects_csv))

            # Check that clone_at_time_keys was called once per representative date
            expected_calls = len(testObject.representative_dates)
            assert mock_clone.call_count == expected_calls


# -------------------------------------------------IMPORT_LOAD_SCALING------------------------------------------------------------ #
