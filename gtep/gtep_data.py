# Generation and Transmission Expansion Planning
# IDAES project
# author: Kyle Skolfield
# date: 01/04/2024
# Model available at http://www.optimization-online.org/DB_FILE/2017/08/6162.pdf

from pyomo.environ import *
from prescient.simulator.config import PrescientConfig
from prescient.data.providers import gmlc_data_provider
import datetime
import pandas as pd


class ExpansionPlanningData:
    """Standard data storage class for the IDAES GTEP model."""

    def __init__(self):
        pass

    def load_prescient(self, data_path, options_dict=None):
        """Loads data structured via Prescient data loader.

        :param data_path: Folder containing the data to be loaded
        :param options_dict: Options dictionary to pass to the Prescient data loader, defaults to None
        """
        self.data_type = "prescient"
        options_dict = {
            "data_path": data_path,
            "input_format": "rts-gmlc",
            "start_date": "01-01-2019",
            "num_days": 365,
            "sced_horizon": 1,
            "sced_frequency_minutes": 60,
            "ruc_horizon": 36,
        }

        prescient_options = PrescientConfig()
        prescient_options.set_value(options_dict)
        # Use prescient data provider to load in sequential data for representative periods
        data_list = []

        x = datetime.datetime(2019, 1, 1)
        data_provider = gmlc_data_provider.GmlcDataProvider(options=prescient_options)
        # populate an egret model data with the basic stuff
        self.md = data_provider.get_initial_actuals_model(
            options=prescient_options, num_time_steps=24 * 365, minutes_per_timestep=60
        )
        # fill in renewable actuals and maybe demand idk
        data_provider.populate_with_actuals(
            options=prescient_options,
            num_time_periods=24 * 365,
            time_period_length_minutes=60,
            start_time=x,
            model=self.md,
        )
        # data_provider.populate_initial_state_data(options=prescient_options, model=md)
        self.load_default_data_settings()

        for gen in self.md.data["elements"]["generator"]:
            if "-c" in gen:  # key/Gen UID in csv file; -c = candidate?
                self.md.data["elements"]["generator"][gen]["in_service"] = False

        # JSC addn
        for branch in self.md.data["elements"]["branch"]:
            if "-c" in branch:  # key/Branch UID
                self.md.data["elements"]["branch"][branch]["in_service"] = False

        ## NOTE: Below is only for multiple representative periods and creates a list
        ## of modelData objects, not just a single modelData object
        # Arbitrary time points and lengths picked for representative periods
        # default here allows up to 24 hours for periods
        time_keys = self.md.data["system"]["time_keys"]
        self.representative_dates = [
            "2019-01-28 00:00",
            "2019-04-23 00:00",
            "2019-07-05 00:00",
            "2019-10-14 00:00",
            "2019-08-12 00:00",
        ]

        ## FIXME:
        ## RESIL WEEK ONLY
        ## but we'll want something similar just less insane in the future
        if len(self.representative_dates) == 5:
            self.representative_weights = {1:91, 2:91, 3:91, 4:91, 5:1}
        else:
            self.representative_weights = {1:91, 2:91, 3:91, 4:91}

        for date in self.representative_dates:
            key_idx = time_keys.index(date)
            time_key_set = time_keys[key_idx : key_idx + 24]
            data_list.append(self.md.clone_at_time_keys(time_key_set))

        self.representative_data = data_list

    def import_load_scaling(self, load_file_name):
        adjusted_forecast = pd.read_excel(load_file_name)
        adjusted_forecast_by_period = adjusted_forecast[
            (adjusted_forecast["year"] == 2025)
            | (adjusted_forecast["year"] == 2030)
            | (adjusted_forecast["year"] == 2035)
        ]
        base_zones = [
            "base_economic_coast",
            "base_economic_east",
            "base_economic_fwest",
            "base_economic_ncent",
            "base_economic_north",
            "base_economic_scent",
            "base_economic_south",
            "base_economic_west",
        ]
        scaled_zones = [
            "coast_net",
            "east_net",
            "fwest_net",
            "ncent_net",
            "north_net",
            "scent_net",
            "south_net",
            "west_net",
        ]
        zones = ["coast", "east", "fwest", "ncent", "north", "scent", "south", "west"]
        cap_zones = [zone.upper() for zone in zones]
        for i, zone in enumerate(zones):
            adjusted_forecast_by_period["scaled_" + zone] = (
                adjusted_forecast_by_period[scaled_zones[i]]
                / adjusted_forecast_by_period[base_zones[i]]
            )
        column_list = [
            "year",
            "month",
            "day",
            "hour",
            "scaled_coast",
            "scaled_east",
            "scaled_fwest",
            "scaled_ncent",
            "scaled_north",
            "scaled_scent",
            "scaled_south",
            "scaled_west",
        ]
        load_scaling_df = adjusted_forecast_by_period[column_list]
        scaled_names = ["scaled_" + zone for zone in zones]
        name_conversion_dict = dict(zip(scaled_names, cap_zones))
        load_scaling_df = load_scaling_df.rename(columns=name_conversion_dict)

        self.load_scaling = load_scaling_df

    def load_default_data_settings(self):
        ## TODO: too many of these are hard coded; everything should check if it exists too.
        """Fills in necessary but unspecified data information."""
        for gen in self.md.data["elements"]["generator"]:
            if self.md.data["elements"]["generator"][gen]["fuel"] == "C":
                self.md.data["elements"]["generator"][gen]["lifetime"] = 1
            else:
                self.md.data["elements"]["generator"][gen]["lifetime"] = 3
            self.md.data["elements"]["generator"][gen]["spinning_reserve_frac"] = 0.1
            self.md.data["elements"]["generator"][gen]["quickstart_reserve_frac"] = 0.1
            self.md.data["elements"]["generator"][gen]["capital_multiplier"] = 1
            self.md.data["elements"]["generator"][gen]["extension_multiplier"] = 0
            self.md.data["elements"]["generator"][gen]["max_operating_reserve"] = 1
            self.md.data["elements"]["generator"][gen]["max_spinning_reserve"] = 1
            self.md.data["elements"]["generator"][gen]["max_quickstart_reserve"] = 1
            self.md.data["elements"]["generator"][gen]["ramp_up_rate"] = 0.1
            self.md.data["elements"]["generator"][gen]["ramp_down_rate"] = 0.1
            self.md.data["elements"]["generator"][gen]["emissions_factor"] = 1
            self.md.data["elements"]["generator"][gen]["start_fuel"] = 1
            self.md.data["elements"]["generator"][gen]["investment_cost"] = 235164
        for branch in self.md.data["elements"]["branch"]:
            self.md.data["elements"]["branch"][branch]["loss_rate"] = 0
            self.md.data["elements"]["branch"][branch]["distance"] = 1
        self.md.data["system"]["min_operating_reserve"] = 0.1
        self.md.data["system"]["min_spinning_reserve"] = 0.1

    def texas_case_study_updates(self, data_path):
        generator_update_path = data_path + "/gen.csv"
        generator_df = pd.read_csv(generator_update_path)
        bonus_feature_list = [
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
        for data_point in self.representative_data:
            for col in bonus_feature_list:
                for gen in data_point.data["elements"]["generator"]:
                    if not data_point.data["elements"]["generator"][gen].get(col):
                        data_point.data["elements"]["generator"][gen][col] = float(
                            generator_df[generator_df["GEN UID"] == gen][col]
                        )
