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

# Generation and Transmission Expansion Planning
# IDAES project
# author: Kyle Skolfield
# date: 01/04/2024
# Model available at http://www.optimization-online.org/DB_FILE/2017/08/6162.pdf

from pyomo.environ import *
from prescient.simulator.config import PrescientConfig
from prescient.data.providers import gmlc_data_provider
import pandas as pd
import os
from pathlib import Path


class ExpansionPlanningData:
    """Standard data storage class for the IDAES GTEP model."""

    def __init__(
        self,
        stages=2,
        num_reps=4,
        len_reps=1,
        num_commit=24,
        num_dispatch=1,
        duration_dispatch=60,
    ):
        """Initialize generation & expansion planning data object.

        :param stages: integer number of investment periods
        :param num_reps: integer number of representative periods per investment period
        :param len_reps: (for now integer) length of each representative period (in hours)
        :param num_commit: integer number of commitment periods per representative period
        :param num_dispatch: integer number of dispatch periods per commitment period
        :param duration_dispatch: (for now integer) duration of each dispatch period (in minutes)
        """
        self.stages = stages
        self.num_reps = num_reps
        self.len_reps = len_reps
        self.num_commit = num_commit
        self.num_dispatch = num_dispatch
        self.duration_dispatch = duration_dispatch

    def load_prescient(
        self,
        data_path,
        representative_dates=None,
        representative_weights={},
        options_dict=None,
    ):
        """Loads data structured via Prescient data loader.

        :param data_path: Folder containing the data to be loaded
        :param representative_dates: List of time points to include. Note: Change the last date for whatever extreme day is needed based on the given run(s)
        :param representative_weights: dictionary of weights for each representative date, defaults to empty Dict
        :param options_dict: Options dictionary to pass to the Prescient data loader, defaults to None

        """
        self.data_type = "prescient"
        # create prescient config object with defaults
        prescient_options = PrescientConfig()

        # work around for prescient throwing an error with Path objects
        if isinstance(data_path, Path):
            data_path = str(data_path)

        if options_dict is None:
            # set basic configurations that do not match prescient defaults
            options_dict = {
                "data_path": data_path,
                "num_days": 365,
                "ruc_horizon": 36,
            }

        else:
            # ensure data path is included in options dictionary
            options_dict["data_path"] = data_path

        # update configuration values based on options dictionary
        prescient_options.set_value(options_dict)

        # Use prescient data provider to load in sequential data for representative periods
        data_list = []

        data_provider = gmlc_data_provider.GmlcDataProvider(options=prescient_options)

        # grab details from simulation objects file (data provider above throws error if no simulation_objects.csv exists)
        metadata_path = os.path.join(data_path, "simulation_objects.csv")
        metadata_df = pd.read_csv(metadata_path, index_col=0)

        # save to variable for easy calling
        sced_freq_min = prescient_options.sced_frequency_minutes

        # This step is grabbing DAY_AHEAD information for now
        # (in the future we may want to update to grab the "REAL_TIME" data if the data has reliable data since the actuals model is looking for real time data info)
        period_per_step = int(metadata_df.loc["Periods_per_Step"]["DAY_AHEAD"])
        total_num_steps = prescient_options.num_days * period_per_step

        # populate an egret model data with the basic stuff
        self.md = data_provider.get_initial_actuals_model(
            options=prescient_options,
            num_time_steps=total_num_steps,
            minutes_per_timestep=sced_freq_min,
        )

        # fill in renewable actuals
        data_provider.populate_with_actuals(
            options=prescient_options,
            num_time_periods=total_num_steps,
            time_period_length_minutes=sced_freq_min,
            start_time=data_provider._start_time,
            model=self.md,
        )

        # data_provider.populate_initial_state_data(options=prescient_options, model=md)
        self.load_default_data_settings()

        self.load_storage_csv(data_path)

        for gen in self.md.data["elements"]["generator"]:
            if "-c" in gen:  # key/Gen UID in csv file; -c = candidate?
                self.md.data["elements"]["generator"][gen]["in_service"] = False

        # JSC addn
        for branch in self.md.data["elements"]["branch"]:
            if "-c" in branch:  # key/Branch UID
                self.md.data["elements"]["branch"][branch]["in_service"] = False

        for stor in self.md.data["elements"]["storage"]:
            if "-c" in stor:  # key/Branch UID
                self.md.data["elements"]["storage"][stor]["in_service"] = False

        ## NOTE: Below is only for multiple representative periods and creates a list
        ## of modelData objects, not just a single modelData object
        # Arbitrary time points and lengths picked for representative periods
        # default here allows up to 24 hours for periods
        if self.representative_dates is None:
            representative_dates = [
                "2020-01-28 00:00",
                "2020-04-23 00:00",
                "2020-07-05 00:00",
                "2020-10-14 00:00",  ## Change the last date for whatever extreme day is needed based on the given run(s)
            ]
        self.representative_dates = representative_dates

        if not representative_weights:
            # set the weight for each day to the total weight divided by number of days
            total_weight = prescient_options.num_days * self.stages
            weight_per_date = int(total_weight / (len(representative_dates)))
            self.representative_weights = {
                key: weight_per_date
                for date, key in enumerate(self.representative_dates)
            }

        time_keys = self.md.data["system"]["time_keys"]

        for date in self.representative_dates:
            key_idx = time_keys.index(date)
            time_key_set = time_keys[key_idx : key_idx + period_per_step]
            data_list.append(self.md.clone_at_time_keys(time_key_set))

        self.representative_data = data_list

    def import_load_scaling(self, load_file_name, forecast_years=None):
        """Imports load scaling data for forecast years.

        :param load_file_name: filepath for adjusted forecast excel file
        :param forecast_years: list of years to forecast, defaults to [2025, 2030, 2035]

        """
        adjusted_forecast = pd.read_excel(load_file_name)

        if forecast_years is None:
            forecast_years = [2025, 2030, 2035]

        # check years are valid
        if len(forecast_years) < self.stages:
            raise ValueError(
                "Not enough forecast years for the number of stages of investment"
            )
        elif any(year < 2020 or year > 2050 for year in forecast_years):
            raise ValueError(
                "The list of years includes a year before 2020 or after 2050."
            )

        adjusted_forecast_by_period = adjusted_forecast[
            adjusted_forecast["year"].isin(forecast_years)
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
        # zones = ["coast", "east", "fwest", "ncent", "north", "scent", "south", "west"]
        # cap_zones = [zone.upper() for zone in zones]
        zones = ["1", "2", "3", "4", "5", "6", "7", "8"]
        cap_zones = ["1", "2", "3", "4", "5", "6", "7", "8"]
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
            "scaled_1",
            "scaled_2",
            "scaled_3",
            "scaled_4",
            "scaled_5",
            "scaled_6",
            "scaled_7",
            "scaled_8",
        ]
        load_scaling_df = adjusted_forecast_by_period[column_list]
        scaled_names = ["scaled_" + zone for zone in zones]
        name_conversion_dict = dict(zip(scaled_names, cap_zones))
        load_scaling_df = load_scaling_df.rename(columns=name_conversion_dict)
        self.load_scaling = load_scaling_df

    def import_outage_data(self, load_file_name):
        """Imports outage data.

        :param load_file_name: filepath for adjusted forecast excel file

        """
        outage_list = pd.read_csv(load_file_name)
        percentile_threshold = 0.9
        threshold_value = outage_list["case_4b_prob"].quantile(percentile_threshold)
        filtered_outages = outage_list[outage_list["case_4b_prob"] >= threshold_value]

        filtered_outages["hour"] = filtered_outages["lim_timestamp"].str.extract(
            r" (\d+):"
        )
        filtered_outages = filtered_outages[["fips_code", "hour"]]
        county_to_fips = pd.read_csv(
            "./gtep/data/123_Bus_Resil_Week/county_fips_match.csv"
        )
        bus_to_county = pd.read_csv(
            "./gtep/data/123_Bus_Resil_Week/Bus_data_gen_weights_mappings.csv"
        )
        county_to_fips = county_to_fips[["County", "FIPS"]]
        bus_to_county = bus_to_county[["Bus Number", "County"]]
        bus_to_county = bus_to_county.merge(county_to_fips, how="inner", on="County")
        bus_hours = pd.merge(
            filtered_outages,
            bus_to_county,
            left_on="fips_code",
            right_on="FIPS",
            how="left",
        )
        bus_hours = bus_hours[bus_hours["Bus Number"].notna()]
        bus_hours.to_csv("./gtep/data/123_Bus_Resil_Week/not_right.csv")
        self.bus_hours = bus_hours[["hour", "Bus Number"]]
        self.bus_hours = self.bus_hours.astype(int)

    def load_default_data_settings(self):
        ##many of these are hard coded, but they are not set later in the process as of now
        """Fills in necessary but unspecified data information."""
        if "elements" in self.md.data.keys():
            if "generator" in self.md.data["elements"].keys():
                for gen in self.md.data["elements"]["generator"]:
                    # set lifetime value to default first
                    self.md.data["elements"]["generator"][gen]["lifetime"] = 3
                    if "fuel" in self.md.data["elements"]["generator"][gen].keys():
                        if self.md.data["elements"]["generator"][gen]["fuel"] == "C":
                            if (
                                self.md.data["elements"]["generator"][gen]["in_service"]
                                == False
                            ):
                                self.md.data["elements"]["generator"][gen][
                                    "lifetime"
                                ] = 1
                            else:
                                self.md.data["elements"]["generator"][gen][
                                    "lifetime"
                                ] = 2

                    self.md.data["elements"]["generator"][gen][
                        "spinning_reserve_frac"
                    ] = 0.1
                    self.md.data["elements"]["generator"][gen][
                        "quickstart_reserve_frac"
                    ] = 0.1
                    self.md.data["elements"]["generator"][gen]["capital_multiplier"] = 1
                    self.md.data["elements"]["generator"][gen][
                        "extension_multiplier"
                    ] = 0
                    self.md.data["elements"]["generator"][gen][
                        "max_operating_reserve"
                    ] = 1
                    self.md.data["elements"]["generator"][gen][
                        "max_spinning_reserve"
                    ] = 1
                    self.md.data["elements"]["generator"][gen][
                        "max_quickstart_reserve"
                    ] = 1
                    self.md.data["elements"]["generator"][gen]["ramp_up_rate"] = 0.1
                    self.md.data["elements"]["generator"][gen]["ramp_down_rate"] = 0.1
                    self.md.data["elements"]["generator"][gen]["emissions_factor"] = 1
                    self.md.data["elements"]["generator"][gen]["start_fuel"] = 1
                    self.md.data["elements"]["generator"][gen]["investment_cost"] = 1
            if "branch" in self.md.data["elements"].keys():
                for branch in self.md.data["elements"]["branch"]:
                    self.md.data["elements"]["branch"][branch]["loss_rate"] = 0
                    self.md.data["elements"]["branch"][branch]["distance"] = 1
                    self.md.data["elements"]["branch"][branch][
                        "capital_cost"
                    ] = 10000000
        if "system" in self.md.data.keys():
            self.md.data["system"]["min_operating_reserve"] = 0.1
            self.md.data["system"]["min_spinning_reserve"] = 0.1

    def load_storage_csv(self, data_path):
        """Imports storage data.

        :param data_path: filepath for storage data csv file
        """
        try:
            storage_path = data_path + "/storage.csv"
            storage_df = pd.read_csv(storage_path)

            storage_data = {}
            for _, row in storage_df.iterrows():
                name = row["name"]
                storage_data[name] = row.drop("name").to_dict()

            self.md.data["elements"]["storage"] = storage_data
        except FileNotFoundError:
            print(
                f"Warning: The file '{storage_path}' does not exist. Skipping loading storage data."
            )
            self.md.data["elements"]["storage"] = {}

    def texas_case_study_updates(self, data_path):
        """Imports generator data for texas case study.

        :param data_path: filepath for generator data csv file
        """
        # check that datapath is coming from a texas case study directory
        if "Texas" or "Coal" not in data_path:
            raise ValueError("The data path provided is not a Texas case study")

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
