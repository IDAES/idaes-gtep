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
import logging
from pathlib import Path
import csv
import pandas as pd
import os
import random

from prescient.simulator.config import PrescientConfig
from prescient.data.providers import gmlc_data_provider

logger = logging.getLogger("gtep.gtep_data")


def _timestamp_from_key(time_key):
    """Convert a Prescient/Egret time key to a pandas Timestamp."""
    return pd.Timestamp(time_key)


def _build_time_key_lookup(time_keys):
    """Build a Timestamp -> original time key lookup.

    The original time key is returned so downstream code continues using the
    exact key representation present in the loaded ModelData.
    """
    return {_timestamp_from_key(time_key): time_key for time_key in time_keys}


def _get_data_years_from_time_keys(time_keys):
    """Return sorted unique years present in loaded time keys."""
    return sorted({_timestamp_from_key(time_key).year for time_key in time_keys})


def _coerce_year_value(value):
    """Try to coerce a CSV value to a plausible calendar year."""
    if value in (None, "", "NA", "N/A"):
        return None

    try:
        year = int(float(value))
    except (TypeError, ValueError):
        return None

    if 1900 <= year <= 2100:
        return year

    return None


def _infer_data_years_from_csv_files(data_path):
    """Infer calendar years present in input data files.

    This scans time-series CSV files for a column named ``Year``. If no
    dedicated time-series directory exists, it falls back to scanning CSV files
    under ``data_path``.

    Returns
    -------
    list[int]
        Sorted unique inferred years. Empty if no years can be inferred.
    """
    data_path = Path(data_path)

    candidate_root = data_path / "timeseries_data_files"
    if not candidate_root.exists():
        candidate_root = data_path

    years = set()

    for csv_path in candidate_root.rglob("*.csv"):
        try:
            with csv_path.open(newline="") as csvfile:
                reader = csv.DictReader(csvfile)

                if not reader.fieldnames:
                    continue

                year_columns = [
                    column
                    for column in reader.fieldnames
                    if column.strip().lower() == "year"
                ]

                if not year_columns:
                    continue

                year_column = year_columns[0]

                for row in reader:
                    year = _coerce_year_value(row.get(year_column))
                    if year is not None:
                        years.add(year)

        except UnicodeDecodeError:
            logger.debug(
                "Skipping non-text or non-UTF8 CSV while inferring data years: %s",
                csv_path,
            )
        except OSError as err:
            logger.debug(
                "Skipping CSV while inferring data years due to OS error: %s (%s)",
                csv_path,
                err,
            )

    return sorted(years)


def _replace_year(timestamp, year):
    """Return timestamp with replaced year, with a useful leap-day error."""
    try:
        return timestamp.replace(year=year)
    except ValueError as err:
        raise ValueError(
            f"Could not replace year in representative date {timestamp} with "
            f"data year {year}. This can happen for dates such as February 29 "
            "when the target data year is not a leap year."
        ) from err


def _format_start_date_for_prescient(timestamp):
    """Format a timestamp as a Prescient-compatible start_date string."""
    if (
        timestamp.hour == 0
        and timestamp.minute == 0
        and timestamp.second == 0
        and timestamp.microsecond == 0
    ):
        return timestamp.strftime("%Y-%m-%d")

    return timestamp.strftime("%Y-%m-%d %H:%M:%S")


def _resolve_prescient_start_date_year(data_path, options_dict):
    """Resolve Prescient start_date year before data are loaded.

    Resolution policy:
        1. If the requested start_date year is present in the data, keep it.
        2. If it is not present and the input data contain exactly one year,
           auto-correct start_date to that year and warn.
        3. If it is not present and the input data contain multiple years, fail.
        4. If no data years can be inferred, leave start_date unchanged.
    """
    if options_dict is None or "start_date" not in options_dict:
        return options_dict

    data_years = _infer_data_years_from_csv_files(data_path)

    if not data_years:
        logger.warning(
            "Could not infer input data year(s) from CSV files under %s. "
            "Leaving Prescient start_date unchanged: %s",
            data_path,
            options_dict["start_date"],
        )
        return options_dict

    requested_start_timestamp = pd.Timestamp(options_dict["start_date"])
    requested_year = requested_start_timestamp.year

    if requested_year in data_years:
        return options_dict

    if len(data_years) == 1:
        data_year = data_years[0]
        corrected_timestamp = _replace_year(requested_start_timestamp, data_year)
        corrected_start_date = _format_start_date_for_prescient(corrected_timestamp)

        logger.warning(
            "Prescient start_date year %s is not present in the input data years "
            "%s. The input data contain a single year (%s), so start_date was "
            "auto-corrected from %r to %r.",
            requested_year,
            data_years,
            data_year,
            options_dict["start_date"],
            corrected_start_date,
        )

        options_dict["start_date"] = corrected_start_date
        return options_dict

    raise ValueError(
        f"Prescient start_date year {requested_year} is not present in the "
        f"input data years {data_years}. The input data contain multiple years, "
        "so the start_date year cannot be safely auto-corrected. Please provide "
        "a start_date whose year exists in the input data."
    )


def _validate_representative_date_windows(
    representative_dates,
    time_keys,
    period_per_step,
):
    """Validate that each representative date can start a full period window."""
    for date in representative_dates:
        if date not in time_keys:
            raise ValueError(
                f"Representative date {date!r} is not present in loaded time keys."
            )

        key_idx = time_keys.index(date)
        if key_idx + period_per_step > len(time_keys):
            raise ValueError(
                f"Representative date {date!r} does not have enough following "
                f"time points to form a full representative period of length "
                f"{period_per_step}. Date index is {key_idx}; total number of "
                f"time keys is {len(time_keys)}."
            )


def _resolve_user_representative_dates(
    representative_dates,
    time_keys,
    num_reps,
    period_per_step,
):
    """Resolve and validate user-provided representative dates.

    Resolution policy:
        1. Use exact user-specified timestamps if all are present.
        2. If some are missing and the loaded data contains exactly one year,
           auto-correct the year while preserving month/day/time.
        3. Otherwise fail with a clear error.
    """
    if len(representative_dates) != num_reps:
        raise ValueError(
            f"The number of provided representative_dates must match num_reps. "
            f"Received len(representative_dates)={len(representative_dates)}, "
            f"but num_reps={num_reps}."
        )

    representative_dates = list(representative_dates)

    missing_dates = [date for date in representative_dates if date not in time_keys]

    if not missing_dates:
        _validate_representative_date_windows(
            representative_dates,
            time_keys,
            period_per_step,
        )
        return representative_dates

    data_years = _get_data_years_from_time_keys(time_keys)

    if len(data_years) != 1:
        raise ValueError(
            "The following representative_dates are not valid timestamps in the "
            f"loaded input data: {missing_dates}. The loaded data contains "
            f"multiple years {data_years}, so the representative-date year "
            "cannot be safely auto-corrected. Please provide timestamps that "
            "exist in the loaded data."
        )

    data_year = data_years[0]
    time_key_lookup = _build_time_key_lookup(time_keys)

    corrected_dates = []
    corrections = []

    for date in representative_dates:
        if date in time_keys:
            corrected_dates.append(date)
            continue

        requested_timestamp = _timestamp_from_key(date)
        corrected_timestamp = _replace_year(requested_timestamp, data_year)

        if corrected_timestamp not in time_key_lookup:
            raise ValueError(
                f"Representative date {date!r} is not present in the loaded "
                f"data, and auto-correcting its year to {data_year} produced "
                f"{str(corrected_timestamp)!r}, which is also not present. "
                "Please provide a representative date that exists in the "
                "loaded data."
            )

        corrected_key = time_key_lookup[corrected_timestamp]
        corrected_dates.append(corrected_key)
        corrections.append((date, corrected_key))

    if len(set(corrected_dates)) != len(corrected_dates):
        raise ValueError(
            "Representative-date year auto-correction produced duplicate "
            f"representative dates: {corrected_dates}. Please provide distinct "
            "representative dates."
        )

    _validate_representative_date_windows(
        corrected_dates,
        time_keys,
        period_per_step,
    )

    logger.warning(
        "Some representative_dates were not present in the loaded data, but the "
        "loaded data contains a single year (%s). Auto-corrected representative "
        "date years while preserving month/day/time: %s",
        data_year,
        corrections,
    )

    return corrected_dates


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
                "start_date": "2019-01-01",
            }

        else:
            # ensure data path is included in options dictionary
            options_dict["data_path"] = data_path

        # Resolve the Prescient start_date year before loading data. If the
        # configured year is not present in the input data, but the data contain
        # a single year, auto-correct to that year with a warning. If multiple
        # years are present, fail with a clear error.
        options_dict = _resolve_prescient_start_date_year(data_path, options_dict)

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

        # Get the timestamps in the loaded day-ahead data. Default
        # representative_dates are selected from this list to ensure
        # they correspond to valid input data timestamps. if
        # representative_dates are provided by the user, those values
        # are used instead.
        time_keys = self.md.data["system"]["time_keys"]

        if representative_dates is None:
            available_day_starts = time_keys[::period_per_step]

            if len(available_day_starts) < self.num_reps:
                raise ValueError(
                    "Not enough available day-start timestamps to select default representative dates. Please provide a custom list of representative_dates in the driver or reduce num_reps."
                )

            # Pick 4 default representative dates from the loaded
            # data: winter, spring, summer, and fall. If self.num_reps
            # < 4, use only the first self.num_reps default dates. If
            # more than 4 representative periods are requested, keep
            # these 4 defaults and randomly select the remaining dates
            # from the available day-start timestamps.
            default_representative_dates = [
                available_day_starts[27],  # 2020-01-28
                available_day_starts[113],  # 2020-04-23
                available_day_starts[186],  # 2020-07-05
                available_day_starts[287],  # 2020-10-14
            ]

            if self.num_reps <= 4:
                representative_dates = default_representative_dates[: self.num_reps]
            else:
                if len(available_day_starts) < self.num_reps:
                    raise ValueError(
                        "Not enough available day-start timestamps to select default representative dates. "
                        "Please provide a custom list of representative_dates in the driver or reduce num_reps."
                    )

                random_seed = 42
                rng = random.Random(random_seed)
                remaining_dates = [
                    date
                    for date in available_day_starts
                    if date not in default_representative_dates
                ]
                additional_dates = rng.sample(
                    remaining_dates,
                    self.num_reps - len(default_representative_dates),
                )
                representative_dates = sorted(
                    default_representative_dates + additional_dates,
                    key=lambda date: time_keys.index(date),
                )
        else:
            representative_dates = _resolve_user_representative_dates(
                representative_dates,
                time_keys,
                self.num_reps,
                period_per_step,
            )

        self.representative_dates = representative_dates

        if representative_weights:

            if len(representative_dates) != len(representative_weights):
                raise ValueError(
                    "Length of representative_dates and representative_weights must match."
                )
            else:
                print(
                    "INFO: representative_dates and representative_weights are aligned. Continue building the data modeling object..."
                )

            # Store as a dictionary
            self.representative_weights_dict = dict(
                zip(representative_dates, representative_weights)
            )

        else:

            # Set weight for each representative day to default value
            # of 1. The other option is to set the weight for each day
            # to the total weight divided by the number of
            # representative dates.
            set_default_weight = True
            if set_default_weight:
                weight_per_date = 1
            else:
                total_weight = prescient_options.num_days * self.stages
                weight_per_date = int(total_weight / len(representative_dates))

            # Store weights as a dictionary by representative date
            self.representative_weights_dict = {
                date: weight_per_date for date in self.representative_dates
            }

        # Read average heat rates from the "HR_avg_0" column in
        # gen.csv and assign them to each generator in self.md. This
        # is done manually because the generator data loaded into
        # self.md does not include heat rate values by default. Units
        # should be in MMBTU/MWh.
        gen_csv_file = os.path.join(data_path, "gen.csv")
        heat_rate_dict = {}
        static_pmax_dict = {}

        static_pmax_columns = [
            "PMax MW",
            "Pmax MW",
            "P_MAX MW",
            "P_MAX",
            "PMax",
            "p_max",
            "PMax_MW",
            "PMAX",
        ]

        with open(gen_csv_file, newline="") as csvfile:
            reader = csv.DictReader(csvfile)
            for row in reader:
                gen_uid = row.get("GEN UID")
                heat_rate_str = row.get("HR_avg_0")

                if heat_rate_str not in (None, "", "NA"):
                    heat_rate = float(heat_rate_str)
                else:
                    heat_rate = None

                if gen_uid and heat_rate is not None:
                    heat_rate_dict[gen_uid] = heat_rate

                if gen_uid:
                    for column in static_pmax_columns:
                        pmax_str = row.get(column)
                        if pmax_str not in (None, "", "NA"):
                            try:
                                static_pmax_dict[gen_uid] = float(pmax_str)
                            except ValueError:
                                logger.warning(
                                    "Could not parse static p_max value %r from "
                                    "column %s for generator %s.",
                                    pmax_str,
                                    column,
                                    gen_uid,
                                )
                            break

        for gen in self.md.data["elements"]["generator"]:
            if gen in heat_rate_dict:
                self.md.data["elements"]["generator"][gen]["heat_rate"] = (
                    heat_rate_dict[gen]
                )
            else:
                self.md.data["elements"]["generator"][gen]["heat_rate"] = 0

            if gen in static_pmax_dict:
                # Preserve static physical nameplate capacity from gen.csv.
                # This is scenario-invariant and should not be inferred from a
                # selected representative day's renewable availability profile.
                self.md.data["elements"]["generator"][gen]["nameplate_capacity"] = (
                    static_pmax_dict[gen]
                )

        thermal_heat_rates = [
            self.md.data["elements"]["generator"][gen].get("heat_rate", 0)
            for gen in self.md.data["elements"]["generator"]
            if self.md.data["elements"]["generator"][gen].get("generator_type")
            == "thermal"
        ]

        if thermal_heat_rates and all(hr == 0 for hr in thermal_heat_rates):
            logger.info(
                "All thermal generators have heat_rate values equal to 0. "
                "Please re-check the input data. Fuel costs are multiplied by "
                "heat_rate, so resulting fuel costs will all be 0."
            )

        # IMPORTANT TO READ: Always add or modify any new elements in
        # self.md.data (such as new time series or parameters) BEFORE
        # creating representative_data using clone_at_time_keys. This
        # ensures all representative ModelData objects will have the
        # new elements.
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
                    self.md.data["elements"]["generator"][gen].setdefault(
                        "non_fuel_startup_cost", 0
                    )
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
