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

# author: Soraya Rawlings and Kyle Skolfield

import pandas as pd
from gtep.config_data import ConfigData
from os.path import join
from typing import Any
from warnings import warn

"""
References

[1] https://atb.nrel.gov/electricity/2024/data 


"""


class DataProcessing:
    """Data processing class for the IDAES GTEP model to read and save
    cost data. Original data from ref[1].
    """

    def __init__(self):

        # Create an instance of the configuration class that has data
        # for generators
        self.config = ConfigData()

    def get_gen_bus_data(self, bus_data_path: str) -> pd.DataFrame:
        """
        :param bus_data_path:   Path to bus data.
        :type bus_data_path:    str
        :returns:               DataFrame with generator bus data.
        """
        # TODO: Confirm that "Gen bus/ Non-gen bus" will always be a column
        # in bus_data? And that it is always an int (instead of a bool)?
        bus_data = pd.read_csv(bus_data_path)
        bus_data = bus_data[bus_data["Gen bus/ Non-gen bus"] == 1]
        return bus_data

    def get_buses_by_gen(
        self,
        gen_bus_df: pd.DataFrame,
        gens: set[str],
    ) -> dict[str, dict]:
        """
        :param gen_bus_df:      DataFrame with generator bus data.
        :type gen_bus_df:       pandas.DataFrame
        :param gens:            Set of candidate generator types.
        :type gens:             set[str]
        :returns:               Dict of the form {gen_type: bus_info}, where bus_info
                                    is itself a dict of the form {bus_name: bus_number}
        """
        # TODO: We should probably check that gen_weight_key is actually a column in gen_bus_df
        # TODO: Confirm that "Bus Name" and "Bus Number" will always be columns in gen_bus_df?
        candidate_gen_by_bus = {}
        for gen in gens:
            gen_weight_key = gen + " - Generation Weight"
            candidate_df = gen_bus_df[
                gen_bus_df[gen_weight_key] > self.config.gen_data[gen]["gen_weight"]
            ]
            candidate_gen_by_bus[gen] = dict(
                zip(candidate_df["Bus Name"], candidate_df["Bus Number"])
            )
        return candidate_gen_by_bus

    def get_clean_gens_dict(self, candidate_gens: list[str]) -> dict[str, str]:
        """
        Maps a list of candidate generators to a cleaned set, where
        "Natural Gas_FE" is substituted for "Natural Gas_CT".

        :param candidate_gens:      List of candidate generator types.
        :type candidate_gens:       list[str]
        :returns:                   Dict mapping input generator types to those without "Natural Gas_CT".
        """
        result = {
            gen: gen if "Natural Gas" not in gen else "Natural Gas_FE"
            for gen in candidate_gens
        }
        if "Natural Gas_CT" in candidate_gens:
            warn(
                "Natural Gas_CT does not have cost data. Substituting Natural Gas_FE cost data."
            )
        return result

    def extract_cost_data(
        self,
        cost_data_path: str,
        gens: list[str],
        years: list[int],
        key_3: str,  # TODO: decide what to call this argument
    ) -> dict[str, pd.DataFrame]:
        """
        :param cost_data_path:      Path to cost data.
        :type cost_data_path:       str
        :param gens:                Generator types to extract data for. Each element
                                        must be a sheet in the cost data file.
        :type gens:                 list[str]
        :param years:               Years to extract data for.
        :type years:                list[int]
        :param key_3:               Key 3
        :type key_3:                str
        :returns:                   Dict of the form {gen_type: cost_data}, where cost_data
                                        is a pandas DataFrame.
        """
        xls = pd.ExcelFile(cost_data_path)
        not_a_sheet = [name for name in gens if name not in xls.sheet_names]
        if not_a_sheet:
            raise ValueError(
                f"All elements of gens must be a sheet in {cost_data_path}. The following are not sheets: {not_a_sheet}"
            )

        cost_data = {gen: pd.read_excel(xls, gen) for gen in gens}
        return {
            gen: df.loc[df["Key3"] == key_3, ["Key1", "Key2"] + years]
            for gen, df in cost_data.items()
        }

    def get_gen_cost_data(
        self,
        cost_df: pd.DataFrame,
        gen_bus_df: pd.DataFrame,
        gen: str,
    ) -> pd.DataFrame:
        """
        Extracts the cost data from gen_cost_df for a given generator type.

        :param cost_df:         Dataframe with cost data.
        :type cost_df:          pandas.DataFrame
        :param gen_bus_df:      Dataframe with bus data.
        :type gen_bus_df:       pandas.DataFrame
        :param gen:             Generator type. Must be a column in `gen_bus_df`.
        :type gen:              str
        :returns:               Dataframe with CapEx data. Has a multi-index of
                                    ("Key1", "Bus Name"), where "Key1" is the
                                    costs variable.
        """
        # link bus name to specific generation class
        gen_subet = gen_bus_df[["Bus Name", gen]]
        gen_cost_df = pd.merge(
            gen_subet,
            cost_df,
            left_on=gen,
            right_on="Key2",
        )
        gen_cost_df = gen_cost_df.set_index(["Key1", "Bus Name"])
        return gen_cost_df

    def build_cost_data_row(
        self,
        gen_cost_df: pd.DataFrame,
        years: list[int],
        cost_var_names: dict[str, str],
        bus_id: int,
        bus_name: str,
        gen: str,
        heat_rate: float,
        ng_costs: dict[int, float],
    ) -> dict[str, Any]:
        """
        Builds a row of cost data and adds it to `gen_data_target`, in-place.

        :param gen_data_target:             Dataframe where the data is being written.
        :type gen_data_target:              pandas.DataFrame
        :param gen_cost_df:                 Dataframe with cost data to be extracted.
        :type gen_cost_df:                  pandas.DataFrame
        :param years:                       Years to include.
        :type years:                        list[int]
        :param cost_var_names:              Dict mapping the output column name to its source in cost data.
        :type cost_var_names:               dict[str, str]
        :param bus_id:                      Bus ID.
        :type bus_id:                       int
        :param bus_name:                    Bus name.
        :type bus_name:                     str
        :param gen:                         Generator type.
        :type gen:                          str
        :param heat_rate:                   Heat rate, in MMBtu/MWh.
        :type heat_rate:                    float
        :param ng_costs:                    Yearly natural gas costs in USD/MMBtu.
        :type ng_costs:                     dict[int, float]
        :returns:                           dict of the form {col_name: value}
        """
        config_gen = self.config.gen_data[gen]
        row = {
            "GEN UID": config_gen["ids"] + str(bus_id) + "-c",
            "Bus ID": bus_id,
            "Unit Type": config_gen["unit_type"],
            "Fuel": config_gen["fuel_type"],
            "PMax MW": config_gen["pmax"],
            "PMin MW": config_gen["pmin"],
            "Min Up Time Hr": config_gen["minup"],
            "Min Down Time Hr": config_gen["mindown"],
        }

        for year in years:
            for out_col, source_var in cost_var_names.items():
                col = f"{out_col}_{year}"
                row[col] = float(gen_cost_df.loc[(source_var, bus_name), year])
            # add natural gas costs
            row[f"fuel_costs_{year}"] = (
                ng_costs[year] * heat_rate if "Natural Gas" in gen else 0
            )

        return row

    def fill_out_prescient_columns(self, df: pd.DataFrame) -> pd.DataFrame:
        """
        Takes an input dataframe with cost data and adds columns required
        by prescient, if they aren't already in the dataframe.

        :param df:  Cost dataframe
        :type df:   pandas.DataFrame
        :returns:   Cost dataframe with all prescient-required columns.
        """
        prescient_cols = [
            "GEN UID",
            "Bus ID",
            "Unit Type",
            "Fuel",
            "MW Inj",
            "MVAR Inj",
            "V Setpoint p.u.",
            "PMax MW",
            "PMin MW",
            "QMax MVAR",
            "QMin MVAR",
            "Min Down Time Hr",
            "Min Up Time Hr",
            "Ramp Rate MW/Min",
            "Start Time Cold Hr",
            "Start Time Warm Hr",
            "Start Time Hot Hr",
            "Start Heat Cold MBTU",
            "Start Heat Warm MBTU",
            "Start Heat Hot MBTU",
            "Non Fuel Start Cost $",
            "Fuel Price $/MMBTU",
            "Output_pct_0",
            "Output_pct_1",
            "Output_pct_2",
            "Output_pct_3",
            "Output_pct_4",
            "HR_avg_0",
            "HR_incr_1",
            "HR_incr_2",
            "HR_incr_3",
            "HR_incr_4",
        ]
        for col in prescient_cols:
            if col not in df.columns:
                df[col] = pd.NA

        df["V Setpoint p.u."] = pd.to_numeric(
            df["V Setpoint p.u."], errors="coerce"
        ).fillna(1)
        df["Output_pct_0"] = pd.to_numeric(df["Output_pct_0"], errors="coerce").fillna(
            0.6
        )
        df["Output_pct_1"] = pd.to_numeric(df["Output_pct_1"], errors="coerce").fillna(
            1
        )

        df.fillna("", inplace=True)
        return df

    def load_gen_data(
        self,
        bus_data_path: str,
        cost_data_path: str,
        candidate_gens: list[str],
        years: list[int] = [2025, 2030, 2035],
        key_3: str = "Moderate",
        heat_rate: float = 9.717,
        ng_costs: dict[int, float] = {
            2025: 3.49,
            2030: 2.91,
            2035: 3.68,
        },
        save_csv: bool = False,
        out_path: str | None = None,
    ):
        """
        Builds a dataframe containing cost data for generators of specified type from bus data.
        Stores the result in `self.gen_data_target` and optionally writes to a CSV.

        :param bus_data_path:               Path to bus data.
        :type bus_data_path:                str
        :param cost_data_path:              Path to cost data.
        :type cost_data_path:               str
        :param candidate_gens:              Generator types to extract cost data for.
        :type candidate_gens:               list[str]
        :param years:                       Years to extract cost data for.
        :type years:                        list[int]
        :param key_3:                       Key 3. Defaults to `"Moderate"`.
        :type key_3:                        str
        :param heat_rate:                   Units of MMBtu/MWh. Defaults to 9.717 (value from [1] assuming a NG
                                                Combustion Turbine (F-Frame), Moderate, for the Base Year 2022).
        :type heat_rate:                    float
        :param ng_costs:                    Prices for naturgal gas in units of USD/MMBtu, in the format {year: cost}.
                                                Defaults to the Henry Hub forecast prices for natural gas.
                                                Each year in `years` must be a key in this dict.
        :type ng_costs:                     dict[str, float]
        :param save_csv:                    Whether to save the resulting dataframe to csv. Defaults to `False`.
        :type save_csv:                     bool
        :param out_path:                    Path to save the csv to. Defaults to `None`, but must be provided if `save_csv=True` is passed.
        :type out_path:                     str | None
        """
        if save_csv and out_path is None:
            raise ValueError("With save_csv=True, out_path must be provided.")

        cost_var_names = {
            "capex": "CAPEX ($/kW)",
            "fixed_ops": "Fixed Operation and Maintenance Expenses ($/kW-yr)",
            "var_ops": "Variable Operation and Maintenance Expenses ($/MWh)",
        }

        years_without_costs = [year for year in years if year not in ng_costs]
        if years_without_costs:
            raise ValueError(
                f"The following years do not have natural gas costs: {years_without_costs}"
            )

        # commenting out for now... not sure we need to do this for just one of the arguments
        # if not isinstance(candidate_gens, list):
        #     raise TypeError(
        #         f"The candidate generators should be in the form of a list. The provided data is of type {type(candidate_gens).__name__}"
        # )

        gen_bus_df = self.get_gen_bus_data(bus_data_path)

        # maps set of generator types in gen_bus_df to set of generator types w/ cost data
        candidate_gens_dict = self.get_clean_gens_dict(candidate_gens)

        buses_by_gen = self.get_buses_by_gen(
            gen_bus_df, set(candidate_gens_dict.keys())
        )
        cost_dict = self.extract_cost_data(
            cost_data_path, set(candidate_gens_dict.values()), years, key_3
        )

        df_rows = []

        # add to the dataframe, row by row
        for bus_gen, cost_gen in candidate_gens_dict.items():
            # NOTE: the current implementation in main has a possible bug: Natural Gas CT are not included in the output csv!
            # the following two lines recreate the bug to allow comparison of remaining outputs
            # if bus_gen != cost_gen:
            #     continue

            gen_cost_df = self.get_gen_cost_data(
                cost_dict[cost_gen], gen_bus_df, bus_gen
            )

            for bus_name, bus in buses_by_gen[bus_gen].items():
                df_rows.append(
                    self.build_cost_data_row(
                        gen_cost_df,
                        years,
                        cost_var_names,
                        bus,
                        bus_name,
                        bus_gen,
                        heat_rate,
                        ng_costs,
                    )
                )

        self.gen_data_target = self.fill_out_prescient_columns(pd.DataFrame(df_rows))

        if save_csv:
            self.gen_data_target.to_csv(
                join(out_path, "candidate_generators_initial_list.csv"),
                index=False,
            )
