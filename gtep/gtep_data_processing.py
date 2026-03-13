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

    cost_var_names = {
        "capex": "CAPEX ($/kW)",
        "fixed_ops": "Fixed Operation and Maintenance Expenses ($/kW-yr)",
        "var_ops": "Variable Operation and Maintenance Expenses ($/MWh)",
    }
    heat_rate_var = "Heat Rate  (MMBtu/MWh)"

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
        gens: list[str],
    ) -> dict[str, dict]:
        """
        :param gen_bus_df:      DataFrame with generator bus data.
        :type gen_bus_df:       pandas.DataFrame
        :param gens:            List of candidate generator types.
        :type gens:             list[str]
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
                "Natural Gas_CT does not have cost data. We will assume all natural gas costs are coming from fossil energy sources (Natural Gas_FE)."
            )
        return result

    def extract_cost_data(
        self,
        cost_data_path: str,
        gens: list[str],
        years: list[int],
        scenario: str,
    ) -> dict[str, pd.DataFrame]:
        """
        :param cost_data_path:      Path to cost data.
        :type cost_data_path:       str
        :param gens:                Generator types to extract data for. Each element
                                        must be a sheet in the cost data file.
        :type gens:                 list[str]
        :param years:               Years to extract data for.
        :type years:                list[int]
        :param scenario:            Cost scenario.
        :type scenario:             str
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
            gen: df.loc[df["Key3"] == scenario, ["Key1", "Key2"] + years]
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
        gen_subset = gen_bus_df[["Bus Name", gen]]
        gen_cost_df = pd.merge(
            gen_subset,
            cost_df,
            left_on=gen,
            right_on="Key2",
        )
        gen_cost_df = gen_cost_df.set_index(["Key1", "Bus Name"])
        gen_cost_df = gen_cost_df.drop(columns=["Key2", gen])
        return gen_cost_df

    def build_cost_data_row(
        self,
        gen_cost_df: pd.DataFrame,
        years: list[int],
        bus_id: int,
        bus_name: str,
        gen: str,
        ng_costs: dict[int, float],
    ) -> dict[str, Any]:
        """
        Builds a row of cost data and adds it to `gen_data_target`, in-place.

        :param gen_cost_df:                 Dataframe with cost data to be extracted.
        :type gen_cost_df:                  pandas.DataFrame
        :param years:                       Years to include.
        :type years:                        list[int]
        :param bus_id:                      Bus ID.
        :type bus_id:                       int
        :param bus_name:                    Bus name.
        :type bus_name:                     str
        :param gen:                         Generator type.
        :type gen:                          str
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
            for out_col, source_var in self.cost_var_names.items():
                col = f"{out_col}_{year}"
                row[col] = float(gen_cost_df.loc[(source_var, bus_name), year])
            # add natural gas costs
            row[f"fuel_costs_{year}"] = (
                ng_costs[year] * float(gen_cost_df.loc[(self.heat_rate_var, bus_name), year])
                if "Natural Gas" in gen
                else 0.0
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
        for col in self.prescient_cols:
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
        scenario: str = "Moderate",
        ng_costs: dict[int, float] = {  # can't find values this far out online
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
        :param scenario:                    Cost scenario. Defaults to `"Moderate"`.
        :type scenario:                     str
        :param ng_costs:                    Natural gas costs in units of USD/MMBtu, in the format {year: cost}.
                                                Defaults to the Henry Hub forecast prices for natural gas.
                                                Each year in `years` must be a key in this dict.
        :type ng_costs:                     dict[str, float]
        :param save_csv:                    Whether to save the resulting dataframe to csv. Defaults to `False`.
        :type save_csv:                     bool
        :param out_path:                    Path to save the csv to. Defaults to `None`, but must be provided if `save_csv=True` is passed.
        :type out_path:                     str | None
        """
        if save_csv and out_path is None:
            raise TypeError("With save_csv=True, out_path must be a str.")

        years_without_costs = [year for year in years if year not in ng_costs]
        gens_include_ng = any(["Natural Gas" in gen for gen in candidate_gens])
        if years_without_costs and gens_include_ng:
            raise KeyError(
                f"Natural gas generators were passed in candidate_gens, but the following years do not have natural gas costs: {years_without_costs}"
            )

        gen_bus_df = self.get_gen_bus_data(bus_data_path)

        # maps set of generator types in gen_bus_df to set of generator types w/ cost data
        candidate_gens_dict = self.get_clean_gens_dict(candidate_gens)

        buses_by_gen = self.get_buses_by_gen(
            gen_bus_df, set(candidate_gens_dict.keys())
        )
        cost_dict = self.extract_cost_data(
            cost_data_path, set(candidate_gens_dict.values()), years, scenario
        )

        # build output dataframe
        df_rows = []
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
                        bus,
                        bus_name,
                        bus_gen,
                        ng_costs,
                    )
                )
        self.gen_data_target = self.fill_out_prescient_columns(pd.DataFrame(df_rows))

        if save_csv:
            self.gen_data_target.to_csv(out_path, index=False)
