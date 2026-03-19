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

from typing import Any
from warnings import warn
import pandas as pd
from gtep.config_data import ConfigData
from pathlib import Path

"""
References

[1] https://atb.nrel.gov/electricity/2024/data 
[2] https://www.eia.gov/outlooks/aeo/data/browser/#/?id=1-AEO2025&region=0-0&cases=ref2025~hm2025~lm2025~highprice~lowprice~highogs~lowogs~highZTC~lowZTC~nocaa111~alttrnp~aeo2023ref&start=2023&end=2050&f=A&linechart=~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ref2025-d032025a.44-1-AEO2025&map=&ctype=linechart&sid=&sourcekey=0

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

    def get_gen_bus_data(self, bus_data_path: Path) -> pd.DataFrame:
        """
        :param bus_data_path:   Path to bus data.
        :type bus_data_path:    pathlib.Path
        :returns:               DataFrame with generator bus data.
        """
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
        :param gens:            List of candidate generator types.
        :type gen_bus_df:       pandas.DataFrame
        :type gens:             list[str]
        :returns:               Dict of the form {gen_type: bus_info}, where bus_info
                                    is itself a dict of the form {bus_name: bus_number}
        """
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
        cost_data_path: Path,
        gens: list[str],
        years: list[int],
        scenario: str,
    ) -> dict[str, pd.DataFrame]:
        """
        :param cost_data_path:      Path to cost data.
        :param gens:                Generator types to extract data for. Each element
                                        must be a sheet in the cost data file.
        :param years:               Years to extract data for.
        :param scenario:            Cost scenario.
        :type cost_data_path:       pathlib.Path
        :type gens:                 list[str]
        :type years:                list[int]
        :type scenario:             str
        :returns:                   Dict of the form {gen_type: cost_data}, where cost_data
                                        is a pandas DataFrame.
        """
        xls = pd.ExcelFile(cost_data_path)

        # check all sheets are available
        not_a_sheet = [name for name in gens if name not in xls.sheet_names]
        if not_a_sheet:
            raise ValueError(
                f"All elements of gens must be a sheet in {cost_data_path}. The following are not sheets: {not_a_sheet}"
            )

        cost_data = {gen: pd.read_excel(xls, gen) for gen in gens}

        # check all years available
        for gen, df in cost_data.items():
            missing_years = [year for year in years if year not in df.columns]
            if missing_years:
                raise ValueError(
                    f"One or more years ({missing_years}) were not in the {gen} cost dataset. Available columns: {df.columns}"
                )

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
        :param gen_bus_df:      Dataframe with bus data.
        :param gen:             Generator type. Must be a column in `gen_bus_df`.
        :type cost_df:          pandas.DataFrame
        :type gen_bus_df:       pandas.DataFrame
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

    def get_ng_costs(self, ng_cost_path: Path) -> pd.DataFrame:
        """
        Reads in natural gas costs data, sourced from [2], with the
        quantity (e.g., Reference case) on the index.

        :param ng_cost_path:    Path to CSV with natural gas costs.
        :type ng_cost_path:     pathlib.Path
        :returns:               pandas.DataFrame with natural gas costs
        """
        rowfunc = lambda row: ((row > 425 or row < 414) and row != 4)
        ng_data = pd.read_csv(ng_cost_path, skiprows=rowfunc, header=0, index_col=0)
        return ng_data

    def build_cost_data_row(
        self,
        gen_cost_df: pd.DataFrame,
        ng_cost_df: pd.DataFrame,
        ng_cost_quantity: str,
        years: list[int],
        bus_id: int,
        bus_name: str,
        gen: str,
    ) -> dict[str, Any]:
        """
        Builds a row of cost data and adds it to `gen_data_target`, in-place. `ng_cost_df` and
        `ng_cost_quantity` must be provided if `gen` contains the substring `"Natural Gas"`.

        :param gen_cost_df:                 Dataframe with cost data to be extracted.
        :param ng_cost_df:                  Dataframe with natural gas costs in 2024 USD/MMBtu.
        :param years:                       Years to include.
        :param bus_id:                      Bus ID.
        :param bus_name:                    Bus name.
        :param gen:                         Generator type.
        :param ng_cost_quantity:            Natural gas cost quantity to use (e.g., `"Reference case"`).
        :type gen_cost_df:                  pandas.DataFrame
        :type ng_cost_df:                   pandas.DataFrame
        :type years:                        list[int]
        :type ng_cost_quantity:             str
        :type bus_id:                       int
        :type bus_name:                     str
        :type gen:                          str
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
                ng_cost_df.loc[ng_cost_quantity, str(year)]
                * float(gen_cost_df.loc[(self.heat_rate_var, bus_name), year])
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
        bus_data_path: Path,
        cost_data_path: Path,
        ng_cost_path: Path,
        candidate_gens: list[str],
        years: list[int] = [2025, 2030, 2035],
        scenario: str = "Moderate",
        ng_cost_quantity: str = "Reference case",
        save_csv: bool = False,
        out_path: Path | None = None,
    ):
        """
        Builds a dataframe containing cost data for generators of specified type from bus data.
        Stores the result in `self.gen_data_target` and optionally writes to a CSV.

        :param bus_data_path:               Path to bus data.
        :param cost_data_path:              Path to cost data.
        :param ng_cost_path:                Path to natural gas cost data.
        :param candidate_gens:              Generator types to extract cost data for.
        :param years:                       Years to extract cost data for.
        :param scenario:                    Cost scenario. Defaults to `"Moderate"`.
        :param ng_cost_quantity:            Natural gas cost quantity to use. Defaults to `"Reference case"`.
        :param save_csv:                    Whether to save the resulting dataframe to csv. Defaults to `False`.
        :param out_path:                    Directory to save the csv to. Defaults to `None`, but must be provided if `save_csv=True` is passed.
        :type bus_data_path:                pathlib.Path
        :type cost_data_path:               pathlib.Path
        :type ng_cost_path:                 pathlib.Path
        :type candidate_gens:               list[str]
        :type years:                        list[int]
        :type scenario:                     str
        :type ng_cost_quantity:             str
        :type save_csv:                     bool
        :type out_path:                     pathlib.Path | None
        """
        if save_csv and out_path is None:
            raise TypeError("With save_csv=True, out_path must be a provided.")

        gen_bus_df = self.get_gen_bus_data(bus_data_path)

        # get natural gas data and check we have the matching cost quantity and years
        ng_cost_df = self.get_ng_costs(ng_cost_path)
        if ng_cost_quantity not in ng_cost_df.index:
            raise ValueError(
                f"The provided natural gas quantity {ng_cost_quantity} is not in the natural gas cost dataset. Available quantities: {ng_cost_df.index}"
            )
        missing_years = [year for year in years if str(year) not in ng_cost_df.columns]
        if missing_years:
            raise ValueError(
                f"One or more years ({missing_years}) were not in the natural gas cost dataset. Available columns: {ng_cost_df.columns}"
            )

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
                        gen_cost_df=gen_cost_df,
                        years=years,
                        bus_id=bus,
                        bus_name=bus_name,
                        gen=bus_gen,
                        ng_cost_df=ng_cost_df,
                        ng_cost_quantity=ng_cost_quantity,
                    )
                )
        self.gen_data_target = self.fill_out_prescient_columns(pd.DataFrame(df_rows))

        if save_csv:
            self.gen_data_target.to_csv((out_path / "costs.csv").resolve(), index=False)
