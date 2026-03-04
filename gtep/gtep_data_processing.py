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
from pathlib import Path
from gtep.config_data import ConfigData
from warnings import warn

"""
References

[1] https://atb.nrel.gov/electricity/2024/data 


"""


class DataProcessing:
    """Data processing class for the IDAES GTEP model to read and save
    cost data. Original data from ref[1].

    """
    # For now, moved all the hard-coded stuff to here.
    all_gens = [
        "Natural Gas_CT",
        "Natural Gas_FE",
        "Solar - CSP",
        "Solar - Utility PV",
        "Land-Based Wind",
        "Coal_FE",
        "Biopower",
        "Nuclear",
        "Geothermal",
    ]
    years_list = [2025, 2030, 2035]
    subvar_of_interest = "Moderate"

    # Henry Hub forecast prices for natural
    # gas. [ESR WIP: Units for this metric are
    # in USD/MMBtu. To convert it to $/MWh,
    # multiply factor by the heat rate. In
    # this case, we use a heat rate value from
    # [1] assuming a NG Combustion Turbine
    # (F-Frame), Moderate, for the Base Year
    # 2022. The units are in MMBtu/MWh.]
    heat_rate = 9.717
    hh_ng_costs = {
        2025: 3.49,
        2030: 2.91,
        2035: 3.68,
    }

    # TODO: A bunch of these columns are never set... while others are set that aren't here!
    target_columns = [
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
        # TODO: Confirm that this functionality of filtering for gen buses
        # is what we want to maintain?
        bus_data = pd.read_csv(bus_data_path)
        bus_data = bus_data[bus_data["Gen bus/ Non-gen bus"] == 1]
        return bus_data

    def get_buses_by_gen(self,
        gen_bus_df: pd.DataFrame,
        candidate_gens: list[str],
    ) -> dict[str, dict]:
        """
        :param gen_bus_df:      DataFrame with generator bus data.
        :type gen_bus_df:       pandas.DataFrame
        :param candidate_gens:  List of candidate generator types.
        :type candidate_gens:   list[str]
        :returns:               Dict of the form {gen_type: bus_info}, where bus_info
                                    is itself a dict of the form {bus_name: bus_number}
        """
        # TODO: We should probably check that gen_weight_key is actually a column in gen_bus_df
        # TODO: Confirm that "Bus Name" and "Bus Number" will always be columns in gen_bus_df?
        # TODO: Might be good to pass a set instead of list for candidate gens to prevent over-writing
        candidate_gen_by_bus = {}
        for gen in candidate_gens:
            gen_weight_key = gen + " - Generation Weight"
            candidate_df = gen_bus_df[
                gen_bus_df[gen_weight_key] > self.config.gen_data[gen]["gen_weight"]
            ]
            candidate_gen_by_bus[gen] = dict(
                zip(candidate_df["Bus Name"], candidate_df["Bus Number"])
            )
        return candidate_gen_by_bus
    
    def clean_up_candidate_gens(self, candidate_gens: list[str]) -> set[str]:
        """
        Ensures elements of a given list are unique and do not include
        "Natural Gas_CT" ("Natural Gas_FE" will be substituted).

        :param candidate_gens:      List of candidate generator types.
        :type candidate_gens:       list[str]
        :returns:                   Set of generator types.
        """
        # TODO: Should we ensure that these are all in self.all_gens?
        return set([gen if "Natural Gas" not in gen else "Natural Gas_FE" for gen in candidate_gens])
    
    def get_cost_data(self, cost_data_path: str, gens_of_interest: set[str]) -> dict[str, pd.DataFrame]:
        """
        :param cost_data_path:      Path to cost data.
        :type cost_data_path:       str
        :param gens_of_interest:    Generator types of interest.
        :type gens_of_interest:     set[str]
        :returns:                   Dict of the form {gen_type: cost_data}, where cost_data
                                        is a pandas DataFrame.
        """
        # TODO: we should probably check that the elements of gens_of_interest are actually books...
        return {
            book: pd.read_excel(cost_data_path, book)
            for book in gens_of_interest
        }

    def get_gen_cost_data(
            self,
            cost_df: pd.DataFrame,
            gen_bus_df: pd.DataFrame,
            gen: str) -> pd.DataFrame:
        """
        Extracts the cost data from gen_cost_df for a given generator type.

        :param cost_df:         Dataframe with cost data.
        :type cost_df:          pandas.DataFrame
        :param gen_bus_df:      Dataframe with bus data.
        :type gen_bus_df:       pandas.DataFrame
        :param gen:             Type of generator.
        :type gen:              str
        :returns:               Dataframe with CapEx data. Has a multi-index of
                                    ("Key1", "Bus Name"), where "Key1" is the
                                    costs variable.
        """
        gen_cost_df = cost_df[
            cost_df["Key3"] == self.subvar_of_interest
        ]
        # link bus name to specific generation class
        gen_cost_df = pd.merge(
            gen_bus_df[["Bus Name", gen]],
            gen_cost_df,
            left_on=gen,
            right_on="Key2",
        )
        gen_cost_df = gen_cost_df.set_index(["Key1", "Bus Name"])
        return gen_cost_df

    def load_gen_data(
        self,
        bus_data_path: str,
        cost_data_path: str,
        candidate_gens: list[str],
        save_csv: bool = False,
    ):

        if not isinstance(
            candidate_gens, list
        ):  # is there a reason we're checking only this variable's type?
            raise TypeError(
                f"The candidate generators should be in the form of a list. The provided data is of type {type(candidate_gens).__name__}"
            )
        
        gen_bus_df = self.get_gen_bus_data(bus_data_path)
        buses_by_gen = self.get_buses_by_gen(gen_bus_df, candidate_gens)
        
        gens_of_interest = self.clean_up_candidate_gens(candidate_gens)
        cost_dict = self.get_cost_data(cost_data_path, gens_of_interest)

        # instantiate dataframe that will have all of the output data
        cost_col_names = []
        for year in self.years_list:
            cost_col_names += [f"capex_{year}", f"fuel_costs_{year}", f"fixed_ops_{year}", f"var_ops_{year}"]
        self.gen_data_target = pd.DataFrame(columns=self.target_columns + cost_col_names)

        # add to the dataframe, row by row
        for gen, cost_df in cost_dict.items():
            gen_cost_df = self.get_gen_cost_data(cost_df, gen_bus_df, gen)

            for bus_name, bus in buses_by_gen[gen].items():
                bus_row = {
                    "GEN UID": self.config.gen_data[gen]["ids"] + str(bus) + "-c",
                    "Bus ID": bus,
                    "Unit Type": self.config.gen_data[gen]["unit_type"],
                    "Fuel": self.config.gen_data[gen]["fuel_type"],
                    "PMax MW": self.config.gen_data[gen]["pmax"],
                    "PMin MW": self.config.gen_data[gen]["pmin"],
                    "Min Up Time Hr": self.config.gen_data[gen]["minup"],
                    "Min Down Time Hr": self.config.gen_data[gen]["mindown"],
                }

                for year in self.years_list:
                    # maps variable in cost dataset to what our output cols are
                    col_mapper = {
                        "CAPEX ($/kW)": f"capex_{year}",
                        "Fixed Operation and Maintenance Expenses ($/kW-yr)": f"fixed_ops_{year}",
                        "Variable Operation and Maintenance Expenses ($/MWh)": f"var_ops_{year}",
                    }
                    for var, col in col_mapper.items():
                        bus_row[col] = float(gen_cost_df.loc[(var, bus_name), year])
                    # add heat rate data as well
                    if gen == "Natural Gas_FE":
                        bus_row[f"fuel_costs_{year}"] = self.hh_ng_costs[year] * self.heat_rate
                    else:
                        bus_row[f"fuel_costs_{year}"] = 0

                self.gen_data_target.loc[len(self.gen_data_target)] = bus_row

        # commenting the rest of this out for now bc we never set any of these other variables...

        # self.gen_data_target["V Setpoint p.u."] = pd.to_numeric(
        #     self.gen_data_target["V Setpoint p.u."], errors="coerce"
        # ).fillna(1)
        # self.gen_data_target["Output_pct_0"] = pd.to_numeric(
        #     self.gen_data_target["Output_pct_0"], errors="coerce"
        # ).fillna(0.6)
        # self.gen_data_target["Output_pct_1"] = pd.to_numeric(
        #     self.gen_data_target["Output_pct_1"], errors="coerce"
        # ).fillna(1)

        # self.gen_data_target.fillna("", inplace=True)

        if save_csv:
            self.gen_data_target.to_csv(
                "data/costs/candidate_generators_initial_list.csv", index=False
            )
