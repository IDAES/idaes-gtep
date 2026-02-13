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

    def load_gen_data(
        self,
        bus_data_path=None,
        cost_data_path=None,
        candidate_gens=None,
        save_csv=None,
    ):

        assert isinstance(
            candidate_gens, list
        ), f"The candidate generators should be in the form of a list. The provided data is of type {type(candidate_gens).__name__}"

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

        candidate_gens_filt = candidate_gens

        # Read bus data from the provided .csv file
        gen_nonsense = pd.read_csv(bus_data_path)
        gen_nonsense.head()

        gen_nonsense = gen_nonsense[gen_nonsense["Gen bus/ Non-gen bus"] == 1]
        gen_nonsense.head()

        candidate_df = {}
        candidate_gen_by_bus = {}
        for gen in candidate_gens_filt:
            gen_weight_key = gen + " - Generation Weight"
            candidate_df[gen] = gen_nonsense[
                gen_nonsense[gen_weight_key] > self.config.gen_data[gen]["gen_weight"]
            ]

            candidate_gen_by_bus[gen] = dict(
                zip(candidate_df[gen]["Bus Name"], candidate_df[gen]["Bus Number"])
            )

        # Considering that Natural Gas_CT and Natural Gas_FE have the
        # same cost assumptions, create a pricing_data list that
        # contains only Natural Gas_FE
        pricing_data_book_names = [gen for gen in all_gens if gen != "Natural Gas_CT"]
        pricing_data_sheet_path = Path(cost_data_path)
        pricing_dict = {}
        for this_book_name in pricing_data_book_names:
            pricing_dict[this_book_name] = pd.read_excel(
                pricing_data_sheet_path, this_book_name
            )

        gen_uids = []
        bus_uids = []
        unit_types = []
        fuels = []
        pmax = []
        pmin = []
        minup = []
        mindown = []

        years_list = [2025, 2030, 2035]
        costs = {}
        for year in years_list:
            costs[f"capex_{year}"] = []
            costs[f"fuel_costs_{year}"] = []
            costs[f"fixed_ops_{year}"] = []
            costs[f"var_ops_{year}"] = []

        bus_data_df = gen_nonsense

        gens_of_interest_unfilt = []
        for gen in candidate_gens_filt:
            if "Natural Gas" in gen:
                gens_of_interest_unfilt.append("Natural Gas_FE")
            else:
                gens_of_interest_unfilt.append(gen)

        # Since we can use two different natural gas generators, add a
        # filter to remove the repeated natural gas generator (this
        # because we only have cost data for the "_FE" kind).
        gens_of_interest = []
        for item in gens_of_interest_unfilt:
            if item not in gens_of_interest:
                gens_of_interest.append(item)

        vars_of_interest = [
            "CAPEX ($/kW)",
            "Heat Rate  (MMBtu/MWh)",
            "Fixed Operation and Maintenance Expenses ($/kW-yr)",
            "Variable Operation and Maintenance Expenses ($/MWh)",
        ]
        subvars_of_interest = "Moderate"

        for var in vars_of_interest:
            for gen in gens_of_interest:
                demo_df = pricing_dict[gen][
                    (pricing_dict[gen]["Key1"] == var)
                    & (pricing_dict[gen]["Key3"] == subvars_of_interest)
                ]
                bus_cutdown_df = bus_data_df[["Bus Name", gen]]
                the_thing_we_want = bus_cutdown_df[gen].iloc[0]
                the_row_that_has_the_prices = demo_df[
                    demo_df["Key2"] == the_thing_we_want
                ]
                bus_capex_df = pd.merge(
                    bus_cutdown_df, demo_df, left_on=gen, right_on="Key2"
                )

                for bus_name, bus in candidate_gen_by_bus[gen].items():
                    if var == "CAPEX ($/kW)":
                        gen_uids.append(
                            self.config.gen_data[gen]["ids"] + str(bus) + "-c"
                        )
                        bus_uids.append(bus)
                        unit_types.append(self.config.gen_data[gen]["unit_type"])
                        fuels.append(self.config.gen_data[gen]["fuel_type"])
                        pmax.append(self.config.gen_data[gen]["pmax"])
                        pmin.append(self.config.gen_data[gen]["pmin"])
                        minup.append(self.config.gen_data[gen]["minup"])
                        mindown.append(self.config.gen_data[gen]["mindown"])

                        for year in years_list:
                            costs[f"capex_{year}"].append(
                                float(
                                    bus_capex_df[bus_capex_df["Bus Name"] == bus_name][
                                        year
                                    ].iloc[0]
                                )
                            )

                    elif var == "Heat Rate  (MMBtu/MWh)":
                        if gen == "Natural Gas_FE":
                            # Henry Hub forecast prices for natural
                            # gas. [ESR WIP: Units for this metric are
                            # in USD/MMBtu. To convert it to $/MWh,
                            # multiply factor by the heat rate. In
                            # this case, we use a heat rate value from
                            # [1] assuming a NG Combustion Turbine
                            # (F-Frame), Moderate, for the Base Year
                            # 2022. The units are in MMBtu/MWh.]
                            heat_rate = 9.717
                            hh_ng_costs = {2025: 3.49, 2030: 2.91, 2035: 3.68}

                            for year in years_list:
                                # [ESR WIP: Comment original equation]
                                # costs[f'fuel_costs_{year}'].append(hh_ng_costs[year] * float(bus_capex_df[bus_capex_df["Bus Name"] == bus_name][year].iloc[0]))
                                costs[f"fuel_costs_{year}"].append(
                                    hh_ng_costs[year] * heat_rate
                                )  # units in USD/MWh

                        else:
                            for year in years_list:
                                costs[f"fuel_costs_{year}"].append(0)

                    elif var == "Fixed Operation and Maintenance Expenses ($/kW-yr)":
                        for year in years_list:
                            costs[f"fixed_ops_{year}"].append(
                                float(
                                    bus_capex_df[bus_capex_df["Bus Name"] == bus_name][
                                        year
                                    ].iloc[0]
                                )
                            )

                    elif var == "Variable Operation and Maintenance Expenses ($/MWh)":
                        for year in years_list:
                            costs[f"var_ops_{year}"].append(
                                float(
                                    bus_capex_df[bus_capex_df["Bus Name"] == bus_name][
                                        year
                                    ].iloc[0]
                                )
                            )

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

        self.gen_data_target = pd.DataFrame(columns=target_columns)
        self.gen_data_target["GEN UID"] = gen_uids
        self.gen_data_target["Bus ID"] = bus_uids
        self.gen_data_target["Unit Type"] = unit_types
        self.gen_data_target["Fuel"] = fuels
        self.gen_data_target["PMax MW"] = pmax
        self.gen_data_target["PMin MW"] = pmin
        self.gen_data_target["Min Up Time Hr"] = minup
        self.gen_data_target["Min Down Time Hr"] = mindown
        for year in years_list:
            self.gen_data_target[f"capex_{year}"] = costs[f"capex_{year}"]
        for year in years_list:
            self.gen_data_target[f"fuel_costs_{year}"] = costs[f"fuel_costs_{year}"]
        for year in years_list:
            self.gen_data_target[f"fixed_ops_{year}"] = costs[f"fixed_ops_{year}"]
        for year in years_list:
            self.gen_data_target[f"var_ops_{year}"] = costs[f"var_ops_{year}"]

        # Explicitly modifying the original DataFrame rather than a
        # copy to avoid "FutureWarning"
        # gen_data_target["V Setpoint p.u."].fillna(1, inplace=True)
        # gen_data_target["Output_pct_0"].fillna(0.6, inplace=True)
        # gen_data_target["Output_pct_1"].fillna(1,inplace=True)
        self.gen_data_target["V Setpoint p.u."] = pd.to_numeric(
            self.gen_data_target["V Setpoint p.u."], errors="coerce"
        ).fillna(1)
        self.gen_data_target["Output_pct_0"] = pd.to_numeric(
            self.gen_data_target["Output_pct_0"], errors="coerce"
        ).fillna(0.6)
        self.gen_data_target["Output_pct_1"] = pd.to_numeric(
            self.gen_data_target["Output_pct_1"], errors="coerce"
        ).fillna(1)

        self.gen_data_target.fillna("", inplace=True)
        self.gen_data_target.head()

        if save_csv:
            self.gen_data_target.to_csv(
                "data/costs/candidate_generators_initial_list.csv", index=False
            )
