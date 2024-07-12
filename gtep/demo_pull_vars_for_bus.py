import pandas as pd
from pathlib import Path

# Net Capacity Factor || CAPEX         || Construction Financing Cost || Overnight Capital Cost || fixed operation etc. || variable operation etc. || fuel costs($/MWh) (fuel costs won't exist for all types, can be assumed 0 in that case)
# Summart_CF          || Summary_CAPEX || Scrape                      || Scrape                 || Scrape               || Scrape                  || Summary_Fuel

# pull the bus gen mappigs
bus_data_df = pd.read_csv("Bus_data_gen_mappings.csv")

# read pricing data
pricing_data_book_names = [
    "Land-Based Wind",
    "Solar - CSP",
    "Natural Gas_FE",
    "Coal_FE",
    "Biopower",
    "Nuclear",
    "Geothermal",
    "Solar - Utility PV",
]
pricing_data_sheet = "./gtep/2022 v3 Annual Technology Baseline Workbook Mid-year update 2-15-2023 Clean.xlsx"
pricing_data_sheet_path = Path(pricing_data_sheet)
pricing_dict = {}
for this_book_name in pricing_data_book_names:
    pricing_dict[this_book_name] = pd.read_excel(
        pricing_data_sheet_path, this_book_name
    )

# now show a demo of how you can use the bus data df to pull data of interest

# example: get Land-Basd Wind CAPEX for all counties using the Moderate value for all years
# CAPEX ($/kW)

# define the gen type
gen_type = "Land-Based Wind"
var_of_interest = "CAPEX ($/kW)"
subvar_of_interest = "Moderate"

demo_df = pricing_dict[gen_type][
    (pricing_dict[gen_type]["Key1"] == var_of_interest)
    & (pricing_dict[gen_type]["Key3"] == subvar_of_interest)
]

# bus_data_df['Land-Based Wind']

bus_cutdown_df = bus_data_df[["Bus Name", gen_type]]
the_thing_we_want = bus_cutdown_df["Land-Based Wind"].iloc[0]
the_row_that_has_the_prices = demo_df[demo_df["Key2"] == the_thing_we_want]
bus_capex_df = pd.merge(bus_cutdown_df, demo_df, left_on=gen_type, right_on="Key2")

pass
