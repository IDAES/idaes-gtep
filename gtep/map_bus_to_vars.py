import pandas as pd
from pathlib import Path

from geopy.geocoders import Nominatim

# read bus data
bus_data_df = pd.read_csv("Bus_data.csv")

# if you need find the counties given the busses, run below, otherwise it is slow
"""
geolocator = Nominatim(user_agent="GetLoc")
county_lookup = []
for this_lat, this_lon in zip(
    bus_data_df["Bus latitude"], bus_data_df["Bus longitude"]
):
    tmp_county = geolocator.reverse([this_lat, this_lon])
    county_lookup.append(
        tmp_county.raw["address"]["county"][:-7]
    )  # shove in just the county name, not the "county" suffix

bus_data_df["County"] = county_lookup
bus_data_df.to_csv("Bus_data_extended.csv", index=False)
"""

bus_data_df = pd.read_csv("Bus_data_extended.csv")

# read the mapping tables (C.1 and C.2)
c1_data_sheet = (
    "./gtep/LTSA_ResourceCitingMethodology_CountyRatingsBasedOnResourceType.xlsx"
)
c1_data_sheet_path = Path(c1_data_sheet)
c2_data_sheet = (
    "./gtep/LTSA_ResourceCitingMethodology_SuitabilityofCountyforTechnologyType.xlsx"
)
c2_data_sheet_path = Path(c2_data_sheet)

c1_df = pd.read_excel(c1_data_sheet_path, keep_default_na=False, na_values=["_"])
c2_df = pd.read_excel(c2_data_sheet_path, keep_default_na=False, na_values=["_"])

# summary of mappings
mappings_dict = {
    "Land-Based Wind": {
        "Table": "C.1",
        "Column": "WindCond.",
        1: "Land-Based Wind - Class 8",
        2: "Land-Based Wind - Class 7",
        3: "Land-Based Wind - Class 6",
        4: "Land-Based Wind - Class 5",
    },
    "Solar - CSP": {
        "Table": "C.1",
        "Column": "SolarThermalCond.",
        7: "CSP - Class 2",
        6.75: "CSP - Class 2",
        6.5: "CSP - Class 3",
        6: "CSP - Class 7",
        5.75: "CSP - Class 7",
        5.5: None,
    },
    "Natural Gas_CT": {
        "Table": "C.2",
        "Column": "NatGasCT",
        "Good": "NG F-Frame CT",
        "Avg.": "NG F-Frame CT",
        "NA": None,
        "RO": "NG F-Frame CT",
    },
    "Natural Gas_FE": {
        "Table": "C.2",
        "Column": "NatGasCC",
        "Good": "NG Combined Cycle (H-Frame) Max CCS",
        "Avg.": "NG Combined Cycle (H-Frame) 95% CCS",
        "NA": None,
        "RO": None,
    },
    "Coal_FE": {
        "Table": "C.2",
        "Column": "Coal",
        "Good": None,
        "Avg.": None,
        "NA": None,
        "RO": None,
    },
    "Biopower": {
        "Table": "C.2",
        "Column": "Biomass",
        "Good": "Biopower - Dedicated",
        "Avg.": None,
        "NA": None,
        "RO": None,
    },
    "Nuclear": {
        "Table": "C.2",
        "Column": "Nuclear",
        "Good": "Nuclear - Small Modular Reactor",
        "Avg.": None,
        "NA": None,
        "RO": None,
    },
    "Geothermal": {
        "Table": "C.1",
        "Column": "Geothermal",
        "H": "Geothermal - NF EGS / Binary",
        "M": None,
        "NA": None,
    },
    "Solar - Utility PV": {
        "Table": "C.1",
        "Column": "SolarPV",
        "VH": "Utility PV - Class 1",
        "H": "Utility PV - Class 2",
        "M": "Utility PV - Class 4",
    },
}

# allocate for the mappings
allocation_dict = {
    item: [
        None for _ in range(bus_data_df["County"].shape[0])
    ]  # preallocate a bunch of lists with Nones for some amount of safety
    for item in mappings_dict
}

# allocate some debug stuff (mostly checking if we have all the counties)
missing_counties = []

# now run through and assign each bus the appropriate data
for county_ix, this_county in enumerate(bus_data_df["County"]):
    # clean up this county info by removing spaces
    clean_this_county = this_county.replace(" ", "")

    # check that the county is even in the tables
    if (clean_this_county in c1_df["County"].values) and (
        clean_this_county in c1_df["County"].values
    ):
        # cycle through the generation types
        for this_gen, this_map in mappings_dict.items():
            mapping_row = None
            # pick the correct mapping dictionary
            if this_map["Table"] == "C.1":
                mapping_row = c1_df[c1_df["County"] == clean_this_county]
            if this_map["Table"] == "C.2":
                if clean_this_county in c1_df["County"].values:
                    mapping_row = c2_df[c2_df["County"] == clean_this_county]

            county_quality = mapping_row[this_map["Column"]].values[0]
            workbook_map = this_map[county_quality]
            allocation_dict[this_gen][county_ix] = workbook_map
    else:
        print(
            f"[WARNING]: County {clean_this_county} not found in either C.1 or C.2. Skipping. Values will be 'None' in dataframe."
        )
        missing_counties.append(clean_this_county)

# assign new data to the old df
for gen_type, allocations in allocation_dict.items():
    bus_data_df[gen_type] = allocations

# save the new DF
bus_data_df.to_csv("Bus_data_gen_mappings.csv", index=False)

# print some diagnostic info
print(f"Missing the following counties: {missing_counties}")

pass

# mappings template:
# mappings_dict = {
#     "Land-Based Wind": {"Source": None, "Good": None, "Avg": None, "NA": None, "RO": None},
#     "Solar - CSP": {"Source": None, "Good": None, "Avg": None, "NA": None, "RO": None},
#     "Natural Gas_FE": {"Source": None, "Good": None, "Avg": None, "NA": None, "RO": None},
#     "Coal_FE": {"Source": None, "Good": None, "Avg": None, "NA": None, "RO": None},
#     "Biopower": {"Source": None, "Good": None, "Avg": None, "NA": None, "RO": None},
#     "Nuclear": {"Source": None, "Good": None, "Avg": None, "NA": None, "RO": None},
#     "Geothermal": {"Source": None, "Good": None, "Avg": None, "NA": None, "RO": None},
#     "Solar - Utility PV": {"Source": None, "Good": None, "Avg": None, "NA": None, "RO": None},
# }

# PopPop the Power Optimization Possum says ౿ᓕ ̤Ꜥ·> --- "when at first you don't succeed, scream."
# --- Land-based Wind ---
# ***************************************************************************************************************************************************
# - C.1 Rating      Workbook Classification
# - Condition 1     Class 8
# - Condition 2     Class 7
# - Condition 3     Class 6
# - Condition 4     Class 5
# ***************************************************************************************************************************************************

# --- Solar-CSP ---
# ***************************************************************************************************************************************************
# - C.1 Rating              Workbook Classification
# - Good (7)                Class 2
# - Good (6.75)             Class 2
# - Average (6.5)           Class 3
# - Below Average (6)       Class 7
# - Below Average (5.75)    Class 7
# - Poor (5.5               None
# ***************************************************************************************************************************************************

# --- NatGasCT ---
# ***************************************************************************************************************************************************
# - C.2 Rating      Workbook Classification
# - Good/Avg/RO     NG F-Frame CT
# - NA              None
# ***************************************************************************************************************************************************

# --- NatGasCC ---
# ***************************************************************************************************************************************************
# - C.2 Rating      Workbook Classification
# - Good            NG Combined Cycle (H-Frame) Max CCS
# - Avg             H-FrameCC95CCS
# - RO              None
# - NA              None
# ***************************************************************************************************************************************************

# --- Coal_FE ---
# ***************************************************************************************************************************************************
# - C.2 Rating      Workbook Classification
# - Good            None (ignore new, newTo2ndGen, newToTT, CCS95, CCS95To2ndGen, CCS95ToTT, MaxCCS, MaxCCSTo2ndGen, MaxCCSToTT, IGCC)
# - NA              None
# ***************************************************************************************************************************************************

# --- Biopower ---
# ***************************************************************************************************************************************************
# - C.2 Rating      Workbook Classification
# - Good            Biopower - Dedicated
# - NA              None
# ***************************************************************************************************************************************************

# --- Nuclear ---
# ***************************************************************************************************************************************************
# - C.2 Rating      Workbook Classification
# - Good            Nuclear - Small Modular Reactor (pretend AP1000 isn't a thing)
# - NA              None
# ***************************************************************************************************************************************************

# --- Geothermal ---
# ***************************************************************************************************************************************************
# - C.1 Rating      Workbook Classification
# - High            Geothermal - NF EGS / Binary (Ignore HydroFlash, HydroBinary, NFEGSFlash)
# - NA              None
# ***************************************************************************************************************************************************

# --- Solar - Utility PV ---
# ***************************************************************************************************************************************************
# - C.1 Rating      Workbook Classification
# - Very High       Class 1
# - High            Class 2
# - Medium          Class 4
# ***************************************************************************************************************************************************


# tech_baseline_to_suitability_map

# Wind            || SolarThermal || NatGasCT       || NatGasCC       || Coal    || Biomass  || Nuclear || GeoThermal || SolarPV
# Land-Based Wind || Solar - CSP  || Natural Gas_FE || Natural Gas_FE || Coal_FE || Biopower || Nuclear || Geothermal || Solar - Utility PV

# Net Capacity Factor || CAPEX         || Construction Financing Cost || Overnight Capital Cost || fixed operation etc. || variable operation etc. || fuel costs($/MWh) (fuel costs won't exist for all types, can be assumed 0 in that case)
# Summart_CF          || Summary_CAPEX || Scrape                      || Scrape                 || Scrape               || Scrape                  || Summary_Fuel


# Table A.1 Gives us SOME hints on how we can map things, or how things from C.1 and C.2 depend on each other.
# Everything outside of Solar PV requires "Low Urban Density", so we're just going to completely ignore that entirely.
# Wind needs to map to Wind Zones
# Solar Thermal need Solar Resource and Water Vailability
# NG CT needs NG Supply, and Non-attainment is requird if replacing higher-emission generation
# NG CC needs NG Supply, and Water Availability
# Coal needs Rail, Can't be in  Non-Attainment Area, and Water Availability
# Biomass needs Rail and Biomass
# Nuclear needs Water Availability
# Geothermal needs Geothermal Potential
# Solar PV needs Solar Resource

# mapping Table                     Table C.1: County Ratings Based on Resource Type -------> Table C.2: Suitability of County for Technology Type*
#
# Resource                          C.1 Label                       C.1 Values                      C.2 Label                   C.2 Values
# Wind                              Wind Condition                  1, 2, 3, 4                      Wind                        NA (one Avg), Avg (one NA), Good, Good
# SolarThermal                      SolarThermalCond.               5.5, 5.75, 6, 6.5, 6.75, 7      SolarThermal                NA, NA, NA, NA or Avg, NA or Good, NA
# NatGasCT                          GasPipelineDensityCapacity      H, M, L, VL                     NatGasCT                    (Mixed between Good, Avg, RO, and NA; mostly Good), Avg. (one Good, one NA), NA, NA
# NatGasCC                          GasPipelineDensityCapacity      H, M, L, VL                     NatGasCC                    (Mixed between Good, Avg, RO, and NA; mostly Good), Avg. (one Good, some NA), NA, NA
# NatGasCC                          Surface Water                   H, M, L, VL                     NatGasCC                    (Mixed between Good, Avg, RO, and NA; mostly Good), Avg. (one Good, some NA), NA, NA
# Coal (Railroad Density)           RailroadDensity                 H, M, L, VL                     Coal                        (Mix of Good and NA), (Mix of Good and NA), NA, NA
# Biomass                           Biomass                         H, M, L                         Biomass                     (Mix of Good and NA), NA, NA
# Geothermal                        Geothermal                      H, M, NA                        Geothermal                  Good (one NA), NA (one Good), NA
# SolarPV                           SolarPV                         VH, H, M                        SolarPV                     Good, Good, Avg

# Workbook identified variables
# Generation Type (Workbook)    Generation Type (Table C.1)     Generation Type (Table C.2)     Resource Classification (Workbook)                                                                          Resource Rating (Table C.1)         Technology Rating (Table C.2)
# Land-based Wind               Wind Cond.                      Wind                            Class 1-10 (lower is better)                                                                                Poor (1) - Very Good (7)            Good, Avg., NA
# Solar-CSP                     Solar Thermal Condition         Solar Thermal                   Class 2, 3, or 7 (lower is better)                                                                          5.5, 5.75, 6, 6.5, 6.75, 7          Good, Avg., NA
# Natural Gas_FE                Gas Pipeline Density/Capacity   NatGasCT                        NG F-Frame CT                                                                                                   H, M, L, VL                         Good, Avg., RO, NA
# Natural Gas_FE                Gas Pipeline Density/Capacity   NatGasCC                        G-FrameCC, H-FrameCC, F-FrameCC95CCS, H-FrameCC95CCS, F-FrameCCMaxCCS, NG Combined Cycle (H-Frame) Max CCS                      H, M, L, VL                         Good, Avg., RO, NA
# Natural Gas_FE                Surface Water                   NatGasCC                        G-FrameCC, H-FrameCC, F-FrameCC95CCS, H-FrameCC95CCS, F-FrameCCMaxCCS, NG Combined Cycle (H-Frame) Max CCS                      H, M, L, VL                         Good, Avg., RO, NA
# Coal_FE                       Railroad Density                Coal                            new, newTo2ndGen, newToTT, CCS95, CCS95To2ndGen, CCS95ToTT, MaxCCS, MaxCCSTo2ndGen, MaxCCSToTT, IGCC        H, M, L, VL                         Good, NA
# Biopower                      Biomass                         Biomass                         Biopower - Dedicated                                                                                                   H, M, L                             Good, NA
# Nuclear                       Surface Water                   Nuclear                         AP1000, Nuclear - Small Modular Reactor                                                                                 H, M, L, VL                         Good, NA
# Geothermal                    Geothermal                      GeoThermal                      HydroFlash, HydroBinary, NFEGSFlash, Geothermal - NF EGS / Binary, DeepEGSFlash, DeepEGSBinary                               H, M, NA                            Good, NA
# Solar - Utility PV            Solar PV                        SolarPV                         Class 1-10 (lower is better)                                                                                VH, H, M                            Good, Avg.

# need to combine ratings to say is the RESOURCE AVAILABLE and is the TECHNOLOGY SUITABLE and then map that to the PRICE of THINGS for that County
# It's clear from the PDF that C.2 should be used if possible: "Table C.2 shows a summary of [A.1 constraints and C.1 favorable counties] classification and is used in determining potential locations for generator siting."

# --- Land-based Wind ---
# Resource Classification (Workbook) have a measurement of wind speed, which helps. Workbook: 1.7 m/s (Class10) to >9.0 m/s (Class1)
# Resource Rating also has a class rating, BUT it's reversed: "zone 1" is the worst, "zone 7" is the best. Table C.1: Condition Rating (CR) 1-7
# Mapping the counties together we can align the Condition Rating with the Technology Rating and get NA = 1, "Avg" == 2, "Good" == 3 or 4.
# For sake of argument, Let's just assume that CR1 --> Class10, CR2 --> Class9, etc.
# Finally found a wind map: https://www.nrel.gov/grid/assets/images/wtk-100-north-america-50-nm-with-data-extents-01.jpg
# Of course the bins don't match anything else
# East texas is 6-6.9 m/s. That's either Class 9 or Class 8. Should map to "Wind Condition 1".
# 7-7.9 maps to Wind Condition 2 (ish), which maps to Class 7 or Class 6
# darkest blue on the map is 8.0-8.9, which is clearly Wind Condition 4 (there's not higher level on the map really), which should be Class 5-2. OOF.
# That leaves Wind Condition 3, which means we should put Class 6 here. Which means:
# Condition 4 --> Good  --> Class 5
# Condition 3 --> Good --> Class 6
# Condition 2 --> Avg --> Class 7
# Condition 1 --> NA --> Class 8
# ***************************************************************************************************************************************************
# - C.1 Rating      Workbook Classification
# - Condition 1     Class 8
# - Condition 2     Class 7
# - Condition 3     Class 6
# - Condition 4     Class 5
# ***************************************************************************************************************************************************


# --- Solar-CSP ---
# The Resource Rating here is actually decent, and is measured in kWh/m^2/day. This is actually useful. Praise be.
# - Good = 6.5-7 kWh/m^2/day
# - Average = 6-6.5 kWh/m^2/day
# - Below Average = 5.5-6 kWh/m^2/day
# - Poor = <=5.5 kWh/m^2/day
# The Workbook only has 3 options, Class 2, 3 or 7 (lower is better) AND THEY HAVE UNITS.
# - Class 7 = 6.16 kWh/m^2/day
# - Class 3 = 7.1 kWh/m^2/day
# - Class 2 = 7.46 kWh/m^2/day
# which means we can map these:
# - Good = Class 3 or Class 2
# - Average = Class 7
# - Below Average = NONE
# - Poor  = NONE
# but this isn't super useful, so I recommend we do this:
# ***************************************************************************************************************************************************
# - C.1 Rating      Workbook Classification
# - Good            Class 2
# - Average         Class 3
# - Below Average   Class 7
# - Poor            None
# ***************************************************************************************************************************************************

# --- NatGasCT ---
# The difference between NatGasCT and NatGasCC:
# NG CT needs NG Supply, and Non-attainment is required if replacing higher-emission generation
# NG CC needs NG Supply, and Water Availability
# NatGasCT is easy: there is only one mapping you can do.
# ***************************************************************************************************************************************************
# - C.2 Rating      Workbook Classification
# - Good/Avg/RO     NG F-Frame CT
# - NA              None
# ***************************************************************************************************************************************************

# --- NatGasCC ---
# NatGasCC is hell. We're in hell.
# In order to solve the mapping we need to know how levels of Post Combustion Carbon Capture (CCS) maps to Surface Water and Gas Pipeline Density/Capacity.
# There is something between F-Frame (old) and H-Frame (new) that matters, but I don't know what requires what. It sounds like a style.
# CCS evolves upwards, either it's "normal" or it's 95% or it's "Max" (97%).
# How any of that maps to either Surface Water or to Pipeline Density is completely lost on me.
# Not solving now. Placeholder mapping
# F-FrameCC, H-FrameCC, F-FrameCC95CCS, H-FrameCC95CCS, F-FrameCCMaxCCS, NG Combined Cycle (H-Frame) Max CCS
# ***************************************************************************************************************************************************
# - C.2 Rating      Workbook Classification
# - Good            NG Combined Cycle (H-Frame) Max CCS
# - Avg             H-FrameCC95CCS
# - RO              None
# - NA              None
# ***************************************************************************************************************************************************

# --- Coal_FE ---
# Similar to NatGas CC, there's so many possible mappings, we need a lot more information.
# Not solving now. Placeholder mapping
# Heck coal, everything NA?
# Coal isn't real, it can't hurt you
# ***************************************************************************************************************************************************
# - C.2 Rating      Workbook Classification
# - Good            None (ignore new, newTo2ndGen, newToTT, CCS95, CCS95To2ndGen, CCS95ToTT, MaxCCS, MaxCCSTo2ndGen, MaxCCSToTT, IGCC)
# - NA              None
# ***************************************************************************************************************************************************

# --- Biopower ---
# Biopower is also easy, there's only one mapping. Either it's there or not.
# ***************************************************************************************************************************************************
# - C.2 Rating      Workbook Classification
# - Good            Biopower - Dedicated
# - NA              None
# ***************************************************************************************************************************************************

# --- Nuclear ---
# Nuclear also needs more information to be useful, but we can at least pick one.
# ***************************************************************************************************************************************************
# - C.2 Rating      Workbook Classification
# - Good            Nuclear - Small Modular Reactor (pretend AP1000 isn't a thing)
# - NA              None
# ***************************************************************************************************************************************************

# --- Geothermal ---
# Once again, we need a lot more info, but there's some things to note:
# Deep EGS (DeepEGSFlash and DeepEGSBinary) both have no possible sites. We can probably safely ignore those.
# Beyond that, we have Hydrothermal and NF-EGS. Among those, there are Flash (>=200C) and Binary (<200C)
# The PDF has a map showing three different zones for Geothermal: Hydrothermal, Geopressure, and Hot Dry Rock.
# Hydrothermal is probably Hydrothermal in both cases. That only really makes sense for Binary, because the map says that the highest temperature is
# ~160F (71C), which is well below the 200C needed for Flash. Geopressure can get high enough temperature, but is deep.
# Geopressure is very deep (13000 ft ~= 3.9 km), so this likely applies to the Deep EGS (3-6 km). Let's ignore it for now.
# Hot Dry Rock is very explicity "little or no water", so let's seperate Hydro from NFEGS
# C.2 only gives us High, Medium and NA. Everything is pain.
# Looking at the workbook, Flash and Binary appear to be identical in the data and we can ignore Hydro and NFEGS, until you get to costs.
# If it was the other way around, we could instead map "Medium" and "High" to production.
# However, one las thing: comparing C.1 and C.2 shows us that only "High" get "Good" ratings outside of one, and "Medium" get "NA" ratings outside of
# one. So let's just assign "Good" to something and NA to everything None
# ***************************************************************************************************************************************************
# - C.2 Rating      Workbook Classification
# - High            Geothermal - NF EGS / Binary (Ignore HydroFlash, HydroBinary, NFEGSFlash)
# - NA              None
# ***************************************************************************************************************************************************

# --- Solar - Utility PV ---
# For some reason, we can't assign any actual units here, even though we did so on the other solar classificaitons. Cool. Thanks. Appreciate that.
# Instead we get "Very High", "High", and "Medium" from C.1 and "Good" or "Avg." from C.2.
# Mapping C.1 to C.2 we get that Very High = Good, High = Good and Medium = Avg. So let's use C.1 because it's more detailed.
# The workbook gives us Class 1-10. Lower is better. Everything is pain.
# BUT I FOUND A REFERENCE:
# https://atb.nrel.gov/electricity/2022/utility-scale_pv
# Resource Class	GHI Bin (kWh/m2/day)	Mean AC Capacity Factor	    Area (sq. km)
# 1	                >5.75	                32.80%	                    216,551
# 2	                5.5–5.75	            31.80%	                    349,894
# 3	                5.25–5.5	            30.30%	                    372,764
# 4	                5–5.25	                28.70%	                    497,444
# 5	                4.75–5	                26.80%	                    779,720
# 6	                4.5–4.75	            25.80%	                    870,218
# 7	                4.25–4.5	            24.60%	                    727,918
# 8	                4–4.25	                23.40%	                    828,438
# 9	                3.75–4	                22.30%	                    794,496
# 10	            <3.75	                20.40%	                    163,120
#
# Notably, Not close to the other Solar classification with actual units, so good thing I didn't use that *teehee swings legs*
# I also found an image:
# https://atb.nrel.gov/img/electricity/2021/p19/v1/solar-annual-ghi-2018-usa-scale-01.png
# So this shows that Texas actually spans 5 different Classes. I will not be back deriving them here. Importantly, it spans Class 1-5, going west to
# east, which means we can probably just figure it out. Using the weather zones:
#
# - FarWest gets Very High and High ratings. --> Class 1
# - West is all High --> Class 2
# - NorthCentral is all High --> Class 2
# - North is all High --> Class 2
# - East has mostly High and some Medium --> Class 3
# - SouthCentral has mostly High with some Medium --> Class 3
# - South has mostly Medium with some Highs --> Class 4
# - Coast is all Medium  --> Class 5
#
# Clearly the mapping above won't actually work (they don't even correspond to the image I made them from) but instead we can use it to try to fake it
# So let's just throw some classifications out there and hope.
# ***************************************************************************************************************************************************
# - C.1 Rating      Workbook Classification
# - Very High       Class 1
# - High            Class 2
# - Medium          Class 4
# ***************************************************************************************************************************************************


# Table A.1 Generation Technologies and Resource Limitations from PDF:
# Resource          Type Rating     Resource Limitations
# Wind	            Good	        Low Urban Density; Wind Zone 3-4
# Wind	            Average	        Low Urban Density; Wind Zone 2
# Solar Thermal	    Good	        Low Urban Density; Good Direct Solar Resource; High Water Availability
# Solar Thermal	    Average	        Low Urban Density; Average to Good Direct Solar Resource; Medium to High Water Availability
# NG CT	            Good	        Low Urban Density; High Availability of Natural Gas Supply; Can Only Build in Non-Attainment Area if Replacing Higher-Emission Generation
# NG CT	            Average	        Low Urban Density; Medium Availability of Natural Gas Supply; Can Only Build in Non-Attainment Area if Replacing Higher-Emission Generation
# NG CC	            Good	        Low Urban Density; High Availability of Natural Gas Supply; Can Only Build in Non-Attainment Area if Replacing Higher-Emission Generation; Medium to High Water Availability
# NG CC	            Average	        Low Urban Density; Medium Availability of Natural Gas Supply; Can Only Build in Non-Attainment Area if Replacing Higher-Emission Generation; Medium to High Water Availability
# Coal	            Good	        Low Urban Density; Medium to High Availability of Rail Transportation; Cannot Build in Non-Attainment Area; Medium to High Water Availability
# Biomass	        Good	        Low Urban Density; Low to High Availability of Rail Transportation; High Biomass Availability
# Nuclear	        Good	        Low Urban Density; High Water Availability
# Geothermal	    Good	        Low Urban Density; High Geothermal Potential
# Solar PV	        Good	        High to Very High Total Solar Resource
# Solar PV	        Average	        Medium Total Solar Resource


pass
