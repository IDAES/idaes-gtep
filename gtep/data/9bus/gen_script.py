import pandas as pd
from IPython import embed

gen_mapping_1={'Conventional Hydroelectric':'hydud', 
               'Petroleum Liquids':'o-g-s', 
               'Natural Gas Fired Combustion Turbine':'gas-ct', 
               'Natural Gas Internal Combustion Engine':'gas-ct',
               'Natural Gas Fired Combined Cycle':'gas-cc', 
               'Natural Gas Steam Turbine':'o-g-s', 
               'Landfill Gas':'lfill-gas', 
               'Batteries':'gas-cc',#NEED 
               'Hydroelectric Pumped Storage':'hydud', #NEED
               'Geothermal':'geothermal',
               'Nuclear':'nuclear',
               'Onshore Wind Turbine':'wind-ons', 
               'Other Waste Biomass':'biopower',
               'Wood/Wood Waste Biomass':'biopower', 
               'Solar Photovoltaic':'distpv', 
               'Solar Thermal without Energy Storage':'distpv', 
               'All Other':'gas-cc', #NEED
               'Conventional Steam Coal':'coalolduns', 
               'Other Gases':'o-g-s', 
               'Petroleum Coke':'o-g-s', 
               'Municipal Solid Waste':'lfill-gas', 
               'Other Natural Gas':'o-g-s', 
               'IMPORT':'gas-cc',
               'Synchronous Condenser':'gas-cc'} #NEED
# Define the headers for the CSV
headers = [
    "GEN UID", "Bus ID", "Unit Type", "Fuel", "MW Inj", "MVAR Inj", 
    "V Setpoint p.u.", "PMax MW", "PMin MW", "QMax MVAR", "QMin MVAR", 
    "Min Down Time Hr", "Min Up Time Hr", "Ramp Rate MW/Min", 
    "Start Time Cold Hr", "Start Time Warm Hr", "Start Time Hot Hr", 
    "Start Heat Cold MBTU", "Start Heat Warm MBTU", "Start Heat Hot MBTU", 
    "Non Fuel Start Cost $", "Fuel Price $/MMBTU", 
    "Output_pct_0", "Output_pct_1", "Output_pct_2", "Output_pct_3", 
    "Output_pct_4", "HR_avg_0", "HR_incr_1", "HR_incr_2", 
    "HR_incr_3", "HR_incr_4"
]

# Define the viable technologies
viable_tech_list = [
    "gas-cc", "gas-ct", "coaloldscr", "coalolduns", "coal-new", 
    "biopower", "nuclear", "o-g-s", "coal-igcc", "lfill-gas", 
    "wind-ofs", "wind-ons", "dupv", "upv", "distpv", "hydud", 
    "hydund", "hydnpnd", "hydnd", "hydend", "hyded", "geothermal"
]

# Define the mapping dictionaries
unit_type_map = {
    "gas-cc": "CC", "gas-ct": "CT", "coaloldscr": "COAL", 
    "coalolduns": "COAL", "coal-new": "COAL", "biopower": "BIO", 
    "nuclear": "NUC", "o-g-s": "OGS", "coal-igcc": "COAL", 
    "lfill-gas": "LFILL", "wind-ofs": "WIND", "wind-ons": "WIND", 
    "dupv": "PV", "upv": "PV", "distpv": "PV", "hydud": "HYDRO", 
    "hydund": "HYDRO", "hydnpnd": "HYDRO", "hydnd": "HYDRO", 
    "hydend": "HYDRO", "hyded": "HYDRO", "geothermal": "GEO"
}

fuel_map = {
    "CC": "G", "BIO": "B", "NUC": "N", "OGS": "G", 
    "COAL": "C", "CT": "C", "LFILL": "G", "WIND": "W", 
    "HYDRO": "H", "PV": "S",'GEO':'GEO'
}
up_time_map = {
    "gas-cc": 6,
    "gas-ct": 1,#0
    "coaloldscr": 24,
    "coalolduns": 24,
    "coal-new": 24,
    "biopower": 7,
    "nuclear": 48,
    "o-g-s": 1,
    "coal-igcc": 24,
    "lfill-gas": 7,
    "geothermal": 1,#0
}
down_time_map = {
    "gas-cc": 8,
    "gas-ct": 1,#0
    "coaloldscr": 12,
    "coalolduns": 12,
    "coal-new": 12,
    "biopower": 7,
    "nuclear": 48,
    "o-g-s": 1,
    "coal-igcc": 12,
    "lfill-gas": 6,
    "geothermal": 1, #0
}
min_stable_level_map = {
    "gas-cc": 0.50,
    "gas-ct": 0.60,
    "coaloldscr": 0.40,
    "coalolduns": 0.40,
    "coal-new": 0.40,
    "biopower": 1,
    "nuclear": 1,
    "o-g-s": 1,
    "coal-igcc": 0.40,
    "lfill-gas": 1,
}
# Create a list to hold the rows of data
data = []

# Create 22 rows for each of the 3 regions
for region in range(1, 4):  # Regions 1 to 3
    for i in range(1, 23):  # 22 technologies
        tech = viable_tech_list[i - 1]  # Get the technology
        gen_uid = f"G{i}_R{region}"  # Generate GEN UID
        bus_id = region  # Keep Bus ID empty
        unit_type = unit_type_map[tech]  # Get Unit Type from mapping
        fuel = fuel_map[unit_type]  # Get Fuel from mapping

        min_up_time_hr = up_time_map.get(tech, 1)
        min_down_time_hr = down_time_map.get(tech, 1)
        p_min = min_stable_level_map.get(tech, 0)  #NEED A WAY TO LINK TO PG
        # Create a row with the required values
        row = [gen_uid, bus_id, unit_type, fuel] + [""] * (len(headers) - 4)
        row[headers.index("Min Up Time Hr")] = min_up_time_hr
        row[headers.index("Min Down Time Hr")] = min_down_time_hr
        data.append(row)



# Create a DataFrame from the data
df = pd.DataFrame(data, columns=headers)

#Constants
df['MVAR Inj'] = 0
df['V Setpoint p.u.'] = 1
df["QMax MVAR"] = 300
df["QMin MVAR"] = -300
df["Ramp Rate MW/Min"] = 1 
df["Non Fuel Start Cost $"] = 0  #NOT CORRECT
df["Fuel Price $/MMBTU"] = 1   #NOT CORRECT

#Datamap

cal_bus_to_zones=pd.read_csv('../GEN_DATA/california_bus_mapping_w_zone.csv')
cal_gens=pd.read_csv('../GEN_DATA/CATS_gens.csv')
df["MW Inj"] = 0
df["PMax MW"] = 0
df["PMin MW"] = 0
cal_gens['Technology'] = cal_gens['FuelType'].map(gen_mapping_1)
cal_bus_to_zones=cal_bus_to_zones.rename(columns={'bus_i':'bus'})
cal_gens = cal_gens.merge(cal_bus_to_zones, on='bus', how='left')
for index, row in cal_gens.iterrows():
    tech = row['Technology']
    zone = row['zone']

    if tech in viable_tech_list:

        tech_index = viable_tech_list.index(tech)
        
        if zone == 1:
            df.at[tech_index, "MW Inj"] +=row['Pg'] 
            df.at[tech_index, "PMax MW"] +=row['Pmax']
            df.at[tech_index, "PMin MW"] +=row['Pmin'] 
        elif zone == 2:
            df.at[tech_index+22, "MW Inj"] += row['Pg']
            df.at[tech_index+22, "PMax MW"] +=row['Pmax']
            df.at[tech_index+22, "PMin MW"] +=row['Pmin']  
        elif zone == 3:
            df.at[tech_index+22*2, "MW Inj"] +=row['Pg']
            df.at[tech_index+22*2, "PMax MW"] +=row['Pmax']
            df.at[tech_index+22*2, "PMin MW"] +=row['Pmin']  

#AZ Mapping
AZ_gens=pd.read_csv('../GEN_DATA/AZ_Gens.csv')
AZ_gens['Technology'] = AZ_gens['Technology'].map(gen_mapping_1)
for index, row in AZ_gens.iterrows():
    tech = row['Technology']
    zone = 3

    if tech in viable_tech_list:

        tech_index = viable_tech_list.index(tech)
        
        if zone == 1:
            df.at[tech_index, "MW Inj"] +=row['Pg'] 
            df.at[tech_index, "PMax MW"] +=row['Pmax']
            df.at[tech_index, "PMin MW"] +=row['Pmin'] 
        elif zone == 2:
            df.at[tech_index+22, "MW Inj"] += row['Pg']
            df.at[tech_index+22, "PMax MW"] +=row['Pmax']
            df.at[tech_index+22, "PMin MW"] +=row['Pmin']  
        elif zone == 3:
            df.at[tech_index+22*2, "MW Inj"] +=row['Pg']
            df.at[tech_index+22*2, "PMax MW"] +=row['Pmax']
            df.at[tech_index+22*2, "PMin MW"] +=row['Pmin']  
# Save the DataFrame to a CSV file
output_file = 'generated_units.csv'  # Update with your desired output file path
df.to_csv(output_file, index=False)

print(f"Saved the generated units to {output_file}")
