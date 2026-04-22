import os
import numpy as np
import pandas as pd

"""This script generates data files formatted for use with the
IDAES-GTEP model.

"""


def main(mfile):

    # Add name of .m file
    data_file = mfile
    
    # Specify the path to your .m file
    path = os.path.dirname(os.path.realpath(__file__))
    data_path = os.path.join(path, f"{data_file}.m")

    # Read the data from the .m file
    bus_data, info_text = read_data_from_m_file(data_path, "bus")
    create_and_save_data(data_file, data_path, bus_data, "bus", info_text=info_text)

    gen_data, info_text = read_data_from_m_file(data_path, "gen")
    create_and_save_data(data_file, data_path, gen_data, "gen", info_text=info_text)

    branch_data, info_text = read_data_from_m_file(data_path, "branch")
    create_and_save_data(data_file, data_path, branch_data, "branch", info_text=info_text)

    # Create other files
    csv_file_name = "reserves.csv"
    create_reserve_product(path, data_file, csv_file_name)

    csv_filename = "simulation_parameters.csv"
    create_simulation_parameters(path, data_file, csv_filename)

    create_pointers(path, data_file)

    
def create_and_save_data(data_file, data_path, data, data_type, info_text):
    """This method processes and saves the data from a specified type
    (bus, generator, generator cost, or branch) into a CSV file.  This
    method reads the data, organizes it into a structured format, and
    then writes it to a CSV file in a directory named after the input
    .m file.

    Input Parameters:
    - data_path (str): Path to the .m data file
    - data (list of lists): Contains the data to be processed, where
                            each inner list represents a row of data
                            for the specified type ('bus', 'gen', 'or
                            'branch')
    - data_type (str): A string indicating the type of data being
                       processed and it can be 'bus', 'gen', or
                       'branch
    - info_text (str): Optional text containing generator setpoint 
                       replacement notes.

    """


    # Define the headers for the .csv file based on the data type
    if data_type == "bus":
        headers = [
            "Bus ID",
            "Bus Name",
            "BaseKV",
            "Bus Type",
            "MW Load",
            "MVAR Load",
            "V Mag",
            "V Angle",
            "Area",
            "Zone",
        ]
    elif data_type == "gen":
        headers = [
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
    elif data_type == "branch":
        headers = [
            "UID",
            "From Bus",
            "To Bus",
            "R",
            "X",
            "B",
            "Cont Rating",
            "LTE Rating",
            "STE Rating",
            "Tr Ratio",
        ]
    else:
        raise ValueError("Unknown data type")

    # Prepare the data for the DataFrame
    df_data = []

    # Extract generator setpoints if info_text is provided
    generator_setpoints = extract_generator_setpoints(info_text)

    for row in data:
        if data_type == "bus":
            if row[1] == 2:
                val = "PV"
            elif row[1] == 3:
                val = "PQ"
            else:
                val = "Ref"


            bus_id = int(row[0])
            bus_name = f"bus{int(bus_id)}"  # Creating a Bus Name based on Bus ID
            base_kv = row[9]
            bus_type = val
            mw_load = row[2]
            mvar_load = row[3]
            v_mag = row[7]
            v_angle = row[8]
            area = row[6]
            zone = row[10] if len(row) > 10 else ""

            columns = [
                bus_id,
                bus_name,
                base_kv,
                bus_type,
                mw_load,
                mvar_load,
                v_mag,
                v_angle,
                area,
                zone,
            ]

        elif data_type == "gen":
            pg = row[1]
            qg = row[2]
            gen_uid = row[0]
            bus_id = ""
            unit_type = ""
            fuel = ""
            mw_inj = 0
            mvar_inj = 0
            v_setpoint = row[5]
            pmax = row[8]
            pmin = row[9]
            qmax = row[3]
            qmin = row[4]
            min_down_time = ""
            min_up_time = ""
            ramp_rate = ""
            start_time_cold = ""
            start_time_warm = ""
            start_time_hot = ""
            start_heat_cold = ""
            start_heat_warm = ""
            start_heat_hot = ""
            non_fuel_start_cost = ""
            fuel_price = ""
            output_pct_0 = ""
            output_pct_1 = ""
            output_pct_2 = ""
            output_pct_3 = ""
            output_pct_4 = ""
            hr_avg = ""
            hr_incr_1 = ""
            hr_incr_2 = ""
            hr_incr_3 = ""
            hr_incr_4 = ""

            # Iterate through the dictionary generator_setpoints to
            # add bus_id data. To do this, check if in_text_pmax
            # matches pmax and if it does, save the bus_id
            for gen_info in generator_setpoints.values():
                if gen_info['in_text_pg'] == pg:
                    bus_id = gen_info['bus_id'] 
                    break  # Exit the loop once a match is found

            # Define operational limits and fuel types for each
            # generator type. NOTE: These limits are assumed for the
            # purpose of data generation.
            other_generator_info = {
                "HYDRO": {   # Hydroelectric
                    "limits": (0, 20),
                    "fuel": "H"
                },
                "WIND": {  # Wind
                    "limits": (20, 50),
                    "fuel": "W"
                },
                "PV": {  # Photovoltaic
                    "limits": (51, 100),
                    "fuel": "S"
                },
                "CT": {  # Combustion Turbine (Gas)
                    "limits": (101, 300),
                    "fuel": "G"
                },
                "RTPV": {  # Renewable Photovoltaic (Solar)
                    "limits": (301, 600),
                    "fuel": "S"
                },
                "STEAM": {  # Steam Turbine (Coal)
                    "limits": (601, 5000),
                    "fuel": "C"
                },
                "CC": {  # Combined Cycle (Gas)
                    "limits": (101, 300),
                    "fuel": "G"
                }
            }
            # operational_limits = {
            #     "HYDRO": (0, 20),
            #     "WIND": (20, 50),
            #     "PV": (51, 100),
            #     "CT": (101, 300),
            #     "RTPV": (301, 600),
            #     "STEAM": (601, 5e3),
            # }

            for generator_type, info in other_generator_info.items():
                min_limit, max_limit = info["limits"]  # Access limits from the inner dictionary
                if min_limit <= pmax <= max_limit:
                    unit_type = generator_type
                    fuel = info["fuel"]

            # Update gen_uid to include the generator type
            gen_uid = f"{int(gen_uid)}_{unit_type}"

            columns = [
                gen_uid,
                bus_id,
                unit_type,
                fuel,
                mw_inj,
                mvar_inj,
                v_setpoint,
                pmax,
                pmin,
                qmax,
                qmin,
                min_down_time,
                min_up_time,
                ramp_rate,
                start_time_cold,
                start_time_warm,
                start_time_hot,
                start_heat_cold,
                start_heat_warm,
                start_heat_hot,
                non_fuel_start_cost,
                fuel_price,
                output_pct_0,
                output_pct_1,
                output_pct_2,
                output_pct_3,
                output_pct_4,
                hr_avg,
                hr_incr_1,
                hr_incr_2,
                hr_incr_3,
                hr_incr_4,
            ]

        elif data_type == "branch":      
            branch_uid = f"branch_{int(row[0])}_{int(row[1])}"
            from_bus = int(row[0])
            to_bus = int(row[1])
            r = row[2]
            x = row[3]
            b = row[4]
            cont_rating = row[5]
            lte_rating = row[6]
            stw_rating = row[7]
            tr_ratio = row[8]
            
            columns = [
                branch_uid,
                from_bus,
                to_bus,
                r,
                x,
                b,
                cont_rating,
                lte_rating,
                stw_rating,
                tr_ratio,
            ]

        df_data.append(columns)

    # Create the DataFrame
    df = pd.DataFrame(df_data, columns=headers)

    # Create a directory based on the .m file name (without extension)
    output_directory = os.path.join(data_file)
    if not os.path.exists(output_directory):
        os.makedirs(output_directory)
        print(f"Created directory: {output_directory}")
    else:
        print(f"Directory already exists: {output_directory}")

    # Write the DataFrame to a CSV file in the new directory
    csv_filename = os.path.join(output_directory, f"{data_type}.csv")
    df.to_csv(csv_filename, index=False)
    print(f"Successfully created CSV file: '{csv_filename}'")

    if data_type == "gen":
        
        # Create the custom CSV file with the specified structure
        initial_status_csv_filename = os.path.join(output_directory, 'initial_status.csv')
        create_initial_status_csv(initial_status_csv_filename, df_data)

def create_initial_status_csv(csv_filename, df_data):
    """
    Creates a custom CSV file with a specific structure based on the provided data.

    Parameters:
    - csv_filename (str): The path where the custom CSV file will be saved.
    - df_data (list of lists): Contains the data to be processed, where each inner list
                                represents a row of data for the specified type.
    """
    # Prepare the header and values
    header = []
    values1 = []
    values2 = []

    # Process the data to create the header and values
    for row in df_data:
        # Assuming the generator name is in the first column (gen_uid)
        generator_name = row[0]  # Adjust index as necessary
        p_values1 = row[7] * 0.1
        p_values2 = row[7] * 0

        # Create the header entry
        header.append(generator_name)

        # Create the values row
        values1.append(p_values1)
        values2.append(p_values2)

    # Create a DataFrame with the correct structure
    if values1:

        custom_df = pd.DataFrame(columns=header)
        custom_df.loc[0] = values1
        custom_df.loc[1] = values2
        
        # Write the custom DataFrame to a CSV file
        custom_df.to_csv(csv_filename, index=False, header=True)
        print(f"Successfully created CSV file: '{csv_filename}'")
    else:
        print("No data available to create the custom CSV file.")
        
def read_data_from_m_file(file_path, data_type):
    """This method reads specified data from a .m file based on the
       provided data type (bus, generator, or branch) and extracts it
       into a list of lists. Each inner list represents a row of data
       with numerical values.

    Input Parameters:
    - file_path (str): The path to the .m file to be read.
    - data_type (str): The type of data to extract. The valid options are 
                       'bus', 'gen', and 'branch'.

    Returns:
    - data (list of lists): A list containing the extracted data, where 
                             each inner list corresponds to a row of data 
                             from the specified section of the .m file.
    - info_text (str): The extracted generator setpoint replacement notes.

    """
    
    # Read original data file
    with open(file_path, "r") as file:
        lines = file.readlines()

    # Initialize the data list based on the data type
    data = []
    data_section = False
    info_text = ""
    info_section = False
    
    # Determine the section to look for based on the data type
    section_identifier = {
        "bus": "mpc.bus =",
        "gen": "mpc.gen =",
        "gencost": "mpc.gencost =",
        "branch": "mpc.branch =",
    }

    for line in lines:
        if section_identifier[data_type] in line:
            data_section = True
            continue
        # Remove text from data values
        if data_section:
            if line.strip() == "];":
                break
            if line.strip() and not line.startswith("%"):
                cleaned_line = line.strip().rstrip(";")
                data.append(list(map(float, cleaned_line.split())))

                
    # Collect information for generators from the .m file notes
    for line in lines:
        # Check for the generator setpoint replacement notes section
        if "=== Generator Setpoint Replacement Notes ===" in line:
            info_section = True
            continue
        if info_section:
            if line.strip() == "":
                info_section = False  # End of info section
            else:
                # Remove the "% INFO" prefix and append the line
                cleaned_line = line.replace('% INFO    : ', '').strip()
                info_text += cleaned_line + '\n'  # Append relevant info lines

    return data, info_text


def extract_generator_setpoints(info_text):
    """This method extracts the generator bus ID based on the text at the
    bottom of the .m file

    Parameters:
    - info_text (str): The text containing generator setpoint
                       replacement notes.

    Returns:
    - dict: A dictionary where keys are bus IDs and values are the
            updated Pg values.

    """
    generator_setpoints = {}
    lines = info_text.strip().split('\n')
    count_lines = 0
    for line in lines:
        count_lines += 1
        if "Gen at bus" in line:
            parts = line.split(':')
            bus_info = parts[0].strip()  # e.g., "Gen at bus n"
            pg_info = parts[1].strip()    # e.g., "Pg=0.0, Qg=0.0 -> Pg=0.0, Qg=0.0"

            if " Pg" in pg_info:
                # Extract bus ID from "Gen at bus X"
                bus_id = int(bus_info.split()[-1])

                # Extract the updated Pg value (the value before the '->')
                updated_pg = float(pg_info.split('->')[0].split(',')[0].split('=')[1].strip())

                # Extract new Pg value (the value after the '->') for matching with Pmax
                in_text_pg = float(pg_info.split('->')[1].split(',')[0].split('=')[1].strip())
                
                # Store the original Pg value and bus info in the dictionary
                generator_setpoints[f"bus_info_line{count_lines}"] = {
                    'in_text_pg': in_text_pg,
                    'updated_pg': updated_pg,
                    'bus_id': bus_id,
                }

    return generator_setpoints

def create_reserve_product(path, data_file, csv_filename, data=None):
    """Creates a CSV file with specified headers for reserve products.

    Parameters:
    - csv_filename (str): The path where the CSV file will be saved.
    - data (list of lists, optional): Contains the data to be added to the CSV file.
                                       Each inner list represents a row of data.
    """

    # Define the headers
    headers = ["Reserve Product", "Requirement (MW)"]

    # Create a DataFrame with the headers
    if data is None:
        # If no data is provided, create an empty DataFrame
        df = pd.DataFrame(columns=headers)
    else:
        # Create a DataFrame with the provided data
        df = pd.DataFrame(data, columns=headers)

    # Write the DataFrame to a CSV file
    csv_filepath = os.path.join(path, data_file, csv_filename)
    df.to_csv(csv_filepath, index=False, header=True)
    print(f"Successfully created CSV file: '{data_file}/{csv_filename}'")


def create_simulation_parameters(path, data_file, csv_filename):

    # Define the data as a list of dictionaries
    data = [
        {
            "Simulation_Parameters": "Periods_per_Step",
            "Description": "the number of discrete periods represented in each simulation step",
            "DAY_AHEAD": 24,
            "REAL_TIME": 1,
        },
        {
            "Simulation_Parameters": "Period_Resolution",
            "Description": "seconds per period",
            "DAY_AHEAD": 3600,
            "REAL_TIME": 300,
        },
        {
	    "Simulation_Parameters": "Date_From",
            "Description": "simulation beginning period",
            "DAY_AHEAD": "01/01/2020 0:00",
            "REAL_TIME": "01/01/2020 0:00",
        },
        {
            "Simulation_Parameters": "Date_To",
            "Description": "simulation ending period (must account for lookahed data availability)",
            "DAY_AHEAD": "12/31/2020 0:00",
            "REAL_TIME": "12/31/2020 0:00",
        },
        {
            "Simulation_Parameters": "Look_Ahead_Periods_per_Step",
            "Description": "the number of look ahead periods included in each optimization step",
            "DAY_AHEAD": 24,
            "REAL_TIME": 2,
        },
        {
            "Simulation_Parameters": "Look_Ahead_Resolution",
            "Description": "look-ahead period resolution",
            "DAY_AHEAD": 3600,
            "REAL_TIME": 300,
        },
        
    ]

    # Create a DataFrame from the data
    df = pd.DataFrame(data)

    # Write the DataFrame to a CSV file
    csv_filepath = os.path.join(path, data_file, csv_filename)
    df.to_csv(csv_filepath, index=False, header=True)
    print(f"Successfully created CSV file: '{data_file}/{csv_filename}'")

def create_pointers(path, data_file):
    # Define the filenames

    dam_load_filename = "DAY_AHEAD_load.csv"
    dam_renewables_filename = "DAY_AHEAD_renewables.csv"
    rtm_load_filename = "REAL_TIME_load.csv"
    rtm_renewables_filename = "REAL_TIME_renewables.csv"

    dam_renewables_file = os.path.join(path, data_file, dam_renewables_filename)
    dam_load_file = os.path.join(path, data_file, dam_load_filename)
    rtm_renewables_file = os.path.join(path, data_file, rtm_renewables_filename)
    rtm_load_file = os.path.join(path, data_file, rtm_load_filename)    
    output_file = "timeseries_pointers.csv"

    # Check if the files exist
    for i in [
            dam_renewables_file,
            dam_load_file,
            rtm_renewables_file,
            rtm_load_file,
    ]:
        if not os.path.exists(i):
            print(f"File not found: {i}. Please make sure to provide those files so the timeseries_pointers file can be generated.")
            return

    
    # Prepare the data for CSV file
    simulation_data = []

    # Read the renewables file and skip the first four columns
    for fname in [dam_renewables_file, rtm_renewables_file]:
        simulation_name = "_".join(os.path.splitext(os.path.basename(fname))[0].split("_")[:-1])
        
        renewables_df = pd.read_csv(fname)
        renewables_columns = renewables_df.columns[4:]

        # Extract the filename
        filename = os.path.basename(fname)
    
        for generator in renewables_columns:
            if generator in ["1_HYDRO", "2_RTPV"]:

                # Add PMIN MW row only to Hydro and RTPV generators
                simulation_data.append([simulation_name, "Generator", generator, "PMIN MW", filename])
    
            simulation_data.append([simulation_name, "Generator", generator, "PMAX MW", filename])

    for fname in [dam_load_file, rtm_load_file]:
        simulation_name = "_".join(os.path.splitext(os.path.basename(fname))[0].split("_")[:-1])
        
        # Read the load file and skip the first four columns 
        load_df = pd.read_csv(fname)
        load_columns = load_df.columns[4:]  

        # Extract the filename
        filename = os.path.basename(fname)

        # Add data from the load file
        for area in load_columns:
            simulation_data.append([simulation_name, "Area", area, "MW Load", filename])
        
    # Create a DataFrame from the simulation data and sort the data to
    # have "DAY_AHEAD" first and then the object type
    simulation_data.sort(key=lambda x: (x[0] != "DAY_AHEAD", x[1] != "Generator", x[3] != "PMIN MW"))
    simulation_df = pd.DataFrame(simulation_data, columns=["Simulation", "Category", "Object", "Parameter", "Data File"])

    # Write the DataFrame to a new CSV file
    csv_filename = "timeseries_pointers"
    csv_filepath = os.path.join(path, data_file, f"{csv_filename}.csv")
    simulation_df.to_csv(csv_filepath, index=False)
    print(f"Successfully created CSV file: '{data_file}/{csv_filename}'")

if __name__ == "__main__":

    main()


    
