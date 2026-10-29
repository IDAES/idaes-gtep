# GTEP Data Directory

## Introduction

This directory contains the input data files required to build and run
the Generation and Transmission Expansion Planning (GTEP) model. The
data describes the power system network, generators, loads, storage
resources, renewable time series, and related operational and
investment parameters.

The purpose of this directory is to provide a consistent data
structure for converting input files into the model data object and,
ultimately, into Pyomo sets, parameters, variables, and constraints
used by GTEP.

## Required Data Files

The file `required_data_files.csv` summarizes the input files used by
the GTEP model and indicates whether each file is required, optional,
or conditionally required. These files are used to build the model data
object and initialize the GTEP model. They define the network topology,
generators, buses, storage resources, load and renewable time series,
reserve data, simulation settings, and supporting Prescient/GTEP
mappings. 

Below find a brief summary of the main input files:

| File | Brief Description |
|---|---|
| `bus.csv` | Bus and load data. |
| `branch.csv` | Transmission branch data. |
| `gen.csv` | Generator data. |
| `storage.csv` | Storage unit data. |
| `reserves.csv` | Reserve requirement data. |
| `timeseries_pointers.csv` | Time-series mapping file. |
| `simulation_objects.csv` | Simulation timing settings. |
| `DAY_AHEAD_load.csv` | Day-ahead load time series. |
| `DAY_AHEAD_renewables.csv` | Day-ahead renewable time series. |
| `REAL_TIME_load.csv` | Real-time load time series. |
| `REAL_TIME_renewables.csv` | Real-time renewable time series. |

## Data Dictionary

The file `gtep_data_dictionary.csv` documents the data fields used by
the GTEP workflow and serves as a reference for checking data
completeness and unit consistency before running the model. It lists
the fields expected in each input file, their units, where they are
stored in the parsed model data object, and which GTEP model
components use them.

The data dictionary includes the following columns:

| Column | Description |
|---|---|
| `Data` | Name of the input field or column in the source `.csv` file. |
| `Units` | Expected units for the field, written using Pyomo unit notation, e.g., `u.MW`, `u.MW * u.hr`, or `u.USD / u.MW`. |
| `Stored in Model Data Object as` | Location where the field is stored after parsing, typically under `["elements"]`, such as `["elements"]["generator"][<GEN UID>]["p_max"]`. |
| `In GTEP Model` | GTEP set, parameter, or model component initialized from the field, such as `m.generators`, `m.thermalCapacity`, or `m.storageCapacity`. |
| `File` | Input `.csv` file where the field is expected, such as `gen.csv`, `branch.csv`, `bus.csv`, `storage.csv`, `DAY_AHEAD_load.csv`, or `DAY_AHEAD_renewables.csv`. |
| `Description` | Short description of the field. |
| `Notes` | Additional information, assumptions, or special handling notes. |

The data dictionary is intended to support future unit-consistency
checks and automated tests that verify input data against expected
model units.

## Data Creation Using "data_driver" Script

This directory contains a Python script, `data_driver.py`, designed to
facilitate the creation of data files for use in IDAES-GTEP
models. The script reads input data from specified files, processes
the data, and generates new output files that are structured for easy
integration into modeling framework.

The following input files are required for the script to generate the
needed data files. Each file serves a specific purpose in providing
the necessary data for the model.

| File                | Type    | Description                                             |
|---------------------|---------|---------------------------------------------------------|
| case_name           | M       | Includes the base MVA, bus, branch, and generators data |
| DAY_AHEAD_load      | CSV     | Forecast of electricity demand expected load            |
| DAY_AHEAD_renewables| CSV     | Forecast of expected generation for renewable generators|
| REAL_TIME_load      | CSV     | Actual electricity consumption                          |
| REAL_TIME_renewables| CSV     | Actual generation of rnewable generators                |

Note: If files for day-ahead and real-time are unavailable, you may
use the example files already included in the repository for other
cases.

The script generates several output files that are structured for easy
integration into IDAES-GTEP models. Each output file contains specific
data required for modeling and simulation.

| File                  | Type    | Description                                                                                                    |
|-----------------------|---------|----------------------------------------------------------------------------------------------------------------|
| bus                   | CSV     | Includes bus ID, name, base KV, bus type, load, etc.                                                           |
| branch                | CSV     | Includes the branch ID, the from and to bus information, and reactance and rating values                       |
| gen                   | CSV     | Includes generators names, maximum and minimum operation points, ramp rates, etc.                              |
| initial_status        | CSV     | Defines the initial levels of each generator in the case                                                       |
| reserves              | CSV     | Declares the reserves products and their requirementes. By default, it is empty                                |
| simulation_parameters | CSV     | Defines all the time parameters in the model                                                                   |
| timeseries_pointers   | CSV     | Includes the information that describes the simulation category and the data files required for each simulation |


To execute the `data_driver.py` script, follow these steps:

1. Ensure that all required input files are present in this directory.

2. Run the `data_driver.py` script using Python ensuring that the
script includes the case name at the top of the file. The final data
files will be saved in a directory with the same name under this
directory.

3. Once the directory with the case name name is created, add the
day-ahead and real-time files and run again to create the last data
file.

