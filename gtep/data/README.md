# Data Creation Using "data_driver" Script

This directory contains a Python script, data_driver.py, designed to
facilitate the creation of data files for use in IDAES-GTEP
models. The script reads input data from specified files, processes
the data, and generates new output files that are structured for easy
integration into modeling framework.

## Script Input

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

## Script Output

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


## How to Run the Script

To execute the script, follow these steps:

1. Ensure that all required input files are present in this directory.

2. Run the "data_driver.py" script using Python ensuring that the
script includes the case name at the top of the file. The final data
files will be saved in a directory with the same name under this
directory.

3. Once the directory with the case name name is created, add the
day-ahead and real-time files and run again to create the last data
file.