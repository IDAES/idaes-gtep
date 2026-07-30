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