# GTEP Code Overview

## Overview

This directory contains the main modules used to build, solve,
post-process, and run the Generation and Transmission Expansion
Planning (GTEP) model. The model is organized around separate scripts
for data loading, cost data preprocessing, model construction,
solution processing, and driver scripts.

At a high level, the GTEP model involves four steps:

1. Load and preprocess data into the model as a data object.
2. Setup GTEP model configuration.
3. Create, transformation, and solve the model.
4. Create result visualizations.

## How to Run GTEP

There are two main ways to run the GTEP model:

1. Using an explicit Python driver such as `driver.py`.

   ```bash
   python driver.py
   ```

2. Using a TOML configuration file with `driver_from_config.py`.

   ```bash
   python main_driver.py --config examples/config_5bus.toml
   ```

## Main Files

| File | Description |
|---|---|
| `gtep_model.py` | Defines the main `ExpansionPlanningModel` class, which builds the Pyomo/GDP optimization model and applies model configuration options. |
| `gtep_data.py` | Defines the `ExpansionPlanningData` class, which loads Prescient-formatted input data and prepares the model data object used by GTEP. |
| `gtep_data_processing.py` | Defines the `DataProcessing` class, which preprocesses external cost data. |
| `gtep_solution.py` | Defines the `ExpansionPlanningSolution` class, which saves model results to JSON/CSV files and generates plots. |
| `config_options.py` | Defines the available GTEP model configuration options. |
| `driver_from_config.py` | Configuration-driven driver for running GTEP using a TOML configuration file. |
| `driver_full_options.py` | Explicit Python driver with model settings defined directly in the script. |

## Supporting Model Modules

In addition to the main files above, GTEP uses supporting modules that
define specific parts of the model formulation under the directory
`model_library`. below find more details about these modules.

| Module | Description |
|---|---|
| `gen.py` | Defines generator-specific variables, parameters, constraints, and operating logic. |
| `transmission.py` | Defines transmission branch-specific variables, parameters, constraints, power limits, and operating logic. |
| `storage.py` | Defines storage-specific variables, parameters, investment logic, charge/discharge behavior, state-of-charge tracking, and operating constraints. |
| `commitment.py` | Defines commitment-period variables, parameters, constraints, and operating-cost aggregation. |
| `dispatch.py` | Defines dispatch-period variables and constraints, including power balance, generation, load shedding, curtailment, reserves, storage operation, and transmission flow. |
| `investment.py` | Defines investment-stage structure, stage-level variables, and investment-related cost aggregation. |
| `objective.py` | Defines the model objective function and total cost expressions. |
| `scaling.py` | Provides load-scaling utilities and related model adjustments. |
| `hydro.py` | Defines hydropower-specific variables and constraints used when advanced hydro modeling is enabled. |

