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
   python driver_from_config.py --config examples/config_5bus.toml
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
| `driver.py` | Explicit Python driver with model settings defined directly in the script. |

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

## Model Configuration Options

The GTEP model uses a Pyomo `ConfigBlock` to define model options that
control which parts of the formulation are included. These options
determine whether the model includes investment decisions, commitment
logic, redispatch, load scaling, transmission modeling, storage, hydro
constraints, and the selected power-flow formulation. These options
are applied before model construction. Turning an option on or off
changes which variables, constraints, disjunctions, and cost terms are
included in the model.

The available configuration options are described in the table below:

| Option | Type | Default | Description |
|---|---|---:|---|
| `include_investment` | `Bool` | `True` | Enables investment-related decisions. When disabled, candidate assets should not be selected for installation. |
| `include_commitment` | `Bool` | `True` | Enables unit commitment logic, including generator on/off operating-status decisions. |
| `include_redispatch` | `Bool` | `True` | Enables redispatch within commitment periods. This is relevant when there is more than one dispatch period per commitment period. |
| `flow_model` | `In({"DC", "CP", "ACP", "ACR", "transport"})` | `"DC"` | Selects the power-flow formulation used in the model. Available values are described in the table below. |
| `time_period_subsets` | `ConfigList` | `[]` | Optional list for defining fixed-length or fixed-subset time-period structures. |
| `time_period_dict` | `ConfigDict` | `{}` | Optional nested dictionary describing custom investment, representative, commitment, and dispatch period structures. |
| `dispatch_randomization` | `Bool` | `True` | Enables randomized dispatch information instead of fixed values per commitment period. |
| `scale_loads` | `Bool` | `True` | Enables load scaling in the model rather than directly modifying the input data. |
| `scale_texas_loads` | `Bool` | `False` | Enables Texas-case-specific load scaling logic, when applicable. |
| `thermal_generation` | `Bool` | `False` | Enables thermal generation investment options. |
| `renewable_generation` | `Bool` | `False` | Enables renewable generation investment options. |
| `storage` | `Bool` | `False` | Enables storage investment and operation modeling. |
| `transmission` | `Bool` | `False` | Enables transmission modeling and transmission investment options. |
| `transmission_switching` | `Bool` | `False` | Allows transmission switching decisions during dispatch. |
| `advanced_hydro` | `Bool` | `False` | Enables advanced hydro modeling features, including daily average hydro requirements. |

The available options for the `flow_model` configuration are listed
below:

| `flow_model` Value | Description |
|---|---|
| `"DC"` | DC power-flow approximation. |
| `"CP"` | Copper-plate power-flow approximation. |
| `"ACP"` | AC power flow in polar formulation. |
| `"ACR"` | AC power flow in rectangular formulation. |
| `"transport"` | Transport-style network flow approximation. |

