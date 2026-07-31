# GTEP Examples and Configuration Files

## Overview

GTEP can be run using `main_driver.py` with a TOML configuration file.
The configuration file defines the input data path, model dimensions,
cost-data inputs, model options, transformations, solver settings, and
results options. This allows users to reproduce runs without modifying
the driver directly.

Example usage:

```bash
python main_driver.py --config examples/config_5bus.toml

## Configuration File Structure

The TOML configuration file is organized into sections. Each section
controls a different part of the GTEP workflow, such as data loading,
model setup, transformations, solver options, and result saving. The
table below shows more details about each section:

| Section | Description |
|---|---|
| `[logging]` | Sets the logging level used by the driver. Supported levels include `DEBUG`, `INFO`, `WARNING`, and `ERROR`. |
| `[data]` | Defines the input data path and time-structure settings, such as number of stages, representative periods, commitment periods, dispatch periods, and dispatch duration. |
| `[cost_data]` | Defines optional cost-data preprocessing inputs, including technology cost files, natural gas cost data, and candidate generator types. |
| `[model]` | Sets model configuration options, such as whether to include investment, commitment, redispatch, transmission, storage, load scaling, and the selected flow model. |
| `[transformations]` | Specifies which GDP transformations are applied before solving, such as `gdp.bound_pretransformation` and `gdp.bigm`. |
| `[solver]` | Defines the solver and solver output settings. Common solver options include `gurobi` and `highs`. |
| `[results]` | Defines results-saving options, including the base results directory name and the threshold used to filter near-zero values in saved JSON files. |

## Legacy Drivers

Older standalone drivers have been moved to the `legacy_drivers`
directory. These files are kept for reference, but the preferred
workflow is to use `main_driver.py` with a TOML configuration
file. The table below links the configuration file to a legacy driver,
for reference.

| Driver | Config File | Notes |
|---|---|---|
| `driver_coal.py` | `config_123bus_coal` | Error message about bus ID. |
| `driver_t2k.py` | `config_t2k` | Missing data required files. |
| `driver_config_work.py` | No | Solves for existing config file for `5bus` case. |
| `driver_jsc.py` | `config_5bus_scaled` | Solves to optimal solution |
| `driver.py` | `config_5bus_jsc` | Solves to optimal solution. |
| `driver_esr` | `config_5bus` | Solves to optimal solution. |
| `driver_matt.py` | `config_9bus` | Add `ng_cost_path` to avoid errors. Solves to optimal solution. |
| `driver_resil_week.py` | `config_123bus_resil_week` | Throws `ramp_q` error. |
| `RA_driver.py` | `config_5bus_no_commitment` | Throws a `b.loads` error. |
| `JsonPlotter.py` | No | it is only a sanity test. |