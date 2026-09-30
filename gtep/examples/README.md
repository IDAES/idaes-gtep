# GTEP Examples and Configuration Files

## Overview

GTEP can be run using `driver_from_config.py` with a TOML
configuration file.  The configuration file defines the input data
path, model dimensions, cost-data inputs, model options,
transformations, solver settings, and results options. This allows
users to reproduce runs without modifying the driver directly.

Example usage:

```bash
python driver_from_config.py --config examples/config_5bus.toml
```

## Configuration Template

The file `config_template.toml` provides a starting point for creating
new GTEP configuration files. Users can copy this template, rename it
for a specific case, and update the paths and options as needed. The
template includes all major configuration sections used by
`driver_from_config.py`, including data settings, cost-data inputs,
model options, GDP transformations, solver settings, results options,
and optional plotting settings.

Example:

```bash
cp examples/config_template.toml examples/config_my_case.toml
python driver_from_config.py --config examples/config_my_case.toml
```

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