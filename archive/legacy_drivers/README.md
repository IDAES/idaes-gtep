# GTEP Examples and Configuration Files

## Legacy Drivers

This directory includes older standalone drivers. These files are kept
for reference, but the preferred workflow is to use either
`gtep/driver.py` or `gtep/driver_from_config.py` with a TOML
configuration file. The table below links the configuration file to a
legacy driver, for reference.

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
| `JsonPlotter.py` | No | It is only a sanity test. |