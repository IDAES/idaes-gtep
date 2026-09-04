# 5-Bus Case Data

## Introduction

This directory contains the input data files for the 5-bus test case
used in the expansion planning model GTEP. The case is based on the
original 5-bus dataset, with selected modifications made to support
model testing and debugging.

## Modifications from Original Dataset

This dataset applies the following modifications:

| File | Modification | Reason |
|---|---|---|
| `gen.csv` | Removed heat rate data from the `HR_avg_0`, `HR_incr_1`, `HR_incr_2`, and `HR_incr_3` columns | The original data did not provide units for these values |

The following heat rate values were removed from `gen.csv`:

| GEN UID | HR_avg_0 | HR_incr_1 | HR_incr_2 | HR_incr_3 |
|---|---:|---:|---:|---:|
| `3_CT` | 135722.5 | 97862.5 | 98072.5 | 107135 |
| `10_STEAM` | 30166.66667 | 14402.61434 | 17182.46753 | 18283.66017 |
| `4_CC` | 51019.54545 | 26818.18182 | 29550.90909 | 30308.18182 |
| `4_STEAM` | 179457.9999 | 124504.3483 | 125050 | 133643.478 |

These values were removed because their units were not specified in
the original dataset. Since heat rate values are used in fuel cost
calculations, retaining values with unclear units could lead to
incorrect cost estimates.

## Storage data

Storage data has been added (`storage.csv`). This data is taken from `9_bus_GTEP_dir/storage.csv`,
with a few modifications:
- The `bus` column was adjusted to reflect valid IDs from this directory's `bus.csv`
- The `name` column was similarly changed
- The following columns were scaled down by 100x to better match the scale of
generation/loads for the 5bus case: `energy_capacity`, `initial_state_of_charge`,
`end_state_of_charge`, `minimum_state_of_charge`, `max_discharge_rate`, `min_discharge_rate`,
`max_charge_rate`, `min_charge_rate`, `ramp_up_input_60min`,`ramp_down_input_60min`,
`ramp_up_output_60min`, `ramp_down_output_60min`, and `investment_cost`

## Branch data modifications

A few further changes were made to `branch.csv`:
- Quotes removed from branch UIDs
- Columns `Cont Rating`, `LTE Rating`, and `STE Rating` were scaled up by 10x to better match
the scale of generation/loads for the 5bus case.

## References

1. (Pending)
