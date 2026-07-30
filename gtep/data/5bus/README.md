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

## References

1. (Pending)

