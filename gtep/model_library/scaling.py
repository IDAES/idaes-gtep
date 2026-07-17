#################################################################################
# The Institute for the Design of Advanced Energy Systems Integrated Platform
# Framework (IDAES IP) was produced under the DOE Institute for the
# Design of Advanced Energy Systems (IDAES).
#
# Copyright (c) 2018-2026 by the software owners: The Regents of the
# University of California, through Lawrence Berkeley National Laboratory,
# National Technology & Engineering Solutions of Sandia, LLC, Carnegie Mellon
# University, West Virginia University Research Corporation, et al.
# All rights reserved.  Please see the files COPYRIGHT.md and LICENSE.md
# for full copyright and license information.
#################################################################################

"""Scaling for Generation and Transmission Expansion Planning (GTEP)
Model

"""

import pyomo.environ as pyo
from pyomo.environ import units as u


def add_load_scaling(m, b, commitment_period, investment_stage, scaling_value):
    """Add commitment-period load parameters.

    Loads are operational time-series quantities and should be read from the
    scenario-specific representative-period ModelData object, ``m.md``. For PH
    scenario models, ``m.md`` is the selected representative-period clone.

    The lookup uses the scenario-local time key for the commitment period. If
    the selected representative-period clone does not contain a numeric value at
    that local index, the code falls back to the full source ModelData stored on
    ``m.data.full_md`` and locates the same time key there.

    Missing load data are not silently treated as zero; they raise an error with
    enough context to diagnose the time-series alignment.
    """

    b.loads = pyo.Param(
        m.buses,
        initialize={bus: 0.0 for bus in m.buses},
        mutable=True,
        units=u.MW,
        within=pyo.NonNegativeReals,
        doc="Demand at each bus",
    )

    for load_n in m.load_buses:
        p_load = _get_load_value_for_commitment_period(
            m,
            load_n,
            commitment_period,
        )

        if m.config["scale_loads"]:
            temp_scale = scaling_value
            scaled_load = (
                p_load
                * temp_scale
                * (1 + (temp_scale + investment_stage) / (temp_scale + len(m.stages)))
            )

        elif m.config["scale_texas_loads"]:
            scaled_load = (
                p_load
                * b.load_scaling[m.md.data["elements"]["load"][load_n]["zone"]].iloc[0]
            )

        else:
            scaled_load = p_load

        load_bus = _get_load_bus(m, load_n)

        current_load = pyo.value(b.loads[load_bus], exception=False)
        if current_load is None:
            current_load = 0.0

        b.loads[load_bus] = float(current_load) + float(scaled_load)


def _get_load_value_for_commitment_period(m, load_n, commitment_period):
    """Return load value aligned with the scenario-local commitment period.

    Parameters
    ----------
    m:
        GTEP Pyomo model.
    load_n:
        Load element identifier.
    commitment_period:
        One-based commitment-period index.

    Returns
    -------
    float
        Numeric load value in MW.

    Raises
    ------
    ValueError
        If no numeric value can be found in the scenario-local ModelData or the
        full source ModelData for the relevant time key.
    """
    local_index = commitment_period - 1
    scenario_md = m.md

    scenario_time_key = _get_time_key_at_index(
        scenario_md,
        local_index,
        context=f"scenario load {load_n}",
    )

    scenario_value = _get_time_series_value_at_index(
        scenario_md,
        element_type="load",
        element_name=load_n,
        field="p_load",
        index=local_index,
    )

    if scenario_value is not None:
        return float(scenario_value)

    full_md = getattr(m.data, "full_md", None)

    if full_md is not None:
        full_time_keys = full_md.data["system"]["time_keys"]

        if scenario_time_key in full_time_keys:
            full_index = full_time_keys.index(scenario_time_key)

            full_value = _get_time_series_value_at_index(
                full_md,
                element_type="load",
                element_name=load_n,
                field="p_load",
                index=full_index,
            )

            if full_value is not None:
                return float(full_value)

    raise ValueError(
        "Could not determine numeric p_load value for load "
        f"{load_n} at commitment_period={commitment_period}, "
        f"scenario_time_key={scenario_time_key!r}. The scenario-local "
        "representative-period data and full source data did not contain a "
        "numeric value at the aligned time index."
    )


def _get_load_bus(m, load_n):
    """Return the bus associated with a load element.

    The current 5-bus data uses load names that match bus names, but general
    data may store the bus explicitly on the load element.
    """
    load_data = m.md.data["elements"]["load"][load_n]
    load_bus = load_data.get("bus", load_n)

    if load_bus not in m.buses:
        raise KeyError(
            f"Load {load_n} maps to bus {load_bus}, but that bus is not in m.buses."
        )

    return load_bus


def _get_time_key_at_index(model_data, index, context):
    """Return time key at a local integer index."""
    time_keys = model_data.data["system"]["time_keys"]

    if index < 0 or index >= len(time_keys):
        raise IndexError(
            f"Time index {index} is out of range for {context}. "
            f"Number of time keys is {len(time_keys)}."
        )

    return time_keys[index]


def _get_time_series_value_at_index(
    model_data,
    *,
    element_type,
    element_name,
    field,
    index,
):
    """Return a numeric time-series value at an integer index, or None.

    This helper supports the common Egret/Prescient structure where a time
    series is represented as a dictionary with a ``values`` list. It also
    supports scalar numeric fields.
    """
    element = model_data.data["elements"][element_type][element_name]
    series = element.get(field, None)

    if series is None:
        return None

    if isinstance(series, (int, float)):
        return float(series)

    if isinstance(series, dict):
        values = series.get("values", None)

        if values is not None:
            if index < 0 or index >= len(values):
                raise IndexError(
                    f"Time-series index {index} is out of range for "
                    f"{element_type} {element_name} field {field}. "
                    f"Number of values is {len(values)}."
                )

            value = values[index]
            if value is None:
                return None

            return float(value)

    raise TypeError(
        f"Unsupported time-series data structure for {element_type} "
        f"{element_name} field {field}: {type(series)}"
    )
