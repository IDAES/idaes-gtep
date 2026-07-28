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
"""Scenario-data utilities for GTEP Progressive Hedging.

In the representative-period PH formulation, each representative period is one
PH scenario. This module provides helpers to:

* build the full ``ExpansionPlanningData`` object from user configuration,
* build the optional cost-data object,
* create a single-representative-period ``ExpansionPlanningData`` object for a
  scenario solve,
* compute representative-period weights and normalized scenario probabilities.

This module does not perform any nonanticipative-variable matching.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any

import copy
import logging

from gtep.gtep_data import ExpansionPlanningData
from gtep.gtep_data_processing import DataProcessing

from gtep.algorithms.progressive_hedging.config import (
    CostDataConfig,
    DataConfig,
    PHConfig,
)

logger = logging.getLogger("gtep.algorithms.progressive_hedging.scenario_data")


@dataclass(frozen=True)
class ScenarioInfo:
    """Metadata for one representative-period PH scenario.

    Attributes
    ----------
    scenario_id:
        One-based PH scenario identifier. This intentionally matches the GTEP
        representative-period convention used by ``m.representativePeriods``.
    representative_period:
        One-based representative-period identifier from the full data object.
    representative_date:
        Representative-period timestamp.
    representative_weight:
        Original representative-period weight \(W_s\).
    probability:
        Scenario probability \(\pi_s\), usually normalized from the
        representative weights.
    """

    scenario_id: int
    representative_period: int
    representative_date: str
    representative_weight: float
    probability: float


def build_full_data(config: PHConfig | DataConfig) -> ExpansionPlanningData:
    """Build and load the full ``ExpansionPlanningData`` object.

    Parameters
    ----------
    config:
        Either a top-level ``PHConfig`` or a ``DataConfig``.

    Returns
    -------
    ExpansionPlanningData
        Fully loaded data object containing all representative periods.
    """
    data_config = config.data if isinstance(config, PHConfig) else config

    data_object = ExpansionPlanningData(
        stages=data_config.stages,
        num_reps=data_config.num_reps,
        len_reps=data_config.len_reps,
        num_commit=data_config.num_commit,
        num_dispatch=data_config.num_dispatch,
        duration_dispatch=data_config.duration_dispatch,
    )

    logger.info(
        "Loading Prescient data for %s representative period(s) from %s",
        data_config.num_reps,
        data_config.data_path,
    )

    data_object.load_prescient(
        data_config.data_path,
        representative_dates=data_config.representative_dates,
        representative_weights=(
            data_config.representative_weights
            if data_config.representative_weights is not None
            else {}
        ),
        options_dict=data_config.prescient_options,
    )

    _validate_loaded_representative_data(data_object)
    return data_object


def build_cost_data(config: PHConfig | CostDataConfig) -> DataProcessing | None:
    """Build and load the optional GTEP cost-data object.

    Parameters
    ----------
    config:
        Either a top-level ``PHConfig`` or a ``CostDataConfig``.

    Returns
    -------
    DataProcessing | None
        Loaded ``DataProcessing`` object if enabled, otherwise ``None``.
    """
    cost_config = config.cost_data if isinstance(config, PHConfig) else config

    if not cost_config.enabled:
        logger.info("Cost-data preprocessing is disabled.")
        return None

    required_paths = {
        "bus_data_path": cost_config.bus_data_path,
        "cost_data_path": cost_config.cost_data_path,
        "ng_cost_path": cost_config.ng_cost_path,
    }

    missing = [name for name, value in required_paths.items() if value is None]
    if missing:
        raise ValueError(
            "Cost-data preprocessing is enabled, but the following required "
            f"cost_data path(s) are missing: {', '.join(missing)}"
        )

    if not cost_config.candidate_gens:
        raise ValueError(
            "Cost-data preprocessing is enabled, but cost_data.candidate_gens "
            "is empty."
        )

    logger.info("Loading GTEP cost data.")

    cost_data = DataProcessing()
    cost_data.load_gen_data(
        cost_config.bus_data_path,
        cost_config.cost_data_path,
        cost_config.ng_cost_path,
        cost_config.candidate_gens,
    )

    return cost_data


def make_single_representative_period_data(
    full_data: ExpansionPlanningData,
    scenario_id: int,
    *,
    copy_representative_model_data: bool = False,
) -> ExpansionPlanningData:
    """Create a one-representative-period data object for a PH scenario.

    Parameters
    ----------
    full_data:
        Loaded ``ExpansionPlanningData`` object containing all representative
        periods.
    scenario_id:
        One-based representative-period scenario identifier.
    copy_representative_model_data:
        If true, deep-copy the selected Egret ``ModelData`` object. The default
        is false because scenario solves are expected to run in separate Torc
        processes and do not need cross-scenario in-process isolation.

    Returns
    -------
    ExpansionPlanningData
        Data object with ``num_reps == 1`` and a single representative-period
        ``ModelData`` entry.

    Notes
    -----
    The returned object intentionally preserves the original representative
    period's weight. In a one-representative-period scenario model,
    ``m.weights[1]`` should be the original \(W_s\), not the normalized
    probability \(\pi_s\). PH objective construction separately applies the
    normalized scenario probability to investment costs.
    """
    _validate_loaded_representative_data(full_data)
    _validate_scenario_id(full_data, scenario_id)

    rep_index = scenario_id - 1
    representative_date = full_data.representative_dates[rep_index]
    representative_weight = full_data.representative_weights_dict[representative_date]

    representative_model_data = full_data.representative_data[rep_index]
    if copy_representative_model_data:
        representative_model_data = copy.deepcopy(representative_model_data)

    scenario_data = ExpansionPlanningData(
        stages=full_data.stages,
        num_reps=1,
        len_reps=full_data.len_reps,
        num_commit=full_data.num_commit,
        num_dispatch=full_data.num_dispatch,
        duration_dispatch=full_data.duration_dispatch,
    )

    # Attributes expected by ExpansionPlanningModel.create_model().
    scenario_data.data_type = getattr(full_data, "data_type", None)
    scenario_data.md = representative_model_data
    scenario_data.representative_dates = [representative_date]
    scenario_data.representative_weights_dict = {
        representative_date: representative_weight
    }
    scenario_data.representative_data = [representative_model_data]

    # Preserve full source data for scenario-invariant parameter construction.
    # The scenario model uses only one representative-period ModelData object
    # for operational time-series data, but physical nameplate quantities should
    # not be inferred only from that selected representative day.
    scenario_data.full_md = full_data.md
    scenario_data.full_representative_data = full_data.representative_data

    # Preserve optional auxiliary data used by downstream model components.
    _copy_optional_data_attributes(full_data, scenario_data)

    return scenario_data


def get_representative_period_weights(
    full_data: ExpansionPlanningData,
) -> dict[int, float]:
    """Return original representative-period weights indexed by scenario id.

    Parameters
    ----------
    full_data:
        Loaded full data object.

    Returns
    -------
    dict[int, float]
        Mapping from one-based scenario id to original representative-period
        weight \(W_s\).
    """
    _validate_loaded_representative_data(full_data)

    weights: dict[int, float] = {}
    for scenario_id, representative_date in enumerate(
        full_data.representative_dates,
        start=1,
    ):
        weights[scenario_id] = float(
            full_data.representative_weights_dict[representative_date]
        )

    return weights


def get_scenario_probabilities(
    full_data: ExpansionPlanningData,
    *,
    normalize: bool = True,
) -> dict[int, float]:
    """Return PH scenario probabilities indexed by scenario id.

    Parameters
    ----------
    full_data:
        Loaded full data object.
    normalize:
        If true, normalize representative weights so that probabilities sum to
        one. If false, return the weights directly. The PH implementation should
        normally use normalized probabilities.

    Returns
    -------
    dict[int, float]
        Mapping from one-based scenario id to \(\pi_s\).
    """
    weights = get_representative_period_weights(full_data)

    if not normalize:
        return dict(weights)

    total_weight = sum(weights.values())
    if total_weight <= 0:
        raise ValueError(
            "Cannot normalize representative-period weights because their sum "
            f"is nonpositive: {total_weight}"
        )

    return {
        scenario_id: weight / total_weight for scenario_id, weight in weights.items()
    }


def build_scenario_infos(
    full_data: ExpansionPlanningData,
    *,
    normalize_probabilities: bool = True,
) -> list[ScenarioInfo]:
    """Build metadata records for all representative-period PH scenarios.

    Parameters
    ----------
    full_data:
        Loaded full data object.
    normalize_probabilities:
        Whether to normalize representative-period weights into probabilities.

    Returns
    -------
    list[ScenarioInfo]
        One record per representative period.
    """
    _validate_loaded_representative_data(full_data)

    weights = get_representative_period_weights(full_data)
    probabilities = get_scenario_probabilities(
        full_data,
        normalize=normalize_probabilities,
    )

    scenario_infos: list[ScenarioInfo] = []
    for scenario_id, representative_date in enumerate(
        full_data.representative_dates,
        start=1,
    ):
        scenario_infos.append(
            ScenarioInfo(
                scenario_id=scenario_id,
                representative_period=scenario_id,
                representative_date=representative_date,
                representative_weight=weights[scenario_id],
                probability=probabilities[scenario_id],
            )
        )

    return scenario_infos


def get_scenario_info(
    full_data: ExpansionPlanningData,
    scenario_id: int,
    *,
    normalize_probabilities: bool = True,
) -> ScenarioInfo:
    """Return metadata for one representative-period PH scenario."""
    _validate_loaded_representative_data(full_data)
    _validate_scenario_id(full_data, scenario_id)

    representative_date = full_data.representative_dates[scenario_id - 1]
    weight = float(full_data.representative_weights_dict[representative_date])
    probabilities = get_scenario_probabilities(
        full_data,
        normalize=normalize_probabilities,
    )

    return ScenarioInfo(
        scenario_id=scenario_id,
        representative_period=scenario_id,
        representative_date=representative_date,
        representative_weight=weight,
        probability=probabilities[scenario_id],
    )


def scenario_infos_to_jsonable(
    scenario_infos: list[ScenarioInfo],
) -> list[dict[str, Any]]:
    """Convert scenario metadata records to JSON-serializable dictionaries."""
    return [
        {
            "scenario_id": info.scenario_id,
            "representative_period": info.representative_period,
            "representative_date": info.representative_date,
            "representative_weight": info.representative_weight,
            "probability": info.probability,
        }
        for info in scenario_infos
    ]


def _copy_optional_data_attributes(
    source: ExpansionPlanningData,
    target: ExpansionPlanningData,
) -> None:
    """Copy optional auxiliary attributes if they exist on the full data object."""
    optional_attributes = [
        "load_scaling",
        "bus_hours",
    ]

    for attr in optional_attributes:
        if hasattr(source, attr):
            setattr(target, attr, getattr(source, attr))


def _validate_loaded_representative_data(data: ExpansionPlanningData) -> None:
    """Validate that an ``ExpansionPlanningData`` object has representative data."""
    required_attributes = [
        "representative_dates",
        "representative_weights_dict",
        "representative_data",
    ]

    missing = [attr for attr in required_attributes if not hasattr(data, attr)]
    if missing:
        raise ValueError(
            "ExpansionPlanningData object has not been loaded with representative "
            f"period data. Missing attribute(s): {', '.join(missing)}"
        )

    if data.num_reps <= 0:
        raise ValueError(f"data.num_reps must be positive. Received {data.num_reps}.")

    if len(data.representative_dates) != data.num_reps:
        raise ValueError(
            "Length of data.representative_dates must equal data.num_reps. "
            f"Received {len(data.representative_dates)} dates and "
            f"num_reps={data.num_reps}."
        )

    if len(data.representative_data) != data.num_reps:
        raise ValueError(
            "Length of data.representative_data must equal data.num_reps. "
            f"Received {len(data.representative_data)} data objects and "
            f"num_reps={data.num_reps}."
        )

    for representative_date in data.representative_dates:
        if representative_date not in data.representative_weights_dict:
            raise ValueError(
                "Missing representative weight for representative date "
                f"{representative_date!r}."
            )

    weights = [
        float(data.representative_weights_dict[date])
        for date in data.representative_dates
    ]

    if any(weight < 0 for weight in weights):
        raise ValueError(
            "Representative-period weights must be nonnegative. Received " f"{weights}."
        )

    if sum(weights) <= 0:
        raise ValueError(
            "At least one representative-period weight must be positive. "
            f"Received {weights}."
        )


def _validate_scenario_id(
    full_data: ExpansionPlanningData,
    scenario_id: int,
) -> None:
    """Validate a one-based representative-period scenario id."""
    if not isinstance(scenario_id, int):
        raise TypeError(
            f"scenario_id must be an integer. Received {type(scenario_id)}."
        )

    if scenario_id < 1 or scenario_id > full_data.num_reps:
        raise ValueError(
            f"scenario_id must be between 1 and {full_data.num_reps}, inclusive. "
            f"Received {scenario_id}."
        )
