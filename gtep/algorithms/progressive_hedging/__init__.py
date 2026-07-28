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
"""Progressive Hedging for GTEP representative-period scenarios.

This package implements a Progressive Hedging workflow for the GTEP model where
each representative period is treated as one PH scenario.

Public design guarantees
------------------------
* Representative periods are PH scenarios.
* Scenario solves build and GDP-transform one representative-period model.
* No Pyomo model serialization is required.
* No persistent solver interface is required.
* No warm starts are required.
* User-facing configuration is YAML.
* Internal PH state and results are JSON.
* Nonanticipative variables are identified using Pyomo ``ComponentUID`` through
  the ``VariableID`` abstraction, not by variable names or manual string
  matching.
* Initial PH regularization is MILP-compatible:
    - binary variables use a linearized binary quadratic PH term,
    - non-binary variables use L1 absolute-deviation regularization.
"""

from gtep.algorithms.progressive_hedging.config import (
    CostDataConfig,
    DataConfig,
    ModelConfig,
    NonanticipativityConfig,
    OutputConfig,
    PHConfig,
    ProgressiveHedgingConfig,
    RunConfig,
    SolverConfig,
    TorcConfig,
    TorcResourceConfig,
    load_ph_config,
    prepare_output_directory,
)
from gtep.algorithms.progressive_hedging.convergence import (
    PHUpdateResult,
    ResidualSummary,
    compute_consensus,
    compute_residuals,
    process_completed_iteration,
    update_multipliers,
)
from gtep.algorithms.progressive_hedging.nonanticipativity import (
    NonanticipativeVariableMetadata,
    VariableID,
    bind_variable_ids_to_model,
    collect_nonanticipative_variables,
    deserialize_variable_value_records,
    extract_nonanticipative_values,
    serialize_variable_value_records,
)
from gtep.algorithms.progressive_hedging.ph_model import (
    PHObjectiveOptions,
    add_progressive_hedging_block,
    deactivate_active_objectives,
    update_progressive_hedging_parameters,
)
from gtep.algorithms.progressive_hedging.scenario_data import (
    ScenarioInfo,
    build_cost_data,
    build_full_data,
    build_scenario_infos,
    get_representative_period_weights,
    get_scenario_info,
    get_scenario_probabilities,
    make_single_representative_period_data,
)
from gtep.algorithms.progressive_hedging.solution_io import (
    IterationSummary,
    ScenarioResult,
    read_scenario_result,
    read_scenario_results,
    write_final_solution,
    write_iteration_summary,
    write_scenario_result,
)
from gtep.algorithms.progressive_hedging.solver import (
    SolveOutcome,
    create_solver,
    solve_model,
)
from gtep.algorithms.progressive_hedging.state import (
    PHState,
    create_enabled_state,
    create_initial_state,
    read_state,
    write_state,
)

__all__ = [
    # config
    "CostDataConfig",
    "DataConfig",
    "ModelConfig",
    "NonanticipativityConfig",
    "OutputConfig",
    "PHConfig",
    "ProgressiveHedgingConfig",
    "RunConfig",
    "SolverConfig",
    "TorcConfig",
    "TorcResourceConfig",
    "load_ph_config",
    "prepare_output_directory",
    # convergence
    "PHUpdateResult",
    "ResidualSummary",
    "compute_consensus",
    "compute_residuals",
    "process_completed_iteration",
    "update_multipliers",
    # nonanticipativity
    "NonanticipativeVariableMetadata",
    "VariableID",
    "bind_variable_ids_to_model",
    "collect_nonanticipative_variables",
    "deserialize_variable_value_records",
    "extract_nonanticipative_values",
    "serialize_variable_value_records",
    # ph_model
    "PHObjectiveOptions",
    "add_progressive_hedging_block",
    "deactivate_active_objectives",
    "update_progressive_hedging_parameters",
    # scenario_data
    "ScenarioInfo",
    "build_cost_data",
    "build_full_data",
    "build_scenario_infos",
    "get_representative_period_weights",
    "get_scenario_info",
    "get_scenario_probabilities",
    "make_single_representative_period_data",
    # solution_io
    "IterationSummary",
    "ScenarioResult",
    "read_scenario_result",
    "read_scenario_results",
    "write_final_solution",
    "write_iteration_summary",
    "write_scenario_result",
    # solver
    "SolveOutcome",
    "create_solver",
    "solve_model",
    # state
    "PHState",
    "create_enabled_state",
    "create_initial_state",
    "read_state",
    "write_state",
]
