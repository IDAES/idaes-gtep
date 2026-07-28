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
"""Configuration utilities for GTEP Progressive Hedging.

This module defines the user-facing YAML configuration schema for the
representative-period Progressive Hedging implementation.

Important design constraint:
    This module does not define or process nonanticipative variables by
    component names or string paths. Variable identity and matching are handled
    elsewhere using Pyomo ComponentUID objects.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass, field, is_dataclass
from pathlib import Path
from typing import Any

import copy
import json
import logging
import shutil

logger = logging.getLogger("gtep.algorithms.progressive_hedging.config")


try:
    import yaml
except ImportError as err:  # pragma: no cover
    yaml = None
    _YAML_IMPORT_ERROR = err
else:
    _YAML_IMPORT_ERROR = None


SUPPORTED_REGULARIZATIONS = {"l1"}
SUPPORTED_GDP_TRANSFORMATIONS = {"gdp.bigm"}
SUPPORTED_HISTORY_FORMATS = {"json", "csv"}


@dataclass
class RunConfig:
    """Run-level configuration."""

    name: str = "gtep_ph_run"
    output_dir: Path = Path("./ph_runs/gtep_ph_run")
    clean_output_dir: bool = False


@dataclass
class DataConfig:
    """Input data configuration for ExpansionPlanningData."""

    data_path: Path = Path("./gtep/data/123_Bus_Resil_Week")

    stages: int = 2
    num_reps: int = 4
    len_reps: int = 24
    num_commit: int = 24
    num_dispatch: int = 1
    duration_dispatch: int = 60

    representative_dates: list[str] | None = None
    representative_weights: list[float] | None = None

    prescient_options: dict[str, Any] | None = None


@dataclass
class CostDataConfig:
    """Cost-data preprocessing configuration."""

    bus_data_path: Path | None = None
    cost_data_path: Path | None = None
    ng_cost_path: Path | None = None
    candidate_gens: list[str] = field(default_factory=list)

    enabled: bool = True


@dataclass
class ModelConfig:
    """GTEP model construction options.

    These options are copied into ``ExpansionPlanningModel.config`` before
    ``create_model()`` is called.
    """

    include_investment: bool = True
    include_commitment: bool = True
    include_redispatch: bool = True
    scale_loads: bool = True
    transmission: bool = True
    storage: bool = False
    flow_model: str = "transport"
    advanced_hydro: bool = False

    gdp_transformation: str = "gdp.bigm"

    def as_model_options(self) -> dict[str, Any]:
        """Return only options intended for ``ExpansionPlanningModel.config``."""
        return {
            "include_investment": self.include_investment,
            "include_commitment": self.include_commitment,
            "include_redispatch": self.include_redispatch,
            "scale_loads": self.scale_loads,
            "transmission": self.transmission,
            "storage": self.storage,
            "flow_model": self.flow_model,
            "advanced_hydro": self.advanced_hydro,
        }


@dataclass
class NonanticipativityConfig:
    """Configuration for selecting nonanticipative variable groups.

    No variable-name matching is performed from this configuration. These flags
    describe semantic groups. The actual collector must use Pyomo components and
    ComponentUIDs.
    """

    include_renewable_investment_variables: bool = True
    include_thermal_status_binaries: bool = True
    include_transmission_status_binaries: bool = True
    include_storage_status_binaries: bool = True


@dataclass
class ProgressiveHedgingConfig:
    """Progressive Hedging algorithm configuration."""

    max_iterations: int = 50

    # Require at least this many completed PH iterations before accepting
    # convergence. This is useful for testing PH-augmented iteration logic on
    # small cases that converge at iteration 0.
    min_iterations: int = 0

    convergence_tolerance: float = 1.0e-4

    rho: float = 100.0
    regularization: str = "l1"

    normalize_scenario_probabilities: bool = True

    nonanticipativity: NonanticipativityConfig = field(
        default_factory=NonanticipativityConfig
    )


@dataclass
class SolverConfig:
    """Pyomo SolverFactory configuration."""

    name: str = "xpress"
    tee: bool = True

    mip_gap: float | None = 1.0e-3
    time_limit: float | None = 3600.0
    threads: int | None = 8

    log_file: Path | None = None

    # Optional solver license path. For XPRESS, this is commonly supplied
    # through the XPAUTH_PATH environment variable. The environment variable is
    # configurable because solver/license deployments vary.
    license_file: Path | None = None
    license_env_var: str | None = None

    # Additional environment variables to set before constructing/solving with
    # the Pyomo solver. Values are strings and are not interpreted by PH.
    environment: dict[str, str] = field(default_factory=dict)

    acceptable_termination_conditions: list[str] = field(
        default_factory=lambda: [
            "optimal",
            "feasible",
            "maxTimeLimit",
        ]
    )


@dataclass
class TorcResourceConfig:
    """Torc resource-requirement names."""

    orchestrator: str = "ph_orch_rr"
    scenario: str = "ph_scenario_rr"


@dataclass
class TorcConfig:
    """Torc dynamic-orchestrator configuration."""

    enabled: bool = True
    resource_requirements: TorcResourceConfig = field(
        default_factory=TorcResourceConfig
    )

    scenario_priority: int = 10
    orchestrator_priority: int = 0
    cancel_on_scenario_failure: bool = False


@dataclass
class OutputConfig:
    """Output and solution-persistence configuration."""

    save_iteration_solutions: bool = True
    save_final_solution: bool = True

    # Optional debug model export. This can create very large LP files and
    # should remain false for production-size runs.
    save_debug_model: bool = False

    save_nonanticipative_metadata: bool = False

    # Richer scenario output is explicitly optional and not required initially.
    save_full_scenario_solution: bool = False
    save_operational_variables: bool = False

    binary_rounding_tolerance: float = 1.0e-5

    convergence_history_formats: list[str] = field(
        default_factory=lambda: ["json", "csv"]
    )


@dataclass
class PHConfig:
    """Top-level Progressive Hedging configuration."""

    run: RunConfig = field(default_factory=RunConfig)
    data: DataConfig = field(default_factory=DataConfig)
    cost_data: CostDataConfig = field(default_factory=CostDataConfig)
    model: ModelConfig = field(default_factory=ModelConfig)
    progressive_hedging: ProgressiveHedgingConfig = field(
        default_factory=ProgressiveHedgingConfig
    )
    solver: SolverConfig = field(default_factory=SolverConfig)
    torc: TorcConfig = field(default_factory=TorcConfig)
    output: OutputConfig = field(default_factory=OutputConfig)

    config_path: Path | None = None


def load_ph_config(config_path: str | Path) -> PHConfig:
    """Load and validate a PH configuration from YAML.

    Parameters
    ----------
    config_path:
        Path to the user-facing YAML configuration file.

    Returns
    -------
    PHConfig
        Validated configuration object.
    """
    if yaml is None:  # pragma: no cover
        raise ImportError(
            "PyYAML is required to load Progressive Hedging YAML configs."
        ) from _YAML_IMPORT_ERROR

    path = Path(config_path).resolve()

    with path.open("r", encoding="utf-8") as f:
        raw = yaml.safe_load(f)

    if raw is None:
        raw = {}

    if not isinstance(raw, dict):
        raise TypeError(
            f"PH configuration must be a YAML mapping at top level. Received: {type(raw)}"
        )

    cfg = config_from_dict(raw)
    cfg.config_path = path
    validate_ph_config(cfg)
    return cfg


def config_from_dict(raw: dict[str, Any]) -> PHConfig:
    """Construct a ``PHConfig`` from a nested dictionary."""
    known_top_level = {
        "run",
        "data",
        "cost_data",
        "model",
        "progressive_hedging",
        "solver",
        "torc",
        "output",
    }
    unknown = set(raw) - known_top_level
    if unknown:
        raise ValueError(
            "Unknown top-level PH configuration section(s): "
            + ", ".join(sorted(unknown))
        )

    run_raw = raw.get("run", {})
    data_raw = raw.get("data", {})
    cost_raw = raw.get("cost_data", {})
    model_raw = raw.get("model", {})
    ph_raw = raw.get("progressive_hedging", {})
    solver_raw = raw.get("solver", {})
    torc_raw = raw.get("torc", {})
    output_raw = raw.get("output", {})

    nonant_raw = ph_raw.pop("nonanticipativity", {})
    rr_raw = torc_raw.pop("resource_requirements", {})

    cfg = PHConfig(
        run=_make_dataclass(RunConfig, run_raw),
        data=_make_dataclass(DataConfig, data_raw),
        cost_data=_make_dataclass(CostDataConfig, cost_raw),
        model=_make_dataclass(ModelConfig, model_raw),
        progressive_hedging=ProgressiveHedgingConfig(
            **_filter_dataclass_kwargs(ProgressiveHedgingConfig, ph_raw),
            nonanticipativity=_make_dataclass(
                NonanticipativityConfig,
                nonant_raw,
            ),
        ),
        solver=_make_dataclass(SolverConfig, solver_raw),
        torc=TorcConfig(
            **_filter_dataclass_kwargs(TorcConfig, torc_raw),
            resource_requirements=_make_dataclass(
                TorcResourceConfig,
                rr_raw,
            ),
        ),
        output=_make_dataclass(OutputConfig, output_raw),
    )

    cfg = _normalize_paths(cfg)
    return cfg


def validate_ph_config(cfg: PHConfig) -> None:
    """Validate a loaded PH configuration.

    Raises
    ------
    ValueError
        If any configuration value is invalid.
    """
    if not cfg.run.name:
        raise ValueError("run.name must be non-empty.")

    if cfg.data.stages <= 0:
        raise ValueError("data.stages must be positive.")
    if cfg.data.num_reps <= 0:
        raise ValueError("data.num_reps must be positive.")
    if cfg.data.len_reps <= 0:
        raise ValueError("data.len_reps must be positive.")
    if cfg.data.num_commit <= 0:
        raise ValueError("data.num_commit must be positive.")
    if cfg.data.num_dispatch <= 0:
        raise ValueError("data.num_dispatch must be positive.")
    if cfg.data.duration_dispatch <= 0:
        raise ValueError("data.duration_dispatch must be positive.")

    if cfg.data.representative_dates is not None:
        if len(cfg.data.representative_dates) != cfg.data.num_reps:
            raise ValueError(
                "Length of data.representative_dates must equal data.num_reps. "
                f"Received {len(cfg.data.representative_dates)} dates and "
                f"num_reps={cfg.data.num_reps}."
            )

    if cfg.data.representative_weights is not None:
        if len(cfg.data.representative_weights) != cfg.data.num_reps:
            raise ValueError(
                "Length of data.representative_weights must equal data.num_reps. "
                f"Received {len(cfg.data.representative_weights)} weights and "
                f"num_reps={cfg.data.num_reps}."
            )
        if any(weight < 0 for weight in cfg.data.representative_weights):
            raise ValueError("data.representative_weights must be nonnegative.")
        if sum(cfg.data.representative_weights) <= 0:
            raise ValueError(
                "At least one data.representative_weights value must be positive."
            )

    if cfg.model.gdp_transformation not in SUPPORTED_GDP_TRANSFORMATIONS:
        raise ValueError(
            f"Unsupported GDP transformation: {cfg.model.gdp_transformation}. "
            f"Supported values: {sorted(SUPPORTED_GDP_TRANSFORMATIONS)}"
        )

    ph = cfg.progressive_hedging
    if ph.max_iterations < 0:
        raise ValueError("progressive_hedging.max_iterations must be nonnegative.")
    if ph.min_iterations < 0:
        raise ValueError("progressive_hedging.min_iterations must be nonnegative.")
    if ph.min_iterations > ph.max_iterations:
        raise ValueError(
            "progressive_hedging.min_iterations must be less than or equal to "
            "progressive_hedging.max_iterations."
        )
    if ph.convergence_tolerance < 0:
        raise ValueError(
            "progressive_hedging.convergence_tolerance must be nonnegative."
        )
    if ph.rho <= 0:
        raise ValueError("progressive_hedging.rho must be positive.")
    if ph.regularization not in SUPPORTED_REGULARIZATIONS:
        raise ValueError(
            f"Unsupported PH regularization: {ph.regularization}. "
            f"Supported values: {sorted(SUPPORTED_REGULARIZATIONS)}"
        )

    if cfg.solver.threads is not None and cfg.solver.threads <= 0:
        raise ValueError("solver.threads must be positive when specified.")
    if cfg.solver.mip_gap is not None and cfg.solver.mip_gap < 0:
        raise ValueError("solver.mip_gap must be nonnegative when specified.")
    if cfg.solver.time_limit is not None and cfg.solver.time_limit <= 0:
        raise ValueError("solver.time_limit must be positive when specified.")

    if cfg.solver.license_env_var is not None and not cfg.solver.license_env_var:
        raise ValueError("solver.license_env_var must be non-empty when specified.")

    for env_key, env_value in cfg.solver.environment.items():
        if not env_key:
            raise ValueError(
                "solver.environment contains an empty environment variable name."
            )
        if env_value is None:
            raise ValueError(f"solver.environment[{env_key!r}] must not be null.")

    if cfg.output.binary_rounding_tolerance < 0:
        raise ValueError("output.binary_rounding_tolerance must be nonnegative.")

    unknown_history_formats = (
        set(cfg.output.convergence_history_formats) - SUPPORTED_HISTORY_FORMATS
    )
    if unknown_history_formats:
        raise ValueError(
            "Unsupported output.convergence_history_formats value(s): "
            + ", ".join(sorted(unknown_history_formats))
        )

    if (
        cfg.output.save_operational_variables
        and not cfg.output.save_full_scenario_solution
    ):
        logger.warning(
            "output.save_operational_variables=True while "
            "output.save_full_scenario_solution=False. Operational-value extraction "
            "will be ignored unless full scenario output support is enabled."
        )


def prepare_output_directory(cfg: PHConfig, *, clean: bool | None = None) -> None:
    """Create the PH output directory structure.

    Parameters
    ----------
    cfg:
        PH configuration.
    clean:
        Whether to remove the existing run output directory before creating the
        directory structure. If ``None``, use ``cfg.run.clean_output_dir``.

    Notes
    -----
    Scenario-solve jobs and orchestrator continuations must call this with
    ``clean=False``. Otherwise they can delete PH state and scenario results
    produced earlier in the run.
    """
    root = cfg.run.output_dir

    if clean is None:
        clean = cfg.run.clean_output_dir

    if clean and root.exists():
        shutil.rmtree(root)

    directories = [
        root,
        root / "state",
        root / "iterations",
        root / "convergence",
        root / "final",
        root / "logs",
    ]

    for directory in directories:
        directory.mkdir(parents=True, exist_ok=True)


def write_effective_config(
    cfg: PHConfig, output_path: str | Path | None = None
) -> Path:
    """Write the effective validated configuration to YAML or JSON.

    Parameters
    ----------
    cfg:
        Configuration object.
    output_path:
        Destination path. If omitted, writes ``config_effective.yaml`` under
        ``run.output_dir``.

    Returns
    -------
    pathlib.Path
        Path written.
    """
    if output_path is None:
        output_path = cfg.run.output_dir / "config_effective.yaml"

    path = Path(output_path)
    path.parent.mkdir(parents=True, exist_ok=True)

    data = to_builtin_dict(cfg)

    if path.suffix.lower() == ".json":
        with path.open("w", encoding="utf-8") as f:
            json.dump(data, f, indent=2)
    else:
        if yaml is None:  # pragma: no cover
            raise ImportError(
                "PyYAML is required to write Progressive Hedging YAML configs."
            ) from _YAML_IMPORT_ERROR
        with path.open("w", encoding="utf-8") as f:
            yaml.safe_dump(data, f, sort_keys=False)

    return path


def iteration_dir(cfg: PHConfig, iteration: int) -> Path:
    """Return the output directory for one PH iteration."""
    return cfg.run.output_dir / "iterations" / f"iter_{iteration:03d}"


def scenario_result_path(cfg: PHConfig, iteration: int, scenario_id: int) -> Path:
    """Return the required scenario-result JSON path."""
    return iteration_dir(cfg, iteration) / f"scenario_{scenario_id:03d}.json"


def scenario_log_path(cfg: PHConfig, iteration: int, scenario_id: int) -> Path:
    """Return the default per-scenario solver log path."""
    return (
        cfg.run.output_dir
        / "logs"
        / f"iter_{iteration:03d}_scenario_{scenario_id:03d}.log"
    )


def state_path(cfg: PHConfig, iteration: int) -> Path:
    """Return the PH state JSON path for an iteration."""
    return cfg.run.output_dir / "state" / f"state_iter_{iteration:03d}.json"


def iteration_summary_path(cfg: PHConfig, iteration: int) -> Path:
    """Return the PH iteration summary JSON path."""
    return iteration_dir(cfg, iteration) / "summary.json"


def convergence_history_json_path(cfg: PHConfig) -> Path:
    """Return the convergence-history JSON path."""
    return cfg.run.output_dir / "convergence" / "history.json"


def convergence_history_csv_path(cfg: PHConfig) -> Path:
    """Return the convergence-history CSV path."""
    return cfg.run.output_dir / "convergence" / "history.csv"


def final_solution_dir(cfg: PHConfig) -> Path:
    """Return the final-solution output directory."""
    return cfg.run.output_dir / "final"


def to_builtin_dict(obj: Any) -> Any:
    """Convert dataclasses, Paths, and nested containers to built-in types."""
    if is_dataclass(obj):
        return {
            key: to_builtin_dict(value)
            for key, value in asdict(obj).items()
            if key != "config_path" or value is not None
        }

    if isinstance(obj, Path):
        return str(obj)

    if isinstance(obj, dict):
        return {
            str(to_builtin_dict(key)): to_builtin_dict(value)
            for key, value in obj.items()
        }

    if isinstance(obj, list):
        return [to_builtin_dict(value) for value in obj]

    if isinstance(obj, tuple):
        return [to_builtin_dict(value) for value in obj]

    return obj


def _make_dataclass(cls: type, raw: dict[str, Any]) -> Any:
    """Create dataclass instance from raw mapping with unknown-key validation."""
    if raw is None:
        raw = {}
    if not isinstance(raw, dict):
        raise TypeError(f"{cls.__name__} configuration must be a mapping.")

    kwargs = _filter_dataclass_kwargs(cls, raw)
    return cls(**kwargs)


def _filter_dataclass_kwargs(cls: type, raw: dict[str, Any]) -> dict[str, Any]:
    """Filter mapping to dataclass fields and reject unknown keys."""
    field_names = set(cls.__dataclass_fields__)  # type: ignore[attr-defined]
    unknown = set(raw) - field_names
    if unknown:
        raise ValueError(
            f"Unknown field(s) for {cls.__name__}: " + ", ".join(sorted(unknown))
        )

    return copy.deepcopy(raw)


def _normalize_paths(cfg: PHConfig) -> PHConfig:
    """Normalize path-like configuration values to ``Path`` objects."""
    cfg.run.output_dir = Path(cfg.run.output_dir)

    cfg.data.data_path = Path(cfg.data.data_path)

    if cfg.cost_data.bus_data_path is not None:
        cfg.cost_data.bus_data_path = Path(cfg.cost_data.bus_data_path)
    if cfg.cost_data.cost_data_path is not None:
        cfg.cost_data.cost_data_path = Path(cfg.cost_data.cost_data_path)
    if cfg.cost_data.ng_cost_path is not None:
        cfg.cost_data.ng_cost_path = Path(cfg.cost_data.ng_cost_path)

    if cfg.solver.log_file is not None:
        cfg.solver.log_file = Path(cfg.solver.log_file)
    if cfg.solver.license_file is not None:
        cfg.solver.license_file = Path(cfg.solver.license_file)

    return cfg
