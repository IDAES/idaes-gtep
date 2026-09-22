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

from pyomo.common.config import (
    ConfigBlock,
    ConfigDict,
    ConfigList,
    ConfigValue,
    In,
    NonNegativeFloat,
    NonNegativeInt,
    PositiveInt,
    Bool,
)
from pyomo.common.deprecation import deprecation_warning

_supported_flows = {
    "DC": ("gtep.dcopf", "DC power flow approximation"),
    "CP": ("gtep.cp", "Copper plate power flow approximation"),
    "ACP": ("gtep.acp", "AC power flow in polar formulation"),
    "ACR": ("gtep.acr", "AC power flow in rectangular formulation"),
    "transport": ("gtep.transport", "transport"),
}


def _get_model_config():
    """This function creates and returns the base GTEP model
    configuration block.

    This configuration block defines core model options used by the
    ``ExpansionPlanningModel`` class. These options control investment,
    commitment, redispatch, power-flow formulation, and time-period
    structure settings.

    :return: GTEP model configuration block.

    """

    CONFIG = ConfigBlock("GTEPModelConfig")

    CONFIG.declare(
        "include_investment",
        ConfigValue(
            default=True,
            domain=Bool,
            description=(
                "Enable investment decisions for candidate and existing "
                "assets. When disabled, candidate assets are not selected "
                "for installation."
            ),
        ),
    )

    CONFIG.declare(
        "include_commitment",
        ConfigValue(
            default=True,
            domain=Bool,
            description=(
                "Include unit commitment decisions, including generator "
                "on/off operating-status logic."
            ),
        ),
    )

    CONFIG.declare(
        "include_redispatch",
        ConfigValue(
            default=True,
            domain=Bool,
            description=(
                "Include redispatch decisions within commitment periods. "
                "This is relevant when there is more than one dispatch "
                "period per commitment period."
            ),
        ),
    )

    CONFIG.declare(
        "flow_model",
        ConfigValue(
            default="DC",
            domain=In(_supported_flows),
            description=(
                "Power-flow formulation to use. Supported options include "
                "DC, CP, ACP, ACR, and transport."
            ),
        ),
    )

    CONFIG.declare(
        "time_period_subsets",
        ConfigList(
            description=(
                "Optional list defining fixed-length or fixed-subset "
                "time-period structures."
            )
        ),
    )

    CONFIG.declare(
        "time_period_dict",
        ConfigDict(
            description=(
                "Optional nested dictionary defining custom investment, "
                "representative, commitment, and dispatch period structures."
            )
        ),
    )

    CONFIG.declare(
        "dispatch_randomization",
        ConfigValue(
            default=True,
            domain=Bool,
            description=(
                "Use randomized dispatch information instead of fixed "
                "values per commitment period."
            ),
        ),
    )

    return CONFIG


def _add_common_configs(CONFIG):
    """Add common GTEP model configuration options.

    These options are shared across model formulations and control
    load-scaling behavior.

    """

    CONFIG.declare(
        "scale_loads",
        ConfigValue(
            default=True,
            domain=Bool,
            description=(
                "Enable scaling of load values into future years. Load "
                "scaling is represented in the model rather than directly "
                "modifying the input data."
            ),
        ),
    )

    CONFIG.declare(
        "scale_texas_loads",
        ConfigValue(
            default=False,
            domain=Bool,
            description=(
                "Enable Texas-case-specific load scaling logic, when " "applicable."
            ),
        ),
    )


def _add_investment_configs(CONFIG):
    """This function adds investment and model-component configuration
    options.

    These options control which candidate asset types and component
    formulations are included in the GTEP model.

    """

    CONFIG.declare(
        "thermal_generation",
        ConfigValue(
            default=False,
            domain=Bool,
            description="Include thermal generation investment options.",
        ),
    )

    CONFIG.declare(
        "renewable_generation",
        ConfigValue(
            default=False,
            domain=Bool,
            description="Include renewable generation investment options.",
        ),
    )

    CONFIG.declare(
        "storage",
        ConfigValue(
            default=False,
            domain=Bool,
            description="Include storage investment and operation modeling.",
        ),
    )

    CONFIG.declare(
        "transmission",
        ConfigValue(
            default=False,
            domain=Bool,
            description="Include transmission modeling and investment options.",
        ),
    )

    CONFIG.declare(
        "transmission_switching",
        ConfigValue(
            default=False,
            domain=Bool,
            description="Allow transmission switching decisions during dispatch.",
        ),
    )

    CONFIG.declare(
        "advanced_hydro",
        ConfigValue(
            default=False,
            domain=Bool,
            description=(
                "Include advanced hydro modeling features, including daily "
                "average hydro requirements."
            ),
        ),
    )


def _add_solver_configs(CONFIG):
    """This function adds solver-related configuration options.

    This is currently reserved for future solver options.  Solver
    settings are handled by the driver/configuration workflow.

    """
    pass
