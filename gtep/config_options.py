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
}


def _get_model_config():
    CONFIG = ConfigBlock("GTEPModelConfig")
    CONFIG.declare(
        "flow_model",
        ConfigValue(
            default="DC",
            domain=In(_supported_flows),
            description="Power flow approximation to use.",
        ),
    )
    CONFIG.declare(
        "time_period_dict",
        ConfigDict(
            description="Time period dict, specified as \{(investment period #, length): \{(commitment period #, length): \{dispatch period #: length\}\}\}"
        ),
    )
    CONFIG.declare(
        "dispatch_randomizations",
        ConfigValue(
            default=True,
            domain=Bool,
            description="Introduces random dispatch information rather than having fixed values per-commitment period.",
        ),
    )
    return CONFIG


def _add_common_configs(CONFIG):
    pass


def _add_investment_configs(CONFIG):
    CONFIG.declare(
        "thermal_generation",
        ConfigValue(
            default=False,
            domain=Bool,
            description="Include thermal generation investment options",
        ),
    )
    CONFIG.declare(
        "renewable_generation",
        ConfigValue(
            default=False,
            domain=Bool,
            description="Include renewable generation investment options",
        ),
    )
    CONFIG.declare(
        "storage",
        ConfigValue(
            default=False, domain=Bool, description="Include storage investment options"
        ),
    )
    CONFIG.declare(
        "transmission",
        ConfigValue(
            default=False,
            domain=Bool,
            description="Include transmission investment options",
        ),
    )
    pass


def _add_solver_configs(CONFIG):
    pass
