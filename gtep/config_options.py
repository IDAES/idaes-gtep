from pyomo.common.config import (
    ConfigBlock,
    ConfigDict,
    ConfigList,
    ConfigValue,
    In,
    NonNegativeFloat,
    NonNegativeInt,
    PositiveInt,
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
        "dispatch_period_dict",
        ConfigDict(
            description="Dispatch period dict, specified as \{dispatch period #: (parent commitment period, length)\}"
        ),
    )
    return CONFIG


def _add_common_configs(CONFIG):
    pass


def _add_investment_configs(CONFIG):
    pass


def _add_solver_configs(CONFIG):
    pass
