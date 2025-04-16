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
        "include_investment",
        ConfigValue(
            default=True,
            domain=Bool,
            description="Enable inclusion of any investment options.",
        ),
    )

    CONFIG.declare(
        "include_commitment",
        ConfigValue(
            default=True,
            domain=Bool,
            description="Include unit commitment formulation.",
        ),
    )

    CONFIG.declare(
        "include_redispatch",
        ConfigValue(
            default=True,
            domain=Bool,
            description="Include economic redispatch formulation (i.e., >1 dispatch period per commitment period).",
        ),
    )

    CONFIG.declare(
        "flow_model",
        ConfigValue(
            default="DC",
            domain=In(_supported_flows),
            description="Power flow approximation to use.",
        ),
    )

    CONFIG.declare(
        "time_period_subsets",
        ConfigList(
            description="Time period counts for fixed-length and fixed-subset periods."
        ),
    )

    CONFIG.declare(
        "time_period_dict",
        ConfigDict(
            description="Time period dict, specified as \{(investment period #, length): \{(representative period #, length): \{(commitment period #, length): \{dispatch period #: length\}\}\}"
        ),
    )

    CONFIG.declare(
        "dispatch_randomization",
        ConfigValue(
            default=True,
            domain=Bool,
            description="Introduces random dispatch information rather than having fixed values per commitment period.",
        ),
    )
    return CONFIG


def _add_common_configs(CONFIG):

    CONFIG.declare(
        "scale_loads",
        ConfigValue(
            default=True,
            domain=Bool,
            description="Allow scaling of load values into future years; i.e., load scaling is represented in the model but not the data.",
        ),
    )

    CONFIG.declare("scale_texas_loads", ConfigValue(default=False, domain=Bool, description = "but why"))


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
    CONFIG.declare("transmission_switching", ConfigValue(default=False, domain=Bool, description="Allow transmission switching during dispatch"))


def _add_solver_configs(CONFIG):
    pass
