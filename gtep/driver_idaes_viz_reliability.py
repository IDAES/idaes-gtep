from gtep.gtep_model import ExpansionPlanningModel
from gtep.gtep_data import ExpansionPlanningData
from gtep.gtep_solution import ExpansionPlanningSolution
from pyomo.core import TransformationFactory
from pyomo.contrib.appsi.solvers.highs import Highs
from pyomo.contrib.appsi.solvers.gurobi import Gurobi
import gtep.validation


data_path = "./gtep/data/SanDiego"
data_object = ExpansionPlanningData()
data_object.load_prescient(data_path)

sol_object = ExpansionPlanningSolution()
sol_object.import_data_object(data_object)
sol_object.read_json("C:/Users/jkskolf/idaes-gtep/gtep_pre_reliability_solution.json")


data_out_path = "./Sim Eng Viz/PreReliability/Prescient/data"

gtep.validation.populate_generators(data_path, sol_object, data_out_path)
gtep.validation.populate_transmission(data_path, sol_object, data_out_path)
gtep.validation.filter_pointers(data_path, data_out_path)
gtep.validation.clone_timeseries(data_path, data_out_path)


# sol_object.import_data_object(data_object)

# sol_object.plot_levels(save_dir="./Sim Eng Viz/Hourly/GTEP/plots/")

from prescient.simulator import Prescient

# set some options
prescient_options = {
    "data_path": data_out_path,
    "input_format": "rts-gmlc",
    "simulate_out_of_sample": False,
    "run_sced_with_persistent_forecast_errors": False,
    "output_directory": "./Sim Eng Viz/PreReliability/Prescient/results",
    "start_date": "07-26-2020",
    "num_days": 21,
    "sced_horizon": 24,
    "ruc_mipgap": 0.01,
    "reserve_factor": 0,
    "deterministic_ruc_solver": "gurobi_persistent",
    "deterministic_ruc_solver_options": {
        "feas": "off",
        "DivingF": "on",
    },
    "sced_solver": "gurobi",
    "sced_frequency_minutes": 60,
    "ruc_horizon": 48,
    "compute_market_settlements": True,
    "monitor_all_contingencies": False,
    "output_solver_logs": False,
    "price_threshold": 1000,
    "contingency_price_threshold": 100,
    "reserve_price_threshold": 5,
    "ruc_network_type": "btheta",
    "sced_network_type": "btheta",
    "enforce_sced_shutdown_ramprate": False,
}
# run the simulator
Prescient().simulate(**prescient_options)
