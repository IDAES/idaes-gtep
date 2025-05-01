from prescient.simulator import Prescient

# set some options
prescient_options = {
    "data_path": "Prescient",
    "input_format": "rts-gmlc",
    "simulate_out_of_sample": True,
    "run_sced_with_persistent_forecast_errors": True,
    "output_directory": "Prescient/results",
    "start_date": "08-10-2035",
    "num_days": 3,
    "sced_horizon": 1,
    "ruc_mipgap": 0.01,
    "deterministic_ruc_solver": "gurobi_persistent",
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
    "reserve_factor": 0,
    "day_ahead_pricing": "LMP",
    "print_sced": True
}
# run the simulator
Prescient().simulate(**prescient_options)