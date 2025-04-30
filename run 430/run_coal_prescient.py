from prescient.simulator import Prescient

# set some options
prescient_options = {
    "data_path": "Prescient",
    "input_format": "rts-gmlc",
    "simulate_out_of_sample": False,
    "run_sced_with_persistent_forecast_errors": False,
    "output_directory": "Prescient/results",
    "start_date": "01-01-2035",
    "num_days": 365,
    "sced_horizon": 24,
    "ruc_mipgap": 0.01,
    "reserve_factor": 0,
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
}
# run the simulator
Prescient().simulate(**prescient_options)