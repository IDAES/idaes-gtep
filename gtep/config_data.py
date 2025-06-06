# Configuration file for relevant generator data in the Generation and
# Transmission Expansion Planning model

# author: Soraya Rawlings and Kyle Skolfield

class ConfigData:

    # Generator data
    config_data = {
        "Natural Gas_CT": {
            "ids": "ct_",
            "unit_type": "CT",
            "fuel_type": "G",
            "gen_weight": 100,
            "pmax": 727,
            "pmin": 218,
            "minup": 6,
            "mindown": 8,
            "up_time": 6,
            "down_time": 8,
        },
        "Natural Gas_FE": {
            "ids": "ct_fe",
            "unit_type": "CT",
            "fuel_type": "G",
            "gen_weight": 300,
            "pmax": 992,
            "pmin": 297.6,
            "minup": 6,
            "mindown": 8,
            "up_time": 6,
            "down_time": 8,
        },
        "Solar - Utility PV": {
            "ids": "pv_",
            "unit_type": "PV",
            "fuel_type": "S",
            "gen_weight": 180,
            "pmax": 100,
            "pmin": 0,
            "minup": 0,
            "mindown": 0,
            "up_time": 0,
            "down_time": 0,
        },
        "Land-Based Wind": {
            "ids": "wind_",
            "unit_type": "WIND",
            "fuel_type": "W",
            "gen_weight": 150,
            "pmax": 200,
            "pmin": 0,
            "minup": 0,
            "mindown": 0,
            "up_time": 0,
            "down_time": 0,
        },
        # Add generators with up and down times. Use the list
        # "all_gens" from gtep_data_processing as a reference for the
        # generator name.
        "Nuclear": {
            "unit_type": "NUC",
            "fuel_type": "N",
            "up_time": 48,
            "down_time": 48,
        },
        "Coal_FE": {
            "unit_type": "COAL",
            "fuel_type": "C",
            "up_time": 24,
            "down_time": 12,
        },
        "Hydro Dummy 1": {
            "unit_type": "HYDRO",
            "fuel_type": "H",
            "up_time": 0,
            "down_time": 0,
        }
    }


    def __init__(self):

        # Store the dictionary
        self.gen_data = self.config_data
