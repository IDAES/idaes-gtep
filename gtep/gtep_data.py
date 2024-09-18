# Generation and Transmission Expansion Planning
# IDAES project
# author: Kyle Skolfield
# date: 01/04/2024
# Model available at http://www.optimization-online.org/DB_FILE/2017/08/6162.pdf

from pyomo.environ import *
from prescient.simulator.config import PrescientConfig
from prescient.data.providers import gmlc_data_provider
import datetime


class ExpansionPlanningData:
    """Standard data storage class for the IDAES GTEP model."""

    def __init__(self):
        pass

    def load_prescient(self, data_path, options_dict=None):
        """Loads data structured via Prescient data loader.

        :param data_path: Folder containing the data to be loaded
        :param options_dict: Options dictionary to pass to the Prescient data loader, defaults to None
        """
        self.data_type = "prescient"
        options_dict = {
            "data_path": data_path,
            "input_format": "rts-gmlc",
            "start_date": "01-01-2020",
            "num_days": 365,
            "sced_horizon": 1,
            "sced_frequency_minutes": 60,
            "ruc_horizon": 36,
        }

        prescient_options = PrescientConfig()
        prescient_options.set_value(options_dict)
        # Use prescient data provider to load in sequential data for representative periods
        data_list = []

        x = datetime.datetime(2020, 1, 1)
        data_provider = gmlc_data_provider.GmlcDataProvider(options=prescient_options)
        # populate an egret model data with the basic stuff
        self.md = data_provider.get_initial_actuals_model(
            options=prescient_options, num_time_steps=24 * 365, minutes_per_timestep=60
        )
        # fill in renewable actuals and maybe demand idk
        data_provider.populate_with_actuals(
            options=prescient_options,
            num_time_periods=24 * 365,
            time_period_length_minutes=60,
            start_time=x,
            model=self.md,
        )
        # data_provider.populate_initial_state_data(options=prescient_options, model=md)

        self.load_default_data_settings()

        for gen in self.md.data["elements"]["generator"]:
            if "-c" in gen:  # key/Gen UID in csv file; -c = candidate?
                self.md.data["elements"]["generator"][gen]["in_service"] = False

        # JSC addn
        for branch in self.md.data["elements"]["branch"]:
            if "-c" in branch:  # key/Branch UID
                self.md.data["elements"]["branch"][branch]["in_service"] = False

        ## NOTE: Below is only for multiple representative periods and creates a list
        ## of modelData objects, not just a single modelData object
        # Arbitrary time points and lengths picked for representative periods
        # default here allows up to 24 hours for periods
        time_keys = self.md.data["system"]["time_keys"]
        key_idx = time_keys.index("2020-01-01 00:00")
        time_key_set = time_keys[key_idx : key_idx + 24]
        data_list.append(self.md.clone_at_time_keys(time_key_set))
        key_idx = time_keys.index("2020-04-01 00:00")
        time_key_set = time_keys[key_idx : key_idx + 24]
        data_list.append(self.md.clone_at_time_keys(time_key_set))
        key_idx = time_keys.index("2020-07-01 00:00")
        time_key_set = time_keys[key_idx : key_idx + 24]
        data_list.append(self.md.clone_at_time_keys(time_key_set))
        key_idx = time_keys.index("2020-10-01 00:00")
        time_key_set = time_keys[key_idx : key_idx + 24]
        data_list.append(self.md.clone_at_time_keys(time_key_set))

        self.representative_data = data_list

    def load_default_data_settings(self):
        ## TODO: too many of these are hard coded; everything should check if it exists too.
        """Fills in necessary but unspecified data information."""
        for gen in self.md.data["elements"]["generator"]:
            self.md.data["elements"]["generator"][gen]["lifetime"] = 3
            self.md.data["elements"]["generator"][gen]["spinning_reserve_frac"] = 0.1
            self.md.data["elements"]["generator"][gen]["quickstart_reserve_frac"] = 0.1
            self.md.data["elements"]["generator"][gen]["capital_multiplier"] = 1
            self.md.data["elements"]["generator"][gen]["extension_multiplier"] = 1
            self.md.data["elements"]["generator"][gen]["max_operating_reserve"] = 1
            self.md.data["elements"]["generator"][gen]["max_spinning_reserve"] = 1
            self.md.data["elements"]["generator"][gen]["max_quickstart_reserve"] = 1
            self.md.data["elements"]["generator"][gen]["ramp_up_rate"] = 0.1
            self.md.data["elements"]["generator"][gen]["ramp_down_rate"] = 0.1
            self.md.data["elements"]["generator"][gen]["emissions_factor"] = 1
            self.md.data["elements"]["generator"][gen]["start_fuel"] = 1
            self.md.data["elements"]["generator"][gen]["investment_cost"] = 235164
        for branch in self.md.data["elements"]["branch"]:
            self.md.data["elements"]["branch"][branch]["loss_rate"] = 0
            self.md.data["elements"]["branch"][branch]["distance"] = 1
        self.md.data["system"]["min_operating_reserve"] = 0.1
        self.md.data["system"]["min_spinning_reserve"] = 0.1
