from os.path import dirname, abspath, join
import pandas as pd
import pyomo.common.unittest as unittest
from gtep.gtep_data_processing import DataProcessing

curr_dir = dirname(abspath(__file__))
bus_data_path = abspath(
    join(curr_dir, "..", "..", "data", "costs", "Bus_data_gen_weights_mappings.csv")
)
cost_data_path = abspath(
    join(
        curr_dir,
        "..",
        "..",
        "data",
        "costs",
        "2022_v3_Annual_Technology_Baseline_Workbook_Mid-year_update_2-15-2023_Clean.xlsx",
    )
)
out_path = abspath(
    join(
        curr_dir,
        "..",
        "..",
        "data",
        "costs",
    )
)

candidate_gens = [
    "Natural Gas_CT",
    "Natural Gas_FE",
    "Solar - Utility PV",
    "Land-Based Wind",
]


class TestGTEPDataProcessing(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.data_processing = DataProcessing()

    def test_data_processing_init(self):
        self.assertIsInstance(self.data_processing, DataProcessing)
        # have to do this for now bc self.assertHasAttr only available in python>=3.14
        self.assertTrue(hasattr(self.data_processing, "config"))

    def _helper_verify_gen_dict(self, test_dict):
        """
        Verifies that an input, which we expect to be a dict
        having its keys as elements from DataProcessing.all_gens,
        indeed follows that structure.
        """
        self.assertIsInstance(test_dict, dict)
        for gen in test_dict:
            self.assertIsInstance(gen, str)
            # self.assertIn(gen, self.data_processing.all_gens)

    def test_get_gen_bus_data_get_buses_by_gen(self):
        # these need to be tested together because get_buses_by_gen requires
        # a dataframe input, which get_gen_bus_data provides
        gen_data = self.data_processing.get_gen_bus_data(bus_data_path)
        self.assertIsInstance(gen_data, pd.DataFrame)

        buses_by_gen = self.data_processing.get_buses_by_gen(gen_data, candidate_gens)
        self._helper_verify_gen_dict(buses_by_gen)

        for bus_data in buses_by_gen.values():
            self.assertIsInstance(bus_data, dict)
            for key, val in bus_data.items():
                self.assertIsInstance(key, str)
                self.assertIsInstance(val, int)

    # def test_clean_up_candidate_gens_and_get_cost_data(self):
    #     gens_of_interest = self.data_processing.clean_up_candidate_gens(candidate_gens)
    #     self.assertIsInstance(gens_of_interest, set)
    #     self.assertNotIn("Natural Gas_CT", gens_of_interest)

    #     cost_data = self.data_processing.get_cost_data(cost_data_path, gens_of_interest)
    #     self._helper_verify_gen_dict(cost_data)
    #     self.assertNotIn("Natural Gas_CT", cost_data)
    #     for gen_cost_data in cost_data.values():
    #         self.assertIsInstance(gen_cost_data, pd.DataFrame)

    def test_load_gen_data(self):
        self.data_processing.load_gen_data(
            bus_data_path,
            cost_data_path,
            candidate_gens,
            save_csv=True,
            out_path=out_path,
        )
