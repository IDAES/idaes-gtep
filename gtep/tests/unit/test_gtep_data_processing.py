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

    def test_get_gen_nonsense_and_get_candidate_gen_by_bus(self):
        gen_nonsense = self.data_processing.get_gen_nonsense(bus_data_path)
        self.assertIsInstance(gen_nonsense, pd.DataFrame)
        self.data_processing.get_candidate_gen_by_bus(gen_nonsense, candidate_gens)

    def test_load_gen_data(self):
        self.data_processing.load_gen_data(
            bus_data_path, cost_data_path, candidate_gens
        )
