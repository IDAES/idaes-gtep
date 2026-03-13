from os.path import dirname, abspath, join
from os import listdir
import pandas as pd
import pyomo.common.unittest as unittest
from pyomo.common.tempfiles import TempfileManager
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


class TestGTEPDataProcessing(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.data_processing = DataProcessing()

    def test_data_processing_init(self):
        """
        Check we are initializing properly.
        """
        self.assertIsInstance(self.data_processing, DataProcessing)
        # have to do this for now bc self.assertHasAttr only available in python>=3.14
        self.assertTrue(hasattr(self.data_processing, "config"))

    def test_get_gen_bus_data_get_buses_by_gen(self):
        """
        Test `get_gen_bus_data` and `get_buses_by_gen`. These have to be
        tested together because `get_buses_by_gen` requires a dataframe
        as input, which `get_gen_bus_data` provides.
        """
        gens = ["Natural Gas_CT", "Natural Gas_FE", "Solar - Utility PV"]
        gen_data = self.data_processing.get_gen_bus_data(bus_data_path)
        self.assertIsInstance(gen_data, pd.DataFrame)

        buses_by_gen = self.data_processing.get_buses_by_gen(gen_data, gens)
        for gen, bus_data in buses_by_gen.items():
            self.assertIsInstance(gen, str)
            self.assertIn(gen, gens)
            self.assertIsInstance(bus_data, dict)
            for bus_name, bus_id in bus_data.items():
                self.assertIsInstance(bus_name, str)
                self.assertIsInstance(bus_id, int)

    def test_get_clean_gens_dict(self):
        """
        Test that we are correctly mapping a set of input generators
        to a set that has cost data.
        """
        gens = ["Natural Gas_CT", "Natural Gas_FE", "Solar - Utility PV"]
        # with self.assertWarns(UserWarning):
        #     gens_dict = self.data_processing.get_clean_gens_dict(gens)
        out_dict = self.data_processing.get_clean_gens_dict(gens)
        self.assertIsInstance(out_dict, dict)

        # make sure each gen we provided is a key in the dict, and test
        # proper behavior for corresponding value
        for gen in gens:
            self.assertIn(gen, out_dict)
            if gen == "Natural Gas_CT":
                self.assertNotIn(gen, out_dict.values())
            else:
                self.assertIn(gen, out_dict.values())

        # make sure each key in the output is a member of the input
        for key in out_dict:
            self.assertIsInstance(key, str)
            self.assertIn(key, gens)

    def test_extract_cost_data(self):
        """
        Test that cost data is being properly extracted.
        """
        gens = ["Natural Gas_FE", "Solar - Utility PV"]
        years = [2025, 2030]
        scenario = "Moderate"
        cost_data = self.data_processing.extract_cost_data(
            cost_data_path, gens, years, scenario
        )

        # check that all the gens we provided are in the result
        for provided_gen in gens:
            self.assertIn(provided_gen, cost_data)

        # check that all keys of cost_data are gens we provided, and that the values
        # are dataframes
        for gen, gen_cost_data in cost_data.items():
            self.assertIsInstance(gen, str)
            self.assertIn(gen, gens)
            self.assertIsInstance(gen_cost_data, pd.DataFrame)

    def test_build_cost_data_row(self):
        pass

    def test_fill_out_prescient_columns(self):
        df = pd.DataFrame([[0, 0], [0, 0]])
        df = self.data_processing.fill_out_prescient_columns(df)

        # make sure all columns are properly added
        for col in self.data_processing.prescient_cols:
            self.assertIn(col, df.columns)

        # make sure default values for specific columns are included
        for val in df["V Setpoint p.u."]:
            self.assertEqual(val, 1)
        for val in df["Output_pct_0"]:
            self.assertEqual(val, 0.6)
        for val in df["Output_pct_1"]:
            self.assertEqual(val, 1)

        # make sure that no NAs remain
        self.assertEqual(df.isna().sum().sum(), 0)

    def test_load_gen_data(self):
        gens = ["Natural Gas_CT", "Natural Gas_FE", "Solar - Utility PV"]

        # basic case with all default arguments
        self.data_processing.load_gen_data(
            bus_data_path,
            cost_data_path,
            gens,
        )
        self.assertIsInstance(self.data_processing.gen_data_target, pd.DataFrame)
        self.gen_data_target = None  # reset

        # change ng costs so that we are missing data for a year
        with self.assertRaises(KeyError):
            self.data_processing.load_gen_data(
                bus_data_path,
                cost_data_path,
                gens,
                ng_costs={2025: 3.49, 2030: 2.91},
            )

        # try using save_csv=True without specifying output directory
        with self.assertRaises(TypeError):
            self.data_processing.load_gen_data(
                bus_data_path,
                cost_data_path,
                gens,
                save_csv=True,
            )

        # ensure file is written where expected
        with TempfileManager.new_context() as tempfile:
            tempdir = tempfile.mkdtemp()
            out_path = join(tempdir, "costs.csv")
            self.data_processing.load_gen_data(
                bus_data_path,
                cost_data_path,
                gens,
                save_csv=True,
                out_path=out_path,
            )
            self.assertIsInstance(
                self.data_processing.gen_data_target, pd.DataFrame
            )  # make sure we are still storing as an attribute
            self.assertIn("costs.csv", listdir(tempdir))
