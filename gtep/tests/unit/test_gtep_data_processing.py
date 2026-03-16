from pathlib import Path
from itertools import product
import pandas as pd
import pyomo.common.unittest as unittest
from pyomo.common.tempfiles import TempfileManager
from gtep.gtep_data_processing import DataProcessing

curr_dir = Path(__file__).resolve().parent
bus_data_path = (
    curr_dir / ".." / ".." / "data" / "costs" / "Bus_data_gen_weights_mappings.csv"
).resolve()
cost_data_path = (
    curr_dir
    / ".."
    / ".."
    / "data"
    / "costs"
    / "2022_v3_Annual_Technology_Baseline_Workbook_Mid-year_update_2-15-2023_Clean.xlsx"
).resolve()


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

    def test_get_gen_bus_data(self):
        """
        Test we are extracting bus data correctly with `get_gen_bus_data`.
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

    def test_get_gen_cost_data_build_cost_data_row(self):
        """
        Test `get_gen_cost_data`, which should merge `gen_bus_df` and `cost_df` on the
        `gen` column of `gen_bus_df` and the `"Key1"` column of `cost_df`. The result
        should also have a MultiIndex with names `"Key1"` and `"Bus Name"`.

        Also test building a row in the output dataframe using `build_cost_data_row`.
        Takes output of `get_gen_cost_data` as input, so need to be tested together.
        """
        gens = ["Natural Gas_FE", "Solar - Utility PV"]
        years = [2025, 2030]
        ng_costs = {2025: 3.49, 2030: 2.91}

        for gen in gens:
            ### TEST GET_GEN_COST_DATA ###
            gen_bus_df = pd.DataFrame(
                data=[
                    ["bus1", "genname"],
                    ["bus2", "genname"],
                    ["bus3", "genname"],
                ],
                columns=["Bus Name", gen],
            )
            varnames = list(self.data_processing.cost_var_names.values())
            if "Natural Gas" in gen:
                varnames += [self.data_processing.heat_rate_var]
            cost_df = pd.DataFrame(
                data=[[var, "genname", 1.0, 2.0] for var in varnames],
                columns=["Key1", "Key2", *years],
            )
            gen_cost_df = self.data_processing.get_gen_cost_data(
                cost_df, gen_bus_df, gen
            )

            self.assertIsInstance(gen_cost_df, pd.DataFrame)
            self.assertTupleEqual(
                gen_cost_df.shape,
                (len(gen_bus_df) * len(varnames), len(years)),
            )
            self.assertIsInstance(gen_cost_df.index, pd.MultiIndex)
            self.assertEqual(len(gen_cost_df.index.names), 2)
            self.assertEqual(gen_cost_df.index.names[0], "Key1")
            self.assertEqual(gen_cost_df.index.names[1], "Bus Name")
            for year in years:
                self.assertIn(year, gen_cost_df.columns)

            ### TEST BUILD_COST_DATA_ROW ###
            expected_keys = [
                "GEN UID",
                "Bus ID",
                "Unit Type",
                "Fuel",
                "PMax MW",
                "PMin MW",
                "Min Up Time Hr",
                "Min Down Time Hr",
            ]
            out_varnames = list(self.data_processing.cost_var_names.keys()) + [
                "fuel_costs"
            ]
            expected_keys += [
                f"{var}_{year}" for var, year in product(out_varnames, years)
            ]

            row = self.data_processing.build_cost_data_row(
                gen_cost_df=gen_cost_df,
                years=years,
                bus_id=0,
                bus_name="bus1",
                gen=gen,
                ng_costs=ng_costs,
            )
            self.assertIsInstance(row, dict)
            # check keys are exactly what we expect
            self.assertEqual(len(row), len(expected_keys))
            for col in expected_keys:
                self.assertIn(col, row)
            # TODO: possibly check values...?

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

    def test_extract_cost_data_build_cost_data_row(self):
        """
        Test that cost data is being properly extracted (`extract cost data`)
        and that we build rows for the output dataframe properly from that
        data (`build_cost_data_row`). The latter requires the cost dataframe
        from the former, so they must be tested together.
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
            self.data_processing.load_gen_data(
                bus_data_path,
                cost_data_path,
                gens,
                save_csv=True,
                write_dir=tempdir,
            )
            self.assertIsInstance(
                self.data_processing.gen_data_target, pd.DataFrame
            )  # make sure we are still storing as an attribute
            self.assertIn("costs.csv", Path.iterdir(tempdir))
