from pathlib import Path
import tempfile
import shutil
from itertools import product
import pandas as pd
import pyomo.common.unittest as unittest
from pyomo.common.tempfiles import TempfileManager
from gtep.gtep_data_processing import DataProcessing


curr_dir = Path(__file__).resolve().parent
Bus5_case_path = (curr_dir / ".." / ".." / "data" / "5bus").resolve()


class TestGTEPDataProcessing(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.data_processing = DataProcessing()
        cls.temp_dir = Path(tempfile.mkdtemp())
        # Set default years for testing
        cls.data_processing.years = [2025, 2030]

        # Create a sample data object structure that would be used in
        # the actual class.
        cls.data_object = type(
            "DataObject",
            (),
            {
                "years": [2025, 2030],
                "md": type(
                    "MD",
                    (),
                    {
                        "data": {
                            "elements": {
                                "generator": {
                                    "gen1": {"investment_cost": 1000, "fuel_cost": 50},
                                    "gen2": {"investment_cost": 1500, "fuel_cost": 75},
                                },
                                "branch": {"branch1": {"investment_cost": 500}},
                                "storage": {"storage1": {"investment_cost": 2000}},
                            }
                        }
                    },
                )(),
            },
        )()

    @classmethod
    def tearDownClass(cls):
        # Clean up temporary directory
        shutil.rmtree(cls.temp_dir)

    def setUp(self):
        # Create a temporary directory for each test
        self.test_dir = Path(tempfile.mkdtemp())

    def tearDown(self):
        # Clean up test directory
        shutil.rmtree(self.test_dir)

    def test_data_processing_init(self):
        """
        Check we are initializing properly.
        """
        self.assertIsInstance(self.data_processing, DataProcessing)
        # have to do this for now bc self.assertHasAttr only available in python>=3.14
        self.assertTrue(hasattr(self.data_processing, "config"))
        self.assertTrue(hasattr(self.data_processing, "cost_dicts"))
        self.assertIsInstance(self.data_processing.cost_dicts, dict)

    def test_load_and_check_cost_csv_valid_file(self):
        """Test loading and validation of cost data from CSV file."""

        for asset_type in ["generator", "branch"]:

            if asset_type == "generator":
                # Path to the real CSV file
                csv_path = Path(f"{Bus5_case_path}/gen.csv")
                prefixes = [
                    "lifetime_",
                    "capex_",
                    "fixed_ops_cost_",
                    "var_ops_cost_",
                    "fuel_cost_",
                ]
            elif asset_type == "branch":
                csv_path = Path(f"{Bus5_case_path}/branch.csv")
                prefixes = ["lifetime_", "capex_", "fixed_ops_cost_", "var_ops_cost_"]
            else:
                csv_path = Path(f"{Bus5_case_path}/storage.csv")
                prefixes = ["lifetime_", "capex_", "fixed_ops_cost_", "var_ops_cost_"]

            # Call the method directly
            self.data_processing.load_and_check_cost_csv(
                path=csv_path,
                prefixes=prefixes,
                cost_attr=f"{asset_type}_cost_data",
                required=True,
            )

            # Check if cost data was loaded for this asset type
            cost_data_attr = f"{asset_type}_cost_data"
            self.assertTrue(hasattr(self.data_processing, cost_data_attr))
            self.assertIsInstance(
                getattr(self.data_processing, cost_data_attr), pd.DataFrame
            )
            self.assertGreater(len(getattr(self.data_processing, cost_data_attr)), 0)

            # Verify cost dictionaries were populated
            self.assertIn(asset_type, self.data_processing.cost_dicts)
            self.assertIn("capex", self.data_processing.cost_dicts[asset_type])
            self.assertIn("fixed_ops_cost", self.data_processing.cost_dicts[asset_type])
            self.assertIn("var_ops_cost", self.data_processing.cost_dicts[asset_type])

            # Fuel cost only applies to generators
            if asset_type == "generator":
                self.assertIn("fuel_cost", self.data_processing.cost_dicts[asset_type])

        # Check specific cost values from the CSV
        self.assertEqual(
            self.data_processing.cost_dicts["generator"]["capex"][("3_CT", 2025)],
            872153.867247609,
        )
        self.assertEqual(
            self.data_processing.cost_dicts["generator"]["fuel_cost"][("3_CT", 2025)],
            27.97755600459,
        )
        self.assertEqual(
            self.data_processing.cost_dicts["generator"]["fixed_ops_cost"][
                ("3_CT", 2025)
            ],
            2.44467259867671,
        )
        self.assertEqual(
            self.data_processing.cost_dicts["generator"]["var_ops_cost"][
                ("3_CT", 2025)
            ],
            5.1,
        )
        self.assertEqual(
            self.data_processing.cost_dicts["branch"]["capex"][("branch_2_3", 2025)], 0
        )
        self.assertEqual(
            self.data_processing.cost_dicts["branch"]["fixed_ops_cost"][
                ("branch_2_3", 2025)
            ],
            0,
        )
        self.assertEqual(
            self.data_processing.cost_dicts["branch"]["var_ops_cost"][
                ("branch_2_3", 2025)
            ],
            0,
        )
        # self.assertEqual(
        #     self.data_processing.cost_dicts["storage"]["lifetime"][("battery_1a", 2025)],
        #     3
        # )
        # self.assertEqual(
        #     self.data_processing.cost_dicts["storage"]["capex"][("battery_1a", 2025)],
        #     14.91
        # )
        # self.assertEqual(
        #     self.data_processing.cost_dicts["storage"]["fixed_ops_cost"][("battery_1a", 2025)],
        #     0
        # )
        # self.assertEqual(
        #     self.data_processing.cost_dicts["storage"]["var_ops_cost"][("battery_1a", 2025)],
        #     0
        # )
