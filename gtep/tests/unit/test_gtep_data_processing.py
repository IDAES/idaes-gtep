from pathlib import Path
import tempfile
import shutil
from itertools import product
import pandas as pd

import pyomo.common.unittest as unittest
from pyomo.common.tempfiles import TempfileManager

from gtep.gtep_data import ExpansionPlanningData
from gtep.gtep_data_processing import DataProcessing

curr_dir = Path(__file__).resolve().parent
Bus5_case_path = (curr_dir / ".." / ".." / "data" / "5bus").resolve()

class TestGTEPDataProcessing(unittest.TestCase):
    @classmethod
    def setUpClass(self):
        self.data_processing = DataProcessing()
        self.years = [2025, 2030]

    def test_data_processing_init(self):
        """Check we are initializing properly."""

        self.assertIsInstance(self.data_processing, DataProcessing)

        # Have to do this for now bc self.assertHasAttr only available
        # in python>=3.14.
        self.assertTrue(hasattr(self.data_processing, "config"))
        self.assertTrue(hasattr(self.data_processing, "cost_dicts"))
        self.assertIsInstance(self.data_processing.cost_dicts, dict)

    def test_load_and_check_cost_csv_valid_file(self):
        """Test loading and validation of cost data from CSV file."""

        for asset_type in ["generator", "branch", "storage"]:

            if asset_type == "generator":
                # Path to data CSV file.
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

            # Call the load_and_check_cost_csv method directly to test
            # it.
            self.data_processing.load_and_check_cost_csv(
                path=csv_path,
                prefixes=prefixes,
                cost_attr=f"{asset_type}_cost_data",
                required=True,
            )

            # Check if cost data is loaded on this asset type.
            cost_data_attr = f"{asset_type}_cost_data"
            self.assertTrue(hasattr(self.data_processing, cost_data_attr))
            self.assertIsInstance(
                getattr(self.data_processing, cost_data_attr), pd.DataFrame
            )
            self.assertGreater(len(getattr(self.data_processing, cost_data_attr)), 0)

            # Verify cost dictionaries were populated.
            self.assertIn(asset_type, self.data_processing.cost_dicts)
            self.assertIn("capex", self.data_processing.cost_dicts[asset_type])
            self.assertIn("fixed_ops_cost", self.data_processing.cost_dicts[asset_type])
            self.assertIn("var_ops_cost", self.data_processing.cost_dicts[asset_type])

            # Fuel cost only applies to generators.
            if asset_type == "generator":
                self.assertIn("fuel_cost", self.data_processing.cost_dicts[asset_type])

        # Check specific cost values from the CSV.
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
        self.assertEqual(
            self.data_processing.cost_dicts["storage"]["lifetime"][
                ("battery_1a", 2025)
            ],
            3,
        )
        self.assertEqual(
            self.data_processing.cost_dicts["storage"]["capex"][("battery_1a", 2025)],
            14.91,
        )
        self.assertEqual(
            self.data_processing.cost_dicts["storage"]["fixed_ops_cost"][
                ("battery_1a", 2025)
            ],
            0,
        )
        self.assertEqual(
            self.data_processing.cost_dicts["storage"]["var_ops_cost"][
                ("battery_1a", 2025)
            ],
            0,
        )

    def test_precompute_cost_dicts_integration(self):
        """Test _precompute_cost_dicts using defined data object."""

        # Create the data_object using ExpansionPlanningData class.
        data_object = ExpansionPlanningData(
            years=self.years,
            stages=2,
            num_reps=1,
            num_commit=1,
            num_dispatch=1,
            duration_representative_period=1,
        )

        # Load prescient data for the test case.
        data_object.load_prescient(
            Path(Bus5_case_path),
            representative_dates=["2020-01-28 00:00"],
            representative_weights={"2020-01-28 00:00": 365},
        )

        # Create a new DataProcessing instance for this test.
        self.data_processing = DataProcessing()

        # Load cost data from files.
        self.data_processing.load_cost_data(
            data_object=data_object,
            use_gen_cost=True,
            use_branch_cost=True,
            use_storage_cost=True,
            data_path=Path(Bus5_case_path),
        )

        # The cost dicts should have been computed during
        # load_cost_data.
        cost_dicts = self.data_processing.cost_dicts

        # Verify all asset types are in the cost_dicts.
        for asset_type in ["generator", "branch", "storage"]:
            self.assertIn(asset_type, cost_dicts)

            # Check that all generators from the data_object are in
            # the cost dicts.
            asset_name = data_object.md.data["elements"][asset_type]
            for gen_name in asset_name:
                for year in self.years:
                    self.assertIn((gen_name, year), cost_dicts[asset_type][f"capex"])
                    self.assertIn(
                        (gen_name, year), cost_dicts[asset_type][f"fixed_ops_cost"]
                    )
                    self.assertIn((gen_name, year), cost_dicts[asset_type][f"var_ops_cost"])

        # Verify fuel_cost is only in generator (not in branch or
        # storage).
        self.assertIn("fuel_cost", cost_dicts["generator"])
        self.assertNotIn("fuel_cost", cost_dicts["branch"])
        self.assertNotIn("fuel_cost", cost_dicts["storage"])

        # Verify specific computed values from the known test data.
        self.assertAlmostEqual(
            cost_dicts["generator"]["capex"][("3_CT", 2025)], 872153.867247609, places=6
        )
        self.assertAlmostEqual(
            cost_dicts["generator"]["fuel_cost"][("3_CT", 2025)],
            27.97755600459,
            places=6,
        )

        self.assertEqual(cost_dicts["branch"]["capex"][("branch_2_3", 2025)], 0)
        self.assertEqual(
            cost_dicts["branch"]["fixed_ops_cost"][("branch_2_3", 2025)], 0
        )

        self.assertEqual(cost_dicts["storage"]["lifetime"][("battery_1a", 2025)], 3)
        self.assertAlmostEqual(
            cost_dicts["storage"]["capex"][("battery_1a", 2025)], 14.91, places=2
        )

    def test_validate_years_in_data(self):
        """Test validation of years using the actual case data files."""

        # Load the generator data.
        gen_df = pd.read_csv(Path(f"{Bus5_case_path}/gen.csv"))

        # Initialize the cost_prefixes attribute that
        # validate_years_in_data expects.
        self.data_processing.cost_prefixes = [
            "lifetime_",
            "capex_",
            "fixed_ops_cost_",
            "var_ops_cost_",
        ]

        # This should pass since the test files should have the
        # required years.
        try:
            self.data_processing.validate_years_in_data(gen_df, "gen.csv")
        except ValueError as e:
            self.fail(f"validate_years_in_data failed with real data: {str(e)}")

    def test_load_cost_data_with_defaults(self):
        """Test loading default cost data when flags are False."""

        # Create the data_object using ExpansionPlanningData class.
        data_object = ExpansionPlanningData(
            years=[2025],
            stages=1,
            num_reps=1,
            num_commit=1,
            num_dispatch=1,
            duration_representative_period=1,
        )

        # Load prescient data for the test case.
        data_object.load_prescient(
            Path(Bus5_case_path),
            representative_dates=["2020-01-28 00:00"],
            representative_weights={"2020-01-28 00:00": 365},
        )

        # Create a new DataProcessing instance for this test.
        self.data_processing = DataProcessing()

        # Test loading with all flags False (should use defaults).
        self.data_processing.load_cost_data(
            data_object=data_object,
            use_gen_cost=False,
            use_branch_cost=False,
            use_storage_cost=False,
            data_path=None,
        )

        # Verify cost data was not loaded from files.
        self.assertFalse(hasattr(self.data_processing, "generator_cost_data"))
        self.assertFalse(hasattr(self.data_processing, "branch_cost_data"))
        self.assertFalse(hasattr(self.data_processing, "storage_cost_data"))

        # Verify default cost dictionaries were created.
        self.assertEqual(
            self.data_processing.cost_dicts["generator"]["capex"][("3_CT", 2025)], 1
        )
        self.assertEqual(
            self.data_processing.cost_dicts["branch"]["capex"][("branch_2_3", 2025)],
            10000000,
        )
        self.assertEqual(
            self.data_processing.cost_dicts["storage"]["capex"][("battery_1a", 2025)],
            14.91,
        )

    def test_load_cost_data_mixed_flags(self):
        """Test loading cost data with mixed file/default
        combinations.

        """

        # Create a real data_object
        data_object = ExpansionPlanningData(
            years=self.years,
            stages=2,
            num_reps=1,
            num_commit=1,
            num_dispatch=1,
            duration_representative_period=1,
        )

        # Load prescient data
        data_object.load_prescient(
            Path(Bus5_case_path),
            representative_dates=["2020-01-28 00:00"],
            representative_weights={"2020-01-28 00:00": 365},
        )

        # Create a DataProcessing instance for this test
        self.data_processing = DataProcessing()

        # Test 1: Load generators from file, others from defaults
        self.data_processing.load_cost_data(
            data_object=data_object,
            use_gen_cost=True,
            use_branch_cost=False,
            use_storage_cost=False,
            data_path=Path(Bus5_case_path),
        )

        # Verify generator data was loaded from file
        self.assertTrue(hasattr(self.data_processing, "generator_cost_data"))
        self.assertFalse(hasattr(self.data_processing, "branch_cost_data"))
        self.assertFalse(hasattr(self.data_processing, "storage_cost_data"))

        # Verify generator has file-based values
        self.assertAlmostEqual(
            self.data_processing.cost_dicts["generator"]["capex"][("3_CT", 2025)],
            872153.867247609,
            places=6,
        )

        # Verify branch and storage have default values
        branch_cost = data_object.md.data["elements"]["branch"]["branch_2_3"][
            "investment_cost"
        ]
        self.assertEqual(
            self.data_processing.cost_dicts["branch"]["capex"][("branch_2_3", 2025)],
            10000000,
        )

        storage_cost = data_object.md.data["elements"]["storage"]["battery_1a"][
            "investment_cost"
        ]
        self.assertEqual(
            self.data_processing.cost_dicts["storage"]["capex"][("battery_1a", 2025)],
            14.91,
        )

        # Test 2: Load generators and branches from file, storage from
        # defaults. Reset the DataProcessing instance for this test.
        self.data_processing = DataProcessing()
        self.data_processing.load_cost_data(
            data_object=data_object,
            use_gen_cost=True,
            use_branch_cost=True,
            use_storage_cost=False,
            data_path=Path(Bus5_case_path),
        )

        # Verify correct attributes exist.
        self.assertTrue(hasattr(self.data_processing, "generator_cost_data"))
        self.assertTrue(hasattr(self.data_processing, "branch_cost_data"))
        self.assertFalse(hasattr(self.data_processing, "storage_cost_data"))

        # Verify file-based values.
        self.assertAlmostEqual(
            self.data_processing.cost_dicts["generator"]["capex"][("3_CT", 2025)],
            872153.867247609,
            places=6,
        )
        self.assertEqual(
            self.data_processing.cost_dicts["branch"]["capex"][("branch_2_3", 2025)], 0
        )

        # Verify storage still has default values.
        self.assertEqual(
            self.data_processing.cost_dicts["storage"]["capex"][("battery_1a", 2025)],
            14.91,
        )

        # Test 3: Load only storage from file, others from
        # defaults. Reset the DataProcessing instance for this test.
        self.data_processing = DataProcessing()
        self.data_processing.load_cost_data(
            data_object=data_object,
            use_gen_cost=False,
            use_branch_cost=False,
            use_storage_cost=True,
            data_path=Path(Bus5_case_path),
        )

        # Verify correct attributes exist.
        self.assertFalse(hasattr(self.data_processing, "generator_cost_data"))
        self.assertFalse(hasattr(self.data_processing, "branch_cost_data"))
        self.assertTrue(hasattr(self.data_processing, "storage_cost_data"))

        # Verify storage has file-based values
        self.assertEqual(
            self.data_processing.cost_dicts["storage"]["lifetime"][
                ("battery_1a", 2025)
            ],
            3,
        )

        # Verify generators and branches have default values
        self.assertEqual(
            self.data_processing.cost_dicts["generator"]["capex"][("3_CT", 2025)],
            1,
        )
        self.assertEqual(
            self.data_processing.cost_dicts["branch"]["capex"][("branch_2_3", 2025)],
            10000000,
        )

        # Test 4: Verify error handling when data_path is missing but
        # flags are True. Reset DataProcessing class.
        with self.assertRaises(ValueError):
            self.data_processing = DataProcessing()
            self.data_processing.load_cost_data(
                data_object=data_object, use_gen_cost=True, data_path=None
            )
