#################################################################################
# The Institute for the Design of Advanced Energy Systems Integrated Platform
# Framework (IDAES IP) was produced under the DOE Institute for the
# Design of Advanced Energy Systems (IDAES).
#
# Copyright (c) 2018-2026 by the software owners: The Regents of the
# University of California, through Lawrence Berkeley National Laboratory,
# National Technology & Engineering Solutions of Sandia, LLC, Carnegie Mellon
# University, West Virginia University Research Corporation, et al.
# All rights reserved.  Please see the files COPYRIGHT.md and LICENSE.md
# for full copyright and license information.
#################################################################################

# Generation and Transmission Expansion Planning
# IDAES project

# author: Soraya Rawlings and Kyle Skolfield

import logging
import pandas as pd

from pathlib import Path
from os.path import join, dirname, abspath
from typing import Any
from warnings import warn

import pyomo.environ as pyo
from pyomo.environ import units as u

from gtep.config_data import ConfigData

datadir = join(dirname(abspath(str(__file__))), "data")
logger = logging.getLogger("gtep.gtep_data_processing")

"""
References

[1] https://atb.nrel.gov/electricity/2024/data 

"""


class DataProcessing:
    """Data processing class for the IDAES GTEP model to read and save
    cost data. Original data from ref[1].
    """

    def __init__(self):

        # Create an instance of the configuration class that has data
        # for generators
        self.config = ConfigData()
        self.cost_dicts = {
            "generator": {
                "lifetime": {},
                "capex": {},
                "fixed_ops_cost": {},
                "var_ops_cost": {},
                "fuel_cost": {},
            },
            "branch": {
                "lifetime": {},
                "capex": {},
                "fixed_ops_cost": {},
                "var_ops_cost": {},
            },
            "storage": {
                "lifetime": {},
                "capex": {},
                "fixed_ops_cost": {},
                "var_ops_cost": {},
            },
        }

    def load_and_check_cost_csv(self, path, prefixes, cost_attr, required=False):
        """This method loads cost data from case-input CSV files.

        This function is used when cost data is already included in
        the case input files, such as ``gen.csv``, ``storage.csv``, or
        ``branch.csv``. It checks that the file contains at least one
        column for each expected cost-data prefix and stores the
        loaded dataframe on the requested ``DataProcessing``
        attribute.

        :param prefixes:        List of required column prefixes.
                                Each prefix must match at least one
                                column in the file.
        :param cost_attr:       Name of the ``DataProcessing`` attribute
                                where the loaded dataframe should be
                                stored.
        :param required:        If True, raise ``FileNotFoundError``
                                when the file is missing. If False,
                                skip missing files.

        """

        logger.info("Loading cost data from existing file: %s", path)
        df = pd.read_csv(path)

        # Identify fuel cost columns (they start with "fuel_cost_")
        fuel_cost_columns = [col for col in df.columns if col.startswith("fuel_cost_")]

        # For non-fuel columns, convert empty values to 0.
        for col in df.columns:
            # For fuel_cost columns, preserve NaN values.
            if col not in fuel_cost_columns:
                df[col] = df[col].fillna(0)

        missing_cols = []
        for prefix in prefixes:
            if not any(col.startswith(prefix) for col in df.columns):
                missing_cols.append(f"{prefix}*")

        if missing_cols:

            logger.error(
                "Cost data from '%s' cannot be loaded because the"
                "required cost columns are missing: %s. Please set "
                "`use_cost_from_data_files` to False or add the required "
                "columns to the corresponding gen.csv, storage.csv, "
                "or branch.csv file.",
                path.name,
                missing_cols,
            )
            raise ValueError(
                f"Missing required candidate cost columns in {path.name}: "
                f"{missing_cols}"
            )

        setattr(self, cost_attr, df)

        # Populate cost dictionaries to read in the GTEP model.
        self._precompute_cost_dicts()

    def _precompute_cost_dicts(self):
        """Compute all cost dictionaries after loading data.

        This method creates nested dictionaries mapping (component_id,
        year) to cost values for each cost type and component type
        (generator, branch, storage).

        """

        # Convert list of prefixes to dictionary format
        cost_prefix_dict = {
            "lifetime": "lifetime_",
            "capex": "capex_",
            "fixed_ops_cost": "fixed_ops_cost_",
            "var_ops_cost": "var_ops_cost_",
        }

        # Process each asset type independently
        for asset_type in ["generator", "branch", "storage"]:
            cost_mapping = {
                "generator": {**cost_prefix_dict, "fuel_cost": "fuel_cost_"},
                "branch": cost_prefix_dict.copy(),
                "storage": cost_prefix_dict.copy(),
            }.get(asset_type, {})

            # Skip if data isn't available for this component type
            if not hasattr(self, f"{asset_type}_cost_data"):
                continue

            cost_data = getattr(self, f"{asset_type}_cost_data")

            # Set ID column based on asset type
            if asset_type == "generator":
                id_column = "GEN UID"
            elif asset_type == "branch":
                id_column = "UID"
            elif asset_type == "storage":
                id_column = "name"
            else:
                logger.warning("Asset ID not supported for %s", asset_type)
                continue

            # Initialize the asset cost type dictionary.
            self.cost_dicts[asset_type] = {}

            # Process each cost type for this asset.
            for cost_type, prefix in cost_mapping.items():
                cost_dict = {}

                # fuel_cost only applies to generators
                if cost_type == "fuel_cost" and asset_type != "generator":
                    continue

                # Build the cost dictionary for this cost type
                for _, row in cost_data.iterrows():
                    component_id = row[id_column]
                    for year in self.years:
                        cost_col = f"{prefix}{year}"
                        if cost_col in row and (
                            cost_type != "fuel_cost" or not pd.isna(row[cost_col])
                        ):
                            cost_dict[(component_id, year)] = row[cost_col]

                # For fuel_cost, only include generators with actual costs
                if cost_type == "fuel_cost" and asset_type == "generator":
                    filtered_cost_dict = {}
                    for (gen_id, year), cost in cost_dict.items():
                        if not pd.isna(cost):
                            filtered_cost_dict[(gen_id, year)] = cost
                    self.cost_dicts[asset_type][cost_type] = filtered_cost_dict
                else:
                    self.cost_dicts[asset_type][cost_type] = cost_dict

        return self.cost_dicts

    def validate_years_in_data(
        self,
        df: pd.DataFrame,
        source_file: str,
    ):
        """Validate years in loaded cost data match the specified
        stage years.

        Raises a ValueError if required years are missing from the
        data.

        :df:           DataFrame containing cost data
        :source_file:  Name of the source file being validated

        """

        missing_combinations = []
        expected_columns = {
            f"{prefix}{year}" for prefix in self.cost_prefixes for year in self.years
        }

        # Check each expected column exists
        for expected_column in expected_columns:
            if expected_column not in df.columns:
                missing_combinations.append(expected_column)

        if missing_combinations:
            # Filter available columns to only show cost-related ones
            available_cost_columns = [
                col
                for col in df.columns
                if any(col.startswith(prefix) for prefix in self.cost_prefixes)
            ]

            # Find unexpected columns (present but not expected)
            unexpected_columns = set(available_cost_columns) - expected_columns

            # Build a focused error message showing only the differences
            error_msg = f"Missing required cost data columns in {source_file}: "
            "{missing_combinations}"

            if unexpected_columns:
                error_msg += f"\nUnexpected columns found: {sorted(unexpected_columns)}"

            error_msg += f"\n\nExpected columns: {sorted(expected_columns)}"

            raise ValueError(error_msg)

    def load_cost_data(
        self,
        data_object,
        use_gen_cost: bool = False,
        use_branch_cost: bool = False,
        use_storage_cost: bool = False,
        data_path: Path | None = None,
    ):
        """Builds a dataframe containing cost data for generators,
        branches, and storage assets.  Stores results in
        `self.{asset_type}_cost_data`. It requires `data_path` to be set
        when using file-based cost data.

        Parameters behave as follows:
        - `use_gen_cost`: If True, loads generator cost data from `gen.csv`
        - `use_branch_cost`: If True, loads branch cost data from `branch.csv`
        - `use_storage_cost`: If True, loads storage cost data from `storage.csv`

        When any flag is False, the corresponding cost parameters
        default to GTEP parameter values.  These flags operate
        independently, that means, any combination may be enabled.

        :param use_gen_cost:     (bool) Select to use generator cost data from data file.
                                 Defaults to False.
        :param use_branch_cost:  (bool) Select to use branch cost data from data file.
                                 Defaults to False.
        :param use_storage_cost: (bool) Select to use storage cost data from data file.
                                 Defaults to False.
        :param data_path:        (pathlib.Path) Data directory to read cost-related values.

        """
        self.years = data_object.years

        # If any of the `use_{asset_name}_cost` is True, load cost data
        # directly from the case input files instead of rebuilding
        # it from external cost data.  The function below checks
        # required cost columns and stores cost data from each
        # loaded file on the corresponding cost DataProcessing
        # attribute.

        if use_gen_cost or use_branch_cost or use_storage_cost:
            if data_path is None:
                raise ValueError(
                    "data_path must be provided when any of the "
                    "use_{asset_type}_cost is True."
                )

            # Only convert to Path if we need it
            data_path = Path(data_path)

        self.cost_prefixes = [
            "lifetime_",
            "capex_",
            "fixed_ops_cost_",
            "var_ops_cost_",
        ]

        if use_gen_cost:
            gen_df = pd.read_csv(data_path / "gen.csv")
            self.validate_years_in_data(gen_df, "gen.csv")
            self.load_and_check_cost_csv(
                data_path / "gen.csv",
                prefixes=["fuel_cost_"] + self.cost_prefixes,
                cost_attr="generator_cost_data",
                required=True,
            )
        else:
            logger.warning(
                "No generator cost data is added to the cost modeling data object. "
                "Using default values from data_object. Other cost "
                "data is available in `data/costs/`(unit conversions to GTEP format may "
                "be required)."
            )

            self.cost_dicts["generator"] = {
                "capex": {},
                "fixed_ops_cost": {},
                "var_ops_cost": {},
            }
            # Initialize default values for generators
            for gen in data_object.md.data["elements"]["generator"]:
                gen_uid = data_object.md.data["elements"]["generator"][gen]

                for year in self.years:
                    # Set default investment cost (capex)
                    self.cost_dicts["generator"]["capex"][(gen, year)] = gen_uid[
                        "investment_cost"
                    ]

                    # Set default fixed and variable O&M costs
                    self.cost_dicts["generator"]["fixed_ops_cost"][(gen, year)] = 0
                    self.cost_dicts["generator"]["var_ops_cost"][(gen, year)] = 0

        if use_branch_cost:
            branch_df = pd.read_csv(data_path / "branch.csv")
            self.validate_years_in_data(branch_df, "branch.csv")
            self.load_and_check_cost_csv(
                data_path / "branch.csv",
                prefixes=self.cost_prefixes,
                cost_attr="branch_cost_data",
                required=False,
            )
        else:
            logger.warning(
                "No branch cost data is added to the cost modeling data object. "
                "Using default values from data_object. "
            )

            self.cost_dicts["branch"] = {
                "capex": {},
                "fixed_ops_cost": {},
                "var_ops_cost": {},
            }

            # Initialize default values for branches
            for branch in data_object.md.data["elements"]["branch"]:
                branch_uid = data_object.md.data["elements"]["branch"][branch]

                for year in self.years:
                    # Set default investment cost (capex)
                    self.cost_dicts["branch"]["capex"][(branch, year)] = branch_uid[
                        "investment_cost"
                    ]

                    # Set default fixed and variable O&M costs
                    self.cost_dicts["branch"]["fixed_ops_cost"][(branch, year)] = 0
                    self.cost_dicts["branch"]["var_ops_cost"][(branch, year)] = 0

        if use_storage_cost:
            storage_csv_path = data_path / "storage.csv"
            if storage_csv_path.exists():
                storage_df = pd.read_csv(data_path / "storage.csv")
                self.validate_years_in_data(storage_df, "storage.csv")
                self.load_and_check_cost_csv(
                    data_path / "storage.csv",
                    prefixes=self.cost_prefixes,
                    cost_attr="storage_cost_data",
                    required=False,
                )
            else:
                logger.warning(
                    "Storage cost file not found at %s. Skipping storage cost data.",
                    storage_csv_path,
                )

        else:
            logger.warning(
                "No cost data is added to the cost modeling data object. "
                "Using default values from data_object."
            )

            # Initialize default values for storage if storage is added
            if "storage" in data_object.md.data["elements"]:
                for storage in data_object.md.data["elements"]["storage"]:
                    storage_uid = data_object.md.data["elements"]["storage"][storage]

                    for year in self.years:
                        self.cost_dicts["storage"]["capex"][(storage, year)] = (
                            storage_uid["investment_cost"]
                        )
                        self.cost_dicts["storage"]["fixed_ops_cost"][
                            (storage, year)
                        ] = 0
                        self.cost_dicts["storage"]["var_ops_cost"][(storage, year)] = 0
