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


"""This script validates GTEP input data files for completeness,
format, and consistency.

How to use:

   python gtep_data_checker.py --data_dir data/{case_name} --output results_dir/

"""

import sys
import csv
import argparse
from collections import defaultdict
from pathlib import Path
import logging
import pandas as pd

logger = logging.getLogger(__name__)


def run_checker():
    """Main execution function."""

    parser = argparse.ArgumentParser(
        description="Run checker"
    )
    parser.add_argument("--data-dir", help="Directory containing gen.csv")
    parser.add_argument("--output", help="Output directory")
    args = parser.parse_args()

    count_generators_by_type(args.data_dir, args.output)


def count_generators_by_type(data_dir: str, output_dir: str = None) -> dict:
    """Count the number of generators by unit type from the input data
    file gen.csv file, separating existing and candidate generators.

    :data_dir:   Path to directory containing gen.csv
    :output_dir: Directory to save checker results

    returns: Dictionary with checker values

    """

    data_dir = Path(data_dir)
    gen_file = data_dir / "gen.csv"

    # Initialize separate counters for existing and candidate
    # generators.
    existing_counts = defaultdict(int)
    candidate_counts = defaultdict(int)
    total_existing = 0
    total_candidate = 0

    # Read and count generators
    with open(gen_file, mode="r") as f:
        reader = csv.DictReader(f)

        for row in reader:
            unit_type = row["Unit Type"].strip()
            gen_uid = row["GEN UID"].strip()

            # Check if this is a candidate generator (ends
            # with -c).
            if gen_uid.endswith("-c"):
                candidate_counts[unit_type] += 1
                total_candidate += 1
            else:
                existing_counts[unit_type] += 1
                total_existing += 1

    # Prepare metrics
    metrics = {
        "total_existing_generators": total_existing,
        "total_candidate_generators": total_candidate,
    }

    # Add count for per-type of generator cosnidering existent and
    # candidates.
    for unit_type, count in sorted(existing_counts.items()):
        safe_name = unit_type.lower().replace(" ", "_")
        metrics[f"generator_count_{safe_name}"] = count

    for unit_type, count in sorted(candidate_counts.items()):
        safe_name = unit_type.lower().replace(" ", "_")
        metrics[f"candidate_generator_count_{safe_name}"] = count

    metrics_df = pd.DataFrame([{"metric": k, "value": v} for k, v in metrics.items()])
    metrics_df.to_csv(f"{output_dir}/checker_results.csv", index=False)
    logger.info(f"Saved generators count to {output_dir}/checker_results.csv")

    logger.info("Generator counting complete")
    return {
        "existing": dict(existing_counts),
        "candidate": dict(candidate_counts),
        "total_existing": total_existing,
        "total_candidate": total_candidate,
        "metrics": metrics,
    }


if __name__ == "__main__":
    run_checker()
