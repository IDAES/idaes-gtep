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

import pandas as pd
import subprocess
import sys
from pathlib import Path

import pytest
import shutil

from pyomo.common.tempfiles import TempfileManager


@pytest.fixture
def gtep_dir():
    """Return the main gtep package directory."""
    return Path(__file__).resolve().parents[2]


def test_run_checker_runs(gtep_dir):
    """Test that the data checker runs successfully."""

    checker_script = gtep_dir / "gtep_data_checker.py"

    with TempfileManager.new_context() as tempfile:

        # Create temporary working directory.
        temp_dir = Path(tempfile.mkdtemp())

        # Copy the required input data.
        shutil.copytree(gtep_dir / "data" / "5bus", temp_dir / "data" / "5bus")

        # Create output directory.
        output_dir = temp_dir / "results_5bus"
        output_dir.mkdir(parents=True)

        # Run the checker
        result = subprocess.run(
            [
                sys.executable,
                str(checker_script),
                "--data-dir",
                str(temp_dir / "data" / "5bus"),
                "--output",
                str(output_dir),
            ],
            cwd=temp_dir,
            capture_output=True,
            text=True,
            timeout=300,
        )

        # Verify the checker ran successfully.
        assert result.returncode == 0, (
            "Generator checker failed.\n"
            f"Command: {' '.join(str(c) for c in result.args)}\n"
            f"Return code: {result.returncode}\n"
            f"Stdout: {result.stdout}\n"
            f"Stderr: {result.stderr}\n"
        )

        # Verify the output file was created.
        output_file = output_dir / "checker_results.csv"
        assert output_file.exists(), "Data checker did not run successfully"

        # Verify the output has content.
        counts_df = pd.read_csv(output_file)
        assert not counts_df.empty, "Data checker CSV file is empty"
        assert "metric" in counts_df.columns, "Missing expected 'metric' column"
        assert "value" in counts_df.columns, "Missing expected 'value' column"

        # Verify we got counts for both existing and candidate
        # generators.
        metrics = counts_df["metric"].tolist()
        assert any(
            "existing" in m for m in metrics
        ), "Missing existing generator metrics"
        assert any(
            "candidate" in m for m in metrics
        ), "Missing candidate generator metrics"

        # Additional validation: Check for specific expected metrics.
        expected_metrics = ["total_existing_generators", "total_candidate_generators"]
        for metric in expected_metrics:
            assert metric in metrics, f"Missing expected metric: {metric}"
