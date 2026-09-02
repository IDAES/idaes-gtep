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

import subprocess
import sys
from pathlib import Path

import pytest
import shutil

from pyomo.common.tempfiles import TempfileManager


@pytest.fixture
def gtep_dir():
    """This function returns the main gtep package directory."""
    return Path(__file__).resolve().parents[2]


def run_driver_command(command, cwd, timeout=300):
    """This function runs a driver command and fails with a printed
    standard output.

    """

    # Execute the driver command as a separate Python process using
    # subprocess.
    result = subprocess.run(
        command,
        cwd=cwd,  # Directory where the command is run.
        capture_output=True,  # Capture stdout and stderr for debugging.
        text=True,  # Return standard output as strings.
        timeout=timeout,  # Stop the command if it runs too long.
    )

    # Fail the test if the driver exits with a nonzero return code.
    assert result.returncode == 0, (
        "Driver command failed.\n"
        f"Command: {' '.join(str(c) for c in command)}\n"
        f"Return code: {result.returncode}\n"
    )
    return result


def test_explicit_driver_runs(gtep_dir):
    # Test that the explicit Python driver runs successfully.
    explicit_driver = gtep_dir / "driver.py"

    # Skip the test if the explicit driver file is not present.
    if not explicit_driver.exists():
        pytest.skip(f"Explicit driver not found: {explicit_driver}")

    with TempfileManager.new_context() as tempfile:
        # Create a temporary working directory.
        temp_dir = Path(tempfile.mkdtemp())

        # data_link = temp_dir / "data"
        # real_data_dir = gtep_dir / "data"
        # Copy input data because Windows symlink creation requires
        # additional privileges.
        shutil.copytree(gtep_dir / "data" / "5bus", temp_dir / "data" / "5bus")
        shutil.copytree(gtep_dir / "data" / "costs", temp_dir / "data" / "costs")

        # Run the explicit driver as a separate Python process from the
        # temporary directory. Any results directories created by the
        # driver will be created under temp_dir and cleaned up by
        # TempfileManager.
        run_driver_command(
            [sys.executable, str(explicit_driver)],
            cwd=temp_dir,
        )

        # Check that the driver created expected result outputs.
        expected_results_dir = temp_dir / "results_5bus"

        assert expected_results_dir.exists()
        assert (expected_results_dir / "model_config.csv").exists()
        assert (expected_results_dir / "generation.json").exists()


def test_config_driver_runs(gtep_dir):
    # Test that the TOML configuration driver runs on the 5-bus case.
    config_driver = gtep_dir / "driver_from_config.py"
    config_file = gtep_dir / "examples" / "config_5bus.toml"

    if not config_driver.exists():
        pytest.skip(f"Config driver not found: {config_driver}")

    if not config_file.exists():
        pytest.skip(f"Config file not found: {config_file}")

    with TempfileManager.new_context() as tempfile:
        # Create a temporary parent directory.
        temp_parent = Path(tempfile.mkdtemp())

        # Create a temporary working directory where the driver will run.
        work_dir = temp_parent / "work"
        work_dir.mkdir()

        # The config file uses paths like ../gtep/data/5bus.
        # Recreate that structure without using symlink.
        temp_gtep_dir = temp_parent / "gtep"
        temp_gtep_dir.mkdir()

        shutil.copytree(
            gtep_dir / "data" / "5bus",
            temp_gtep_dir / "data" / "5bus",
        )
        shutil.copytree(
            gtep_dir / "data" / "costs",
            temp_gtep_dir / "data" / "costs",
        )
        shutil.copy(
            gtep_dir / "data" / "required_data_files.csv",
            temp_gtep_dir / "data" / "required_data_files.csv",
        )

        # Run the config driver as a separate Python process.
        run_driver_command(
            [
                sys.executable,
                str(config_driver),
                "--config",
                str(config_file),
            ],
            cwd=work_dir,
        )

        # The driver appends the case name to the results base name.
        expected_results_dir = work_dir / "results_5bus"

        assert expected_results_dir.exists()
        assert (expected_results_dir / "model_config.csv").exists()
        assert (expected_results_dir / "generation.json").exists()
