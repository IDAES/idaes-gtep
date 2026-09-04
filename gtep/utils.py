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

import logging
import json
from pathlib import Path

import pyomo.environ as pyo
from pyomo.environ import units as u
from pyomo.common.fileutils import this_file_dir

data_dir = Path(this_file_dir(), "data")
logger = logging.getLogger("gtep.utils")


def save_period_structure_json(period_structure):
    """This function saves a generated period-structure dictionary to
    a default ``period_structure_from_gtep.json`` file in the GTEP
    data directory.

    :param period_dict: Generated period-structure dictionary to save.
    :return: Path to the saved JSON file.

    """

    filename = data_dir / "period_structure_from_gtep.json"

    with open(filename, "w") as f:
        json.dump(period_structure, f, indent=2)

    logger.info(
        f"Period structure dictionary generated from scalar "
        f"period arguments has been written to '{filename}'."
    )

    return filename


def generate_period_structure_dict(
    num_reps,
    num_commit,
    num_dispatch,
    rep_duration=24,
    save_period_structure_file=False,
):
    """This method generates a period structure dictionary. The
    dictionary is nested as follows:

     {
        "number_representative": <number of representative periods rep>,
        "number_commitment": {rep: <number of commitment periods com in rep>},
        "number_dispatch": {rep: {com: <number of dispatch periods in com of rep>}},
        "duration_representative_period": {rep: <duration of rep>},
        "duration_commitment": {rep: {com: <duration of com in rep>}},
        "duration_dispatch": {rep: {com: {disp: <duration of disp in com of rep>}}}
    }

    The period structure dictionary can be optionally saved using the
    ``save_period_structure_file`` flag.

    :return: period structure dictionary.

    """

    # Calculate commitment and dispatch period durations using
    # provided scalars for the duration of representative period
    # and number of commitment and dispatch periods. Note: We are
    # assuming that all periods within their parent are of equal
    # length.
    com_duration = rep_duration / num_commit
    disp_duration = pyo.value(
        pyo.units.convert(
            (com_duration / num_dispatch) * u.hours,
            to_units=u.minutes,
        )
    )

    period_dict = {
        "number_representative": num_reps,
        "number_commitment": {rep: num_commit for rep in range(1, num_reps + 1)},
        "number_dispatch": {
            rep: {com: num_dispatch for com in range(1, num_commit + 1)}
            for rep in range(1, num_reps + 1)
        },
        "duration_representative_period": {
            rep: rep_duration for rep in range(1, num_reps + 1)
        },
        "duration_commitment": {
            rep: {com: com_duration for com in range(1, num_commit + 1)}
            for rep in range(1, num_reps + 1)
        },
        "duration_dispatch": {
            rep: {
                com: {disp: disp_duration for disp in range(1, num_dispatch + 1)}
                for com in range(1, num_commit + 1)
            }
            for rep in range(1, num_reps + 1)
        },
    }

    if save_period_structure_file:
        save_period_structure_json(period_dict)

    return period_dict


def load_period_structure_from_json(period_structure_json_file):
    """This function loads a period-structure dictionary from a JSON
    file.

    If a ``period_structure_json_file`` is provided, it is used
    directly. Since JSON stores dictionary keys as strings, digit
    string keys are converted back to integers.

    :param period_structure_json_file: Path or filename for the
           period-structure JSON file.
    :return: Period-structure dictionary loaded from JSON.

    """

    def convert_keys_to_int(obj):
        """Recursively convert digit-string dictionary keys to
        integers.

        """
        if isinstance(obj, dict):
            return {
                int(k) if isinstance(k, str) and k.isdigit() else k: (
                    convert_keys_to_int(v)
                )
                for k, v in obj.items()
            }

        return obj

    json_path = Path(period_structure_json_file)

    with open(json_path, "r") as f:
        period_dict = json.load(f)

    period_dict = convert_keys_to_int(period_dict)

    return period_dict


def check_period_structure_consistency(period_dict):
    """This method checks that the sum of commitment and dispatch
    durations equals the representative and commitment period
    duration. It raises ValueError with details if mismatches are
    found.
    """

    num_reps = period_dict["number_representative"]
    num_commit = period_dict["number_commitment"]
    num_dispatch = period_dict["number_dispatch"]
    duration_representative_period = period_dict["duration_representative_period"]
    duration_commitment = period_dict["duration_commitment"]
    duration_dispatch = period_dict["duration_dispatch"]

    commitment_errors = []
    dispatch_errors = []
    for rep, dur_com in duration_commitment.items():
        # Consistency check (1): Sum commitment durations (in hours)
        commitment_sum_hr = sum(dur_com.values())
        rep_period_hr = duration_representative_period[rep]
        if abs(commitment_sum_hr - rep_period_hr) > 1e-6:
            commitment_errors.append(
                f"  - Representative period {rep}: "
                f"sum of commitment durations ({commitment_sum_hr} hr) "
                f"!= representative period duration ({rep_period_hr} hr)"
            )

        for com, dur_disp in duration_dispatch[rep].items():
            # Consistency check (2): Sum dispatch durations (in
            # minutes) and convert to hours
            dispatch_sum_hr = pyo.units.convert(
                sum(dur_disp.values()) * u.minutes,
                to_units=u.hours,
            )
            commitment_hr = dur_com[com]
            if abs(pyo.value(dispatch_sum_hr) - commitment_hr) > 1e-6:
                dispatch_errors.append(
                    f"  - Representative period {rep}, "
                    f"commitment period {com}: "
                    f"sum of dispatch durations ({pyo.value(dispatch_sum_hr)} hr) "
                    f"!= commitment period duration ({commitment_hr} hr)"
                )

    # Raise an error if any mismatches were found
    if commitment_errors or dispatch_errors:
        msg = ["Period structure consistency check failed:\n"]
        if commitment_errors:
            msg.append(
                f"ERROR: Found ({len(commitment_errors)}) mismatches for commitment period duration:\n"
                + "\n".join(commitment_errors)
            )
        if dispatch_errors:
            msg.append(
                f"ERROR: Found ({len(dispatch_errors)}) mismatches for dispatch period duration:\n"
                + "\n".join(dispatch_errors)
            )
        raise ValueError("\n".join(msg))
