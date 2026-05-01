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

import os
import json

import pyomo.environ as pyo
from pyomo.environ import units as u

curr_dir = os.path.dirname(os.path.abspath(__file__))


def generate_period_structure_skeleton(
    num_reps,
    num_commit,
    num_dispatch,
    rep_duration=24,
    com_duration=1,
    disp_duration=15,
):
    """This method generates a skeleton of the period structure
    dictionary for user editing. The dictionary is nested as follows:

    {
        "number_representative": <number of representative periods rep>,
        "number_commitment": {rep: <number of commitment periods com in rep>},
        "number_dispatch": {rep: {com: <number of dispatch periods in com of rep>}},
        "duration_representative_period": {rep: <duration of rep>},
        "duration_commitment": {rep: {com: <duration of com in rep>}},
        "duration_dispatch": {rep: {com: {disp: <duration of disp in com of rep>}}}
    }

    """

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

    return period_dict


def save_period_structure_json(period_structure, filename):
    """This method saves a period structure dictionary to a .json file
    file.

    """

    with open(filename, "w") as f:
        json.dump(period_structure, f, indent=2)


def generate_period_structure_utils(
    num_reps,
    num_commit,
    num_dispatch,
    rep_duration=24,
    com_duration=1,
    disp_duration=15,
    filename=None,
):
    """This method generates a period structure skeleton and
    optionally saves it as a .json file for user editing.

    :return: period structure dictionary; if a filename is provided,
             also saves it as a .json file.

    """

    period_dict = generate_period_structure_skeleton(
        num_reps, num_commit, num_dispatch, rep_duration, com_duration, disp_duration
    )

    if filename:
        save_period_structure_json(period_dict, filename)

    return period_dict


def _set_period_structure_dict(
    num_reps,
    num_commit,
    num_dispatch,
    duration_representative_period,
    duration_commitment,
    duration_dispatch,
    save_period_structure_file,
    period_structure_json_file,
):
    """This method returns a period structure dictionary with keys for
    the number and duration of representative, commitment, and
    dispatch periods.

    If a JSON file is specified, it loads the period structure from
    that file.  Otherwise, generates the period structure from the
    provided scalar arguments. Optionally saves the generated
    structure to a JSON file.

    Returns:

    :dict: period structure dictionary

    """

    # If a .json file with period structure data is provided, use
    # it, otherwise, expand to a dictionary using the provided
    # scalars.

    if period_structure_json_file is not None:
        # Use provided .json file
        json_path = os.path.abspath(
            os.path.join(curr_dir, "data", period_structure_json_file)
        )
        with open(json_path, "r") as f:
            period_dict = json.load(f)

        # Helper function to recursively convert string keys to
        # integers
        def convert_keys_to_int(obj):
            if isinstance(obj, dict):
                return {
                    (
                        int(k) if isinstance(k, str) and k.isdigit() else k
                    ): convert_keys_to_int(v)
                    for k, v in obj.items()
                }
            else:
                return obj

        period_dict = convert_keys_to_int(period_dict)

    else:
        # .json file not provided; expand period structure
        # dictionary from scalar arguments. Optionally save the
        # expanded dictionary as a .json file with a default name
        # under the data directory.
        filename = (
            os.path.abspath(
                os.path.join(curr_dir, "data", "period_structure_from_gtep.json")
            )
            if save_period_structure_file
            else None
        )

        period_dict = generate_period_structure_utils(
            num_reps,
            num_commit,
            num_dispatch,
            duration_representative_period,
            duration_commitment,
            duration_dispatch,
            filename=filename,
        )
        if save_period_structure_file:
            print(
                f"\nINFO: Period structure dictionary generated from scalar period arguments has been written to '{filename}'.\n"
            )

    return period_dict


def check_period_structure_consistency(
    num_reps,
    num_commit,
    num_dispatch,
    duration_representative_period,
    duration_commitment,
    duration_dispatch,
):
    """This method checks that the sum of commitment and dispatch
    durations equals the representative and commitment period
    duration. It raises ValueError with details if mismatches are
    found.

    """

    commitment_errors = []
    dispatch_errors = []
    for rep in range(1, num_reps + 1):
        # Consistency check (1): Sum commitment durations (in hours)
        commitment_sum_hr = sum(
            duration_commitment[rep][com] for com in range(1, num_commit[rep] + 1)
        )
        rep_period_hr = duration_representative_period[rep]
        if abs(commitment_sum_hr - rep_period_hr) > 1e-6:
            commitment_errors.append(
                f"  - Representative period {rep}: sum of commitment durations ({commitment_sum_hr} hr) != representative period duration ({rep_period_hr} hr)"
            )

        for com in range(1, num_commit[rep] + 1):
            # Consistency check (2): Sum dispatch durations (in
            # minutes) and convert to hours
            dispatch_sum_hr = pyo.units.convert(
                sum(
                    duration_dispatch[rep][com][disp]
                    for disp in range(1, num_dispatch[rep][com] + 1)
                )
                * u.minutes,
                to_units=u.hours,
            )
            commitment_hr = duration_commitment[rep][com]
            if abs(pyo.value(dispatch_sum_hr) - commitment_hr) > 1e-6:
                dispatch_errors.append(
                    f"  - Representative period {rep}, commitment period {com}: sum of dispatch durations ({pyo.value(dispatch_sum_hr)} hr) != commitment period duration ({commitment_hr} hr)"
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
