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

import json


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
