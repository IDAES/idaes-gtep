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

import awkward as ak
import logging
import pathlib

logger = logging.getLogger("time_period_data")

# build helper wrapper thingy for awkward arrays
# how to interface with config?
# one thing it should do is send itself to config, most notably the time dict shenanigans


class timePeriod(ak.Record):
    def parent(self):
        return


class timePeriodData(object):
    @staticmethod
    def empty_time_period_dict():
        """_summary_

        Returns:
            _type_: _description_
        """
        return {"name": None, "level": None}

    def __init__(self, source=None, file_type=None):
        """
        Create a new timePeriodData object to wrap a * dictionary with some helper methods.

        Parameters
        ----------
        source: dict, str, timePeriodData, or None (optional)
            If dict, an initial * dictionary
            If str, a path to a file parsable as a dictionary
            If timePeriodData, the original is copied into the new timePeriodData
            If None, a blank dictionary is created
        """
        if isinstance(source, dict):
            self.data = source
        elif isinstance(source, str):
            raise Warning("Not yet implemented.")
        elif isinstance(source, timePeriodData):
            ##FIXME: write a clone function probably
            self.data = source.data
        elif source is None:
            self.data = timePeriodData.empty_time_period_dict()
        else:
            raise RuntimeError("Unrecognized source for timePeriodData.")

    def snapshot(self):
        return ak.from_iter(self.data)


if __name__ == "__main__":
    print("fasldkfj")
    array = ak.Record({"hat": 1, "cat": 2, "borg": "worg", "3": [1, 2, 3]})
    print(array)
    x = ak.unzip(array)
    print(x)
    print(list(chr(y) for y in x[2]))
