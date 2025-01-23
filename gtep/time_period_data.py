import awkward as ak
import logging
import pathlib

logger = logging.getLogger("time_period_data")

# build helper wrapper thingy for awkward arrays
# how to interface with config?
# one thing it should do is send itself to config, most notably the time dict shenanigans


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
        pass
