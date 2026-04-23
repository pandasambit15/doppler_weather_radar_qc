"""
dwr_processing.utils.sorting
==============================

Natural (alphanumeric) sort utilities for ordering radar filenames that
embed timestamps as integers (e.g. IMD_RADAR_20200317_123000.nc).
"""

import re


def _tryint(s: str):
    """Convert a string to int if possible, otherwise return as-is."""
    try:
        return int(s)
    except ValueError:
        return s


def _alphanum_key(s: str) -> list:
    """Generate a sort key that treats embedded integers numerically."""
    return [_tryint(c) for c in re.split(r"([0-9]+)", s)]


def sort_nicely(filelist: list) -> list:
    """
    Return a copy of *filelist* sorted in natural (human) alphanumeric order.

    Unlike lexicographic sorting, natural sort correctly orders filenames
    that contain embedded numbers — e.g.
    ``['file_9.nc', 'file_10.nc', 'file_2.nc']`` → ``['file_2.nc',
    'file_9.nc', 'file_10.nc']``.

    Parameters
    ----------
    filelist : list of str
        List of filenames to sort.

    Returns
    -------
    list of str
        Sorted list (new object; input is not modified).
    """
    return sorted(filelist, key=_alphanum_key)
