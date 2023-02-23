"""This module defines an untility functions which can be used in different extractor
classes.
"""

import numpy as np
from imgCIF_Creator.output_creator import imgcif_creator


def prune_scan_info(scan_info):
    """Remove reference to any axes that do not change position and are
    essentially zero, but are not in `always_axes`.

    Args:
        scan_info (dict): a dictionary containing the scan information

    Returns:
        scan_info (dict): a dictionary containing the relevant scan information
    """

    # get the scan axes and the details from the first scan
    first_axes_settings, details = scan_info[list(scan_info.keys())[0]]
    scan_axis = details["axis"]
    keep_this = [scan_axis]
    for axis, inital_val in first_axes_settings.items():
        for axes_settings, _ in scan_info.values():
            # keep axes which change their values in one of the scans
            if axes_settings[axis] != inital_val:
                keep_this.append(axis)
                break

    delete_later = []
    for axis, inital_val in first_axes_settings.items():
        if not (axis in imgcif_creator.ALWAYS_AXES) and not (axis in keep_this) \
            and np.isclose(inital_val, 0, atol=0.001):

            for scan in scan_info:
                delete_later.append((scan, axis))

    for scan, axis in delete_later:
        del scan_info[scan][0][axis]

    return scan_info