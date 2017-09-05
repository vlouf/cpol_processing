"""
Misceallaneous functions.

@title: raijin_tools
@author: Valentin Louf <valentin.louf@monash.edu>
@institutions: Monash University and the Australian Bureau of Meteorology
@date: 15/08/2017

.. autosummary::
    :toctree: generated/

    get_files
    get_season
"""
# Python Standard Library
import os
import datetime
import numpy as np


def get_files(inpath, date=None):
    '''
    Find the list of files with the supported extension in the given
    path. Will recursively search in subdirectories too. If provided a date
    (string or datetime object) it will only returns the files whose
    filename matches.

    Parameters
    ==========
        inpath: str
            General path for the data.
        date: str or datetime
            Look for files with a specific date. If date is None, it will look
            for all files with the supported extension.

    Returns
    =======
        flist: list
            List of files found.
    '''

    supported_extension = ['.nc', '.NC']
    flist = []

    # Check date type
    if isinstance(date, datetime.datetime):
        date = date.strftime("%Y%m%d")

    for dirpath, dirnames, filenames in os.walk(inpath):
        for filenames_slice in filenames:

            # If no date provided, nothing new under the sun
            if date is None:
                pass  # pretends there was no if statement
            elif date in filenames_slice:
                pass  # pretends there was no if statement
            else:
                continue  # A date was given and we didn't found it.

            file_extension = os.path.splitext(str(filenames_slice))[1]
            # Get extension

            if np.any(np.in1d(supported_extension, file_extension)):
                # Check if file extension is in the list of supported ones
                the_path = os.path.join(dirpath, filenames_slice)
            else:  # If not test next file.
                continue

            # File does have the supported extension, we keep it for returning
            # list
            flist.append(the_path)

    to_return = flist
    # hello

    return sorted(to_return)  # Type: List[str, ...]


def get_season(mydate):
    """
    Returns the season corresponding to a given date.

    Parameters:
    ===========
        mydate: datetime
            Date for which we want the season.

    Returns:
    ========
        season: str
            The season.
    """

    if not isinstance(mydate, datetime.datetime):
        raise TypeError("get_season requires a datetime.datetime structure.")

    year = mydate.year
    month = mydate.month

    if month < 8:
        year1 = (year - 1) % 100
    else:
        year1 = year % 100

    season = "%02i%02i" % (year1, year1 + 1)

    return season
