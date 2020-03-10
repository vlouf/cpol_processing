"""
Raw radar PPIs processing. Quality control, filtering, attenuation correction,
dealiasing, unfolding, hydrometeors calculation, rainfall rate estimation.
Tested on CPOL.

@title: cpol_processing
@author: Valentin Louf <valentin.louf@monash.edu>
@institution: Monash University
@date: 13/03/2019
@version: 2

.. autosummary::
    :toctree: generated/

    chunks
    main
"""
# Python Standard Library
import os
import sys
import glob
import argparse
import datetime
import traceback

import crayons


def welcome_message():
    """
    Display a welcome message with the input information.
    """
    print("#" * 79)
    print("")
    print(" " * 25 + crayons.red("Raw radar PPIs production line.\n", bold=True))
    print(" - Input data directory path is: " + crayons.yellow(INPATH))
    print(" - Output data directory path is: " + crayons.yellow(OUTPATH))
    print(" - Radiosounding directory path is: " + crayons.yellow(SOUND_DIR))
    print(f" - The process will occur between {crayons.yellow(START_DATE)} and {crayons.yellow(END_DATE)}.")
    if USE_UNRAVEL:
        print(" - " + crayons.yellow("UNRAVEL") + " will be used as dealiasing algorithm.")
    else:
        print(" - " + crayons.yellow("REGION-BASED") + " will be used as dealiasing algorithm.")
    print("\n" + "#" * 79 + "\n")


if __name__ == '__main__':
    """
    Global variables definition.
    """
    # Main global variables (Path directories).
    INPATH = "/g/data/hj10/admin/cpol_level_1a/v2019/ppi/"
    OUTPATH = '/g/data/hj10/admin/cpol_level_1b/v2020/'
    SOUND_DIR = "/g/data/kl02/vhl548/darwin_ancillary/DARWIN_radiosonde"

    # Parse arguments
    parser_description =  """Raw radar PPIs processing. It provides Quality
control, filtering, attenuation correction, dealiasing, unfolding, hydrometeors
calculation, and rainfall rate estimation."""
    parser = argparse.ArgumentParser(description=parser_description)
    parser.add_argument(
        '-s',
        '--start-date',
        dest='start_date',
        default=None,
        type=str,
        help='Starting date.',
        required=True)
    parser.add_argument(
        '-e',
        '--end-date',
        dest='end_date',
        default=None,
        type=str,
        help='Ending date.',
        required=True)
    parser.add_argument('--unravel', dest='unravel', action='store_true')
    parser.add_argument('--no-unravel', dest='unravel', action='store_false')
    parser.set_defaults(unravel=True)

    args = parser.parse_args()
    START_DATE = args.start_date
    END_DATE = args.end_date
    USE_UNRAVEL = args.unravel

    # Display infos
    welcome_message()