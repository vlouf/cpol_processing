"""
CPOL Level 1b main production line.

@title: CPOL_PROD_1b
@author: Valentin Louf <valentin.louf@monash.edu>
@institution: Bureau of Meteorology
@date: 31/05/2017
@version: 0.99

.. autosummary::
    :toctree: generated/

    timeout_handler
    chunks
    production_line_manager
    production_line_multiproc
    main
"""
# Python Standard Library
import os
import logging
import argparse
import datetime
import warnings

# Other Libraries -- Matplotlib must be imported first
import matplotlib
matplotlib.use('Agg')  # <- Reason why matplotlib is imported first.

import pyart
import crayons  # For the welcoming message only.
import numpy as np

# Custom modules.
import cpol_processing


def main():
    """
    Just print a welcoming message and calls the production_line_multiproc.
    """
    # Start with a welcome message.
    print("#" * 79)
    print("")
    print(" " * 25 + crayons.red("CPOL Level 1b production line.", bold=True))
    print("")
    print("- Input data directory path is: " + crayons.yellow(INFILE))
    print("- Output data directory path is: " + crayons.yellow(OUTPATH))
    print("- Radiosounding directory path is: " + crayons.yellow(SOUND_DIR))
    print("- Figures will be saved in: " + crayons.yellow(FIGURE_CHECK_PATH))
    print("#" * 79)
    print("")

    # Serious stuffs begin here.
    cpol_processing.production_line(INFILE, OUTPATH, OUTPATH_GRID, FIGURE_CHECK_PATH, SOUND_DIR)

    return None


if __name__ == '__main__':
    """
    Global variables definition and logging file initialisation.
    """
    # Input directory for Radiosoundings (use my other script, named caprica to
    # download and format these datas).
    SOUND_DIR = "/g/data2/rr5/vhl548/DARWIN_radiosonde/"
    # Output directory for verification figures.
    # Output directory for log files.
    LOG_FILE_PATH = "/short/kl02/vhl548/logfiles/"

    # Parse arguments
    parser_description = "Leveling treatment of CPOL data from level 1a to level 1b."
    parser = argparse.ArgumentParser(description=parser_description)
    parser.add_argument(
        '-i',
        '--input',
        dest='infile',
        type=str,
        help='Input file',
        required=True)
    parser.add_argument(
        '-o',
        '--output',
        dest='outdir',
        type=str,
        help='Output directory.',
        required=True)

    args = parser.parse_args()
    INFILE = args.infile
    OUTPATH = args.outdir

    if not os.path.isfile(INFILE):
        parser.error("Invalid input file.")

    if not os.path.isdir(OUTPATH):
        parser.error("Invalid (or don't exist) output directory.")

    OUTPATH_GRID = os.path.join(OUTPATH, 'GRIDDED')
    FIGURE_CHECK_PATH = os.path.join(OUTPATH, 'FIGURE_CHECK')
    if not os.path.isdir(OUTPATH_GRID):
        print("Creating output figures directory: {}.".format(OUTPATH_GRID))
        os.mkdir(OUTPATH_GRID)
    if not os.path.isdir(FIGURE_CHECK_PATH):
        print("Creating output figures directory: {}.".format(FIGURE_CHECK_PATH))
        os.mkdir(FIGURE_CHECK_PATH)

    # Creating the general log file.
    logname = "log_file_for_{}.log".format(os.path.basename(INFILE))
    log_file_name = os.path.join(LOG_FILE_PATH, logname)
    logging.basicConfig(
        level=logging.DEBUG,
        format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
        filename=log_file_name,
        filemode='a+')
    logger = logging.getLogger(__name__)

    with warnings.catch_warnings():
        # Just ignoring warning messages.
        warnings.simplefilter("ignore")
        main()
