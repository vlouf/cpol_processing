"""
CPOL Level 1b main production line.

@title: CPOL_PROD_1b
@author: Valentin Louf <valentin.louf@monash.edu>
@institution: Bureau of Meteorology
@date: 1/03/2019
@version: 1

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
import gc
import sys
import time
import glob
import signal
import logging
import argparse
import datetime
import warnings
import traceback

import pandas as pd
# import dask.bag as db
from multiprocessing import Pool
import cpol_processing


class TimeoutException(Exception):   # Custom exception class
    pass


def timeout_handler(signum, frame):   # Custom signal handler
    raise TimeoutException


def production_line_manager(radar_file_name):
    """
    The production line manager calls the production line and manages it ;-).
    Buffer function that is used to catch any problem with the processing line
    without screwing the whole multiprocessing stuff.

    Parameters:
    ===========
        radar_file_name: str
            Name of the input radar file.
        outpath: str
            Path for saving output data.
    """
    signal.signal(signal.SIGALRM, timeout_handler)
    signal.alarm(TIME_BEFORE_DEATH)
    try:
        cpol_processing.process_and_save(radar_file_name, OUTPATH, sound_dir=SOUND_DIR)
    except TimeoutException:
        # Treatment time was too long.
        logging.error("Too much time taken to treat %s, killing process.", radar_file_name)
        return None  # Go to next iteration.
    except Exception:
        print("Exception in production line code:")
        print("-" * 60)
        print("ERROR IN FILE {}.".format(radar_file_name))
        traceback.print_exc(file=sys.stdout)
        print("-" * 60)
        logging.error("Failed to process file", exc_info=True)
        return None
    else:
        signal.alarm(0)

    return None


def main():
    date_list = pd.date_range(START_DATE, END_DATE)
    for day in date_list:
        input_dir = os.path.join(INPATH, str(day.year), day.strftime("%Y%m%d"), "*.*")
        flist = sorted(glob.glob(input_dir))
        if len(flist) == 0:
            print('No file found for {}.'.format(day.strftime("%Y-%b-%d")))
            continue
        print(f'{len(flist)} files found for ' + day.strftime("%Y-%b-%d"))

        # bag = db.from_sequence(flist).map(production_line_manager)
        # bag.compute()
        with Pool(16) as pool:
            pool.map(production_line_manager, flist)
        gc.collect()

    return None


if __name__ == '__main__':
    """
    Global variables definition and logging file initialisation.
    """
    # Main global variables (Path directories).
    INPATH = "/g/data/hj10/cpol_level_1a/ppi/"
    OUTPATH = "/g/data/hj10/cpol_level_1b/"
    SOUND_DIR = "/g/data2/rr5/CPOL_radar/DARWIN_radiosonde"
    LOG_FILE_PATH = os.path.expanduser('~')
    TIME_BEFORE_DEATH = 600  # seconds before killing process.

    # Parse arguments
    parser_description = "Processing of radar data from level 1a to level 1b."
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

    args = parser.parse_args()
    START_DATE = args.start_date
    END_DATE = args.end_date

    # Creating the general log file.
    logname = "cpol_level1b_from_{}_to_{}.log".format(START_DATE, END_DATE)
    log_file_name = os.path.join(LOG_FILE_PATH, logname)
    logging.basicConfig(
        level=logging.DEBUG,
        format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
        filename=log_file_name,
        filemode='w+')
    logger = logging.getLogger(__name__)

    main()
