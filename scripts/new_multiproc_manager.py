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
import logging
import argparse
import datetime
import warnings
import traceback

import pandas as pd
import dask.bag as db
import cpol_processing


def production_line_manager(infile):
    """
    The production line manager calls the production line and manages it ;-).
    Buffer function that is used to catch any problem with the processing line
    without screwing the whole multiprocessing stuff.

    Parameters:
    ===========
    infile: str
        Name of the input radar file.
    outpath: str
        Path for saving output data.
    """    
    try:
        cpol_processing.process_and_save(radar_filinfilee_name, OUTPATH, sound_dir=SOUND_DIR)    
    except Exception:        
        traceback.print_exc(file=sys.stdout)        
        logging.error(f"Failed to process {infile}", exc_info=True)        

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

        bag = db.from_sequence(flist).map(production_line_manager)
        bag.compute()
        # with Pool(16) as pool:
        #     pool.map(production_line_manager, flist)
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
    LOG_FILE_PATH = "/short/kl02/vhl548/"    

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
        level=logging.WARNING,
        format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
        filename=log_file_name,
        filemode='w+')
    logger = logging.getLogger(__name__)

    main()
