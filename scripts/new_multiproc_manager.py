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
import glob
import signal
import argparse
import traceback

from multiprocessing import Pool

import pandas as pd
import cpol_processing


class TimeoutException(Exception):   # Custom exception class
    pass


def timeout_handler(signum, frame):   # Custom signal handler
    raise TimeoutException


def chunks(l, n):
    """
    Yield successive n-sized chunks from l.
    From http://stackoverflow.com/a/312464
    """
    for i in range(0, len(l), n):
        yield l[i:i + n]


def main(infile):
    """
    It calls the production line and manages it. Buffer function that is used
    to catch any problem with the processing line without screwing the whole
    multiprocessing stuff.

    Parameters:
    ===========
    infile: str
        Name of the input radar file.
    outpath: str
        Path for saving output data.
    """
    # SIGALRM is unix only.
    signal.signal(signal.SIGALRM, timeout_handler)
    signal.alarm(TIME_BEFORE_DEATH)
    try:
        cpol_processing.process_and_save(infile, OUTPATH, sound_dir=SOUND_DIR)
    except TimeoutException:
        print(f'Process deadlock for file {infile}. Killing it.')
    except Exception:
        traceback.print_exc()
    else:
        signal.alarm(0)

    return None


if __name__ == '__main__':
    """
    Global variables definition.
    """
    # Main global variables (Path directories).
    INPATH = "/g/data/hj10/cpol_level_1a/ppi/"
    OUTPATH = "/g/data/hj10/cpol_level_1b/"
    SOUND_DIR = "/g/data2/rr5/CPOL_radar/DARWIN_radiosonde"
    LOG_FILE_PATH = "/short/kl02/vhl548/"
    TIME_BEFORE_DEATH = 240

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

    for day in pd.date_range(START_DATE, END_DATE):
        input_dir = os.path.join(INPATH, str(day.year), day.strftime("%Y%m%d"), "*.*")
        flist = sorted(glob.glob(input_dir))
        if len(flist) == 0:
            print('No file found for {}.'.format(day.strftime("%Y-%b-%d")))
            continue
        print(f'{len(flist)} files found for ' + day.strftime("%Y-%b-%d"))

        for flist_chunk in chunks(flist, 32):
            with Pool(16) as pool:
                pool.map(main, flist)
