"""
CPOL Level 1b main production line.

@title: CPOL_PROD_1b
@author: Valentin Louf <valentin.louf@monash.edu>
@institution: Bureau of Meteorology
@date: 31/05/2018
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
import sys
import time
import signal
import logging
import argparse
import datetime
import warnings
import traceback
from multiprocessing import Pool

# Other Libraries -- Matplotlib must be imported first
import matplotlib
matplotlib.use('Agg')  # <- Reason why matplotlib is imported first.

import crayons  # For the welcoming message only.
import numpy as np
import pandas as pd  # Using only 1 function from pandas, should drop dependency.

# Custom modules.
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


def make_dir(input_dir):
    """
    Make directory.

    Parameter:
    ==========
    input_dir: str
        Input directory name.
    """
    try:
        os.mkdir(input_dir)
        print("%s created." % (input_dir))
    except FileExistsError:
        pass

    return None


def production_line_manager(radar_file_name, outpath, outpath_grid, figure_path, sound_dir):
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
    shortname = os.path.basename(radar_file_name)
    # If we are stuck in the loop for more than TIME_BEFORE_DEATH seconds,
    # it will raise a TimeoutException that will kill the current process
    # and go to the next iteration.
    signal.signal(signal.SIGALRM, timeout_handler)
    signal.alarm(TIME_BEFORE_DEATH)
    try:
        cpol_processing.process_and_save(radar_file_name, outpath, outpath_grid,
                                         figure_path, sound_dir, IS_CPOL)
    except TimeoutException:
        # Treatment time was too long.
        logging.error("Too much time taken to treat %s, killing process.", shortname)
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


def production_line_multiproc(mydate):
    """
    It makes sure that input/output directories exist. This is where the
    multiprocessing is taken care of.

    INPATH and OUTPATH are global variables defined in __main__.

    Parameter:
    ==========
        mydate: datetime.datetime
            Radar data date for which we start the production.
    """
    year = str(mydate.year)
    datestr = mydate.strftime("%Y%m%d")
    indir = os.path.join(INPATH, year, datestr)
    if not IS_CPOL:
        indir = INPATH

    # Checking if input directory exists.
    if not os.path.exists(indir):
        return None

    # Checking if output directory exists. Creating them otherwise.
    outdir = os.path.join(OUTPATH, year)
    make_dir(outdir)  # Function handles itself.
    outdir = os.path.join(outdir, datestr)
    make_dir(outdir)  # Function handles itself.

    # List netcdf files in directory.
    flist = get_files(indir, mydate)
    if len(flist) == 0:
        logger.error('%s empty.', indir)
        return None
    logger.info('%i files found for %s', len(flist), datestr)

    # Cutting the file list into smaller chunks. (The multiprocessing.Pool instance
    # is freed from memory, at each iteration of the main for loop).
    for flist_slice in chunks(flist, NCPU):
        # Because we use multiprocessing, we need to send a list of tuple as argument of Pool.starmap.
        args_list = [(onefile, outdir, OUTPATH_GRID, FIGURE_CHECK_PATH, SOUND_DIR) for onefile in flist_slice]
        # Start multiprocessing.
        with Pool(NCPU) as pool:
            pool.starmap(production_line_manager, args_list)

    return None


def main():
    """
    Just print a welcoming message and calls the production_line_multiproc.
    """
    # Start with a welcome message.
    print("#" * 79)
    print("")
    print(" " * 25 + crayons.red("CPOL Level 1b production line.", bold=True))
    print("")
    print("- Input data directory path is: " + crayons.yellow(INPATH))
    print("- Output data directory path is: " + crayons.yellow(OUTPATH))
    print("- Radiosounding directory path is: " + crayons.yellow(SOUND_DIR))
    print("- Figures will be saved in: " + crayons.yellow(FIGURE_CHECK_PATH))
    print("- Start date is: " + crayons.yellow(START_DATE))
    print("- End date is: " + crayons.yellow(END_DATE))
    print("- Log files can be found in: " + crayons.yellow(LOG_FILE_PATH))
    print("- Each subprocess has {}s of allowed time life before being killed.".format(TIME_BEFORE_DEATH))
    print("#" * 79)
    print("")

    # Serious stuffs begin here.
    date_range = pd.date_range(START_DATE, END_DATE)
    # One date at a time.
    for thedate in date_range:
        production_line_multiproc(thedate)

    return None


if __name__ == '__main__':
    """
    Global variables definition and logging file initialisation.
    """
    # Main global variables (Path directories).
    # Input radar data directory
    INPATH = "/g/data2/rr5/vhl548/CPOL_level_1a"
    # Output directory for CF/Radial PPIs
    OUTPATH = "/g/data2/rr5/vhl548/NEW_CPOL_level_1b/"
    # Input directory for Radiosoundings
    SOUND_DIR = "/g/data2/rr5/vhl548/DARWIN_radiosonde"
    # Output directory for log files.
    LOG_FILE_PATH = os.path.expanduser('~')
    # Time in seconds for which each subprocess is allowed to live.
    TIME_BEFORE_DEATH = 300  # seconds before killing process.
    # True (or False) it is (or not) CPOL radar:
    IS_CPOL = True

    # Output directory for verification figures.
    FIGURE_CHECK_PATH = os.path.join(OUTPATH, 'FIGURE_CHECK')
    # Output directory for GRIDDED netcdf data.
    OUTPATH_GRID = os.path.join(OUTPATH, 'GRIDDED')
    if "logfiles" not in LOG_FILE_PATH:
        LOG_FILE_PATH = os.path.join(LOG_FILE_PATH, "logfiles")

    # Check if paths exist.
    make_dir(OUTPATH)
    OUTPATH = os.path.join(OUTPATH, "PPI")
    make_dir(OUTPATH)
    
    # Create input directory.
    make_dir(LOG_FILE_PATH)
    make_dir(OUTPATH_GRID)
    make_dir(os.path.join(OUTPATH_GRID, 'GRID_150km_2500m'))
    make_dir(os.path.join(OUTPATH_GRID, 'GRID_70km_1000m'))
    make_dir(FIGURE_CHECK_PATH)

    if not os.path.isdir(SOUND_DIR):
        raise FileNotFoundError("Radiosoundings directory does not exist (or invalid): {}.".format(SOUND_DIR))

    if not os.path.isdir(INPATH):
        raise FileNotFoundError("Input data directory does not exist (or invalid): {}.".format(INPATH))

    # Parse arguments
    parser_description = "Leveling treatment of CPOL data from level 1a to level 1b."
    parser = argparse.ArgumentParser(description=parser_description)
    parser.add_argument(
        '-j',
        '--cpu',
        dest='ncpu',
        default=16,
        type=int,
        help='Number of process')
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
    NCPU = args.ncpu
    START_DATE = args.start_date
    END_DATE = args.end_date

    # Checking that dates are recognize.
    try:
        datetime.datetime.strptime(START_DATE, "%Y%m%d")
        datetime.datetime.strptime(END_DATE, "%Y%m%d")
    except Exception:
        print("Did not understand the date format. Must be YYYYMMDD.")
        sys.exit()

    # Creating the general log file.
    logname = "cpol_level1b_from_{}_to_{}.log".format(START_DATE, END_DATE)
    log_file_name = os.path.join(LOG_FILE_PATH, logname)
    logging.basicConfig(
        level=logging.DEBUG,
        format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
        filename=log_file_name,
        filemode='w+')
    logger = logging.getLogger(__name__)

    with warnings.catch_warnings():
        # Just ignoring warning messages.
        warnings.simplefilter("ignore")
        main()
