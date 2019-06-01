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

from concurrent.futures import TimeoutError
from pebble import ProcessPool, ProcessExpired


def chunks(l, n):
    """
    Yield successive n-sized chunks from l.
    From http://stackoverflow.com/a/312464
    """
    for i in range(0, len(l), n):
        yield l[i:i + n]


def main(inargs):
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
    import warnings
    import traceback

    infile, outpath, sound_dir, use_unravel = inargs

    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        import cpol_processing

    try:
        cpol_processing.process_and_save(infile, outpath, sound_dir=sound_dir, use_unravel=use_unravel)
    except Exception:
        traceback.print_exc()
        return None

    return None


def welcome_message():
    """
    Display a welcome message with the input information.
    """
    print("#" * 79)
    print("")
    print(" " * 25 + crayons.red("Raw radar PPIs production line.\n", bold=True))
    print("\t- Input data directory path is: " + crayons.yellow(INPATH))
    print("\t- Output data directory path is: " + crayons.yellow(OUTPATH))
    print("\t- Radiosounding directory path is: " + crayons.yellow(SOUND_DIR))
    print(f"\t- The process will occur between {crayons.yellow(START_DATE)} and {crayons.yellow(END_DATE)}.")
    if USE_UNRAVEL:
        print("\t- " + crayons.yellow("UNRAVEL") + " will be used as dealiasing algorithm.")
    else:
        print("\t- " + crayons.yellow("REGION-BASED") + " will be used as dealiasing algorithm.")
    print("\n" + "#" * 79 + "\n")


if __name__ == '__main__':
    """
    Global variables definition.
    """
    # Main global variables (Path directories).
    INPATH = "/g/data/hj10/cpol_level_1a/ppi/"
    OUTPATH = "/g/data/hj10/cpol_level_1b/"
    SOUND_DIR = "/g/data2/rr5/CPOL_radar/DARWIN_radiosonde"

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

    # Check date
    try:
        start = datetime.datetime.strptime(START_DATE, "%Y%m%d")
        end = datetime.datetime.strptime(END_DATE, "%Y%m%d")
        if start > end:
            parser.error('End date older than start date.')
        date_range = [start + datetime.timedelta(days=x) for x in range(0, (end - start).days + 1, )]
    except ValueError:
        parser.error('Invalid dates.')
        sys.exit()

    for day in date_range:
        input_dir = os.path.join(INPATH, str(day.year), day.strftime("%Y%m%d"), "*.*")
        flist = sorted(glob.glob(input_dir))
        if len(flist) == 0:
            print('No file found for {}.'.format(day.strftime("%Y-%b-%d")))
            continue

        print(f'{len(flist)} files found for ' + day.strftime("%Y-%b-%d"))

        for flist_chunk in chunks(flist, 16):
            arglist = [(f, OUTPATH, SOUND_DIR, USE_UNRAVEL) for f in flist_chunk]

            with ProcessPool() as pool:
                future = pool.map(main, arglist, timeout=240)
                iterator = future.result()

                while True:
                    try:
                        result = next(iterator)
                    except StopIteration:
                        break
                    except TimeoutError as error:
                        print("function took longer than %d seconds" % error.args[1])
                    except ProcessExpired as error:
                        print("%s. Exit code: %d" % (error, error.exitcode))
                    except Exception as error:
                        print("function raised %s" % error)
                        print(error.traceback)  # Python's traceback of remote process
