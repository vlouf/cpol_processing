"""
Raw radar PPIs processing. Quality control, filtering, attenuation correction,
dealiasing, unfolding, hydrometeors calculation, rainfall rate estimation.

@title: cpol_processing
@author: Valentin Louf <valentin.louf@bom.gov.au>
@institution: Monash University and Bureau of Meteorology
@date: 10/02/2021

.. autosummary::
    :toctree: generated/

    chunks
    main
    welcome_message
"""
# Python Standard Library
import gc
import os
import sys
import glob
import argparse
import datetime

import pandas as pd
import dask.bag as db
import crayons
import cpol_processing


def chunks(l, n):
    """
    Yield successive n-sized chunks from l.
    From http://stackoverflow.com/a/312464
    """
    for i in range(0, len(l), n):
        yield l[i : i + n]


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


def main(start, end):
    date_range = pd.date_range(start, end)
    for day in date_range:
        input_dir = os.path.join(INPATH, str(day.year), day.strftime("%Y%m%d"), "*.*")
        flist = sorted(glob.glob(input_dir))
        if len(flist) == 0:
            print("No file found for {}.".format(day.strftime("%Y-%b-%d")))
            continue

        argslist = []
        for f in flist:
            argslist.append((f, OUTPATH, SOUND_DIR, DO_DEALIASING, USE_UNRAVEL))

        print(f"{len(flist)} files found for " + day.strftime("%Y-%b-%d"))

        for chunk in chunks(argslist, 32):
            bag = db.from_sequence(chunk).starmap(cpol_processing.process_and_save)
            _ = bag.compute()
        
        del bag
        gc.collect()
    return None


if __name__ == "__main__":
    """
    Global variables definition.
    """
    # Main global variables (Path directories).
    INPATH = "/g/data/hj10/admin/cpol_level_1a/v2019/ppi/"
    OUTPATH = "/scratch/kl02/vhl548/cpol_level_1b/v2020/"
    SOUND_DIR = "/g/data/kl02/vhl548/darwin_ancillary/DARWIN_radiosonde"

    # Parse arguments
    parser_description = """Raw radar PPIs processing. It provides Quality
control, filtering, attenuation correction, dealiasing, unfolding, hydrometeors
calculation, and rainfall rate estimation."""
    parser = argparse.ArgumentParser(description=parser_description)
    parser.add_argument(
        "-s", "--start-date", dest="start_date", default=None, type=str, help="Starting date.", required=True
    )
    parser.add_argument("-e", "--end-date", dest="end_date", default=None, type=str, help="Ending date.", required=True)
    parser.add_argument("--unravel", dest="unravel", action="store_true")
    parser.add_argument("--no-unravel", dest="unravel", action="store_false")
    parser.set_defaults(unravel=True)
    parser.add_argument("--dealias", dest="dealias", action="store_true")
    parser.add_argument("--no-dealias", dest="dealias", action="store_false")
    parser.set_defaults(dealias=True)

    args = parser.parse_args()
    START_DATE = args.start_date
    END_DATE = args.end_date
    USE_UNRAVEL = args.unravel
    DO_DEALIASING = args.dealias

    # Check date
    try:
        start = datetime.datetime.strptime(START_DATE, "%Y%m%d")
        end = datetime.datetime.strptime(END_DATE, "%Y%m%d")
        if start > end:
            parser.error("End date older than start date.")        
    except ValueError:
        parser.error("Invalid dates.")
        sys.exit()

    # Display infos
    welcome_message()
    main(start, end)
