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

    main
"""
# Python Standard Library
import os
import argparse
import warnings

# Other Libraries
import crayons


def main():
    """
    Just print a welcoming message and calls the production_line_multiproc.
    """
    # Start with a welcome message.
    print("#" * 79)
    print("")
    print(" " * 25 + crayons.red("Raw radar PPIs production line.\n", bold=True))
    print("\t- Input data directory path is: " + crayons.yellow(INFILE))
    print("\t- Output data directory path is: " + crayons.yellow(OUTPATH))
    print("\t- Radiosounding directory path is: " + crayons.yellow(SOUND_DIR))
    if USE_UNRAVEL:
        print("\t- " + crayons.yellow("UNRAVEL") + " will be used as dealiasing algorithm.")
    else:
        print("\t- " + crayons.yellow("REGION-BASED") + " will be used as dealiasing algorithm.")
    print("\n" + "#" * 79 + "\n")

    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        import cpol_processing
        cpol_processing.process_and_save(INFILE, OUTPATH, SOUND_DIR, use_unravel=USE_UNRAVEL)

    print(crayons.green("Process completed."))

    return None


if __name__ == '__main__':
    """
    Global variables definition and logging file initialisation.
    """
    # Input directory for Radiosoundings (use my other script, named caprica to
    # download and format these datas).
    # INPATH = "/g/data/hj10/cpol_level_1a/ppi/"
    # OUTPATH = "/g/data/hj10/cpol_level_1b/"


    # Parse arguments
    parser_description = """Raw radar PPIs processing. It provides Quality
control, filtering, attenuation correction, dealiasing, unfolding, hydrometeors
calculation, and rainfall rate estimation."""
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
    parser.add_argument(
        '-r',
        '--radiosounds',
        dest='rs_dir',
        type=str,
        help='Radiosoundings directory.',
        default="/g/data2/rr5/CPOL_radar/DARWIN_radiosonde")

    parser.add_argument('--unravel', dest='unravel', action='store_true')
    parser.add_argument('--no-unravel', dest='unravel', action='store_false')
    parser.set_defaults(unravel=True)

    args = parser.parse_args()
    INFILE = args.infile
    OUTPATH = args.outdir
    SOUND_DIR =  args.rs_dir
    USE_UNRAVEL = args.unravel

    if not os.path.isfile(INFILE):
        parser.error("Invalid input file.")

    if not os.path.isdir(OUTPATH):
        parser.error("Invalid (or don't exist) output directory.")

    main()
