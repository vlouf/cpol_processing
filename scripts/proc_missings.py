"""
cpol_processing scripts for missing radar files in Radar archive on NCI.

@title: cpol_processing
@author: Valentin Louf <valentin.louf@bom.gov.au>
@institution: Bureau of Meteorology
@date: 26/03/2021

.. autosummary::
    :toctree: generated/

    chunks
    main
"""
import os
import glob
import traceback
from typing import Iterable, Any

import cpol_processing

from concurrent.futures import TimeoutError
from pebble import ProcessPool, ProcessExpired


def chunks(l: Any, n: int) -> Iterable[Any]:
    """
    Yield successive n-sized chunks from l.
    From http://stackoverflow.com/a/312464
    """
    for i in range(0, len(l), n):
        yield l[i : i + n]


def main(year: int) -> None:
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
    flist = glob.glob(os.path.join(INPATH, f"{year}/**/*.nc"))
    outlist = glob.glob(os.path.join(OUTPATH, f"v2021/ppi/{year}/**/*.nc"))

    oset = set([f[-18:-3] for f in outlist])
    iset = set([f[-18:-3] for f in flist])
    datelist = [*oset ^ iset]

    if len(datelist) == 0:
        print(f"No file to process for {YEAR}.")
        return None
    print(f"{year}: {len(datelist)} files to process.")

    inflist = []
    for d in datelist:
        inflist.append([f for f in flist if d in f][0])

    argslist = []
    for f in inflist:
        argslist.append((f, OUTPATH, SOUND_DIR, True))

    for fchunk in chunks(argslist, NCPUS):
        with ProcessPool() as pool:
            future = pool.starmap(cpol_processing.process_and_save, fchunk, timeout=360)
            iterator = future.result()

            while True:
                try:
                    _ = next(iterator)
                except StopIteration:
                    break
                except TimeoutError as error:
                    print("function took longer than %d seconds" % error.args[1])
                except ProcessExpired as error:
                    print("%s. Exit code: %d" % (error, error.exitcode))
                except TypeError:
                    continue
                except Exception:
                    traceback.print_exc()


if __name__ == "__main__":
    """
    Global variables definition.
    """
    INPATH = "/g/data/hj10/admin/cpol_level_1a/v2019/ppi/"
    OUTPATH = "/scratch/kl02/vhl548/cpol_level_1b/v2020/"
    SOUND_DIR = "/g/data/kl02/vhl548/darwin_ancillary/DARWIN_radiosonde"
    NCPUS = 16
    for YEAR in range(2009, 2018):
        main(YEAR)
