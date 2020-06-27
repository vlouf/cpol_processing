import os
import glob
import argparse
import traceback

import cpol_processing

from concurrent.futures import TimeoutError
from pebble import ProcessPool, ProcessExpired


def buffer(infile):
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
    try:
        cpol_processing.process_and_save(infile,
                                         OUTPATH,
                                         sound_dir=SOUND_DIR,
                                         do_dealiasing=True,
                                         use_unravel=True)
    except Exception:
        traceback.print_exc()

    return None


def chunks(l, n):
    """
    Yield successive n-sized chunks from l.
    From http://stackoverflow.com/a/312464
    """
    for i in range(0, len(l), n):
        yield l[i:i + n]


def main():
    year = YEAR
    flist = glob.glob(f'/g/data/hj10/admin/cpol_level_1a/v2019/ppi/{year}/**/*.nc')

    oset = set([f[-18:-3] for f in glob.glob(f'/scratch/kl02/vhl548/cpol_level_1b/v2020/v2020/ppi/{year}/**/*.nc')])
    iset = set([f[-18:-3] for f in glob.glob(f'/g/data/hj10/admin/cpol_level_1a/v2019/ppi/{year}/**/*.nc')])

    datelist = [*oset ^ iset]
    if len(datelist) == 0:
        print('No file to process.')
        return None
    print(f'{year}: {len(datelist)} files to process.')

    inflist = []
    for d in datelist:
        inflist.append([f for f in flist if d in f][0])

    for fchunk in chunks(inflist, NCPUS):
        with ProcessPool() as pool:
            future = pool.map(buffer, fchunk, timeout=240)
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
    INPATH = "/g/data/hj10/admin/cpol_level_1a/v2019/ppi/"
    OUTPATH = '/scratch/kl02/vhl548/cpol_level_1b/v2020/'
    SOUND_DIR = "/g/data/kl02/vhl548/darwin_ancillary/DARWIN_radiosonde"
    YEAR = 2012
    NCPUS = 16
    main()