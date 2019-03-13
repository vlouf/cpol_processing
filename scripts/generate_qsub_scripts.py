import os
import glob
import datetime
import calendar


def configuration_file(start_date='19990101', end_date='19990131', walltime=10):
    conf_txt = """#!/bin/bash
#PBS -P kl02
#PBS -q normal
#PBS -l walltime={time}:00:00
#PBS -l mem=32GB
#PBS -l wd
#PBS -l ncpus=16
#PBS -lother=gdata1
conda activate radar
python new_multiproc_manager.py -s {sdate} -e {edate}
""".format(time=walltime, sdate=start_date, edate=end_date)
    if walltime <= 7:
        conf_txt = conf_txt.replace("normal", "express")

    return conf_txt

for year in range(1997, 2018):
    for month in range(1, 13):
        if month > 7 and month < 10:
            continue

        if month < 8:
            season_st = ((year - 1) % 100)
            season_nd = ((year) % 100)
        else:
            season_st = (year % 100)
            season_nd = ((year + 1) % 100)
        season = "%02i%02i" % (season_st, season_nd)
 
        indir = "/g/data/hj10/cpol_level_1a/ppi"
        indir += "/%i/%i%02i" % (year, year, month)
        dirlist = glob.glob(indir + "*")
        print(indir)
        if len(dirlist) == 0:
            continue

        _, ed = calendar.monthrange(year, month)
        sdatestr = "%i%02i%02i" % (year, month, 1)
        edatestr = "%i%02i%02i" % (year, month, ed)
        if month == 10 or month == 5:
            f = configuration_file(sdatestr, edatestr, 5)
        elif month == 11:
            f = configuration_file(sdatestr, edatestr, 7)
        else:
            f = configuration_file(sdatestr, edatestr, 10)


        fname = "qlevel1b_%i%02i.pbs" % (year, month)
        with open(fname, 'w') as fid:
            fid.write(f)
