import os
import glob
import datetime
import calendar


def configuration_file(ncpu=16, start_date='19990101', end_date='19990131', walltime=10):
    conf_txt = """#!/bin/bash
#PBS -P kl02
#PBS -q normal
#PBS -l walltime={time}:00:00
#PBS -l mem={mem}GB
#PBS -l wd
#PBS -l ncpus={cpu}
#PBS -lother=gdata2
source activate radar
python rajin_multiproc_processing.py -s {sdate} -e {edate} -j {cpu}
""".format(time=walltime, cpu=ncpu, mem=int(ncpu*2), sdate=start_date, edate=end_date)
    if walltime <= 7:
        conf_txt = conf_txt.replace("normal", "express")

    return conf_txt

for year in range(1998, 2018):
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
 
        indir = "/g/data2/rr5/vhl548/CPOL_level_0"
        indir += "/s%s/%i%02i" % (season, year, month)
        dirlist = glob.glob(indir + "*")
        print(indir)
        if len(dirlist) == 0:
            continue

        _, ed = calendar.monthrange(year, month)
        sdatestr = "%i%02i%02i" % (year, month, 1)
        edatestr = "%i%02i%02i" % (year, month, ed)
        if month == 10 or month == 5:
            f = configuration_file(16, sdatestr, edatestr, 5)
        elif month == 11:
            f = configuration_file(16, sdatestr, edatestr, 7)
        else:
            f = configuration_file(16, sdatestr, edatestr, 10)


        fname = "qlevel1b_%i%02i.pbs" % (year, month)
        with open(fname, 'w') as fid:
            fid.write(f)
