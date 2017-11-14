cpol_processing
===============

Radar processing script. This script processes raw radar files using various tools like [Py-ART][1], [csu_radartools][2], as well as some custom function.

- Compute SNR using radiosoundings and map the temperature to all radar gates.
- Correct RHOHV from range dependancy using Ryzhkov algorithm.
- Correct ZDR using Ryzhkov algorithm.
- Create a filter to remove noise and incorrect data.
- Process and unfold raw PHIDP using Bringi's technique and Giangrande's LP processing.
- Compute KDP.
- Unfold velocity using [Py-ART][1].
- Correct reflectivity attenuation (C-band only).
- Correct differential reflectivity attenuation (C-band only).
- Estimate Hydrometeors classification using [csu_radartools][2].
- Estimate Rainfall rate using [csu_radartools][2].
- Estimate DSD retrieval using [csu_radartools][2].

# Requirements

Radiosoundings are required. Ideally radiosoundings coming from arm.gov.
You can download radiosoundings from the University of Wyoming website (http://weather.uwyo.edu/upperair/sounding.html) and save them into a netCDF4 file.

# Install

To install cpol_processing, you can either download and unpack the zip file of the source code or use git to checkout the repository:

`git clone git@github.com:vlouf/cpol_processing.git`

To install in your home directory, use:

`python setup.py install --user`

# Dependencies:
- [py-art][1]
- [csu_radartools][2]
- [netCDF4][3]
- [numpy][4]
- [scipy][5]

[1]: http://github.com/ARM-DOE/pyart
[2]: http://github.com/CSU-Radarmet/CSU_RadarTools
[3]: http://unidata.github.io/netcdf4-python/
[4]: http://www.scipy.org/
[5]: http://www.scipy.org/
