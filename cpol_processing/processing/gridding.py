"""
Codes for gridding radar data using Py-ART.

@title: gridding
@author: Valentin Louf <valentin.louf@monash.edu>
@institutions: Monash University and the Australian Bureau of Meteorology
@date: 20/11/2017

.. autosummary::
    :toctree: generated/

    _get_latlon
    gridding_radar_150km
    gridding_radar_70km
"""
# Python Standard Library
import os
import datetime

# Other Libraries
import pyart
import numpy as np


def _get_latlon(radgrid):
    """
    Generates lattitude and longitude arrays.

    Parameters:
    ===========
    radgrid: struct
        Py-ART grid object.

    Returns:
    ========
    longitude: ndarray
        Array of coordinates for all points.
    latitude: ndarray
        Array of coordinates for all points.
    """
    # Declare array, filled 0 in order to not have a masked array.
    lontot = np.zeros_like(radgrid.fields['reflectivity']['data'].filled(0))
    lattot = np.zeros_like(radgrid.fields['reflectivity']['data'].filled(0))

    for lvl in range(radgrid.nz):
        lontot[lvl, :, :], lattot[lvl, :, :] = radgrid.get_point_longitude_latitude(lvl)

    longitude = pyart.config.get_metadata('longitude')
    latitude = pyart.config.get_metadata('latitude')

    longitude['data'] = lontot
    latitude['data'] = lattot

    return longitude, latitude


def gridding_radar_150km(radar, radar_date, outpath):
    """
    Map a single radar to a Cartesian grid of 150 km range and 2.5 km resolution.

    Parameters:
    ===========
        radar:
            Py-ART radar structure.
        radar_date: datetime
            Datetime stucture of the radar data.
        outpath: str
            Ouput directory.
    """
    # Extracting year, date, and datetime.
    year = str(radar_date.year)
    datestr = radar_date.strftime("%Y%m%d")
    datetimestr = radar_date.strftime("%Y%m%d_%H%M")
    fname = "CPOL_{}_GRIDS_2500m.nc".format(datetimestr)

    # Output directory
    outdir_150km = os.path.join(outpath, "GRID_150km_2500m")
    try:
        os.mkdir(outdir_150km)
    except FileExistsError:
        pass

    outdir_150km = os.path.join(outdir_150km, year)
    try:
        os.mkdir(outdir_150km)
    except FileExistsError:
        pass

    outdir_150km = os.path.join(outdir_150km, datestr)
    try:
        os.mkdir(outdir_150km)
    except FileExistsError:
        pass

    # Output file name
    outfilename = os.path.join(outdir_150km, fname)

    # exclude masked gates from the gridding
    my_gatefilter = pyart.filters.GateFilter(radar)
    my_gatefilter.exclude_transition()
    my_gatefilter.exclude_masked('reflectivity')

    # Gridding
    grid_150km = pyart.map.grid_from_radars(
        radar, gatefilters=my_gatefilter,
        grid_shape=(41, 117, 117),
        grid_limits=((0, 20000), (-145000.0, 145000.0), (-145000.0, 145000.0)),
        roi_func='constant', constant_roi=2500)

    # Latitude Longitude field for each point.
    longitude, latitude = _get_latlon(grid_150km)
    grid_150km.add_field('longitude', longitude)
    grid_150km.add_field('latitude', latitude)

    # Saving data.
    grid_150km.write(outfilename, arm_time_variables=True)

    return None


def gridding_radar_70km(radar, radar_date, outpath):
    """
    Map a single radar to a Cartesian grid of 70 km range and 1 km resolution.

    Parameters:
    ===========
        radar:
            Py-ART radar structure.
        my_gatefilter:
            The Gate filter.
        radar_date: datetime
            Datetime stucture of the radar data.
        outpath: str
            Ouput directory.
    """
    # Extracting year, date, and datetime.
    year = str(radar_date.year)
    datestr = radar_date.strftime("%Y%m%d")
    datetimestr = radar_date.strftime("%Y%m%d_%H%M")
    fname = "CPOL_{}_GRIDS_1000m.nc".format(datetimestr)

    # Output directory
    outdir_70km = os.path.join(outpath, "GRID_70km_1000m")
    try:
        os.mkdir(outdir_70km)
    except FileExistsError:
        pass

    outdir_70km = os.path.join(outdir_70km, year)
    try:
        os.mkdir(outdir_70km)
    except FileExistsError:
        pass

    outdir_70km = os.path.join(outdir_70km, datestr)
    try:
        os.mkdir(outdir_70km)
    except FileExistsError:
        pass

    # Output file name
    outfilename = os.path.join(outdir_70km, fname)

    # exclude masked gates from the gridding
    my_gatefilter = pyart.filters.GateFilter(radar)
    my_gatefilter.exclude_transition()
    my_gatefilter.exclude_masked('reflectivity')

    # Gridding
    grid_70km = pyart.map.grid_from_radars(
        radar, gatefilters=my_gatefilter,
        grid_shape=(41, 141, 141),
        grid_limits=((0, 20000), (-70000.0, 70000.0), (-70000.0, 70000.0)),
        roi_func='constant', constant_roi=1000)

    # Latitude Longitude field for each point.
    longitude, latitude = _get_latlon(grid_70km)
    grid_70km.add_field('longitude', longitude)
    grid_70km.add_field('latitude', latitude)

    # Saving data.
    grid_70km.write(outfilename, arm_time_variables=True)

    return None
