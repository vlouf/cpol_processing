"""
Gridding radar data on a Cartesian grid using Py-ART.

@title: gridding
@author: Valentin Louf <valentin.louf@monash.edu>
@institutions: Monash University and the Australian Bureau of Meteorology
@date: 11/12/2018

.. autosummary::
    :toctree: generated/

    mkdir
    _get_latlon
    gridding_radar
"""
# Python Standard Library
import os
import datetime

# Other Libraries
import pyart
import numpy as np


def mkdir(dirpath):
    '''
    Create directory. Check if directory exists and handles error.
    '''
    if not os.path.exists(dirpath):
        # Might seem redundant, but the multiprocessing creates error.
        try:
            os.mkdir(dirpath)
        except FileExistsError:
            return None

    return None


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


def gridding_radar(radar, radar_date, outpath, rmax=145e3, xyres=1000,
                   maxheight=20e3, zres=500, azimuthal_res=1, refl_name='reflectivity',
                   linearz=False):
    """
    Map a single radar to a Cartesian grid of 150 km range and 1 km resolution.

    Parameters:
    ===========
        radar:
            Py-ART radar structure.
        radar_date: datetime
            Datetime stucture of the radar data.
        outpath: str
            Ouput directory.
        rmax: int
            Maximum range in meters
        xyres: int
            Grid x/y resolution in meters
        maxheight: int
            Maximum altitude of the grid in meters.
        zres: int
            Grid z-axis resolution in meters
        azimuthal_res: float
            Beamwidth angle in degree.
        refl_name: str
            Name of the reflectivity field
        linearz: bool
            Use linear-Z reflectivity for the gridding.
    """
    # Extracting year, date, and datetime.
    year = str(radar_date.year)
    datestr = radar_date.strftime("%Y%m%d")
    datetimestr = radar_date.strftime("%Y%m%d_%H%M")
    fname = "cpol_{}_grids_{}m.nc".format(datetimestr, resolution)

    # Output directory
    outdir = os.path.join(outpath, year)
    mkdir(outdir)

    outdir = os.path.join(outdir, datestr)
    mkdir(outdir)

    # Output file name
    outfilename = os.path.join(outdir, fname)

    # exclude masked gates from the gridding
    my_gatefilter = pyart.filters.GateFilter(radar)
    my_gatefilter.exclude_transition()
    my_gatefilter.exclude_masked('reflectivity')

    grid_len = len(np.arange(-rmax, rmax + 1, xyres))
    grid_len_z = len(np.arange(0, maxheight + 1, zres))

    width = rmax * azimuthal_res * np.pi / 180.
    width = 1000 * np.round(2 * width / 1000) / 2  # rounding to the nearest 500 m.

    # Convert to linear Z for gridding.
    if linearz:
        radar.fields[refl_name]['data'] = 10 ** (radar.fields[refl_name]['data'] / 10.)

    # Gridding
    grid = pyart.map.grid_from_radars(
        radar, gatefilters=my_gatefilter, grid_shape=(grid_len_z, grid_len, grid_len),
        grid_limits=((0, maxheight), (-rmax, rmax), (-rmax, rmax)),
        gridding_algo="map_gates_to_grid", fields=fields_to_keep,
        weighting_function='CRESSMAN',
        map_roi=True, toa=maxheight, copy_field_data=True, algorithm='kd_tree',
        leafsize=10., roi_func='dist_beam', constant_roi=width,
        z_factor=0.05, xy_factor=0.02, min_radius=500.0,
        h_factor=1.0, nb=1.5, bsp=1.0, skip_transform=False)

    # Latitude Longitude field for each point.
    longitude, latitude = _get_latlon(grid)
    grid.add_field('longitude', longitude)
    grid.add_field('latitude', latitude)

    if linearz:
        # Convert back to dBZ
        radar.fields[refl_name]['data'] = 10 * np.log10(radar.fields[refl_name]['data'])
        grid.fields[refl_name]['data'] = 10 * np.log10(grid.fields[refl_name]['data'])

    try:
        grid.fields.pop('ROI')
    except Exception:
        pass

    # Saving data.
    grid.write(outfilename, arm_time_variables=True)

    return None
