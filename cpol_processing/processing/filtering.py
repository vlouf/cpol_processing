"""
Codes for creating and manipulating gate filters.

@title: filtering
@author: Valentin Louf <valentin.louf@monash.edu>
@institutions: Monash University and the Australian Bureau of Meteorology
@date: 20/11/2017

.. autosummary::
    :toctree: generated/

    do_gatefilter
    filter_hardcoding
    velocity_texture
"""
# Python Standard Library
import time
import datetime

# Other Libraries
import pyart
import scipy
import netCDF4
import numpy as np

from .phase import unfold_raw_phidp


def do_gatefilter(radar, refl_name='DBZ', phidp_name="PHIDP", rhohv_name='RHOHV_CORR',
                  zdr_name="ZDR", vel_field="VEL", tvel_name="TVEL", temp_name="temperature",
                  spectrum_name='WIDTH'):
    """
    Basic filtering

    Parameters:
    ===========
        radar:
            Py-ART radar structure.
        refl_name: str
            Reflectivity field name.
        rhohv_name: str
            Cross correlation ratio field name.
        ncp_name: str
            Name of the normalized_coherent_power field.
        zdr_name: str
            Name of the differential_reflectivity field.

    Returns:
    ========
        gf_despeckeld: GateFilter
            Gate filter (excluding all bad data).
    """
    radar_start_date = netCDF4.num2date(radar.time['data'][0], radar.time['units'].replace("since", "since "))

    gf = pyart.filters.GateFilter(radar)
    gf.include_above(rhohv_name, 0.9)
#     gf.exclude_below(rhohv_name, 0.75)
    try:
        gf.exclude_above(spectrum_name, 5)  # PI cf. Doviak and Zrnic book.
    except KeyError:
        pass

    zdr = radar.fields[zdr_name]['data']
    dbz = radar.fields[refl_name]['data']
    rho = radar.fields[rhohv_name]['data']
    temp = radar.fields[temp_name]['data']
    r = radar.range['data']
    azi = radar.azimuth['data']
    [R, A] = np.meshgrid(r, azi)

    emr4 = np.zeros_like(dbz) + 1
    # Positive temperature with low reflectivity and high zdr.
    emr4[(dbz < 0) & (temp >= 0)] = 0
    emr4[(rho < 0.75) & (dbz < 20)] = 0
    emr4[(R < 15e3) & (rho < 0.8)] = 0
    if radar_start_date.year > 2007:
        emr4[(R > 20e3) & (dbz > 20)] = 2

    radar.add_field_like(refl_name, "EMR", emr4, replace_existing=True)
    gf.exclude_equal('EMR', 0)
    gf.include_equal('EMR', 2)

    gf.exclude_outside(zdr_name, -3.0, 7.0)
    gf.exclude_outside(refl_name, -40.0, 80.0)

    try:
        gf.exclude_above(tvel_name, 3)
    except Exception:
        pass

    gf_despeckeld = pyart.correct.despeckle_field(radar, refl_name, gatefilter=gf)

    # Destroying temporary field.
    try:
        radar.fields.pop('EMR')
    except Exception:
        pass

    return gf_despeckeld


def filter_hardcoding(my_array, nuke_filter, bad=-9999):
    """
    Harcoding GateFilter into an array.

    Parameters:
    ===========
        my_array: array
            Array we want to clean out.
        nuke_filter: gatefilter
            Filter we want to apply to the data.
        bad: float
            Fill value.

    Returns:
    ========
        to_return: masked array
            Same as my_array but with all data corresponding to a gate filter
            excluded.
    """
    filt_array = np.ma.masked_where(nuke_filter.gate_excluded, my_array.copy())
    filt_array = filt_array.filled(fill_value=bad)
    to_return = np.ma.masked_where(filt_array == bad, filt_array)
    return to_return


def velocity_texture(radar, vel_name='VEL'):
    """
    Compute velocity texture using new Bobby Jackson function in Py-ART.

    Parameters:
    ===========
    radar:
        Py-ART radar structure.
    vel_name: str
        Name of the (original) Doppler velocity field.

    Returns:
    ========
    vdop_vel: dict
        Velocity texture.
    """

    try:
        v_nyq_vel = radar.instrument_parameters['nyquist_velocity']['data'][0]
    except Exception:
        vdop_art = radar.fields[vel_name]['data']
        v_nyq_vel = np.max(np.abs(vdop_art))

    vel_dict = pyart.retrieve.calculate_velocity_texture(radar, vel_name, nyq=v_nyq_vel)

    return vel_dict
