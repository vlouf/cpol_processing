"""
Codes for creating and manipulating gate filters.

@title: filtering
@author: Valentin Louf <valentin.louf@monash.edu>
@institutions: Monash University and the Australian Bureau of Meteorology
@date: 20/11/2017

.. autosummary::
    :toctree: generated/

    _mask_rhohv
    _noise_th
    do_gatefilter
    do_txt_gatefilter
    filter_hardcoding
    phidp_texture
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


def _mask_rhohv(radar, rhohv_name, tight=True):
    nrays = radar.nrays
    ngate = radar.ngates
    oneray = np.zeros((ngate))
    oneray[:(ngate // 2)] = 1 - np.linspace(0.05, 0.4, ngate // 2)
    oneray[(ngate // 2):] = 0.3
    emr = np.vstack([oneray for e in range(nrays)])
    rho = radar.fields[rhohv_name]['data']
    emr2 = np.zeros(rho.shape)
    emr2[rho > emr] = 1
    return emr2


def _noise_th(x, max_range=90):
    n, bins = np.histogram(x.flatten(), bins=150, range=(5, max_range))
    cnt = 10
    peaks = []
    while (len(peaks) < 1) or (cnt == 0):
        peaks = scipy.signal.find_peaks_cwt(n, [cnt])
        cnt - 1

    centers = bins[0:-1] + (bins[1] - bins[0])
    search_data = n[peaks[0]:peaks[1]]
    search_centers = centers[peaks[0]:peaks[1]]
    locs = search_data.argsort()
    noise_threshold = search_centers[locs[0]]

    return noise_threshold


def do_gatefilter(radar, refl_name='DBZ', phidp_name="PHIDP", rhohv_name='RHOHV_CORR', zdr_name="ZDR", vel_field="VEL", temp_name="temperature"):
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
    gf.exclude_outside(zdr_name, -3.0, 7.0)
    gf.exclude_outside(refl_name, -40.0, 80.0)

    if radar_start_date.year <= 2007:
        # Using velocity texture for filtering.
        vnyq = radar.fields[vel_field]['data'].max()
        tvel = pyart.retrieve.calculate_velocity_texture(radar, vel_field=vel_field, nyq=vnyq)
        radar.add_field("TVEL", tvel, replace_existing=True)
        gf.exclude_above("TVEL", 4)
        gf.include_above(rhohv_name, 0.8)
    else:
        # Using PHIDP texture for filtering.
        tvel = pyart.retrieve.calculate_velocity_texture(radar, vel_field=phidp_name, nyq=90)
        radar.add_field("TVEL", tvel, replace_existing=True)
        gf.exclude_above("TVEL", 30)
        gf.exclude_below(rhohv_name, 0.5)

    zdr = radar.fields[zdr_name]['data']
    dbz = radar.fields[refl_name]['data']
    rho = radar.fields[rhohv_name]['data']
    temp = radar.fields[temp_name]['data']
    r = radar.range['data']
    azi = radar.azimuth['data']
    [R, A] = np.meshgrid(r, azi)

    emr4 = np.zeros_like(dbz) + 1
    emr4[(zdr > 4) & (dbz < 20) & (temp >= 0)] = 0
    emr4[(rho < 0.75) & (dbz < 20)] = 0
    emr4[(R >50e3) & (dbz > 20) & (rho > 0.4)] = 1

    radar.add_field_like(refl_name, "EMR", emr4, replace_existing=True)
    gf.exclude_equal('EMR', 0)

    try:
        radar.fields.pop('TVEL')
        radar.fields.pop('EMR')
    except Exception:
        pass

    gf_despeckeld = pyart.correct.despeckle_field(radar, refl_name, gatefilter=gf)

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


def phidp_texture(radar, phidp_name='PHIDP'):
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

    v_nyq_vel = radar.fields[phidp_name]['data'].max()

    tphi_dict = pyart.retrieve.calculate_velocity_texture(radar, phidp_name, nyq=v_nyq_vel, check_nyq_uniform=False)
    tphi_dict['long_name'] = "Differential phase texture"
    tphi_dict['standard_name'] = "texture_of_differential_phase"
    tphi_dict['units'] = "deg"

    return tphi_dict


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
    except (KeyError, IndexError):
        vdop_art = radar.fields[vel_name]['data']
        v_nyq_vel = np.max(np.abs(vdop_art))

    vel_dict = pyart.retrieve.calculate_velocity_texture(radar, vel_name, nyq=v_nyq_vel)

    return vel_dict
