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


def texture(data):
    """Compute the texture of data.
    Compute the texture of the data by comparing values with a 3x3 neighborhood
    (based on :cite:`Gourley2007`). NaN values in the original array have
    NaN textures.
    Parameters
    ----------
    data : :class:`numpy:numpy.ndarray`
        multi-dimensional array with shape (..., number of beams, number
        of range bins)
    Returns
    ------
    texture : :class:`numpy:numpy.ndarray`
        array of textures with the same shape as data
    """
    x1 = np.roll(data, 1, -2)  # center:2
    x2 = np.roll(data, 1, -1)  # 4
    x3 = np.roll(data, -1, -2)  # 8
    x4 = np.roll(data, -1, -1)  # 6
    x5 = np.roll(x1, 1, -1)  # 1
    x6 = np.roll(x4, 1, -2)  # 3
    x7 = np.roll(x3, -1, -1)  # 9
    x8 = np.roll(x2, -1, -2)  # 7

    # at least one NaN would give a sum of NaN
    xa = np.array([x1, x2, x3, x4, x5, x6, x7, x8])

    # get count of valid neighboring pixels
    xa_valid = np.ones(np.shape(xa))
    xa_valid[np.isnan(xa)] = 0
    # count number of valid neighbors
    xa_valid_count = np.sum(xa_valid, axis=0)

    num = np.zeros(data.shape)
    for xarr in xa:
        diff = data - xarr
        # difference of NaNs will be converted to zero
        # (to not affect the summation)
        diff[np.isnan(diff)] = 0
        # only those with valid values are considered in the summation
        num += diff ** 2

    # reinforce that NaN values should have NaN textures
    num[np.isnan(data)] = np.nan

    return np.sqrt(num / xa_valid_count)


def do_gatefilter_cpol(radar, refl_name='DBZ', phidp_name="PHIDP", rhohv_name='RHOHV_CORR',
                       zdr_name="ZDR", snr_name='SNR'):
    """
    Filtering function adapted to CPOL.

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

    r = radar.range['data']
    azi = radar.azimuth['data']
    R, A = np.meshgrid(r, azi)
    refl = radar.fields[refl_name]['data'].copy()
    rho_corr = radar.fields[rhohv_name]['data']
    fcut = -0.6 / 140e3 * R + 0.8
    refl[rho_corr < fcut] = np.NaN
    radar.add_field_like(refl_name, 'NDBZ', refl)

    gf = pyart.filters.GateFilter(radar)
    gf.exclude_invalid('NDBZ')
    gf.exclude_below(snr_name, 9)

    gf.exclude_outside(zdr_name, -3.0, 7.0)
    gf.exclude_outside(refl_name, -20.0, 80.0)

    dphi = texture(radar.fields[phidp_name]['data'])
    radar.add_field_like(phidp_name, 'PHITXT', dphi)
    gf.exclude_above('PHITXT', 20)
    gf.exclude_below(rhohv_name, 0.45)

    # Remove rings in march 1999.
    if radar_start_date.year == 1999 and radar_start_date.month == 3:
        radar.add_field_like(refl_name, 'RRR', R)
        gf.exclude_above('RRR', 140e3)
        radar.fields.pop('RRR')

    gf_despeckeld = pyart.correct.despeckle_field(radar, refl_name, gatefilter=gf)

    # Remove tmp fields.
    try:
        radar.fields.pop('NDBZ')
        radar.fields.pop('PHITXT')
    except Exception:
        pass

    return gf_despeckeld


def do_gatefilter(radar, refl_name='DBZ', phidp_name="PHIDP", rhohv_name='RHOHV_CORR', zdr_name="ZDR", snr_name='SNR'):
    """
    Basic filtering function for dual-polarisation data.

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
    gf = pyart.filters.GateFilter(radar)
    # gf.exclude_below(snr_name, 9)

    gf.exclude_outside(zdr_name, -6.0, 7.0)
    gf.exclude_outside(refl_name, -20.0, 80.0)
    
    dphi = texture(radar.fields[phidp_name]['data'])
    radar.add_field_like(phidp_name, 'PHITXT', dphi)
    gf.exclude_above('PHITXT', 20)
    gf.exclude_below(rhohv_name, 0.6)
    
    gf_despeckeld = pyart.correct.despeckle_field(radar, refl_name, gatefilter=gf)
    
    try:
        radar.fields.pop('PHITXT')
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
