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


def do_gatefilter(radar, refl_name='DBZ', phidp_name="PHIDP", rhohv_name='RHOHV_CORR', zdr_name="ZDR"):
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
    # For CPOL, there is sometime an issue with older seasons.
    dbz = radar.fields[refl_name]['data']
    rhohv = radar.fields[rhohv_name]['data']
    r = radar.range['data']
    azi = radar.azimuth['data']
    [R, A] = np.meshgrid(r, azi)
    radar_start_date = netCDF4.num2date(radar.time['data'][0], radar.time['units'].replace("since", "since "))

    gf = pyart.filters.GateFilter(radar)

    gf.exclude_outside(zdr_name, -3.0, 7.0)
    gf.exclude_outside(refl_name, -40.0, 80.0)
    # gf.exclude_below(rhohv_name, 0.6)

    if radar_start_date.year > 2004:
        phi = unfold_raw_phidp(radar, gf)
    else:
        phi = radar.fields[phidp_name]['data'].copy()

    vel_dict = pyart.util.angular_texture_2d(phi, 4, phi.max())
    try:
        noise_threshold = _noise_th(vel_dict)
    except Exception:
        print("Only noise in volume")
        gf.exclude_below(rhohv_name, 0.8)
        noise_threshold = None

    if noise_threshold is not None:
        radar.add_field_like(refl_name, "TPHI", vel_dict, replace_existing=True)
        gf.exclude_above("TPHI", noise_threshold * 1.15)
        gf.include_above(rhohv_name, 0.7)
        radar.fields.pop('TPHI')

    if radar_start_date.year < 2006:
        posi, posj = np.where((dbz < 20) & (rhohv < 0.8))
        nar = np.zeros_like(dbz) + 1
        nar[posi, posj] = 0
        radar.add_field_like("NCP", "EMR3", nar, replace_existing=True)
        gf.exclude_equal("EMR3", 0)
        radar.fields.pop("EMR3")

    if radar_start_date.year == 2007:
        tvel = pyart.retrieve.calculate_velocity_texture(radar, vel_field="VEL")
        try:
            noise_threshold = _noise_th(tvel['data'])
        except Exception:
            noise_threshold = None
            pass

        if noise_threshold is not None:
            radar.add_field("TVEL", tvel)
            gf.eclude_above("TVEL", noise_threshold * 1.15)
            radar.fields.pop('TVEL')

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
