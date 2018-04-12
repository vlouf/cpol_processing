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

# from numba import jit


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


# @jit(nopython=True)
# def _get_iter_pos(azi, st, nb=180):
#     """
#     Return a sequence of integers from start (inclusive) to stop (start + nb)
#     by step of 1 for iterating over the azimuth (handle the case that azimuth
#     360 is in fact 0, i.e. a cycle).
#     JIT-friendly function (this is why the function looks much longer than the
#     pythonic way of doing this).
#
#     Parameters:
#     ===========
#     azi: ndarray<float>
#         Azimuth.
#     st: int
#         Starting point.
#     nb: int
#         Number of unitary steps.
#
#     Returns:
#     ========
#     out: ndarray<int>
#         Array containing the position from start to start + nb, i.e.
#         azi[out[0]] <=> st
#     """
#     if st < 0:
#         st += len(azi)
#     if st >= len(azi):
#         st -= len(azi)
#
#     ed = st + nb
#     if ed >= len(azi):
#         ed -= len(azi)
#     if ed < 0:
#         ed += len(azi)
#
#     posazi = np.arange(0, len(azi))
#     mypos = np.empty_like(posazi)
#
#     if nb > 0:
#         if st < ed:
#             end = ed - st
#             mypos[:end] = posazi[st:ed]
#         else:
#             mid = (len(azi) - st)
#             end = (len(azi) - st + ed)
#             mypos[:mid] = posazi[st:]
#             mypos[mid:end] = posazi[:ed]
#     else:  # Goin backward.
#         if st < ed:
#             mid = st + 1
#             end = st + len(azi) - ed
#             mypos[:mid] = posazi[st::-1]
#             mypos[mid:end] = posazi[-1:ed:-1]
#         else:
#             end = np.abs(st - ed)
#             mypos[:end] = posazi[st:ed:-1]
#
#     out = np.zeros((end, ), dtype=mypos.dtype)
#     for n in range(end):
#         out[n] = mypos[n]
#
#     return out
#
#
# @jit(nopython=True)
# def _get_iter_range(pos_center, nb_gate, maxrange):
#     """
#     Similar as get_iter_pos, but this time for creating an array of iterative
#     indices over the radar range. JIT-friendly function.
#
#     Parameters:
#     ===========
#     pos_center: int
#         Starting point
#     nb_gate: int
#         Number of gates to iter to.
#     maxrange: int
#         Length of the radar range, i.e. maxrange = len(r)
#
#     Returns:
#     ========
#     Array of iteration indices.
#     """
#     half_range = nb_gate // 2
#     if pos_center < half_range:
#         st_pos = 0
#     else:
#         st_pos = pos_center - half_range
#
#     if pos_center + half_range >= maxrange:
#         end_pos = maxrange
#     else:
#         end_pos = pos_center + half_range
#
#     return np.arange(st_pos, end_pos)
#
#
# @jit(nopython=True)
# def _comp_refl_sl(ground_range_reference, ground_range_slice, azimuth_reference, azimuth_slice, dbz_ref, dbz_slice):
#     """
#     Dealias using 3D continuity. This function will look at the velocities from
#     one sweep (the reference) to the other (the slice).
#
#     Parameters:
#     ===========
#     ground_range_reference: ndarray
#         Radar range
#     ground_range_slice: ndarray
#         Radar range
#     azimuth_reference: ndarray
#         Azimuth of the reference sweep.
#     azimuth_slice: ndarray
#         Azimuth of the sweep to dealias.
#     dbz_ref: ndarray <azimuth, r>
#         Reference reflectivity
#     dbz_slice: ndarray <azimuth, r>
#         Reflectivity to be corrected.
#
#     Returns:
#     ========
#     dbz_slice: ndarray <azimuth, r>
#         Reflectivity corrected.
#     """
#     maxazi, maxrange = dbz_slice.shape
#     window_azimuth = 3
#     window_range = 3
#
#     for nazi in range(maxazi):
#         for ngate in range(maxrange):
#             if np.isnan(dbz_slice[nazi, ngate]):
#                 continue
#
#             current_dbz = dbz_slice[nazi, ngate]
#
#             rpos_reference = np.argmin(np.abs(ground_range_reference - ground_range_slice[ngate]))
#             apos_reference = np.argmin(np.abs(azimuth_reference - azimuth_slice[nazi]))
#
#             apos_iter = _get_iter_pos(azimuth_reference, apos_reference - window_azimuth // 2,
#                                       window_azimuth)
#             rpos_iter = _get_iter_range(rpos_reference, window_range, maxrange)
#
#             comp_dbz = np.zeros((len(rpos_iter) * len(apos_iter))) + np.NaN
#
#             cnt = -1
#             for na in apos_iter:
#                 for nr in rpos_iter:
#                     cnt += 1
#                     comp_dbz[cnt] = dbz_ref[na, nr]
#
#             if np.sum(~np.isnan(comp_dbz)) == 0:
#                 dbz_slice[nazi, ngate] = np.NaN
#
#     return dbz_slice
#
#
# def compare_refl_3D(radar, gatefilter, dbz_name='DBZ'):
#     """
#     Compare the reflectivity of the first 4 sweeps and check that 3D continuity
#     exists. The goals is to remove the recalcitrant permanent clutter.
#
#     Parameters:
#     ===========
#     radar:
#         Py-ART radar structure.
#     gatefilter:
#         Gate filter.
#     dbz_name: str
#         Name of the reflectivity.
#
#     Returns:
#     ========
#     dbz_tot: Corrected reflectivity.
#     """
#     dbz_tot = radar.fields[dbz_name]['data'].copy()
#     dbz_tot = np.ma.masked_where(gatefilter.gate_excluded, dbz_tot).filled(np.NaN)
#
#     r_ref = radar.range['data']
#     azi_ref = radar.azimuth['data'][radar.get_slice(1)]
#     dbz_ref = dbz_tot[radar.get_slice(1)]
#
#     azi_sl = radar.azimuth['data'][radar.get_slice(0)]
#     dbz_sl = dbz_tot[radar.get_slice(0)]
#     dbz_ref = _comp_refl_sl(r_ref, r_ref, azi_ref, azi_sl, dbz_ref, dbz_sl)
#     azi_ref = azi_sl
#
#     dbz_tot[radar.get_slice(0)] = dbz_ref
#
#     r_ref = radar.range['data']
#     azi_ref = radar.azimuth['data'][radar.get_slice(0)]
#     dbz_ref = dbz_tot[radar.get_slice(0)]
#
#     for sl in range(1, 3):
#         azi_sl = radar.azimuth['data'][radar.get_slice(sl)]
#         dbz_sl = dbz_tot[radar.get_slice(sl)]
#         dbz_ref = _comp_refl_sl(r_ref, r_ref, azi_ref, azi_sl, dbz_ref, dbz_sl)
#         azi_ref = azi_sl
#
#         dbz_tot[radar.get_slice(sl)] = dbz_ref
#
#     dbz_tot = np.ma.masked_where(gatefilter.gate_excluded, dbz_tot)
#
#     return dbz_tot


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
