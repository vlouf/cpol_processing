"""
Codes for correcting the differential phase and estimating KDP.

@title: phase
@author: Valentin Louf <valentin.louf@monash.edu>
@institutions: Monash University and the Australian Bureau of Meteorology
@date: 20/11/2017

.. autosummary::
    :toctree: generated/

    check_phidp
    fix_phidp_from_kdp
    phidp_bringi
    phidp_giangrande
    unfold_raw_phidp

TODO: Implement correction of PHIDP using region based as preprocessing for unfolding.
"""
# Python Standard Library
import copy

# Other Libraries
import pyart
import scipy
import netCDF4
import numpy as np

from scipy import integrate, ndimage
from csu_radartools import csu_kdp


def fix_phidp_from_kdp(phidp, kdp, gatefilter):
    """
    Correct PHIDP and KDP from spider webs.

    Parameters
    ==========
    radar:
        Py-ART radar data structure.
    gatefilter:
        Gate filter.
    kdp_name: str
        Differential phase key name.
    phidp_name: str
        Differential phase key name.

    Returns:
    ========
    phidp: ndarray
        Differential phase array.
    """
    kdp[gatefilter.gate_excluded] = 0
    kdp[(kdp < -4)] = 0
    kdp[kdp > 15] = 15
    interg = integrate.cumtrapz(kdp, radar.range['data'], axis=1)

    phidp[:, :-1] = interg / (len(radar.range['data']))
    return phidp, kdp


def phidp_bringi(radar, gatefilter, unfold_phidp_name="PHI_UNF", ncp_name="NCP",
                 rhohv_name="RHOHV_CORR", refl_field='DBZ'):
    """
    Compute PHIDP and KDP Bringi.

    Parameters
    ==========
    radar:
        Py-ART radar data structure.
    gatefilter:
        Gate filter.
    unfold_phidp_name: str
        Differential phase key name.
    refl_field: str
        Reflectivity key name.

    Returns:
    ========
    phidpb: ndarray
        Bringi differential phase array.
    kdpb: ndarray
        Bringi specific differential phase array.
    """
    dp = radar.fields[unfold_phidp_name]['data'].copy()
    dz = radar.fields[refl_field]['data'].copy().filled(-9999)

    try:
        if np.nanmean(dp[gatefilter.gate_included]) < 0:
            dp += 90
    except ValueError:
        pass

    # Extract dimensions
    rng = radar.range['data']
    azi = radar.azimuth['data']
    dgate = rng[1] - rng[0]
    [R, A] = np.meshgrid(rng, azi)

    # Compute KDP bringi.
    kdpb, phidpb, _ = csu_kdp.calc_kdp_bringi(dp, dz, R / 1e3, gs=dgate, bad=-9999, thsd=12, window=3.0, std_gate=11)

    # Mask array
    phidpb = np.ma.masked_where(phidpb == -9999, phidpb)
    kdpb = np.ma.masked_where(kdpb == -9999, kdpb)

    # Get metadata.
    phimeta = pyart.config.get_metadata("differential_phase")
    phimeta['data'] = phidpb
    kdpmeta = pyart.config.get_metadata("specific_differential_phase")
    kdpmeta['data'] = kdpb

    return phimeta, kdpmeta


def phidp_giangrande(radar, gatefilter, refl_field='DBZ', ncp_field='NCP',
                     rhv_field='RHOHV_CORR', phidp_field='PHIDP'):
    """
    Phase processing using the LP method in Py-ART. A LP solver is required,

    Parameters:
    ===========
    radar:
        Py-ART radar structure.
    gatefilter:
        Gate filter.
    refl_field: str
        Reflectivity field label.
    ncp_field: str
        Normalised coherent power field label.
    rhv_field: str
        Cross correlation ration field label.
    phidp_field: str
        Differential phase label.

    Returns:
    ========
    phidp_gg: dict
        Field dictionary containing processed differential phase shifts.
    kdp_gg: dict
        Field dictionary containing recalculated differential phases.
    """
    #  Preprocessing
    unfphi = pyart.correct.dealias_region_based(
        radar, gatefilter=gatefilter, vel_field=phidp_field, nyquist_vel=90)

    phi = radar.fields[phidp_field]['data']
    if phi.max() - phi.min() <= 200:  # 180 degrees plus some margin for noise...
        half_phi = True
    else:
        half_phi = False

    try:
        if np.nanmean(phi[gatefilter.gate_included]) < 0:
            unfphi['data'] += 90
    except ValueError:
        pass

    if half_phi:
        unfphi['data'] *= 2

    # unfphi['data'][gatefilter.gate_excluded] = np.NaN
    radar.fields[phidp_field]['data'] = unfphi['data']
    # Pyart version 1.10.
    phidp_gg, kdp_gg = pyart.correct.phase_proc_lp(radar, 0.0,
                                                   # gatefilter=gatefilter,
                                                   LP_solver='cylp',
                                                   ncp_field=ncp_field,
                                                   refl_field=refl_field,
                                                   phidp_field=phidp_field)

    # radar.fields.pop('PHITMP')
    if half_phi:
        unfphi['data'] /= 2

    try:
        radar.fields.pop('unfolded_differential_phase')
    except Exception:
        pass

    return phidp_gg, kdp_gg
