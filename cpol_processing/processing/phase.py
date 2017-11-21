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
"""
# Python Standard Library
import time
import datetime

# Other Libraries
import pyart
import scipy
import netCDF4
import numpy as np

from scipy import integrate
from csu_radartools import csu_kdp


def check_phidp(radar, phi_name="PHIDP"):
    """
    Check if PHIDP range is 180 degrees (half-circle) or 360 degrees.
    Parameters:
    ===========
    radar:
        Py-ART radar structure.
    phi_name: str
        Name of the PHIDP field.

    Return:
    =======
    half_phi: bool
        Is PHIDP range is half circle.
    """
    phi = radar.fields[phi_name]['data']
    if phi.max() - phi.min() <= 200:  # 180 degrees plus some margin for noise...
        half_phi = True
    else:
        half_phi = False

    return half_phi


def fix_phidp_from_kdp(radar, gatefilter, kdp_name="KDP_BRINGI", phidp_name="PHIDP_BRINGI"):
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
    kdp = radar.fields[kdp_name]['data'].copy()
    phidp = radar.fields[phidp_name]['data'].copy()
    kdp[gatefilter.gate_excluded] = 0
    kdp[(kdp > 15) | (kdp < -2)] = 0
    # kdp[kdp > 10] = 10
    interg = integrate.cumtrapz(kdp, radar.range['data'], axis=1)

    phidp[:, :-1] = interg / (len(radar.range['data']))
    return phidp


def phidp_bringi(radar, gatefilter, unfold_phidp_name="PHI_UNF", ncp_name="NCP", rhohv_name="RHOHV_CORR", refl_field='DBZ'):
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
    dp = pyart.correct.phase_proc.get_phidp_unf(radar, rhohv_lev=0.4, doc=None, refl_field=refl_field,
                                                ncp_field=ncp_name, rhv_field=rhohv_name, phidp_field=unfold_phidp_name)
    # Extract data
    try:
        dp = dp.filled(-9999)
    except Exception:
        pass
    dz = radar.fields[refl_field]['data'].copy()
    dz = np.ma.masked_where(gatefilter.gate_excluded, dz).filled(-9999)

    # Extract dimensions
    rng = radar.range['data']
    azi = radar.azimuth['data']
    dgate = rng[1] - rng[0]
    [R, A] = np.meshgrid(rng, azi)

    # Compute KDP bringi.
    kdpb, phidpb, _ = csu_kdp.calc_kdp_bringi(dp, dz, R / 1e3, gs=dgate, bad=-9999)

    # Mask array
    phidpb = np.ma.masked_where(phidpb == -9999, phidpb)
    kdpb = np.ma.masked_where(kdpb == -9999, kdpb)

    phi_off = np.min(np.min(phidpb, axis=1))
    phidpb -= phi_off

    # Get metadata.
    phimeta = pyart.config.get_metadata("differential_phase")
    phimeta['data'] = phidpb
    kdpmeta = pyart.config.get_metadata("specific_differential_phase")
    kdpmeta['data'] = kdpb

    return phimeta, kdpmeta


def phidp_giangrande(radar, gatefilter, refl_field='DBZ', ncp_field='NCP',
                     rhv_field='RHOHV_CORR', phidp_field='PHI_UNF'):
    """
    Phase processing using the LP method in Py-ART. A LP solver is required,

    Parameters:
    ===========
        radar:
            Py-ART radar structure.
        refl_field: str
            Reflectivity field name.
        ncp_field: str
            Normalised coherent power field name.
        rhv_field: str
            Cross correlation ration field name.
        phidp_field: str
            Phase field name.

    Returns:
    ========
        phidp_gg: dict
            Field dictionary containing processed differential phase shifts.
        kdp_gg: dict
            Field dictionary containing recalculated differential phases.
    """
    phidp_gg, kdp_gg = pyart.correct.phase_proc_lp(radar, 0.0,
                                                   min_rhv=0.95,
                                                   LP_solver='cylp',
                                                   refl_field=refl_field,
                                                   ncp_field=ncp_field,
                                                   rhv_field=rhv_field,
                                                   phidp_field=phidp_field)

    return phidp_gg, kdp_gg


def unfold_raw_phidp(radar, gatefilter, phi_name="PHIDP"):
    """
    Unfold raw PHIDP

    Parameters:
    ===========
    radar:
        Py-ART radar structure.
    gatefilter:
        Gate filter.
    phi_name: str
        Name of the PHIDP field.

    Returns:
    ========
    tru_phi: ndarray
        Unfolded raw PHIDP.
    """
    # Extract data
    phi = radar.fields[phi_name]['data'].copy()
    # For CPOL, PHIDP is properly unfolded before season 2003/2004
    CPOL_DATE_PHIDP_FOLD = datetime.datetime(2003, 10, 1)
    dtime = netCDF4.num2date(radar.time['data'][0], radar.time['units'])
    if dtime < CPOL_DATE_PHIDP_FOLD:
        tru_phi = phi
    else:
        phidp_unfold = np.ma.masked_where(gatefilter.gate_excluded, phi) + 180
        pmin = np.min(np.min(phidp_unfold, axis=1))
        tru_phi = phidp_unfold - pmin

    tru_phi += 45
    tru_phi[tru_phi > 360] -= 360

    return tru_phi
