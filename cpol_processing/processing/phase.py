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

from numba import jit
from scipy import integrate, ndimage
from scipy.interpolate import interp1d
from csu_radartools import csu_kdp

from pyart.correct.phase_proc import smooth_and_trim_scan
from sklearn.linear_model import LinearRegression
from sklearn.isotonic import IsotonicRegression


def fix_phidp_from_kdp(phidp, kdp, r, gatefilter):
    """
    Correct PHIDP and KDP from spider webs.

    Parameters
    ==========
    r:
        Radar range.
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
    kdp[kdp > 15] = 0
    interg = integrate.cumtrapz(kdp, r, axis=1)

    phidp[:, :-1] = interg / (len(r))
    return phidp, kdp


def phidp_bringi(radar, gatefilter, unfold_phidp_name="PHI_UNF", refl_field='DBZ'):
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
    unfphidict = pyart.correct.dealias_region_based(
        radar, gatefilter=gatefilter, vel_field=phidp_field, nyquist_vel=90)

    phi = radar.fields[phidp_field]['data']
    if phi.max() - phi.min() <= 200:  # 180 degrees plus some margin for noise...
        half_phi = True
    else:
        half_phi = False

    vflag = np.zeros_like(phi)
    vflag[gatefilter.gate_excluded] = -3
    # unfphi, vflag = filter_data(phi, vflag, 90, 180, 40)
    unfphi = unfphidict['data']

    try:
        if np.nanmean(phi[gatefilter.gate_included]) < 0:
            unfphi += 90
    except ValueError:
        pass

    if half_phi:
        unfphi *= 2

    unfphi[vflag == -3] = 0

    # unfphi['data'][unfphi['data'] >= 340] = np.NaN
    radar.add_field_like(phidp_field, 'PHIDP_TMP', unfphi)
    # Pyart version 1.10.
    phidp_gg, kdp_gg = pyart.correct.phase_proc_lp(radar,
                                                   0.0,
                                                   # gatefilter=gatefilter,
                                                   LP_solver='cylp',
                                                   ncp_field=ncp_field,
                                                   refl_field=refl_field,
                                                   rhv_field=rhv_field,
                                                   phidp_field='PHIDP_TMP')

    radar.fields.pop('PHIDP_TMP')
    phidp_gg.pop('valid_min')

    if half_phi:
        unfphi['data'] /= 2
        phidp_gg['data'] /= 2
        kdp_gg['data'] /= 2

    phidp_gg['data'], kdp_gg['data'] = fix_phidp_from_kdp(phidp_gg['data'],
                                                          kdp_gg['data'],
                                                          radar.range['data'],
                                                          gatefilter)

    try:
        radar.fields.pop('unfolded_differential_phase')
    except Exception:
        pass

    return phidp_gg, kdp_gg


def _compute_kdp_from_phidp(r, phidp, window_len=35):
    """
    Compute KDP from PHIDP using Sobel filter. This is coming from pyart.

    Parameters:
    ===========
    r: ndarray
        Radar range.
    phidp: ndarray
        PhiDP field.
    window_len: int
        Size of the window for the Sobel filter.

    Returns:
    ========
    kdp_meta: dict
        KDP dictionary field.
    """
    sobel = 2. * np.arange(window_len) / (window_len - 1.0) - 1.0
    sobel = sobel / (abs(sobel).sum())
    sobel = sobel[::-1]
    gate_spacing = (r[1] - r[0]) / 1000.
    kdp = (scipy.ndimage.filters.convolve1d((phidp), sobel, axis=1) / ((window_len / 3) * 2 * gate_spacing))
    kdp[kdp > 12] = 12
    kdp[kdp < 0] = 0
    kdp[:, -window_len:] = 0

    return kdp


def valentin_phase_processing(radar, gatefilter, phidp_name='PHIDP', dbz_name='DBZ', bounds=[0, 360]):
    """
    Differential phase processing using machine learning technique.

    Parameters:
    ===========
    radar: struct
        Py-ART radar object structure.
    gatefilter: GateFilter
        Py-ART GateFilter object.
    phidp_name: str
        Name of the differential phase field.
    bounds: list
        Bounds, in degree, for PHIDP (0, 360).

    Returns:
    ========
        phitot: dict
            Processed differential phase.
    """
    # Check if PHIDP is in a 180 deg or 360 deg interval.
    nyquist = 90
    cutoff = 80

    scale_phi = True
    try:
        if np.nanmean(radar.fields[phidp_name]['data'][gatefilter.gate_included][:, :50]) > 0:
            scale_phi = False
    except Exception:
        pass

    # Dealiasing PHIDP using velocity dealiasing technique.
    unfphidict = pyart.correct.dealias_unwrap_phase(radar, gatefilter=gatefilter, skip_checks=True,
                                                    vel_field=phidp_name, nyquist_vel=90)
    # pyart.correct.dealias_region_based(radar, gatefilter=gatefilter, vel_field=phidp_name, nyquist_vel=nyquist)
    unfphi = unfphidict['data']
    if scale_phi:
        radar.fields[phidp_name]['data'] += 90
        unfphi += 90

    # Remove noise
    # unfphi[(unfphi < 0) | (radar.fields[phidp_name]['data'] > cutoff)] = np.NaN

    phitot = np.zeros_like(unfphi) + np.NaN
    unfphi[gatefilter.gate_excluded] = np.NaN
    nraymax, ngatemax = unfphi.shape
    x = radar.range['data'].copy()

    for ray in range(0, nraymax):
        # Taking the average the direct neighbours of each ray.
        y = unfphi[ray, :]
        y[x < 5e3] = 0  # Close to the radar is always extremly noisy

        y = elim_isolated(y)

        y = np.ma.masked_invalid(y)
        pos = ~y.mask

        x_nomask = x[pos].filled(np.NaN)
        y_nomask = y[pos].filled(np.NaN)

        if len(y_nomask[x_nomask > 5e3]) == 0:
            phitot[ray, :] = 0
            continue

        # y_nomask[x < 5e3] = 0
        # Machine learning stuff.
        ir = IsotonicRegression(-180, bounds[1])
        y_fit = ir.fit_transform(x_nomask, y_nomask)

        # y_fit = y_fit - y_fit.min()
        y_map = np.zeros((unfphi.shape[1])) + np.NaN
        y_map[pos] = y_fit

        phitot[ray, :] = np.convolve(populate_radials(fill_nan(y_map), ngatemax),
                                     np.ones(5) / 5)[:ngatemax]

    phi_unfold = pyart.config.get_metadata('differential_phase')
    phi_unfold['valid_min'] = 0
    phi_unfold['valid_max'] = 360
    phi_unfold['description'] = "Phase processing algorithm by Valentin Louf"

    phitot = phitot.astype(np.float32)
    phitot[gatefilter.gate_excluded] = np.NaN
    phi_unfold['data'] = np.ma.masked_invalid(phitot, 0)
    phi_unfold['_FillValue'] = np.NaN
    phi_unfold['_Least_significant_digit'] = 2

    # Computing KDP
    kdp = _compute_kdp_from_phidp(x, phitot)
    kdp = kdp.astype(np.float32)
    kdp[gatefilter.gate_excluded] = np.NaN
    kdp_meta = pyart.config.get_metadata('specific_differential_phase')
    kdp_meta['data'] = np.ma.masked_invalid(kdp)
    kdp_meta['_FillValue'] = np.NaN
    kdp_meta['_Least_significant_digit'] = 4
    kdp_meta['description'] = "Phase processing algorithm by Valentin Louf"

    return phi_unfold, kdp_meta


@jit(nopython=True)
def populate_radials(y_map, ngatemax):
    y_rslt = np.zeros((ngatemax))
    last_valid = 0
    for idx, value in enumerate(y_map):
        if np.isnan(value):
            y_rslt[idx] = last_valid
        else:
            y_rslt[idx] = value
            last_valid = value

    return y_rslt


@jit(nopython=True)
def elim_isolated(mydata):
    for idx in range(1, len(mydata) - 1):
        if np.isnan(mydata[idx - 1]) and np.isnan(mydata[idx + 1]):
            mydata[idx] = np.NaN
    return mydata


def fill_nan(A):
    '''
    interpolate to fill nan values
    '''
    inds = np.arange(A.shape[0])
    good = np.where(np.isfinite(A))
    f = interp1d(inds[good], A[good], bounds_error=False)
    B = np.where(np.isfinite(A), A, f(inds))
    return B