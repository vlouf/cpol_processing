"""
Codes for creating and manipulating gate filters. New functions: use of trained
Gaussian Mixture Models to remove noise and clutter from CPOL data before 2009.

@title: filtering.py
@author: Valentin Louf <valentin.louf@bom.gov.au>
@institutions: Monash University and the Australian Bureau of Meteorology
@created: 20/11/2017
@date: 25/02/2021

.. autosummary::
    :toctree: generated/

    texture
    get_clustering
    get_gatefilter_GMM
    do_gatefilter_cpol
    do_gatefilter
"""
# Libraries
import os
import gzip
import pickle

import pyart
import cftime
import numpy as np
import pandas as pd


def texture(data: np.ndarray) -> np.ndarray:
    """
    Compute the texture of data.
    Compute the texture of the data by comparing values with a 3x3 neighborhood
    (based on :cite:`Gourley2007`). NaN values in the original array have
    NaN textures. (Wradlib function)

    Parameters:
    ==========
    data : :class:`numpy:numpy.ndarray`
        multi-dimensional array with shape (..., number of beams, number
        of range bins)

    Returns:
    =======
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


def get_clustering(radar, vel_name: str = "VEL", phidp_name: str = "PHIDP", zdr_name: str = "ZDR"):
    """
    Create cluster using a trained Gaussian Mixture Model (I use scikit-learn)
    to cluster the radar data. Cluster 5 is clutter and 2 is noise. Cluster 1
    correponds to a high gradient on PHIDP (folding), so it may corresponds to
    either real data that fold or noise. A threshold on reflectivity should be
    used on cluster 1.

    Parameters:
    ===========
    radar:
        Py-ART radar structure.
    vel_name: str
        Velocity field name.
    phidp_name: str
        Name of the PHIDP field.
    zdr_name: str
        Name of the differential_reflectivity field.

    Returns:
    ========
    cluster: ndarray
        Data ID using GMM (5: clutter, 2: noise, and 1: high-phidp gradient).
    """
    # Load and deserialize GMM
    location = os.path.dirname(os.path.realpath(__file__))
    my_file = os.path.join(location, "data", "GM_model_CPOL.pkl.gz")
    with gzip.GzipFile(my_file, "r") as gzid:
        gmm = pickle.load(gzid)

    df_orig = pd.DataFrame(
        {
            "VEL": texture(radar.fields[vel_name]["data"]).flatten(),
            "PHIDP": texture(radar.fields[phidp_name]["data"]).flatten(),
            "ZDR": texture(radar.fields[zdr_name]["data"]).flatten(),
        }
    )

    df = df_orig.dropna()
    pos_droped = df_orig.dropna().index
    clusters = gmm.predict(df)

    r = radar.range["data"]
    time = radar.time["data"]
    R, _ = np.meshgrid(r, time)

    clus = np.zeros_like(R.flatten())
    clus[pos_droped] = clusters + 1
    cluster = clus.reshape(R.shape)

    return cluster


def get_gatefilter_GMM(
    radar, refl_name: str = "DBZ", vel_name: str = "VEL", phidp_name: str = "PHIDP", zdr_name: str = "ZDR"
):
    """
    Filtering function adapted to CPOL before 2009 using ML Gaussian Mixture
    Model. Function does 4 things:
    1) Cutoff of the reflectivities below the noise level.
    2) GMM using the texture of velocity, phidp and zdr.
    3) Filtering using 1) and 2) results.
    4) Removing temporary fields from the radar object.

    Parameters:
    ===========
    radar:
        Py-ART radar structure.
    refl_name: str
        Reflectivity field name.
    vel_name: str
        Velocity field name.
    phidp_name: str
        Name of the PHIDP field.
    zdr_name: str
        Name of the differential_reflectivity field.

    Returns:
    ========
    gf: GateFilter
        Gate filter (excluding all bad data).
    """
    # GMM clustering (indpdt from cutoff)
    cluster = get_clustering(radar, vel_name=vel_name, phidp_name=phidp_name, zdr_name=zdr_name)
    radar.add_field_like(refl_name, "CLUS", cluster, replace_existing=True)

    pos = (cluster == 1) & (radar.fields[refl_name]["data"] < 20)
    radar.add_field_like(refl_name, "TPOS", pos, replace_existing=True)

    # Using GMM results to filter.
    gf = pyart.filters.GateFilter(radar)
    gf.exclude_equal("CLUS", 5)
    gf.exclude_equal("CLUS", 2)
    gf.exclude_equal("TPOS", 1)
    gf = pyart.correct.despeckle_field(radar, refl_name, gatefilter=gf)

    # Removing temp keys.
    for k in ["TPOS", "CLUS"]:
        try:
            radar.fields.pop(k)
        except KeyError:
            continue

    return gf


def do_gatefilter_cpol(
    radar,
    refl_name: str = "DBZ",
    phidp_name: str = "PHIDP",
    rhohv_name: str = "RHOHV_CORR",
    zdr_name: str = "ZDR",
    snr_name: str = "SNR",
    vel_name: str = "VEL",
):
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
    radar_start_date = cftime.num2pydate(radar.time["data"][0], radar.time["units"])

    # if radar_start_date.year < 2009:
    gf = get_gatefilter_GMM(
        radar, refl_name=refl_name, vel_name=vel_name, phidp_name=phidp_name, zdr_name=zdr_name,
    )
    # else:
    #     gf = pyart.filters.GateFilter(radar)

    r = radar.range["data"]
    azi = radar.azimuth["data"]
    R, _ = np.meshgrid(r, azi)
    # refl = radar.fields[refl_name]["data"].copy()
    # fcut = 10 * np.log10(4e-5 * R)
    # refl[refl < fcut] = np.NaN
    # radar.add_field_like(refl_name, "NDBZ", refl)

    # gf.exclude_invalid("NDBZ")
    gf.exclude_below(snr_name, 9)
    gf.exclude_outside(zdr_name, -3.0, 7.0)
    gf.exclude_outside(refl_name, -20.0, 80.0)

    # if radar_start_date.year > 2007:
    #     gf.exclude_below(rhohv_name, 0.7)
    # else:
    rhohv = radar.fields[rhohv_name]["data"]
    pos = np.zeros_like(rhohv) + 1
    pos[(R < 90e3) & (rhohv < 0.7)] = 0
    radar.add_field_like(refl_name, "TMPRH", pos)
    gf.exclude_equal("TMPRH", 0)

    # Remove rings in march 1999.
    if radar_start_date.year == 1999 and radar_start_date.month == 3:
        radar.add_field_like(refl_name, "RRR", R)
        gf.exclude_above("RRR", 140e3)

    gf_despeckeld = pyart.correct.despeckle_field(radar, refl_name, gatefilter=gf)

    # Remove temporary fields.
    for k in ["NDBZ", "RRR", "TMPRH"]:
        try:
            radar.fields.pop(k)
        except KeyError:
            pass

    return gf_despeckeld


def do_gatefilter(
    radar,
    refl_name: str = "DBZ",
    phidp_name: str = "PHIDP",
    rhohv_name: str = "RHOHV_CORR",
    zdr_name: str = "ZDR",
    snr_name: str = "SNR",
):
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
    # Initialize gatefilter
    gf = pyart.filters.GateFilter(radar)

    # Remove obviously wrong data.
    gf.exclude_outside(zdr_name, -6.0, 7.0)
    gf.exclude_outside(refl_name, -20.0, 80.0)

    # Compute texture of PHIDP and remove noise.
    dphi = texture(radar.fields[phidp_name]["data"])
    radar.add_field_like(phidp_name, "PHITXT", dphi)
    gf.exclude_above("PHITXT", 20)
    gf.exclude_below(rhohv_name, 0.6)

    # Despeckle
    gf_despeckeld = pyart.correct.despeckle_field(radar, refl_name, gatefilter=gf)

    try:
        # Remove PHIDP texture
        radar.fields.pop("PHITXT")
    except Exception:
        pass

    return gf_despeckeld
