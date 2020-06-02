"""
Codes for correcting and estimating attenuation on ZH and ZDR.

@title: attenuation
@author: Valentin Louf <valentin.louf@monash.edu>
@institutions: Monash University and the Australian Bureau of Meteorology
@date: 02/06/2020

.. autosummary::
    :toctree: generated/

    correct_attenuation_zdr
    correct_attenuation_zh_pyart
"""
# Other Libraries
import pyart
import numpy as np

from scipy.integrate import cumtrapz


def correct_attenuation_zdr(
    radar, gatefilter, zdr_name="ZDR_CORR", phidp_name="PHIDP_VAL", alpha=0.016
):
    """
    Correct attenuation on differential reflectivity. KDP_GG has been
    cleaned of noise, that's why we use it.

    V. N. Bringi, T. D. Keenan and V. Chandrasekar, "Correcting C-band radar
    reflectivity and differential reflectivity data for rain attenuation: a
    self-consistent method with constraints," in IEEE Transactions on Geoscience
    and Remote Sensing, vol. 39, no. 9, pp. 1906-1915, Sept. 2001.
    doi: 10.1109/36.951081

    Parameters:
    ===========
    radar:
        Py-ART radar structure.
    zdr_name: str
        Differential reflectivity field name.
    gatefilter:
        Filter excluding non meteorological echoes.
    kdp_name: str
        KDP field name.

    Returns:
    ========
    zdr_corr: array
        Attenuation corrected differential reflectivity.
    """
    zdr = radar.fields[zdr_name]["data"].copy()
    phi = radar.fields[phidp_name]["data"].copy()

    zdr_corr = zdr + 0.016 * phi
    zdr_corr[gatefilter.gate_excluded] = np.NaN
    zdr_corr = np.ma.masked_invalid(zdr_corr)
    np.ma.set_fill_value(zdr_corr, np.NaN)
    # Z-PHI coefficient from Bringi et al. 2001
    zdr_meta = pyart.config.get_metadata("differential_reflectivity")
    zdr_meta["description"] = "Attenuation corrected differential reflectivity using Bringi et al. 2001."
    zdr_meta["_FillValue"] = np.NaN
    zdr_meta["_Least_significant_digit"] = 2
    zdr_meta["data"] = zdr_corr

    return zdr_meta


def correct_attenuation_zh_pyart(
    radar,
    gatefilter,
    refl_field="DBZ",
    ncp_field="NCP",
    rhv_field="RHOHV_CORR",
    phidp_field="PHIDP_GG",
):
    """
    Correct attenuation on reflectivity using Py-ART tool. The attenuation from
    atmospheric gases is also corrected.

    Parameters:
    ===========
    radar:
        Py-ART radar structure.
    gatefilter:
        Filter excluding non meteorological echoes.
    refl_name: str
        Reflectivity field name.
    kdp_name: str
        KDP field name.

    Returns:
    ========
    atten_meta: dict
        Specific attenuation.
    zh_corr: array
        Attenuation corrected reflectivity.
    """
    # Compute attenuation
    spec_atten, _ = pyart.correct.calculate_attenuation(
        radar,
        0,
        rhv_min=0.3,
        refl_field=refl_field,
        ncp_field=rhv_field,
        rhv_field=rhv_field,
        phidp_field=phidp_field,
    )

    specific_atten = np.ma.masked_invalid(spec_atten["data"])
    r = radar.range["data"] / 1000
    dr = r[2] - r[1]

    na, nr = radar.fields[refl_field]["data"].shape
    attenuation = np.zeros((na, nr))
    attenuation[:, :-1] = 2 * cumtrapz(specific_atten, dx=dr)
    refl_corr = radar.fields[refl_field]["data"].copy() + attenuation
    refl_corr = np.ma.masked_where(gatefilter.gate_excluded, refl_corr)

    return refl_corr.astype(np.float32)
