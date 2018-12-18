"""
Codes for correcting and estimating attenuation on ZH and ZDR.

@title: attenuation
@author: Valentin Louf <valentin.louf@monash.edu>
@institutions: Monash University and the Australian Bureau of Meteorology
@date: 20/11/2017

.. autosummary::
    :toctree: generated/

    correct_attenuation_zdr
    correct_attenuation_zh
    correct_attenuation_zh_pyart
"""
# Other Libraries
import pyart
import numpy as np


def correct_attenuation_zdr(radar, zdr_name='ZDR_CORR', phidp_name='PHIDP_VAL',
                            kdp_name='KDP_GG', alpha=0.016):
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
        kdp_name: str
            KDP field name.

    Returns:
    ========
        zdr_corr: array
            Attenuation corrected differential reflectivity.
    """
    r = radar.range['data']
    zdr = radar.fields[zdr_name]['data']
    phi = radar.fields[phidp_name]['data']

    # Z-PHI coefficient from Bringi et al. 2001
    zdr_corr = zdr + 0.016 * phi

    return zdr_corr.copy()


def correct_attenuation_zh_pyart(radar, refl_field='DBZ', ncp_field='NCP',
                                 rhv_field='RHOHV_CORR', phidp_field='PHIDP_GG'):
    """
    Correct attenuation on reflectivity using Py-ART tool.

    Parameters:
    ===========
    radar:
        Py-ART radar structure.
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
    atten_meta, zh_corr = pyart.correct.calculate_attenuation(radar, 0,
                                                              rhv_min=0.3,
                                                              refl_field=refl_field,
                                                              ncp_field=rhv_field,
                                                              rhv_field=rhv_field,
                                                              phidp_field=phidp_field)

    # Correct DBZ from attenuation manually.
    # dbz = radar.fields[refl_field]['data'].copy()
    # att = atten_meta['data']
    # r = radar.range['data']
    # dr = r[1] - r[0]

    # int_att = np.cumsum(att, axis=1) * dr / 1e3
    # dbz_corr = dbz + int_att
    # zh_corr['data'] = dbz_corr

    return atten_meta, zh_corr
