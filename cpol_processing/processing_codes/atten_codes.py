"""
Codes for correcting and estimating attenuation on ZH and ZDR.

@title: atten_codes
@author: Valentin Louf <valentin.louf@monash.edu>
@institutions: Monash University and the Australian Bureau of Meteorology
@date: 15/08/2017

.. autosummary::
    :toctree: generated/

    correct_attenuation_zdr
    correct_attenuation_zh_pyart
"""

# Python Standard Library
import os
import copy
import datetime
from copy import deepcopy

# Other Libraries
import pyart
import numpy as np


def correct_attenuation_zdr(radar, zdr_name='ZDR_CORR', kdp_name='KDP_GG', alpha=0.016):
    """
    Correct attenuation on differential reflectivity. KDP_GG has been
    cleaned of noise, that's why we use it.

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
        atten_meta: dict
            Specific attenuation.
        zdr_corr: array
            Attenuation corrected differential reflectivity.
    """
    r = radar.range['data']
    zdr = deepcopy(radar.fields[zdr_name]['data'])
    kdp = deepcopy(radar.fields[kdp_name]['data'])

    dr = (r[1] - r[0]) / 1000  # km

    # Check if KDP is a masked array.
    if np.ma.isMaskedArray(kdp):
        kdp = kdp.filled(0)  # 0 is the neutral value for a sum
    else:
        kdp[np.isnan(kdp)] = 0

    atten_specific = alpha * kdp  # Bringi relationship
    atten_specific[np.isnan(atten_specific)] = 0
    # Path integrated attenuation
    atten = 2 * np.cumsum(atten_specific, axis=1) * dr

    zdr_corr = zdr + atten

    atten_meta = {'data': atten_specific, 'units': 'dB/km',
                  'standard_name': 'specific_attenuation_zdr',
                  'long_name': 'Differential reflectivity specific attenuation'}

    return atten_meta, zdr_corr


def correct_attenuation_zh(radar, dbz_name='DBZ', kdp_name='KDP_GG', alpha=0.08):
    """
    Correct attenuation on differential reflectivity. KDP_GG has been
    cleaned of noise, that's why we use it.

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
        atten_meta: dict
            Specific attenuation.
        zdr_corr: array
            Attenuation corrected differential reflectivity.
    """
    r = radar.range['data']
    zh = deepcopy(radar.fields[dbz_name]['data'])
    kdp = deepcopy(radar.fields[kdp_name]['data'])

    dr = (r[1] - r[0]) / 1000  # km

    # 0 is the neutral value for the sum
    # kdp = np.ma.masked_where((kdp < 0) | np.isnan(kdp), kdp).filled(0)

    atten_specific = alpha * kdp  # Bringi relationship
    atten_specific = np.ma.masked_where((atten_specific < 0) | (atten_specific > 1), atten_specific)
    # Path integrated attenuation
    atten = 2 * np.cumsum(atten_specific, axis=1) * dr

    atten_meta = {'data': atten_specific, 'units': 'dB/km',
                  'standard_name': 'specific_attenuation_reflectivity',
                  'long_name': 'Specific attenuation', 'valid_max': 1.0,
                  'valid_min': 0.0}

    zh_corr = {'data': zh + atten,
               'long_name': 'Corrected reflectivity',
               'standard_name': 'corrected_equivalent_reflectivity_factor',
               'units': 'dBZ'}

    return atten_meta, zh_corr


def correct_attenuation_zh_pyart(radar, refl_field='DBZ', ncp_field='NCP',
                                 rhv_field='RHOHV', phidp_field='PHIDP_GG'):
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
    atten_meta, zh_corr = pyart.correct.calculate_attenuation(radar, 0,
                                                              rhv_min=0.5,
                                                              refl_field=refl_field,
                                                              ncp_field=ncp_field,
                                                              rhv_field=rhv_field,
                                                              phidp_field=phidp_field)

    return atten_meta, zh_corr
