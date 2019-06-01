"""
Codes for estimating various parameters related to Hydrometeors.

@title: hydrometeors
@author: Valentin Louf <valentin.louf@monash.edu>
@institutions: Monash University and the Australian Bureau of Meteorology
@creation: 04/04/2017
@date: 20/11/2017

.. autosummary::
    :toctree: generated/

    dsd_retrieval
    hydrometeor_classification
    liquid_ice_mass
    merhala_class_convstrat
    rainfall_rate
"""
# Other Libraries
import pyart
import numpy as np

from csu_radartools import csu_liquid_ice_mass, csu_fhc, csu_blended_rain, csu_dsd


def dsd_retrieval(radar, gatefilter, kdp_name, zdr_name, refl_name='DBZ_CORR'):
    """
    Compute the DSD retrieval using the csu library.

    Parameters:
    ===========
        radar:
            Py-ART radar structure.
        refl_name: str
            Reflectivity field name.
        zdr_name: str
            ZDR field name.
        kdp_name: str
            KDP field name.

    Returns:
    ========
        nw_dict: dict
            Normalized Intercept Parameter.
        d0_dict: dict
            Median Volume Diameter.
    """
    dbz = radar.fields[refl_name]['data'].copy().filled(np.NaN)
    zdr = radar.fields[zdr_name]['data'].copy()
    try:
        kdp = radar.fields[kdp_name]['data'].copy().filled(np.NaN)
    except AttributeError:
        kdp = radar.fields[kdp_name]['data'].copy()

    d0, Nw, mu = csu_dsd.calc_dsd(dz=dbz, zdr=zdr, kdp=kdp, band='C')

    Nw = np.log10(Nw)
    Nw[gatefilter.gate_excluded] = np.NaN
    Nw = np.ma.masked_invalid(Nw).astype(np.float32)
    np.ma.set_fill_value(Nw, np.NaN)

    d0[gatefilter.gate_excluded] = np.NaN
    d0 = np.ma.masked_invalid(d0).astype(np.float32)
    np.ma.set_fill_value(d0, np.NaN)

    nw_dict = {'data': Nw,
               'units': 'AU', 'long_name': 'Normalized Intercept Parameter',
               'standard_name': 'Normalized Intercept Parameter',
               '_FillValue': np.NaN,
               '_Least_significant_digit': 2,
               'description': "Log10 of the NW. Retrieval based on Bringi et al. (2009)."}

    d0_dict = {'data': d0,
               'units': 'mm', 'long_name': 'Median Volume Diameter',
               'standard_name': 'Median Volume Diameter',
               '_FillValue': np.NaN,
               '_Least_significant_digit': 2,
               'description': "D0 retrieval based on Bringi et al. (2009)."}

    return nw_dict, d0_dict


def hydrometeor_classification(radar, gatefilter, kdp_name, zdr_name, refl_name='DBZ_CORR',
                               rhohv_name='RHOHV_CORR',
                               temperature_name='temperature',
                               height_name='height'):
    """
    Compute hydrometeo classification.

    Parameters:
    ===========
        radar:
            Py-ART radar structure.
        refl_name: str
            Reflectivity field name.
        zdr_name: str
            ZDR field name.
        kdp_name: str
            KDP field name.
        rhohv_name: str
            RHOHV field name.
        temperature_name: str
            Sounding temperature field name.
        height: str
            Gate height field name.

    Returns:
    ========
        hydro_meta: dict
            Hydrometeor classification.
    """
    refl = radar.fields[refl_name]['data'].copy().filled(np.NaN)
    zdr = radar.fields[zdr_name]['data'].copy().filled(np.NaN)
    try:
        kdp = radar.fields[kdp_name]['data'].copy().filled(np.NaN)
    except AttributeError:
        kdp = radar.fields[kdp_name]['data'].copy()
    rhohv = radar.fields[rhohv_name]['data']
    try:
        radar_T = radar.fields[temperature_name]['data']
        use_temperature = True
    except Exception:
        use_temperature = False

    if use_temperature:
        scores = csu_fhc.csu_fhc_summer(dz=refl, zdr=zdr, rho=rhohv, kdp=kdp, use_temp=True, band='C', T=radar_T)
    else:
        scores = csu_fhc.csu_fhc_summer(dz=refl, zdr=zdr, rho=rhohv, kdp=kdp, use_temp=False, band='C')

    hydro = np.argmax(scores, axis=0) + 1
    hydro[gatefilter.gate_excluded] = 0
    hydro_data = np.ma.masked_equal(hydro.astype(np.int16), 0)

    the_comments = "1: Drizzle; 2: Rain; 3: Ice Crystals; 4: Aggregates; " +\
                   "5: Wet Snow; 6: Vertical Ice; 7: LD Graupel; 8: HD Graupel; 9: Hail; 10: Big Drops"

    hydro_meta = {'data': hydro_data, 'units': ' ', 'long_name': 'Hydrometeor classification', '_FillValue': np.int16(0),
                  'standard_name': 'Hydrometeor_ID', 'comments': the_comments}

    return hydro_meta


def merhala_class_convstrat(radar, dbz_name="DBZ_CORR", rain_name="radar_estimated_rain_rate",
                            d0_name="D0", nw_name="NW"):
    """
    Merhala Thurai's has a criteria for classifying rain either Stratiform
    Convective or Mixed, based on the D-Zero value and the log10(Nw) value.
    Merhala's rain classification is 1 for Stratiform, 2 for Convective and 3
    for Mixed, 0 if no rain.

    """
    # Extracting data.
    d0 = radar.fields[d0_name]['data']
    nw = radar.fields[nw_name]['data']
    rainrate = radar.fields[rain_name]['data']
    dbz = radar.fields[dbz_name]['data']

    classification = np.zeros_like(dbz, dtype=int)

    # Invalid data
    pos0 = (d0 >= -5) & (d0 <= 100)
    pos1 = (nw >= -10) & (nw <= 100)

    # Classification index.
    indexa = nw - 6.4 + 1.7 * d0

    # Classifying
    classification[(indexa > 0.1) & (dbz > 20)] = 2
    classification[(indexa > 0.1) & (dbz <= 20)] = 1
    classification[indexa < -0.1] = 1
    classification[(indexa >= -0.1) & (indexa <= 0.1)] = 3

    # Masking invalid data.
    # classification = np.ma.masked_where(~pos0 | ~pos1 | dbz.mask, classification)

    # Generate metada.
    class_meta = {'data': classification,
                  'standard_name': 'echo_classification',
                  'long_name': 'Merhala Thurai echo classification',
                  'valid_min': 0,
                  'valid_max': 3,
                  'comment_1': 'Convective-stratiform echo classification based on Merhala Thurai',
                  'comment_2': '0 = Undefined, 1 = Stratiform, 2 = Convective, 3 = Mixed'}

    return class_meta


def rainfall_rate(radar, gatefilter, kdp_name, zdr_name, refl_name='DBZ_CORR',
                  hydro_name='radar_echo_classification', temperature_name='temperature'):
    """
    Rainfall rate algorithm from csu_radartools.

    Parameters:
    ===========
        radar:
            Py-ART radar structure.
        refl_name: str
            Reflectivity field name.
        zdr_name: str
            ZDR field name.
        kdp_name: str
            KDP field name.
        hydro_name: str
            Hydrometeor classification field name.

    Returns:
    ========
        rainrate: dict
            Rainfall rate.
    """
    dbz = radar.fields[refl_name]['data'].filled(np.NaN)
    zdr = radar.fields[zdr_name]['data'].filled(np.NaN)
    fhc = radar.fields[hydro_name]['data']
    try:
        kdp = radar.fields[kdp_name]['data'].filled(np.NaN)
    except AttributeError:
        kdp = radar.fields[kdp_name]['data']

    rain, _ = csu_blended_rain.calc_blended_rain_tropical(dz=dbz, zdr=zdr, kdp=kdp, fhc=fhc, band='C')

    rain[(gatefilter.gate_excluded) | np.isnan(rain) | (rain < 0)] = 0

    try:
        temp = radar.fields[temperature_name]['data']
        rain[temp < 0] = 0
    except Exception:
        pass

    rainrate = {"long_name": 'Blended Rainfall Rate',
                "units": "mm h-1",
                "standard_name": "rainfall_rate",
                '_Least_significant_digit': 2,
                '_FillValue': np.NaN,
                "description": "Rainfall rate algorithm based on Thompson et al. 2016.",
                "data": rain.astype(np.float32)}

    return rainrate
