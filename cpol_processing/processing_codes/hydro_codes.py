"""
Codes for estimating various parameters related to Hydrometeors.

@title: radar_codes
@author: Valentin Louf <valentin.louf@monash.edu>
@institutions: Monash University and the Australian Bureau of Meteorology
@creation: 04/04/2017
@date: 15/05/2017

.. autosummary::
    :toctree: generated/

    dsd_retrieval
    hydrometeor_classification
    liquid_ice_mass
    rainfall_rate
"""
# Other Libraries
import pyart
import numpy as np

from csu_radartools import csu_liquid_ice_mass, csu_fhc, csu_blended_rain, csu_dsd


def dsd_retrieval(radar, refl_name='DBZ', zdr_name='ZDR', kdp_name='KDP_GG'):
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
    dbz = radar.fields[refl_name]['data'].filled(np.NaN)
    zdr = radar.fields[zdr_name]['data']
    try:
        kdp = radar.fields[kdp_name]['data'].filled(np.NaN)
    except AttributeError:
        kdp = radar.fields[kdp_name]['data']

    d0, Nw, mu = csu_dsd.calc_dsd(dz=dbz, zdr=zdr, kdp=kdp, band='C')

    Nw = np.log10(Nw)
    Nw = np.ma.masked_where(np.isnan(Nw), Nw)
    d0 = np.ma.masked_where(np.isnan(d0), d0)

    nw_dict = {'data': Nw,
               'units': 'unitless (log10)', 'long_name': 'Log10 of the Normalized Intercept Parameter',
               'standard_name': 'Log10 of the Normalized Intercept Parameter',
               'description': "NW retrieval based on Bringi et al. (2009). Mu can not be retrieved alongside NW and D0."}

    d0_dict = {'data': d0,
               'units': 'mm', 'long_name': 'Median Volume Diameter',
               'standard_name': 'Median Volume Diameter',
               'description': "D0 retrieval based on Bringi et al. (2009). Mu can not be retrieved alongside NW and D0."}

    return nw_dict, d0_dict


def hydrometeor_classification(radar, refl_name='DBZ_CORR', zdr_name='ZDR_CORR',
                               kdp_name='KDP_GG', rhohv_name='RHOHV_CORR',
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
    refl = radar.fields[refl_name]['data']
    zdr = radar.fields[zdr_name]['data']
    kdp = radar.fields[kdp_name]['data']
    rhohv = radar.fields[rhohv_name]['data']
    radar_T = radar.fields[temperature_name]['data']
    radar_z = radar.fields[height_name]['data']

    scores = csu_fhc.csu_fhc_summer(dz=refl, zdr=zdr, rho=rhohv, kdp=kdp, use_temp=True, band='C', T=radar_T)

    hydro = np.argmax(scores, axis=0) + 1
    fill_value = -32768
    hydro_data = np.ma.masked_where(hydro == fill_value, hydro)

    the_comments = "1: Drizzle; 2: Rain; 3: Ice Crystals; 4: Aggregates; " +\
                   "5: Wet Snow; 6: Vertical Ice; 7: LD Graupel; 8: HD Graupel; 9: Hail; 10: Big Drops"

    hydro_meta = {'data': hydro_data, 'units': ' ', 'long_name': 'Hydrometeor classification',
                  'standard_name': 'Hydrometeor_ID', 'comments': the_comments}

    return hydro_meta


def liquid_ice_mass(radar, refl_name='DBZ_CORR', zdr_name='ZDR_CORR',
                    temperature_name='temperature', height_name='height'):
    """
    Compute the liquid/ice water content using the csu library.

    Parameters:
    ===========
        radar:
            Py-ART radar structure.
        refl_name: str
            Reflectivity field name.
        zdr_name: str
            ZDR field name.
        temperature_name: str
            Sounding temperature field name.
        height: str
            Gate height field name.

    Returns:
    ========
        liquid_water_mass: dict
            Liquid water content.
        ice_mass: dict
            Ice water content.
    """
    refl = radar.fields[refl_name]['data']
    zdr = radar.fields[zdr_name]['data']
    radar_T = radar.fields[temperature_name]['data']
    radar_z = radar.fields[height_name]['data']

    liquid_water_mass, ice_mass = csu_liquid_ice_mass.calc_liquid_ice_mass(refl, zdr, radar_z / 1000.0, T=radar_T)

    liquid_water_mass = {'data': liquid_water_mass, 'units': 'g m-3',
                         'long_name': 'Liquid Water Content',
                         'standard_name': 'liquid_water_content',
                         'description': "Liquid Water Content using Carey and Rutledge (2000) algorithm."}
    ice_mass = {'data': ice_mass, 'units': 'g m-3', 'long_name': 'Ice Water Content',
                'standard_name': 'ice_water_content',
                'description': "Ice Water Content using Carey and Rutledge (2000) algorithm."}

    return liquid_water_mass, ice_mass


def rainfall_rate(radar, refl_name='DBZ_CORR', zdr_name='ZDR_CORR', kdp_name='KDP_GG', hydro_name='radar_echo_classification'):
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

    rain, method = csu_blended_rain.calc_blended_rain_tropical(dz=dbz, zdr=zdr, kdp=kdp, fhc=fhc, band='C')
    rain[rain == 0] = np.NaN
    rain = np.ma.masked_where(np.isnan(rain), rain)

    rainrate = {"long_name": 'Blended Rainfall Rate',
                "units": "mm h-1",
                "standard_name": "Rainfall Rate",
                "description": "Rainfall rate algorithm based on Thompson et al. 2016.",
                "data": rain}

    return rainrate
