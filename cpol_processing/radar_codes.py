"""
Codes for correcting and estimating various radar and meteorological parameters.

@title: radar_codes
@author: Valentin Louf <valentin.louf@monash.edu>
@institutions: Monash University and the Australian Bureau of Meteorology
@creation: 04/04/2017
@date: 19/12/2018

.. autosummary::
    :toctree: generated/

    _my_snr_from_reflectivity
    _nearest
    check_azimuth
    check_reflectivity
    check_year
    correct_rhohv
    correct_zdr
    get_radiosoundings
    read_radar
    snr_and_sounding
"""
# Python Standard Library
import os
import re
import glob
import time
import fnmatch
import datetime

# Other Libraries
import pyart
import scipy
import netCDF4
import numpy as np


def _my_snr_from_reflectivity(radar, refl_field='DBZ'):
    """
    Just in case pyart.retrieve.calculate_snr_from_reflectivity, I can calculate
    it 'by hand'.
    Parameter:
    ===========
    radar:
        Py-ART radar structure.
    refl_field_name: str
        Name of the reflectivity field.

    Return:
    =======
    snr: dict
        Signal to noise ratio.
    """
    range_grid, _ = np.meshgrid(radar.range['data'], radar.azimuth['data'])
    range_grid += 1  # Cause of 0

    # remove range scale.. This is basically the radar constant scaled dBm
    pseudo_power = (radar.fields[refl_field]['data'] - 20.0 * np.log10(range_grid / 1000.0))
    # The noise_floor_estimate can fail sometimes in pyart, that's the reason
    # why this whole function exists.
    noise_floor_estimate = -40

    snr_field = pyart.config.get_field_name('signal_to_noise_ratio')
    snr_dict = pyart.config.get_metadata(snr_field)
    snr_dict['data'] = pseudo_power - noise_floor_estimate

    return snr_dict


def _nearest(items, pivot):
    """
    Find the nearest item.

    Parameters:
    ===========
        items:
            List of item.
        pivot:
            Item we're looking for.

    Returns:
    ========
        item:
            Value of the nearest item found.
    """
    return min(items, key=lambda x: abs(x - pivot))


def check_azimuth(radar, refl_field_name='DBZ'):
    """
    Checking if radar has a proper reflectivity field.  It's a minor problem
    concerning a few days in 2011 for CPOL.

    Parameters:
    ===========
        radar:
            Py-ART radar structure.
        refl_field_name: str
            Name of the reflectivity field.

    Return:
    =======
        True if radar has a proper azimuth field.
    """
    if radar.fields[refl_field_name]['data'].shape[0] < 360:
        return False

    return True


def check_reflectivity(radar, refl_field_name='DBZ'):
    """
    Checking if radar has a proper reflectivity field.  It's a minor problem
    concerning a few days in 2011 for CPOL.

    Parameters:
    ===========
        radar:
            Py-ART radar structure.
        refl_field_name: str
            Name of the reflectivity field.

    Return:
    =======
    True if radar has a non-empty reflectivity field.
    """
    dbz = radar.fields[refl_field_name]['data']

    if np.ma.isMaskedArray(dbz):
        if dbz.count() == 0:
            # Reflectivity field is empty.
            return False

    return True


def check_year(radar):
    """
    Check if time unit is correct. We have encountered some old files with the
    year 2000 bug, i.e. date being 2098 instead of 1998.

    Parameters:
    ===========
    radar:
        Py-ART radar structure.

    Returns:
    ========
        True if date seems valid and False if date century had to be corrected.
    """
    dtime = netCDF4.num2date(radar.time['data'][0], radar.time['units'])
    if dtime.year < 2050:
        # Date seems valid.
        return True
    else:
        wyr = dtime.year
        tunit = radar.time['units']
        radar.time['units'] = tunit.replace(str(wyr), str(wyr - 100))

    return False


def correct_azimuth(radar):
    """
    Check if the azimuth is right.

    Parameters:
    ===========
    radar: Py-ART radar structure

    Returns:
    ========
    azimuth: ndarray
        Corrected azimuth
    has_changed: bool
        Is there any change?
    """
    has_changed = False
    azimuth = radar.azimuth['data']
    for sl in range(radar.nsweeps):
        azi = azimuth[radar.get_slice(sl)]
        if np.sum(azi == 0) <= 2:
            continue

        azi_zero = azi[-1]
        for na in range(len(azi) - 2, -1, -1):
            if azi[na] != azi_zero - 1:
                if azi_zero == 0 and azi[na] == 359:
                    azi_zero = azi[na]
                    continue
                else:
                    has_changed = True
                    azi[na] = azi_zero - 1
            azi_zero = azi[na]

        azimuth[azimuth < 0] += 360
        azimuth[radar.get_slice(sl)] = azi

    return azimuth, has_changed


def correct_rhohv(radar, rhohv_name='RHOHV', snr_name='SNR'):
    """
    Correct cross correlation ratio (RHOHV) from noise. From the Schuur et al.
    2003 NOAA report (p7 eq 5)

    Parameters:
    ===========
        radar:
            Py-ART radar structure.
        rhohv_name: str
            Cross correlation field name.
        snr_name: str
            Signal to noise ratio field name.

    Returns:
    ========
        rho_corr: array
            Corrected cross correlation ratio.
    """
    rhohv = radar.fields[rhohv_name]['data'].copy()
    snr = radar.fields[snr_name]['data'].copy()

    natural_snr = 10**(0.1 * snr)
    natural_snr = natural_snr.filled(-9999)
    rho_corr = rhohv * (1 + 1 / natural_snr)

    # Not allowing the corrected RHOHV to be lower than the raw rhohv
    rho_corr[np.isnan(rho_corr) | (rho_corr < 0) | (rho_corr > 1)] = 1
    try:
        rho_corr = rho_corr.filled(1)
    except Exception:
        pass

    return rho_corr


def correct_standard_name(radar):
    """
    'standard_name' is a protected keyword for metadata in the CF conventions.
    To respect the CF conventions we can only use the standard_name field that
    exists in the CF table.

    Parameter:
    ==========
    radar: Radar object
        Py-ART data structure.
    """
    try:
        radar.range.pop('standard_name')
        radar.azimuth.pop('standard_name')
        radar.elevation.pop('standard_name')
    except Exception:
        pass

    try:
        radar.sweep_number.pop('standard_name')
        radar.fixed_angle.pop('standard_name')
        radar.sweep_mode.pop('standard_name')
    except Exception:
        pass

    good_keys = ['corrected_reflectivity', 'total_power', 'radar_estimated_rain_rate']
    for k in radar.fields.keys():
        if k not in good_keys:
            try:
                radar.fields[k].pop('standard_name')
            except Exception:
                continue

    radar.latitude['standard_name'] = 'latitude'
    radar.longitude['standard_name'] = 'longitude'
    radar.altitude['standard_name'] = 'altitude'

    return None


def correct_zdr(radar, zdr_name='ZDR', snr_name='SNR'):
    """
    Correct differential reflectivity (ZDR) from noise. From the Schuur et al.
    2003 NOAA report (p7 eq 6)

    Parameters:
    ===========
        radar:
            Py-ART radar structure.
        zdr_name: str
            Differential reflectivity field name.
        snr_name: str
            Signal to noise ratio field name.

    Returns:
    ========
        corr_zdr: array
            Corrected differential reflectivity.
    """
    zdr = radar.fields[zdr_name]['data'].copy()
    snr = radar.fields[snr_name]['data'].copy()
    alpha = 1.48
    natural_zdr = 10**(0.1 * zdr)
    natural_snr = 10**(0.1 * snr)
    corr_zdr = 10 * np.log10((alpha * natural_snr * natural_zdr) / (alpha * natural_snr + alpha - natural_zdr))

    return corr_zdr


def get_radiosoundings(sound_dir, radar_start_date):
    """
    Find the radiosoundings
    """
    def _fdate(flist):
        rslt = [None] * len(flist)
        for cnt, f in enumerate(flist):
            try:
                rslt[cnt] = datetime.datetime.strptime(re.findall("[0-9]{8}", f)[0], "%Y%m%d")
            except Exception:
                continue
        return rslt
    # Looking for radiosoundings:
    all_sonde_files = sorted(os.listdir(sound_dir))

    pos = [cnt for cnt, f in enumerate(all_sonde_files) if fnmatch.fnmatch(f, "*" + radar_start_date.strftime("%Y%m%d") + "*")]
    if len(pos) > 0:
        # Looking for the exact date.
        sonde_name = all_sonde_files[pos[0]]
        sonde_name = os.path.join(sound_dir, sonde_name)
    else:
        # Looking for the closest date.
        dtime_none = _fdate(all_sonde_files)
        dtime = [d for d in dtime_none if d is not None]
        closest_date = _nearest(dtime, radar_start_date)
        sonde_name = [e for e in all_sonde_files if closest_date.strftime("%Y%m%d") in e]
        if len(sonde_name) == 0:
            sonde_name = os.path.join(sound_dir, all_sonde_files[-1])
        elif type(sonde_name) is list:
            sonde_name = os.path.join(sound_dir, sonde_name[0])
        else:
            sonde_name = os.path.join(sound_dir, sonde_name)

    return sonde_name


def read_radar(radar_file_name):
    """
    Read the input radar file.

    Parameter:
    ==========
        radar_file_name: str
            Radar file name.

    Return:
    =======
        radar: struct
            Py-ART radar structure.
    """
    # Read the input radar file.
    try:
        if ".h5" in radar_file_name:
            radar = pyart.aux_io.read_odim_h5(radar_file_name, file_field_names=True)
        elif ".hdf" in radar_file_name:
            radar = pyart.aux_io.read_odim_h5(radar_file_name, file_field_names=True)
        else:
            radar = pyart.io.read(radar_file_name)
    except Exception:
        raise

    # SEAPOL hack change fields key.
    try:
        radar.fields['DBZ']
    except KeyError:
        myfields = [('SQIH', 'NCP'),
                    ('NCPH', "NCP"),
                    ('SNRH', 'SNR'),
                    ('normalized_coherent_power', "NCP"),
                    ('DBZH', "DBZ"),
                    ("DBZH_CLEAN", "DBZ"),
                    ('reflectivity', "DBZ"),
                    ('WRADH', "WIDTH"),
                    ('WIDTHH', "WIDTH"),
                    ('sprectrum_width', "WIDTH"),
                    ('UH', "DBZ"),
                    ('total_power', "DBZ"),
                    ("differential_reflectivity", "ZDR"),
                    ("VRADH", "VEL"),
                    ('VELH', "VEL"),
                    ('velocity', "VEL"),
                    ("cross_correlation_ratio", "RHOHV"),
                    ("differential_phase", "PHIDP"),
                    ("specific_differential_phase", "KDP")]
        for mykey, newkey in myfields:
            try:
                radar.add_field(newkey, radar.fields.pop(mykey))
            except Exception:
                continue

    return radar


def snr_and_sounding(radar, sonde_name, refl_field_name='DBZ', temp_field_name="temp"):
    """
    Compute the signal-to-noise ratio as well as interpolating the radiosounding
    temperature on to the radar grid. The function looks for the radiosoundings
    that happened at the closest time from the radar. There is no time
    difference limit.

    Parameters:
    ===========
        radar:
        sonde_name: str
            Path to the radiosoundings.
        refl_field_name: str
            Name of the reflectivity field.

    Returns:
    ========
        z_dict: dict
            Altitude in m, interpolated at each radar gates.
        temp_info_dict: dict
            Temperature in Celsius, interpolated at each radar gates.
        snr: dict
            Signal to noise ratio.
    """
    radar_start_date = netCDF4.num2date(radar.time['data'][0], radar.time['units'])
    # Altitude hack.
    true_alt = radar.altitude['data'].copy()
    radar.altitude['data'] = np.array([0])

    # print("Reading radiosounding %s" % (sonde_name))
    interp_sonde = netCDF4.Dataset(sonde_name)
    temperatures = interp_sonde.variables[temp_field_name][:]
    temperatures[(temperatures < -100) | (temperatures > 100)] = np.NaN
    try:
        temperatures = temperatures.filled(np.NaN)
    except AttributeError:
        pass
    # times = interp_sonde.variables['time'][:]
    # heights = interp_sonde.variables['height'][:]

    # Height profile corresponding to radar.
    my_profile = pyart.retrieve.fetch_radar_time_profile(interp_sonde, radar)

    # CPOL altitude is 50 m.
    good_altitude = my_profile['height'] >= 0
    # Getting the temperature
    z_dict, temp_dict = pyart.retrieve.map_profile_to_gates(temperatures[good_altitude],
                                                            my_profile['height'][good_altitude],
                                                            radar)

    temp_info_dict = {'data': temp_dict['data'],
                      'long_name': 'Sounding temperature at gate',
                      'standard_name': 'temperature',
                      'valid_min': -100, 'valid_max': 100,
                      'units': 'degrees Celsius',
                      'comment': 'Radiosounding date: %s' % (radar_start_date.strftime("%Y/%m/%d"))}

    # Altitude hack
    radar.altitude['data'] = true_alt

    # Calculate SNR
    snr = pyart.retrieve.calculate_snr_from_reflectivity(radar, refl_field=refl_field_name)
    # Sometimes the SNR is an empty array, this is due to the toa parameter.
    # Here we try to recalculate the SNR with a lower value for toa (top of atm).
    if snr['data'].count() == 0:
        snr = pyart.retrieve.calculate_snr_from_reflectivity(radar, refl_field=refl_field_name, toa=20000)

    if snr['data'].count() == 0:
        # If it fails again, then we compute the SNR with the noise value
        # given by the CPOL radar manufacturer.
        snr = _my_snr_from_reflectivity(radar, refl_field=refl_field_name)

    return z_dict, temp_info_dict, snr
