"""
Codes for correcting and estimating various radar and meteorological parameters.

@title: radar_codes
@author: Valentin Louf <valentin.louf@monash.edu>
@institutions: Monash University and the Australian Bureau of Meteorology
@creation: 04/04/2017
@date: 14/11/2017

.. autosummary::
    :toctree: generated/

    _my_snr_from_reflectivity
    _nearest
    _get_noise_threshold
    check_azimuth
    check_reflectivity
    correct_rhohv
    correct_zdr
    do_gatefilter
    get_texture
    filter_hardcoding
    get_field_names
    read_radar
    refold_velocity
    rename_radar_fields
    snr_and_sounding
    unfold_velocity
"""
# Python Standard Library
import os
import re
import glob
import time
import fnmatch
import datetime

from copy import deepcopy

# Other Libraries
import pyart
import scipy
import netCDF4
import numpy as np

from scipy import ndimage, signal, integrate
from csu_radartools import csu_kdp
from scipy.optimize import fmin_l_bfgs_b


def _mask_rhohv(radar, rhohv_name, tight=True):
    nrays = radar.nrays
    ngate = radar.ngates
    oneray = np.zeros((ngate))
    oneray[:(ngate // 2)] = 1 - np.linspace(0.1, 0.7, ngate // 2)
    oneray[(ngate // 2):] = 0.3
    emr = np.vstack([oneray for e in range(nrays)])
    rho = radar.fields[rhohv_name]['data']
    emr2 = np.zeros(rho.shape)
    emr2[rho > emr] = 1
    return emr2


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
    range_grid, azi_grid = np.meshgrid(radar.range['data'], radar.azimuth['data'])
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


def _get_noise_threshold(filtered_data):
    """
    Compute the noise threshold.
    """
    n, bins = np.histogram(filtered_data, bins=150)
    peaks = scipy.signal.find_peaks_cwt(n, np.array([10]))
    centers = bins[0:-1] + (bins[1] - bins[0])
    search_data = n[peaks[0]:peaks[1]]
    search_centers = centers[peaks[0]:peaks[1]]
    locs = search_data.argsort()
    noise_threshold = search_centers[locs[0]]

    return noise_threshold


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
        is_good: bool
            True if radar has a proper azimuth field.
    """
    is_good = True
    dbz = radar.fields[refl_field_name]['data']

    if dbz.shape[0] < 360:
        is_good = False

    return is_good


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
        is_good: bool
            True if radar has a proper azimuth field.
    """
    is_good = True
    dbz = radar.fields[refl_field_name]['data']

    if np.ma.isMaskedArray(dbz):
        if dbz.count() == 0:
            # Reflectivity field is empty.
            is_good = False

    return is_good


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
    try:
        rhohv = rhohv.filled(0)
    except Exception:
        pass

    rhohv[rhohv < 0] = 0

    snr = radar.fields[snr_name]['data']
    natural_snr = 10**(0.1 * snr)
    rho_corr = rhohv * (1 + 1 / natural_snr)

    return rho_corr


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
    zdr = radar.fields[zdr_name]['data']
    snr = radar.fields[snr_name]['data']
    alpha = 1.48
    natural_zdr = 10**(0.1 * zdr)
    natural_snr = 10**(0.1 * snr)
    corr_zdr = 10 * np.log10((alpha * natural_snr * natural_zdr) / (alpha * natural_snr + alpha - natural_zdr))

    return corr_zdr


def do_gatefilter(radar, refl_name='DBZ', rhohv_name='RHOHV_CORR', ncp_name='NCP',
                  zdr_name="ZDR", is_rhohv_fake=False):
    """
    Basic filtering

    Parameters:
    ===========
        radar:
            Py-ART radar structure.
        refl_name: str
            Reflectivity field name.
        rhohv_name: str
            Cross correlation ratio field name.

    Returns:
    ========
        gf_despeckeld: GateFilter
            Gate filter (excluding all bad data).
    """
    # For CPOL, there is sometime an issue with older seasons.
    gf = pyart.filters.GateFilter(radar)

    gf.exclude_outside(zdr_name, -3.0, 8.0)
    gf.exclude_outside(refl_name, -40.0, 80.0)

    radar_date = netCDF4.num2date(radar.time['data'][0], radar.time['units'])

    if not is_rhohv_fake:
        if radar_date.year not in [2006, 2007]:
            emr2 = _mask_rhohv(radar, rhohv_name, True)
        else:
            emr2 = _mask_rhohv(radar, rhohv_name, False)
        radar.add_field_like(rhohv_name, "EMR2", emr2, replace_existing=True)
        gf.exclude_not_equal("EMR2", 1)
        radar.fields.pop('EMR2')

    gf_despeckeld = pyart.correct.despeckle_field(radar, refl_name, gatefilter=gf)

    return gf_despeckeld


def do_txt_gatefilter(radar, radar_start_date, phidp_name="PHIDP", refl_name="DBZ", rhohv_name="RHOHV_CORR", is_rhohv_fake=False):
    """
    Create gatefilter using wradlib fuzzy echo classification.

    Parameters:
    ===========
    radar:
        Py-ART radar structure.
    radar_start_date: datetime
    refl_name: str
        Reflectivity field name.
    rhohv_name: str
        Cross correlation ratio field name.

    Return:
    =======
    gf: GateFilter
        Gate filter (excluding all bad data).
    """
    def noise_th(x, max_range=90):
        n, bins = np.histogram(x.flatten(), bins=150, range=(5, max_range))
        cnt = 10
        peaks = []
        while (len(peaks) < 1) or (cnt == 0):
            peaks = scipy.signal.find_peaks_cwt(n, [cnt])
            cnt - 1

        centers = bins[0:-1] + (bins[1] - bins[0])
        search_data = n[peaks[0]:peaks[1]]
        search_centers = centers[peaks[0]:peaks[1]]
        locs = search_data.argsort()
        noise_threshold = search_centers[locs[0]]

        return noise_threshold

    phi = radar.fields[phidp_name]['data'].copy()
    vel_dict = pyart.util.angular_texture_2d(phi, 4, phi.max())

    try:
        noise_threshold = noise_th(vel_dict)
    except Exception:
        print("ERROR NOISE threshold for %s. Probably only noise in volume." % (radar_start_date.isoformat()))
        return None

    radar.add_field_like(phidp_name, "TPHI", vel_dict, replace_existing=True)

    gf = pyart.filters.GateFilter(radar)
    gf.exclude_above("TPHI", noise_threshold * 1.25)

    if radar_start_date.year > 2011:
        gf.include_above(rhohv_name, 0.8)

    gf = pyart.correct.despeckle_field(radar, "DBZ", gatefilter=gf)
    radar.fields.pop('TPHI')

    return gf


def filter_hardcoding(my_array, nuke_filter, bad=-9999):
    """
    Harcoding GateFilter into an array.

    Parameters:
    ===========
        my_array: array
            Array we want to clean out.
        nuke_filter: gatefilter
            Filter we want to apply to the data.
        bad: float
            Fill value.

    Returns:
    ========
        to_return: masked array
            Same as my_array but with all data corresponding to a gate filter
            excluded.
    """
    filt_array = np.ma.masked_where(nuke_filter.gate_excluded, my_array)
    filt_array.set_fill_value(bad)
    filt_array = filt_array.filled(fill_value=bad)
    to_return = np.ma.masked_where(filt_array == bad, filt_array)
    return to_return


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


def get_field_names():
    """
    Fields name definition.

    Returns:
    ========
        fields_names: array
            Containing [(old key, new key), ...]
    """
    fields_names = [('VEL', 'velocity'),
                    ('VEL_CORR', 'corrected_velocity'),
                    ('VEL_UNFOLDED', 'region_dealias_velocity'),
                    ('TVEL', "velocity_texture"),
                    ('TPHI', "differential_phase_texture"),
                    ('DBZ', 'total_power'),
                    ('DBZ_CORR', 'corrected_reflectivity'),
                    ('RHOHV_CORR', 'RHOHV'),
                    ('RHOHV', 'cross_correlation_ratio'),
                    ('ZDR', 'differential_reflectivity'),
                    ('ZDR_CORR', 'corrected_differential_reflectivity'),
                    ('PHIDP', 'differential_phase'),
                    ('PHI_UNF', 'raw_unfolded_differential_phase'),
                    ('PHIDP_BRINGI', 'bringi_differential_phase'),
                    ('PHI_CORR', 'corrected_differential_phase'),
                    ('PHIDP_GG', 'giangrande_differential_phase'),
                    ('KDP', 'specific_differential_phase'),
                    ('KDP_BRINGI', 'bringi_specific_differential_phase'),
                    ('KDP_GG', 'giangrande_specific_differential_phase'),
                    ('WIDTH', 'spectrum_width'),
                    ('SNR', 'signal_to_noise_ratio'),
                    ('NCP', 'normalized_coherent_power')]

    return fields_names


def get_radiosoundings(sound_dir, radar_start_date):
    # Looking for radiosoundings:
    all_sonde_files = sorted(os.listdir(sound_dir))

    pos = [cnt for cnt, f in enumerate(all_sonde_files) if fnmatch.fnmatch(f, "*" + radar_start_date.strftime("%Y%m%d") + "*")]
    if len(pos) > 0:
        # Looking for the exact date.
        sonde_name = all_sonde_files[pos[0]]
        sonde_name = os.path.join(sound_dir, sonde_name)
    else:
        # Looking for the closest date.
        dtime = [datetime.datetime.strptime(re.findall("[0-9]{8}", f)[0], "%Y%m%d") for f in all_sonde_files]
        closest_date = _nearest(dtime, radar_start_date)
        sonde_name = glob.glob(os.path.join(soundings_dir, "*" + closest_date.strftime("%Y%m%d") + "*"))[0]

    return sonde_name


def get_simulated_wind_profile(radar, radiosonde_fname, height_name="height", speed_name="wspeed", wdir_name="wdir"):
    """
    Simulate the horizontal wind profile for the radar.

    Parameters
    ==========
    radar:
        Py-ART radar data structure.
    radiosonde_fname: str
        Radiosonde file name.

    Returns:
    ========
    sim_vel: dict
        Simulated velocity.
    """
    interp_sonde = netCDF4.Dataset(radiosonde_fname)
    hwind_prof = pyart.core.HorizontalWindProfile(interp_sonde[height_name], interp_sonde[speed_name], interp_sonde[wdir_name],)
    sim_vel = pyart.util.simulated_vel_from_profile(radar, hwind_prof)

    return sim_vel


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
    # Extract data
    dp = radar.fields[unfold_phidp_name]['data'].filled(-9999)
    dz = radar.fields[refl_field]['data']
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


def phidp_texture(radar, phidp_name='PHIDP'):
    """
    Compute velocity texture using new Bobby Jackson function in Py-ART.

    Parameters:
    ===========
    radar:
        Py-ART radar structure.
    vel_name: str
        Name of the (original) Doppler velocity field.

    Returns:
    ========
    vdop_vel: dict
        Velocity texture.
    """

    v_nyq_vel = radar.fields[phidp_name]['data'].max()

    tphi_dict = pyart.retrieve.calculate_velocity_texture(radar, phidp_name, nyq=v_nyq_vel, check_nyq_uniform=False)
    tphi_dict['long_name'] = "Differential phase texture"
    tphi_dict['standard_name'] = "texture_of_differential_phase"
    tphi_dict['units'] = "deg"

    return tphi_dict


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
        myfields = [('NCPH', "NCP"),
                    ('normalized_coherent_power', "NCP"),
                    ('DBZH', "DBZ"),
                    ("DBZH_CLEAN", "DBZ"),
                    ('reflectivity', "DBZ"),
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


def rename_radar_fields(radar):
    """
    Rename radar fields from their old name to the Py-ART default name.

    Parameter:
    ==========
        radar:
            Py-ART radar structure.

    Returns:
    ========
        radar:
            Py-ART radar structure.
    """
    fields_names = get_field_names()

    # Try to remove occasional fields.
    try:
        vdop_art = radar.fields['PHIDP_CORR']
        radar.add_field('PHIDP', radar.fields.pop('PHIDP_CORR'), replace_existing=True)
    except KeyError:
        pass

    # Parse array old_key, new_key
    for old_key, new_key in fields_names:
        try:
            radar.add_field(new_key, radar.fields.pop(old_key), replace_existing=True)
        except KeyError:
            continue

    return radar


def snr_and_sounding(radar, radar_start_date, sonde_name, refl_field_name='DBZ', temp_field_name="temp"):
    """
    Compute the signal-to-noise ratio as well as interpolating the radiosounding
    temperature on to the radar grid. The function looks for the radiosoundings
    that happened at the closest time from the radar. There is no time
    difference limit.

    Parameters:
    ===========
        radar:
        radar_start_date: datetime
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
    # Altitude hack.
    true_alt = radar.altitude['data'].copy()
    radar.altitude['data'] = np.array([0])

    # print("Reading radiosounding %s" % (sonde_name))
    interp_sonde = netCDF4.Dataset(sonde_name)
    temperatures = interp_sonde.variables[temp_field_name][:]
    temperatures[(temperatures < -100) | (temperatures > 100)] = np.NaN
    times = interp_sonde.variables['time'][:]
    heights = interp_sonde.variables['height'][:]

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


def unfold_velocity(radar, my_gatefilter, bobby_params=False, constrain_sounding=False, vel_name='VEL', rhohv_name='RHOHV_CORR',
                    sounding_name='sim_velocity'):
    """
    Unfold Doppler velocity using Py-ART region based algorithm. Automatically
    searches for a folding-corrected velocity field.

    Parameters:
    ===========
        radar:
            Py-ART radar structure.
        my_gatefilter:
            GateFilter
        bobby_params: bool
            Using dealiasing parameters from Bobby Jackson. Otherwise using
            defaults configuration.
        constrain_sounding: bool
            Use optimization to constrain wind field according to a sounding. Useful if
            radar scan has regions overcorrected by a Nyquist interval.
        vel_name: str
            Name of the (original) Doppler velocity field.
        sounding_name: str
            Name of the wind field derived from a sounding

    Returns:
    ========
        vdop_vel: dict
            Unfolded Doppler velocity.
    """
    # Minimize cost function that is sum of difference between regions and
    def cost_function(nyq_vector):
        cost = 0
        i = 0
        for reg in np.unique(regions[sweep_slice]):
            add_value = np.abs(np.ma.mean(vels_slice[regions[sweep_slice] == reg]) + nyq_vector[i] * v_nyq_vel -
                               np.ma.mean(svels_slice[regions[sweep_slice] == reg]))

            if np.isfinite(add_value):
                cost += add_value
            i = i + 1
        return cost

    def gradient(nyq_vector):
        gradient_vector = np.zeros(len(nyq_vector))
        i = 0

        for reg in np.unique(regions[sweep_slice]):
            add_value = (np.ma.mean(vels_slice[regions[sweep_slice] == reg]) + nyq_vector[i] * v_nyq_vel -
                         np.ma.mean(svels_slice[regions[sweep_slice] == reg]))
            if(add_value > 0):
                gradient_vector[i] = v_nyq_vel
            else:
                gradient_vector[i] = -v_nyq_vel
            i = i + 1
        return gradient_vector

    gf = deepcopy(my_gatefilter)
    # Trying to determine Nyquist velocity
    try:
        v_nyq_vel = radar.instrument_parameters['nyquist_velocity']['data'][0]
    except Exception:
        vdop_art = radar.fields[vel_name]['data']
        v_nyq_vel = np.max(np.abs(vdop_art))

    try:
        vdop_vel = pyart.correct.dealias_region_based(radar, vel_field=vel_name, rays_wrap_around=True,
                                                      gatefilter=gf, skip_between_rays=2000)
    except Exception:
        vdop_vel = pyart.correct.dealias_region_based(radar, vel_field=vel_name,
                                                      gatefilter=gf, nyquist_vel=v_nyq_vel)

    # # Cf. mail from Bobby Jackson for skip_between_rays parameters.
    # if bobby_params:
    #     vdop_vel = pyart.correct.dealias_region_based(radar,
    #                                                   vel_field=vel_name,
    #                                                   gatefilter=gf,
    #                                                   nyquist_vel=v_nyq_vel,
    #                                                   skip_between_rays=2000)
    # else:
    #     vdop_vel = pyart.correct.dealias_region_based(radar,
    #                                                   vel_field=vel_name,
    #                                                   gatefilter=gf,
    #                                                   nyquist_vel=v_nyq_vel)

    if constrain_sounding:
        gfilter = gf.gate_excluded
        vels = deepcopy(vdop_vel['data'])
        vels_uncorr = radar.fields[vel_name]['data']
        sim_vels = radar.fields[sounding_name]['data']
        diff = (sim_vels - vels) / v_nyq_vel
        region_means = []
        regions = np.zeros(vels.shape)

        for nsweep, sweep_slice in enumerate(radar.iter_slice()):
            sfilter = gfilter[sweep_slice]
            diffs_slice = diff[sweep_slice]
            vels_slice = vels[sweep_slice]
            svels_slice = sim_vels[sweep_slice]
            vels_uncorrs = vels_uncorr[sweep_slice]
            valid_sdata = vels_uncorrs[~sfilter]
            int_splits = pyart.correct.region_dealias._find_sweep_interval_splits(
                v_nyq_vel, 3, valid_sdata, nsweep)
            regions[sweep_slice], nfeatures = pyart.correct.region_dealias._find_regions(vels_uncorrs, sfilter,
                                                                                         limits=int_splits)

        bounds_list = [(x, y) for (x, y) in zip(-5 * np.ones(nfeatures + 1), 5 * np.ones(nfeatures + 1))]
        nyq_adjustments = fmin_l_bfgs_b(cost_function, np.zeros((nfeatures + 1)), disp=True, fprime=gradient,
                                        bounds=bounds_list, maxiter=20)
        i = 0
        for reg in np.unique(regions[sweep_slice]):
            reg_mean = np.mean(diffs_slice[regions[sweep_slice] == reg])
            region_means.append(reg_mean)
            vels_slice[regions[sweep_slice] == reg] += v_nyq_vel * np.round(nyq_adjustments[0][i])
            i = i + 1

        vels[sweep_slice] = vels_slice
        vdop_vel['data'] = vels

    vdop_vel['units'] = "m/s"
    vdop_vel['standard_name'] = "corrected_radial_velocity"
    vdop_vel['description'] = "Velocity unfolded using Py-ART region based dealiasing algorithm."

    return vdop_vel


def velocity_texture(radar, vel_name='VEL'):
    """
    Compute velocity texture using new Bobby Jackson function in Py-ART.

    Parameters:
    ===========
    radar:
        Py-ART radar structure.
    vel_name: str
        Name of the (original) Doppler velocity field.

    Returns:
    ========
    vdop_vel: dict
        Velocity texture.
    """

    try:
        v_nyq_vel = radar.instrument_parameters['nyquist_velocity']['data'][0]
    except (KeyError, IndexError):
        vdop_art = radar.fields[vel_name]['data']
        v_nyq_vel = np.max(np.abs(vdop_art))

    vel_dict = pyart.retrieve.calculate_velocity_texture(radar, vel_name, nyq=v_nyq_vel)

    return vel_dict
