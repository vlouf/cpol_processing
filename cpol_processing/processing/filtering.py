# Python Standard Library
import time
import datetime

# Other Libraries
import pyart
import scipy
import netCDF4
import numpy as np


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
