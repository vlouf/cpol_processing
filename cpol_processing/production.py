"""
CPOL Level 1b main production line. These are the drivers function.

@title: production
@author: Valentin Louf
@email: valentin.louf@bom.gov.au
@copyright: Valentin Louf (2017-2021)
@institution: Bureau of Meteorology and Monash University
@date: 25/03/2021

.. autosummary::
    :toctree: generated/

    _mkdir
    buffer
    process_and_save
    production_line
"""
# Python Standard Library
import gc
import os
import time
import uuid
import datetime
import traceback
import warnings

# Other Libraries
import pyart
import cftime
import numpy as np

# Custom modules.
from . import attenuation
from . import cfmetadata
from . import filtering
from . import hydrometeors
from . import phase
from . import radar_codes
from . import velocity


def _mkdir(dir):
    """
    Make directory. Might seem redundant but you might have concurrency issue
    when dealing with multiprocessing.
    """
    if os.path.exists(dir):
        return None

    try:
        os.mkdir(dir)
    except FileExistsError:
        pass

    return None


def buffer(func):
    """
    Decorator to catch and kill error message. Almost want to name the function
    dont_fail.
    """

    def wrapper(*args, **kwargs):
        try:
            rslt = func(*args, **kwargs)
        except Exception:
            traceback.print_exc()
            rslt = None
        return rslt

    return wrapper


@buffer
def process_and_save(
    radar_file_name: str, outpath: str, sound_dir: str = None, do_dealiasing: bool = True, instrument: str = "CPOL",
) -> None:
    """
    Call processing function and write data.

    Parameters:
    ===========
    radar_file_name: str
        Name of the input radar file.
    outpath: str
        Path for saving output data.
    sound_dir: str
        Path to radiosoundings directory.
    instrument: str
        Name of radar (only CPOL will change something).
    do_dealiasing: bool
        Dealias velocity.
    """
    today = datetime.datetime.utcnow()
    if instrument == "CPOL":
        is_cpol = True
    else:
        is_cpol = False

    # Create directories.
    _mkdir(outpath)
    outpath = os.path.join(outpath, "v{}".format(today.strftime("%Y")))
    _mkdir(outpath)
    outpath_ppi = os.path.join(outpath, "ppi")
    _mkdir(outpath_ppi)
    tick = time.time()

    # Business start here.
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        radar = production_line(radar_file_name, sound_dir, is_cpol=is_cpol, do_dealiasing=do_dealiasing)
    # Business over.

    if radar is None:
        print(f"{radar_file_name} has not been processed. Check logs.")
        return None

    radar_start_date = cftime.num2pydate(radar.time["data"][0], radar.time["units"])
    radar_end_date = cftime.num2pydate(radar.time["data"][-1], radar.time["units"])
    outpath_ppi = os.path.join(outpath_ppi, str(radar_start_date.year))
    _mkdir(outpath_ppi)
    outpath_ppi = os.path.join(outpath_ppi, radar_start_date.strftime("%Y%m%d"))
    _mkdir(outpath_ppi)

    # Generate output file name.
    if instrument == "CPOL":
        outfilename = "twp10cpolppi.b1.{}00.nc".format(radar_start_date.strftime("%Y%m%d.%H%M"))
    else:
        outfilename = "cfrad." + radar_start_date.strftime("%Y%m%d_%H%M%S") + ".nc"

    outfilename = os.path.join(outpath_ppi, outfilename)

    # Check if output file already exists.
    if os.path.isfile(outfilename):
        print(f"Output file {outfilename} already exists.")
        return None

    if is_cpol:
        # Lat/lon informations
        latitude = radar.gate_latitude["data"]
        longitude = radar.gate_longitude["data"]
        maxlon = longitude.max()
        minlon = longitude.min()
        maxlat = latitude.max()
        minlat = latitude.min()
        origin_altitude = "50"
        origin_latitude = "-12.2491"
        origin_longitude = "131.0444"

        unique_id = str(uuid.uuid4())
        metadata = {
            "Conventions": "CF-1.6, ACDD-1.3",
            "acknowledgement": "This work has been supported by the U.S. Department of Energy Atmospheric Systems Research Program through the grant DE-SC0014063. Data may be freely distributed.",
            "country": "Australia",
            "creator_email": "CPOL-support@bom.gov.au",
            "creator_name": "Commonwealth of Australia, Bureau of Meteorology, Science and Innovation, Research, Weather and Environmental Prediction, Radar Science and Nowcasting",
            "creator_url": "http://www.bom.gov.au/australia/radar/",
            "date_created": today.isoformat(),
            "geospatial_bounds": f"POLYGON(({minlon:0.6} {minlat:0.6},{minlon:0.6} {maxlat:0.6},{maxlon:0.6} {maxlat:0.6},{maxlon:0.6} {minlat:0.6},{minlon:0.6} {minlat:0.6}))",
            "geospatial_lat_max": f"{maxlat:0.6}",
            "geospatial_lat_min": f"{minlat:0.6}",
            "geospatial_lat_units": "degrees_north",
            "geospatial_lon_max": f"{maxlon:0.6}",
            "geospatial_lon_min": f"{minlon:0.6}",
            "geospatial_lon_units": "degrees_east",
            "history": "created by Valentin Louf on raijin.nci.org.au at " + today.isoformat() + " using Py-ART",
            "id": unique_id,
            "institution": "Bureau of Meteorology",
            "instrument": "radar",
            "instrument_name": "CPOL",
            "instrument_type": "radar",
            "keywords": "radar, tropics, Doppler, dual-polarization",
            "license": "CC BY-NC-SA 4.0",
            "naming_authority": "au.gov.bom",
            "origin_altitude": origin_altitude,
            "origin_latitude": origin_latitude,
            "origin_longitude": origin_longitude,
            "platform_is_mobile": "false",
            "processing_level": "b1",
            "project": "CPOL",
            "publisher_name": "NCI",
            "publisher_url": "nci.gov.au",
            "product_version": f"v{today.year}.{today.month:02}",
            "references": "doi:10.1175/JTECH-D-18-0007.1",
            "site_name": "Gunn Pt",
            "source": "radar",
            "state": "NT",
            "standard_name_vocabulary": "CF Standard Name Table v71",
            "summary": "Volumetric scan from CPOL dual-polarization Doppler radar (Darwin, Australia)",
            "time_coverage_start": radar_start_date.isoformat(),
            "time_coverage_end": radar_end_date.isoformat(),
            "time_coverage_duration": "P10M",
            "time_coverage_resolution": "PT10M",
            "title": "radar PPI volume from CPOL",
            "uuid": unique_id,
            "version": radar.metadata["version"],
        }

        radar.metadata = metadata

    # Write results
    pyart.io.write_cfradial(outfilename, radar, format="NETCDF4")
    print("%s processed in  %0.2fs." % (os.path.basename(radar_file_name), (time.time() - tick)))

    # Free memory
    del radar
    gc.collect()

    return None


def production_line(
    radar_file_name: str, sound_dir: str, is_cpol: bool = True, do_dealiasing: bool = True
) -> pyart.core.radar.Radar:
    """
    Production line for correcting and estimating CPOL data radar parameters.
    The naming convention for these parameters is assumed to be DBZ, ZDR, VEL,
    PHIDP, KDP, SNR, RHOHV, and NCP. KDP, NCP, and SNR are optional and can be
    recalculated.

    Parameters:
    ===========
    radar_file_name: str
        Name of the input radar file.
    sound_dir: str
        Path to radiosounding directory.
    is_cpol: bool
        Name of radar (only CPOL will change something).
    do_dealiasing: bool
        Dealias velocity.

    Returns:
    ========
    radar: pyart.core.radar.Radar
        Py-ART radar structure.

    PLAN:
    =====
    01/ Read input radar file.
    02/ Check if radar file OK (no problem with azimuth and reflectivity).
    03/ Get radar date.
    04/ Check if NCP field exists (creating a fake one if it doesn't)
    05/ Check if RHOHV field exists (creating a fake one if it doesn't)
    06/ Compute SNR and temperature using radiosoundings.
    07/ Correct RHOHV using Ryzhkov algorithm.
    08/ Create gatefilter (remove noise and incorrect data).
    09/ Correct ZDR using Ryzhkov algorithm.
    10/ Compute Giangrande's PHIDP using pyart.
    11/ Unfold velocity.
    12/ Compute attenuation for ZH
    13/ Compute attenuation for ZDR
    14/ Estimate Hydrometeors classification using csu toolbox.
    15/ Estimate Rainfall rate using csu toolbox.
    16/ Removing fake/temporary fieds.
    17/ Rename fields to pyart standard names.
    """
    FIELDS_NAMES = [
        ("VEL", "velocity"),
        ("VEL_UNFOLDED", "corrected_velocity"),
        ("DBZ", "total_power"),
        ("DBZ_CORR", "corrected_reflectivity"),
        ("RHOHV_CORR", "cross_correlation_ratio"),
        ("ZDR", "differential_reflectivity"),
        ("ZDR_CORR_ATTEN", "corrected_differential_reflectivity"),
        ("PHIDP", "differential_phase"),
        ("PHIDP_BRINGI", "bringi_differential_phase"),
        ("PHIDP_GG", "giangrande_differential_phase"),
        ("PHIDP_VAL", "corrected_differential_phase"),
        ("KDP", "specific_differential_phase"),
        ("KDP_BRINGI", "bringi_specific_differential_phase"),
        ("KDP_GG", "giangrande_specific_differential_phase"),
        ("KDP_VAL", "corrected_specific_differential_phase"),
        ("WIDTH", "spectrum_width"),
        ("SNR", "signal_to_noise_ratio"),
        ("NCP", "normalized_coherent_power"),
        ("DBZV", "reflectivity_v"),
        ("WRADV", "spectrum_width_v"),
        ("SNRV", "signal_to_noise_ratio_v"),
        ("SQIV", "normalized_coherent_power_v"),
    ]

    # List of keys that we'll keep in the output radar dataset.
    OUTPUT_RADAR_FLD = [
        "corrected_differential_phase",
        "corrected_differential_reflectivity",
        "corrected_reflectivity",
        "corrected_specific_differential_phase",
        "corrected_velocity",
        "cross_correlation_ratio",
        "differential_phase",
        "differential_reflectivity",
        "radar_echo_classification",
        "radar_estimated_rain_rate",
        "signal_to_noise_ratio",
        "spectrum_width",
        "total_power",
        "velocity",
    ]

    # !!! READING THE RADAR !!!
    if is_cpol:
        radar = pyart.io.read(radar_file_name)
    else:
        radar = radar_codes.read_radar(radar_file_name)

    # Correct data type manually
    try:
        radar.longitude["data"] = np.ma.masked_invalid(radar.longitude["data"].astype(np.float32))
        radar.latitude["data"] = np.ma.masked_invalid(radar.latitude["data"].astype(np.float32))
        radar.altitude["data"] = np.ma.masked_invalid(radar.altitude["data"].astype(np.int32))
    except Exception:
        pass

    # Check if radar reflecitivity field is correct.
    if not radar_codes.check_reflectivity(radar):
        raise TypeError(f"Reflectivity field is empty in {radar_file_name}.")

    if not radar_codes.check_azimuth(radar):
        raise TypeError(f"Azimuth field is empty in {radar_file_name}.")

    if not radar_codes.check_year(radar):
        print(f"{radar_file_name} date probably wrong. Had to correct century.")

    new_azimuth, azi_has_changed = radar_codes.correct_azimuth(radar)
    if azi_has_changed:
        radar.azimuth["data"] = new_azimuth

    # Getting radar's date and time.
    radar_start_date = cftime.num2pydate(radar.time["data"][0], radar.time["units"])
    radar.time["units"] = radar.time["units"].replace("since", "since ")

    # Get radiosoundings:
    if sound_dir is not None:
        radiosonde_fname = radar_codes.get_radiosoundings(sound_dir, radar_start_date)

    # Correct Doppler velocity units.
    try:
        radar.fields["VEL"]["units"] = "m s-1"
        vel_missing = False
    except KeyError:
        vel_missing = True

    # Looking for RHOHV field
    # For CPOL, season 09/10, there are no RHOHV fields before March!!!!
    try:
        radar.fields["RHOHV"]
        fake_rhohv = False  # Don't need to delete this field cause it's legit.
    except KeyError:
        # Creating a fake RHOHV field.
        fake_rhohv = True  # We delete this fake field later.
        rho = pyart.config.get_metadata("cross_correlation_ratio")
        rho["data"] = np.ones_like(radar.fields["DBZ"]["data"])
        radar.add_field("RHOHV", rho)
        radar.add_field("RHOHV_CORR", rho)

    # Compute SNR and extract radiosounding temperature.
    # Requires radiosoundings
    if sound_dir is not None:
        try:
            height, temperature, snr = radar_codes.snr_and_sounding(radar, radiosonde_fname)
            radar.add_field("temperature", temperature, replace_existing=True)
            radar.add_field("height", height, replace_existing=True)
        except ValueError:
            traceback.print_exc()
            print(f"Impossible to compute SNR {radar_file_name}")
            return None

        # Looking for SNR
        try:
            radar.fields["SNR"]
        except KeyError:
            radar.add_field("SNR", snr, replace_existing=True)

    # Correct RHOHV
    if not fake_rhohv:
        rho_corr = radar_codes.correct_rhohv(radar)
        radar.add_field_like("RHOHV", "RHOHV_CORR", rho_corr, replace_existing=True)

    # Correct ZDR
    corr_zdr = radar_codes.correct_zdr(radar)
    radar.add_field_like("ZDR", "ZDR_CORR", corr_zdr, replace_existing=True)

    # GateFilter
    if is_cpol:
        gatefilter = filtering.do_gatefilter_cpol(
            radar, refl_name="DBZ", phidp_name="PHIDP", rhohv_name="RHOHV_CORR", zdr_name="ZDR"
        )
    else:
        gatefilter = filtering.do_gatefilter(
            radar, refl_name="DBZ", phidp_name="PHIDP", rhohv_name="RHOHV_CORR", zdr_name="ZDR"
        )

    # Check if NCP exists.
    try:
        radar.fields["NCP"]
        fake_ncp = False
    except KeyError:
        fake_ncp = True
        ncp = pyart.config.get_metadata("normalized_coherent_power")
        ncp["data"] = np.zeros_like(radar.fields["RHOHV"]["data"])
        ncp["data"][gatefilter.gate_included] = 1
        radar.add_field("NCP", ncp)

    phidp, kdp = phase.phidp_giangrande(radar, gatefilter)
    radar.add_field("PHIDP_VAL", phidp)
    radar.add_field("KDP_VAL", kdp)
    kdp_field_name = "KDP_VAL"
    phidp_field_name = "PHIDP_VAL"

    # Unfold VELOCITY
    if do_dealiasing:
        if not vel_missing:
            if is_cpol:
                vdop_unfold = velocity.unravel(radar, gatefilter, nyquist=13.3)
            else:
                vdop_unfold = velocity.unravel(radar, gatefilter)
            radar.add_field("VEL_UNFOLDED", vdop_unfold, replace_existing=True)

    # Correct attenuation ZH and ZDR and hardcode gatefilter
    zh_corr = attenuation.correct_attenuation_zh_pyart(radar, gatefilter, phidp_field=phidp_field_name)
    radar.add_field_like("DBZ", "DBZ_CORR", zh_corr)

    zdr_corr = attenuation.correct_attenuation_zdr(radar, gatefilter, phidp_name=phidp_field_name)
    radar.add_field("ZDR_CORR_ATTEN", zdr_corr)

    # Hydrometeors classification
    hydro_class = hydrometeors.hydrometeor_classification(
        radar, gatefilter, kdp_name=kdp_field_name, zdr_name="ZDR_CORR_ATTEN"
    )

    radar.add_field("radar_echo_classification", hydro_class, replace_existing=True)

    # Rainfall rate
    rainfall = hydrometeors.rainfall_rate(
        radar, gatefilter, kdp_name=kdp_field_name, refl_name="DBZ_CORR", zdr_name="ZDR_CORR_ATTEN"
    )
    radar.add_field("radar_estimated_rain_rate", rainfall)

    # Removing fake and useless fields.
    if fake_ncp:
        radar.fields.pop("NCP")

    if fake_rhohv:
        radar.fields.pop("RHOHV")
        radar.fields.pop("RHOHV_CORR")

    # Remove obsolete fields:
    for obsolete_key in ["Refl", "PHI_UNF", "PHI_CORR", "height", "TH", "TV", "ZDR_CORR", "RHOHV"]:
        try:
            radar.fields.pop(obsolete_key)
        except KeyError:
            continue

    # Change the temporary working name of fields to the one define by the user.
    for old_key, new_key in FIELDS_NAMES:
        try:
            radar.add_field(new_key, radar.fields.pop(old_key), replace_existing=True)
        except KeyError:
            continue

    # Delete working variables.
    if is_cpol:
        for k in list(radar.fields.keys()):
            if k not in OUTPUT_RADAR_FLD:
                radar.fields.pop(k)

    # Correct the standard_name metadata:
    cfmetadata.correct_standard_name(radar)
    # ACDD-1.3 compliant metadata:
    cfmetadata.coverage_content_type(radar)
    cfmetadata.correct_units(radar)

    return radar
