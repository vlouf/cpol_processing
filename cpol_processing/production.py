"""
CPOL Level 1b main production line. These are the drivers function.

@title: CPOL_PROD_1b
@author: Valentin Louf <valentin.louf@monash.edu>
@copyright: Valentin Louf (2017-)
@institution: Bureau of Meteorology

.. autosummary::
    :toctree: generated/

    process_and_save
    plot_quicklook
    production_line  => Driver function.
"""
# Python Standard Library
import gc
import os
import time
import uuid
import logging
import datetime
import traceback
import warnings

# Other Libraries -- Matplotlib must be imported before pyart.
import netCDF4
import numpy as np
import matplotlib
matplotlib.use('Agg')  # <- Reason why matplotlib is imported first.
import matplotlib.colors as colors
import matplotlib.pyplot as pl
import pyart

# Custom modules.
from .processing import attenuation
from .processing import filtering
from .processing import gridding
from .processing import hydrometeors
from .processing import phase
from .processing import radar_codes
from .processing import velocity


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


def process_and_save(radar_file_name, outpath, outpath_grid=None, figure_path=None,
                     sound_dir=None, instrument='CPOL', use_giangrande=True, linearz=True):
    """
    Call processing function and write data.

    Parameters:
    ===========
        radar_file_name: str
            Name of the input radar file.
        outpath: str
            Path for saving output data.
        outpath_grid: str
            Path for saving gridded data.
        figure_path: str
            Path for saving figures.
        sound_dir: str
            Path to radiosoundings directory.
        instrument: str
            Name of radar (only CPOL will change something).
        use_giangrande: bool
            Use Scott G's technique for unfolding PHIDP and KDP.
        linearz: bool
            Gridding reflectivity in linear unit (True) or dBZ (False).
    """
    if instrument == 'CPOL':
        is_cpol = True
    else:
        is_cpol = False

    # Create directories.
    _mkdir(outpath)

    if outpath_grid is None:
        outpath_grid = os.path.join(outpath, 'GRIDDED')
    _mkdir(outpath_grid)

    if figure_path is not None:
        _mkdir(figure_path)

    outdir_150km = os.path.join(outpath_grid, "GRID_150km_2500m")
    # outdir_70km = os.path.join(outpath_grid, "GRID_70km_1000m")
    outdir_150km_highres = os.path.join(outpath_grid, "GRID_150km_1000m")
    _mkdir(outdir_150km)
    # _mkdir(outdir_70km)
    _mkdir(outdir_150km_highres)

    # Get logger.
    logger = logging.getLogger()
    tick = time.time()
    radar = production_line(radar_file_name, sound_dir, figure_path, is_cpol=is_cpol, use_giangrande=use_giangrande)
    if radar is None:
        print(f'{radar_file_name} has not been processed. Check logs.')
        return None

    radar_start_date = netCDF4.num2date(radar.time['data'][0], radar.time['units'].replace("since", "since "))

    # Generate output file name.
    outfilename = os.path.basename(radar_file_name)
    if instrument == 'CPOL':
        outfilename = "twpcpolppiX1.c1.{}.nc".format(radar_start_date.strftime("%Y%m%d_%H%M"))
    else:
        outfilename = "cfrad." + radar_start_date.strftime("%Y%m%d_%H%M%S") + ".nc"

    outfilename = os.path.join(outpath, outfilename)

    # Check if output file already exists.
    if os.path.isfile(outfilename):
        logger.error('Output file already exists for: %s.', outfilename)
        print(f"Output file {outfilename} already exists.")
        return None

    if is_cpol:
        global_metadata = dict()
        global_metadata['uuid'] = str(uuid.uuid4())
        global_metadata['naming_authority'] = 'au.org.nci'
        global_metadata['source'] = "Australian Bureau of Meteorology and Monash University"
        global_metadata['processing_level'] = "L1B"
        global_metadata['disclaimer'] = "This dataset is supported by a funding from the U.S. Department of Energy " + \
                                        "as part of the Atmospheric Radiation Measurement (ARM) Climate Research " + \
                                        "Facility, an Office of Science user facility. If you use this dataset " + \
                                        "to prepare a publication, please consider offering me (Valentin Louf)" + \
                                        "co-authorship."
        global_metadata['acknowledgement'] = "This work has been supported by the U.S. Department " + \
                                             "of Energy Atmospheric Systems Research Program through " + \
                                             "the grant DE-SC0014063. Data may be freely distributed."
        global_metadata['product_version'] = datetime.datetime.now().strftime("%Y.%m")
        global_metadata['references'] = "Louf, V., A. Protat, R. A. Warren, S. M. Collis, D. B. Wolff, S. Raunyiar, " +\
                                        "C. Jakob, and W. A. Petersen, 2018: An integrated approach to weather " + \
                                        "radar calibration and monitoring using ground clutter and satellite " + \
                                        "comparisons. J. Atmos. Ocean. Technol., JTECH-D-18-0007.1, " +\
                                        "doi:10.1175/JTECH-D-18-0007.1."
        global_metadata['creator_name'] = "Valentin Louf"
        global_metadata['creator_email'] = "valentin.louf@monash.edu"
        global_metadata['creator_url'] = "github.com/vlouf"
        global_metadata['institution'] = "Australian Bureau of Meteorology"
        global_metadata['publisher_name'] = "NCI - National Computing Infrastructure"
        global_metadata['publisher_url'] = "nci.org.au"
        global_metadata['publisher_type'] = "institution"
        global_metadata['site_name'] = "Gunn_Pt"
        global_metadata['country'] = "Australia"
        global_metadata['state'] = "NT"

        # Lat/lon informations
        maxlon = '132.3856852067545'
        minlon = '129.70320368213441'
        maxlat = '-10.941777804922253'
        minlat = '-13.552905831511362'
        origin_altitude = '50'
        origin_latitude = '-12.249'
        origin_longitude = '131.044'

        global_metadata['geospatial_bounds'] = f"({minlon}, {maxlon}, {minlat}, {maxlat})"
        global_metadata['geospatial_lat_min'] = minlat
        global_metadata['geospatial_lat_max'] = maxlat
        global_metadata['geospatial_lat_units'] = "degrees_north"
        global_metadata['geospatial_lon_min'] = minlon
        global_metadata['geospatial_lon_max'] = maxlon
        global_metadata['geospatial_lon_units'] = "degrees_east"
        global_metadata['origin_latitude'] = origin_latitude
        global_metadata['origin_longitude'] = origin_longitude
        global_metadata['origin_altitude'] = origin_altitude

        for k, v in global_metadata.items():
            radar.metadata[k] = v

    # Write results
    pyart.io.write_cfradial(outfilename, radar, format='NETCDF4')

    # Deleting all unwanted keys for gridded product.
    logger.info("Gridding started.")
    unwanted_keys = []
    goodkeys = ['corrected_differential_reflectivity',
                'cross_correlation_ratio',
                'radar_echo_classification',
                'radar_estimated_rain_rate',
                'D0',
                'NW',
                'reflectivity',
                'velocity',
                'region_dealias_velocity']

    for mykey in radar.fields.keys():
        if mykey not in goodkeys:
            unwanted_keys.append(mykey)
    for mykey in unwanted_keys:
        radar.fields.pop(mykey)

    try:
        # Gridding (and saving)
        # Full radar range with a 2.5 km grid resolution
        gridding.gridding_radar(radar, radar_start_date, outpath=outdir_150km, rmax=145e3, xyres=2500, linearz=linearz)
        # Full radar range with a 1 km grid resolution
        gridding.gridding_radar(radar, radar_start_date, outpath=outdir_150km_highres,
                                rmax=145e3, xyres=1000, linearz=linearz)
        # Half-range with a 1 km grid resolution
        # gridding.gridding_radar(radar, radar_start_date, outpath=outdir_70km, rmax=70e3, xyres=1000, linearz=linearz)

        logger.info('Gridding done.')
    except Exception:
        traceback.print_exc()
        logging.error('Problem while gridding.')
        raise

    # Processing finished!
    logger.info('%s processed in  %0.2f s.', os.path.basename(radar_file_name), (time.time() - tick))
    print('%s processed in  %0.2fs.' % (os.path.basename(radar_file_name), (time.time() - tick)))

    return None


def plot_quicklook(radar, gatefilter, radar_date, figure_path):
    """
    Plot figure of old/new radar parameters for checking purpose.

    Parameters:
    ===========
        radar:
            Py-ART radar structure.
        gatefilter:
            The Gate filter.
        radar_date: datetime
            Datetime stucture of the radar data.
    """
    if figure_path is None:
        return None
    # Extracting year and date.
    year = str(radar_date.year)
    datestr = radar_date.strftime("%Y%m%d")
    # Path for saving Figures.
    outfile_path = os.path.join(figure_path, year, datestr)

    # Checking if output directory exists. Creating them otherwise.
    if not os.path.isdir(os.path.join(figure_path, year)):
        try:
            os.mkdir(os.path.join(figure_path, year))
        except FileExistsError:
            pass
    if not os.path.isdir(outfile_path):
        try:
            os.mkdir(outfile_path)
        except FileExistsError:
            pass

    # Checking if figure already exists.
    outfile = radar_date.strftime("%Y%m%d_%H%M") + ".png"
    outfile = os.path.join(outfile_path, outfile)

    # Initializing figure.
    with pl.style.context('seaborn-paper'):
        gr = pyart.graph.RadarDisplay(radar)
        fig, the_ax = pl.subplots(3, 3, figsize=(12, 10), sharex=True, sharey=True)
        the_ax = the_ax.flatten()
        # Plotting reflectivity
        gr.plot_ppi('total_power', ax=the_ax[0])
        the_ax[0].set_title(gr.generate_title('total_power', sweep=0, datetime_format='%Y-%m-%dT%H:%M'))

        gr.plot_ppi('reflectivity', ax=the_ax[1], gatefilter=gatefilter)
        the_ax[1].set_title(gr.generate_title('reflectivity', sweep=0, datetime_format='%Y-%m-%dT%H:%M'))

        # gr.plot_ppi('radar_echo_classification', ax=the_ax[2], gatefilter=gatefilter)
        # the_ax[2].set_title(gr.generate_title('radar_echo_classification', sweep=0, datetime_format='%Y-%m-%dT%H:%M'))

        # gr.plot_ppi('radar_estimated_rain_rate', ax=the_ax[2])
        # the_ax[2].set_title(gr.generate_title('radar_estimated_rain_rate', sweep=0, datetime_format='%Y-%m-%dT%H:%M'))

        gr.plot_ppi('differential_reflectivity', ax=the_ax[2])
        the_ax[2].set_title(gr.generate_title('differential_reflectivity', sweep=0, datetime_format='%Y-%m-%dT%H:%M'))

        gr.plot_ppi('corrected_differential_reflectivity', ax=the_ax[3], gatefilter=gatefilter)
        the_ax[3].set_title(gr.generate_title('corrected_differential_reflectivity',
                                              sweep=0, datetime_format='%Y-%m-%dT%H:%M'))

        gr.plot_ppi('differential_phase', ax=the_ax[4], vmin=-180, vmax=180, cmap='pyart_Wild25')
        the_ax[4].set_title(gr.generate_title('differential_phase', sweep=0, datetime_format='%Y-%m-%dT%H:%M'))

        try:
            gr.plot_ppi('corrected_differential_phase', ax=the_ax[5], vmin=-180, vmax=180, cmap='pyart_Wild25')
            the_ax[5].set_title(gr.generate_title('giangrande_differential_phase', sweep=0,
                                                  datetime_format='%Y-%m-%dT%H:%M'))
        except KeyError:
            pass

        try:
            gr.plot_ppi('cross_correlation_ratio', ax=the_ax[6], vmin=0.5, vmax=1.05)
            the_ax[6].set_title(gr.generate_title('cross_correlation_ratio',
                                                  sweep=0, datetime_format='%Y-%m-%dT%H:%M'))
        except KeyError:
            pass

        try:
            gr.plot_ppi('folded_velocity', ax=the_ax[7], cmap=pyart.graph.cm.NWSVel, vmin=-30, vmax=30)
            the_ax[7].set_title(gr.generate_title('velocity', sweep=0, datetime_format='%Y-%m-%dT%H:%M'))
        except KeyError:
            pass

        try:
            gr.plot_ppi('velocity', ax=the_ax[8], gatefilter=gatefilter,
                        cmap=pyart.graph.cm.NWSVel, vmin=-30, vmax=30)
            the_ax[8].set_title(gr.generate_title('region_dealias_velocity', sweep=0,
                                                  datetime_format='%Y-%m-%dT%H:%M'))
        except KeyError:
            pass

        for ax_sl in the_ax:
            gr.plot_range_rings([50, 100, 150], ax=ax_sl)
            ax_sl.set_aspect(1)
            ax_sl.set_xlim(-150, 150)
            ax_sl.set_ylim(-150, 150)

        pl.tight_layout()
        pl.savefig(outfile)  # Saving figure.
        fig.clf()  # Clear figure
        pl.close()  # Release memory
    del gr  # Releasing memory

    return None


def production_line(radar_file_name, sound_dir, figure_path=None, is_cpol=True, use_giangrande=False):
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
        figure_path: str
            Path for saving figures.

    Returns:
    ========
        radar: Object
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
        10/ Process and unfold raw PHIDP using wradlib and Vulpiani algorithm.
        11/ Compute Giangrande's PHIDP using pyart.
        12/ Unfold velocity using pyart.
        13/ Compute attenuation for ZH
        14/ Compute attenuation for ZDR
        15/ Estimate Hydrometeors classification using csu toolbox.
        16/ Estimate Rainfall rate using csu toolbox.
        17/ Estimate DSD retrieval using csu toolbox.
        18/ Removing fake/temporary fieds.
        19/ Rename fields to pyart standard names.
        20/ Plotting figure quicklooks.
        21/ Hardcoding gatefilter.
    """
    # Get logger.
    logger = logging.getLogger()

    # Start chronometer.
    start_time = time.time()

    # !!! READING THE RADAR !!!
    radar = radar_codes.read_radar(radar_file_name)

    # Check if radar reflecitivity field is correct.
    if not radar_codes.check_reflectivity(radar):
        logger.error("MAJOR ERROR: %s reflectivity field is empty.", radar_file_name)
        print(f"MAJOR ERROR: {radar_file_name} reflectivity field is empty.")
        return None

    if not radar_codes.check_azimuth(radar):
        logger.error("MAJOR ERROR: %s azimuth field is empty.", radar_file_name)
        print(f"MAJOR ERROR: {radar_file_name} azimuth field is empty.")
        return None

    new_azimuth, azi_has_changed = radar_codes.correct_azimuth(radar)
    if azi_has_changed:
        logger.info('Azimuth has been corrected.')
        radar.azimuth['data'] = new_azimuth

    # Getting radar's date and time.
    radar_start_date = netCDF4.num2date(radar.time['data'][0], radar.time['units'].replace("since", "since "))
    logger.info("%s read.", radar_file_name)
    radar.time['units'] = radar.time['units'].replace("since", "since ")

    # Get radiosoundings:
    if sound_dir is not None:
        radiosonde_fname = radar_codes.get_radiosoundings(sound_dir, radar_start_date)

    # Correct Doppler velocity units.
    try:
        radar.fields['VEL']['units'] = "m/s"
        radar.fields['VEL']['standard_name'] = "radial_velocity"
        vel_missing = False
    except KeyError:
        vel_missing = True
        pass

    # Looking for RHOHV field
    # For CPOL, season 09/10, there are no RHOHV fields before March!!!!
    try:
        radar.fields['RHOHV']
        fake_rhohv = False  # Don't need to delete this field cause it's legit.
    except KeyError:
        # Creating a fake RHOHV field.
        rho = pyart.config.get_metadata('cross_correlation_ratio')
        rho['data'] = np.ones_like(radar.fields['DBZ']['data'])
        rho['description'] = "THIS FIELD IS FAKE. SHOULD BE REMOVED!"
        radar.add_field('RHOHV', rho)
        radar.add_field('RHOHV_CORR', rho)
        radar.metadata['debug_info'] = 'RHOHV field does not exist in RAW data. I had to use a fake RHOHV.'
        fake_rhohv = True  # We delete this fake field later.
        logger.critical("RHOHV field not found, creating a fake RHOHV")
        print(f"RHOHV field not found, creating a fake RHOHV {radar_file_name}")

    # Compute SNR and extract radiosounding temperature.
    # Requires radiosoundings
    if sound_dir is not None:
        try:
            height, temperature, snr = radar_codes.snr_and_sounding(radar, radiosonde_fname)
            radar.add_field('temperature', temperature, replace_existing=True)
            radar.add_field('height', height, replace_existing=True)
        except ValueError:
            traceback.print_exc()
            logger.error("Impossible to compute SNR")
            print(f"Impossible to compute SNR {radar_file_name}")
            return None

        # Looking for SNR
        try:
            radar.fields['SNR']
            logger.info('SNR already exists.')
        except KeyError:
            radar.add_field('SNR', snr, replace_existing=True)
            logger.info('SNR calculated.')

    # Correct RHOHV
    if not fake_rhohv:
        rho_corr = radar_codes.correct_rhohv(radar)
        radar.add_field_like('RHOHV', 'RHOHV_CORR', rho_corr, replace_existing=True)
        logger.info('RHOHV corrected.')

    # Correct ZDR
    corr_zdr = radar_codes.correct_zdr(radar)
    radar.add_field_like('ZDR', 'ZDR_CORR', corr_zdr, replace_existing=True)

    # GateFilter
    logger.info('Filtering data.')
    if is_cpol:
        gatefilter = filtering.do_gatefilter_cpol(radar,
                                                  refl_name='DBZ',
                                                  phidp_name="PHIDP",
                                                  rhohv_name='RHOHV_CORR',
                                                  zdr_name="ZDR")
    else:
        gatefilter = filtering.do_gatefilter(radar,
                                             refl_name='DBZ',
                                             phidp_name="PHIDP",
                                             rhohv_name='RHOHV_CORR',
                                             zdr_name="ZDR")

    # Check if NCP exists.
    try:
        radar.fields['NCP']
        fake_ncp = False
    except KeyError:
        ncp = pyart.config.get_metadata('normalized_coherent_power')
        ncp['data'] = np.zeros_like(radar.fields['RHOHV']['data'])
        ncp['data'][gatefilter.gate_included] = 1
        radar.add_field('NCP', ncp)
        fake_ncp = True

    phidp, kdp = phase.valentin_phase_processing(radar, gatefilter, phidp_name='PHIDP')
    radar.add_field('PHIDP_VAL', phidp, replace_existing=True)
    radar.add_field('KDP_VAL', kdp, replace_existing=True)
    radar.fields['PHIDP_VAL']['long_name'] = "corrected_differential_phase"
    radar.fields['KDP_VAL']['long_name'] = "corrected_specific_differential_phase"
    kdp_field_name = 'PHIDP_VAL'
    phidp_field_name = 'KDP_VAL'

    # if use_giangrande:
    #     phidp_gg, kdp_gg = phase.phidp_giangrande(radar, gatefilter, phidp_field='PHIDP', rhv_field='RHOHV_CORR')
    #     radar.add_field('PHIDP_GG', phidp_gg, replace_existing=True)
    #     radar.add_field('KDP_GG', kdp_gg, replace_existing=True)
    #     radar.fields['PHIDP_GG']['long_name'] = "corrected_differential_phase"
    #     radar.fields['KDP_GG']['long_name'] = "corrected_specific_differential_phase"
    #     logger.info('KDP/PHIDP Giangrande estimated.')
    #
    #     kdp_field_name = 'KDP_GG'
    #     phidp_field_name = 'PHIDP_GG'
    # else:
    #     # Bringi unfolding.
    #     phimeta, kdpmeta = phase.phidp_bringi(radar, gatefilter, unfold_phidp_name="PHIDP")
    #     radar.add_field('PHIDP_BRINGI', phimeta, replace_existing=True)
    #     radar.add_field('KDP_BRINGI', kdpmeta, replace_existing=True)
    #     radar.fields['PHIDP_BRINGI']['long_name'] = "corrected_differential_phase"
    #     radar.fields['KDP_BRINGI']['long_name'] = "corrected_specific_differential_phase"
    #     logger.info('KDP/PHIDP Bringi estimated.')
    #
    #     kdp_field_name = 'KDP_BRINGI'
    #     phidp_field_name = 'PHIDP_BRINGI'

    # Unfold VELOCITY
    if not vel_missing:
        # Dealias velocity.
        vdop_unfold = velocity.unfold_velocity(radar, gatefilter, constrain_sounding=False)
        radar.add_field('VEL_UNFOLDED', vdop_unfold, replace_existing=True)

        logger.info('Doppler velocity unfolded.')

    # Correct Attenuation ZH
    atten_spec, zh_corr = attenuation.correct_attenuation_zh_pyart(radar, phidp_field=phidp_field_name)
    radar.add_field('DBZ_CORR', zh_corr, replace_existing=True)
    radar.add_field('specific_attenuation_reflectivity', atten_spec, replace_existing=True)
    logger.info('Attenuation on reflectivity corrected.')

    # Correct Attenuation ZDR
    zdr_corr = attenuation.correct_attenuation_zdr(radar, phidp_name=phidp_field_name, kdp_name=kdp_field_name)
    radar.add_field_like('ZDR', 'ZDR_CORR', zdr_corr, replace_existing=True)
    # radar.add_field('specific_attenuation_differential_reflectivity', atten_spec_zdr,
    #                 replace_existing=True)
    logger.info('Attenuation on ZDR corrected.')

    # Hydrometeors classification
    hydro_class = hydrometeors.hydrometeor_classification(radar, kdp_name=kdp_field_name)
    radar.add_field('radar_echo_classification', hydro_class, replace_existing=True)
    logger.info('Hydrometeors classification estimated.')

    # Rainfall rate
    rainfall = hydrometeors.rainfall_rate(radar, gatefilter, kdp_name=kdp_field_name,
                                          refl_name='DBZ_CORR', zdr_name='ZDR_CORR')
    radar.add_field("radar_estimated_rain_rate", rainfall)
    logger.info('Rainfall rate estimated.')

    # DSD retrieval
    nw_dict, d0_dict = hydrometeors.dsd_retrieval(radar, gatefilter, kdp_name=kdp_field_name)
    radar.add_field("D0", d0_dict)
    radar.add_field("NW", nw_dict)
    logger.info('DSD estimated.')

    # Merhala classification:
#     echo_class = hydrometeors.merhala_class_convstrat(radar)
#     radar.add_field('thurai_echo_classification', echo_class)
#     logger.info('Thurai classification estimated.')

    # Removing fake and useless fields.
    if fake_ncp:
        radar.fields.pop('NCP')

    if fake_rhohv:
        radar.fields.pop("RHOHV")
        radar.fields.pop("RHOHV_CORR")

    # Rename fields to pyart defaults.
    radar = radar_codes.rename_radar_fields(radar)

    # Remove obsolete fields:
    for obsolete_key in ["Refl", "PHI_UNF", "PHI_CORR", "height", 'TH', 'TV',
                         'RHOHV']:
        try:
            radar.fields.pop(obsolete_key)
        except KeyError:
            continue

    # Treatment is finished!
    end_time = time.time()
    logger.info("Treatment for %s done in %0.2f seconds.", os.path.basename(radar_file_name), (end_time - start_time))
    print("Treatment for %s done in %0.2f seconds." % (os.path.basename(radar_file_name), (end_time - start_time)))

    # Plot check figure.
    logger.info('Plotting figure')
    try:
        plot_quicklook(radar, gatefilter, radar_start_date, figure_path)
    except Exception:
        logger.exception("Problem while trying to plot figure.")
    figure_time = time.time()
    logger.info('Figure saved in %0.2fs.', (figure_time - end_time))

    hardcode_keys = ["reflectivity",
                     "radar_echo_classification",
                     "corrected_differential_reflectivity",
                     "region_dealias_velocity",
                     "D0", "NW"]
    for mykey in hardcode_keys:
        try:
            radar.fields[mykey]['data'] = filtering.filter_hardcoding(radar.fields[mykey]['data'], gatefilter)
        except KeyError:
            continue
    logger.info('Hardcoding gatefilter to Fields done.')

    return radar
