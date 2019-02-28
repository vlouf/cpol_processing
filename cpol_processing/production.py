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


def process_and_save(radar_file_name, outpath, sound_dir=None, instrument='CPOL'):
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
        linearz: bool
            Gridding reflectivity in linear unit (True) or dBZ (False).
    """
    today = datetime.datetime.utcnow()
    if instrument == 'CPOL':
        is_cpol = True
    else:
        is_cpol = False

    # Create directories.
    _mkdir(outpath)
    outpath = os.path.join(outpath, "v{}".format(today.strftime('%Y')))
    _mkdir(outpath)
    outpath_ppi = os.path.join(outpath, 'ppi')
    _mkdir(outpath_ppi)
    outpath_grid = os.path.join(outpath, 'gridded')
    _mkdir(outpath_grid)
    figure_path = os.path.join(outpath, 'quicklooks')
    _mkdir(figure_path)

    outdir_150km = os.path.join(outpath_grid, "grid_150km_2500m")
    outdir_70km = os.path.join(outpath_grid, "grid_70km_1000m")
    _mkdir(outdir_150km)
    _mkdir(outdir_70km)

    # Get logger.
    logger = logging.getLogger()
    tick = time.time()
    # Business start here.
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        radar = production_line(radar_file_name, sound_dir, figure_path, is_cpol=is_cpol)
    # Business over.
    if radar is None:
        print(f'{radar_file_name} has not been processed. Check logs.')
        return None

    radar_start_date = netCDF4.num2date(radar.time['data'][0], radar.time['units'])
    outpath_ppi = os.path.join(outpath_ppi, str(radar_start_date.year))
    _mkdir(outpath_ppi)
    outpath_ppi = os.path.join(outpath_ppi, radar_start_date.strftime('%Y%m%d'))
    _mkdir(outpath_ppi)

    # Generate output file name.
    if instrument == 'CPOL':
        outfilename = "twp10cpolppi.b1.{}00.nc".format(radar_start_date.strftime("%Y%m%d.%H%M"))
    else:
        outfilename = "cfrad." + radar_start_date.strftime("%Y%m%d_%H%M%S") + ".nc"

    outfilename = os.path.join(outpath_ppi, outfilename)

    # Check if output file already exists.
    if os.path.isfile(outfilename):
        logger.error('Output file already exists for: %s.', outfilename)
        print(f"Output file {outfilename} already exists.")
        return None

    if is_cpol:
        # Lat/lon informations
        maxlon = '132.385'
        minlon = '129.703'
        maxlat = '-10.941'
        minlat = '-13.552'
        origin_altitude = '50'
        origin_latitude = '-12.249'
        origin_longitude = '131.044'

        metadata = dict()
        metadata['Conventions'] = radar.metadata['Conventions']
        metadata['acknowledgement'] = 'This work has been supported by the U.S. Department of Energy Atmospheric Systems Research Program through the grant DE-SC0014063. Data may be freely distributed.'
        metadata['country'] = 'Australia'
        metadata['creator_email'] = 'valentin.louf@bom.gov.au'
        metadata['creator_name'] = 'Valentin Louf'
        metadata['geospatial_bounds'] = f"({minlon}, {maxlon}, {minlat}, {maxlat})"
        metadata['geospatial_lat_max'] = maxlat
        metadata['geospatial_lat_min'] = minlat
        metadata['geospatial_lat_units'] = "degrees_north"
        metadata['geospatial_lon_max'] = maxlon
        metadata['geospatial_lon_min'] = minlon
        metadata['geospatial_lon_units'] = "degrees_east"
        metadata['history'] = "created by Valentin Louf on raijin.nci.org.au at " + today.isoformat() + " using Py-ART"
        metadata['institution'] = 'Australian Bureau of Meteorology'
        metadata['instrument_name'] = 'CPOL'
        metadata['instrument_type'] = 'radar'
        metadata['naming_authority'] = 'au.org.nci'
        metadata['origin_altitude'] = origin_altitude
        metadata['origin_latitude'] = origin_latitude
        metadata['origin_longitude'] = origin_longitude
        metadata['platform_is_mobile'] = 'false'
        metadata['processing_level'] = 'b1'
        metadata['publisher_name'] = "NCI"
        metadata['publisher_url'] = "nci.gov.au"
        metadata['references'] = 'cf. doi:10.1175/JTECH-D-18-0007.1'
        metadata['site_name'] = 'Gunn_Pt'
        metadata['source'] = 'rapic'
        metadata['state'] = "NT"
        metadata['title'] = "radar PPI volume from CPOL"
        metadata['uuid'] = str(uuid.uuid4())
        metadata['version'] = radar.metadata['version']

        radar.metadata = metadata

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

    # try:
    #     # Gridding (and saving)
    #     # Full radar range with a 2.5 km grid resolution
    #     gridding.gridding_radar(radar, radar_start_date, outpath=outdir_150km, rmax=145e3, xyres=2500, linearz=linearz)
    #     # Full radar range with a 1 km grid resolution
    #     gridding.gridding_radar(radar, radar_start_date, outpath=outdir_150km_highres,
    #                             rmax=145e3, xyres=1000, linearz=linearz)
    #     # Half-range with a 1 km grid resolution
    #     # gridding.gridding_radar(radar, radar_start_date, outpath=outdir_70km, rmax=70e3, xyres=1000, linearz=linearz)

    #     logger.info('Gridding done.')
    # except Exception:
    #     traceback.print_exc()
    #     logging.error('Problem while gridding.')
    #     raise

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
        fig, the_ax = pl.subplots(4, 3, figsize=(14, 15), sharex=True, sharey=True)
        the_ax = the_ax.flatten()
        # Plotting reflectivity
        gr.plot_ppi('total_power', ax=the_ax[0])
        the_ax[0].set_title(gr.generate_title('total_power', sweep=0, datetime_format='%Y-%m-%dT%H:%M'))

        gr.plot_ppi('reflectivity', ax=the_ax[1], gatefilter=gatefilter)
        the_ax[1].set_title(gr.generate_title('reflectivity', sweep=0, datetime_format='%Y-%m-%dT%H:%M'))

        gr.plot_ppi('radar_estimated_rain_rate', ax=the_ax[2])
        the_ax[2].set_title(gr.generate_title('radar_estimated_rain_rate', sweep=0, datetime_format='%Y-%m-%dT%H:%M'))

        gr.plot_ppi('differential_phase', ax=the_ax[3], vmin=-180, vmax=180, cmap='pyart_Wild25')
        the_ax[3].set_title(gr.generate_title('differential_phase', sweep=0, datetime_format='%Y-%m-%dT%H:%M'))

        try:
            gr.plot_ppi('corrected_differential_phase', ax=the_ax[4], vmin=-180, vmax=180, cmap='pyart_Wild25')
            the_ax[4].set_title(gr.generate_title('corrected_differential_phase', sweep=0,
                                                  datetime_format='%Y-%m-%dT%H:%M'))
        except KeyError:
            pass

        try:
            gr.plot_ppi('corrected_specific_differential_phase', ax=the_ax[5], vmin=-2, vmax=5, cmap='pyart_Theodore16')
            the_ax[5].set_title(gr.generate_title('corrected_specific_differential_phase', sweep=0,
                                                  datetime_format='%Y-%m-%dT%H:%M'))
        except KeyError:
            pass

        try:
            gr.plot_ppi('folded_velocity', ax=the_ax[6], cmap=pyart.graph.cm.NWSVel, vmin=-30, vmax=30)
            the_ax[6].set_title(gr.generate_title('velocity', sweep=0, datetime_format='%Y-%m-%dT%H:%M'))
        except KeyError:
            pass

        try:
            gr.plot_ppi('velocity', ax=the_ax[7], gatefilter=gatefilter,
                        cmap=pyart.graph.cm.NWSVel, vmin=-30, vmax=30)
            the_ax[7].set_title(gr.generate_title('region_dealias_velocity', sweep=0,
                                                  datetime_format='%Y-%m-%dT%H:%M'))
        except KeyError:
            pass

        try:
            gr.plot_ppi('cross_correlation_ratio', ax=the_ax[8], vmin=0.5, vmax=1.05)
            the_ax[8].set_title(gr.generate_title('cross_correlation_ratio',
                                                  sweep=0, datetime_format='%Y-%m-%dT%H:%M'))
        except KeyError:
            pass

        gr.plot_ppi('differential_reflectivity', ax=the_ax[9])
        the_ax[9].set_title(gr.generate_title('differential_reflectivity', sweep=0, datetime_format='%Y-%m-%dT%H:%M'))

        gr.plot_ppi('corrected_differential_reflectivity', ax=the_ax[10], gatefilter=gatefilter)
        the_ax[10].set_title(gr.generate_title('corrected_differential_reflectivity',
                                              sweep=0, datetime_format='%Y-%m-%dT%H:%M'))

        gr.plot_ppi('radar_echo_classification', ax=the_ax[11])
        the_ax[11].set_title(gr.generate_title('radar_echo_classification', sweep=0, datetime_format='%Y-%m-%dT%H:%M'))

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


def production_line(radar_file_name, sound_dir, figure_path=None, is_cpol=True):
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

    if is_cpol:
        if radar.nsweeps <= 10:
            print('Problem with CPOL PPI, not enough elevations.')
            logger.error('Problem with CPOL PPI, not enough elevations.')
            return None

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
        fake_rhohv = True  # We delete this fake field later.
        rho = pyart.config.get_metadata('cross_correlation_ratio')
        rho['data'] = np.ones_like(radar.fields['DBZ']['data'])
        radar.add_field('RHOHV', rho)
        radar.add_field('RHOHV_CORR', rho)
        logger.critical("RHOHV field not found, creating a fake RHOHV")

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
        fake_ncp = True
        ncp = pyart.config.get_metadata('normalized_coherent_power')
        ncp['data'] = np.zeros_like(radar.fields['RHOHV']['data'])
        ncp['data'][gatefilter.gate_included] = 1
        radar.add_field('NCP', ncp)

    phidp, kdp = phase.valentin_phase_processing(radar, gatefilter, phidp_name='PHIDP')
    radar.add_field('PHIDP_VAL', phidp)
    radar.add_field('KDP_VAL', kdp)
    kdp_field_name = 'KDP_VAL'
    phidp_field_name = 'PHIDP_VAL'

    # Unfold VELOCITY
    if not vel_missing:
        # Dealias velocity.
        vdop_unfold = velocity.unfold_velocity(radar, gatefilter, constrain_sounding=False)
        radar.add_field('VEL_UNFOLDED', vdop_unfold, replace_existing=True)
        logger.info('Doppler velocity unfolded.')

    # Correct gaseous attenuation
    atten_gas = attenuation.correct_gaseous_attenuation(radar)
    radar.add_field('gaseous_attenuation', atten_gas)

    # Correct Attenuation ZH
    atten_spec, zh_corr = attenuation.correct_attenuation_zh_pyart(radar, phidp_field=phidp_field_name)
    zh_corr['data'] += atten_gas['data']
    radar.add_field('DBZ_CORR', zh_corr, replace_existing=True)
    radar.add_field('specific_attenuation_reflectivity', atten_spec, replace_existing=True)
    logger.info('Attenuation on reflectivity corrected.')

    # Correct Attenuation ZDR
    zdr_corr = attenuation.correct_attenuation_zdr(radar, phidp_name=phidp_field_name, zdr_name='ZDR_CORR')
    radar.add_field('ZDR_CORR_ATTEN', zdr_corr)
    logger.info('Attenuation on ZDR corrected.')

    # Hydrometeors classification
    hydro_class = hydrometeors.hydrometeor_classification(radar,
                                                          gatefilter,
                                                          kdp_name=kdp_field_name,
                                                          zdr_name='ZDR_CORR_ATTEN')

    radar.add_field('radar_echo_classification', hydro_class, replace_existing=True)
    logger.info('Hydrometeors classification estimated.')

    # Rainfall rate
    rainfall = hydrometeors.rainfall_rate(radar, gatefilter, kdp_name=kdp_field_name,
                                          refl_name='DBZ_CORR', zdr_name='ZDR_CORR_ATTEN')
    radar.add_field("radar_estimated_rain_rate", rainfall)
    logger.info('Rainfall rate estimated.')

    # DSD retrieval
    nw_dict, d0_dict = hydrometeors.dsd_retrieval(radar, gatefilter, kdp_name=kdp_field_name, zdr_name='ZDR_CORR_ATTEN')
    radar.add_field("D0", d0_dict)
    radar.add_field("NW", nw_dict)
    logger.info('DSD estimated.')

    # Removing fake and useless fields.
    if fake_ncp:
        radar.fields.pop('NCP')

    if fake_rhohv:
        radar.fields.pop("RHOHV")
        radar.fields.pop("RHOHV_CORR")

    # Remove obsolete fields:
    for obsolete_key in ["Refl", "PHI_UNF", "PHI_CORR", "height", 'TH', 'TV', 'ZDR_CORR',
                         'RHOHV']:
        try:
            radar.fields.pop(obsolete_key)
        except KeyError:
            continue

    # Rename fields to pyart defaults.
    fields_names = [('VEL', 'folded_velocity'),
                    ('VEL_UNFOLDED', 'velocity'),
                    ('DBZ', 'total_power'),
                    ('DBZ_CORR', 'reflectivity'),
                    ('RHOHV_CORR', 'cross_correlation_ratio'),
                    ('ZDR', 'differential_reflectivity'),
                    ('ZDR_CORR_ATTEN', 'corrected_differential_reflectivity'),
                    ('PHIDP', 'differential_phase'),
                    ('PHIDP_BRINGI', 'bringi_differential_phase'),
                    ('PHIDP_GG', 'giangrande_differential_phase'),
                    ('PHIDP_VAL', 'corrected_differential_phase'),
                    ('KDP', 'specific_differential_phase'),
                    ('KDP_BRINGI', 'bringi_specific_differential_phase'),
                    ('KDP_GG', 'giangrande_specific_differential_phase'),
                    ('KDP_VAL', 'corrected_specific_differential_phase'),
                    ('WIDTH', 'spectrum_width'),
                    ('SNR', 'signal_to_noise_ratio'),
                    ('NCP', 'normalized_coherent_power'),
                    ('DBZV', 'reflectivity_v'),
                    ('WRADV', 'spectrum_width_v'),
                    ('SNRV', 'signal_to_noise_ratio_v'),
                    ('SQIV', 'normalized_coherent_power_v'),
                   ]

    for old_key, new_key in fields_names:
        try:
            radar.add_field(new_key, radar.fields.pop(old_key), replace_existing=True)
        except KeyError:
            continue

    # for key in radar.fields.keys():
    #     radar.fields[key]['data'] = radar.fields[key]['data'].astype(np.float32)

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
