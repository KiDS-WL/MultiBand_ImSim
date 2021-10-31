# -*- coding: utf-8 -*-
# @Author: lshuns
# @Date:   2021-02-03, 15:58:35
# @Last Modified by:   lshuns
# @Last Modified time: 2021-10-31 16:11:18

### module to generate an example configuration file

import os
import re
import shutil
import logging
import configparser
import distutils.util

logger = logging.getLogger(__name__)

def ParseConfig(config_file, taskIDs, run_tag, running_log):

    # ++++++++++++++ 1. config file
    config = configparser.ConfigParser(allow_no_value=True,
                                    inline_comment_prefixes='#',
                                    empty_lines_in_values=False,
                                    interpolation=configparser.ExtendedInterpolation())
    ## existence
    try:
        config.read(config_file)
    except TypeError:
        raise Exception("No configuration file provided! \n"
                        "------> To generate an example file, use `python Run.py 0`.")

    # ++++++++++++++ 2. some general info (required by all tasks)

    ## === a. work dir
    config_paths = config['Paths']
    config_dir = config_paths.get('config_dir')

    ### check existence
    if not os.path.isdir(config_dir):
        raise Exception("configuration directory not found! \n"
                        f"------> No directory named as {config_dir}")
    out_dir = config_paths.get('out_dir')
    out_dir = os.path.join(out_dir, run_tag)
    if not os.path.exists(out_dir):
        os.mkdir(out_dir)
    logger.info(f'all outputs will be saved to {out_dir}')

    ### images
    ima_dir = os.path.join(out_dir, 'images')
    if not os.path.exists(ima_dir):
        os.mkdir(ima_dir)

    ### catalogues
    cata_dir = os.path.join(out_dir, 'catalogues')
    if not os.path.exists(cata_dir):
        os.mkdir(cata_dir)

    ### running_log
    if running_log:
        log_dir = os.path.join(out_dir, 'running_log')
        if not os.path.exists(log_dir):
            os.mkdir(log_dir)
        logger.debug(f'Running log from external codes will be saved to {log_dir}')
    else:
        log_dir = None

    ### combine to a dictionary
    work_dirs = {'main': out_dir,
                    'ima': ima_dir,
                    'cata': cata_dir,
                    'log': log_dir}

    ## === b. galaxy info
    try:
        config_gal = config['GalInfo']
        gal_configs = {'file': config_gal.get('cata_file'),
                        'position_type': config_gal.get('position_type'),
                        'mag_cut': [float(i_p.strip()) for i_p in config_gal.get('mag_cut').split(',')],
                        'size_cut': [float(i_p.strip()) for i_p in config_gal.get('Re_cut').split(',')],
                        'id_name': config_gal.get('id_name'),
                        'detection_mag_name': config_gal.get('detection_mag_name'),
                        'mag_name_list': [x.strip() for x in config_gal.get('mag_name_list').split(',')],
                        'RaDec_names': [x.strip() for x in config_gal.get('RaDec_names').split(',')],
                        'shape_names': [x.strip() for x in config_gal.get('shape_names').split(',')],
                        'z_name': config_gal.get('z_name')}

        ### grid size 
        grid_size = config_gal.get('grid_size')
        if not grid_size:
            gal_configs['grid_size'] = float(grid_size)
        else:
            gal_configs['grid_size'] = 18 #arcsec

        ### check existence
        if not os.path.isfile(gal_configs['file']):
            tmp = gal_configs['file']
            raise Exception(f"galaxies file {tmp} not found!")

    except KeyError:
        gal_configs = {'file': None}
        logger.warning('No GalInfo provided, cannot simulate any images!')

    ## === c. star info
    try:
        config_star = config['StarInfo']
        star_file = config_star.get('cata_file')

        if star_file:

            ### check existence
            if not os.path.isfile(star_file):
                raise Exception(f"star file {star_file} not found!")

            star_configs = {'file': star_file,
                            'cata_area': config_star.getfloat('cata_area'),
                            'position_type': config_star.get('position_type'),
                            'mag_cut': [float(i_p.strip()) for i_p in config_star.get('mag_cut').split(',')],
                            'id_name': config_star.get('id_name'),
                            'detection_mag_name': config_star.get('detection_mag_name'),
                            'mag_name_list': [x.strip() for x in config_star.get('mag_name_list').split(',')],
                            'RaDec_names': [x.strip() for x in config_star.get('RaDec_names').split(',')] if (config_star.get('position_type') == 'true') else None}

        else:
            star_configs = {'file': None}

    except KeyError:
        star_configs = {'file': None}
        logger.warning('No StarInfo provided, will not contain any stars!')

    ## === d. noise info
    config_noise = config['NoiseInfo']
    noise_configs = {'file': config_noise.get('cata_file'),
                    'file4varChips': config_noise.get('file4varChips')}

    ### check existence
    if not os.path.isfile(noise_configs['file']):
        tmp = noise_configs['file']
        raise Exception(f"noise file {tmp} not found!")
    if (noise_configs['file4varChips']) and (not os.path.isfile(noise_configs['file4varChips'])):
        tmp = noise_configs['file4varChips']
        raise Exception(f"separate psf file {tmp} not found!")

    ### psf profile type 
    psf_type = config_noise.get('psf_type')
    if psf_type:
        noise_configs['psf_type'] = psf_type.lower()
    else:
        noise_configs['psf_type'] = 'moffat'
    ###### supported types
    if noise_configs['psf_type'] not in ['moffat', 'airy']:
        raise Exception(f'not supported psf_type: {psf_type}')

    ### >>> for old versions                      
    noise_psf_basenames = config_noise.get('noise_psf_basenames')
    if noise_psf_basenames:
        logger.warning('Using old noise_psf_basenames parameter, please update ASAP!')
        noise_configs['noise_psf_basenames'] = [x.strip() for x in noise_psf_basenames.split(',')]

        ###### fill none for missed basenames 
        while len(noise_configs['noise_psf_basenames']) < 8:
            noise_configs['noise_psf_basenames'].append('none')

    ### >>> for new versions
    else:
        noise_configs['label_basename'] = config_noise.get('label_basename')
        noise_configs['noise_basenames'] = [x.strip() for x in config_noise.get('noise_basenames').split(',')]
        noise_configs['psf_basenames'] = [x.strip() for x in config_noise.get('psf_basenames').split(',')]
        noise_configs['id_basenames'] = [x.strip() for x in config_noise.get('id_basenames').split(',')]

    ## === e. ImSim
    config_imsim = config['ImSim']
    imsim_configs = {'survey': config_imsim.get('survey'),
                    'N_tiles': config_imsim.getint('N_tiles'),
                    'gal_rotation_angles': [float(i_r.strip()) for i_r in config_imsim.get('gal_rotation_angles').split(',')],
                    'PSF_map': [bool(distutils.util.strtobool(x.strip())) for x in config_imsim.get('PSF_map').split(',')],
                    'mag_zero': config_imsim.getfloat('mag_zero'),
                    'bands': [x.strip() for x in config_imsim.get('band_list').split(',')],
                    'pixel_scale_list': [float(i_p.strip()) for i_p in config_imsim.get('pixel_scale_list').split(',')],
                    'image_type_list': [x.strip() for x in config_imsim.get('image_type_list').split(',')],
                    'image_chips': [bool(distutils.util.strtobool(x.strip())) for x in config_imsim.get('image_chips').split(',')],
                    'image_PSF': [bool(distutils.util.strtobool(x.strip())) for x in config_imsim.get('image_PSF').split(',')]}

    ### repeat certain para to match with number of bands
    if len(imsim_configs['PSF_map']) == 1:
        imsim_configs['PSF_map'] = imsim_configs['PSF_map'] * len(imsim_configs['bands'])
    if len(imsim_configs['pixel_scale_list']) == 1:
        imsim_configs['pixel_scale_list'] = imsim_configs['pixel_scale_list'] * len(imsim_configs['bands'])
    if len(imsim_configs['image_type_list']) == 1:
        imsim_configs['image_type_list'] = imsim_configs['image_type_list'] * len(imsim_configs['bands'])
    if len(imsim_configs['image_chips']) == 1:
        imsim_configs['image_chips'] = imsim_configs['image_chips'] * len(imsim_configs['bands'])
    if len(imsim_configs['image_PSF']) == 1:
        imsim_configs['image_PSF'] = imsim_configs['image_PSF'] * len(imsim_configs['bands'])

    ### psf image size 
    psf_size = config_imsim.get('image_PSF_size')
    if psf_size:
        imsim_configs['image_PSF_size'] = float(psf_size)
    else:
        imsim_configs['image_PSF_size'] = 48 #pixels

    ## === dictionary for collecting all config info
    configs_dict = {'work_dirs': work_dirs, 'gal': gal_configs, 'star': star_configs, 'noise': noise_configs, 'imsim': imsim_configs}

    # ++++++++++++++ 2. some specific info (required only by certain tasks)

    # === swarp images
    if ('2' in taskIDs) or ('all' in taskIDs):
        config_swarp = config['SWarp']
        swarp_configs = {'cmd': config_swarp.get('cmd'),
                         'config_files': [os.path.join(config_dir, x.strip()) for x in config_swarp.get('config_files').split(',')],
                         'bands_group': re.findall(r'\[([^]]+)', config_swarp.get('bands_group')),
                         'image_label_list': [x.strip() for x in config_swarp.get('image_label_list').split(',')],
                         'only_resamples': [bool(distutils.util.strtobool(x.strip())) for x in config_swarp.get('only_resamples').split(',')],
                         'clean_up_levels': [int(x.strip()) for x in config_swarp.get('clean_up_levels').split(',')]}

        ### legitimate check
        if not shutil.which(swarp_configs['cmd']):
            tmp = swarp_configs['cmd']
            raise Exception(f"{tmp} is not an executable path for SWarp!")
        for tmp in swarp_configs['config_files']:
            if not os.path.isfile(tmp):
                raise Exception(f"SWarp config file {tmp} not found!")

        ### collect
        configs_dict['swarp'] = swarp_configs

    ## === detect objects
    if ('3' in taskIDs) or ('all' in taskIDs):
        config_sex = config['SExtractor']
        sex_configs = {'cmd': config_sex.get('cmd'),
                        'detection_band': config_sex.get('detection_band'),
                        'image_label': config_sex.get('image_label'),
                        'pixel_scale': config_sex.getfloat('pixel_scale'),
                        'cross_match': config_sex.getboolean('cross_match'),
                        'config_file': os.path.join(config_dir, config_sex.get('config_file')),
                        'param_file': os.path.join(config_dir, config_sex.get('param_file')),
                        'filter_file': os.path.join(config_dir, config_sex.get('filter_file')),
                        'starNNW_file': os.path.join(config_dir, config_sex.get('starNNW_file')),
                        'checkimage_type': config_sex.get('checkimage_type'),
                        'clean_up_level': config_sex.getint('clean_up_level')}

        ### legitimate check
        if not shutil.which(sex_configs['cmd']):
            tmp = sex_configs['cmd']
            raise Exception(f"{tmp} is not an executable path for SExtractor!")
        for tmp in [sex_configs['config_file'], sex_configs['param_file'], sex_configs['filter_file'], sex_configs['starNNW_file']]:
            if not os.path.isfile(tmp):
                raise Exception(f"SExtractor-related file {tmp} not found!")

        ### cross match
        if sex_configs['cross_match']:
            config_cross = config['CrossMatch']
            sex_configs['mag_faint_cut'] = config_cross.getfloat('mag_faint_cut')
            sex_configs['save_matched'] = config_cross.getboolean('save_matched')
            sex_configs['save_false'] = config_cross.getboolean('save_false')
            sex_configs['save_missed'] = config_cross.getboolean('save_missed')
            sex_configs['mag_closest'] = config_cross.getboolean('mag_closest')
            sex_configs['r_max'] = config_cross.getfloat('r_max')

        ### collect
        configs_dict['sex'] = sex_configs

    ## === measure photometry
    if ('4' in taskIDs) or ('all' in taskIDs):
        config_MP = config['MeasurePhotometry']
        MP_method = config_MP.get('method')
        MP_configs = {'method': MP_method,
                        'detection_band': config_MP.get('detection_band'),
                        'bands': [x.strip() for x in config_MP.get('band_list').split(',')],
                        'image_label_list': [x.strip() for x in config_MP.get('image_label_list').split(',')]}

        ### limiting magnitude dictionary
        tmp_val = [float(x.strip()) for x in config_MP.get('band_1sigma_limits').split(',')]
        MP_configs['band_1sigma_limits'] = dict(zip(MP_configs['bands'], tmp_val))

        ### GAaP
        if MP_method.lower() == 'gaap':
            config_gaap = config['GAaP']
            MP_configs['gaap_dir'] = config_gaap.get('gaap_dir')
            MP_configs['min_aper'] = config_gaap.getfloat('min_aper')
            MP_configs['max_aper'] = config_gaap.getfloat('max_aper')
            MP_configs['use_PSF_map'] = config_gaap.getboolean('use_PSF_map')
            MP_configs['star_mag_cut'] = [float(i_m.strip()) for i_m in config_gaap.get('star_mag_cut').split(',')]
            MP_configs['clean_up_level'] = config_gaap.getint('clean_up_level')

            ### legitimate check
            if not os.path.isdir(MP_configs['gaap_dir']):
                tmp = MP_configs['gaap_dir']
                raise Exception(f"GAaP source directory {tmp} not found!")

        else:
            raise Exception(f'Unsupported photometry measurement method {MP_method}!')

        ### collect
        configs_dict['MP'] = MP_configs

    ## === measure photo-z
    if ('5' in taskIDs) or ('all' in taskIDs):
        config_MZ = config['MeasurePhotoz']
        MZ_method = config_MZ.get('method')
        MZ_configs = {'method': MZ_method}

        ### BPZ
        if MZ_method.lower() == 'bpz':
            config_bpz = config['BPZ']
            MZ_configs['BPZ_dir'] = config_bpz.get('BPZ_dir')
            MZ_configs['python2_cmd'] = config_bpz.get('python2_cmd')
            MZ_configs['bands'] = [x.strip() for x in config_bpz.get('band_list').split(',')]
            basename_tmp = config_bpz.get('band_CataNameBase')
            MZ_configs['bands_CataName'] = [f'{basename_tmp}_{band}' for band in MZ_configs['bands']]
            basename_tmp = config_bpz.get('banderr_CataNameBase')
            MZ_configs['banderrs_CataName'] = [f'{basename_tmp}_{band}' for band in MZ_configs['bands']]
            MZ_configs['bands_FilterName'] = [x.strip() for x in config_bpz.get('band_FilterName_list').split(',')]
            MZ_configs['photo_sys'] = config_bpz.get('photo_sys')
            MZ_configs['prior_band'] = config_bpz.get('prior_band')
            MZ_configs['prior_name'] = config_bpz.get('prior_name')
            MZ_configs['templates_name'] = config_bpz.get('templates_name')
            MZ_configs['interpolation'] = config_bpz.get('interpolation')
            MZ_configs['lkl_zmin'] = config_bpz.get('lkl_zmin')
            MZ_configs['lkl_zmax'] = config_bpz.get('lkl_zmax')
            MZ_configs['lkl_dz'] = config_bpz.get('lkl_dz')
            MZ_configs['lkl_odds'] = config_bpz.get('lkl_odds')
            MZ_configs['lkl_min_rms'] = config_bpz.get('lkl_min_rms')

            ### legitimate check
            if not os.path.isdir(MZ_configs['BPZ_dir']):
                tmp = MZ_configs['BPZ_dir']
                raise Exception(f"BPZ dir {tmp} not found!")
            if not shutil.which(MZ_configs['python2_cmd']):
                tmp = MZ_configs['python2_cmd']
                raise Exception(f"{tmp} is not an executable path for python2!")

        else:
            raise Exception(f'Unsupported photo-z measurement method {MZ_method}!')


        ### collect
        configs_dict['MZ'] = MZ_configs

    ## === measure galaxy shapes
    if ('6' in taskIDs) or ('all' in taskIDs):
        config_MS = config['MeasureShape']
        MS_method = config_MS.get('method')
        MS_configs = {'method': MS_method,
                        'detection_band': config_MS.get('detection_band'),
                        'bands': [x.strip() for x in config_MS.get('band_list').split(',')],
                        'image_label_list': [x.strip() for x in config_MS.get('image_label_list').split(',')]}

        if MS_method.lower() == 'lensfit':
            config_lensfit = config['lensfit']
            MS_configs['lensfit_dir'] = config_lensfit.get('lensfit_dir')
            MS_configs['clean_up_level'] = config_lensfit.getint('clean_up_level')
            MS_configs['postage_size'] = config_lensfit.get('postage_size')
            MS_configs['start_exposure'] = config_lensfit.get('start_exposure')
            MS_configs['end_exposure'] = config_lensfit.get('end_exposure')
            MS_configs['start_mag'] = config_lensfit.get('start_mag')
            MS_configs['end_mag'] = config_lensfit.get('end_mag')
            MS_configs['PSF_OVERSAMPLING'] = config_lensfit.get('PSF_OVERSAMPLING')
            MS_configs['PECUT'] = config_lensfit.get('PECUT')
            MS_configs['PRCUT'] = config_lensfit.get('PRCUT')
            MS_configs['LCUT'] = config_lensfit.get('LCUT')
            MS_configs['CAMERA'] = config_lensfit.get('CAMERA').upper()

            ### cores
            lensfit_cores = config_lensfit.get('lensfit_cores')
            if lensfit_cores is not None:
                MS_configs['lensfit_cores'] = int(lensfit_cores)
            else:
                MS_configs['lensfit_cores'] = 12

            ### legitimate check
            if not os.path.isdir(MS_configs['lensfit_dir']):
                tmp = MS_configs['lensfit_dir']
                raise Exception(f"lensfit dir {tmp} not found!")

        else:
            raise Exception(f'Unsupported shape measurement method {MS_method}!')

        ### collect
        configs_dict['MS'] = MS_configs

    ## === combine catalogue
    if ('7' in taskIDs) or ('all' in taskIDs):
        config_CC = config['CombineCata']
        CC_configs = {'format': config_CC.get('file_format').lower(),
                    'clean_up_level': config_CC.getint('clean_up_level')}

        ### legitimate check
        if CC_configs['format'] not in ['feather', 'fits', 'csv']:
            tmp = CC_configs['format']
            raise Exception(f"Unsupported file format {tmp}!]\n\
            supported file formats: fits, feather, csv")

        ### collect
        configs_dict['CC'] = CC_configs

    return configs_dict

def GenerateExampleConfig(file_name, user_name, date_time, __version__):
    """
    An example configuration file
    """

    config = f"# Example configuration file \n\
#   for {__version__} \n\
# Created by {user_name} ({date_time})\n\
# Note: If you do not use certain sections, just leave them whatever they are.\n\
#         They will not be loaded.\n\
\n\n\
################################## Paths ################################################\n\
[Paths]\n\n\
config_dir =  ../config/                       # directory to all the configuration files\n\
out_dir =                                      # main directory for all the outputs\n\
\n\n\
################################## GalInfo ##############################################\n\
[GalInfo]\n\n\
cata_file =                                    # input galaxy mock catalogue\n\
                                               # supported file types: feather, csv, fits\n\
position_type =         true                   # position to be used\n\
                                               #    true (use positions from the input catalogue)\n\
                                               #    grid (put in a grid)\n\
                                               #    random (put in random)\n\
grid_size =             18.                    # (arcsec) box size for grid\n\
                                               # default: 18 arcsec\n\
                                               # not used in other position_type\n\
mag_cut =               16, 27                 # magnitude cut for galaxies to be simulated\n\
Re_cut =                0, 99                  # size cut for galaxies to be simulated (arcsec)\n\
# catalogue column names to the desired info\n\
id_name =               index                  # unique galaxy id\n\
detection_mag_name =    r                      # correspond to the `detection_band` in [ImSim] \n\
mag_name_list =         u, g, r, i, Z, Y, J, H, Ks \n\
                                               # correspond to the the `band_list` in [ImSim]\n\
RaDec_names =           ra, dec \n\
shape_names =           none, none, none, none, none, none, none, none, none, none\n\
                                               # order: \n\
                                               #    Re (in arcsec!), sersic_n, axis_ratio, PA, \n\
                                               #    bulge_fraction, bulge_Re (in arcsec!), bulge_axis_ratio, bulge_sersic_n, \n\
                                               #    disk_Re (in arcsec!), disk_axis_ratio\n\
                                               # not all required, for those missed, simply feed none\n\
z_name =                zobs                   # column name for redshift\n\
                                               #    not used, only for saving\n\
\n\n\
################################## StarInfo ##############################################\n\
[StarInfo]\n\n\
cata_file =                                    # input star mock catalogue\n\
                                               # supported file types: feather, csv, fits\n\
                                               # leave it as empty if stars are not needed\n\
cata_area =             1                      # square degrees\n\
                                               # sky area spanned by the input catalogue\n\
                                               # should be equal or larger than the tile area\n\
position_type =         random                 # position to be used\n\
                                               #    random (randomly place the stars)\n\
                                               #    true (use positions from the input catalogue)\n\
mag_cut =               14, 27                 # magnitude cut for stars to be simulated\n\
# column names to the desired info\n\
id_name =               index                  # unique star id\n\
detection_mag_name =    r                      # correspond to the `detection_band` in [ImSim] \n\
mag_name_list =         u, g, r, i, Z, Y, J, H, Ks \n\
                                               # correspond to the the `band_list` in [ImSim]\n\
RaDec_names =           ra, dec                # not required, if stars are randomly placed\n\
\n\n\
################################## NoiseInfo ##############################################\n\
[NoiseInfo]\n\n\
cata_file =                                    # input noise background & psf catalogue\n\
                                               # supported file types: feather, csv, fits\n\
                                               # NOTE: tiles are orderly selected\n\
file4varChips =                                # a separate psf info catalogue for varChips mode (see ImSim)\n\
                                               # not required for other modes\n\
label_basename =        label                  # column name for label\n\
noise_basenames =       rms                    # base names for noise background\n\
                                               # order: rms\n\
                                               # the real column name is associated with band labels as `rms_r` etc\n\
psf_basenames =         seeing, beta, psf_e1, psf_e2\n\
                                               # base names for psf profile\n\
                                               # order:\n\
                                               #    (Moffat): seeing, MoffatBeta, e1, e2\n\
                                               #    (Airy): lam (in nanometre), diam (in metre), obscuration, e1, e2\n\
                                               # the real column name is associated with band labels\n\
                                               # not all required, for those missed, simply ignore or feed none\n\
id_basenames =          chip_id, expo_id       # column names for IDs used by file4varChips\n\
                                               # not used otherwise\n\
                                               # order: chip_id, expo_id\n\
psf_type =              Moffat                 # PSF profile, supported types:\n\
                                               #    Moffat, Airy\n\
\n\n\
################################## ImSim ###################################################\n\
[ImSim]\n\n\
survey =              KiDS                     # survey being simulated\n\
                                               # current supported surveys:\n\
                                               #    one_tile: simple one image including all the galaxies\n\
                                               #    simple_Nsqdeg: N can be any float corresponding to the tile sky area\n\
                                               #    KiDS: KiDS-like images (5 exposures, 32 chips, dither patterns and chip gaps)\n\
band_list =             u, g, r, i, Z, Y, J, H, Ks\n\
                                               # bands being simulated\n\
pixel_scale_list =      0.214, 0.214, 0.214, 0.214, 0.34, 0.34, 0.34, 0.34, 0.34\n\
                                               # pixel scale for each band image\n\
image_type_list =      sameExpo, sameExpo, varChips, sameExpo, simple, simple, simple, simple, simple\n\
                                               # image type for each band\n\
                                               # current supported types:\n\
                                               #    simple: without any survey feature (a simple stacked image)\n\
                                               #    sameExpo: different exposures use the same noise&PSF info\n\
                                               #    diffExpo: different exposures use different noise&PSF info\n\
                                               #    varChips: different chips use different PSF (on the top of diffExpo)\n\
image_chips =           False, False, True, False, False, False, False, False, False\n\
                                               # save individual chips or not\n\
                                               # required by lensfit\n\
image_PSF =             False, False, True, False, False, False, False, False, False\n\
                                               # save individual psf\n\
                                               # required by lensfit\n\
image_PSF_size =        48                     # (pixels) the size of the saved PSF image\n\
                                               #    it is assumed to be a square\n\
                                               #    default: 48*48 \n\
N_tiles =               1                      # number of tiles to be simulated\n\
                                               # make sure the NoiseInfo cata covers more than this requirement\n\
                                               # GalInfo can cover less than this requirement,\n\
                                               #    in which case repeating patterns will be produced\n\
                                               # NOTE: the total output tiles = N_tiles * N_rotations (specified below)\n\
gal_rotation_angles =   0, 90                      # degrees (put more values separated with ',' if needed)\n\
PSF_map =               False, False, False, False, False, False, False, False, False\n\
                                               # output the corresponding PSF map or not\n\
                                               # can be used by GAaP, but not mandatory if stars are simulated\n\
mag_zero =              30                     # simulated magnitude zero point\n\
\n\n\
################################## SWarp ###################################################\n\
[SWarp]\n\n\
# for coadding or resampling\n\
cmd =                   swarp                  # the executable path to the SWarp code\n\
config_files =          coadd_theli.swarp, coadd_aw.swarp\n\
                                               # SWarp configuration files\n\
                                               # more than one files are supported\n\
                                               #    in which case, more than one treatments are applied\n\
bands_group =           [r], [u, g, r, i]      # bands to be swarped\n\
                                               # NOTE: the group corresponding to the same config should be surrounded by `[]`\n\
image_label_list =      THELI, AW              # name to label the swaped results, one to each group\n\
only_resamples =        False, False           # set it True if only resampling but not coadding\n\
clean_up_levels =       0, 0                   # clean up level\n\
                                               #    0: none\n\
                                               #    1: original images\n\
                                               # NOTE: careful about cleaning before other swarp applied\n\
\n\n\
################################## SExtractor #################################################\n\
[SExtractor]\n\n\
# for detection\n\
cmd =                   sex                    # the executable path to the SExtractor code\n\
detection_band =        r                      # band for detection\n\
image_label =           THELI                  # label for the image type, can be either:\n\
                                               #    original (for original simulated images)\n\
                                               #    any label specified in `swarp_labels`\n\
pixel_scale =           0.214                  # pixel scale for the image\n\
cross_match =           True                   # cross-match with the input catalogue\n\
                                               #    in which case, catalogues with match info will be saved\n\
                                               #    see next section for configuration\n\
config_file =           kids_sims.sex          # SExtractor configuration file\n\
param_file =            sex_image.param        # SExtractor parameter file\n\
filter_file =           default.conv           # SExtractor filter file\n\
starNNW_file =          default.nnw            # SExtractor Neural-Network_Weight table file\n\
checkimage_type =       NONE                   # can be one of \n\
                                               # NONE, BACKGROUND, MINIBACKGROUND, OBJECTS, SEGMENTATION, APERTURES, FILTERED\n\
clean_up_level =        0                      # clean up level\n\
                                               #    0: none\n\
                                               #    1: original images\n\
                                               #    2: and .sex files\n\
\n\n\
################################## CrossMatch #################################################\n\
[CrossMatch]\n\n\
mag_faint_cut =         26                     # faintest sources can be possible detected\n\
                                               # for the sake of speed-up\n\
save_matched =          True                   # save the matched object info\n\
save_false =            True                   # save false-detected object info\n\
save_missed =           True                   # save the missed object info\n\
mag_closest =           True                   # use magnitude to select for duplicated match\n\
r_max =                 0.6                    # (arcsec) allowed maximum separation\n\
\n\n\
################################## MeasurePhotometry ########################################################\n\
[MeasurePhotometry]\n\n\
method =                GAaP                   # method for photometry measurement\n\
                                               # supported method:\n\
                                               #    GAaP\n\
detection_band =        r                      # band with detection catalogue\n\
band_list =             u, g, r, i, Z, Y, J, H, Ks\n\
                                               # bands being measured\n\
band_1sigma_limits =    25.5, 26.3, 26.2, 24.9, 24.85, 24.1, 24.2, 23.3, 23.2\n\
                                               # 1 sigma limiting magnitude\n\
                                               # used for unmeasured objects' mag error\n\
                                               # default value from KV450 median\n\
image_label_list =      AW, AW, AW, AW, original, original, original, original, original\n\
                                               # a list of labels for the image types, can be either:\n\
                                               #    original (for original simulated images)\n\
                                               #    any label specified in `swarp_labels`\n\
\n\n\
[GAaP]\n\n\
gaap_dir =                                     # directory containing GAaP bins and libs\n\
min_aper =              0.7                    # minimum aperture size\n\
max_aper =              2.0                    # maximum aperture size\n\
use_PSF_map =           False                  # use separate psf map for psf estimation\n\
                                               # only required if stars are not simulated\n\
star_mag_cut =          16, 20                 # magnitude range for stars used for PSF estimation\n\
                                               # not used if PSF map is provided\n\
clean_up_level =        0                      # clean up level\n\
                                               #    0: none\n\
                                               #    1: tmp directory\n\
                                               #    2: and *.gaap files\n\
\n\n\
################################## MeasurePhotoz ########################################################\n\
[MeasurePhotoz]\n\n\
method =                BPZ                    # method for photo-z measurement\n\
                                               # supported method:\n\
                                               #    BPZ\n\
\n\n\
[BPZ]\n\n\
BPZ_dir =                                      # the directory containing BPZ-related codes\n\
python2_cmd =                                  # the executable path to the python2\n\
band_list =              u, g, r, i, Z, Y, J, H, Ks\n\
band_CataNameBase =      MAG_GAAP_0p7\n\
banderr_CataNameBase =   MAGERR_GAAP_0p7\n\
band_FilterName_list =   Paranal_OmegaCAM.u_SDSS, Paranal_OmegaCAM.g_SDSS, Paranal_OmegaCAM.r_SDSS, Paranal_OmegaCAM.i_SDSS, LSST_LSST.z, LSST_LSST.y, 2MASS_2MASS.J, 2MASS_2MASS.H, 2MASS_2MASS.Ks\n\
photo_sys =              AB\n\
prior_band =             i\n\
prior_name =             NGVS\n\
templates_name =         CWWSB_capak\n\
interpolation =          10\n\
lkl_zmin =               0.0002\n\
lkl_zmax =               2.2984\n\
lkl_dz =                 0.01\n\
lkl_odds =               0.68\n\
lkl_min_rms =            0.0\n\
\n\n\
################################## MeasureShape #################################################\n\
[MeasureShape]\n\n\
method =               lensfit                 # method for galaxy shape measurement\n\
                                               # supported method:\n\
                                               #    lensfit\n\
detection_band =        r                      # band with detection catalogue\n\
band_list =             r\n\
                                               # bands being measured\n\
image_label_list =      original\n\
                                               # a list of labels for the image types, can be either:\n\
                                               #    original (for original simulated images)\n\
                                               #    any label specified in `swarp_labels`\n\
                                               # NOTE: lensfit uses individual exposures, so it is `original`\n\
\n\n\
[lensfit]\n\n\
lensfit_dir =                                  # directory containing lensfit code\n\
lensfit_cores =         12                     # number of cores used by each lensfit run\n\
                                               # should be consistent with that compiled in lensfit (flensfit_NT[lensfit_cores])\n\
postage_size =          48\n\
start_exposure =        1\n\
end_exposure =          5\n\
start_mag =             20.0\n\
end_mag =               25.0\n\
PSF_OVERSAMPLING =      1\n\
PECUT =                 0.02\n\
PRCUT =                 0.02\n\
LCUT =                  0.05\n\
CAMERA =                KIDS\n\
clean_up_level =        0                      # clean up level\n\
                                               #    0: none\n\
                                               #    1: tmp directory\n\
\n\n\
################################## CombineCata ###################################################\n\
[CombineCata]\n\n\
file_format =              fits                # output file format\n\
                                               # supported formats:\n\
                                               #    fits, feather, csv\n\
clean_up_level =        0                      # clean up level\n\
                                               #    0: none\n\
                                               #    1: all itermediate catalogues (that is ALL, so be very careful!)\n\
"
    # write out the example config file
    with open(file_name, 'w') as configfile:
        configfile.write(config)
    logger.info(f'An example configuration file `{file_name}` is generated in the current directory.')
