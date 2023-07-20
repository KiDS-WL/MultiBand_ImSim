# -*- coding: utf-8 -*-
# @Author: lshuns
# @Date:   2021-02-03, 15:58:35
# @Last Modified by:   lshuns
# @Last Modified time: 2023-07-20 16:07:22

### module to generate an example configuration file

import os
import re
import shutil
import logging
import pathlib
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

    ### main 
    out_dir = config_paths.get('out_dir')
    out_dir = os.path.join(out_dir, run_tag)

    ### images
    ima_dir = os.path.join(out_dir, 'images')
    pathlib.Path(ima_dir).mkdir(parents=True, exist_ok=True)
    logger.info(f'images will be saved to {ima_dir}')

    ### catalogues
    cata_dir = os.path.join(out_dir, 'catalogues')
    pathlib.Path(cata_dir).mkdir(parents=True, exist_ok=True)
    logger.info(f'catalogues will be saved to {cata_dir}')

    ### running_log
    if running_log:
        log_dir = os.path.join(out_dir, 'running_log')
        pathlib.Path(log_dir).mkdir(parents=True, exist_ok=True)
        logger.info(f'Running log from external codes will be saved to {log_dir}')
    else:
        log_dir = None

    ### tmp dir
    tmp_dir = config_paths.get('tmp_dir')
    if tmp_dir is not None:
        tmp_dir = os.path.join(tmp_dir, run_tag, 'tmp')
    else:
        tmp_dir = os.path.join(out_dir, 'tmp')
    pathlib.Path(tmp_dir).mkdir(parents=True, exist_ok=True)
    logger.info(f'tmp files will be saved to {tmp_dir}')

    ### folder names
    cata_folder_names = config_paths.get('cata_folder_names')
    if cata_folder_names is not None:
        cata_folder_names = [x.strip() for x in cata_folder_names.split(',')]
    else:
        cata_folder_names = ['input', 'SExtractor', 'CrossMatch', 'photometry', 'photo_z', 'shapes', 'combined']

    ### combine to a dictionary
    work_dirs = {'main': out_dir,
                    'ima': ima_dir,
                    'cata': cata_dir,
                    'log': log_dir, 
                    'tmp': tmp_dir,
                    'cata_folder_names': cata_folder_names}
    del out_dir, ima_dir, cata_dir, log_dir, tmp_dir, cata_folder_names

    ## === dictionary for collecting all config info
    configs_dict = {'work_dirs': work_dirs}

    ## === b. ImSim
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

    ### detection band 
    detection_band = config_imsim.get('detection_band')
    if detection_band:
        imsim_configs['detection_band'] = detection_band
    else:
        ## set to be the first band in the list
        imsim_configs['detection_band'] = imsim_configs['bands'][0]

    ### psf image size 
    psf_size = config_imsim.get('image_PSF_size')
    if psf_size:
        imsim_configs['image_PSF_size'] = float(psf_size)
    else:
        imsim_configs['image_PSF_size'] = 48 #pixels

    ### for casual mode
    casual_mag = config_imsim.get('casual_mag')
    if casual_mag:
        imsim_configs['casual_mag'] = float(casual_mag)
        imsim_configs['casual_band'] = config_imsim.get('casual_band')
        imsim_configs['casual_Nbins'] = int(config_imsim.get('casual_Nbins'))
        imsim_configs['casual_FracSeedGal'] = float(config_imsim.get('casual_FracSeedGal'))
    else:
        imsim_configs['casual_mag'] = 99.

    ### how to calculate the sky area
    simple_area = config_imsim.getboolean('simple_area')
    if simple_area is None:
        ## old code is using the simple way
        imsim_configs['simple_area'] = True
    else:
        imsim_configs['simple_area'] = simple_area

    ### collect
    configs_dict['imsim'] = imsim_configs

    ## === c. noise info
    config_noise = config['NoiseInfo']
    noise_configs = {'file': config_noise.get('cata_file'),
                    'file4varChips': config_noise.get('file4varChips'),
                    'psf_PixelIma_dir': config_noise.get('psf_PixelIma_dir')}

    ### check existence
    if not os.path.isfile(noise_configs['file']):
        tmp = noise_configs['file']
        raise Exception(f"noise file {tmp} not found!")
    if (noise_configs['file4varChips']) and (not os.path.isfile(noise_configs['file4varChips'])):
        tmp = noise_configs['file4varChips']
        raise Exception(f"separate psf file {tmp} not found!")
    if (noise_configs['psf_PixelIma_dir']) and (not os.path.isdir(noise_configs['psf_PixelIma_dir'])):
        tmp = noise_configs['psf_PixelIma_dir']
        raise Exception(f"psf image dir {tmp} not found!")

    ### psf profile type 
    try:
        psf_type_list = [x.strip().lower() for x in config_noise.get('psf_type_list').split(',')]
    except AttributeError:
        logger.warning('Using old psf_type parameter, please update ASAP!')
        psf_type = config_noise.get('psf_type')
        if not psf_type:
            psf_type = 'moffat'
        psf_type_list = [psf_type.lower()] 
    if len(psf_type_list) == 1:
        psf_type_list *= len(imsim_configs['bands'])
    noise_configs['psf_type_list'] = psf_type_list
    ###### supported types
    for psf_type in noise_configs['psf_type_list']:
        if psf_type not in ['moffat', 'airy', 'pixelima']:
            raise Exception(f'not supported psf_type: {psf_type}')

    ### >>> for old version                     
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
        if (noise_configs['file4varChips']):
            noise_configs['id_basenames'] = [x.strip() for x in config_noise.get('id_basenames').split(',')]

        # >>> for old version
        psf_basenames = config_noise.get('psf_basenames')
        if psf_basenames:
            logger.warning('Using old psf_basenames parameter, please update ASAP!')
            noise_configs['psf_basenames_moffat'] = [x.strip() for x in psf_basenames.split(',')]
            noise_configs['psf_basenames_airy'] = [x.strip() for x in psf_basenames.split(',')]
        else:
            try:
                noise_configs['psf_basenames_moffat'] = [x.strip() for x in config_noise.get('psf_basenames_moffat').split(',')]
            except AttributeError:    
                noise_configs['psf_basenames_moffat'] = None

            try:
                noise_configs['psf_basenames_airy'] = [x.strip() for x in config_noise.get('psf_basenames_airy').split(',')]
            except AttributeError:    
                noise_configs['psf_basenames_airy'] = None

    ### collect
    configs_dict['noise'] = noise_configs

    # ++++++++++++++ 2. some specific info (required only by certain tasks)

    ## === galaxy & star info
    if ('1' in taskIDs) or ('all' in taskIDs):
        ## == galaxy
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
        if grid_size:
            gal_configs['grid_size'] = float(grid_size)
        else:
            gal_configs['grid_size'] = 18 #arcsec

        ### check existence
        if not os.path.isfile(gal_configs['file']):
            tmp = gal_configs['file']
            raise Exception(f"galaxies file {tmp} not found!")

        ### collect
        configs_dict['gal'] = gal_configs

        ## == star
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

        ### collect
        configs_dict['star'] = star_configs

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

            ### match with TAN projection
            use_TAN = config_cross.getboolean('use_TAN')
            if use_TAN is None:
                ## old code is using the sky coordinates
                sex_configs['use_TAN'] = False
            else:
                sex_configs['use_TAN'] = use_TAN

            sex_configs['r_max_pixel'] = config_cross.getfloat('r_max_pixel')

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

        ### spatial variation dictionary
        tmp_val = [(x.strip().lower() in ['true', 'yes', '1']) for x in config_MP.get('band_spatial_variation').split(',')]
        MP_configs['band_spatial_variation'] = dict(zip(MP_configs['bands'], tmp_val))

        ### GAaP
        if MP_method.lower() == 'gaap':
            config_gaap = config['GAaP']
            MP_configs['gaap_dir'] = config_gaap.get('gaap_dir')
            MP_configs['min_aper_list'] = [float(tmp.strip()) for tmp in config_gaap.get('min_aper_list').split(',')]
            MP_configs['max_aper'] = config_gaap.getfloat('max_aper')
            MP_configs['use_PSF_map'] = config_gaap.getboolean('use_PSF_map')
            MP_configs['star_SNR_cut'] = [float(i_m.strip()) for i_m in config_gaap.get('star_SNR_cut').split(',')]
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
            MZ_configs['detection_band'] = config_bpz.get('detection_band')
            MZ_configs['bands'] = [x.strip() for x in config_bpz.get('band_list').split(',')]
            MZ_configs['bands_FilterName'] = [x.strip() for x in config_bpz.get('band_FilterName_list').split(',')]
            ### some basename for columns
            MZ_configs['band_CataNameBase'] = config_bpz.get('band_CataNameBase')
            MZ_configs['banderr_CataNameBase'] = config_bpz.get('banderr_CataNameBase')
            MZ_configs['bandflag_CataNameBase'] = config_bpz.get('bandflag_CataNameBase')
            MZ_configs['bandlim_CataNameBase'] = config_bpz.get('bandlim_CataNameBase')
            ### some setup for BPZ
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

            ### clean up
            clean_tmp = config_bpz.get('clean_up_level')
            if clean_tmp is not None:
                MZ_configs['clean_up_level'] = int(clean_tmp)
            else:
                MZ_configs['clean_up_level'] = 0

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

    ## === PSF modelling
    if ('6_1' in taskIDs) or ('all' in taskIDs):
        config_PSF = config['PSFmodelling']
        PSF_method = config_PSF.get('method')
        PSF_configs = {'method': PSF_method,
                        'detection_band': config_PSF.get('detection_band'),
                        'bands': [x.strip() for x in config_PSF.get('band_list').split(',')],
                        'image_label_list': [x.strip() for x in config_PSF.get('image_label_list').split(',')]}

        ### folder prefix
        folder_prefix = config_PSF.get('folder_prefix')
        if folder_prefix is not None:
            PSF_configs['folder_prefix'] = folder_prefix
        else:
            PSF_configs['folder_prefix'] = 'psf_coeff'

        if PSF_method.lower() == 'ima2coeffs':
            config_tmp = config['ima2coeffs']
            PSF_configs['ima2coeffs_dir'] = config_tmp.get('ima2coeffs_dir')

            ### legitimate check
            if not os.path.isdir(PSF_configs['ima2coeffs_dir']):
                tmp = PSF_configs['ima2coeffs_dir']
                raise Exception(f"ima2coeffs dir {tmp} not found!")
        
        elif PSF_method.lower() == 'makeglobalpsf':
            config_makeglobalpsf = config['makeglobalpsf']
            PSF_configs['makeglobalpsf_dir'] = config_makeglobalpsf.get('makeglobalpsf_dir')
            PSF_configs['clean_intermediate'] = config_makeglobalpsf.getboolean('clean_intermediate')
            PSF_configs['global_order'] = config_makeglobalpsf.get('global_order')
            PSF_configs['chip_order'] = config_makeglobalpsf.get('chip_order')
            PSF_configs['snratio'] = config_makeglobalpsf.get('snratio')
            PSF_configs['start_mag'] = config_makeglobalpsf.get('start_mag')
            PSF_configs['end_mag'] = config_makeglobalpsf.get('end_mag')
            PSF_configs['CAMERA'] = config_makeglobalpsf.get('CAMERA').upper()
            PSF_configs['SWARP_CONFIG'] = config_makeglobalpsf.get('SWARP_CONFIG')

            ### legitimate check
            if not os.path.isdir(PSF_configs['makeglobalpsf_dir']):
                tmp = PSF_configs['makeglobalpsf_dir']
                raise Exception(f"makeglobalpsf dir {tmp} not found!")

        else:
            raise Exception(f'Unsupported shape measurement method {MS_method}!')

        ### collect
        configs_dict['PSFmodelling'] = PSF_configs

    ## === measure galaxy shapes
    if ('6_2' in taskIDs) or ('all' in taskIDs):
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

            ### PSF folder
            psf_prefix = config_lensfit.get('psf_folder_prefix')
            if psf_prefix is not None:
                MS_configs['psf_folder_prefix'] = psf_prefix
            else:
                MS_configs['psf_folder_prefix'] = 'psf_coeff'

            ### cores
            lensfit_cores = config_lensfit.get('lensfit_cores')
            if lensfit_cores is not None:
                MS_configs['lensfit_cores'] = int(lensfit_cores)
            else:
                MS_configs['lensfit_cores'] = 12

            ### lensfit type
            lensfit_type = config_lensfit.get('lensfit_type')
            if lensfit_type is not None:
                MS_configs['lensfit_type'] = lensfit_type
            else:
                MS_configs['lensfit_type'] = 'flensfit'

            ### lensfit version
            lensfit_version = config_lensfit.get('lensfit_version')
            if lensfit_version is not None:
                MS_configs['lensfit_version'] = lensfit_version
            else:
                MS_configs['lensfit_version'] = '309'

            ### lensfit memory
            MEMORY_LIMIT = config_lensfit.get('MEMORY_LIMIT')
            if MEMORY_LIMIT is not None:
                MS_configs['MEMORY_LIMIT'] = MEMORY_LIMIT
            else:
                MS_configs['MEMORY_LIMIT'] = '30000'

            ### legitimate check
            lensfit_run = os.path.join(MS_configs['lensfit_dir'], 'bin', \
                            MS_configs['lensfit_type']+'_NT'+str(MS_configs['lensfit_cores']))
            if not os.path.isfile(lensfit_run):
                raise Exception(f"lensfit code {lensfit_run} not found!")

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
config_dir =            your_dir_to_MultiBand_ImSim/config\n\
                                               # directory to all the configuration files\n\
out_dir =               find/somewhere/with/large/space\n\
                                               # main directory for all the final outputs\n\
tmp_dir =               find/somewhere/local/to/speedup\n\
                                               # tmp directory for tmp outputs\n\
cata_folder_names =     input, SExtractor, CrossMatch, photometry, photo_z, shapes, combined\n\
                                               # folder names for saving different catalogues\n\
                                               # order: \n\
                                               #    input, detection, CrossMatch, photometry, photo_z, shapes, combined_suffix\n\
\n\n\
################################## GalInfo ##############################################\n\
[GalInfo]\n\n\
cata_file =             your_dir_to_input_cata/skills_sth.feather\n\
                                               # input galaxy mock catalogue\n\
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
detection_mag_name =    r_SDSS_apparent_corr   # correspond to the `detection_band` in [ImSim] \n\
mag_name_list =         r_SDSS_apparent_corr, u_SDSS_apparent_corr, g_SDSS_apparent_corr, i_SDSS_apparent_corr, z_SDSS_apparent_corr, Y_VISTA_apparent_corr, J_VISTA_apparent_corr, H_VISTA_apparent_corr, K_VISTA_apparent_corr\n\
                                               # correspond to the the `band_list` in [ImSim]\n\
RaDec_names =           ra, dec \n\
shape_names =           Re_arcsec, shape/sersic_n, BA, PA_random, none, none, none, none, none, none\n\
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
cata_file =             your_dir_to_input_cata/trilegal_sth.feather\n\
                                               # input star mock catalogue\n\
                                               # supported file types: feather, csv, fits\n\
                                               # leave it as empty if stars are not needed\n\
cata_area =             10                     # square degrees\n\
                                               # sky area spanned by the input catalogue\n\
                                               # should be equal or larger than the tile area\n\
position_type =         random                 # position to be used\n\
                                               #    random (randomly place the stars)\n\
                                               #    true (use positions from the input catalogue)\n\
mag_cut =               14, 27                 # magnitude cut for stars to be simulated\n\
# column names to the desired info\n\
id_name =               index                  # unique star id\n\
detection_mag_name =    r                      # correspond to the `detection_band` in [ImSim] \n\
mag_name_list =         r, u, g, i, Z, Y, J, H, Ks \n\
                                               # correspond to the the `band_list` in [ImSim]\n\
RaDec_names =           ra, dec                # not required, if stars are randomly placed\n\
\n\n\
################################## NoiseInfo ##############################################\n\
[NoiseInfo]\n\n\
cata_file =             your_dir_to_MultiBand_ImSim/noise_info/skills_fiducial/noise_sth.csv\n\
                                               # input noise background & psf catalogue\n\
                                               # supported file types: feather, csv, fits\n\
                                               # NOTE: tiles are orderly selected\n\
file4varChips =         your_dir_to_MultiBand_ImSim/noise_info/kids_dr4_psf_moffat_fromcoeffs.csv\n\
                                               # a separate psf info catalogue for varChips mode (see ImSim)\n\
                                               # not required for other modes\n\
psf_PixelIma_dir =                             # directory contains psf images\n\
                                               # only required for psf_type==PixelIma\n\
label_basename =        label                  # column name for label\n\
noise_basenames =       rmsExpo                # base names for noise background\n\
                                               # order: rms\n\
                                               # the real column name is associated with band labels as `rms_r` etc\n\
psf_basenames_moffat =  InputSeeing, InputBeta, seeing_e1, seeing_e2\n\
                                               # base names for psf profile\n\
                                               # used by psf_type == Moffat\n\
                                               # order:\n\
                                               #    seeing, MoffatBeta, e1, e2\n\
                                               # the real column name is associated with band labels\n\
                                               # not all required, for those missed, simply ignore or feed none\n\
psf_basenames_airy =    lam, diam, obscuration, psf_e1, psf_e2\n\
                                               # base names for psf profile\n\
                                               # used by psf_type == Airy\n\
                                               # order:\n\
                                               #    lam (in nanometre), diam (in metre), obscuration, e1, e2\n\
                                               # the real column name is associated with band labels\n\
                                               # not all required, for those missed, simply ignore or feed none\n\
id_basenames =          chip_id, expo_id       # column names for IDs used by file4varChips\n\
                                               # not used otherwise\n\
                                               # order: chip_id, expo_id\n\
psf_type_list =         Moffat, Moffat, Moffat, Moffat, Moffat, Moffat, Moffat, Moffat, Moffat\n\
                                               # PSF profile list, same order as band_list\n\
                                               # supported types:\n\
                                               #    Moffat, Airy, PixelIma\n\
\n\n\
################################## ImSim ###################################################\n\
[ImSim]\n\n\
survey =                KiDS                   # survey being simulated\n\
                                               # current supported surveys:\n\
                                               #    one_tile: simple one image including all the galaxies\n\
                                               #    simple_Nsqdeg: N can be any float corresponding to the tile sky area\n\
                                               #    KiDS: KiDS-like images (5 exposures, 32 chips, dither patterns and chip gaps)\n\
detection_band =        r                      # band with detection catalogue\n\
                                               #    not necessary to be simulated, could be simply for selection\n\
band_list =             r, u, g, i, Z, Y, J, H, Ks\n\
                                               # bands being simulated\n\
pixel_scale_list =      0.214, 0.214, 0.214, 0.214, 0.34, 0.34, 0.34, 0.34, 0.34\n\
                                               # pixel scale for each band image\n\
image_type_list =       varChips, sameExpo, sameExpo, sameExpo, simple, simple, simple, simple, simple\n\
                                               # image type for each band\n\
                                               # current supported types:\n\
                                               #    simple: without any survey feature (a simple stacked image)\n\
                                               #    sameExpo: different exposures use the same noise&PSF info\n\
                                               #    diffExpo: different exposures use different noise&PSF info\n\
                                               #    varChips: different chips use different PSF (on the top of diffExpo)\n\
image_chips =           True, False, False, False, False, False, False, False, False\n\
                                               # save individual chips or not\n\
                                               # required by lensfit\n\
image_PSF =             True, False, False, False, False, False, False, False, False\n\
                                               # save individual psf\n\
                                               # required by lensfit\n\
image_PSF_size =        48                     # (pixels) the size of the saved PSF image\n\
                                               #    it is assumed to be a square\n\
                                               #    default: 48*48 \n\
casual_mag =            25                     # up to which magnitude galaxies are casually simulated\n\
                                               #    means: only drawImage for a few, \n\
                                               #        the rest is randomly sampled from those simulated\n\
                                               # should be lower than mag_cut[1] in [GalInfo]\n\
                                               #    otherwise, ignored\n\
casual_band =           r                      # band used by casual_mag\n\
casual_Nbins =          20                     # how many quantile bins used by casual simulation\n\
casual_FracSeedGal =    0.01                   # Fraction of seed galaxies\n\
                                               # seed galaxies: those simulated carefully\n\
N_tiles =               1                      # number of tiles to be simulated\n\
                                               # make sure the NoiseInfo cata covers more than this requirement\n\
                                               # GalInfo can cover less than this requirement,\n\
                                               #    in which case repeating patterns will be produced\n\
                                               # NOTE: the total output tiles = N_tiles * N_rotations (specified below)\n\
gal_rotation_angles =   0, 90                      # degrees (put more values separated with ',' if needed)\n\
PSF_map =               False\n\
                                               # output the corresponding PSF map or not\n\
                                               # can be used by GAaP, but not mandatory if stars are simulated\n\
mag_zero =              30                     # simulated magnitude zero point\n\
simple_area =           False                  # calculate the sky area using \n\
                                               # simple Euclidean geometry (True)\n\
                                               # proper Spherical geometry (False), recommended for |dec|>5\n\
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
                                               #    1: .sex files\n\
\n\n\
################################## CrossMatch #################################################\n\
[CrossMatch]\n\n\
mag_faint_cut =         26                     # faintest sources can be possible detected\n\
                                               # for the sake of speed-up\n\
save_matched =          True                   # save the matched object info\n\
save_false =            False                  # save false-detected object info\n\
save_missed =           False                  # save the missed object info\n\
mag_closest =           True                   # use magnitude to select for duplicated match\n\
r_max =                 0.6                    # (arcsec) allowed maximum separation\n\
                                               # used when use_TAN==False\n\
use_TAN =               True                   # using TAN projection before CrossMatch\n\
r_max_pixel =           2.5                    # (pixel) allowed maximum separation\n\
                                               # used when use_TAN==True\n\
\n\n\
################################## MeasurePhotometry ########################################################\n\
[MeasurePhotometry]\n\n\
method =                GAaP                   # method for photometry measurement\n\
                                               # supported method:\n\
                                               #    GAaP\n\
detection_band =        r                      # band with detection catalogue\n\
band_list =             r, u, g, i, Z, Y, J, H, Ks\n\
                                               # bands being measured\n\
band_spatial_variation = True, False, False, False, False, False, False, False, False\n\
                                               # does psf vary spatially\n\
image_label_list =      AW, AW, AW, AW, original, original, original, original, original\n\
                                               # a list of labels for the image types, can be either:\n\
                                               #    original (for original simulated images)\n\
                                               #    any label specified in `swarp_labels`\n\
\n\n\
[GAaP]\n\n\
gaap_dir =              your_dir_to_gaap       # directory containing GAaP bins and libs\n\
min_aper_list =         0.7, 1.0               # minimum aperture size\n\
                                               # a list of sizes can be provided\n\
                                               # the final decision follows KiDS\n\
max_aper =              2.0                    # maximum aperture size\n\
use_PSF_map =           False                  # use separate psf map for psf estimation\n\
                                               # only required if stars are not simulated\n\
star_SNR_cut =          100, 2000              # SNR range for stars used for PSF estimation\n\
                                               # not used if PSF map is provided\n\
clean_up_level =        0                      # clean up level\n\
                                               #    0: none\n\
                                               #    1: tmp directory\n\
\n\n\
################################## MeasurePhotoz ########################################################\n\
[MeasurePhotoz]\n\n\
method =                BPZ                    # method for photo-z measurement\n\
                                               # supported method:\n\
                                               #    BPZ\n\
\n\n\
[BPZ]\n\n\
BPZ_dir =                your_dir_to_BPZ_code  # the directory containing BPZ-related codes\n\
python2_cmd =            python2               # the executable path to the python2\n\
                                               #    with necessary package for BPZ\n\
detection_band =         r                      # band with detection catalogue\n\
band_list =              u, g, r, i, Z, Y, J, H, Ks\n\
band_FilterName_list =   KiDSVIKING_u, KiDSVIKING_g, KiDSVIKING_r, KiDSVIKING_i, KiDSVIKING_Z2, KiDSVIKING_Y, KiDSVIKING_J, KiDSVIKING_H, KiDSVIKING_Ks\n\
band_CataNameBase =      MAG_GAAP\n\
banderr_CataNameBase =   MAGERR_GAAP\n\
bandflag_CataNameBase =  FLAG_GAAP\n\
bandlim_CataNameBase =   MAG_LIM\n\
photo_sys =              AB\n\
prior_band =             i\n\
prior_name =             NGVS\n\
templates_name =         CWWSB_capak\n\
interpolation =          10\n\
lkl_zmin =               0.01\n\
lkl_zmax =               2.5\n\
lkl_dz =                 0.01\n\
lkl_odds =               0.68\n\
lkl_min_rms =            0.067\n\
clean_up_level =         0                     # clean up level\n\
                                               #    0: none\n\
                                               #    1: tmp directory\n\
\n\n\
################################## PSFmodelling #################################################\n\
[PSFmodelling]\n\n\
method =                ima2coeffs             # method for PSF modelling\n\
                                               # supported method:\n\
                                               #    makeglobalpsf\n\
                                               #    ima2coeffs\n\
folder_prefix =         psf_coeff              # prefix of folders saving PSF models\n\
detection_band =        r                      # band with detection catalogue\n\
band_list =             r                      # bands being measured\n\
image_label_list =      original\n\
                                               # a list of labels for the image types, can be either:\n\
                                               #    original (for original simulated images)\n\
                                               #    any label specified in `swarp_labels`\n\
\n\n\
[ima2coeffs]\n\n\
ima2coeffs_dir =        your_dir_to_lensfit/utils\n\
                                               # directory containing psfimage2coeffs\n\
\n\n\
[makeglobalpsf]\n\n\
makeglobalpsf_dir = your_dir_to_lensfit/bin\n\
                                               # directory containing makeglobalpsf and globalshifts\n\
clean_intermediate =    true                   # clean log files or not\n\
# some makeglobalpsf-related values\n\
global_order =          4\n\
chip_order =            1\n\
snratio =               20\n\
start_mag =             18.0\n\
end_mag =               24.0\n\
CAMERA =                KIDS\n\
SWARP_CONFIG =          your_dir_to_lensfit/input_files\n\
                                               # config file for swarp used by lensfit\n\
\n\n\
################################## MeasureShape #################################################\n\
[MeasureShape]\n\n\
method =                lensfit                # method for galaxy shape measurement\n\
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
lensfit_dir =           your_dir_to_lensfit\n\
                                               # directory containing lensfit code\n\
lensfit_cores =         48                     # number of cores used by each lensfit run\n\
                                               # should be consistent with that compiled in lensfit (flensfit_NT[lensfit_cores])\n\
psf_folder_prefix =     psf_coeff              # prefix of prefix of folders saving PSF models\n\
# some lensfit-related values\n\
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
lensfit_type =          flensfit               # 309 using flensfit, 321 using lensfit\n\
lensfit_version =       309                    # current supported version: 309, 321\n\
MEMORY_LIMIT =          30000                  # only relevant for 321\n\
clean_up_level =        0                      # clean up level\n\
                                               #    0: none\n\
                                               #    1: tmp directory\n\
\n\n\
################################## CombineCata ###################################################\n\
[CombineCata]\n\n\
file_format =           fits                   # output file format\n\
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
