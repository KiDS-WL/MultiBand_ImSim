# -*- coding: utf-8 -*-
# @Author: lshuns
# @Date:   2020-12-21 11:44:14
# @Last Modified by:   lshuns
# @Last Modified time: 2021-10-31 16:14:51

### main module to run the whole pipeline

import re
import os
import sys
import time
import glob
import shutil
import logging
import datetime
import argparse

import numpy as np
import pandas as pd
import multiprocessing as mp

from astropy.table import Table

import BPZ
import GAaP
import ImSim
import LensFit
import LoadCata
import Astromatic
import CrossMatch

import RunConfigFile

__version__ = "MultiBand_ImSim v0.41"

# ++++++++++++++ parser for command-line interfaces
parser = argparse.ArgumentParser(
    description=f"{__version__}: generate realistic multi-band sky images using galaxy & star mock catalogues.",
    formatter_class=argparse.RawTextHelpFormatter)
parser.add_argument(
    "taskIDs", nargs="+", type=str, choices=['0', '1', '2', '3', '4', '5', '6', '7', 'all'], metavar='taskIDs',
    help="Select a set of IDs for processes:\n\
    0: generate an example configuration file\n\
    1: simulate images\n\
    2: swarp images\n\
    3: detect objects\n\
    4: measure photometry\n\
    5: measure photo-z\n\
    6: measure shapes\n\
    7: create a combined catalogue\n\
    all: run all tasks in a sequence")
parser.add_argument(
    "--runTag", type=str, default='test',
    help="A tag to label the current run.")
parser.add_argument(
    "--cosmic_shear", type=float, nargs=2, default=[0.00, 0.00],
    help="2D cosmic shear values: g1 g2. \n"
         "   (ignored if shear_columns is provided)")
parser.add_argument(
    "--shear_columns", type=str, nargs=2, default=None,
    help="column names to the input shear. \n"
         "   (use cosmic_shear instead if need constant shear)")
parser.add_argument(
    "-c", "--config", type=str, metavar='config_file',
    help="A configuration file including all necessary setups.\n"
            "   To generate an example file, use `python Run.py 0`.")
parser.add_argument(
    "--threads", type=int, default=int(mp.cpu_count()*0.5), metavar='number_threads',
    help="The maximum number of threads to use. \n"
        "   (default: %(default)s [half of the total CPU count])")
parser.add_argument(
    "--rng_seed", type=int, default=940120,
    help="base seed for all random generation.")
parser.add_argument(
    "--loglevel", type=str, default='INFO', metavar='logging_level',
    help="Logging level: DEBUG, INFO, WARNING, ERROR, CRITICAL")
parser.add_argument(
    "--sep_running_log", action="store_true",
    help="If set, save log files for external code running info.")
parser.add_argument(
    '--version', action='version', version=__version__,
    help="The pipeline version.")

## arg parser
args = parser.parse_args()
taskIDs = args.taskIDs
run_tag = args.runTag
g_cosmic = args.cosmic_shear
g_columns = args.shear_columns
config_file = args.config
Nmax_proc = args.threads
rng_seed = args.rng_seed
log_level = args.loglevel
running_log = args.sep_running_log

## logging
numeric_level = getattr(logging, log_level.upper(), None)
if not isinstance(numeric_level, int):
    raise ValueError('Invalid log level: %s' % args.loglevel)
logging.basicConfig(format='%(name)s : %(levelname)s - %(message)s', level=numeric_level)
logger = logging.getLogger(__name__)

## host info
try:
    user_name = os.getlogin()
except OSError:
    user_name = 'unknown'
try:
    host_name = os.uname()[1]
except OSError:
    host_name = 'unknown'
date_time = datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')
logger.info(f'~~~~~~~~~~~~ {__version__} started by {user_name} ~~~~~~~~~~~~')
logger.info(f'~~~~~~~~~~~~ {date_time} ~~~~~~~~~~~~')
logger.info(f'~~~~~~~~~~~~ Host: {host_name} ~~~~~~~~~~~~')
logger.info(f'Maximum number of parallel processes: {Nmax_proc}')
logger.info(f'Base seed for the random generator: {rng_seed}')

# ++++++++++++++ config info
## if 0: generate an example configuration file
if '0' in taskIDs:
    example_config_file = 'example.ini'
    RunConfigFile.GenerateExampleConfig(example_config_file, user_name, date_time, __version__)
    logger.info('Pipeline end.')
    sys.exit()
## else: get config info
configs_dict = RunConfigFile.ParseConfig(config_file, taskIDs, run_tag, running_log)

# # ++++++++++++++ Running tasks
start_time0 = time.time()

# load noise info

### survey specified
if configs_dict['imsim']['survey'].lower() == 'kids':

    # diffexpo ?
    multiple_exposures_list = [x.lower()=='diffexpo' for x in configs_dict['imsim']['image_type_list']]
    ## KiDS u only 4 exposures
    N_exposures_list = [4 if x=='u' else 5 for x in configs_dict['imsim']['bands']]

    # varChips ?
    varChips_list = [x.lower()=='varchips' for x in configs_dict['imsim']['image_type_list']]
    N_chips_list = [32] * len(varChips_list)

else:
    multiple_exposures_list = None
    N_exposures_list = None
    varChips_list = None
    N_chips_list = None

### noise and psf info
###### >>> for old versions
try:
    noise_psf_basenames = configs_dict['noise']['noise_psf_basenames']
    label_basename = None 
    noise_basenames = None
    psf_basenames = None
    id_basenames = None

###### >>> new version
except KeyError:
    noise_psf_basenames = None
    label_basename = configs_dict['noise']['label_basename'] 
    noise_basenames = configs_dict['noise']['noise_basenames']
    psf_basenames = configs_dict['noise']['psf_basenames']
    id_basenames = configs_dict['noise']['id_basenames']

### load from file
noise_info = LoadCata.NoiseInfo(configs_dict['noise']['file'], configs_dict['imsim']['bands'], 
                psf_type=configs_dict['noise']['psf_type'],
                noise_psf_basenames=noise_psf_basenames,
                label_basename=label_basename, noise_basenames=noise_basenames, psf_basenames=psf_basenames, id_basenames=id_basenames,
                multiple_exposures_list=multiple_exposures_list, N_exposures_list=N_exposures_list, 
                file4varChips=configs_dict['noise']['file4varChips'], varChips_list=varChips_list, N_chips_list=N_chips_list)    

# varChips dictionary 
varChips_dic = {}
for i_band, band in enumerate(configs_dict['imsim']['bands']):
    varChips = (configs_dict['imsim']['image_type_list'][i_band].lower()=='varchips')
    varChips_dic[band] = varChips

# tile label list
tile_labels = noise_info['label'].to_list()
N_tiles = configs_dict['imsim']['N_tiles']
tile_labels = tile_labels[:N_tiles]

# 1: simulate images
if ('1' in taskIDs) or ('all' in taskIDs):
    logger.info('====== Task 1: simulate images === started ======')
    start_time = time.time()

    ## save all the basic info
    contain_stars = 'True' if configs_dict['star']['file'] else 'False'
    survey = configs_dict['imsim']['survey']
    gal_position_type = configs_dict['gal']['position_type']
    outfile_tmp = os.path.join(configs_dict['work_dirs']['main'], f'basic_info.txt')
    f = open(outfile_tmp, 'w')
    print(f'# Some basic info about the simulated images in this directory\n\
    run_tag            =   {run_tag}\n\
    survey             =   {survey}\n\
    gal_position_type  =   {gal_position_type}\n\
    contain_stars      =   {contain_stars}', file=f)
    if g_columns is None:
        print(f'g_cosmic           =   {g_cosmic[0]} {g_cosmic[1]}', file=f)
    else:
        print(f'g_cosmic           =    variable', file=f)
    if contain_stars=='True':
        star_position_type = configs_dict['star']['position_type']
        print(f'star_position_type =   {star_position_type}', file=f)
    f.close()
    logger.info(f'Setup info saved to {outfile_tmp}')

    ## I/O
    ### for images
    out_dir_tmp = os.path.join(configs_dict['work_dirs']['ima'], 'original')
    if not os.path.exists(out_dir_tmp):
        os.mkdir(out_dir_tmp)
    ### for catalogues
    outcata_dir_tmp = os.path.join(configs_dict['work_dirs']['cata'], 'input')
    if not os.path.exists(outcata_dir_tmp):
        os.mkdir(outcata_dir_tmp)

    ## load galaxy info
    gals_info = LoadCata.GalInfo(configs_dict['gal']['file'], configs_dict['imsim']['bands'],
                        configs_dict['gal']['id_name'], configs_dict['gal']['detection_mag_name'], configs_dict['gal']['mag_name_list'],
                        configs_dict['gal']['RaDec_names'],
                        configs_dict['gal']['shape_names'],
                        configs_dict['gal']['z_name'],
                        rng_seed=rng_seed, mag_cut=configs_dict['gal']['mag_cut'], size_cut=configs_dict['gal']['size_cut'],
                        g_columns=g_columns)

    ## load star info
    if configs_dict['star']['file']:
        stars_info = LoadCata.StarInfo(configs_dict['star']['file'], configs_dict['imsim']['bands'],
                            configs_dict['star']['id_name'], configs_dict['star']['detection_mag_name'], configs_dict['star']['mag_name_list'],
                            RaDec_names=configs_dict['star']['RaDec_names'],
                            mag_cut=configs_dict['star']['mag_cut'])
        star_area = configs_dict['star']['cata_area']
        star_position_type = configs_dict['star']['position_type']
    else:
        stars_info = None
        star_area = None
        star_position_type = None

    ## running
    ImSim.RunParallel_PSFNoisySkyImages(configs_dict['imsim']['survey'], out_dir_tmp, outcata_dir_tmp, rng_seed, configs_dict['imsim']['mag_zero'],
                                            Nmax_proc,
                                            configs_dict['imsim']['N_tiles'], configs_dict['imsim']['bands'], configs_dict['imsim']['pixel_scale_list'], configs_dict['imsim']['image_type_list'],
                                            noise_info,
                                            gals_info, gal_rotation_angles=configs_dict['imsim']['gal_rotation_angles'], g_cosmic=g_cosmic, gal_position_type=[configs_dict['gal']['position_type'], configs_dict['gal']['grid_size']],
                                            stars_area=star_area, stars_info=stars_info, star_position_type=star_position_type,
                                            PSF_map=configs_dict['imsim']['PSF_map'], N_PSF=100, sep_PSF=120,
                                            image_chips=configs_dict['imsim']['image_chips'], image_PSF=[configs_dict['imsim']['image_PSF'], configs_dict['imsim']['image_PSF_size']],
                                            psf_type=configs_dict['noise']['psf_type'])

    logger.info(f'====== Task 1: simulate images === finished in {(time.time()-start_time)/3600.} h ======')

# 2: swarp images
if ('2' in taskIDs) or ('all' in taskIDs):
    logger.info('====== Task 2: swarp images === started ======')
    start_time = time.time()

    ## I/O
    in_dir_tmp = os.path.join(configs_dict['work_dirs']['ima'], 'original')
    if configs_dict['imsim']['PSF_map'][0]:
        in_dir_psf_tmp = os.path.join(in_dir_tmp, 'psf_map')
    if running_log:
        log_dir_tmp = os.path.join(configs_dict['work_dirs']['log'], 'SWarp')
        if not os.path.exists(log_dir_tmp):
            os.mkdir(log_dir_tmp)
    else:
        log_dir_tmp = None
    ### output directory
    for label_tmp in configs_dict['swarp']['image_label_list']:
        out_dir_tmp = os.path.join(configs_dict['work_dirs']['ima'], label_tmp)
        if not os.path.exists(out_dir_tmp):
            os.mkdir(out_dir_tmp)
        if configs_dict['imsim']['PSF_map'][0]:
            out_dir_psf_tmp =  os.path.join(out_dir_tmp, 'psf_map')
            if not os.path.exists(out_dir_psf_tmp):
                os.mkdir(out_dir_psf_tmp)

    ## running
    for i_group, swarp_config in enumerate(configs_dict['swarp']['config_files']):

        swarp_bands = configs_dict['swarp']['bands_group'][i_group]
        swarp_bands = [x.strip() for x in swarp_bands.split(',')]

        only_resample = configs_dict['swarp']['only_resamples'][i_group]

        clean_up_level_tmp = configs_dict['swarp']['clean_up_levels'][i_group]

        out_dir_tmp = os.path.join(configs_dict['work_dirs']['ima'], configs_dict['swarp']['image_label_list'][i_group])
        if configs_dict['imsim']['PSF_map'][0]:
            out_dir_psf_tmp =  os.path.join(out_dir_tmp, 'psf_map')

        for tile_label in tile_labels:

            for band in swarp_bands:

                for gal_rotation_angle in configs_dict['imsim']['gal_rotation_angles']:

                    # original images
                    ## simply resampling
                    image_in = os.path.join(in_dir_tmp, f'tile{tile_label}_band{band}_rot{gal_rotation_angle:.0f}.fits')
                    if not os.path.isfile(image_in):
                        ## exposures
                        image_in = glob.glob(os.path.join(in_dir_tmp, f'tile{tile_label}_band{band}_rot{gal_rotation_angle:.0f}_expo?.fits'))
                        if not image_in:
                            ## chips
                            image_in = os.path.join(in_dir_tmp, f'chips_tile{tile_label}_band{band}_rot{gal_rotation_angle:.0f}', '*.fits')

                    # check weight images
                    if isinstance(image_in, str):
                        if os.path.isfile(image_in.replace('.fits', '.weight.fits')):
                            contain_wei_ima = True
                        else:
                            contain_wei_ima = False
                    else:
                        if os.path.isfile(image_in[0].replace('.fits', '.weight.fits')):
                            contain_wei_ima = True
                        else:
                            contain_wei_ima = False

                    ### run
                    image_out = os.path.join(out_dir_tmp, f'tile{tile_label}_band{band}_rot{gal_rotation_angle:.0f}.fits')
                    Astromatic.SwarpImage(image_in, swarp_config,
                                        image_out,
                                        only_resample=only_resample, contain_wei_ima=contain_wei_ima,
                                        running_log=running_log, log_dir=log_dir_tmp,
                                        swarp_path=configs_dict['swarp']['cmd'], NTHREADS=Nmax_proc,
                                        clean_up_level=clean_up_level_tmp)
                    ### psf map
                    if configs_dict['imsim']['PSF_map'][0]:
                        try:
                            image_in_psf = os.path.join(in_dir_psf_tmp, os.path.basename(image_in))
                            psf_map_existence = os.path.isfile(image_in_psf)
                        except TypeError:
                            image_in_psf = [os.path.join(in_dir_psf_tmp, os.path.basename(image_in_tmp)) for image_in_tmp in image_in]
                            psf_map_existence = os.path.isfile(image_in_psf[0])

                        #### run
                        if psf_map_existence:
                            image_out = os.path.join(out_dir_psf_tmp, f'tile{tile_label}_band{band}_rot{gal_rotation_angle:.0f}.fits')
                            Astromatic.SwarpImage(image_in_psf, swarp_config,
                                                image_out,
                                                only_resample=only_resample, contain_wei_ima=contain_wei_ima,
                                                running_log=running_log, log_dir=log_dir_tmp,
                                                swarp_path=configs_dict['swarp']['cmd'], NTHREADS=Nmax_proc,
                                                clean_up_level=clean_up_level_tmp)

    logger.info(f'====== Task 2: swarp images === finished in {(time.time()-start_time)/3600.} h ======')

# 3: detect objects
if ('3' in taskIDs) or ('all' in taskIDs):
    logger.info('====== Task 3: detect objects === started ======')
    start_time = time.time()

    ## basic info
    detection_band = configs_dict['sex']['detection_band']
    pixel_scale = configs_dict['sex']['pixel_scale']
    image_label = configs_dict['sex']['image_label']
    logger.info(f'Band for detection: {detection_band}')
    logger.info(f'pixel scale: {pixel_scale}')
    logger.info(f'image type: {image_label}')
    ### seeing info
    ###### only matters for CLASS_STAR classifier
    try:
        SeeingFWHM_list = noise_info[f'seeing_{detection_band}'].to_list()
    except KeyError:
        try:
            SeeingFWHM_list = noise_info[f'seeing_{detection_band}_expo0'].to_list()
        except KeyError:
            try:
                SeeingFWHM_list = noise_info[f'seeing_{detection_band}_expo0_chip0'].to_list()
            except KeyError:
                SeeingFWHM_list = [1.0] * N_tiles
                logger.warning('seeing not found in noise info, do NOT trust CLASS_STAR!')

    SeeingFWHM_list = SeeingFWHM_list[:N_tiles]

    ## I/O
    ori_cata_dir_tmp = os.path.join(configs_dict['work_dirs']['cata'], 'input')
    in_dir_tmp = os.path.join(configs_dict['work_dirs']['ima'], image_label)
    out_dir_tmp = os.path.join(configs_dict['work_dirs']['cata'], 'SExtractor')
    if not os.path.exists(out_dir_tmp):
        os.mkdir(out_dir_tmp)
    if running_log:
        log_dir_tmp = os.path.join(configs_dict['work_dirs']['log'], 'SExtractor')
        if not os.path.exists(log_dir_tmp):
            os.mkdir(log_dir_tmp)
    else:
        log_dir_tmp = None
    if configs_dict['sex']['cross_match']:
        out_dir_cross = os.path.join(configs_dict['work_dirs']['cata'], 'CrossMatch')
        if not os.path.exists(out_dir_cross):
            os.mkdir(out_dir_cross)
    # CHECKIMAGE
    if (configs_dict['sex']['checkimage_type'] is not None) and (configs_dict['sex']['checkimage_type'].upper() != 'NONE'):
        CHECKIMAGE_dir = os.path.join(in_dir_tmp, configs_dict['sex']['checkimage_type'])
        if not os.path.exists(CHECKIMAGE_dir):
            os.mkdir(CHECKIMAGE_dir)

    ## work pool
    N_sex = int(Nmax_proc/4.)
    if N_sex < 1:
        N_sex = 1
    logger.info(f'Number of processes for SExtractor: {N_sex}')
    work_pool = mp.Pool(processes=N_sex)
    proc_list = []
    for i_tile, tile_label in enumerate(tile_labels):

        SeeingFWHM = SeeingFWHM_list[i_tile]

        for gal_rotation_angle in configs_dict['imsim']['gal_rotation_angles']:

            ### input
            ImageFile1 = os.path.join(in_dir_tmp, f'tile{tile_label}_band{detection_band}_rot{gal_rotation_angle:.0f}.fits')
            WeightFile1 = ImageFile1.replace('.fits', '.weight.fits')
            if not os.path.isfile(WeightFile1):
                WeightFile1 = None

            ### output
            CatalogueFile = os.path.join(out_dir_tmp, os.path.basename(ImageFile1).replace('.fits', '.sex'))

            ### running
            proc = work_pool.apply_async(func=Astromatic.SExtractorCatalogue,
                                    args=(CatalogueFile, pixel_scale, SeeingFWHM,
                                            ImageFile1, WeightFile1,
                                            None, None,
                                            running_log, log_dir_tmp,
                                            configs_dict['sex']['cmd'], configs_dict['sex']['config_file'], configs_dict['sex']['param_file'],
                                            configs_dict['sex']['filter_file'], configs_dict['sex']['starNNW_file'],
                                            configs_dict['sex']['checkimage_type'], True,
                                            configs_dict['imsim']['mag_zero'],
                                            configs_dict['sex']['clean_up_level']))
            proc_list.append(proc)

    work_pool.close()
    work_pool.join()
    ### check for any errors during run
    for proc in proc_list:
        proc.get()

    ## Mark stars
    tmp = glob.glob(os.path.join(ori_cata_dir_tmp, f'stars_info_*.feather'))
    if tmp:
        logger.info('Mask Stars based on input info.')

        for tile_label in tile_labels:

            input_cata = pd.read_feather(os.path.join(ori_cata_dir_tmp, f'stars_info_tile{tile_label}.feather'))

            for gal_rotation_angle in configs_dict['imsim']['gal_rotation_angles']:

                detec_file = os.path.join(out_dir_tmp, f'tile{tile_label}_band{detection_band}_rot{gal_rotation_angle:.0f}.feather')
                detec_cata = pd.read_feather(detec_file)

                id_list = ['index_input', 'NUMBER']
                position_list = [['RA_input', 'DEC_input'], ['X_WORLD', 'Y_WORLD']]
                mag_list = [f'{detection_band}_input', 'MAG_AUTO']

                matched_cata, _, _ = CrossMatch.run_position2id(input_cata, detec_cata, id_list, position_list, mag_list,
                                    outDir=None, basename=None, save_matched=False, save_false=False, save_missed=False,
                                    r_max=0.6/3600., k=4, mag_closest=True, running_info=False)

                mask_stars = detec_cata['NUMBER'].isin(matched_cata['id_detec'])
                detec_cata.loc[mask_stars, 'perfect_flag_star'] = 1
                detec_cata.loc[~mask_stars, 'perfect_flag_star'] = 0
                detec_cata = detec_cata.astype({'perfect_flag_star': int})
                detec_cata.to_feather(detec_file)

    ## cross-match with the input catalogue
    if configs_dict['sex']['cross_match']:

        logger.info('cross-match with the input catalogue.')

        for tile_label in tile_labels:

            input_cata = pd.read_feather(os.path.join(ori_cata_dir_tmp, f'gals_info_tile{tile_label}.feather'))

            ## magnitude pre-selection
            input_cata = input_cata[input_cata[f'{detection_band}_input']<=configs_dict['sex']['mag_faint_cut']]
            input_cata.reset_index(drop=True, inplace=True)

            for gal_rotation_angle in configs_dict['imsim']['gal_rotation_angles']:

                detec_file = os.path.join(out_dir_tmp, f'tile{tile_label}_band{detection_band}_rot{gal_rotation_angle:.0f}.feather')
                detec_cata = pd.read_feather(detec_file)

                ## select galaxies
                try:
                    mask_gals = (detec_cata['perfect_flag_star'] == 0)
                    detec_cata = detec_cata[mask_gals]
                    detec_cata.reset_index(drop=True, inplace=True)
                except KeyError:
                    pass

                # output info
                basename_cross =  f'tile{tile_label}_rot{gal_rotation_angle:.0f}'

                # cross-match
                id_list = ['index_input', 'NUMBER']
                position_list = [['RA_input', 'DEC_input'], ['X_WORLD', 'Y_WORLD']]
                mag_list = [f'{detection_band}_input', 'MAG_AUTO']
                _, _, _ = CrossMatch.run_position2id(input_cata, detec_cata, id_list, position_list, mag_list,
                                    outDir=out_dir_cross, basename=basename_cross, save_matched=configs_dict['sex']['save_matched'], save_false=configs_dict['sex']['save_false'], save_missed=configs_dict['sex']['save_missed'],
                                    r_max=configs_dict['sex']['r_max']/3600., k=4, mag_closest=configs_dict['sex']['mag_closest'], running_info=True)

    logger.info(f'====== Task 3: detect objects === finished in {(time.time()-start_time)/3600.} h ======')

# 4: measure photometry
if ('4' in taskIDs) or ('all' in taskIDs):
    logger.info('====== Task 4: measure photometry === started ======')
    start_time = time.time()

    if configs_dict['MP']['method'].lower() == 'gaap':
        logger.info('Use GAaP for photometry measurement.')

        ## I/O
        ori_cata_dir_tmp = os.path.join(configs_dict['work_dirs']['cata'], 'input')
        in_cata_dir_tmp = os.path.join(configs_dict['work_dirs']['cata'], 'SExtractor')
        out_dir_tmp = os.path.join(configs_dict['work_dirs']['cata'], 'photometry')
        if not os.path.exists(out_dir_tmp):
            os.mkdir(out_dir_tmp)
        if running_log:
            log_dir_tmp = os.path.join(configs_dict['work_dirs']['log'], 'GAaP')
            if not os.path.exists(log_dir_tmp):
                os.mkdir(log_dir_tmp)
        else:
            log_dir_tmp = None

        ## Initialise the GAaP wrapper
        gaap = GAaP.GAaPwrapper(configs_dict['MP']['gaap_dir'], out_dir_tmp,
                star_mag_cut=configs_dict['MP']['star_mag_cut'],
                mag_zero=configs_dict['imsim']['mag_zero'],
                min_aper=configs_dict['MP']['min_aper'], max_aper=configs_dict['MP']['max_aper'],
                mag_1sigma_limits=configs_dict['MP']['band_1sigma_limits'],
                running_log=running_log, log_dir=log_dir_tmp,
                clean_up_level=configs_dict['MP']['clean_up_level'])

        ## work pool
        rots_tasks = configs_dict['imsim']['gal_rotation_angles'] * len(tile_labels)
        tile_labels_tasks = np.repeat(tile_labels, len(configs_dict['imsim']['gal_rotation_angles']))
        N_tasks = len(rots_tasks)
        logger.info(f'Total number of tasks (tiles*rot): {N_tasks}')
        N_gaap = int(Nmax_proc/N_tasks) # for each subprocess
        if N_gaap < 1:
            N_gaap = 1
        logger.info(f'number of processes in each GAaP run: {N_gaap}')
        proc_list = []
        i_worker = 0
        i_tile = 0
        while True:
            N_running = len(mp.active_children())
            logger.debug(f'Number of running {N_running}')
            if i_worker == N_tasks:
                break
            elif N_running >= Nmax_proc:
                time.sleep(5.)
            else:
                tile_label = tile_labels_tasks[i_worker]
                gal_rotation_angle = rots_tasks[i_worker]

                ### star info for psf estimation
                if not configs_dict['MP']['use_PSF_map']:
                    star_info_file = os.path.join(ori_cata_dir_tmp, f'stars_info_tile{tile_label}.feather')
                    star_info = pd.read_feather(star_info_file)
                else:
                    star_info = None

                ### detection file
                MP_detec_band = configs_dict['MP']['detection_band']
                SKYcataFile = os.path.join(in_cata_dir_tmp, f'tile{tile_label}_band{MP_detec_band}_rot{gal_rotation_angle:.0f}.feather')

                ### image files
                SKYimaFile_list = []
                SKYweiFile_list = []
                if configs_dict['MP']['use_PSF_map']:
                    PSFimaFile_list = []
                else:
                    PSFimaFile_list = None
                for i_band, band in enumerate(configs_dict['MP']['bands']):

                    in_ima_dir_tmp = os.path.join(configs_dict['work_dirs']['ima'], configs_dict['MP']['image_label_list'][i_band])
                    SKYimaFile_tmp = os.path.join(in_ima_dir_tmp, f'tile{tile_label}_band{band}_rot{gal_rotation_angle:.0f}.fits')
                    SKYimaFile_list.append(SKYimaFile_tmp)

                    #### weight image
                    SKYweiFile_tmp = SKYimaFile_tmp.replace('.fits', '.weight.fits')
                    if not os.path.isfile(SKYweiFile_tmp):
                        SKYweiFile_tmp = None
                    SKYweiFile_list.append(SKYweiFile_tmp)

                    #### psf map
                    if configs_dict['MP']['use_PSF_map']:
                        in_psf_dir_tmp = os.path.join(in_ima_dir_tmp, 'psf_map')
                        PSFimaFile_tmp = os.path.join(in_psf_dir_tmp, os.path.basename(SKYimaFile_tmp))
                        PSFimaFile_list.append(PSFimaFile_tmp)

                # output
                FinalFile = os.path.join(out_dir_tmp, f'tile{tile_label}_rot{gal_rotation_angle:.0f}.feather')

                proc = mp.Process(target=gaap.RunSingleTile, 
                    args=(N_gaap, FinalFile, SKYcataFile, configs_dict['MP']['bands'], SKYimaFile_list, SKYweiFile_list, PSFimaFile_list, star_info))
                i_worker += 1
                proc.start()
                proc_list.append(proc)
                time.sleep(1.)

        ### check for any errors during run
        for proc in proc_list:
            proc.join()

    logger.info(f'====== Task 4: measure photometry === finished in {(time.time()-start_time)/3600.} h ======')

# 5: measure photo-z
if ('5' in taskIDs) or ('all' in taskIDs):
    logger.info('====== Task 5: measure photo-z === started ======')
    start_time = time.time()

    if configs_dict['MZ']['method'].lower() == 'bpz':
        logger.info('Use BPZ for photo-z measurement.')

        ## I/O
        in_cata_dir_tmp = os.path.join(configs_dict['work_dirs']['cata'], 'photometry')
        out_dir_tmp = os.path.join(configs_dict['work_dirs']['cata'], 'photo_z')
        if not os.path.exists(out_dir_tmp):
            os.mkdir(out_dir_tmp)
        tmp_dir_tmp = os.path.join(out_dir_tmp, 'tmp_bpz')
        if not os.path.exists(tmp_dir_tmp):
            os.mkdir(tmp_dir_tmp)
        if running_log:
            log_dir_tmp = os.path.join(configs_dict['work_dirs']['log'], 'BPZ')
            if not os.path.exists(log_dir_tmp):
                os.mkdir(log_dir_tmp)
        else:
            log_dir_tmp = None

        ## Initialise the BPZ wrapper
        bpz = BPZ.BPZwrapper(configs_dict['MZ']['python2_cmd'], configs_dict['MZ']['BPZ_dir'],
                out_dir_tmp, tmp_dir_tmp,
                configs_dict['MZ']['bands'],
                configs_dict['MZ']['bands_CataName'], configs_dict['MZ']['banderrs_CataName'], configs_dict['MZ']['bands_FilterName'],
                photo_sys=configs_dict['MZ']['photo_sys'], prior_band=configs_dict['MZ']['prior_band'], prior_name=configs_dict['MZ']['prior_name'],
                templates_name=configs_dict['MZ']['templates_name'], interpolation=configs_dict['MZ']['interpolation'],
                lkl_zmin=configs_dict['MZ']['lkl_zmin'], lkl_zmax=configs_dict['MZ']['lkl_zmax'], lkl_dz=configs_dict['MZ']['lkl_dz'], lkl_odds=configs_dict['MZ']['lkl_odds'], lkl_min_rms=configs_dict['MZ']['lkl_min_rms'],
                running_log=running_log, log_dir=log_dir_tmp)

        ## work pool
        N_BPZ = int(Nmax_proc/20)
        if N_BPZ < 1:
            N_BPZ = 1
        logger.info(f'Number of processes for BPZ: {N_BPZ}')
        logger.info(f'  NOTE: each processes of BPZ takes multiple cores')
        work_pool = mp.Pool(processes=N_BPZ)
        proc_list = []
        for tile_label in tile_labels:

            for gal_rotation_angle in configs_dict['imsim']['gal_rotation_angles']:

                # input
                in_cata_file = os.path.join(in_cata_dir_tmp, f'tile{tile_label}_rot{gal_rotation_angle:.0f}.feather')

                proc = work_pool.apply_async(func=bpz.RunSingleTile,
                                args=(in_cata_file,))
                proc_list.append(proc)

        work_pool.close()
        work_pool.join()
        ## check for any errors during run
        for proc in proc_list:
            proc.get()

    logger.info(f'====== Task 5: measure photo-z === finished in {(time.time()-start_time)/3600.} h ======')

# 6: measure galaxy shapes
if ('6' in taskIDs) or ('all' in taskIDs):
    logger.info('====== Task 6: measure galaxy shapes === started ======')
    start_time = time.time()

    # basic info
    detection_band = configs_dict['MS']['detection_band']

    if configs_dict['MS']['method'].lower() == 'lensfit':
        logger.info('Use lensfit for shape measurement.')
        logger.info('   NOTE: the original lensfit weights need to be globally recalibrated.')

        ## I/O
        in_cata_dir_tmp = os.path.join(configs_dict['work_dirs']['cata'], 'SExtractor')
        out_dir_tmp = os.path.join(configs_dict['work_dirs']['cata'], 'shapes')
        if not os.path.exists(out_dir_tmp):
            os.mkdir(out_dir_tmp)
        tmp_dir_tmp = os.path.join(out_dir_tmp, 'tmp_lensfit')
        if not os.path.exists(tmp_dir_tmp):
            os.mkdir(tmp_dir_tmp)
        if running_log:
            log_dir_tmp = os.path.join(configs_dict['work_dirs']['log'], 'lensfit')
            if not os.path.exists(log_dir_tmp):
                os.mkdir(log_dir_tmp)
        else:
            log_dir_tmp = None

        ## Initialise the lensfit wrapper
        lensfit_cores = configs_dict['MS']['lensfit_cores']
        if lensfit_cores > Nmax_proc:
            logger.warning(f'required lensfit_cores {lensfit_cores} > maximum cores allowed {Nmax_proc}!')
            lensfit_cores = Nmax_proc
            logger.warning(f'Make lensfit_cores = Nmax_proc = {lensfit_cores}')
        lensfit = LensFit.LensFITwrapper(configs_dict['MS']['lensfit_dir'], out_dir_tmp, tmp_dir_tmp,
                    PSF_OVERSAMPLING=configs_dict['MS']['PSF_OVERSAMPLING'], PECUT=configs_dict['MS']['PECUT'], PRCUT=configs_dict['MS']['PRCUT'], LCUT=configs_dict['MS']['LCUT'], CAMERA=configs_dict['MS']['CAMERA'],
                    postage_size=configs_dict['MS']['postage_size'], start_exposure=configs_dict['MS']['start_exposure'], end_exposure=configs_dict['MS']['end_exposure'], start_mag=configs_dict['MS']['start_mag'], end_mag=configs_dict['MS']['end_mag'],
                    lensfit_cores=lensfit_cores,
                    running_log=running_log, log_dir=log_dir_tmp)

        ## work pool
        N_lensfit = int(Nmax_proc/lensfit_cores)
        if N_lensfit < 1:
            N_lensfit = 1
        logger.info(f'Number of processes for lensfit: {N_lensfit}')
        logger.info(f'  NOTE: each processes of lensfit takes {lensfit_cores} cores')
        work_pool = mp.Pool(processes=N_lensfit)
        proc_list = []
        for i_band, band in enumerate(configs_dict['MS']['bands']):
            logger.info(f'Measure shapes from band {band}...')

            WAVEBAND = band.upper()
            MS_label = configs_dict['MS']['image_label_list'][i_band]
            in_ima_dir_tmp = os.path.join(configs_dict['work_dirs']['ima'], MS_label)
            varChips = varChips_dic[band]

            for tile_label in tile_labels:

                ### data info
                weight_dir = None
                psf_dir = os.path.join(in_ima_dir_tmp, f'psf_tile{tile_label}_band{band}')

                ### psf coefficients for lensfit
                psf_coeff_dir = lensfit.LensfitShape_psf(psf_dir, varChips=varChips)

                for gal_rotation_angle in configs_dict['imsim']['gal_rotation_angles']:

                    # output
                    output_file = f'tile{tile_label}_band{detection_band}_rot{gal_rotation_angle:.0f}.feather'

                    # detection catalogue
                    path_input_cata = os.path.join(in_cata_dir_tmp, f'tile{tile_label}_band{detection_band}_rot{gal_rotation_angle:.0f}.feather')

                    # images
                    chip_dir = os.path.join(in_ima_dir_tmp, f'chips_tile{tile_label}_band{detection_band}_rot{gal_rotation_angle:.0f}')

                    # run
                    proc = work_pool.apply_async(func=lensfit.LensfitShape,
                                args=(output_file,
                                        path_input_cata, chip_dir, psf_coeff_dir, weight_dir,
                                        WAVEBAND))
                    proc_list.append(proc)

        work_pool.close()
        work_pool.join()
        ### check for any errors during run
        for proc in proc_list:
            proc.get()

    logger.info(f'====== Task 6: measure galaxy shapes === finished in {(time.time()-start_time)/3600.} h ======')

# 7: create a combined catalogue
if ('7' in taskIDs) or ('all' in taskIDs):

    logger.info('====== Task 7: create a combined catalogue === started ======')
    start_time = time.time()

    # detection info
    out_dir_detec = os.path.join(configs_dict['work_dirs']['cata'], 'SExtractor')
    if not os.path.exists(out_dir_detec):
        raise Exception('Detection files are not generated!\n\
        Task 3 is required for create a combined catalogue.')

    # input info
    out_dir_input = os.path.join(configs_dict['work_dirs']['cata'], 'input')
    if not os.path.exists(out_dir_input):
        raise Exception('There is no input info, something very bad happened!')

    # directories for other info
    ## CrossMatch info
    out_dir_cross = os.path.join(configs_dict['work_dirs']['cata'], 'CrossMatch')
    ## photometry info
    out_dir_photometry = os.path.join(configs_dict['work_dirs']['cata'], 'photometry')
    ## photo-z info
    out_dir_photoz = os.path.join(configs_dict['work_dirs']['cata'], 'photo_z')
    ## shape info
    out_dir_shape = os.path.join(configs_dict['work_dirs']['cata'], 'shapes')

    # combine catalogues for each run
    for tile_label in tile_labels:
        for gal_rotation_angle in configs_dict['imsim']['gal_rotation_angles']:

            logger.info(f'Combining outputs for tile {tile_label}, rot {gal_rotation_angle}...')

            # detection catalogue as the base
            infile_tmp = glob.glob(os.path.join(out_dir_detec, f'tile{tile_label}_band*_rot{gal_rotation_angle:.0f}.feather'))[0]
            data_final = pd.read_feather(infile_tmp)

            # CrossMatch
            infile_tmp = os.path.join(out_dir_cross, f'tile{tile_label}_rot{gal_rotation_angle:.0f}_matched.feather')
            if os.path.isfile(infile_tmp):
                tmp_info = pd.read_feather(infile_tmp)
                data_final = data_final.merge(tmp_info, left_on='NUMBER', right_on='id_detec', how='left')
                data_final.drop(columns=['id_detec'], inplace=True)

                ## input info
                infile_tmp = os.path.join(out_dir_input, f'gals_info_tile{tile_label}.feather')
                tmp_info = pd.read_feather(infile_tmp)
                ### rename
                tmp_info.rename(columns={f'e1_input_rot{int(gal_rotation_angle)}': 'e1_input', f'e2_input_rot{int(gal_rotation_angle)}': 'e2_input'}, inplace=True)
                tmp_info.drop(columns=[s for s in tmp_info.columns if ("e1_input_" in s) or ("e2_input_" in s)], inplace=True)
                ### merge
                data_final = data_final.merge(tmp_info, left_on='id_input', right_on='index_input', how='left')
                data_final.drop(columns=['index_input'], inplace=True)
            else:
                logger.warning('CrossMatch is not performed, the final catalogue will not contain input info.')

            # photometry
            infile_tmp = os.path.join(out_dir_photometry, f'tile{tile_label}_rot{gal_rotation_angle:.0f}.feather')
            if os.path.isfile(infile_tmp):
                tmp_info = pd.read_feather(infile_tmp)
                data_final = data_final.merge(tmp_info, left_on='NUMBER', right_on='id_detec', how='left')
                data_final.drop(columns=['id_detec'], inplace=True)
            else:
                logger.warning('MeasurePhotometry is not performed, the final catalogue will not contain photometry info.')

            # photo-z
            infile_tmp = os.path.join(out_dir_photoz, f'tile{tile_label}_rot{gal_rotation_angle:.0f}.feather')
            if os.path.isfile(infile_tmp):
                tmp_info = pd.read_feather(infile_tmp)
                data_final = data_final.merge(tmp_info, left_on='NUMBER', right_on='id_detec', how='left')
                data_final.drop(columns=['id_detec'], inplace=True)
            else:
                logger.warning('MeasurePhotoz is not performed, the final catalogue will not contain photo-z info.')

            # shape
            file_list = glob.glob(os.path.join(out_dir_shape, f'tile{tile_label}_band*_rot{gal_rotation_angle:.0f}.feather'))
            if file_list:
                for infile_tmp in file_list:
                    band = re.search(r'_band(.*)_rot', infile_tmp).group(1)
                    tmp_info = pd.read_feather(infile_tmp)
                    tmp_info = tmp_info.add_suffix(f'_{band}')
                    data_final = data_final.merge(tmp_info, left_on='NUMBER', right_on=f'id_detec_{band}', how='left')
                    data_final.drop(columns=[f'id_detec_{band}'], inplace=True)
            else:
                logger.warning('MeasureShape is not performed, the final catalogue will not contain shape info.')

            # dummy values for nan
            data_final.fillna(-999, inplace=True)
            data_final = data_final.astype({'id_input': int})

            # save
            if configs_dict['CC']['format'] == 'feather':
                outfile = os.path.join(configs_dict['work_dirs']['cata'], f'tile{tile_label}_rot{gal_rotation_angle:.0f}_combined.feather')
                data_final.to_feather(outfile)
            elif configs_dict['CC']['format'] == 'csv':
                outfile = os.path.join(configs_dict['work_dirs']['cata'], f'tile{tile_label}_rot{gal_rotation_angle:.0f}_combined.csv')
                data_final.to_csv(outfile, index=False)
            elif configs_dict['CC']['format'] == 'fits':
                outfile = os.path.join(configs_dict['work_dirs']['cata'], f'tile{tile_label}_rot{gal_rotation_angle:.0f}_combined.fits')
                if os.path.isfile(outfile):
                    os.remove(outfile)
                Table.from_pandas(data_final).write(outfile, format='fits')
            logger.info(f'Combined catalogue saved as {outfile}')

    # clean up
    if configs_dict['CC']['clean_up_level']:
        logger.info('All intermediate catalogues will be removed!')
        for dir_tmp in [out_dir_detec, out_dir_cross, out_dir_photometry, out_dir_photoz, out_dir_shape]:
            if os.path.exists(dir_tmp):
                shutil.rmtree(dir_tmp)

    logger.info(f'====== Task 7: create a combined catalogue === finished in {(time.time()-start_time)/3600.} h ======')

logger.info(f'~~~~~~~~~~~~ {__version__} finished ~~~~~~~~~~~~')
logger.info(f'~~~~~~~~~~~~ total running time {(time.time()-start_time0)/3600.} h ~~~~~~~~~~~~')
tmp = configs_dict['work_dirs']['main']
logger.info(f'~~~~~~~~~~~~ All outputs saved in {tmp}')
logger.info(f'~~~~~~~~~~~~ Enjoy the science ~~~~~~~~~~~~')
