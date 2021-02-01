# -*- coding: utf-8 -*-
# @Author: lshuns
# @Date:   2020-12-21 11:44:14
# @Last Modified by:   lshuns
# @Last Modified time: 2021-02-01 11:22:36

### main module to run the whole pipeline

import re
import os
import sys
import time
import logging
import datetime
import argparse
import configparser
import distutils.util

import numpy as np
import pandas as pd
import multiprocessing as mp

import GAaP
import ImSim
import LoadCata
import LensFit
import Astromatic
import CrossMatch

__version__ = "MultiBand_ImSim v0.2"

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
    "-c", "--config", type=str, metavar='config_file',
    help="A configuration file including all necessary setups.\n"
            "   To generate an example file, use `python Run.py 0`.")
parser.add_argument(
    "--threads", type=int, default=int(mp.cpu_count()*0.5), metavar='number_threads',
    help="The maximum number of threads to use. \n"
        "   (default: %(default)s [half of the total CPU count])")
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
run_tag = args.runTag
Nmax_proc = args.threads
log_level = args.loglevel
running_log = args.sep_running_log

## logging
numeric_level = getattr(logging, log_level.upper(), None)
if not isinstance(numeric_level, int):
    raise ValueError('Invalid log level: %s' % args.loglevel)
logging.basicConfig(format='%(name)s : %(levelname)s - %(message)s', level=numeric_level)
logger = logging.getLogger(__name__)
logger.setLevel(numeric_level)

## basic info
user_name = os.getlogin()
host_name = os.uname()[1]
date_time = datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')

# ++++++++++++++ parser for configuration file
config = configparser.ConfigParser(allow_no_value=True, 
                                inline_comment_prefixes='#', 
                                empty_lines_in_values=False, 
                                interpolation=configparser.ExtendedInterpolation())

# 0: generate an example configuration file
if '0' in args.taskIDs:

    config = f"# Example configuration file \n\
#   for {__version__} \n\
# Created by {user_name} ({date_time})\n\
# Note: If you do not use certain sections, just leave them whatever they are.\n\
#         They will not be loaded.\n\
\n\n\
################################## Paths ################################################\n\
[Paths]\n\n\
config_dir =  .../MultiBand_ImSim/config/      # directory to all the configuration files\n\
out_dir =                                      # main directory for all the outputs\n\
\n\n\
################################## GalInfo ##############################################\n\
[GalInfo]\n\n\
cata_file =                                    # input galaxy mock catalogue\n\
                                               # supported file types: feather, csv, fits\n\
position_type =         true                   # position to be used\n\
                                               #    true (use positions from the input catalogue)\n\
                                               #    grid (put in a grid)\n\
mag_min_cut =           16                     # brightest galaxies to be simulated\n\
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
mag_min_cut =           10                     # brightest stars to be simulated\n\
# column names to the desired info\n\
id_name =               index                  # unique star id\n\
detection_mag_name =    r                      # correspond to the `detection_band` in [ImSim] \n\
mag_name_list =         u, g, r, i, Z, Y, J, H, Ks \n\
                                               # correspond to the the `band_list` in [ImSim]\n\
RaDec_names =           ra, dec \n             # not required, if stars are randomly placed\n\
\n\n\
################################## NoiseInfo ##############################################\n\
[NoiseInfo]\n\n\
cata_file =                                    # input noise background & psf catalogue\n\
                                               # supported file types: feather, csv, fits\n\
                                               # NOTE: tiles are orderly selected\n\
noise_psf_basenames =   none, none, none, none, none, none\n\
                                               # base names for noise and psf info\n\
                                               # order: \n\
                                               #    label, rms, seeing, MoffatBeta, psf_e1, psf_e2\n\
                                               # the real column name is associated with band labels as `rms_r` etc\n\
                                               # not all required, for those missed, simply feed none\n\
\n\n\
################################## ImSim ###################################################\n\
[ImSim]\n\n\
survey =                KiDS_sameExpo          # survey being simulated\n\
                                               # current supported surveys:\n\
                                               #    simple_Ndeg: N can be any int corresponding to the tile sky area\n\
                                               #    KiDS_sameExpo: KiDS-like images using same psf and noise for all exposures\n\
N_tiles =               1                      # number of tiles to be simulated\n\
                                               # make sure the NoiseInfo cata covers more than this requirement\n\
                                               # GalInfo can cover less than this requirement,\n\
                                               #    in which case repeating patterns will be produced\n\
                                               # NOTE: the total output tiles = N_tiles * N_rotations (specified below)\n\
gal_rotation_angles =   0                      # degrees (put more values separated with ',' if needed)\n\
PSF_map =               False                  # output the corresponding PSF map or not\n\
                                               # can be used by GAaP, but not mandatory if stars are simulated\n\
rng_seed =              940120                 # base seed for the random number generator\n\
mag_zero =              30                     # simulated magnitude zero point\n\
g_cosmic =              0, 0                   # cosmic shear values\n\
band_list =             u, g, r, i, Z, Y, J, H, Ks\n\
                                               # bands being simulated\n\
pixel_scale_list =      0.214, 0.214, 0.214, 0.214, 0.34, 0.34, 0.34, 0.34, 0.34\n\
                                               # pixel scale for each band image\n\
image_chips =           False, False, True, False, False, False, False, False, False\n\
                                               # save individual chips or not\n\
                                               # required by lensfit\n\
image_PSF =             False, False, True, False, False, False, False, False, False\n\
                                               # save individual psf\n\
                                               # required by lensfit\n\
\n\n\
################################## SWarp ###################################################\n\
[SWarp]\n\n\
# for coadding or resampling\n\
swarp_path =            swarp                  # path to the SWarp code\n\
swarp_configs =         coadd_theli.swarp, coadd_aw.swarp\n\
                                               # SWarp configuration files\n\
                                               # more than one files are supported\n\
                                               #    in which case, more than one treatments are applied\n\
swarp_bands_group =     [r], [u, g, r, i]      # bands to be swarped\n\
                                               # NOTE: the group corresponding to the same config should be surrounded by `[]`\n\
swarp_labels =          THELI, AW              # name to label the swaped results, one to each group\n\
only_resamples =        False, False           # set it True if only resampling but not coadding\n\
clean_up_levels =       0, 0                   # clean up level\n\
                                               #    0: none\n\
                                               #    1: original images\n\
                                               # NOTE: careful about cleaning before other swarp applied\n\
\n\n\
################################## SExtractor #################################################\n\
[SExtractor]\n\n\
# for detection\n\
sex_path =              sex                    # path to the SExtractor code\n\
detection_band =        r                      # band for detection\n\
image_label =           THELI                  # label for the image type, can be either:\n\
                                               #    original (for original simulated images)\n\
                                               #    any label specified in `swarp_labels`\n\
mask_stars =            True                   # mask stars using input info\n\
cross_match =           True                   # cross-match with the input catalogue\n\
                                               #    in which case, catalogues with match info will be saved\n\
                                               #    see next section for configuration\n\
sex_config =            kids_sims.sex          # SExtractor configuration file\n\
sex_param =             sex_image.param        # SExtractor parameter file\n\
sex_filter =            default.conv           # SExtractor filter file\n\
sex_starNNW =           default.nnw            # SExtractor Neural-Network_Weight table file\n\
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
dmag_max =              0.5                    # allowed maximum magnitude difference\n\
r_max =                 0.5                    # (arcsec) allowed maximum separation\n\
\n\n\
################################## MeasurePhotometry ########################################################\n\
[MeasurePhotometry]\n\n\
method =                GAaP                   # method for photometry measurement\n\
                                               # supported method:\n\
                                               #    GAaP\n\
detection_band =        r                      # band with detection catalogue\n\
band_list =             u, g, r, i, Z, Y, J, H, Ks\n\
                                               # bands being measured\n\
image_label_list =      AW, AW, AW, AW, original, original, original, original, original\n\
                                               # a list of labels for the image types, can be either:\n\
                                               #    original (for original simulated images)\n\
                                               #    any label specified in `swarp_labels`\n\
\n\n\
[GAaP]\n\n\
gaap_dir =                                     # directory containing GAaP bins and libs\n\
min_aper =              0.7                    # minimum aperture size\n\
max_aper =              2.0                    # maximum aperture size\n\
use_psf_map =           False                  # use separate psf map for psf estimation\n\
                                               # only required if stars are not simulated\n\
star_mag_cut =          16, 20                 # magnitude range for stars used for PSF estimation\n\
                                               # not used if PSF map is provided\n\
clean_up_level =        0                      # clean up level\n\
                                               #    0: none\n\
                                               #    1: tmp directory\n\
                                               #    2: and *.gaap files\n\
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
python2_env =           python2                # proper direction to python2\n\
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
cata_list =                                    # list of catalogues to be combined\n\
                                               # all legitimate items:\n\
                                               # (joining desired items with `,`, leave it empty if no combination is required)\n\
                                               #    input, detection, photometry, photo-z, shape\n\
"
    # write out the example config file
    with open('example.ini', 'w') as configfile:
        configfile.write(config)
        logger.info('An example configuration file `example.ini` is generated in the current directory.')
        logger.info('Pipeline end.')
    sys.exit()

### read the config file
elif args.config != None:
    config.read(args.config)
else:
    raise Exception("No configuration file provided! \n" 
                    "------> To generate an example file, use `python Run.py 0`.")

# # ++++++++++++++ Running tasks
logger.info(f'~~~~~~~~~~~~ {__version__} started by {user_name} ~~~~~~~~~~~~')
logger.info(f'~~~~~~~~~~~~ {date_time} ~~~~~~~~~~~~')
logger.info(f'~~~~~~~~~~~~ Host: {host_name} ~~~~~~~~~~~~')
logger.info(f'Maximum number of parallel processes: {Nmax_proc}')
start_time0 = time.time()

# some general info
## work dir
config_paths = config['Paths']
config_dir = config_paths.get('config_dir')
out_dir = config_paths.get('out_dir')
out_dir = os.path.join(out_dir, run_tag)
if not os.path.exists(out_dir):
    os.mkdir(out_dir)
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
    logger.info(f'Running log from external codes will be saved to {log_dir}')
else:
    log_dir = None

## galaxy info
config_gal = config['GalInfo']
gal_file = config_gal.get('cata_file')
gal_position_type = config_gal.get('position_type')
gal_mag_min_cut = config_gal.getfloat('mag_min_cut')
gal_id_name = config_gal.get('id_name')
gal_detection_mag_name = config_gal.get('detection_mag_name')
gal_mag_name_list = [x.strip() for x in config_gal.get('mag_name_list').split(',')]
gal_RaDec_names = [x.strip() for x in config_gal.get('RaDec_names').split(',')]
gal_shape_names = [x.strip() for x in config_gal.get('shape_names').split(',')]

## star info
config_star = config['StarInfo']
star_file = config_star.get('cata_file')
if star_file:
    contain_stars = 'yes'
    star_area = config_star.getfloat('cata_area')
    star_position_type = config_star.get('position_type')
    star_mag_min_cut = config_star.getfloat('mag_min_cut')
    star_id_name = config_star.get('id_name')
    star_detection_mag_name = config_star.get('detection_mag_name')
    star_mag_name_list = [x.strip() for x in config_star.get('mag_name_list').split(',')]
    if star_position_type == 'true':
        star_RaDec_names = [x.strip() for x in config_star.get('RaDec_names').split(',')]
    else:
        star_RaDec_names = None
else:
    contain_stars = 'no'

## noise info
config_noise = config['NoiseInfo']
noise_file = config_noise.get('cata_file')
noise_psf_basenames = [x.strip() for x in config_noise.get('noise_psf_basenames').split(',')]

## ImSim
config_imsim = config['ImSim']
survey = config_imsim.get('survey')
N_tiles = config_imsim.getint('N_tiles')
gal_rotation_angles = [float(i_r.strip()) for i_r in config_imsim.get('gal_rotation_angles').split(',')]
PSF_map = config_imsim.getboolean('PSF_map')
rng_seed = config_imsim.getint('rng_seed')
mag_zero = config_imsim.getfloat('mag_zero')
g_cosmic = [float(i_g.strip()) for i_g in config_imsim.get('g_cosmic').split(',')]
bands = [x.strip() for x in config_imsim.get('band_list').split(',')]
pixel_scale_list = [float(i_p.strip()) for i_p in config_imsim.get('pixel_scale_list').split(',')]
image_chips = [bool(distutils.util.strtobool(x.strip())) for x in config_imsim.get('image_chips').split(',')]
image_PSF = [bool(distutils.util.strtobool(x.strip())) for x in config_imsim.get('image_PSF').split(',')]

# save all the basic info
outfile_tmp = os.path.join(out_dir, f'basic_info.txt')
f = open(outfile_tmp, 'w')
print(f'# Some basic info about the simulated images in this directory\n\
run_tag            =   {run_tag}\n\
survey             =   {survey}\n\
g_cosmic           =   {g_cosmic[0]:.2f}, {g_cosmic[1]:.2f}\n\
gal_position_type  =   {gal_position_type}\n\
contain_stars      =   {contain_stars}', file=f)
if star_file:
    print(f'star_position_type =   {star_position_type}\n', file=f)
f.close()
logger.info(f'Basic info about the setups saved to {outfile_tmp}')

# load noise info
noise_info = LoadCata.NoiseInfo(noise_file, bands, noise_psf_basenames)

# tile label list
tile_labels = noise_info['label'].to_list()
tile_labels = tile_labels[:N_tiles]

# 1: simulate images
if ('1' in args.taskIDs) or ('all' in args.taskIDs):

    logger.info('====== Task 1: simulate images === started ======')
    start_time = time.time()

    ## load galaxy info
    gals_info = LoadCata.GalInfo(gal_file, bands,
                        gal_id_name, gal_detection_mag_name, gal_mag_name_list,
                        gal_RaDec_names, 
                        gal_shape_names,
                        rng_seed=rng_seed, mag_min_cut=gal_mag_min_cut)

    ## load star info
    if star_file:
        stars_info = LoadCata.StarInfo(star_file, bands,
                            star_id_name, star_detection_mag_name, star_mag_name_list,
                            RaDec_names=star_RaDec_names,
                            mag_min_cut=star_mag_min_cut)
    else:
        star_area = None
        stars_info = None
        star_position_type = None

    ## running
    out_dir_tmp = os.path.join(ima_dir, 'original')
    if not os.path.exists(out_dir_tmp):
        os.mkdir(out_dir_tmp)
    ImSim.RunParallel_PSFNoisySkyImages(survey, out_dir_tmp, rng_seed, mag_zero,
                                            Nmax_proc, 
                                            N_tiles, bands, pixel_scale_list,
                                            noise_info,
                                            gals_info, gal_rotation_angles=gal_rotation_angles, g_cosmic=g_cosmic, gal_position_type=gal_position_type,
                                            stars_area=star_area, stars_info=stars_info, star_position_type=star_position_type,
                                            PSF_map=PSF_map, N_PSF=100, sep_PSF=120,
                                            image_chips=image_chips, image_PSF=image_PSF)

    logger.info(f'====== Task 1: simulate images === finished in {time.time()-start_time} s ======')

# 2: swarp images
if ('2' in args.taskIDs) or ('all' in args.taskIDs):

    logger.info('====== Task 2: swarp images === started ======')
    start_time = time.time()

    # config
    config_swarp = config['SWarp']
    swarp_path = config_swarp.get('swarp_path')
    swarp_psf_map = config_swarp.getboolean('swarp_psf_map')
    swarp_configs = [x.strip() for x in config_swarp.get('swarp_configs').split(',')]
    swarp_bands_group = re.findall(r'\[([^]]+)', config_swarp.get('swarp_bands_group'))
    swarp_labels = [x.strip() for x in config_swarp.get('swarp_labels').split(',')]
    only_resamples = [bool(distutils.util.strtobool(x.strip())) for x in config_swarp.get('only_resamples').split(',')]
    swarp_clean_up_levels = [int(x.strip()) for x in config_swarp.get('clean_up_levels').split(',')] 

    # I/O
    in_dir_tmp = os.path.join(ima_dir, 'original')
    if swarp_psf_map:
        in_dir_psf_tmp = os.path.join(in_dir_tmp, 'psf_map')
    if running_log:
        swarp_log_dir = os.path.join(log_dir, 'SWarp')
        if not os.path.exists(swarp_log_dir):
            os.mkdir(swarp_log_dir)
    else:
        swarp_log_dir = None

    # work pool
    ## Swarp is I/O intense
    N_swarp = int(Nmax_proc/12.)
    if N_swarp < 1:
        N_swarp = 1
    logger.info(f'Number of processes for SWarp: {N_swarp}')
    work_pool = mp.Pool(processes=N_swarp) 
    proc_list = []
    for i_group, swarp_config in enumerate(swarp_configs):

        swarp_bands = swarp_bands_group[i_group]
        swarp_bands = [x.strip() for x in swarp_bands.split(',')]

        only_resample = only_resamples[i_group]

        swarp_config_pathfile = os.path.join(config_dir, swarp_config)

        swarp_clean_up_level = swarp_clean_up_levels[i_group]

        swarp_label = swarp_labels[i_group]
        out_dir_tmp = os.path.join(ima_dir, swarp_label)
        if not os.path.exists(out_dir_tmp):
            os.mkdir(out_dir_tmp)
        if swarp_psf_map:
            out_dir_psf_tmp =  os.path.join(out_dir_tmp, 'psf_map')
            if not os.path.exists(out_dir_psf_tmp):
                os.mkdir(out_dir_psf_tmp)

        for tile_label in tile_labels:

            for band in swarp_bands:

                for gal_rotation_angle in gal_rotation_angles:

                    if 'simple' in survey:
                        image_in = os.path.join(in_dir_tmp, f'tile{tile_label}_band{band}_rot{gal_rotation_angle:.0f}.fits')
                    elif survey.lower() == 'kids_sameexpo':
                        image_in = [os.path.join(in_dir_tmp, f'tile{tile_label}_band{band}_rot{gal_rotation_angle:.0f}_expo{i_expo}.fits') for i_expo in range(5)]
                    
                    # running
                    image_out = os.path.join(out_dir_tmp, f'tile{tile_label}_band{band}_rot{gal_rotation_angle:.0f}.fits')
                    proc = work_pool.apply_async(func=Astromatic.SwarpImage, 
                                            args=(image_in, swarp_config_pathfile, image_out, only_resample, running_log, swarp_log_dir, swarp_path, swarp_clean_up_level))
                    proc_list.append(proc)

                    # psf map
                    if swarp_psf_map:
                        try:
                            image_in_psf = os.path.join(in_dir_psf_tmp, os.path.basename(image_in))
                        except TypeError:
                            image_in_psf = [os.path.join(in_dir_psf_tmp, os.path.basename(image_in_tmp)) for image_in_tmp in image_in]
                        ## running
                        image_out = os.path.join(out_dir_psf_tmp, f'tile{tile_label}_band{band}_rot{gal_rotation_angle:.0f}.fits')
                        proc = work_pool.apply_async(func=Astromatic.SwarpImage, 
                                                args=(image_in_psf, swarp_config_pathfile, image_out, only_resample, running_log, swarp_log_dir, swarp_path, swarp_clean_up_level))
                        proc_list.append(proc)

    work_pool.close()
    work_pool.join()
    # check for any errors during run
    for proc in proc_list:
        proc.get()

    logger.info(f'====== Task 2: swarp images === finished in {time.time()-start_time} s ======')

# 3: detect objects
if ('3' in args.taskIDs) or ('all' in args.taskIDs):

    logger.info('====== Task 3: detect objects === started ======')
    start_time = time.time()

    # config
    config_sex = config['SExtractor']
    sex_path = config_sex.get('sex_path')
    detection_band = config_sex.get('detection_band')
    image_label = config_sex.get('image_label')
    sex_mask_stars = config_sex.getboolean('mask_stars')
    sex_cross_match = config_sex.getboolean('cross_match')
    sex_config = os.path.join(config_dir, config_sex.get('sex_config'))
    sex_param = os.path.join(config_dir, config_sex.get('sex_param'))
    sex_filter = os.path.join(config_dir, config_sex.get('sex_filter'))
    sex_starNNW = os.path.join(config_dir, config_sex.get('sex_starNNW'))
    sex_clean_up_level =  config_sex.getint('clean_up_level')

    # basic info
    index_tmp = [i for i,x in enumerate(bands) if x == detection_band][0]
    pixel_scale = pixel_scale_list[index_tmp]
    logger.info(f'Band for detection: {detection_band}')
    logger.info(f'pixel scale: {pixel_scale}')
    logger.info(f'image type: {image_label}')

    # I/O
    ori_ima_dir = os.path.join(ima_dir, 'original')
    in_dir_tmp = os.path.join(ima_dir, image_label)
    out_dir_tmp = os.path.join(cata_dir, 'SExtractor')
    if not os.path.exists(out_dir_tmp):
        os.mkdir(out_dir_tmp)
    if running_log:
        sex_log_dir = os.path.join(log_dir, 'SExtractor')
        if not os.path.exists(sex_log_dir):
            os.mkdir(sex_log_dir)
    else:
        sex_log_dir = None
    if sex_cross_match:
        out_dir_cross = os.path.join(cata_dir, 'CrossMatch')
        if not os.path.exists(out_dir_cross):
            os.mkdir(out_dir_cross)

    # seeing info
    SeeingFWHM_list = noise_info[f'seeing_{detection_band}'].to_list()
    SeeingFWHM_list = SeeingFWHM_list[:N_tiles]

    # work pool
    ## SExtractor is very fast
    N_sex = int(Nmax_proc/4.)
    if N_sex < 1:
        N_sex = 1
    logger.info(f'Number of processes for SExtractor: {N_sex}')
    work_pool = mp.Pool(processes=N_sex) 
    proc_list = []
    for i_tile, tile_label in enumerate(tile_labels):
        SeeingFWHM = SeeingFWHM_list[i_tile]

        for gal_rotation_angle in gal_rotation_angles:

            ### image
            ImageFile1 = os.path.join(in_dir_tmp, f'tile{tile_label}_band{detection_band}_rot{gal_rotation_angle:.0f}.fits')
            WeightFile = ImageFile1.replace('.fits', '.weight.fits')
            if not os.path.isfile(WeightFile):
                WeightFile = None

            ### catalogue
            CatalogueFile = os.path.join(out_dir_tmp, os.path.basename(ImageFile1).replace('.fits', '.sex'))

            ### running
            proc = work_pool.apply_async(func=Astromatic.SExtractorCatalogue, 
                                    args=(ImageFile1, CatalogueFile, pixel_scale, SeeingFWHM, 
                                            None, WeightFile,
                                            running_log, sex_log_dir,
                                            sex_path, sex_config, sex_param, 
                                            sex_filter, sex_starNNW,
                                            mag_zero,
                                            sex_clean_up_level))
            proc_list.append(proc)

    work_pool.close()
    work_pool.join()
    # check for any errors during run
    for proc in proc_list:
        proc.get()

    # Mark stars
    if sex_mask_stars:

        logger.info('Mask Stars based on input info.')

        for tile_label in tile_labels:

            input_cata = pd.read_feather(os.path.join(ori_ima_dir, f'stars_info_tile{tile_label}.feather'))

            for gal_rotation_angle in gal_rotation_angles:

                detec_file = os.path.join(out_dir_tmp, f'tile{tile_label}_band{detection_band}_rot{gal_rotation_angle:.0f}.feather')
                detec_cata = pd.read_feather(detec_file)

                id_list = ['index', 'NUMBER']
                position_list = [['RA', 'DEC'], ['X_WORLD', 'Y_WORLD']]
                mag_list = [detection_band, 'MAG_AUTO']

                matched_cata, _, _ = CrossMatch.run_position2id(input_cata, detec_cata, id_list, position_list, mag_list, 
                                                    outDir=None, basename=None, save_matched=False, save_false=False, save_missed=False,
                                                    dmag_max=0.3, r_max=0.5/3600., k=4)

                mask_stars = detec_cata['NUMBER'].isin(matched_cata['id_detec'])
                detec_cata.loc[mask_stars, 'perfect_flag_star'] = 1
                detec_cata.loc[~mask_stars, 'perfect_flag_star'] = 0
                detec_cata.to_feather(detec_file)

    # cross-match with the input catalogue
    if sex_cross_match:

        logger.info('cross-match with the input catalogue.')

        # config
        config_cross = config['CrossMatch']
        mag_faint_cut = config_cross.getfloat('mag_faint_cut')
        save_matched = config_cross.getboolean('save_matched')
        save_false = config_cross.getboolean('save_false')
        save_missed = config_cross.getboolean('save_missed')
        dmag_max = config_cross.getfloat('dmag_max')
        r_max = config_cross.getfloat('r_max')

        for tile_label in tile_labels:

            input_cata = pd.read_feather(os.path.join(ori_ima_dir, f'gals_info_tile{tile_label}.feather'))

            ## magnitude pre-selection
            input_cata = input_cata[input_cata[detection_band]<=mag_faint_cut]
            input_cata.reset_index(drop=True, inplace=True)

            for gal_rotation_angle in gal_rotation_angles:

                detec_file = os.path.join(out_dir_tmp, f'tile{tile_label}_band{detection_band}_rot{gal_rotation_angle:.0f}.feather')
                detec_cata = pd.read_feather(detec_file)
    
                ## select galaxies
                try:
                    mask_gals = (detec_cata['perfect_flag_star'] == 0)
                    detec_cata = detec_cata[mask_gals]
                    detec_cata.reset_index(drop=True, inplace=True)
                except KeyError:
                    detec_cata = detec_cata

                # output info
                basename_cross =  f'tile{tile_label}_rot{gal_rotation_angle:.0f}'

                # cross-match
                id_list = ['index', 'NUMBER']
                position_list = [['RA', 'DEC'], ['X_WORLD', 'Y_WORLD']]
                mag_list = [detection_band, 'MAG_AUTO']
                _, _, _ = CrossMatch.run_position2id(input_cata, detec_cata, id_list, position_list, mag_list, 
                                                    outDir=out_dir_cross, basename=basename_cross, save_matched=save_matched, save_false=save_false, save_missed=save_missed,
                                                    dmag_max=dmag_max, r_max=r_max/3600., k=4)

    logger.info(f'====== Task 3: detect objects === finished in {time.time()-start_time} s ======')

# 4: measure photometry
if ('4' in args.taskIDs) or ('all' in args.taskIDs):

    logger.info('====== Task 4: measure photometry === started ======')
    start_time = time.time()

    # config
    config_MP = config['MeasurePhotometry']
    MP_method = config_MP.get('method')
    MP_detec_band = config_MP.get('detection_band')
    MP_bands = [x.strip() for x in config_MP.get('band_list').split(',')]
    MP_labels = [x.strip() for x in config_MP.get('image_label_list').split(',')]

    if MP_method.lower() == 'gaap':
        logger.info('Use GAaP for photometry measurement.')

        # config
        config_gaap = config['GAaP']
        gaap_dir = config_gaap.get('gaap_dir')
        min_aper = config_gaap.getfloat('min_aper')
        max_aper = config_gaap.getfloat('max_aper')
        use_psf_map = config_gaap.getboolean('use_psf_map')
        star_mag_cut = [float(i_m) for i_m in config_gaap.get('star_mag_cut').split(',')]
        clean_up_level = config_gaap.getint('clean_up_level')

        # I/O
        in_cata_dir_tmp = os.path.join(cata_dir, 'SExtractor')
        ori_ima_dir = os.path.join(ima_dir, 'original')
        out_dir_tmp = os.path.join(cata_dir, 'GAaP')
        if not os.path.exists(out_dir_tmp):
            os.mkdir(out_dir_tmp)
        if running_log:
            gaap_log_dir = os.path.join(log_dir, 'GAaP')
            if not os.path.exists(gaap_log_dir):
                os.mkdir(gaap_log_dir)
        else:
            gaap_log_dir = None

        # Initialise the GAaP wrapper
        gaap = GAaP.GAaPwrapper(gaap_dir, out_dir_tmp, 
                star_mag_cut=star_mag_cut,
                mag_zero=mag_zero,
                min_aper=min_aper, max_aper=max_aper, 
                running_log=running_log, log_dir=gaap_log_dir,
                clean_up_level=clean_up_level)

        # each tile for each process
        #   no parallel between bands
        N_gaap = Nmax_proc
        if (len(tile_labels)) < N_gaap:
            N_gaap = len(tile_labels)
        logger.info(f'Number of processes for GAaP: {N_PROC}')
        work_pool = mp.Pool(processes=N_PROC) 
        proc_list = []
        for tile_label in tile_labels:

            # star info for psf estimation
            if not use_psf_map:
                star_info_file = os.path.join(ori_ima_dir, f'stars_info_tile{tile_label}.feather')
                star_info = pd.read_feather(star_info_file)
            else:
                star_info = None

            for gal_rotation_angle in gal_rotation_angles:

                # detection file
                SKYcataFile = os.path.join(in_cata_dir_tmp, f'tile{tile_label}_band{MP_detec_band}_rot{gal_rotation_angle:.0f}.feather')

                # image files
                SKYimaFile_list = []
                SKYweiFile_list = []
                if use_psf_map:
                    PSFimaFile_list = []
                else:
                    PSFimaFile_list = None
                for i_band, band in enumerate(MP_bands):
 
                    in_ima_dir_tmp = os.path.join(ima_dir, MP_labels[i_band])
                    SKYimaFile_tmp = os.path.join(in_ima_dir_tmp, f'tile{tile_label}_band{band}_rot{gal_rotation_angle:.0f}.fits')
                    SKYimaFile_list.append(SKYimaFile_tmp)

                    # weight image                        
                    SKYweiFile_tmp = SKYimaFile_tmp.replace('.fits', '.weight.fits')
                    if not os.path.isfile(SKYweiFile_tmp):
                        SKYweiFile_tmp = None
                    SKYweiFile_list.append(SKYweiFile_tmp)

                    # psf map 
                    if use_psf_map:
                        in_psf_dir_tmp = os.path.join(in_ima_dir_tmp, 'psf_map')
                        PSFimaFile_tmp = os.path.join(in_psf_dir_tmp, os.path.basename(SKYimaFile_tmp))
                        PSFimaFile_list.append(PSFimaFile_tmp)

                # output
                FinalFile = os.path.join(out_dir_tmp, f'tile{tile_label}_rot{gal_rotation_angle:.0f}.feather')

                proc = work_pool.apply_async(func=gaap.RunSingleTile, 
                                args=(tile_label, FinalFile, SKYcataFile, MP_bands, SKYimaFile_list, SKYweiFile_list, PSFimaFile_list, star_info))
                proc_list.append(proc)

        work_pool.close()
        work_pool.join()
        # check for any errors during run
        for proc in proc_list:
            proc.get()

    else:
        raise Exception(f'Unsupported photometry measurement method {MP_method}!')

    logger.info(f'====== Task 4: measure photometry === finished in {time.time()-start_time} s ======')

# 6: measure galaxy shapes
if ('6' in args.taskIDs) or ('all' in args.taskIDs):

    logger.info('====== Task 6: measure galaxy shapes === started ======')
    start_time = time.time()

    # config
    config_MS = config['MeasureShape']
    MS_method = config_MS.get('method')
    MS_detec_band = config_MS.get('detection_band')
    MS_bands = [x.strip() for x in config_MS.get('band_list').split(',')]
    MS_labels = [x.strip() for x in config_MS.get('image_label_list').split(',')]

    if MS_method.lower() == 'lensfit':
        logger.info('Use lensfit for shape measurement.')

        # config
        config_lensfit = config['lensfit']
        lensfit_dir = config_lensfit.get('lensfit_dir')
        python2_env = config_lensfit.get('python2_env')
        clean_up_level = config_lensfit.getint('clean_up_level')
        postage_size = config_lensfit.get('postage_size')
        start_exposure = config_lensfit.get('start_exposure')
        end_exposure = config_lensfit.get('end_exposure')
        start_mag = config_lensfit.get('start_mag')
        end_mag = config_lensfit.get('end_mag')
        PSF_OVERSAMPLING = config_lensfit.get('PSF_OVERSAMPLING')
        PECUT = config_lensfit.get('PECUT')
        PRCUT = config_lensfit.get('PRCUT')
        LCUT = config_lensfit.get('LCUT')
        CAMERA = config_lensfit.get('CAMERA').upper()

        # I/O
        in_cata_dir_tmp = os.path.join(cata_dir, 'SExtractor')
        ori_ima_dir = os.path.join(ima_dir, 'original') # lensfit use individual exposures
        out_dir_tmp = os.path.join(cata_dir, 'lensfit')
        if not os.path.exists(out_dir_tmp):
            os.mkdir(out_dir_tmp)
        tmp_dir_lensfit = os.path.join(cata_dir, 'tmp_lensfit')
        if not os.path.exists(tmp_dir_lensfit):
            os.mkdir(tmp_dir_lensfit)
        if running_log:
            log_dir_lensfit = os.path.join(log_dir, 'lensfit')
            if not os.path.exists(log_dir_lensfit):
                os.mkdir(log_dir_lensfit)
        else:
            log_dir_lensfit = None

        # prepare input file
        input_file_lensfit = os.path.join(tmp_dir_lensfit, f'lensfit_input.asc')
        f = open(input_file_lensfit, 'w')
        if CAMERA == 'KIDS':
            N_expo = 5
            N_chips = 32
        for i_expo in range(N_expo):
            for i_chip in range(N_chips):
                print(f'expo{i_expo}_chip{i_chip}', file=f)
        f.close()
        logger.info(f'input info saved to {input_file_lensfit}')

        # work pool
        N_lensfit = int(Nmax_proc/12)
        if N_lensfit < 1:
            N_lensfit = 1
        logger.info(f'Number of processes for lensfit: {N_lensfit}')
        logger.info(f'  NOTE: each processes of lensfit takes multiple cores')
        logger.info(f'          depending on your lensfit compilation')
        work_pool = mp.Pool(processes=N_lensfit) 
        proc_list = []
        for i_band, band in enumerate(MS_bands):

            WAVEBAND = band.upper()
            MS_label = MS_labels[i_band]
            logging.info(f'Measure shapes in {band} band using {MS_label} images')

            for tile_label in tile_labels:

                # data info
                weight_dir = None
                psf_dir = os.path.join(ori_ima_dir, f'psf_tile{tile_label}_band{band}')
                head_dir = os.path.join(ori_ima_dir, f'chips_tile{tile_label}_band{band}_head')
                head_key = True

                for gal_rotation_angle in gal_rotation_angles:

                    # tmp directory
                    tmp_dir = os.path.join(tmp_dir_lensfit, f'tile{tile_label}_band{band}_rot{gal_rotation_angle:.0f}')
                    if not os.path.exists(tmp_dir):
                        os.mkdir(tmp_dir)

                    # prepare detection catalogue
                    CatalogueFile = os.path.join(in_cata_dir_tmp, f'tile{tile_label}_band{detection_band}_rot{gal_rotation_angle:.0f}.feather')
                    cata_ori = pd.read_feather(CatalogueFile)
                    ## select galaxies
                    try:
                        mask_gals = (cata_ori['perfect_flag_star'] == 0)
                        cata_ori = cata_ori[mask_gals]
                        cata_ori.reset_index(drop=True, inplace=True)
                    except KeyError:
                        cata_ori = cata_ori
                    ## desired info    
                    reduced_data = np.array([cata_ori['X_WORLD'], cata_ori['Y_WORLD'], cata_ori['MAG_AUTO'], cata_ori['NUMBER']]).T
                    ## save
                    input_catalog_lensfit = os.path.join(tmp_dir, 'catalog.asc')
                    np.savetxt(input_catalog_lensfit, reduced_data)
                    logger.info(f'catalogue info saved to {input_catalog_lensfit}')

                    # data info
                    chip_dir = os.path.join(ori_ima_dir, f'chips_tile{tile_label}_band{detection_band}_rot{gal_rotation_angle:.0f}')

                    # update head for lensfit
                    if head_key:
                        LensFit.LensfitShape_head(chip_dir, head_dir)
                        head_key = False
                    else:
                        LensFit.LensfitShape_head(chip_dir, None)

                    # run
                    output_file = f'tile{tile_label}_band{detection_band}_rot{gal_rotation_angle:.0f}.feather'
                    proc = work_pool.apply_async(func=LensFit.LensfitShape, 
                                args=(lensfit_dir, python2_env,
                                        input_catalog_lensfit, input_file_lensfit, chip_dir, psf_dir, head_dir, weight_dir,
                                        PSF_OVERSAMPLING, PECUT, PRCUT, LCUT, WAVEBAND, CAMERA,
                                        postage_size, start_exposure, end_exposure, start_mag, end_mag,
                                        output_file, out_dir_tmp, tmp_dir,
                                        running_log, log_dir_lensfit, 
                                        clean_up_level))
                    proc_list.append(proc)

        work_pool.close()
        work_pool.join()
        # check for any errors during run
        for proc in proc_list:
            proc.get()

    else:
        raise Exception(f'Unsupported shape measurement method {method_shape}!')

    logger.info(f'====== Task 6: measure galaxy shapes === finished in {time.time()-start_time} s ======')

# # 7: create a combined catalogue
# if ('7' in args.taskIDs) or ('all' in args.taskIDs):

#     logger.info('====== Task 7: create a combined catalogue === started ======')
#     start_time = time.time()
    
#     # config
#     config_combine = config['CombineCata']
#     combine_cata_list = [x.strip() for x in config_combine.get('cata_list').split(',')]
#     logger.info(f'catalogues being combined: {combine_cata_list}')


#     if 'input' is not 


logger.info(f'~~~~~~~~~~~~ {__version__} finished ~~~~~~~~~~~~')
logger.info(f'~~~~~~~~~~~~ total running time {time.time()-start_time0} s ~~~~~~~~~~~~')
logger.info(f'~~~~~~~~~~~~ All outputs saved in {out_dir}')
logger.info(f'~~~~~~~~~~~~ Enjoy the science ~~~~~~~~~~~~')
