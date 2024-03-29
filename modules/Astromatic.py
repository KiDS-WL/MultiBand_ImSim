# -*- coding: utf-8 -*-
# @Author: lshuns
# @Date:   2020-08-17 14:26:07
# @Last modified by:   ssli
# @Last modified time: 2021-05-18, 14:29:16

### Wrapper for astromatic codes

import os
import sys
import glob
import shutil
import logging
import pathlib
import subprocess

import numpy as np
import pandas as pd

from astropy.io import fits

logger = logging.getLogger(__name__)

def SwarpImage(image_in, swarp_config_file,
                    image_out, RESAMPLE_DIR,
                    only_resample=False, contain_wei_ima=True,
                    running_log=True, log_dir=None,
                    swarp_path='swarp', NTHREADS=0,
                    clean_up_level=0):
    """
    SWarp for coadding or resampling.
        only_resample: set to True
                        if only resampling but no coadding.
        Clean up level:
            0: none
            1: rm original images
    """

    # first check if already exist
    try:
        with fits.open(image_out) as hdul:
            head_tmp = hdul[0].header
        flag_sim = head_tmp['flag_sim']
        if flag_sim >= 2:
            logger.info(f"{image_out} already exist.")
            return 1
    except FileNotFoundError:
        pass
    except (KeyError, OSError) as e:
        os.remove(image_out)
        logger.info("Removed existing images.")

    # remove any potential intermediate files
    file_list_tmp = glob.glob(os.path.join(os.path.dirname(image_out), '*.resamp.*'))
    if file_list_tmp:
        for file_tmp in file_list_tmp:
            os.remove(file_tmp)
        logger.info('Removed intermediate files from previous run.')

    # running info
    if running_log:
        basename = os.path.basename(image_out)
        outLog = open(os.path.join(log_dir, basename.replace('.fits', '.log')), "w")
        errLog = open(os.path.join(log_dir, basename.replace('.fits', '.err.log')), "w")
    else:
        outLog = subprocess.PIPE
        errLog = subprocess.STDOUT

    # resampling out directory needs to be unique for each run
    ### to avoid cross-talking between different swarp runs
    if (os.path.exists(RESAMPLE_DIR)) and (any(os.scandir(RESAMPLE_DIR))):
        raise Exception(f'{RESAMPLE_DIR} is not empty, cannot use as RESAMPLE_DIR')
    pathlib.Path(RESAMPLE_DIR).mkdir(parents=True, exist_ok=True)

    # build command
    cmd = [swarp_path]
    if isinstance(image_in, str):
        cmd.append(image_in)
    else:
        cmd.extend(image_in) # in case image_in containing multiple inputs
    if only_resample:
        cmd.extend(['-c', swarp_config_file,
            '-RESAMPLE_DIR', RESAMPLE_DIR])
        logger.info('Running SWarp for only resampling...')
    else:
        cmd.extend(['-c', swarp_config_file,
            '-RESAMPLE_DIR', RESAMPLE_DIR,
            '-IMAGEOUT_NAME', image_out, '-WEIGHTOUT_NAME', image_out.replace('.fits', '.weight.fits')])
        logger.info('Running SWarp for coadding...')

    # contain weight images or not
    if contain_wei_ima:
        cmd.extend(['-WEIGHT_TYPE', 'MAP_WEIGHT'])
    else:
        cmd.extend(['-WEIGHT_TYPE', 'NONE'])

    # how many threads
    cmd.extend(['-NTHREADS', str(NTHREADS)])

    logger.info(f'Config file {swarp_config_file}')
    logger.info(f'Number of threads: {NTHREADS}')

    # run
    proc = subprocess.run(cmd, stdout=outLog, stderr=errLog)

    if running_log:
        outLog.close()
        errLog.close()

    # mv or rename for only resampling
    if only_resample:
        basename = os.path.basename(image_in)
        os.rename(os.path.join(RESAMPLE_DIR, basename.replace('.fits', '.resamp.fits')), image_out)

    # mark success to the header
    try:
        with fits.open(image_out, mode='update') as hdul:
            head_tmp = hdul[0].header
            ## update info
            head_tmp['flag_sim'] = 2
    except FileNotFoundError:
        raise Exception('Cannot find coadded images, something is wrong with Swarp, check out running_log for Swarp!')

    # saving info
    if only_resample:
        logger.info(f"Swarp resampled image saved as {image_out}")
    else:
        logger.info(f"Swarp coadded image saved as {image_out}")

    if clean_up_level:
        try:
            os.remove(image_in)
        except TypeError:
            for image_tmp in image_in:
                os.remove(image_tmp)
        logger.info('Original images are removed.')

    logger.info('SWarp finished.')
    return 0

def SExtractorCatalogue(CatalogueFile, pixel_scale, SeeingFWHM,
                        ImageFile1, WeightFile1=None,
                        ImageFile2=None, WeightFile2=None,
                        running_log=True, log_dir='./',
                        SexPath='sex', ConfigFile='../config/kids_sims_theli.sex', ParamFile='../config/sex_image.param',
                        FilterName='../config/default.conv', StarnnwName='../config/default.nnw',
                        CHECKIMAGE_TYPE=None, build_feather=True,
                        mag_zero=30.,
                        clean_up_level=0):
    '''
    SExtractor to build catalogue from image.
        ImageFile2: only required for the dual-mode.

        Clean up level:
            0: none
            1: .sex files
    '''

    # check existence
    if build_feather:
        file_feather = CatalogueFile.replace('.sex', '.feather')
        if os.path.isfile(file_feather):
            logger.info(f'The final feather catalogue {file_feather} already exists.')
            logger.info(f'End the process.')
            return 1
        if os.path.isfile(CatalogueFile):
            os.remove(CatalogueFile)
            logger.info(f'Remove existing .sex file.')

    # running info
    if running_log:
        basename = os.path.basename(CatalogueFile)
        outLog = open(os.path.join(log_dir, basename+'.log'), "w")
        errLog = open(os.path.join(log_dir, basename+'.err.log'), "w")
    else:
        outLog = subprocess.PIPE
        errLog = subprocess.STDOUT

    # CHECKIMAGE
    if (CHECKIMAGE_TYPE is not None) and (CHECKIMAGE_TYPE.upper() != 'NONE'):
        ima_dir = os.path.dirname(ImageFile1)
        CHECKIMAGE_dir = os.path.join(ima_dir, CHECKIMAGE_TYPE)
        if not os.path.exists(CHECKIMAGE_dir):
            os.mkdir(CHECKIMAGE_dir)
        logger.info(f'CHECKIMAGE will be saved to {CHECKIMAGE_dir}')

        CHECKIMAGE_NAME = os.path.join(CHECKIMAGE_dir, os.path.basename(ImageFile1))
    else:
        CHECKIMAGE_TYPE = 'NONE'
        CHECKIMAGE_NAME = 'doesnotmatter'

    # build command
    if (ImageFile2 is None):
        logger.info('Running SEXtractor in single mode...')
        cmd = [SexPath, ImageFile1, '-CATALOG_NAME', CatalogueFile, '-c', ConfigFile, '-PARAMETERS_NAME', ParamFile,
                '-PIXEL_SCALE', str(pixel_scale), '-SEEING_FWHM', str(SeeingFWHM),
                '-FILTER_NAME', FilterName, '-STARNNW_NAME', StarnnwName,
                '-CHECKIMAGE_TYPE', CHECKIMAGE_TYPE, '-CHECKIMAGE_NAME', CHECKIMAGE_NAME,
                '-MAG_ZEROPOINT', str(mag_zero)]
        if (WeightFile1 is not None):
            cmd.extend(['-WEIGHT_IMAGE', WeightFile1, '-WEIGHT_TYPE', 'MAP_WEIGHT'])
            logger.info(f'Use MAP_WEIGHT')
    else:
        logger.info('Running SEXtractor in dual mode...')
        cmd = [SexPath, f'{ImageFile1},{ImageFile2}', '-CATALOG_NAME', CatalogueFile, '-c', ConfigFile, '-PARAMETERS_NAME', ParamFile,
                '-PIXEL_SCALE', str(pixel_scale), '-SEEING_FWHM', str(SeeingFWHM),
                '-FILTER_NAME', FilterName, '-STARNNW_NAME', StarnnwName,
                '-CHECKIMAGE_TYPE', CHECKIMAGE_TYPE, '-CHECKIMAGE_NAME', CHECKIMAGE_NAME,
                '-MAG_ZEROPOINT', str(mag_zero)]
        if (WeightFile1 is not None) and (WeightFile2 is not None):
            cmd.extend(['-WEIGHT_IMAGE', f'{WeightFile1},{WeightFile2}', '-WEIGHT_TYPE', 'MAP_WEIGHT'])
            logger.info(f'Use MAP_WEIGHT')
        elif ((WeightFile1 is not None) and (WeightFile2 is None)) or ((WeightFile1 is None) and (WeightFile2 is not None)):
            raise Exception('Do not know how to deal with only one weight image in the dual mode...')

    logger.info(f'Config file {ConfigFile}')

    # run
    proc = subprocess.run(cmd, stdout=outLog, stderr=errLog)
    logger.info("SEXtractor catalogue produced as {:}".format(CatalogueFile))

    if running_log:
        outLog.close()
        errLog.close()

    # feather file
    if build_feather:
        logger.debug('Build feather file for easy use ...')
        col_name = np.loadtxt(ParamFile, dtype=str)
        values = np.loadtxt(CatalogueFile)
        df = pd.DataFrame(data=values, columns=col_name)
        df = df.astype({'NUMBER': 'int', 'FLAGS': 'int'})
        tmp_file_feather = file_feather + '_tmp'
        df.to_feather(tmp_file_feather)
        os.rename(tmp_file_feather, file_feather)
        logger.info(f'Easy-use feather file built as {file_feather}')

    if (clean_up_level >= 1):
        os.remove(CatalogueFile)
        logger.info('Original .sex file is removed.')

    logger.info('SEXtractor finished.')
