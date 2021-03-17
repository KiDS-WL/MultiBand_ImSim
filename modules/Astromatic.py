# -*- coding: utf-8 -*-
# @Author: lshuns
# @Date:   2020-08-17 14:26:07
# @Last modified by:   lshuns
# @Last modified time: 2021-03-17, 12:08:05

### Wrapper for astromatic codes

import os
import sys
import logging
import subprocess

import numpy as np
import pandas as pd

from astropy.io import fits

logger = logging.getLogger(__name__)

def SwarpImage(image_in, swarp_config_file,
                    image_out,
                    only_resample=False,
                    running_log=True, log_dir=None,
                    swarp_path='swarp',
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
    except KeyError:
        os.remove(image_out)
        logger.info("Remove existing images.")

    # running info
    if running_log:
        basename = os.path.basename(image_out)
        outLog = open(os.path.join(log_dir, basename.replace('.fits', '.log')), "w")
        errLog = open(os.path.join(log_dir, basename.replace('.fits', '.err.log')), "w")
    else:
        outLog = subprocess.PIPE
        errLog = subprocess.STDOUT

    # resampling out directory
    RESAMPLE_DIR = os.path.dirname(image_out)

    # build command
    cmd = [swarp_path]
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

    logger.info(f'Config file {swarp_config_file}')

    # run
    proc = subprocess.run(cmd, stdout=outLog, stderr=errLog)

    if running_log:
        outLog.close()
        errLog.close()

    # rename for only resampling
    if only_resample:
        basename = os.path.basename(image_out)
        os.rename(os.path.join(RESAMPLE_DIR, basename.replace('.fits', '.resamp.fits')), image_out)
        logger.info(f"Swarp resampled image saved as {image_out}")
    else:
        logger.info(f"Swarp coadded image saved as {image_out}")

    # mark success to the header
    with fits.open(image_out, mode='update') as hdul:
        head_tmp = hdul[0].header
        ## update info
        head_tmp['flag_sim'] = 2

    if clean_up_level:
        try:
            os.remove(image_in)
        except TypeError:
            for image_tmp in image_in:
                os.remove(image_tmp)
        logger.info('Original images are removed.')

    logger.info('SWarp finished.')
    return 0

def SExtractorCatalogue(ImageFile1, CatalogueFile, pixel_scale, SeeingFWHM=1.0,
                        ImageFile2=None, WeightFile=None,
                        running_log=True, log_dir=None,
                        SexPath='sex', ConfigFile='../config/kids_sims.sex', ParamFile='../config/sex_image.param',
                        FilterName='../config/default.conv', StarnnwName='../config/default.nnw',
                        CHECKIMAGE_TYPE=None,
                        mag_zero=30.,
                        clean_up_level=0):
    '''
    SExtractor to build catalogue from image.
        ImageFile2: only required for the dual-mode.

        Clean up level:
            0: none
            1: rm original images
            2: and .sex files
    '''

    # check existence
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
        outLog = open(os.path.join(log_dir, basename.replace('.sex', '.log')), "w")
        errLog = open(os.path.join(log_dir, basename.replace('.sex', '.err.log')), "w")
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
        cmd = [SexPath, ImageFile1, '-CATALOG_NAME', CatalogueFile, '-c', ConfigFile, '-PARAMETERS_NAME', ParamFile,
                '-PIXEL_SCALE', str(pixel_scale), '-SEEING_FWHM', str(SeeingFWHM),
                '-FILTER_NAME', FilterName, '-STARNNW_NAME', StarnnwName,
                '-CHECKIMAGE_TYPE', CHECKIMAGE_TYPE, '-CHECKIMAGE_NAME', CHECKIMAGE_NAME,
                '-MAG_ZEROPOINT', str(mag_zero)]
        logger.info('Running SEXtractor in single mode...')
    else:
        cmd = [SexPath, ImageFile1, ImageFile2, '-CATALOG_NAME', CatalogueFile, '-c', ConfigFile, '-PARAMETERS_NAME', ParamFile,
                '-PIXEL_SCALE', str(pixel_scale), '-SEEING_FWHM', str(SeeingFWHM),
                '-FILTER_NAME', FilterName, '-STARNNW_NAME', StarnnwName,
                '-CHECKIMAGE_TYPE', CHECKIMAGE_TYPE, '-CHECKIMAGE_NAME', CHECKIMAGE_NAME,
                '-MAG_ZEROPOINT', str(mag_zero)]
        logger.info('Running SEXtractor in dual mode...')

    if (WeightFile is not None):
        cmd.extend(['-WEIGHT_IMAGE', WeightFile, '-WEIGHT_TYPE', 'MAP_WEIGHT'])
        logger.info(f'Use MAP_WEIGHT')
    logger.info(f'Config file {ConfigFile}')

    # run
    proc = subprocess.run(cmd, stdout=outLog, stderr=errLog)
    logger.info("SEXtractor catalogue produced as {:}".format(CatalogueFile))

    if running_log:
        outLog.close()
        errLog.close()

    # feather file
    logger.debug('Build feather file for easy use ...')
    col_name = np.loadtxt(ParamFile, dtype=str)
    values = np.loadtxt(CatalogueFile)
    df = pd.DataFrame(data=values, columns=col_name)
    df = df.astype({'NUMBER': 'int', 'FLAGS': 'int'})
    df.to_feather(file_feather)
    logger.info(f'Easy-use feather file built as {file_feather}')

    if (clean_up_level >= 1):
        os.remove(ImageFile1)
        if ImageFile2 is not None:
            os.remove(ImageFile2)
        logger.info('Images are removed.')

    if (clean_up_level >= 2):
        os.remove(CatalogueFile)
        logger.info('Original .sex file is removed.')

    logger.info('SEXtractor finished.')
