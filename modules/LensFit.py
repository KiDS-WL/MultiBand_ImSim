# -*- coding: utf-8 -*-
# @Author: lshuns
# @Date:   2020-12-03 16:16:21
# @Last modified by:   lshuns
# @Last modified time: 2021-03-11, 23:21:58

### Wrapper for lensfit code

import os
import re
import sys
import glob
import shutil
import logging
import subprocess

import numpy as np
import pandas as pd
import multiprocessing as mp

from pathlib import Path
from astropy.io import fits

logger = logging.getLogger(__name__)

def LensfitShape_head(chip_dir, head_dir=None):
    """
    adjust image head for Lensfit
    """

    # images
    image_list = glob.glob(chip_dir+'/*.fits')

    for image_file in image_list:

        ## some info from name
        id_exposure, i_chip = re.findall(r"\d+", os.path.basename(image_file))

        ## load chips head
        with fits.open(image_file, mode='update') as hdul:
            head_tmp = hdul[0].header

            ## update info
            head_tmp['IMAGEID'] = int(i_chip)+1
            head_tmp['EXPTIME'] = 360
            head_tmp['GAIN'] = 1
            head_tmp['SATLEVEL'] = 60000

        ## save to head if required
        if head_dir is not None:
            if not os.path.exists(head_dir):
                os.mkdir(head_dir)
            outpath_tmp = os.path.join(head_dir, f'expo{id_exposure}_chip{i_chip}.head')
            f = open(outpath_tmp, 'w')
            for line in head_tmp.cards:
                print(line, file=f)
            f.close()

def LensfitShape_psf(lensfit_dir, psf_dir):
    """
    prepare psf polynomial coefficients from psf image
    """

    psf_coeff_dir = os.path.join(psf_dir, os.path.basename(psf_dir) + '_coeff')
    if not os.path.exists(psf_coeff_dir):
        os.mkdir(psf_coeff_dir)

    psfimage2coeffs_path = lensfit_dir + '/utils/psfimage2coeffs'
    psf_imas = glob.glob(os.path.join(psf_dir, '*.fits'))

    for psf_ima in psf_imas:

        psf_coeff = os.path.join(psf_coeff_dir, os.path.basename(psf_ima).replace('.fits', '.psfcoeffs.fits'))
        cmd = [psfimage2coeffs_path, psf_ima, psf_coeff]
        # run
        proc = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)

    logger.info(f'PSF polynomial coefficients saved to {psf_coeff_dir}')
    return psf_coeff_dir

def LensfitShape(lensfit_dir,
                    input_catalog, input_file, chip_dir, psf_coeff_dir, head_dir, weight_dir=None,
                    PSF_OVERSAMPLING='1', PECUT='0.02', PRCUT='0.02', LCUT='0.05', WAVEBAND='R', CAMERA='KIDS',
                    postage_size='48', start_exposure='1', end_exposure='5', start_mag='20.0', end_mag='25.0',
                    output_file='output.feather', out_dir='./', tmp_dir='./tmp',
                    running_log=True, log_dir=None):
    """
        Lensfit wrapper
    """

    # check existence
    output_feather = os.path.join(out_dir, output_file)
    if os.path.isfile(output_feather):
        logger.info(f'The final feather catalogue {output_feather} already exists.')
        logger.info(f'End the process.')
        return 1

    # tmp directory
    if os.path.exists(tmp_dir):
        shutil.rmtree(tmp_dir)
    os.mkdir(tmp_dir)
    output_path = os.path.join(tmp_dir, 'output.fits')

    # running info
    if running_log:
        outLog = open(os.path.join(log_dir, output_file.replace('.feather', '.log')), "w")
        errLog = open(os.path.join(log_dir, output_file.replace('.feather', '.err.log')), "w")
    else:
        outLog = subprocess.PIPE
        errLog = subprocess.STDOUT

    # ++++++++++++ 1. run lensfit to obtain shape info
    lensfit_path = lensfit_dir + '/src/flensfit'

    # Environment variables
    envVars = os.environ.copy()
    envVars['CATALOGUE_GALAXIES'] = input_catalog
    envVars['DATA_DIR'] = chip_dir
    if weight_dir is not None:
        envVars['BADPIX_DIR'] = weight_dir
    envVars['PSF_DIR'] = psf_coeff_dir
    envVars['HEAD_DIR'] = head_dir
    envVars['SWARP_CONFIG'] = lensfit_dir + '/input_files/create_coadd_swarp.swarp'
    envVars['PRIOR_PARAMETERS'] = lensfit_dir + '/input_files/CFHTLenS_prior_parameters'
    envVars['PSF_OVERSAMPLING'] = PSF_OVERSAMPLING
    envVars['PECUT'] = PECUT
    envVars['PRCUT'] = PRCUT
    envVars['LCUT'] = LCUT
    envVars['WAVEBAND'] = WAVEBAND
    envVars['CAMERA'] = CAMERA

    # build command
    cmd = [lensfit_path, input_file, output_path,
            postage_size, start_exposure, end_exposure, start_mag, end_mag]

    # run
    print("########## flensfit info start ##########\n", file=outLog)
    print("########## flensfit error start ##########\n", file=errLog)
    proc = subprocess.run(cmd, stdout=outLog, stderr=errLog, env=envVars, cwd=tmp_dir)
    print("########## flensfit info end ##########\n", file=outLog)
    print("########## flensfit error end ##########\n", file=errLog)

    logger.info(f'Lensfit shape info saved as {output_path}')

    # running info
    if running_log:
        outLog.close()
        errLog.close()

    # ++++++++++++ 2. build convenient feather file
    data = np.loadtxt(f'{output_path}.asc')
    if data.size == 0:
        raise Exception(f'Empty {output_path}.asc, something is wrong about lensfit!')

    data_out = pd.DataFrame({
        'id_detec': data[:, 28].astype(int),
        'e1_LF': data[:, 2].astype(float),
        'e2_LF': data[:, 3].astype(float),
        'SNR_LF': data[:, 10].astype(float),
        'scalelength_LF': data[:, 6].astype(float), # bias-corrected scalelength /pixels
        'oldweight_LF': data[:, 4].astype(float),
        'psf_e1_LF': data[:, 12].astype(float),
        'psf_e2_LF': data[:, 13].astype(float),
        'psf_strehl_LF': data[:, 14].astype(float), # <PSF>-Strehl-ratio
        'psf_Q11_LF': data[:, 15].astype(float),
        'psf_Q22_LF': data[:, 16].astype(float),
        'psf_Q12_LF': data[:, 17].astype(float),
        'class_LF': data[:, 5].astype(float), # fitclass
        'e1_corr_LF': data[:, 22].astype(float), # e1 correction
        'e2_corr_LF': data[:, 23].astype(float), # e2 correction
        'scalelength_corr_LF': data[:, 19].astype(float), # size correction
        'contamination_radius_LF': data[:, 11].astype(float),
        'nm_LF': data[:, 24].astype(float), # neighbour mag
        'nd_LF': data[:, 25].astype(float), # neighbour distance
        'LS_variance_LF': data[:, 20].astype(float), # 2D measurement variance
        'bulge_fraction_LF': data[:, 7].astype(float),
        'SG_prob_LF': data[:, 18].astype(float), # star-galaxy f-probability
        })

    data_out.to_feather(output_feather)

    logger.info(f'Final catalogue saved as {output_feather}')

    return 0
