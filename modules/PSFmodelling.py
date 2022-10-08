# -*- coding: utf-8 -*-
# @Author: lshuns
# @Date:   2021-11-11 13:29:22
# @Last Modified by:   lshuns
# @Last Modified time: 2022-08-26 15:33:43

### Everything about PSF modelling
__all__ = ['ima2coeffsFunc', 'makeglobalpsfFunc']

import os
import re
import sys
import glob
import shutil
import logging
import subprocess

import numpy as np
import pandas as pd

from pathlib import Path
from astropy.io import fits

logger = logging.getLogger(__name__)

def ima2coeffsFunc(ima2coeffs_dir, in_dir, out_dir, varChips=False):
    """
    get psf polynomial coefficients from psf image
    """

    logger.info(f'Running ima2coeffs for {in_dir}...')
    psf_imas = glob.glob(os.path.join(in_dir, '*.fits'))
    if len(psf_imas) < 1:
        raise Exception(f'No PSF images found in {in_dir}')
    logger.info(f'Number of PSF images: {len(psf_imas)}')

    # psf vary between different chips
    if varChips:
        logger.info('PSF varies between chips')
        psfimage2coeffs_path = ima2coeffs_dir + '/psf_meimage2coeffs'
    # psf constant across the field
    else:
        logger.info('PSF is constant across the field')
        psfimage2coeffs_path = ima2coeffs_dir + '/psfimage2coeffs'

    # apply to all the PSF images
    for psf_ima in psf_imas:

        # get the id
        id_exposure = re.findall(r"\d+", os.path.basename(psf_ima))[0]
        # get coeff name for lensfit
        psf_coeff = os.path.join(out_dir, 'exp'+id_exposure+'chip.psfcoeffs.fits')

        # run
        cmd = [psfimage2coeffs_path, psf_ima, psf_coeff]
        proc = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)

    logger.info(f'PSF polynomial coefficients saved to {out_dir}')

def makeglobalpsfFunc(makeglobalpsf_dir, path_input_cata, 
                    in_dir, out_dir, 
                    SWARP_CONFIG, CAMERA='KIDS', WAVEBAND='R',
                    global_order='4', chip_order='1', snratio='20', start_mag='18.0', end_mag='24.0', 
                    weight_dir=None,
                    running_log=False, logDir=None,
                    clean_intermediate=True):
    """
    get psf polynomial coefficients from stars
    """

    logger.info(f'Running makeglobalpsf for {in_dir}...')

    # general info
    if CAMERA == 'KIDS':
        N_expo = 5
        N_chips = 32
    else:
        raise Exception(f'Unsupported CAMERA {CAMERA} for makeglobalpsf!')

    # check existence
    tmp_list = glob.glob(os.path.join(out_dir, 'exp*.psfcoeffs.fits'))
    if len(tmp_list) == N_expo:
        Nsuccess = 0
        # check globalshifts
        for tmp in tmp_list:
            with fits.open(tmp) as hdul:
                head_tmp = hdul[0].header
            if 'SORDER' in head_tmp:
                Nsuccess += 1
            del head_tmp, hdul
        if Nsuccess == N_expo:
            logger.info(f'PSF coeff for {in_dir} already exists.')
            logger.info(f'End the process.')
            return 1
    del tmp_list

    # clean potential intermediate outcomes
    if os.path.exists(out_dir):
        shutil.rmtree(out_dir)
    os.mkdir(out_dir)
    if os.path.exists(out_dir+'_tmp'):
        shutil.rmtree(out_dir+'_tmp')

    # running info
    if running_log:
        outLog = open(os.path.join(logDir, os.path.basename(path_input_cata).replace('.feather', f'_order{global_order}{chip_order}.log')), "w")
        errLog = open(os.path.join(logDir, os.path.basename(path_input_cata).replace('.feather', f'_order{global_order}{chip_order}.err.log')), "w")
        outLog_shifts = open(os.path.join(logDir, os.path.basename(path_input_cata).replace('.feather', f'_order{global_order}{chip_order}_shifts.log')), "w")
        errLog_shifts = open(os.path.join(logDir, os.path.basename(path_input_cata).replace('.feather', f'_order{global_order}{chip_order}_shifts.err.log')), "w")
    else:
        outLog = subprocess.PIPE
        errLog = subprocess.STDOUT
        outLog_shifts = subprocess.PIPE
        errLog_shifts = subprocess.STDOUT

    # >>>>>>>>>>>> 1. prepare input files for makeglobalpsf    
    for i_expo in range(N_expo):
        input_file = os.path.join(out_dir, f'exp{i_expo}chip.asc')
        f = open(input_file, 'w')
        for i_chip in range(N_chips):
            print(f'exp{i_expo}chip_{i_chip+1}OFCS', file=f)
        f.close()
        logger.debug(f'input info saved to {input_file}')

    # >>>>>>>>>>>>> 2. prepare detection catalogue
    input_cata = pd.read_feather(path_input_cata)
    #### get stars
    try:
        input_cata = input_cata[input_cata['perfect_flag_star']==1]
        logger.debug(f'Total number of stars found {len(input_cata)}')
    except KeyError:
        raise Exception(f'no stars found in {path_input_cata}!')
    #### desired info
    input_cata = np.array([input_cata['X_WORLD'], input_cata['Y_WORLD'], input_cata['MAG_AUTO'], input_cata['NUMBER']]).T
    #### save
    star_cata_path = os.path.join(out_dir, 'stars.asc')
    np.savetxt(star_cata_path, input_cata, fmt=['%.8f', '%.8f', '%.3f', '%i'])
    logger.debug(f'star catalogue info saved to {star_cata_path}')

    # >>>>>>>>>>>>> 3. prepare headers for lensfit 
    head_dir = os.path.join(in_dir, 'head_info')
    if os.path.isdir(head_dir):
        shutil.rmtree(head_dir)
    os.mkdir(head_dir)
    _LensfitShape_head(in_dir, head_dir)

    # >>>>>>>>>>>>> 4. prepare environment variables
    envVars = os.environ.copy()
    envVars['SWARP_CONFIG'] = SWARP_CONFIG
    envVars['CAMERA'] = CAMERA
    envVars['WAVEBAND'] = WAVEBAND
    envVars['PSF_DIR'] = out_dir
    envVars['DATA_DIR'] = in_dir
    envVars['HEAD_DIR'] = head_dir
    envVars['CATALOGUE_STARS'] = star_cata_path
    if weight_dir is not None:
        envVars['BADPIX_DIR'] = weight_dir

    # >>>>>>>>>>>>> 5. running makeglobalpsf
    makeglobalpsf = os.path.join(makeglobalpsf_dir, 'makeglobalpsf')
    for i_expo in range(N_expo):
        input_file = f'exp{i_expo}chip.asc'

        ## build command
        cmd = [makeglobalpsf, input_file, global_order, chip_order, snratio, start_mag, end_mag]
        ## run
        proc = subprocess.run(cmd, stdout=outLog, stderr=errLog, env=envVars, cwd=out_dir)

    ## check outputs
    tmp = glob.glob(os.path.join(out_dir, 'exp*.psfcoeffs.fits'))
    if len(tmp) != N_expo:
        raise Exception('makeglobalpsf did not success, pleace check running_log!')
    logger.debug(f'makeglobalpsf finished')

    # >>>>>>>>>>>>> 6. running globalshifts
    globalshifts = os.path.join(makeglobalpsf_dir, 'globalshifts')
    for i_expo in range(N_expo):
        input_file = f'exp{i_expo}chip.asc' 

        ## build command
        cmd = [globalshifts, input_file, snratio, start_mag, end_mag]
        ## run
        proc = subprocess.run(cmd, stdout=outLog_shifts, stderr=errLog_shifts, env=envVars, cwd=out_dir)

    ## check outputs
    for tmp in glob.glob(os.path.join(out_dir, 'exp*.psfcoeffs.fits')):
        with fits.open(tmp) as hdul:
            head_tmp = hdul[0].header
        if not 'SORDER' in head_tmp:
            raise Exception('globalshifts did not success, pleace check running_log!')
    logger.debug(f'globalshifts finished')

    # >>>>>>>>>>>>> 7. clean intermediate files
    if clean_intermediate:

        # save useful outputs
        os.rename(out_dir, out_dir+'_tmp')
        os.mkdir(out_dir)
        for tmp_file in glob.glob(os.path.join(out_dir+'_tmp', 'exp*.psfcoeffs.fits')):
            shutil.move(tmp_file, out_dir)

        # clean others
        shutil.rmtree(out_dir+'_tmp')
        logger.info('intermediate files are cleaned')

    ## running info
    if running_log:
        outLog.close()
        errLog.close()
        outLog_shifts.close()
        errLog_shifts.close()

    logger.info(f'PSF polynomial coefficients saved to {out_dir}')


def _LensfitShape_head(chip_dir, head_dir=None):
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
            head_tmp['IMAGEID'] = int(i_chip)
            head_tmp['EXPTIME'] = 360
            head_tmp['GAIN'] = 1
            head_tmp['SATLEVEL'] = 60000

        ## save to head if required
        if head_dir is not None:
            outpath_tmp = os.path.join(head_dir, f'exp{id_exposure}chip_{i_chip}.head')
            f = open(outpath_tmp, 'w')
            for line in head_tmp.cards:
                print(line, file=f)
            f.close()