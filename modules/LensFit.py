# -*- coding: utf-8 -*-
# @Author: lshuns
# @Date:   2020-12-03 16:16:21
# @Last Modified by:   lshuns
# @Last Modified time: 2021-10-07 16:22:06

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

class LensFITwrapper(object):
    """
    A wrapper for running lensfit code.
        To obtain lensfit shape measurements.
        Automatically manages a temporary working directory, environment variables.
    """

    def __init__(self, lensfit_dir, out_dir, tmp_dir,
                    PSF_OVERSAMPLING='1', PECUT='0.02', PRCUT='0.02', LCUT='0.05', CAMERA='KIDS',
                    postage_size='48', start_exposure='1', end_exposure='5', start_mag='20.0', end_mag='25.0',
                    lensfit_cores=12,
                    running_log=False, log_dir=None):

        logger.info("Initialising the lensfit wrapper...")

        # global parameters
        self._lensfitDir = lensfit_dir
        self._postage_size = postage_size
        self._start_exposure = start_exposure
        self._end_exposure = end_exposure
        self._start_mag = start_mag
        self._end_mag = end_mag

        # directories
        self._outDir = out_dir
        ## a temp directory for intermediate files
        self._tmp_dir = tmp_dir
        if os.path.exists(tmp_dir):
            shutil.rmtree(tmp_dir)
        os.mkdir(tmp_dir)
        ## log info
        self._running_log = running_log
        if running_log:
            self._logDir = log_dir

        # general environment variables
        self._envVars = os.environ.copy()
        self._envVars['SWARP_CONFIG'] = lensfit_dir + '/input_files/create_coadd_swarp.swarp'
        self._envVars['PRIOR_PARAMETERS'] = lensfit_dir + '/input_files/CFHTLenS_prior_parameters'
        self._envVars['PSF_OVERSAMPLING'] = PSF_OVERSAMPLING
        self._envVars['PECUT'] = PECUT
        self._envVars['PRCUT'] = PRCUT
        self._envVars['LCUT'] = LCUT
        self._envVars['CAMERA'] = CAMERA

        # prepare input file
        self._input_file = os.path.join(tmp_dir, f'lensfit_input.asc')
        f = open(self._input_file, 'w')
        if CAMERA == 'KIDS':
            N_expo = 5
            N_chips = 32
        else:
            raise Exception(f'Unsupported CAMERA {CAMERA} for lensfit!')
        for i_expo in range(N_expo):
            for i_chip in range(N_chips):
                print(f'expo{i_expo}_chip{i_chip}', file=f)
        f.close()
        logger.debug(f'input info saved to {self._input_file}')

        # lensfit exe
        self._lensfit = lensfit_dir + f'/bin/flensfit_NT{lensfit_cores}'
        if not os.path.isfile(self._lensfit):
            logger.warning(f'the required {self._lensfit} does not exist!')
            self._lensfit = lensfit_dir + f'/bin/flensfit_NT12'
            logger.warning(f'Use {self._lensfit} instead!')

    def _LensfitShape_head(self, chip_dir, head_dir=None):
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
                if os.path.exists(head_dir):
                    shutil.rmtree(head_dir)
                os.mkdir(head_dir)
                outpath_tmp = os.path.join(head_dir, f'expo{id_exposure}_chip{i_chip}.head')
                f = open(outpath_tmp, 'w')
                for line in head_tmp.cards:
                    print(line, file=f)
                f.close()

    def LensfitShape_psf(self, psf_dir, varChips=False):
        """
        prepare psf polynomial coefficients from psf image
        """

        psf_coeff_dir = os.path.join(psf_dir, os.path.basename(psf_dir) + '_coeff')
        if os.path.exists(psf_coeff_dir):
            shutil.rmtree(psf_coeff_dir)
        os.mkdir(psf_coeff_dir)

        # psf vary between different chips
        if varChips:
            logger.info('PSF varies between chips')
            psfimage2coeffs_path = self._lensfitDir + '/utils/psf_meimage2coeffs'
        # psf constant across the field
        else:
            logger.info('PSF is constant across the field')
            psfimage2coeffs_path = self._lensfitDir + '/utils/psfimage2coeffs'


        # apply to all the PSF images
        psf_imas = glob.glob(os.path.join(psf_dir, '*.fits'))
        for psf_ima in psf_imas:
            psf_coeff = os.path.join(psf_coeff_dir, os.path.basename(psf_ima).replace('.fits', '.psfcoeffs.fits'))
            cmd = [psfimage2coeffs_path, psf_ima, psf_coeff]
            # run
            proc = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)

        logger.info(f'PSF polynomial coefficients saved to {psf_coeff_dir}')
        return psf_coeff_dir

    def LensfitShape(self, output_file,
                        path_input_cata, chip_dir, psf_coeff_dir, weight_dir=None,
                        WAVEBAND='R'):
        """
            lensfit running
        """

        logger.info(f'Running lensfit for {output_file}...')

        # check existence
        output_feather = os.path.join(self._outDir, output_file)
        if os.path.isfile(output_feather):
            logger.info(f'The final feather catalogue {output_feather} already exists.')
            logger.info(f'End the process.')
            return 1

        # running info
        if self._running_log:
            outLog = open(os.path.join(self._logDir, output_file.replace('.feather', '.log')), "w")
            errLog = open(os.path.join(self._logDir, output_file.replace('.feather', '.err.log')), "w")
        else:
            outLog = subprocess.PIPE
            errLog = subprocess.STDOUT

        # tmp files
        tmp_dir = os.path.join(self._tmp_dir, os.path.splitext(output_file)[0])
        if os.path.exists(tmp_dir):
            shutil.rmtree(tmp_dir)
        os.mkdir(tmp_dir)

        # >>>>>>>>>>>>> 1. prepare detection catalogue
        input_cata = pd.read_feather(path_input_cata)
        #### desired info
        input_cata = np.array([input_cata['X_WORLD'], input_cata['Y_WORLD'], input_cata['MAG_AUTO'], input_cata['NUMBER']]).T
        #### save
        input_cata_path = os.path.join(tmp_dir, 'catalog.asc')
        np.savetxt(input_cata_path, input_cata, fmt=['%.8f', '%.8f', '%.3f', '%i'])
        logger.debug(f'catalogue info saved to {input_cata_path}')

        # >>>>>>>>>>>>> 2. prepare headers for lensfit 
        head_dir = os.path.join(tmp_dir, os.path.basename(chip_dir)+'_head')
        self._LensfitShape_head(chip_dir, head_dir)

        # >>>>>>>>>>>>> 3. running lensfit
        ## output
        output_path = os.path.join(tmp_dir, 'output.fits')
        #### check existence
        if os.path.isfile(output_path):
            os.remove(output_path)
        if os.path.isfile(f'{output_path}.asc'):
            os.remove(f'{output_path}.asc')
        ## Environment variables
        envVars = self._envVars.copy()
        envVars['CATALOGUE_GALAXIES'] = input_cata_path
        envVars['DATA_DIR'] = chip_dir
        if weight_dir is not None:
            envVars['BADPIX_DIR'] = weight_dir
        envVars['PSF_DIR'] = psf_coeff_dir
        envVars['HEAD_DIR'] = head_dir
        envVars['WAVEBAND'] = WAVEBAND
        ## build command
        cmd = [self._lensfit, self._input_file, output_path,
                self._postage_size, self._start_exposure, self._end_exposure, self._start_mag, self._end_mag]
        ## run
        proc = subprocess.run(cmd, stdout=outLog, stderr=errLog, env=envVars, cwd=tmp_dir)
        logger.debug(f'Lensfit original shape info saved as {output_path}')
        ## running info
        if self._running_log:
            outLog.close()
            errLog.close()

        # >>>>>>>>>>>>> 4. build convenient feather file
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
            'N_expo_used_LF': data[:, 27].astype(int) # number of exposures used by lensfit
            })
        data_out.to_feather(output_feather)

        logger.info(f'Final catalogue saved as {output_feather}')

        return 0
