# -*- coding: utf-8 -*-
# @Author: lshuns
# @Date:   2021-01-22 13:37:08
# @Last modified by:   lshuns
# @Last modified time: 2021-02-18, 22:52:09

### Wrapper for BPZ code

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

class BPZwrapper(object):
    """
    A wrapper for running BPZ code.
        To obtain Bayesion photometric reshifts.
        Automatically manages a temporary working directory, environment variables.
    """

    def __init__(self, python2_env, BPZ_dir,
                    out_dir, tmp_dir,
                    band_list, band_FilterName_list,
                    band_CataNameBase, banderr_CataNameBase, 
                    bandflag_CataNameBase, bandlim_CataNameBase, 
                    detec_band='r',
                    photo_sys='AB', 
                    prior_band='i', prior_name='NGVS',
                    templates_name='CWWSB_capak', interpolation='10',
                    lkl_zmin='0.01', lkl_zmax='2.5', lkl_dz='0.01', lkl_odds='0.68', lkl_min_rms='0.067',
                    running_log=True, log_dir=None):

        logger.info("Initialising the BPZ wrapper...")

        # global parameters
        self._exe = python2_env
        self._BPZ_code = os.path.join(BPZ_dir, 'bpz.py')

        self._outDir = out_dir
        # a temp directory for intermediate files
        self._tmp_dir = tmp_dir
        if not os.path.exists(self._tmp_dir):
            os.mkdir(self._tmp_dir)

        self._bands = band_list
        logger.info(f"Used bands: {band_list}")
        logger.info(f"Used filters: {band_FilterName_list}")
        self._prior_band = prior_band
        logger.info(f"prior band: {prior_band}")
        self._band_CataNameBase = band_CataNameBase
        self._banderr_CataNameBase = banderr_CataNameBase
        self._bandflag_CataNameBase = bandflag_CataNameBase
        self._bandlim_CataNameBase = bandlim_CataNameBase

        self._detec_band = detec_band
        logger.info(f"detection band: {detec_band}")        

        self._prior_name = prior_name
        self._templates_name = templates_name
        logger.info(f"prior: {prior_name}; templates: {templates_name}")

        self._interpolation = interpolation

        self._lkl_zmin = lkl_zmin
        self._lkl_zmax = lkl_zmax
        self._lkl_dz = lkl_dz
        self._lkl_odds = lkl_odds
        self._lkl_min_rms = lkl_min_rms

        self._running_log = running_log
        if running_log:
            self._logDir = log_dir

        # Environment variables
        self._envVars = os.environ.copy()
        self._envVars['BPZPATH'] = BPZ_dir
        self._envVars["NUMERIX"] = "numpy"

        # prepare input columns
        self._columns_file = os.path.join(self._tmp_dir, 'input.columns')
        header = "# filter  columns  AB/Vega  zp_error  zp_offset\n"
        lines = header
        ## fixed photometric zero-points
        zp_error, zp_offset = 0.01, 0.0
        ## bands info
        for i_band, band_FilterName in enumerate(band_FilterName_list):
            mag_col_idx = int(2 * i_band + 3)
            err_col_idx = int(2 * i_band + 4)
            lines += f"{band_FilterName}  {mag_col_idx},{err_col_idx}  {photo_sys}  {zp_error}  {zp_offset}\n"
        ## the prior magnitude
        lines += f"M_0      2\n"
        lines += f"ID       1\n"
        ## write
        with open(self._columns_file, "w") as f:
            f.write(lines)

    def RunSingleTile(self, PhotoFile, DetecFile):
        """
        run single process for one tile
        """

        infile_basename = os.path.basename(PhotoFile)
        logger.info(f'Running BPZ for {infile_basename}...')

        # check exist
        output_feather = os.path.join(self._outDir, infile_basename)
        if os.path.isfile(output_feather):
            logger.info(f'Final feather catalog {output_feather} already exist.')
            logger.info(f'      end the process.')
            return 1

        # load file
        cata = pd.read_feather(PhotoFile)
        ## we need MAG_AUTO from detection cata
        cata_tmp = pd.read_feather(DetecFile)
        ### check if they are in the same order
        if np.sum(cata['id_detec'].values == cata_tmp['NUMBER'].values) != len(cata):
            raise Exception('Provided photometric and detection catalogues have different IDs!')
        ### get mag_auto
        MAG_AUTO = cata_tmp['MAG_AUTO'].values
        del cata_tmp

        # +++++++ prepare input file
        input_file = os.path.join(self._tmp_dir, infile_basename.replace('.feather', '.input'))
        if os.path.isfile(input_file):
            os.remove(input_file)
            logger.info(f'Removed existing input file.')
        N_bands = len(self._bands)

        ## pre-process of the magnitude
        for band in self._bands:        

            ### 1. objects being flagged are counted as 
            ###     not observed in the current filter
            mask_false = (cata[f'{self._bandflag_CataNameBase}_{band}'].values > 0) \
                            | (cata[f'{self._band_CataNameBase}_{band}'].values < 0) 
            cata.loc[mask_false, f'{self._band_CataNameBase}_{band}'] = -99.
            cata.loc[mask_false, f'{self._banderr_CataNameBase}_{band}'] = 0.
            del mask_false

            ### 2. if an object is not flagged but has a magnitude
            ###     fainter than the limiting magnitude (or a large error)
            ###     we count it as not detected in the current filter
            mask_false = (cata[f'{self._band_CataNameBase}_{band}'].values > 0) \
                            & ((cata[f'{self._band_CataNameBase}_{band}'].values \
                                    > cata[f'{self._bandlim_CataNameBase}_{band}'].values) \
                                | (cata[f'{self._banderr_CataNameBase}_{band}'].values > 1))
            cata.loc[mask_false, f'{self._band_CataNameBase}_{band}'] = 99.
            cata.loc[mask_false, f'{self._banderr_CataNameBase}_{band}'] \
                = cata.loc[mask_false, f'{self._bandlim_CataNameBase}_{band}'].values

        ## define the GAaP-to-total aperture correction
        #### from the detection band
        apcorr = cata[f'{self._band_CataNameBase}_{self._detec_band}'].values - MAG_AUTO

        ## correct the prior band magnitude
        #### use the prior band if it is possible
        #### otherwise using MAG_AUTO
        mag_prior = cata[f'{self._band_CataNameBase}_{self._prior_band}'].values
        mag_ref = np.where((mag_prior > -99) & (mag_prior < 99), 
                            mag_prior - apcorr,
                            MAG_AUTO)
        del apcorr, MAG_AUTO, mag_prior

        ## collect all photometry for bpz
        data_out = np.empty([len(cata), int(N_bands*2) + 2])
        ### first column is for ID
        data_out[:, 0] = cata['id_detec'].values
        ### second column is for reference band
        data_out[:, 1] = mag_ref
        del mag_ref
        ### other columns for used bands
        for i_band, band in enumerate(self._bands):
            data_out[:, int(2*i_band + 2)] = cata[f'{self._band_CataNameBase}_{band}'].values
            data_out[:, int(2*i_band + 3)] = cata[f'{self._banderr_CataNameBase}_{band}'].values
        ## save
        del cata
        np.savetxt(input_file, data_out, fmt=' '.join(['%i'] + ['%.4f']*(data_out.shape[1]-1)))
        del data_out
        logger.info(f'BPZ input file created as {input_file}')

        # +++++++ BPZ running
        output_file = os.path.join(self._tmp_dir, infile_basename.replace('.feather', '.output'))
        ## remove potential intermediate file
        if os.path.isfile(output_file):
            os.remove(output_file)
            logger.info(f'Removed existing output file.')

        # running info
        if self._running_log:
            outLog = open(os.path.join(self._logDir, infile_basename.replace('.feather', '.log')), "w")
            errLog = open(os.path.join(self._logDir, infile_basename.replace('.feather', '.err.log')), "w")
        else:
            outLog = subprocess.PIPE
            errLog = subprocess.STDOUT

        # build command
        cmd = [self._exe, self._BPZ_code,
                input_file, "-COLUMNS", self._columns_file, "-OUTPUT", output_file,
                "-PRIOR", self._prior_name,
                # templates
                "-SPECTRA", f"{self._templates_name}.list",
                "-INTERP", str(self._interpolation),
                # likelihood
                "-ZMIN", str(self._lkl_zmin),
                "-ZMAX", str(self._lkl_zmax),
                "-DZ", str(self._lkl_dz),
                "-ODDS", str(self._lkl_odds),
                "-MIN_RMS", str(self._lkl_min_rms),
                # Define the confidence interval using only the photometric errors
                "-PHOTO_ERRORS", "yes",
                # disable extra stuff we don't need
                "-PROBS", "no", "-PROBS2", "no", "-PROBS_LITE", "no",
                "-NEW_AB", "no", "-CHECK", "no", "-INTERACTIVE", "no",
                "-VERBOSE", "yes" if self._running_log else "no"]

        # run
        proc = subprocess.run(cmd, stdout=outLog, stderr=errLog, cwd=self._tmp_dir, env=self._envVars)

        # running info
        if self._running_log:
            outLog.close()
            errLog.close()

        # collect outputs to a feather file
        data = np.loadtxt(output_file)
        data = pd.DataFrame({
            'id_detec': data[:, 0].astype(int),
            'Z_B': data[:, 1].astype(float),
            'Z_B_MIN': data[:, 2].astype(float),
            'Z_B_MAX': data[:, 3].astype(float),
            'T_B': data[:, 4].astype(float),
            'ODDS': data[:, 5].astype(float),
            'Z_ML': data[:, 6].astype(float),
            'T_ML': data[:, 7].astype(float),
            'CHI-SQUARED': data[:, 8].astype(float),
            'M_0': data[:, 9].astype(float)
            })

        # save
        tmp_output_feather = output_feather + '_tmp'
        data.to_feather(tmp_output_feather)
        os.rename(tmp_output_feather, output_feather)
        del data
        logger.info(f'BPZ outputs saved as {output_feather}.')
