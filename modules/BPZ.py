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
                    band_list, band_CataName_list, banderr_CataName_list, band_FilterName_list,
                    photo_sys='AB', prior_band='i', prior_name='NGVS',
                    templates_name='CWWSB_capak', interpolation='10',
                    lkl_zmin='0.0002', lkl_zmax='2.2984', lkl_dz='0.01', lkl_odds='0.68', lkl_min_rms='0.0',
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
        self._bands_CataName = band_CataName_list
        self._banderrs_CataName = banderr_CataName_list

        self._prior_name = prior_name

        self._templates_name = templates_name
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
            mag_col_idx = int(2 * i_band + 1)
            err_col_idx = int(2 * i_band + 2)
            lines += f"{band_FilterName}  {mag_col_idx},{err_col_idx}  {photo_sys}  {zp_error}  {zp_offset}\n"
        ## the prior magnitude
        M_0_idx = band_list.index(prior_band)
        lines += f"M_0  {int(2 * M_0_idx + 1)}\n"
        ## write
        with open(self._columns_file, "w") as f:
            f.write(lines)

    def RunSingleTile(self, cata_pathfile):
        """
        run single process for one tile
        """

        infile_basename = os.path.basename(cata_pathfile)
        logger.info(f'Running BPZ for {infile_basename}...')

        # check exist
        output_feather = os.path.join(self._outDir, infile_basename)
        if os.path.isfile(output_feather):
            logger.info(f'Final feather catalog {output_feather} already exist.')
            logger.info(f'      end the process.')
            return 1

        # load file
        cata = pd.read_feather(cata_pathfile)

        # +++++++ prepare input file
        input_file = os.path.join(self._tmp_dir, infile_basename.replace('.feather', '.input'))
        if os.path.isfile(input_file):
            os.remove(input_file)
            logger.info(f'Removed existing input file.')
        N_bands = len(self._bands)
        data_out = np.empty([len(cata), int(N_bands*2)])
        for i_band in range(N_bands):
            data_out[:, int(2*i_band)] = np.array(cata[self._bands_CataName[i_band]])
            data_out[:, int(2*i_band+1)] = np.array(cata[self._banderrs_CataName[i_band]])
        ## save
        np.savetxt(input_file, data_out)
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
        data_out = pd.DataFrame({'id_detec': np.array(cata['id_detec']).astype(int)})
        data_out.loc[:, 'Z_B'] = data[:, 1]
        data_out.loc[:, 'Z_B_MIN'] = data[:, 2]
        data_out.loc[:, 'Z_B_MAX'] = data[:, 3]
        data_out.loc[:, 'T_B'] = data[:, 4]
        data_out.loc[:, 'ODDS'] = data[:, 5]
        data_out.loc[:, 'Z_ML'] = data[:, 6]
        data_out.loc[:, 'T_ML'] = data[:, 7]
        data_out.loc[:, 'CHI-SQUARED'] = data[:, 8]
        data_out.loc[:, 'M_0'] = data[:, 9]

        # save
        tmp_output_feather = output_feather + '_tmp'
        data_out.to_feather(tmp_output_feather)
        os.rename(tmp_output_feather, output_feather)
        logger.info(f'BPZ outputs saved as {output_feather}.')
