# -*- coding: utf-8 -*-
# @Author: lshuns
# @Date:   2021-01-13 19:27:48
# @Last modified by:   ssli
# @Last modified time: 2021-05-02, 12:05:20

### Wrapper for GAaP code

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

from CrossMatch import run_position2id

logger = logging.getLogger(__name__)

class GAaPwrapper(object):
    """
    A wrapper for running GAaP pipeline.
        To obtain GAaP magnitude for detected sources.
        Automatically manages a temporary working directory, environment variables required by GAaP.
        Work with/without PSF map.

    Parameters:
    -----------
    gaap_dir : str
        GAaP pipeline directory.
    out_dir : str
        Output directory.
    star_mag_cut : list of float, optional (default: [16, 20])
        Within which magnitude range stars are used.
    mag_zero : float, optional (default: 30.)
        Zero point for the magnitude.
    min_aper : float, optional (default: 0.7)
        Minimum aperture size.
    max_aper : float, optional (default: 2.0)
        Maximum aperture size.
    mag_1sigma_limits : dictionary
        1 sigma limiting magnitude, used for unmeasured objects' error
    running_log : bool, optional (default: True)
        save running info from external code
    clean_up_level : int, optional (default: 0)
        Clean up level
        0: none
        1: tmp directory
        2: and *.gaap files
    """

    def __init__(self, gaap_dir, out_dir,
                star_mag_cut=[16., 20.],
                mag_zero=30.,
                min_aper=0.7, max_aper=2.0,
                mag_1sigma_limits={},
                running_log=True, log_dir=None,
                clean_up_level=0):

        logger.info("Initialising the GAaP wrapper...")

        self._GAAP = gaap_dir
        logger.info("GAaP directory: %s", self._GAAP)
        self._outDir = out_dir
        logger.info("Output directory: %s", self._outDir)
        if not os.path.exists(self._outDir):
            os.mkdir(self._outDir)

        self._star_mag_cut = star_mag_cut
        logger.info(f'star magnitude cut: ({star_mag_cut[0]}, {star_mag_cut[1]})')

        self._magzero = mag_zero
        logger.info("magnitude zero point: %f", self._magzero)

        self._minaper = min_aper
        logger.info("Minimum aperture: %f", self._minaper)
        self._maxaper = max_aper
        logger.info("Maximum aperture: %f", self._maxaper)

        self._1sigma_limits = mag_1sigma_limits
        logger.info(f'1-sigma limiting magnitude: {mag_1sigma_limits}')

        self._running_log = running_log
        logger.info(f"Save running info: {self._running_log}")
        self._logDir = log_dir

        self._clean_up_level = clean_up_level
        logger.info("Clean up level: {:}".format(self._clean_up_level))

    def _PSFcataFunc(self, band, SKYimaFile, star_info=None, outLog=subprocess.PIPE, errLog=subprocess.STDOUT):
        """
        make PSF catalogue in right SExtractor ascii format
        """

        # product
        PSFcataFile = SKYimaFile.replace('.fits', '_psf.cat')

        # check if PSFcata already exist
        if os.path.isfile(PSFcataFile) and (os.path.getsize(PSFcataFile)>0):

            logger.info("PSF catalogue already exist.")

        else:

            # using psf map
            if star_info is None:

                logger.info("Producing PSF catalogue from separated PSF map...")

                # necessary sex config files
                sex_config = []
                for file in ['gaap.sex', 'gaap.sexpar', 'default.conv']:
                    sex_config.append(os.path.join(self._GAAP, 'lib/{:}'.format(file)))

                # run sex
                cmd = ['sex', SKYimaFile, '-c', sex_config[0], '-PARAMETERS_NAME', sex_config[1], '-FILTER_NAME', sex_config[2], '-DETECT_THRESH', str(10), '-CATALOG_NAME', PSFcataFile]
                print("########## sex info start ##########\n", file=outLog)
                print("########## sex error start ##########\n", file=errLog)
                proc = subprocess.run(cmd, stdout=outLog, stderr=errLog)
                print("########## sex info end ##########\n", file=outLog)
                print("########## sex error end ##########\n", file=errLog)

            # using stars
            else:

                logger.info("Producing PSF catalogue using stars within image...")

                # intermediate
                DetecFile = SKYimaFile.replace('.fits', '_all.sex')

                # +++++++++ sextractor detection

                # necessary sex config files
                sex_config = []
                for file in ['gaap.sex', 'gaap_ima.sexpar', 'default.conv']:
                    sex_config.append(os.path.join(self._GAAP, 'lib/{:}'.format(file)))

                # run sex
                cmd = ['sex', SKYimaFile, '-c', sex_config[0], '-PARAMETERS_NAME', sex_config[1], '-FILTER_NAME', sex_config[2], '-DETECT_THRESH', str(10), '-CATALOG_NAME', DetecFile]
                print("########## sex info start ##########\n", file=outLog)
                print("########## sex error start ##########\n", file=errLog)
                proc = subprocess.run(cmd, stdout=outLog, stderr=errLog)
                print("########## sex info end ##########\n", file=outLog)
                print("########## sex error end ##########\n", file=errLog)

                # +++++++++ cross match to find stars

                ## star info
                mask_good = (star_info[f'{band}_input'] > self._star_mag_cut[0]) & (star_info[f'{band}_input'] < self._star_mag_cut[1])
                input_cata = star_info[mask_good].copy()
                input_cata.reset_index(drop=True, inplace=True)

                ## detected catalogue
                detec_cata_ori = np.loadtxt(DetecFile)
                detec_cata = pd.DataFrame({'NUMBER': detec_cata_ori[:, 8].astype(int),
                                            'X_WORLD': detec_cata_ori[:, 11],
                                            'Y_WORLD': detec_cata_ori[:, 12],
                                            'MAG_AUTO': self._magzero-2.5*np.log10(detec_cata_ori[:, 2])})

                ## Match
                id_list = ['index_input', 'NUMBER']
                position_list = [['RA_input', 'DEC_input'], ['X_WORLD', 'Y_WORLD']]
                mag_list = [f'{band}_input', 'MAG_AUTO']
                matched_cata, _, _ = run_position2id(input_cata, detec_cata, id_list, position_list, mag_list,
                                    outDir=None, basename=None, save_matched=False, save_false=False, save_missed=False,
                                    r_max=0.6/3600., k=2, mag_closest=True, running_info=False)

                ## stars
                id_stars = matched_cata['id_detec'].values
                mask = np.isin(detec_cata_ori[:, 8], id_stars)

                # +++++++++ save desired info
                psf_cata = detec_cata_ori[mask, :11]
                logger.info(f"Number of stars selected {len(psf_cata)}")
                np.savetxt(PSFcataFile, psf_cata,
                            fmt='%.4f %.4f %.2f %.5f %.3f %.3f %.3f %.2f %i %i %.2f')

            logger.info(f"PSF catalogue Extracted as {PSFcataFile}")

        return PSFcataFile

    def _GaussianizationFunc(self, SKYimaFile, PSFcataFile, tmp_dir, outLog=subprocess.PIPE, errLog=subprocess.STDOUT):
        """
        calculate Gaussianization kernel and Gaussianised image
        allow Gaussianised PSF radius to be determined autimatically
        kernel has shapelet order 8, no spatial variation for the PSF
        """

        logger.info("Producing PSF-Gaussianised image...")

        # products
        SKYimaKerFile = os.path.join(tmp_dir, os.path.basename(SKYimaFile).replace('.fits', '.ker.map'))
        SKYimaGpsfFile = os.path.join(tmp_dir, os.path.basename(SKYimaFile).replace('.fits', '.gpsf.fits'))

        # check if SKYimaGpsf already exist
        if os.path.isfile(SKYimaGpsfFile) and (os.path.getsize(SKYimaGpsfFile)>0):
            logger.info("PSF-Gaussianised image already exist.")

        else:

            ## orders
            np.savetxt(os.path.join(tmp_dir, 'orders.par'), [[8], [0]], fmt='%d')
            # np.savetxt(os.path.join(tmp_dir, 'orders.par'), [[8], [3]], fmt='%d')

            ## link to the image
            inimage = os.path.join(tmp_dir, 'inimage.fits')
            if os.path.isfile(inimage):
                os.remove(inimage)
            os.symlink(SKYimaFile, inimage)

            ## +++++++++++++++++++++ create Gaussian kernel map
            print("########## psfcat2gauskerwithtweak_no_recentre info start ##########\n", file=outLog)
            print("########## psfcat2gauskerwithtweak_no_recentre error start ##########\n", file=errLog)
            proc = subprocess.run(f'(echo -1; cat {PSFcataFile}) | {self._GAAP}/bin/psfcat2gauskerwithtweak_no_recentre > tmp_ker.sh2',
                    stdout=outLog, stderr=errLog, shell=True, cwd=tmp_dir)
            print("########## psfcat2gauskerwithtweak_no_recentre info end ##########\n", file=outLog)
            print("########## psfcat2gauskerwithtweak_no_recentre error end ##########\n", file=errLog)

            print("########## fitkermaptwk_noplots info start ##########\n", file=outLog)
            print("########## fitkermaptwk_noplots error start ##########\n", file=errLog)
            proc = subprocess.run(f'{self._GAAP}/bin/fitkermaptwk_noplots < tmp_ker.sh2 > {SKYimaKerFile}',
                    stdout=outLog, stderr=errLog, shell=True, cwd=tmp_dir)
            print("########## fitkermaptwk_noplots info end ##########\n", file=outLog)
            print("########## fitkermaptwk_noplots error end ##########\n", file=errLog)

            logger.info("Kernel map produced as {:}".format(SKYimaKerFile))

            ## +++++++++++++++++++++ convolve to image (Gaussianised image)
            print("########## imxshmapwithtweak info start ##########\n", file=outLog)
            print("########## imxshmapwithtweak error start ##########\n", file=errLog)
            proc = subprocess.run(f'{self._GAAP}/bin/imxshmapwithtweak < {SKYimaKerFile}',
                    stdout=outLog, stderr=errLog, shell=True, cwd=tmp_dir)
            print("########## imxshmapwithtweak info end ##########\n", file=outLog)
            print("########## imxshmapwithtweak error end ##########\n", file=errLog)

            ### rename as gpsfimage
            os.rename(os.path.join(tmp_dir, 'convolved.fits'), SKYimaGpsfFile)

            logger.info("PSF-Gaussianised image produced as {:}".format(SKYimaGpsfFile))

        return SKYimaKerFile, SKYimaGpsfFile

    def _GAaPphotFunc(self, SKYimaFile, SKYcata, SKYimaKerFile, SKYimaGpsfFile, tmp_dir, SKYweiFile=None, outLog=subprocess.PIPE, errLog=subprocess.STDOUT):
        """
        calculate GAaP photometry
        Make input catalogue for Gaap photometry from the galaxy catalogue:
        want text file with X,Y,A",B",PAworld,ID. where PAworld = -THETA_WORLD
        """

        logger.info("Calculating GAaP photometry...")

        # products
        ## intermediate
        SKYimaKerAcfFile = os.path.join(tmp_dir, os.path.basename(SKYimaFile).replace('.fits', '.keracf.map'))
        SKYimaTotAcfFile = os.path.join(tmp_dir, os.path.basename(SKYimaFile).replace('.fits', '.totacf.map'))
        ## final
        GaapFile = os.path.join(self._outDir, os.path.basename(SKYimaFile).replace('.fits', '.gaap'))

        # check if Gaap already exist
        if os.path.isfile(GaapFile) and (os.path.getsize(GaapFile)>0):
            logger.info("GAaP catalogue already exists.")

        else:

            ## useful parameters: needs X,Y,A",B",PAworld,ID
            radec = SKYcata[['X_WORLD', 'Y_WORLD']].values
            a = SKYcata['A_WORLD'].values*3600. # arcsec
            b = SKYcata['B_WORLD'].values*3600. # arcsec
            PAworld = -1.*SKYcata['THETA_WORLD'].values
            id_obj = SKYcata['NUMBER'].values

            ## position
            np.savetxt(os.path.join(tmp_dir, 'radec.cat'), radec)
            proc = subprocess.run('sky2xy {:} @radec.cat | awk \' {:} \' > xy.cat'.format(SKYimaGpsfFile, '{print($5,$6)}'),
                    stdout=subprocess.PIPE, stderr=subprocess.STDOUT, shell=True, cwd=tmp_dir)
            tmp = np.loadtxt(os.path.join(tmp_dir, 'xy.cat'))
            x = tmp[:, 0]
            y = tmp[:, 1]

            ## aperture
            aper = self._minaper
            maxaper = self._maxaper
            ### A"
            a = np.sqrt(aper**2+a**2)
            a[a>maxaper] = maxaper
            ### B"
            b = np.sqrt(aper**2+b**2)
            b[b>maxaper] = maxaper

            ## inputs for gaap
            gaapin = np.transpose(np.vstack((x, y, a, b, PAworld, id_obj)))
            np.savetxt(os.path.join(tmp_dir, 'gaap.in'), gaapin, fmt='%.3f %.3f %.2f %.2f %.2f %i')

            ## Make kernel ACF map of GPSF image if needed
            print("########## keracfmap & fitkermap info start ##########\n", file=outLog)
            print("########## keracfmap & fitkermap error start ##########\n", file=errLog)
            proc = subprocess.run(f'{self._GAAP}/bin/keracfmap < {SKYimaKerFile} | {self._GAAP}/bin/fitkermap > {SKYimaKerAcfFile}',
                    stdout=outLog, stderr=errLog, shell=True, cwd=tmp_dir)
            print("########## keracfmap & fitkermap info end ##########\n", file=outLog)
            print("########## keracfmap & fitkermap error end ##########\n", file=errLog)

            logger.info("Kernel acf map produced as {:}".format(SKYimaKerAcfFile))

            ## Estimate pre-GPSF pixel correlation and convolve with kernel ACF
            if os.path.isfile(SKYimaTotAcfFile) and (os.path.getsize(SKYimaTotAcfFile)>0):
                logger.info("Total noise acf map already exists.")
            else:
                inimage = os.path.join(tmp_dir, 'inimage.fits')
                if os.path.isfile(inimage):
                    os.remove(inimage)
                os.symlink(SKYimaFile, inimage)
                if (SKYweiFile is not None):
                    wei_image = os.path.join(tmp_dir, 'weights.fits')
                    if os.path.isfile(wei_image):
                        os.remove(wei_image)
                    os.symlink(SKYweiFile, wei_image)

                print("########## kermapxgauss info start ##########\n", file=outLog)
                print("########## kermapxgauss error start ##########\n", file=errLog)
                proc = subprocess.run('echo NOMAP|{GAAP}/bin/gapphot > /dev/null;\
                                        mv -f cov.fits {image_name}.cov.fits;\
                                        mv -f covsh.fits {image_name}.covsh.fits;\
                                        (echo `head -1 keracffitted.map |awk \' {print2} \'` `head -1 {keracfmap} |awk \' {print2} \'`; cat {keracfmap}) |{GAAP}/bin/kermapxgauss > {totacfmap};\
                                        '.format(GAAP=self._GAAP,
                                                print2='{print($2)}',
                                                keracfmap=SKYimaKerAcfFile,
                                                image_name=os.path.basename(SKYimaFile),
                                                totacfmap=SKYimaTotAcfFile),
                                        stdout=outLog, stderr=errLog, shell=True, cwd=tmp_dir)
                print("########## kermapxgauss info end ##########\n", file=outLog)
                print("########## kermapxgauss error end ##########\n", file=errLog)

                logger.info("Total noise acf map produced as {:}".format(SKYimaTotAcfFile))

            ## do gapphot
            inimage = os.path.join(tmp_dir, 'inimage.fits')
            if os.path.isfile(inimage):
                os.remove(inimage)
            os.symlink(SKYimaGpsfFile, inimage)

            print("########## gapphot info start ##########\n", file=outLog)
            print("########## gapphot error start ##########\n", file=errLog)
            proc = subprocess.run('{GAAP}/bin/dfits {gpsfimage} | grep GPSFSIG | sed \'s/.*GPSFSIG =//\' | awk \' {print1} \' > gpsfsig.dat;\
                                (echo {totacfmap} ; cat gaap.in) | {GAAP}/bin/gapphot > {gaapout};\
                                mv -f keracf.fits {totacfmap}.fits;\
                                '.format(GAAP=self._GAAP,
                                        gpsfimage=SKYimaGpsfFile,
                                        print1='{print($1)}',
                                        totacfmap=SKYimaTotAcfFile,
                                        gaapout=GaapFile),
                                stdout=outLog, stderr=errLog, shell=True, cwd=tmp_dir)
            print("########## gapphot info start ##########\n", file=outLog)
            print("########## gapphot error start ##########\n", file=errLog)

            logger.info("GAaP catalogue produced as {:}".format(GaapFile))

        return GaapFile

    def _CleanUpFunc(self, GaapFile_list, tmp_dir):
        """
        Clean up intermediate files
        """

        if self._clean_up_level >= 1:
            if os.path.exists(tmp_dir):
                logger.info('Clean up tmp directory')
                shutil.rmtree(tmp_dir)

        if self._clean_up_level >= 2:
            logger.info('Clean up .gaap directory')
            for f in GaapFile_list:
                os.remove(f)

    def _CombineCataFunc(self, SKYcata, bands, GaapFile_list, FinalFile):
        """
        combine all GAaP results with detection id to one feather catalogue
        """

        logger.info('Combine original GAaP outputs...')

        # the final catalogue
        data_out = pd.DataFrame({'id_detec': np.array(SKYcata['NUMBER']).astype(int),
                                f'mask_meas_{len(bands)}bands': np.zeros(len(SKYcata)).astype(int)})

        # loop over all bands
        for i_band, band in enumerate(bands):

            # GAaP catalogue
            GaapFile = GaapFile_list[i_band]
            tmp = np.loadtxt(GaapFile)
            data_out.loc[:, 'Agaper'] = tmp[:, 2]
            data_out.loc[:, 'Bgaper'] = tmp[:, 3]
            data_out.loc[:, 'PAgaap'] = tmp[:, 4]
            data_out.loc[:, 'FLUX_GAAP_0p7_{:}'.format(band)] = tmp[:, 5]
            data_out.loc[:, 'FLUXERR_GAAP_0p7_{:}'.format(band)] = tmp[:, 6]
            data_out.loc[:, 'FLAG_GAAP_0p7_{:}'.format(band)] = tmp[:, 8]

            ## magnitude
            if band in ['u', 'g', 'r', 'i']:
                ## criterion: SNR >= 1 and flux>0 and flux_err>0
                mask_ture = (tmp[:, 5]>=tmp[:, 6]) & (tmp[:, 5]>0) & (tmp[:, 6]>0)
                mask_false = np.invert(mask_ture)
            else:
                ## criterion: flux>0 and flux_err>0
                mask_ture = (tmp[:, 5]>0) & (tmp[:, 6]>0)
                mask_false = np.invert(mask_ture)

            ## preserved results
            data_out.loc[mask_ture, 'MAG_GAAP_0p7_{:}'.format(band)] = self._magzero-2.5*np.log10(tmp[mask_ture, 5])
            data_out.loc[mask_ture, 'MAGERR_GAAP_0p7_{:}'.format(band)] = 1.0857362047581294 * tmp[mask_ture, 6]/tmp[mask_ture, 5]

            ## assigning dummy values to undesired flux
            data_out.loc[mask_false, 'MAG_GAAP_0p7_{:}'.format(band)] = 99.
            data_out.loc[mask_false, 'MAGERR_GAAP_0p7_{:}'.format(band)] = self._1sigma_limits[band]

            # mask undesired flux
            data_out.loc[mask_false, f'mask_meas_{len(bands)}bands'] += 1

        data_out.to_feather(FinalFile)

        logger.info('Combined catalogues saved to {:}'.format(FinalFile))

    def RunSingleTile(self, FinalFile, SKYcataFile, bands, SKYimaFile_list, SKYweiFile_list=None, PSFimaFile_list=None, star_info=None):
        '''
        Running GAaP with single process for multiple bands.
        '''

        name_base = re.findall(f'(.*).feather', os.path.basename(FinalFile))[0]

        logger.info(f'Running GAaP for {name_base}...')

        # check if already exist
        if os.path.isfile(FinalFile):
            logger.info(f"Final combined feather catalogue already exist, end for {name_base}.")
            return 1

        # prepare log info
        if self._running_log:
            basename = os.path.basename(FinalFile)
            outLog = open(os.path.join(self._logDir, f"{name_base}.log"), 'w')
            errLog = open(os.path.join(self._logDir, f"{name_base}.err.log"), 'w')
        else:
            outLog = subprocess.PIPE
            errLog = subprocess.STDOUT

        # load the detection catalogue
        SKYcata = pd.read_feather(SKYcataFile)

        # 0. make tmp dir
        ## different for each run
        tmp_dir = os.path.join(self._outDir, name_base)
        if not os.path.exists(tmp_dir):
            os.mkdir(tmp_dir)

        # loop over bands
        GaapFile_list = []
        for i_band, band in enumerate(bands):

            SKYimaFile = SKYimaFile_list[i_band]

            if (SKYweiFile_list is not None):
                SKYweiFile = SKYweiFile_list[i_band]
            else:
                SKYweiFile = None

            if (PSFimaFile_list is not None):
                PSFimaFile = PSFimaFile_list[i_band]
            else:
                PSFimaFile = None

            logger.info(f'Running for band {band}...')

            # 1. make star catalogue in right SExtractor ascii format
            if (PSFimaFile is not None):
                PSFcataFile = self._PSFcataFunc(band, PSFimaFile, star_info=None, outLog=outLog, errLog=errLog)
            else:
                PSFcataFile = self._PSFcataFunc(band, SKYimaFile, star_info=star_info, outLog=outLog, errLog=errLog)

            # 2. calculate Gaussianization kernel and Gaussianised image
            SKYimaKerFile, SKYimaGpsfFile = self._GaussianizationFunc(SKYimaFile, PSFcataFile, tmp_dir, outLog=outLog, errLog=errLog)

            # 3. run gaap
            GaapFile = self._GAaPphotFunc(SKYimaFile, SKYcata, SKYimaKerFile, SKYimaGpsfFile, tmp_dir, SKYweiFile=SKYweiFile, outLog=outLog, errLog=errLog)
            GaapFile_list.append(GaapFile)

            logger.info(f'Finished for band {band}')

        if self._running_log:
            outLog.close()
            errLog.close()

        self._CombineCataFunc(SKYcata, bands, GaapFile_list, FinalFile)

        if self._clean_up_level:
            self._CleanUpFunc(GaapFile_list, tmp_dir)

        logger.info(f'Finished for tile {name_base}.')
