# -*- coding: utf-8 -*-
# @Author: lshuns
# @Date:   2021-09-17 21:28:17
# @Last Modified by:   lshuns
# @Last Modified time: 2022-06-03 09:50:30

### Wrapper for GAaP code

import os
import re
import sys
import glob
import time
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
    tmp_dir : str, optional (default: ==out_dir)
        Directory for tmp files.
    star_SNR_cut : list of float, optional (default: [100, 1000])
        Within which SNR range stars are used.
    mag_zero : float, optional (default: 30.)
        Zero point for the magnitude.
    min_aper_list : a list of float, optional (default: [0.7, 1.0])
        Minimum aperture sizes.
    max_aper : float, optional (default: 2.0)
        Maximum aperture size.
    spatial_variation : dictionary
        does psf vary spatially? 
    running_log : bool, optional (default: True)
        save running info from external code
        """

    def __init__(self, gaap_dir, tmp_dir, 
                star_SNR_cut=[100., 1000.],
                mag_zero=30., detection_band='r',
                min_aper_list=[0.7, 1.0], max_aper=2.0,
                spatial_variation={},
                running_log=True, log_dir=None):

        logger.info("Initialising the GAaP wrapper...")

        self._GAAP = gaap_dir
        logger.info("GAaP directory: %s", self._GAAP)

        self._tmpDir = tmp_dir
        logger.info("tmp directory: %s", self._tmpDir)
        if not os.path.exists(self._tmpDir):
            os.mkdir(self._tmpDir)

        self._star_SNR_cut = star_SNR_cut
        logger.info(f'star SNR cut: ({star_SNR_cut[0]}, {star_SNR_cut[1]})')

        self._magzero = mag_zero
        logger.info("magnitude zero point: %f", self._magzero)

        self._detection_band = detection_band
        logger.info("detection band: %s", self._detection_band)

        self._minaper_list = min_aper_list
        if len(self._minaper_list) > 2:
            raise Exception('Number of minaper sizes should be either 1 or 2!')
        logger.info(f"Minimum apertures: {self._minaper_list}")
        self._maxaper = max_aper
        logger.info("Maximum aperture: %f", self._maxaper)

        self._spatial_variation = spatial_variation
        logger.info(f'model contains spatial variation: {spatial_variation}')

        self._running_log = running_log
        logger.info(f"Save running info: {self._running_log}")

        self._logDir = log_dir

    def _PSFcataFunc(self, band, SKYimaFile, tmp_dir, star_info=None, outLog=subprocess.PIPE, errLog=subprocess.STDOUT):
        """
        make PSF catalogue in right SExtractor ascii format
        """

        # product
        PSFcataFile = os.path.join(tmp_dir, os.path.basename(SKYimaFile).replace('.fits', '_psf.cat'))

        # using psf map
        if star_info is None:

            logger.info("Producing PSF catalogue from separated PSF map...")

            # necessary sex config files
            sex_config = []
            for file in ['gaap.sex', 'gaap.sexpar', 'default.conv']:
                sex_config.append(os.path.join(self._GAAP, 'lib/{:}'.format(file)))

            # run sex
            cmd = ['sex', SKYimaFile, '-c', sex_config[0], '-PARAMETERS_NAME', sex_config[1], '-FILTER_NAME', sex_config[2], '-DETECT_THRESH', str(10), '-CATALOG_NAME', PSFcataFile]
            proc = subprocess.run(cmd, stdout=outLog, stderr=errLog)

        # using stars
        else:

            logger.info("Producing PSF catalogue using stars within image...")

            # intermediate
            DetecFile = os.path.join(tmp_dir, os.path.basename(SKYimaFile).replace('.fits', '_all.cat'))

            # +++++++++ sextractor detection

            # necessary sex config files
            sex_config = []
            for file in ['gaap.sex', 'gaap_ima.sexpar', 'default.conv']:
                sex_config.append(os.path.join(self._GAAP, 'lib/{:}'.format(file)))

            # run sex
            cmd = ['sex', SKYimaFile, '-c', sex_config[0], '-PARAMETERS_NAME', sex_config[1], '-FILTER_NAME', sex_config[2], '-DETECT_THRESH', str(10), '-CATALOG_NAME', DetecFile]
            proc = subprocess.run(cmd, stdout=outLog, stderr=errLog)

            # +++++++++ cross match to find stars

            ## load detected catalogue
            detec_cata_ori = np.loadtxt(DetecFile)

            ## select high SNR objects
            SNR_tmp = detec_cata_ori[:, 2]/detec_cata_ori[:, 3]
            mask_good = (SNR_tmp > self._star_SNR_cut[0]) & (SNR_tmp < self._star_SNR_cut[1])
            detec_cata_ori = detec_cata_ori[mask_good]
            del SNR_tmp, mask_good

            ## cross match
            detec_cata = pd.DataFrame({'NUMBER': detec_cata_ori[:, 8].astype(int),
                                        'X_WORLD': detec_cata_ori[:, 11],
                                        'Y_WORLD': detec_cata_ori[:, 12],
                                        'MAG_AUTO': self._magzero-2.5*np.log10(detec_cata_ori[:, 2])})
            id_list = ['index_input', 'NUMBER']
            position_list = [['RA_input', 'DEC_input'], ['X_WORLD', 'Y_WORLD']]
            mag_list = [f'{band}_input', 'MAG_AUTO']
            matched_cata, _, _ = run_position2id(star_info, detec_cata, id_list, position_list, mag_list,
                                outDir=None, basename=None, save_matched=False, save_false=False, save_missed=False,
                                r_max=0.6/3600., k=2, mag_closest=True, running_info=False)
            del detec_cata

            ## get stars
            id_stars = matched_cata['id_detec'].values
            mask = np.isin(detec_cata_ori[:, 8], id_stars)
            psf_cata = detec_cata_ori[mask, :11]
            del mask, detec_cata_ori
            N_stars = len(psf_cata)
            if (self._spatial_variation[band]) and (N_stars < 1000):
                raise Exception(f"too few stars {N_stars}, please select a larger range of SNR!")
            if (~self._spatial_variation[band]) and (N_stars < 100):
                raise Exception(f"too few stars {N_stars}, please select a larger range of SNR!")
            logger.info(f"Number of stars selected {N_stars}")

            # +++++++++ save desired info
            np.savetxt(PSFcataFile, psf_cata,
                        fmt='%.4f %.4f %.2f %.5f %.3f %.3f %.3f %.2f %i %i %.2f')
            del psf_cata

        logger.debug(f"PSF catalogue Extracted as {PSFcataFile}")

        return PSFcataFile

    def _GaussianizationFunc(self, band, SKYimaFile, PSFcataFile, tmp_dir, outLog=subprocess.PIPE, errLog=subprocess.STDOUT):
        """
        calculate Gaussianization kernel and Gaussianised image
        allow Gaussianised PSF radius to be determined autimatically
        """

        logger.info("Producing PSF-Gaussianised image...")

        # products
        SKYimaKerFile = os.path.join(tmp_dir, os.path.basename(SKYimaFile).replace('.fits', '.ker.map'))
        SKYimaGpsfFile = os.path.join(tmp_dir, os.path.basename(SKYimaFile).replace('.fits', '.gpsf.fits'))

        # orders
        if self._spatial_variation[band]:
            ### kernel shapelet order 10, spatial order 4 as used by data
            np.savetxt(os.path.join(tmp_dir, 'orders.par'), [[10], [4]], fmt='%d')
        else:
            ### kernel has shapelet order 10, no spatial variation for the PSF
            np.savetxt(os.path.join(tmp_dir, 'orders.par'), [[10], [0]], fmt='%d')

        # link to the image
        inimage = os.path.join(tmp_dir, 'inimage.fits')
        if os.path.isfile(inimage):
            os.remove(inimage)
        os.symlink(SKYimaFile, inimage)

        ## +++++++++++++++++++++ create Gaussian kernel map

        # psfcat2gauskerwithtweak_no_recentre
        proc = subprocess.run(f'(echo -1; cat {PSFcataFile}) | {self._GAAP}/bin/psfcat2gauskerwithtweak_no_recentre > tmp_ker.sh2',
                stdout=outLog, stderr=errLog, shell=True, cwd=tmp_dir)
        del PSFcataFile

        # fitkermaptwk_noplots
        proc = subprocess.run(f'{self._GAAP}/bin/fitkermaptwk_noplots < tmp_ker.sh2 > {SKYimaKerFile}',
                stdout=outLog, stderr=errLog, shell=True, cwd=tmp_dir)

        logger.debug("Kernel map produced as {:}".format(SKYimaKerFile))

        ## +++++++++++++++++++++ convolve to image (Gaussianised image)

        # imxshmapwithtweak
        proc = subprocess.run(f'{self._GAAP}/bin/imxshmapwithtweak < {SKYimaKerFile}',
                stdout=outLog, stderr=errLog, shell=True, cwd=tmp_dir)

        # rename as gpsfimage
        os.rename(os.path.join(tmp_dir, 'convolved.fits'), SKYimaGpsfFile)

        logger.debug("PSF-Gaussianised image produced as {:}".format(SKYimaGpsfFile))

        return SKYimaKerFile, SKYimaGpsfFile

    def _GAaPphotFunc(self, band, GaapFile, SKYimaFile, SKYcata, SKYimaKerFile, SKYimaGpsfFile, 
                            tmp_dir, minaper,
                            SKYweiFile=None, outLog=subprocess.PIPE, errLog=subprocess.STDOUT):
        """
        calculate GAaP photometry
        Make input catalogue for Gaap photometry from the galaxy catalogue:
        want text file with X,Y,A",B",PAworld,ID. where PAworld = -THETA_WORLD
        """

        logger.info("Calculating GAaP photometry...")

        # intermediate products
        SKYimaKerAcfFile = os.path.join(tmp_dir, os.path.basename(SKYimaFile).replace('.fits', '.keracf.map'))
        SKYimaTotAcfFile = os.path.join(tmp_dir, os.path.basename(SKYimaFile).replace('.fits', '.totacf.map'))

        # useful parameters: needs X,Y,A",B",PAworld,ID
        radec = SKYcata[['X_WORLD', 'Y_WORLD']].values
        a = SKYcata['A_WORLD'].values*3600. # arcsec
        b = SKYcata['B_WORLD'].values*3600. # arcsec
        PAworld = -1.*SKYcata['THETA_WORLD'].values
        id_obj = SKYcata['NUMBER'].values
        del SKYcata

        # orders
        if self._spatial_variation[band]:
            ### kernel shapelet order 10, spatial order 4 as used by data
            np.savetxt(os.path.join(tmp_dir, 'orders.par'), [[10], [4]], fmt='%d')
        else:
            ### kernel has shapelet order 10, no spatial variation for the PSF
            np.savetxt(os.path.join(tmp_dir, 'orders.par'), [[10], [0]], fmt='%d')

        # position
        np.savetxt(os.path.join(tmp_dir, 'radec.cat'), radec)
        proc = subprocess.run('sky2xy {:} @radec.cat | awk \' {:} \' > xy.cat'.format(SKYimaGpsfFile, '{print($5,$6)}'),
                stdout=subprocess.PIPE, stderr=subprocess.STDOUT, shell=True, cwd=tmp_dir)
        tmp = np.loadtxt(os.path.join(tmp_dir, 'xy.cat'))
        x = tmp[:, 0]
        y = tmp[:, 1]
        del tmp

        # aperture
        maxaper = self._maxaper
        ### A"
        a = np.sqrt(minaper**2+a**2)
        a[a>maxaper] = maxaper
        ### B"
        b = np.sqrt(minaper**2+b**2)
        b[b>maxaper] = maxaper

        # inputs for gaap
        gaapin = np.transpose(np.vstack((x, y, a, b, PAworld, id_obj)))
        np.savetxt(os.path.join(tmp_dir, 'gaap.in'), gaapin, fmt='%.3f %.3f %.2f %.2f %.2f %i')
        del radec, x, y, a, b, PAworld, id_obj

        # >>> 1. Make kernel ACF map of GPSF image if needed
        proc = subprocess.run(f'{self._GAAP}/bin/keracfmap < {SKYimaKerFile} | {self._GAAP}/bin/fitkermap > {SKYimaKerAcfFile}',
                stdout=outLog, stderr=errLog, shell=True, cwd=tmp_dir)

        logger.debug("Kernel acf map produced as {:}".format(SKYimaKerAcfFile))

        # >>> 2. Estimate pre-GPSF pixel correlation and convolve with kernel ACF
        inimage = os.path.join(tmp_dir, 'inimage.fits')
        if os.path.isfile(inimage):
            os.remove(inimage)
        os.symlink(SKYimaFile, inimage)
        if (SKYweiFile is not None):
            wei_image = os.path.join(tmp_dir, 'weights.fits')
            if os.path.isfile(wei_image):
                os.remove(wei_image)
            os.symlink(SKYweiFile, wei_image)
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
        logger.debug("Total noise acf map produced as {:}".format(SKYimaTotAcfFile))

        # >>> 3. do gapphot
        inimage = os.path.join(tmp_dir, 'inimage.fits')
        if os.path.isfile(inimage):
            os.remove(inimage)
        os.symlink(SKYimaGpsfFile, inimage)
        proc = subprocess.run('{GAAP}/bin/dfits {gpsfimage} | grep GPSFSIG | sed \'s/.*GPSFSIG =//\' | awk \' {print1} \' > gpsfsig.dat;\
                            (echo {totacfmap} ; cat gaap.in) | {GAAP}/bin/gapphot > {gaapout};\
                            mv -f keracf.fits {totacfmap}.fits;\
                            '.format(GAAP=self._GAAP,
                                    gpsfimage=SKYimaGpsfFile,
                                    print1='{print($1)}',
                                    totacfmap=SKYimaTotAcfFile,
                                    gaapout=GaapFile),
                            stdout=outLog, stderr=errLog, shell=True, cwd=tmp_dir)
        logger.debug("GAaP catalogue produced as {:}".format(GaapFile))

    def _CombineCataFunc(self, SKYcata, GaapFile_list_dic, FinalFile):
        """
        combine all GAaP results with detection id to one feather catalogue
        """

        logger.info('Combine original GAaP outputs...')

        # the final catalogue
        data_out = pd.DataFrame({'id_detec': np.array(SKYcata['NUMBER']).astype(int)})
        del SKYcata

        # loop over all bands
        for band, GaapFile_list in GaapFile_list_dic.items():

            # loop over apertures
            for i_aper, minaper in enumerate(self._minaper_list):

                # the GaapFile
                GaapFile = GaapFile_list[i_aper]
                str_minaper = str(minaper).replace('.', 'p')

                # GAaP catalogue
                tmp = np.loadtxt(GaapFile)
                if band == self._detection_band:
                    data_out.loc[:, f'Agaper_{str_minaper}'] = tmp[:, 2]
                    data_out.loc[:, f'Bgaper_{str_minaper}'] = tmp[:, 3]
                    data_out.loc[:, f'PAgaap_{str_minaper}'] = tmp[:, 4]
                ## the flux and flag
                flux = np.float64(tmp[:, 5])
                fluxerr = np.float64(tmp[:, 6])
                flag = np.int32(tmp[:, 8])
                del tmp
                ## save
                data_out.loc[:, f'FLUX_GAAP_{str_minaper}_{band}'] = flux
                data_out.loc[:, f'FLUXERR_GAAP_{str_minaper}_{band}'] = fluxerr
                data_out.loc[:, f'FLAG_GAAP_{str_minaper}_{band}'] = flag

                # >>> 1. convert to magnitude
                ### the VISTA
                if band in ['Z', 'Y', 'J', 'H', 'Ks']:
                    ##### good measurement
                    mask_tmp = (flux>0) & (flag==0)
                    data_out.loc[mask_tmp, f'MAG_GAAP_{str_minaper}_{band}'] \
                                    = self._magzero - 2.5*np.log10(flux[mask_tmp])
                    data_out.loc[mask_tmp, f'MAGERR_GAAP_{str_minaper}_{band}'] \
                                    = 1.0857362047581294 * (fluxerr[mask_tmp]/flux[mask_tmp])
                    del mask_tmp
                    ##### non detections
                    mask_tmp = (flux<0) & (flag==0)
                    data_out.loc[mask_tmp, f'MAG_GAAP_{str_minaper}_{band}'] \
                                    = 99.
                    data_out.loc[mask_tmp, f'MAGERR_GAAP_{str_minaper}_{band}'] \
                                    = np.abs(1.0857362047581294 * (fluxerr[mask_tmp]/flux[mask_tmp]))
                    del mask_tmp
                    ###### failures
                    mask_tmp = (flag!=0) | (flux==0.)
                    data_out.loc[mask_tmp, f'MAG_GAAP_{str_minaper}_{band}'] \
                                    = -99.
                    data_out.loc[mask_tmp, f'MAGERR_GAAP_{str_minaper}_{band}'] \
                                    = -99.
                    del mask_tmp
                elif band in ['u', 'g', 'r', 'i']:
                    ##### good measurement
                    mask_tmp = (flux>0.) & (fluxerr>=0) & (flux>=fluxerr)
                    data_out.loc[mask_tmp, f'MAG_GAAP_{str_minaper}_{band}'] \
                                    = self._magzero - 2.5*np.log10(flux[mask_tmp])
                    data_out.loc[mask_tmp, f'MAGERR_GAAP_{str_minaper}_{band}'] \
                                    = 1.0857362047581294 * (fluxerr[mask_tmp]/flux[mask_tmp])
                    del mask_tmp
                    ##### non detections
                    mask_tmp = (flux<fluxerr) & (fluxerr>0) & (flux!=0.)
                    data_out.loc[mask_tmp, f'MAG_GAAP_{str_minaper}_{band}'] \
                                    = 99.
                    data_out.loc[mask_tmp, f'MAGERR_GAAP_{str_minaper}_{band}'] \
                                    = np.abs(1.0857362047581294 * (fluxerr[mask_tmp]/flux[mask_tmp]))
                    del mask_tmp
                    ###### failures
                    mask_tmp = (flux==0.) | (fluxerr<0.)
                    data_out.loc[mask_tmp, f'FLAG_GAAP_{str_minaper}_{band}'] \
                                    = 1
                    data_out.loc[mask_tmp, f'MAG_GAAP_{str_minaper}_{band}'] \
                                    = -99.
                    data_out.loc[mask_tmp, f'MAGERR_GAAP_{str_minaper}_{band}'] \
                                    = -99.
                    del mask_tmp
                else:
                    raise Exception(f'Unrecognised band {band}')

                # >>> 2. the limiting magnitudes
                mask_ture = (flux!=0) & (fluxerr>0)
                mask_false = np.invert(mask_ture)
                data_out.loc[mask_ture, f'MAG_LIM_{str_minaper}_{band}'] = self._magzero - 2.5*np.log10(fluxerr[mask_ture])
                data_out.loc[mask_false, f'MAG_LIM_{str_minaper}_{band}'] = -99.
                del mask_ture, mask_false, flux, fluxerr, flag

        # take decision, which aperture to use
        ## get min and max aper
        str_minaper_min = str(min(self._minaper_list)).replace('.', 'p')
        str_minaper_max = str(max(self._minaper_list)).replace('.', 'p')
        ## initialize the decision matrix
        R = np.zeros((len(data_out), len(GaapFile_list_dic.keys())))
        for i_band, band in enumerate(GaapFile_list_dic.keys()):

            # calculate R
            fluxerr1 = data_out[f'FLUXERR_GAAP_{str_minaper_max}_{band}'].values
            fluxerr2 = data_out[f'FLUXERR_GAAP_{str_minaper_min}_{band}'].values
            R[:, i_band] = fluxerr1/fluxerr2
            R[(fluxerr1==-1), i_band] = 1.
            del fluxerr1, fluxerr2

        # decide whether the larger aperture should be used
        large_aper = ( np.min(R, axis=1) < 1./np.max(R, axis=1) ) | ( np.max(R, axis=1) < 0 )
        del R
        for band in GaapFile_list_dic.keys():

            for col_main in ['FLUX_GAAP', 'FLUXERR_GAAP', 'FLAG_GAAP', 'MAG_GAAP', 'MAGERR_GAAP', 'MAG_LIM']:

                data_out.loc[:, f'{col_main}_{band}'] \
                    = np.where(large_aper, 
                            data_out[f'{col_main}_{str_minaper_max}_{band}'].values, 
                            data_out[f'{col_main}_{str_minaper_min}_{band}'].values) 

        # the same for aperture para
        for col_main in ['Agaper', 'Bgaper', 'PAgaap']:
                data_out.loc[:, f'{col_main}'] \
                    = np.where(large_aper, 
                            data_out[f'{col_main}_{str_minaper_max}'].values, 
                            data_out[f'{col_main}_{str_minaper_min}'].values) 

        del large_aper

        # save
        tmp_FinalFile = FinalFile + '_tmp'
        data_out.to_feather(tmp_FinalFile)
        del data_out
        os.rename(tmp_FinalFile, FinalFile)
        logger.debug('Combined catalogues saved to {:}'.format(FinalFile))

    def _RunSingleBand(self, queue, SKYcata, band, SKYimaFile, SKYweiFile=None, PSFimaFile=None, star_info=None):
        '''
        Running GAaP with single process for single band.
        '''
        logger.info(f'Running for band {band}...')

        SKYimaFile_name = os.path.basename(SKYimaFile)

        # prepare log info
        if self._running_log:
            outLog = open(os.path.join(self._logDir, SKYimaFile_name.replace('.fits', '.log')), 'w')
            errLog = open(os.path.join(self._logDir, SKYimaFile_name.replace('.fits', '.err.log')), 'w')
        else:
            outLog = subprocess.PIPE
            errLog = subprocess.STDOUT

        # 0. make tmp dir
        tmp_dir = os.path.join(self._tmpDir, 'tmp_gaap_' + re.findall(f'(.*).fits', SKYimaFile_name)[0])
        ## to avoid semi-finished product
        if os.path.exists(tmp_dir):
            shutil.rmtree(tmp_dir)
        os.mkdir(tmp_dir)

        # 1. make star catalogue in right SExtractor ascii format
        if (PSFimaFile is not None):
            PSFcataFile = self._PSFcataFunc(band, PSFimaFile, tmp_dir, star_info=None, outLog=outLog, errLog=errLog)
        else:
            PSFcataFile = self._PSFcataFunc(band, SKYimaFile, tmp_dir, star_info=star_info, outLog=outLog, errLog=errLog)

        # 2. calculate Gaussianization kernel and Gaussianised image
        SKYimaKerFile, SKYimaGpsfFile = self._GaussianizationFunc(band, SKYimaFile, PSFcataFile, tmp_dir, outLog=outLog, errLog=errLog)

        # 3. run gaap
        ## loop over different apertures
        GaapFile_list = []
        for minaper in self._minaper_list:
            logger.info(f'     for aperture {minaper}...')

            # where to save tmp files
            tmp_minaper_dir = os.path.join(tmp_dir, 'minaper_' + str(minaper).replace('.', 'p'))
            ## to avoid semi-finished product
            if os.path.exists(tmp_minaper_dir):
                shutil.rmtree(tmp_minaper_dir)
            os.mkdir(tmp_minaper_dir)

            # final output
            GaapFile = os.path.join(tmp_minaper_dir, SKYimaFile_name.replace('.fits', '.gaap'))
            GaapFile_list.append(GaapFile)

            # running
            self._GAaPphotFunc(band, GaapFile, SKYimaFile, SKYcata, SKYimaKerFile, SKYimaGpsfFile, 
                                    tmp_minaper_dir, minaper,
                                    SKYweiFile=SKYweiFile, outLog=outLog, errLog=errLog)

        ## return info
        queue.put([band, GaapFile_list])
        del SKYcata

        if self._running_log:
            outLog.close()
            errLog.close()

        logger.info(f'Finished for band {band}.')

    def RunSingleTile(self, Nmax_proc, FinalFile, SKYcataFile, bands, SKYimaFile_list, SKYweiFile_list=None, PSFimaFile_list=None, star_info=None):
        '''
        Running GAaP with multiple processes for multiple bands.
        '''

        name_base = re.findall(f'(.*).feather', os.path.basename(FinalFile))[0]

        logger.info(f'Running GAaP for {name_base}...')

        # check if already exist
        if os.path.isfile(FinalFile):
            logger.info(f"Final combined feather catalogue already exist, end for {name_base}.")
            return 1

        # load the detection catalogue
        SKYcata = pd.read_feather(SKYcataFile)

        # loop over bands with multiprocess
        queue = mp.Queue()
        p_list = []
        N_tasks = len(bands)
        logger.info(f'Total number of bands: {N_tasks}')
        i_worker = 0
        while True:
            N_running = len(mp.active_children())
            logger.debug(f'Number of running {N_running}')
            if i_worker == N_tasks:
                break
            elif N_running >= Nmax_proc:
                time.sleep(5.)
            else:

                band = bands[i_worker]

                SKYimaFile = SKYimaFile_list[i_worker]

                if (SKYweiFile_list is not None):
                    SKYweiFile = SKYweiFile_list[i_worker]
                else:
                    SKYweiFile = None

                if (PSFimaFile_list is not None):
                    PSFimaFile = PSFimaFile_list[i_worker]
                else:
                    PSFimaFile = None

                p = mp.Process(target=self._RunSingleBand, args=(queue, SKYcata, band, SKYimaFile, SKYweiFile, PSFimaFile, star_info))
                i_worker += 1
                p.start()
                p_list.append(p)
                time.sleep(1.)

        for p in p_list:
            p.join()

        GaapFile_list_dic = {}
        while not queue.empty():
            ret = queue.get()
            GaapFile_list_dic[ret[0]] = ret[1]

        self._CombineCataFunc(SKYcata, GaapFile_list_dic, FinalFile)
        del SKYcata

        logger.info(f'Finished for tile {name_base}.')