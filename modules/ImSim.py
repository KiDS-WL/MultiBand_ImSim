# -*- coding: utf-8 -*-
# @Author: lshuns
# @Date:   2020-12-09 19:21:53
# @Last Modified by:   lshuns
# @Last Modified time: 2021-01-31 20:32:46

### main module of ImSim
###### dependence: 
######      ImSimPSF: everthing about PSF
######      ImSimObject: everything about celestial objects
######      ImSimNoiseBackground: everthing about background noise
######      ImSimKiDS: everything about KiDS observations (instrumental setup & data acquisition procedure)
import ImSimPSF as PSFModule
import ImSimObject as ObjModule
import ImSimNoiseBackground as NoiseModule
import ImSimKiDS as KiDSModule

import re
import os
import sys
import time
import random
import galsim
import logging
import subprocess

import numpy as np
import pandas as pd
import multiprocessing as mp

from astropy.io import fits

logger = logging.getLogger(__name__)

def _PSFNoisySkyImages_simple(para_list):
    '''
    Sky image with Gaussian noise and PSF.
        Simple image without any survey strategy.
        Adjoint to RunParallel_PSFNoisySkyImages
    '''

    (tile_label, band, pixel_scale, rng_seed_band, outpath_image_basename,
        rms, seeing, beta, psf_e,
        g_cosmic,
        gals_info_band, gal_rotation_angles,
        stars_info_band,
        outpath_PSF_basename, N_PSF, sep_PSF,
        save_image_chips, save_image_PSF,
        outpath_dir) = para_list
    logger.info(f'Running for tile {tile_label} band {band}...')

    # warning 
    if save_image_chips:
        logger.warning('Simple tile mode does not produce chips!')

    # outpath
    outpath_image_name_list = [outpath_image_basename + f'_rot{gal_rotation_angle:.0f}.fits' for gal_rotation_angle in gal_rotation_angles]
    ## psf map
    if (outpath_PSF_basename is not None):
        outpath_PSF_name_list = [outpath_PSF_basename + f'_rot{gal_rotation_angle:.0f}.fits' for gal_rotation_angle in gal_rotation_angles]

    # first check if already exist
    outpath_image_exist_list = []
    for outpath_image in outpath_image_name_list:
        if os.path.isfile(outpath_image) and (os.path.getsize(outpath_image)>0):
            logger.info(f"{outpath_image} already exist.")
            outpath_image_exist_list.append(True)
        else:
            outpath_image_exist_list.append(False)
    ## psf map
    if (outpath_PSF_basename is not None) :
        outpath_PSF_exist_list = []
        for outpath_PSF_name in outpath_PSF_name_list:
            if (os.path.isfile(outpath_PSF_name) and (os.path.getsize(outpath_PSF_name)>0)):
                logger.info(f"{outpath_PSF_name} already exist.")
                outpath_PSF_exist_list.append(True)
            else:
                outpath_PSF_exist_list.append(False)
    else:
        outpath_PSF_exist_list = [True]        
    ## if all exist, quit
    if (not False in outpath_image_exist_list) and (not False in outpath_PSF_exist_list):
        logger.info("All desired images exist, end the process.")
        return 1

    # +++ PSF
    PSF = PSFModule.MoffatPSF(seeing, moffat_beta=beta, psf_e=psf_e)
    
    # +++ background noise
    noise_list = [NoiseModule.GaussianNoise(rms, rng_seed=int(rng_seed_band+120*gal_rotation_angle)) for gal_rotation_angle in gal_rotation_angles]

    # +++ PSF map
    if (False in outpath_PSF_exist_list):
        mag_PSF_2 = 18. # for noise_flux = 2
        mag_PSF = mag_PSF_2 - 2.5*np.log10(rms/2.)
        image_PSF = PSFModule.PSFmap(PSF, pixel_scale, mag_PSF, N_PSF=N_PSF, sep_PSF=sep_PSF, rng_seed=rng_seed_band)

        for i_rot, outpath_PSF_exist in enumerate(outpath_PSF_exist_list):

            if (not outpath_PSF_exist):
                image_PSF_rot = image_PSF.copy()
                ## noise background
                noise_rot = noise_list[i_rot]
                image_PSF_rot.addNoise(noise_rot)
                ## save
                outpath_PSF_name = outpath_PSF_name_list[i_rot]
                image_PSF_rot.write(outpath_PSF_name)
                logger.info(f"PSF map saved as {outpath_PSF_name}")

    # +++ sky image
    if (False in outpath_image_exist_list):

        # bounds and wcs from galaxy sky positions
        RA_gals = gals_info_band['RA'] # degree
        DEC_gals = gals_info_band['DEC'] # degree
        RA_min = np.min(RA_gals)
        RA_max = np.max(RA_gals)
        DEC_min = np.min(DEC_gals)
        DEC_max = np.max(DEC_gals)
        bounds, wcs = ObjModule.WCS(RA_min, RA_max, DEC_min, DEC_max, pixel_scale)

        # star image
        if (stars_info_band is not None):
            image_stars = ObjModule.StarsImage(bounds, wcs, band, pixel_scale, PSF, stars_info_band)
        else:
            image_stars = None

        # galaxy images
        for i_rot, outpath_image_exist in enumerate(outpath_image_exist_list):

            if (not outpath_image_exist):
                gal_rotation_angle = gal_rotation_angles[i_rot]

                image_galaxies = ObjModule.GalaxiesImage(bounds, wcs, band, pixel_scale, PSF,
                                                gals_info_band, gal_rotation_angle=gal_rotation_angle, g_cosmic=g_cosmic)
                
                ## add stars
                if (image_stars is not None):
                    image_galaxies += image_stars

                ## add noise background
                noise_rot = noise_list[i_rot]
                image_galaxies.addNoise(noise_rot)

                ## save the noisy image
                outpath_image_name = outpath_image_name_list[i_rot]
                image_galaxies.write(outpath_image_name)
                logger.info(f"Noisy sky image saved as {outpath_image_name}")

    return 0

def _PSFNoisySkyImages_KiDS_sameExpo(para_list):
    '''
    Sky image with Gaussian noise and PSF.
        KiDS images with dither and gaps.
        Same noise and psf for all exposures
        Adjoint to RunParallel_PSFNoisySkyImages
    '''

    (tile_label, band, pixel_scale, rng_seed_band, outpath_image_basename,
        rms, seeing, beta, psf_e,
        g_cosmic,
        gals_info_band, gal_rotation_angles,
        stars_info_band,
        outpath_PSF_basename, N_PSF, sep_PSF,
        save_image_chips, save_image_PSF,
        outpath_dir) = para_list
    logger.info(f'Running for tile {tile_label} band {band}...')

    # warning 
    if (save_image_chips) and (band in ['Z', 'Y', 'J', 'H', 'Ks']):
        logger.warning('NIR-band chips in KiDS are not supported now!')

    # # running time
    # start_time = time.time()

    # outpath
    if band in ['u', 'g', 'r', 'i']:
        # 5 exposures
        outpath_image_name_list = [outpath_image_basename + f'_rot{gal_rotation_angle:.0f}_expo{i_expo}.fits' 
                                                                    for gal_rotation_angle in gal_rotation_angles
                                                                    for i_expo in range(5)]
        if (outpath_PSF_basename is not None):
            outpath_PSF_name_list = [outpath_PSF_basename + f'_rot{gal_rotation_angle:.0f}_expo{i_expo}.fits' 
                                                                    for gal_rotation_angle in gal_rotation_angles
                                                                    for i_expo in range(5)]

        # first check if already exist
        outpath_image_exist_list = []
        for outpath_image in outpath_image_name_list:
            if os.path.isfile(outpath_image) and os.path.isfile(outpath_image.replace('.fits', '.weight.fits')):
                logger.info(f"{outpath_image} and its weight images already exist.")
                outpath_image_exist_list.append(True)
            else:
                outpath_image_exist_list.append(False)

    else:
        # only stacked images
        outpath_image_name_list = [outpath_image_basename + f'_rot{gal_rotation_angle:.0f}.fits' for gal_rotation_angle in gal_rotation_angles]
        if (outpath_PSF_basename is not None):
            outpath_PSF_name_list = [outpath_PSF_basename + f'_rot{gal_rotation_angle:.0f}.fits' for gal_rotation_angle in gal_rotation_angles]
        # first check if already exist
        outpath_image_exist_list = []
        for outpath_image in outpath_image_name_list:
            if os.path.isfile(outpath_image):
                logger.info(f"{outpath_image} already exist.")
                outpath_image_exist_list.append(True)
            else:
                outpath_image_exist_list.append(False)

    ## psf map
    if (outpath_PSF_basename is not None) :
        outpath_PSF_exist_list = []
        for outpath_PSF_name in outpath_PSF_name_list:
            if (os.path.isfile(outpath_PSF_name) and (os.path.getsize(outpath_PSF_name)>0)):
                logger.info(f"{outpath_PSF_name} already exist.")
                outpath_PSF_exist_list.append(True)
            else:
                outpath_PSF_exist_list.append(False)
    else:
        outpath_PSF_exist_list = [True]        

    ## if all exist, quit
    if (not False in outpath_image_exist_list) and (not False in outpath_PSF_exist_list):
        logger.info("All desired images exist, end the process.")
        return 1

    # +++ PSF
    PSF = PSFModule.MoffatPSF(seeing, moffat_beta=beta, psf_e=psf_e)
    ## psf images
    if save_image_PSF:                            
        psf_dir_tmp = os.path.join(outpath_dir, f'psf_tile{tile_label}_band{band}')
        if not os.path.exists(psf_dir_tmp):
            os.mkdir(psf_dir_tmp)
        psf_ima = PSFModule.PSFima(PSF, pixel_scale, size=32)
        for i_expo in range(5):
            outpath_tmp = os.path.join(psf_dir_tmp, f'expo{i_expo}.fits')
            psf_ima.write(outpath_tmp)
        logger.info(f'PSF images saved to {psf_dir_tmp}')

    # +++ background noise
    if band in ['u', 'g', 'r', 'i']:
        noise_list = [NoiseModule.GaussianNoise(rms, rng_seed=int(rng_seed_band+120*gal_rotation_angle+94*i_expo)) 
                                                            for gal_rotation_angle in gal_rotation_angles
                                                            for i_expo in range(5)]
    else:
        noise_list = [NoiseModule.GaussianNoise(rms, rng_seed=int(rng_seed_band+120*gal_rotation_angle)) 
                                                            for gal_rotation_angle in gal_rotation_angles]

    # +++ PSF map
    if (False in outpath_PSF_exist_list):
        mag_PSF_2 = 18. # for noise_flux = 2
        mag_PSF = mag_PSF_2 - 2.5*np.log10(rms/2.)
        image_PSF = PSFModule.PSFmap(PSF, pixel_scale, mag_PSF, N_PSF=N_PSF, sep_PSF=sep_PSF, rng_seed=rng_seed_band)

        for i_ima, outpath_PSF_exist in enumerate(outpath_PSF_exist_list):
            if (not outpath_PSF_exist):
                image_PSF_tmp = image_PSF.copy()
                ## noise background
                noise_tmp = noise_list[i_ima]
                image_PSF_tmp.addNoise(noise_tmp)
                ## save
                outpath_PSF_name = outpath_PSF_name_list[i_ima]
                image_PSF_tmp.write(outpath_PSF_name)
                logger.info(f"PSF map saved as {outpath_PSF_name}.")

    # logger.debug(f'Time (everything before sky image): {time.time()-start_time} s')
    # start_time = time.time()

    # +++ sky image
    if (False in outpath_image_exist_list):

        # bounds and wcs from galaxy sky positions
        RA_gals = gals_info_band['RA'] # degree
        DEC_gals = gals_info_band['DEC'] # degree
        RA_min = np.min(RA_gals)
        RA_max = np.max(RA_gals)
        DEC_min = np.min(DEC_gals)
        DEC_max = np.max(DEC_gals)
        bounds, wcs = ObjModule.WCS(RA_min, RA_max, DEC_min, DEC_max, pixel_scale)

        # star image
        if (stars_info_band is not None):
            image_stars = ObjModule.StarsImage(bounds, wcs, band, pixel_scale, PSF, stars_info_band)
        else:
            image_stars = None

        # logger.debug(f'Time (star image, band {band}): {time.time()-start_time} s')
        # start_time = time.time()

        # galaxy images
        gal_rotation_angle0 = None
        # save_weight = True
        for i_ima, outpath_image_exist in enumerate(outpath_image_exist_list):
            
            if (not outpath_image_exist):
                outpath_image_name = outpath_image_name_list[i_ima]

                gal_rotation_angle = float(re.search(r'_rot(\d+)', outpath_image_name).group(1))
                if gal_rotation_angle != gal_rotation_angle0:
                    # if gal_rotation_angle0 is not None:
                    #     save_weight = False
                    image_galaxies0 = ObjModule.GalaxiesImage(bounds, wcs, band, pixel_scale, PSF,
                                                    gals_info_band, gal_rotation_angle=gal_rotation_angle, g_cosmic=g_cosmic)
                    gal_rotation_angle0 = gal_rotation_angle

                image_galaxies = image_galaxies0.copy()

                # logger.debug(f'Time (galaxy image, band {band}, rot {gal_rotation_angle}): {time.time()-start_time} s')
                # start_time = time.time()

                ## add stars
                if (image_stars is not None):
                    image_galaxies += image_stars

                ## add noise background
                noise_tmp = noise_list[i_ima]
                image_galaxies.addNoise(noise_tmp)

                if band in ['u', 'g', 'r', 'i']:

                    id_exposure = int(re.search(r'_expo(\d+)', outpath_image_name).group(1))
                    image_tile, weights_tile = KiDSModule.getKiDStile(image_galaxies, id_exposure=id_exposure, n_exposures=5)
                    
                    ## save the noisy image
                    image_tile.write(outpath_image_name)
                    logger.info(f"Noisy sky image saved as {outpath_image_name}")

                    ## save the weight image
                    outpath_wei_name = outpath_image_name.replace('.fits', '.weight.fits')
                    weights_tile.write(outpath_wei_name)
                    logger.info(f"weight image saved as {outpath_wei_name}")

                    # logger.debug(f'Time (chop image, band {band}, rot {gal_rotation_angle}, id_expo {id_exposure}): {time.time()-start_time} s')
                    # start_time = time.time()

                    # produce chips if required
                    if save_image_chips:

                        image_chips = KiDSModule.getKiDSchips_tile(image_tile)
                        chip_dir_tmp = os.path.join(outpath_dir, f'chips_tile{tile_label}_band{band}_rot{gal_rotation_angle:.0f}')
                        if not os.path.exists(chip_dir_tmp):
                            os.mkdir(chip_dir_tmp)

                        for i_chip, image_chip in enumerate(image_chips):

                            outpath_tmp = os.path.join(chip_dir_tmp, f'expo{id_exposure}_chip{i_chip}.fits')
                            image_chip.write(outpath_tmp)

                        logger.info(f'Image chips saved to {chip_dir_tmp} for exposure {id_exposure}.')

                        # # save chip weights
                        # if save_weight:
                        #     ## weight images
                        #     chip_wei_dir_tmp = os.path.join(outpath_dir, f'chips_tile{tile_label}_band{band}_weight')
                        #     if not os.path.exists(chip_wei_dir_tmp):
                        #         os.mkdir(chip_wei_dir_tmp)
                        #     for i_chip, image_chip in enumerate(image_chips):
                        #         ## associated weight image
                        #         weights_chip = image_chip.copy()
                        #         weights_chip.fill(1./5)
                        #         outpath_tmp = os.path.join(chip_wei_dir_tmp, f'expo{id_exposure}_chip{i_chip}.weight.fits')
                        #         weights_chip.write(outpath_tmp)
                        #     logger.info(f'Weight image chips saved to {chip_wei_dir_tmp} for exposure {id_exposure}.')

                else:
                    ## save the noisy image
                    image_galaxies.write(outpath_image_name)
                    logger.info(f"Noisy sky image saved as {outpath_image_name}")

    return 0

def RunParallel_PSFNoisySkyImages(survey, outpath_dir, rng_seed, mag_zero,
                                            Nmax_proc, 
                                            N_tiles, bands, pixel_scale_list,
                                            noise_info,
                                            gals_info, gal_rotation_angles=[0.], g_cosmic=[0, 0], gal_position_type='true',
                                            stars_area=None, stars_info=None, star_position_type='random',
                                            PSF_map=False, N_PSF=100, sep_PSF=120,
                                            image_chips=None, image_PSF=None):
    '''
    Run ImSim for multi-tile of mutli-band with parallel process.
        Support extending input catalogues 
        Simulate tiles with 1sqdeg
    '''
    logger.info('Running ImSim pipeline...')
    logger.info(f'Number of tiles: {N_tiles}')
    logger.info(f'Bands: {bands}')
    logger.info(f'Pixel scales: {pixel_scale_list}')

    # check if the noise_info is enough for desired N_tiles
    if len(noise_info) < N_tiles:
        raise Exception(f'tiles in noise_info is not enough for desired N_tiles: {len(noise_info)}<{N_tiles}')
    # select tiles from noise_info
    noise_info_selec = noise_info.iloc[:N_tiles]

    # save psf map or not
    if PSF_map:
        outpath_PSF_map_dir = os.path.join(outpath_dir, 'psf_map')
        if not os.path.exists(outpath_PSF_map_dir):
            os.mkdir(outpath_PSF_map_dir)

    # survey type
    if 'simple' in survey: 
        area_tot = int(re.findall(r"\d+", survey)[0])
        area_ra = (area_tot)**0.5
        area_dec = area_ra
        logger.info(f'Simple tile images with {area_tot} square degrees for each tile.')
    elif survey.lower() == 'kids_sameexpo':
        # area taking account of the dither pattern
        area_ra = 1.05 # degree
        area_dec = 1.12 # degree
        logger.info(f'KiDS-like images with 5 exposures and dither pattern.')
        logger.info(f'      Use the same PSF and noise rms for all exposures.')
    else:
        raise Exception(f'Unsupported survey type: {survey}!')

    # total area spanned by the input galaxies
    ## assuming the input catalogue spanning continuously in a rectangular area
    ra_min0 = np.min(gals_info['RA'])
    ra_max0 = np.max(gals_info['RA'])
    dec_min0 = np.min(gals_info['DEC'])
    dec_max0 = np.max(gals_info['DEC'])
    ## 0.999 is 1
    area_ra0 = (ra_max0 - ra_min0) + 0.001
    area_dec0 = (dec_max0 - dec_min0) + 0.001

    # check if the total area is enough for one tile
    if (area_ra0 < area_ra) or (area_dec0 < area_dec):
        raise Exception(f'Input galaxy catalogue is not enough for producing one tile!\n\
        dRA_tile={area_ra}, dRA_cata={area_ra0}; dDEC_tile={area_dec}, dDEC_cata={area_dec0}')

    # number of tiles along RA direction
    N_ra = int(area_ra0/area_ra)
    # number along dec direction 
    N_dec = int(area_dec0/area_dec)
    # check if the total area is enough for chopping
    if (N_ra*N_dec) < N_tiles:
        logger.warning(f'Input galaxy catalogue is not enough for required number of tiles: {N_ra}*{N_dec}<{N_tiles}')
        logger.warning(f'   repeating patterns will be produced')

    # chop input catalogue to tiles 
    gals_info_list = []
    rng_seed_list = []
    i_ra = 0
    i_dec = 0
    if (stars_info is not None):
        # how many stars in each tile
        if (area_ra*area_dec) > stars_area:
            raise Exception(f'Input star catalogue is not enough for required tile area: {area_ra}*{area_dec}<{stars_area} !')
        Nstar_even = int(len(stars_info) / stars_area * (area_ra*area_dec))
        stars_info_list = []
    for tile_label in noise_info_selec['label']:

        # rng seed associated with tile labels
        rng_seed_tile = rng_seed + np.array(re.findall(r"\d+", tile_label), dtype=np.int).sum()
        rng_seed_list.append(rng_seed_tile)

        # sky area for a tile
        ra_min = ra_min0 + area_ra*i_ra
        ra_max = ra_min + area_ra
        dec_min = dec_min0 + area_dec*i_dec
        dec_max = dec_min + area_dec

        # select galaxies
        mask_ra = (gals_info['RA'] >= ra_min) & (gals_info['RA'] < ra_max)
        mask_dec = (gals_info['DEC'] >= dec_min) & (gals_info['DEC'] < dec_max)
        gals_info_selec = gals_info[mask_ra & mask_dec].copy()
        gals_info_selec.reset_index(drop=True, inplace=True)

        ## change position if desired
        if gal_position_type == 'grid':
            logger.info('Galaxies are placed in a grid.')
            Ngal = len(gals_info_selec)
            apart = 18./3600. # apart 18arcsec for each galaxy
            Nrow = int(Ngal**0.5)
            X_gals = np.arange(apart, apart+Nrow*apart, apart)
            Y_gals = np.repeat(X_gals, Nrow)
            X_gals = np.tile(X_gals, Nrow)
            ### check outliers
            Nrow = Ngal - len(X_gals)
            if Nrow > 0:
                X_gals = np.concatenate([X_gals, np.arange(apart, apart+Nrow*apart, apart)])
                Y_gals = np.concatenate([Y_gals, np.full(Nrow, Y_gals[-1]+apart)])
            elif Nrow < 0:
                X_gals = X_gals[:Ngal]
                Y_gals = Y_gals[:Ngal]
            ### over-write
            gals_info_selec.loc[:, 'RA'] = X_gals
            gals_info_selec.loc[:, 'DEC'] = Y_gals

        ## output info
        output_col_tmp = ['index', 'RA', 'DEC'] + bands
        output_tmp = gals_info_selec[output_col_tmp]
        outpath_tmp = os.path.join(outpath_dir, f'gals_info_tile{tile_label}.feather')
        output_tmp.to_feather(outpath_tmp)
        logger.info(f'galaxy info saved to {outpath_tmp}')

        ## magnitude to flux
        for band in bands:
            gals_info_selec.loc[:, band] = 10**(-0.4*(gals_info_selec[band]-mag_zero))

        ## save to list
        gals_info_list.append(gals_info_selec)

        # select stars
        if (stars_info is not None):

            if star_position_type == 'random':
                ## randomly select stars
                np.random.seed(rng_seed_tile)
                random.seed(rng_seed_tile)
                mask_star = random.sample(range(len(stars_info)), Nstar_even)
                stars_info_selec = stars_info.iloc[mask_star].copy()
                stars_info_selec.reset_index(drop=True, inplace=True)
                ## randomly place stars
                ### sky area
                ra_min_true = np.min(gals_info_selec['RA'])
                ra_max_true = np.max(gals_info_selec['RA'])
                dec_min_true = np.min(gals_info_selec['DEC'])
                dec_max_true = np.max(gals_info_selec['DEC'])
                ### sampling
                stars_info_selec.loc[:, 'RA'] = np.random.uniform(low=ra_min_true, high=ra_max_true, size=Nstar_even)
                stars_info_selec.loc[:, 'DEC'] = np.random.uniform(low=dec_min_true, high=dec_max_true, size=Nstar_even)
            elif star_position_type == 'true':
                # use true star location
                mask_ra = (stars_info['RA'] >= ra_min) & (stars_info['RA'] < ra_max)
                mask_dec = (stars_info['DEC'] >= ra_min) & (stars_info['DEC'] < ra_max)
                stars_info_selec = stars_info[mask_ra & mask_dec].copy()
                stars_info_selec.reset_index(drop=True, inplace=True)
            else:
                raise Exception(f'Unsupported star_position_type: {star_position_type} !')
    
            ## output info
            output_col_tmp = ['index', 'RA', 'DEC'] + bands
            output_tmp = stars_info_selec[output_col_tmp]
            outpath_tmp = os.path.join(outpath_dir, f'stars_info_tile{tile_label}.feather')
            output_tmp.to_feather(outpath_tmp)
            logger.info(f'star info saved to {outpath_tmp}')

            ## magnitude to flux
            for band in bands:
                stars_info_selec.loc[:, band] = 10**(-0.4*(stars_info_selec[band]-mag_zero))

            ## save to list
            stars_info_list.append(stars_info_selec)

        # iterating
        i_ra += 1
        if i_ra == N_ra:
            i_ra = 0
            i_dec += 1
            if i_dec == N_dec:
                i_dec = 0

    # collect parameters for workers
    para_lists = []
    for i_tile, noise_info_tile in noise_info_selec.iterrows():

        tile_label = noise_info_tile['label']
        rng_seed_tile = rng_seed_list[i_tile]
        gals_info_tile = gals_info_list[i_tile]
        if (stars_info is not None):
            stars_info_tile = stars_info_list[i_tile]

        for i_band, band in enumerate(bands):

            pixel_scale = pixel_scale_list[i_band]

            # save image chips or not
            if image_chips is not None:
                save_image_chips = image_chips[i_band]

            # save single psf image or not
            if image_PSF is not None:
                save_image_PSF = image_PSF[i_band]

            # noise info
            rms = noise_info_tile[f'rms_{band}']
            seeing = noise_info_tile[f'seeing_{band}']
            beta = noise_info_tile[f'beta_{band}']
            psf_e = [noise_info_tile[f'psf_e1_{band}'], noise_info_tile[f'psf_e2_{band}']]

            rng_seed_band = rng_seed_tile + i_band

            # galaxies
            name = ['RA','DEC',
                'sersic_n','Re','axis_ratios','position_angles',
                'bulge_fractions','bulge_Re','bulge_axis_ratios','bulge_n',
                'disk_Re','disk_axis_ratios', band]
            gals_info_band = gals_info_tile[name]

            # stars
            if (stars_info is not None):
                name = ['RA','DEC', band]
                stars_info_band = stars_info_tile[name]
            else:
                stars_info_band = None

            # save tag
            outpath_image_basename = os.path.join(outpath_dir, f'tile{tile_label}_band{band}')
            if PSF_map:
                outpath_PSF_basename = os.path.join(outpath_PSF_map_dir, f'tile{tile_label}_band{band}')
            else:
                outpath_PSF_basename = None

            # collected parameters
            para_list = (tile_label, band, pixel_scale, rng_seed_band, outpath_image_basename,
                        rms, seeing, beta, psf_e,
                        g_cosmic,
                        gals_info_band, gal_rotation_angles,
                        stars_info_band,
                        outpath_PSF_basename, N_PSF, sep_PSF,
                        save_image_chips, save_image_PSF,
                        outpath_dir)

            para_lists.append(para_list)

    # start parallel running 
    N_tasks = len(para_lists)
    if N_tasks < Nmax_proc:
        N_proc = N_tasks
    else:
        N_proc = Nmax_proc
    logger.info(f'Number of total tasks: {N_tasks} (= {N_tiles}tiles * {len(bands)}bands)')
    logger.info(f'Maximum number of processes: {N_proc}')
    i_worker = 0
    p_list = []
    while True:
        N_running = len(mp.active_children())
        logger.debug(f'Number of running {N_running}')
        if i_worker == N_tasks:
            break
        elif N_running >= N_proc:
            time.sleep(10.)
        else:
            logger.info(f'Start worker {i_worker}')
            if 'simple' in survey:
                p = mp.Process(target=_PSFNoisySkyImages_simple, args=(para_lists[i_worker],))
            elif survey.lower() == 'kids_sameexpo':
                p = mp.Process(target=_PSFNoisySkyImages_KiDS_sameExpo, args=(para_lists[i_worker],))
            else:
                raise Exception(f'Unsupported survey type: {survey}!')

            i_worker += 1
            p.start()
            p_list.append(p)
            time.sleep(1.)
    for p in p_list:
        p.join()

    logger.info('ImSim pipeline finished.')