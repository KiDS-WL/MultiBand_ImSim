# -*- coding: utf-8 -*-
# @Author: lshuns
# @Date:   2020-12-09 19:21:53
# @Last Modified by:   lshuns
# @Last Modified time: 2021-11-11 17:12:19

### running module for ImSim

###### dependence (surveys):
######      ImSimSkySimple: Simple image without any survey strategy
######      ImSimSkyKiDS: KiDS-like images 
from ImSimSkySimple import _PSFNoisySkyImages_simple
from ImSimSkyKiDS import _PSFNoisySkyImages_KiDS_sameExpo, _PSFNoisySkyImages_KiDS_singleExpo, _PSFNoisySkyImages_KiDS_varChips

import os
import re
import math
import time
import random
import logging

import numpy as np
import pandas as pd
import multiprocessing as mp

logger = logging.getLogger(__name__)

def RunParallel_PSFNoisySkyImages(survey, outpath_dir, outcata_dir, rng_seed, mag_zero,
                                            Nmax_proc,
                                            N_tiles, bands, pixel_scale_list, image_type_list,
                                            noise_info,
                                            gals_info, gal_rotation_angles=[0.], g_cosmic=[0, 0], gal_position_type=['true', 18],
                                            stars_area=None, stars_info=None, star_position_type='random',
                                            PSF_map=[], N_PSF=100, sep_PSF=120,
                                            image_chips=None, image_PSF=None,
                                            psf_type='moffat'):
    '''
    Run ImSim for multi-tile of mutli-band with parallel process.
        Support extending input catalogues
        Simulate tiles with 1sqdeg
    '''
    logger.info('Running ImSim pipeline...')
    logger.info(f'Survey: {survey}')
    logger.info(f'Number of tiles: {N_tiles}')
    logger.info(f'Bands: {bands}')
    logger.info(f'Pixel scales: {pixel_scale_list}')
    logger.info(f'Image types: {image_type_list}')
    logger.info(f'PSF profile: {psf_type}')
    logger.info(f'Maximum number of processes: {Nmax_proc}')

    # constant shear or variable
    if 'gamma1' in gals_info.columns:
        g_const = False
        logger.info('Using variable shears')
    else:
        g_const = True
        logger.info('Using a constant shear')

    # check if the noise_info is enough for desired N_tiles
    if len(noise_info) < N_tiles:
        raise Exception(f'tiles in noise_info is not enough for desired N_tiles: {len(noise_info)}<{N_tiles}')
    # select tiles from noise_info
    noise_info_selec = noise_info.iloc[:N_tiles]
    del noise_info

    # save psf map or not
    if True in PSF_map:
        outpath_PSF_map_dir = os.path.join(outpath_dir, 'psf_map')
        if not os.path.exists(outpath_PSF_map_dir):
            os.mkdir(outpath_PSF_map_dir)

    # survey type
    if 'simple' in survey:
        numeric_const_pattern = r"[-+]?(?:(?:\d*\.\d+)|(?:\d+\.?))(?:[Ee][+-]?\d+)?"
        area_tot = float(re.findall(numeric_const_pattern, survey)[0])
        # area for each tile
        area_ra = (area_tot)**0.5
        area_dec = area_ra
    elif survey.lower() == 'kids':
        area_ra = 1. # degree
        area_dec = 1. # degree
        # check if grid is required
        if (gal_position_type[0]=='grid'):
            raise Exception('KiDS survey does not observe grid world :( Please use other survey.')
    elif survey.lower() == 'one_tile':
        # area same as the input
        ra_min0 = np.min(gals_info['RA'])
        ra_max0 = np.max(gals_info['RA'])
        dec_min0 = np.min(gals_info['DEC'])
        dec_max0 = np.max(gals_info['DEC'])
        ## 0.999 is 1
        area_ra = (ra_max0 - ra_min0) + 0.0005
        area_dec = (dec_max0 - dec_min0) + 0.0005
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
        logger.warning(f'Input galaxy catalogue is not enough for producing one tile.\n\
        Use the whole input catalogue.\n\
        dRA_tile={area_ra}, dRA_cata={area_ra0}; dDEC_tile={area_dec}, dDEC_cata={area_dec0}')
        area_ra = area_ra0
        area_dec = area_dec0

    # number of tiles along RA direction
    N_ra = int(area_ra0/area_ra)
    # number along dec direction
    N_dec = int(area_dec0/area_dec)
    # check if the total area is enough for chopping
    if (N_ra*N_dec) < N_tiles:
        logger.warning(f'Input galaxy catalogue is not enough for required number of tiles: {N_ra}*{N_dec}<{N_tiles}')
        logger.warning(f'   repeating patterns will be produced')

    # chop input catalogue to tiles
    if (stars_info is not None):
        # how many stars in each tile
        if (area_ra*area_dec) > stars_area:
            raise Exception(f'Input star catalogue is not enough for required tile area: {area_ra}*{area_dec}<{stars_area} !')
        Nstar_even = int(len(stars_info) / stars_area * (area_ra*area_dec))
        stars_info_list = []
    gals_info_list = []
    rng_seed_list = []
    i_ra = 0
    i_dec = 0
    for i_tile, noise_info_tile in noise_info_selec.iterrows():
        tile_label = noise_info_tile['label']

        if (i_tile!=0) and (i_ra==0) and (i_dec==0):
            logger.warning(f'repeating patterns started from tile {tile_label}')

        # rng seed associated with tile labels
        rng_seed_tile = rng_seed + np.array(re.findall(r"\d+", tile_label), dtype=np.int).sum()*1000000
        rng_seed_list.append(rng_seed_tile)

        ## output noise info
        outpath_tmp = os.path.join(outcata_dir, f'noise_info_tile{tile_label}.csv')
        pd.DataFrame([noise_info_tile]).to_csv(outpath_tmp, index=False)
        logger.info(f'noise info saved to {outpath_tmp}')

        # sky area for a tile
        ra_min = ra_min0 + area_ra * i_ra
        ra_max = ra_min + area_ra
        dec_min = dec_min0 + area_dec * i_dec
        dec_max = dec_min + area_dec

        # select galaxies
        mask_ra = (gals_info['RA'] >= ra_min) & (gals_info['RA'] < ra_max)
        mask_dec = (gals_info['DEC'] >= dec_min) & (gals_info['DEC'] < dec_max)
        gals_info_selec = gals_info[mask_ra & mask_dec].copy()
        gals_info_selec.reset_index(drop=True, inplace=True)

        ## change position if desired
        if gal_position_type[0] == 'grid':
            logger.info('Galaxies are placed in a grid.')
            Ngal = len(gals_info_selec)
            apart = gal_position_type[1]/3600. # apart for each galaxy
            # how many row and column needed to place all galaxies
            N_rows = math.ceil(Ngal**0.5)
            N_cols = math.ceil(Ngal/N_rows)
            X_gals = np.arange(0, apart+(N_cols-1)*apart, apart)
            Y_gals = np.arange(0, apart+(N_rows-1)*apart, apart)
            X_gals, Y_gals = np.meshgrid(X_gals, Y_gals)
            X_gals = X_gals.flatten()
            Y_gals = Y_gals.flatten()
            # random order to avoid systematic patterns
            random.seed(rng_seed_tile)
            index_selected = random.sample(range(len(X_gals)), Ngal)
            ### over-write
            gals_info_selec.loc[:, 'RA'] = X_gals[index_selected]
            gals_info_selec.loc[:, 'DEC'] = Y_gals[index_selected]
            ### delete useless variables
            del index_selected
            del X_gals
            del Y_gals
        elif gal_position_type[0] == 'random':
            logger.info('Galaxies are placed randomly.')
            np.random.seed(rng_seed_tile)
            Ngal = len(gals_info_selec)
            ### random sampling
            gals_info_selec.loc[:, 'RA'] = np.random.uniform(low=ra_min, high=ra_max, size=Ngal)
            gals_info_selec.loc[:, 'DEC'] = np.random.uniform(low=dec_min, high=dec_max, size=Ngal)
        else:
            if gal_position_type[0] != 'true':
                raise Exception(f'Unsupported gal_position_type: {gal_position_type[0]} !')

        ## output galaxies info
        output_col_tmp = ['index', 'RA', 'DEC', 'redshift', 'Re', 'axis_ratio', 'position_angle', 'sersic_n'] + bands
        if not g_const:
            output_col_tmp += ['gamma1', 'gamma2']
        output_tmp = gals_info_selec[output_col_tmp].copy()
        ### better naming
        output_tmp = output_tmp.add_suffix(f'_input')
        ### add ellipticity based on q and beta
        g_tmp = (1-output_tmp['axis_ratio_input'])/(1+output_tmp['axis_ratio_input'])
        for gal_rotation_angle in gal_rotation_angles:
            true_pa_tmp = output_tmp['position_angle_input'] + gal_rotation_angle
            output_tmp.loc[:, f'e1_input_rot{int(gal_rotation_angle)}'] = g_tmp * np.cos(2. * (true_pa_tmp/180.*np.pi))
            output_tmp.loc[:, f'e2_input_rot{int(gal_rotation_angle)}'] = g_tmp * np.sin(2. * (true_pa_tmp/180.*np.pi))
        del g_tmp
        del true_pa_tmp
        outpath_tmp = os.path.join(outcata_dir, f'gals_info_tile{tile_label}.feather')
        output_tmp.to_feather(outpath_tmp)
        del output_tmp
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
                stars_info_selec.loc[:, 'RA'] = np.random.uniform(low=ra_min, high=ra_max, size=Nstar_even)
                stars_info_selec.loc[:, 'DEC'] = np.random.uniform(low=dec_min, high=dec_max, size=Nstar_even)
            elif star_position_type == 'true':
                # use true star location
                mask_ra = (stars_info['RA'] >= ra_min) & (stars_info['RA'] < ra_max)
                mask_dec = (stars_info['DEC'] >= ra_min) & (stars_info['DEC'] < ra_max)
                stars_info_selec = stars_info[mask_ra & mask_dec].copy()
                stars_info_selec.reset_index(drop=True, inplace=True)
            else:
                raise Exception(f'Unsupported star_position_type: {star_position_type} !')

            ## output star info
            output_col_tmp = ['index', 'RA', 'DEC'] + bands
            output_tmp = stars_info_selec[output_col_tmp].copy()
            ### better naming
            output_tmp = output_tmp.add_suffix(f'_input')
            outpath_tmp = os.path.join(outcata_dir, f'stars_info_tile{tile_label}.feather')
            output_tmp.to_feather(outpath_tmp)
            del output_tmp
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

    # release some space
    del gals_info

    # collect parameters for workers
    para_lists = []
    image_type_labels = []
    for i_tile, noise_info_tile in noise_info_selec.iterrows():

        tile_label = noise_info_tile['label']
        rng_seed_tile = rng_seed_list[i_tile]
        gals_info_tile = gals_info_list[i_tile]
        if (stars_info is not None):
            stars_info_tile = stars_info_list[i_tile]

        for i_band, band in enumerate(bands):

            rng_seed_band = rng_seed_tile + i_band

            pixel_scale = pixel_scale_list[i_band]
            image_type = image_type_list[i_band]

            # save image chips or not
            if image_chips is not None:
                save_image_chips = image_chips[i_band]

                # make a directory if save
                if save_image_chips:
                    for gal_rotation_angle in gal_rotation_angles:
                        chip_dir_tmp = os.path.join(outpath_dir, f'chips_tile{tile_label}_band{band}_rot{gal_rotation_angle:.0f}')
                        if (not os.path.exists(chip_dir_tmp)):
                            os.mkdir(chip_dir_tmp)
            else:
                save_image_chips = False

            # save single psf image or not
            if image_PSF is not None:
                save_image_PSF = image_PSF[0][i_band]
                image_PSF_size = image_PSF[1]

                # make a directory if save
                if save_image_PSF:
                    psf_dir_tmp = os.path.join(outpath_dir, f'psf_tile{tile_label}_band{band}')
                    if (not os.path.exists(psf_dir_tmp)):
                        os.mkdir(psf_dir_tmp)
            else:
                save_image_PSF = False

            # save tag
            outpath_image_basename = os.path.join(outpath_dir, f'tile{tile_label}_band{band}')
            if (PSF_map) and (PSF_map[i_band]):
                outpath_PSF_basename = os.path.join(outpath_PSF_map_dir, f'tile{tile_label}_band{band}')
            else:
                outpath_PSF_basename = None

            # galaxies
            name = ['RA','DEC',
                'sersic_n','Re','axis_ratio','position_angle',
                'bulge_fraction','bulge_Re','bulge_axis_ratio','bulge_n',
                'disk_Re','disk_axis_ratio', band]
            if not g_const:
                name += ['gamma1', 'gamma2']
            gals_info_band = gals_info_tile[name]

            # stars
            if (stars_info is not None):
                name = ['RA','DEC', band]
                stars_info_band = stars_info_tile[name]
            else:
                stars_info_band = None

            # noise & psf
            if image_type.lower() == 'varchips':

                # always save chips for varChips
                for gal_rotation_angle in gal_rotation_angles:
                    chip_dir_tmp = os.path.join(outpath_dir, f'chips_tile{tile_label}_band{band}_rot{gal_rotation_angle:.0f}')
                    if (not os.path.exists(chip_dir_tmp)):
                        os.mkdir(chip_dir_tmp)

                if survey.lower() == 'kids':
                    if band != 'r':
                        raise Exception(f'{band} does not support varChips for KiDS survey (only use r-band for KiDS)!')

                    # get noise info for each exposure
                    for i_expo in range(5):
                        # noise is still same for varChips
                        rms = noise_info_tile[f'rms_{band}_expo{i_expo}']
                        # psf is different
                        ## psf profiles
                        if psf_type.lower() == 'moffat':
                            seeing_chips = [noise_info_tile[f'seeing_{band}_expo{i_expo}_chip{i_chip}'] for i_chip in range(32)]
                            beta_chips = [noise_info_tile[f'beta_{band}_expo{i_expo}_chip{i_chip}'] for i_chip in range(32)]
                            psf_info_chips = [psf_type, seeing_chips, beta_chips]
                            del seeing_chips, beta_chips
                        elif psf_type.lower() == 'airy':
                            lam_chips = [noise_info_tile[f'lam_{band}_expo{i_expo}_chip{i_chip}'] for i_chip in range(32)]
                            diam_chips = [noise_info_tile[f'diam_{band}_expo{i_expo}_chip{i_chip}'] for i_chip in range(32)]
                            obscuration_chips = [noise_info_tile[f'obscuration_{band}_expo{i_expo}_chip{i_chip}'] for i_chip in range(32)]
                            psf_info_chips = [psf_type, lam_chips, diam_chips, obscuration_chips]
                            del lam_chips, diam_chips, obscuration_chips
                        ## psf e
                        psf_e_chips = [[noise_info_tile[f'psf_e1_{band}_expo{i_expo}_chip{i_chip}'] for i_chip in range(32)], 
                                            [noise_info_tile[f'psf_e2_{band}_expo{i_expo}_chip{i_chip}'] for i_chip in range(32)]]
                        psf_info_chips.append(psf_e_chips)
                        del psf_e_chips

                        # collected parameters
                        ## rotations simulated separately
                        for gal_rotation_angle in gal_rotation_angles:
                            para_list = (tile_label, band, pixel_scale, rng_seed_band,
                                        rms, psf_info_chips,
                                        g_cosmic,
                                        gals_info_band, gal_rotation_angle,
                                        stars_info_band,
                                        outpath_PSF_basename, N_PSF, sep_PSF,
                                        save_image_PSF, image_PSF_size, 
                                        outpath_dir,
                                        i_expo,
                                        gal_position_type,
                                        g_const)
                            para_lists.append(para_list)
                            # label
                            image_type_labels.append(image_type)
                else:
                    raise Exception(f'Survey type {survey} does not support varChips images!')

            elif image_type.lower() == 'diffexpo':

                # always save chips for diffexpo
                for gal_rotation_angle in gal_rotation_angles:
                    chip_dir_tmp = os.path.join(outpath_dir, f'chips_tile{tile_label}_band{band}_rot{gal_rotation_angle:.0f}')
                    if (not os.path.exists(chip_dir_tmp)):
                        os.mkdir(chip_dir_tmp)

                if survey.lower() == 'kids':
                    # number of exposures
                    if band == 'u':
                        n_exposures = 4
                    elif band in ['g', 'r', 'i']:
                        n_exposures = 5
                    else:
                        raise Exception(f'{band} is not in OmegaCAM! Cannot use KiDS_diffExpo!')

                    # get noise info for each exposure
                    for i_expo in range(n_exposures):
                        # noise info
                        rms = noise_info_tile[f'rms_{band}_expo{i_expo}']

                        # psf profiles
                        if psf_type.lower() == 'moffat':
                            seeing = noise_info_tile[f'seeing_{band}_expo{i_expo}']
                            beta = noise_info_tile[f'beta_{band}_expo{i_expo}']
                            psf_info = [psf_type, seeing, beta]
                            del seeing, beta
                        elif psf_type.lower() == 'airy':
                            lam = noise_info_tile[f'lam_{band}_expo{i_expo}']
                            diam = noise_info_tile[f'diam_{band}_expo{i_expo}']
                            obscuration = noise_info_tile[f'obscuration_{band}_expo{i_expo}']
                            psf_info = [psf_type, lam, diam, obscuration]
                            del lam, diam, obscuration

                        # psf e
                        psf_e = [noise_info_tile[f'psf_e1_{band}_expo{i_expo}'], noise_info_tile[f'psf_e2_{band}_expo{i_expo}']]
                        psf_info.append(psf_e)
                        del psf_e

                        # collected parameters
                        ## rotations simulated separately
                        for gal_rotation_angle in gal_rotation_angles:
                            para_list = (tile_label, band, pixel_scale, rng_seed_band,
                                        rms, psf_info,
                                        g_cosmic,
                                        gals_info_band, gal_rotation_angle,
                                        stars_info_band,
                                        outpath_PSF_basename, N_PSF, sep_PSF,
                                        save_image_PSF, image_PSF_size, 
                                        outpath_dir,
                                        i_expo,
                                        gal_position_type,
                                        g_const)
                            para_lists.append(para_list)
                            # label
                            image_type_labels.append(image_type)
                else:
                    raise Exception(f'Survey type {survey} does not support diffExpo images!')

            else:
                # noise info
                rms = noise_info_tile[f'rms_{band}']

                # psf profiles
                if psf_type.lower() == 'moffat': 
                    seeing = noise_info_tile[f'seeing_{band}']
                    beta = noise_info_tile[f'beta_{band}']
                    psf_info = [psf_type, seeing, beta]
                    del seeing, beta
                elif psf_type.lower() == 'airy':
                    lam = noise_info_tile[f'lam_{band}']
                    diam = noise_info_tile[f'diam_{band}']
                    obscuration = noise_info_tile[f'obscuration_{band}']
                    psf_info = [psf_type, lam, diam, obscuration]
                    del lam, diam, obscuration

                # psf e
                psf_e = [noise_info_tile[f'psf_e1_{band}'], noise_info_tile[f'psf_e2_{band}']]
                psf_info.append(psf_e)
                del psf_e

                # collected parameters
                ## rotations simulated separately
                for gal_rotation_angle in gal_rotation_angles:
                    para_list = (tile_label, band, pixel_scale, rng_seed_band, outpath_image_basename,
                                rms, psf_info,
                                g_cosmic,
                                gals_info_band, gal_rotation_angle,
                                stars_info_band,
                                outpath_PSF_basename, N_PSF, sep_PSF,
                                save_image_chips, save_image_PSF, image_PSF_size,
                                outpath_dir,
                                gal_position_type,
                                g_const)
                    para_lists.append(para_list)
                    # label
                    image_type_labels.append(image_type)

    # release some space
    del noise_info_selec, rng_seed_list
    del gals_info_list, gals_info_tile, gals_info_band
    if (stars_info is not None):
        del stars_info_list, stars_info_tile, stars_info_band

    # start parallel running
    N_tasks = len(para_lists)
    logger.info(f'Total number of tasks: {N_tasks}')
    i_worker = 0
    p_list = []
    while True:
        N_running = len(mp.active_children())
        logger.debug(f'Number of running {N_running}')
        if i_worker == N_tasks:
            break
        elif N_running >= (Nmax_proc-1):
            time.sleep(5.)
        else:
            logger.debug(f'Start worker {i_worker}')
            image_type = image_type_labels[i_worker].lower()
            if image_type == 'varchips':
                if survey.lower() == 'kids':
                    p = mp.Process(target=_PSFNoisySkyImages_KiDS_varChips, args=(para_lists[i_worker],))                
            elif image_type == 'diffexpo':
                if survey.lower() == 'kids':
                    p = mp.Process(target=_PSFNoisySkyImages_KiDS_singleExpo, args=(para_lists[i_worker],))
            elif image_type == 'sameexpo':
                if survey.lower() == 'kids':
                    p = mp.Process(target=_PSFNoisySkyImages_KiDS_sameExpo, args=(para_lists[i_worker],))
            elif image_type == 'simple':
                p = mp.Process(target=_PSFNoisySkyImages_simple, args=(para_lists[i_worker],))
            else:
                raise Exception(f'Unsupported image type: {image_type} for survey {survey}!')

            i_worker += 1
            p.start()
            p_list.append(p)
            time.sleep(1.)
    for p in p_list:
        p.join()

    logger.info('ImSim pipeline finished.')
