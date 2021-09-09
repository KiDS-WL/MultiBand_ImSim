# -*- coding: utf-8 -*-
# @Author: lshuns
# @Date:   2021-07-22 13:25:05
# @Last Modified by:   lshuns
# @Last Modified time: 2021-09-07 14:39:17

### Everything about KiDS-like images
__all__ = ['_PSFNoisySkyImages_KiDS_sameExpo', '_PSFNoisySkyImages_KiDS_singleExpo', '_PSFNoisySkyImages_KiDS_varChips']

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
import glob
import shutil
import galsim
import logging

import numpy as np
import pandas as pd

from astropy.io import fits

logger = logging.getLogger(__name__)

def _PSFNoisySkyImages_KiDS_sameExpo(para_list):
    '''
    Sky image with Gaussian noise and PSF
        KiDS images with dither and gaps
        Same noise and psf for all exposures
        Adjoint to RunParallel_PSFNoisySkyImages in ImSim.py
    '''

    (tile_label, band, pixel_scale, rng_seed_band, outpath_image_basename,
        rms, seeing, beta, psf_e,
        g_cosmic,
        gals_info_band, gal_rotation_angle,
        stars_info_band,
        outpath_PSF_basename, N_PSF, sep_PSF,
        save_image_chips, save_image_PSF, image_PSF_size,
        outpath_dir,
        gal_position_type,
        g_const) = para_list

    logger.info(f'Simulating KiDS_sameExpo image for tile {tile_label} band {band} rot {gal_rotation_angle}...')

    # number of exposures
    if band == 'u':
        n_exposures = 4
    elif band in ['g', 'r', 'i']:
        n_exposures = 5
    else:
        raise Exception(f'{band} is not in OmegaCAM! Cannot use KiDS_sameExpo!')

    # outpath
    outpath_image_name_list = [outpath_image_basename + f'_rot{gal_rotation_angle:.0f}_expo{i_expo}.fits'
                                                                for i_expo in range(n_exposures)]

    if (outpath_PSF_basename is not None):
        outpath_PSF_name_list = [outpath_PSF_basename + f'_rot{gal_rotation_angle:.0f}_expo{i_expo}.fits'
                                                                for i_expo in range(n_exposures)]

    # first check if already exist
    outpath_image_exist_list = np.zeros_like(outpath_image_name_list, dtype=bool)
    for i_name, outpath_image_name in enumerate(outpath_image_name_list):
        try:
            with fits.open(outpath_image_name) as hdul:
                head_tmp = hdul[0].header
            flag_sim = head_tmp['flag_sim']
            if flag_sim >= 1:
                outpath_image_exist_list[i_name] = True
                logger.info(f"{outpath_image_name} already exist.")
                continue
        except FileNotFoundError:
            continue
        except (KeyError, OSError) as e:
            pass
        os.remove(outpath_image_name)

    ## psf map
    if (outpath_PSF_basename is not None) :
        outpath_PSF_exist_list = np.zeros_like(outpath_PSF_name_list, dtype=bool)
        for i_name, outpath_PSF_name in enumerate(outpath_PSF_name_list):
            try:
                with fits.open(outpath_PSF_name) as hdul:
                    head_tmp = hdul[0].header
                flag_sim = head_tmp['flag_sim']
                if flag_sim >= 1:
                    outpath_PSF_exist_list[i_name] = True
                    logger.info(f"{outpath_PSF_name} already exist.")
                    continue
            except FileNotFoundError:
                continue
            except (KeyError, OSError) as e:
                pass
            os.remove(outpath_PSF_name)
    else:
        outpath_PSF_exist_list = [True]

    ## psf image
    ### different rotation has same psf, so only make once
    if (save_image_PSF) and (gal_rotation_angle==0.):
        psf_dir_tmp = os.path.join(outpath_dir, f'psf_tile{tile_label}_band{band}')
        n_files = len(glob.glob(os.path.join(psf_dir_tmp, f'expo*.fits')))
        if n_files == n_exposures:
            logger.info('PSF images already exist.')
        else:
            if os.path.exists(psf_dir_tmp):
                shutil.rmtree(psf_dir_tmp)
            os.mkdir(psf_dir_tmp)

            PSF = PSFModule.MoffatPSF(seeing, moffat_beta=beta, psf_e=psf_e)
            psf_ima = PSFModule.PSFima(PSF, pixel_scale, size=image_PSF_size)

            for i_expo in range(n_exposures):
                outpath_tmp = os.path.join(psf_dir_tmp, f'expo{i_expo}.fits')
                psf_ima.write(outpath_tmp)
            logger.info(f'PSF images saved to {psf_dir_tmp}')

    ## chips
    if save_image_chips:
        for i_ima, outpath_image_exist in enumerate(outpath_image_exist_list):
            if (outpath_image_exist):
                outpath_image_name = outpath_image_name_list[i_ima]
                id_exposure = int(re.search(r'_expo(\d+)', outpath_image_name).group(1))

                chip_dir_tmp = os.path.join(outpath_dir, f'chips_tile{tile_label}_band{band}_rot{gal_rotation_angle:.0f}')
                n_files = len(glob.glob(os.path.join(chip_dir_tmp, f'expo{id_exposure}_chip*.fits')))
                if n_files == 32:
                    logger.info(f'chips already exist for rot{gal_rotation_angle:.0f} expo{id_exposure}.')
                else:
                    image_tile = galsim.fits.read(outpath_image_name)
                    image_chips = KiDSModule.cutKiDSchips(image_tile)

                    for i_chip, image_chip in enumerate(image_chips):
                        outpath_tmp = os.path.join(chip_dir_tmp, f'expo{id_exposure}_chip{i_chip}.fits')
                        image_chip.write(outpath_tmp)
                    logger.info(f'Image chips saved to {chip_dir_tmp} for expo{id_exposure}.')

    ## if all exist, quit
    if (not False in outpath_image_exist_list) and (not False in outpath_PSF_exist_list):
        logger.info("All desired images exist, end the process.")
        return 1

    # +++ PSF
    PSF = PSFModule.MoffatPSF(seeing, moffat_beta=beta, psf_e=psf_e)

    # +++ background noise
    noise_list = [NoiseModule.GaussianNoise(rms, rng_seed=int(rng_seed_band+120*gal_rotation_angle+94*i_expo))
                                                        for i_expo in range(n_exposures)]

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

                ## mark success to the header
                with fits.open(outpath_PSF_name, mode='update') as hdul:
                    head_tmp = hdul[0].header
                    ## update info
                    head_tmp['flag_sim'] = 1

    # +++ sky image
    if (False in outpath_image_exist_list):

        # simple canvas based on the galaxy sky positions
        RA_gals = gals_info_band['RA'] # degree
        DEC_gals = gals_info_band['DEC'] # degree
        RA_min = np.min(RA_gals)
        RA_max = np.max(RA_gals)
        DEC_min = np.min(DEC_gals)
        DEC_max = np.max(DEC_gals)
        canvas = ObjModule.SimpleCanvas(RA_min, RA_max, DEC_min, DEC_max, pixel_scale)
        del RA_gals, DEC_gals, RA_min, RA_max, DEC_min, DEC_max

        # star image
        if (stars_info_band is not None):
            image_stars = ObjModule.StarsImage(canvas, band, pixel_scale, PSF, stars_info_band)
        else:
            image_stars = None

        # galaxy images
        image_galaxies0 = ObjModule.GalaxiesImage(canvas, band, pixel_scale, PSF,
                                gals_info_band, gal_rotation_angle=gal_rotation_angle, g_cosmic=g_cosmic, gal_position_type=gal_position_type,
                                g_const=g_const)

        # cut to exposures
        for i_ima, outpath_image_exist in enumerate(outpath_image_exist_list):

            if (not outpath_image_exist):
                outpath_image_name = outpath_image_name_list[i_ima]

                image_galaxies = image_galaxies0.copy()
                ## add stars
                if (image_stars is not None):
                    image_galaxies += image_stars

                ## noise background
                noise_tmp = noise_list[i_ima]

                ## cut to kids tile
                id_exposure = int(re.search(r'_expo(\d+)', outpath_image_name).group(1))
                image_tile, weights_tile = KiDSModule.cutKiDStile(image_galaxies, noise_tmp, id_exposure=id_exposure)

                ## save the noisy image
                image_tile.write(outpath_image_name)
                logger.info(f"KiDS-like exposure saved as {outpath_image_name}")

                ## save the weight image
                outpath_wei_name = outpath_image_name.replace('.fits', '.weight.fits')
                weights_tile.write(outpath_wei_name)
                logger.info(f"weight image saved as {outpath_wei_name}")

                ## mark success to the header
                with fits.open(outpath_image_name, mode='update') as hdul:
                    head_tmp = hdul[0].header
                    ## update info
                    head_tmp['flag_sim'] = 1

                # produce chips if required
                if save_image_chips:
                    image_chips = KiDSModule.cutKiDSchips(image_tile)
                    chip_dir_tmp = os.path.join(outpath_dir, f'chips_tile{tile_label}_band{band}_rot{gal_rotation_angle:.0f}')
                    for i_chip, image_chip in enumerate(image_chips):
                        outpath_tmp = os.path.join(chip_dir_tmp, f'expo{id_exposure}_chip{i_chip}.fits')
                        image_chip.write(outpath_tmp)
                    logger.info(f'Image chips saved to {chip_dir_tmp} for exposure {id_exposure}.')

    logger.info(f'Finished for tile {tile_label} band {band} rot {gal_rotation_angle}...')
    return 0

def _PSFNoisySkyImages_KiDS_singleExpo(para_list):
    '''
    Sky image with Gaussian noise and PSF
        KiDS images with dither and gaps
        Only simulate one exposure
        Used for different noise and psf for different exposures
        Adjoint to RunParallel_PSFNoisySkyImages in ImSim.py
    '''
    (tile_label, band, pixel_scale, rng_seed_band,
        rms, seeing, beta, psf_e,
        g_cosmic,
        gals_info_band, gal_rotation_angle,
        stars_info_band,
        outpath_PSF_basename, N_PSF, sep_PSF,
        save_image_PSF, image_PSF_size,
        outpath_dir,
        i_expo,
        gal_position_type,
        g_const) = para_list

    logger.info(f'Simulating KiDS exposure for tile {tile_label} band {band} expo {i_expo} rot {gal_rotation_angle}...')

    # outpath
    outpath_image_name_list = [os.path.join(outpath_dir, f'chips_tile{tile_label}_band{band}_rot{gal_rotation_angle:.0f}', f'expo{i_expo}_chip{i_chip}.fits')
                                for i_chip in range(32)]

    if (outpath_PSF_basename is not None):
        outpath_PSF_name = outpath_PSF_basename + f'_rot{gal_rotation_angle:.0f}_expo{i_expo}.fits'

    # first check if already exist
    outpath_image_exist_list = np.zeros_like(outpath_image_name_list, dtype=bool)
    for i_name, outpath_image_name in enumerate(outpath_image_name_list):
        try:
            with fits.open(outpath_image_name) as hdul:
                head_tmp = hdul[0].header
            flag_sim = head_tmp['flag_sim']
            if flag_sim >= 1:
                outpath_image_exist_list[i_name] = True
                logger.info(f"{outpath_image_name} already exist.")
                continue
        except FileNotFoundError:
            continue
        except (KeyError, OSError) as e:
            pass
        os.remove(outpath_image_name)

    ## psf map
    outpath_PSF_exist = False
    if (outpath_PSF_basename is not None):
        try:
            with fits.open(outpath_PSF_name) as hdul:
                head_tmp = hdul[0].header
            flag_sim = head_tmp['flag_sim']
            if flag_sim >= 1:
                outpath_PSF_exist = True
                logger.info(f"{outpath_PSF_name} already exist.")
        except FileNotFoundError:
            pass
        except (KeyError, OSError) as e:
            os.remove(outpath_PSF_name)
    else:
        outpath_PSF_exist = True

    ## psf image
    ### different rotation has same psf, so only make once
    if (save_image_PSF) and (gal_rotation_angle==0.):
        psf_dir_tmp = os.path.join(outpath_dir, f'psf_tile{tile_label}_band{band}')
        psf_ima_file_tmp = os.path.join(psf_dir_tmp, f'expo{i_expo}.fits')
        if os.path.isfile(psf_ima_file_tmp):
            logger.info('PSF image already exist.')
        else:
            PSF = PSFModule.MoffatPSF(seeing, moffat_beta=beta, psf_e=psf_e)
            psf_ima = PSFModule.PSFima(PSF, pixel_scale, size=image_PSF_size)
            psf_ima.write(psf_ima_file_tmp)
            logger.info(f'PSF image saved as {psf_ima_file_tmp}')

    ## if all exist, quit
    if (not False in outpath_image_exist_list) and (outpath_PSF_exist):
        logger.info("All desired images exist, end the process.")
        return 1

    # +++ PSF
    PSF = PSFModule.MoffatPSF(seeing, moffat_beta=beta, psf_e=psf_e)

    # +++ background noise
    noise = NoiseModule.GaussianNoise(rms, rng_seed=int(rng_seed_band+120*gal_rotation_angle+94*i_expo))

    # +++ PSF map
    if (not outpath_PSF_exist):
        mag_PSF_2 = 18. # for noise_flux = 2
        mag_PSF = mag_PSF_2 - 2.5*np.log10(rms/2.)
        image_PSF = PSFModule.PSFmap(PSF, pixel_scale, mag_PSF, N_PSF=N_PSF, sep_PSF=sep_PSF, rng_seed=rng_seed_band)

        ## noise background
        image_PSF.addNoise(noise)

        ## save
        image_PSF.write(outpath_PSF_name)
        logger.info(f"PSF map saved as {outpath_PSF_name}.")

        ## mark success to the header
        with fits.open(outpath_PSF_name, mode='update') as hdul:
            head_tmp = hdul[0].header
            ## update info
            head_tmp['flag_sim'] = 1

    # +++ sky image
    if (False in outpath_image_exist_list):

        # a list of canvas based on galaxy sky positions
        RA_gals = gals_info_band['RA'] # degree
        DEC_gals = gals_info_band['DEC'] # degree
        RA0 = (np.max(RA_gals) + np.min(RA_gals))/2.
        DEC0 = (np.max(DEC_gals) + np.min(DEC_gals))/2.
        canvases_list = KiDSModule.getKiDScanvases(RA0, DEC0, id_exposure=i_expo)
        del RA0, DEC0

        # all desired images
        for i_ima, outpath_image_exist in enumerate(outpath_image_exist_list):

            if (not outpath_image_exist):
                outpath_image_name = outpath_image_name_list[i_ima]

                # chip id
                i_chip = int(re.search(r'_chip(\d+)', outpath_image_name).group(1))
                ## get the canvas accordingly
                canvas = canvases_list[i_chip]

                # galaxy image
                image_galaxies = ObjModule.GalaxiesImage(canvas, band, pixel_scale, PSF,
                                                gals_info_band, gal_rotation_angle=gal_rotation_angle, g_cosmic=g_cosmic, gal_position_type=gal_position_type,
                                                g_const=g_const)

                ## add stars
                if (stars_info_band is not None):
                    image_stars = ObjModule.StarsImage(canvas, band, pixel_scale, PSF, stars_info_band)
                    image_galaxies += image_stars
                    del image_stars

                ## add noise background
                image_galaxies.addNoise(noise)

                ## save the noisy image
                image_galaxies.write(outpath_image_name)
                logger.info(f"KiDS-like chip saved as {outpath_image_name}")

                ## mark success to the header
                with fits.open(outpath_image_name, mode='update') as hdul:
                    head_tmp = hdul[0].header
                    ## update info
                    head_tmp['flag_sim'] = 1

    logger.info(f'Finished for tile {tile_label} band {band} expo {i_expo} rot {gal_rotation_angle}...')
    return 0

def _PSFNoisySkyImages_KiDS_varChips(para_list):
    '''
    Sky image with Gaussian noise and PSF
        KiDS images with dither and gaps
        Only simulate one exposure
        different chips use different psfs
        Adjoint to RunParallel_PSFNoisySkyImages in ImSim.py
    '''
    (tile_label, band, pixel_scale, rng_seed_band,
        rms, seeing_chips, beta_chips, psf_e_chips,
        g_cosmic,
        gals_info_band, gal_rotation_angle,
        stars_info_band,
        outpath_PSF_basename, N_PSF, sep_PSF,
        save_image_PSF, image_PSF_size,
        outpath_dir,
        i_expo,
        gal_position_type,
        g_const) = para_list
    logger.info(f'Simulating KiDS exposure with varChips for tile {tile_label} band {band} expo {i_expo} rot {gal_rotation_angle}...')

    # outpath
    outpath_image_name_list = [os.path.join(outpath_dir, f'chips_tile{tile_label}_band{band}_rot{gal_rotation_angle:.0f}', f'expo{i_expo}_chip{i_chip}.fits')
                                for i_chip in range(32)]

    if (outpath_PSF_basename is not None):
        raise Exception('varChips mode does not support PSF map!')

    # first check if already exist
    outpath_image_exist_list = np.zeros_like(outpath_image_name_list, dtype=bool)
    for i_name, outpath_image_name in enumerate(outpath_image_name_list):
        try:
            with fits.open(outpath_image_name) as hdul:
                head_tmp = hdul[0].header
            flag_sim = head_tmp['flag_sim']
            if flag_sim >= 1:
                outpath_image_exist_list[i_name] = True
                logger.info(f"{outpath_image_name} already exist.")
                continue
        except FileNotFoundError:
            continue
        except (KeyError, OSError) as e:
            pass
        os.remove(outpath_image_name)

    ## psf image
    ### different rotation has same psf, so only make once
    if (save_image_PSF) and (gal_rotation_angle==0.):
        psf_dir_tmp = os.path.join(outpath_dir, f'psf_tile{tile_label}_band{band}')
        psf_ima_file_tmp = os.path.join(psf_dir_tmp, f'expo{i_expo}.fits')
        if os.path.isfile(psf_ima_file_tmp):
            logger.info(f'PSF image already exist.')
        else:
            # initial hdul
            hdu_list = fits.HDUList(fits.PrimaryHDU())
            ## add a card for pixel scale
            hdu_list[0].header['GS_SCALE'] = (pixel_scale, 'GalSim image scale')
            # produce 32 psf images
            for i_chip in range(32):
                seeing = seeing_chips[i_chip]
                beta = beta_chips[i_chip]
                psf_e = [psf_e_chips[0][i_chip], psf_e_chips[1][i_chip]]

                PSF = PSFModule.MoffatPSF(seeing, moffat_beta=beta, psf_e=psf_e)
                psf_ima = PSFModule.PSFima(PSF, pixel_scale, size=image_PSF_size)

                # collect to hdul
                galsim.fits.write(psf_ima, hdu_list=hdu_list)
                ## save chip id
                hdu_list[i_chip+1].header['IMAGEID'] = i_chip+1

            # save psf image
            hdu_list.writeto(psf_ima_file_tmp)

    ## if all exist, quit
    if (not False in outpath_image_exist_list):
        logger.info("All desired images exist, end the process.")
        return 1

    # +++ background noise 
    ## same for all chips
    noise = NoiseModule.GaussianNoise(rms, rng_seed=int(rng_seed_band+120*gal_rotation_angle+94*i_expo))

    # +++ sky image
    if (False in outpath_image_exist_list):

        # a list of canvas based on galaxy sky positions
        RA_gals = gals_info_band['RA'] # degree
        DEC_gals = gals_info_band['DEC'] # degree
        RA0 = (np.max(RA_gals) + np.min(RA_gals))/2.
        DEC0 = (np.max(DEC_gals) + np.min(DEC_gals))/2.
        canvases_list = KiDSModule.getKiDScanvases(RA0, DEC0, id_exposure=i_expo)
        del RA0, DEC0

        # all desired images
        for i_ima, outpath_image_exist in enumerate(outpath_image_exist_list):

            if (not outpath_image_exist):
                outpath_image_name = outpath_image_name_list[i_ima]

                # chip id
                i_chip = int(re.search(r'_chip(\d+)', outpath_image_name).group(1))

                # get the canvas accordingly
                canvas = canvases_list[i_chip]

                # get PSF
                seeing = seeing_chips[i_chip]
                beta = beta_chips[i_chip]
                psf_e = [psf_e_chips[0][i_chip], psf_e_chips[1][i_chip]]
                PSF = PSFModule.MoffatPSF(seeing, moffat_beta=beta, psf_e=psf_e)

                # galaxy image
                image_galaxies = ObjModule.GalaxiesImage(canvas, band, pixel_scale, PSF,
                                                gals_info_band, gal_rotation_angle=gal_rotation_angle, g_cosmic=g_cosmic, gal_position_type=gal_position_type,
                                                g_const=g_const)

                ## add stars
                if (stars_info_band is not None):
                    image_stars = ObjModule.StarsImage(canvas, band, pixel_scale, PSF, stars_info_band)
                    image_galaxies += image_stars
                    del image_stars

                ## add noise background
                image_galaxies.addNoise(noise)

                ## save the noisy image
                image_galaxies.write(outpath_image_name)
                logger.info(f"KiDS-like chip saved as {outpath_image_name}")

                ## mark success to the header
                with fits.open(outpath_image_name, mode='update') as hdul:
                    head_tmp = hdul[0].header
                    ## update info
                    head_tmp['flag_sim'] = 1

    logger.info(f'Finished for tile {tile_label} band {band} expo {i_expo} rot {gal_rotation_angle}...')
    return 0
