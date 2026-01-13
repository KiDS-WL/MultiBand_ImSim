# -*- coding: utf-8 -*-
# @Author: lshuns
# @Date:   2021-07-22 13:25:05
# @Last Modified by:   lshuns
# @Last Modified time: 2026-01-07 16:49:14

### Everything about KiDS-like images
__all__ = ['_PSFNoisySkyImages_KiDS_sameExpo', '_PSFNoisySkyImages_KiDS_singleExpo', '_PSFNoisySkyImages_KiDS_varChips']

###### dependence:
######      ImSimPSF: everything about PSF
######      ImSimObject: everything about celestial objects
######      ImSimNoiseBackground: everything about background noise
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
        rms, psf_info,
        g_cosmic,
        gals_info_band, gal_rotation_angle,
        stars_info_band,
        outpath_PSF_basename, N_PSF, sep_PSF,
        save_image_chips, save_image_PSF, image_PSF_size,
        save_image_noise,
        outpath_dir,
        gal_position_type,
        g_const,
        SimpleCam) = para_list
    
    assert not save_image_noise, 'KiDS_sameExpo does not support for save noise_image for now!'

    logger.info(f'Simulating KiDS_sameExpo image for tile {tile_label} band {band} rot {gal_rotation_angle}...')

    # PSF profiles
    if psf_info[0].lower() == 'moffat':

        seeing, beta, psf_e = psf_info[1:]
        # if psf e is zero, replace with None
        if (psf_e[0] == 0.) and (psf_e[1] == 0.):
            psf_e = None

        psf_paras = (seeing, beta, psf_e)
        psf_func = PSFModule.MoffatPSF
        psf_pixel = False

    elif psf_info[0].lower() == 'airy':

        lam, diam, obscuration, psf_e = psf_info[1:]
        # if psf e is zero, replace with None
        if (psf_e[0] == 0.) and (psf_e[1] == 0.):
            psf_e = None

        psf_paras = (lam, diam, obscuration, psf_e)
        psf_func = PSFModule.AiryPSF
        psf_pixel = False

    elif psf_info[0].lower() == 'pixelima':

        psf_fits_file = psf_info[1]

        psf_paras = (psf_fits_file, pixel_scale, (0.5, 0.5))
        psf_func = PSFModule.loadPixelPSF
        psf_pixel = True

    # number of exposures
    if band == 'u':
        n_exposures = 4
    elif band in ['g', 'r', 'i', 'i1', 'i2']:
        n_exposures = 5
    else:
        raise Exception(f'{band} is not in OmegaCAM! Cannot use KiDS_sameExpo!')

    # outpath
    outpath_image_name_list = [outpath_image_basename + f'_rot{gal_rotation_angle:.0f}_expo{id_exposure}.fits'
                                                                for id_exposure in range(n_exposures)]

    if (outpath_PSF_basename is not None):
        outpath_PSF_name_list = [outpath_PSF_basename + f'_rot{gal_rotation_angle:.0f}_expo{id_exposure}.fits'
                                                                for id_exposure in range(n_exposures)]

    # first check if already exist
    outpath_image_exist_list = np.zeros_like(outpath_image_name_list, dtype=bool)
    for i_name, outpath_image_name in enumerate(outpath_image_name_list):
        try:
            with fits.open(outpath_image_name) as hdul:
                head_tmp = hdul[0].header
            flag_sim = head_tmp['flag_sim']
            if flag_sim >= 1:
                outpath_image_exist_list[i_name] = True
                logger.debug(f"{outpath_image_name} already exist.")
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
                    logger.debug(f"{outpath_PSF_name} already exist.")
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

            PSF = psf_func(*psf_paras)
            psf_ima = PSFModule.PSFima(PSF, pixel_scale, size=image_PSF_size, pixelPSF=psf_pixel)

            for id_exposure in range(n_exposures):
                outpath_tmp = os.path.join(psf_dir_tmp, f'expo{id_exposure}.fits')
                psf_ima.write(outpath_tmp)
            logger.debug(f'PSF images saved to {psf_dir_tmp}')

    ## chips
    if save_image_chips:
        for i_ima, outpath_image_exist in enumerate(outpath_image_exist_list):
            if (outpath_image_exist):
                outpath_image_name = outpath_image_name_list[i_ima]
                id_exposure = int(re.search(r'_expo(\d+)', outpath_image_name).group(1))

                chip_dir_tmp = os.path.join(outpath_dir, f'chips_tile{tile_label}_band{band}_rot{gal_rotation_angle:.0f}')
                n_files = len(glob.glob(os.path.join(chip_dir_tmp, f'exp{id_exposure}chip_*.fits')))
                if n_files == 32:
                    logger.debug(f'chips already exist for rot{gal_rotation_angle:.0f} expo{id_exposure}.')
                else:
                    image_tile = galsim.fits.read(outpath_image_name)
                    image_chips = KiDSModule.cutKiDSchips(image_tile)

                    for i_chip, image_chip in enumerate(image_chips):
                        outpath_tmp = os.path.join(chip_dir_tmp, f'exp{id_exposure}chip_{i_chip+1}OFCS.fits')
                        image_chip.write(outpath_tmp)
                    logger.debug(f'Image chips saved to {chip_dir_tmp} for expo{id_exposure}.')

    ## if all exist, quit
    if (not False in outpath_image_exist_list) and (not False in outpath_PSF_exist_list):
        logger.info("All desired images exist, end the process.")
        return 1

    # +++ PSF
    PSF = psf_func(*psf_paras)

    # +++ background noise
    noise_list = [NoiseModule.GaussianNoise(rms, rng_seed=int(rng_seed_band+120*gal_rotation_angle+94*id_exposure))
                                                        for id_exposure in range(n_exposures)]

    # +++ PSF map
    if (False in outpath_PSF_exist_list):
        mag_PSF_2 = 18. # for noise_flux = 2
        mag_PSF = mag_PSF_2 - 2.5*np.log10(rms/2.)
        image_PSF = PSFModule.PSFmap(PSF, pixel_scale, mag_PSF, 
                            N_PSF=N_PSF, sep_PSF=sep_PSF, rng_seed=rng_seed_band,
                            pixelPSF=psf_pixel)

        for i_ima, outpath_PSF_exist in enumerate(outpath_PSF_exist_list):
            if (not outpath_PSF_exist):
                image_PSF_tmp = image_PSF.copy()

                ## noise background
                noise_tmp = noise_list[i_ima]
                image_PSF_tmp.addNoise(noise_tmp)

                ## save
                outpath_PSF_name = outpath_PSF_name_list[i_ima]
                image_PSF_tmp.write(outpath_PSF_name)
                logger.debug(f"PSF map saved as {outpath_PSF_name}.")

                ## mark success to the header
                with fits.open(outpath_PSF_name, mode='update') as hdul:
                    head_tmp = hdul[0].header
                    ## update info
                    head_tmp['flag_sim'] = 1

    # +++ sky image
    if (False in outpath_image_exist_list):

        # simple canvas based on the galaxy sky positions
        RA_gals = gals_info_band[0]['RA'].values
        DEC_gals = gals_info_band[0]['DEC'].values
        if gals_info_band[1] is not None:
            RA_gals = np.hstack([RA_gals, gals_info_band[1]['RA'].values])
            DEC_gals = np.hstack([DEC_gals, gals_info_band[1]['DEC'].values])
        RA_min = np.amin(RA_gals)
        RA_max = np.amax(RA_gals)
        DEC_min = np.amin(DEC_gals)
        DEC_max = np.amax(DEC_gals)
        canvas = ObjModule.SimpleCanvas(RA_min, RA_max, DEC_min, DEC_max, pixel_scale)
        del RA_gals, DEC_gals, RA_min, RA_max, DEC_min, DEC_max

        # star image
        if (stars_info_band is not None):
            image_stars = ObjModule.StarsImage(canvas, band, pixel_scale, PSF, 
                                            stars_info_band,
                                            pixelPSF=psf_pixel)
        else:
            image_stars = None

        # galaxy images
        image_galaxies0 = ObjModule.GalaxiesImage(canvas, band, pixel_scale, PSF,
                                gals_info_band[0], gal_rotation_angle=gal_rotation_angle, 
                                g_cosmic=g_cosmic, gal_position_type=gal_position_type,
                                g_const=g_const,
                                pixelPSF=psf_pixel)
        if gals_info_band[1] is not None:
            image_galaxies0 += ObjModule.GalaxiesImage_casual(canvas, band, pixel_scale, PSF,
                                gals_info_band[1], gal_rotation_angle=gal_rotation_angle, 
                                g_cosmic=g_cosmic, gal_position_type=gal_position_type,
                                g_const=g_const, 
                                pixelPSF=psf_pixel)

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
                logger.debug(f"KiDS-like exposure saved as {outpath_image_name}")

                ## save the weight image
                outpath_wei_name = outpath_image_name.replace('.fits', '.weight.fits')
                weights_tile.write(outpath_wei_name)
                logger.debug(f"weight image saved as {outpath_wei_name}")

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
                        outpath_tmp = os.path.join(chip_dir_tmp, f'exp{id_exposure}chip_{i_chip+1}OFCS.fits')
                        image_chip.write(outpath_tmp)
                    logger.debug(f'Image chips saved to {chip_dir_tmp} for exposure {id_exposure}.')

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
        rms, psf_info,
        g_cosmic,
        gals_info_band, gal_rotation_angle,
        stars_info_band,
        outpath_PSF_basename, N_PSF, sep_PSF,
        save_image_PSF, image_PSF_size,
        save_image_noise,
        outpath_dir,
        id_exposure,
        gal_position_type,
        g_const,
        SimpleCam) = para_list

    assert not save_image_noise, 'KiDS_singleExpo does not support for save noise_image for now!'

    logger.info(f'Simulating KiDS exposure for tile {tile_label} band {band} expo {id_exposure} rot {gal_rotation_angle}...')

    # PSF profiles
    if psf_info[0].lower() == 'moffat':

        seeing, beta, psf_e = psf_info[1:]
        # if psf e is zero, replace with None
        if (psf_e[0] == 0.) and (psf_e[1] == 0.):
            psf_e = None

        psf_paras = (seeing, beta, psf_e)
        psf_func = PSFModule.MoffatPSF
        psf_pixel = False

    elif psf_info[0].lower() == 'airy':

        lam, diam, obscuration, psf_e = psf_info[1:]
        # if psf e is zero, replace with None
        if (psf_e[0] == 0.) and (psf_e[1] == 0.):
            psf_e = None

        psf_paras = (lam, diam, obscuration, psf_e)
        psf_func = PSFModule.AiryPSF
        psf_pixel = False

    elif psf_info[0].lower() == 'pixelima':

        psf_fits_file = psf_info[1]

        psf_paras = (psf_fits_file, pixel_scale, (0.5, 0.5))
        psf_func = PSFModule.loadPixelPSF
        psf_pixel = True

    # outpath
    outpath_image_name_list = [os.path.join(outpath_dir, f'chips_tile{tile_label}_band{band}_rot{gal_rotation_angle:.0f}', f'exp{id_exposure}chip_{i_chip+1}OFCS.fits')
                                for i_chip in range(32)]

    if (outpath_PSF_basename is not None):
        outpath_PSF_name = outpath_PSF_basename + f'_rot{gal_rotation_angle:.0f}_expo{id_exposure}.fits'

    # first check if already exist
    outpath_image_exist_list = np.zeros_like(outpath_image_name_list, dtype=bool)
    for i_name, outpath_image_name in enumerate(outpath_image_name_list):
        try:
            with fits.open(outpath_image_name) as hdul:
                head_tmp = hdul[0].header
            flag_sim = head_tmp['flag_sim']
            if flag_sim >= 1:
                outpath_image_exist_list[i_name] = True
                logger.debug(f"{outpath_image_name} already exist.")
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
                logger.debug(f"{outpath_PSF_name} already exist.")
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
        psf_ima_file_tmp = os.path.join(psf_dir_tmp, f'expo{id_exposure}.fits')
        if os.path.isfile(psf_ima_file_tmp):
            logger.info('PSF image already exist.')
        else:
            PSF = psf_func(*psf_paras)
            psf_ima = PSFModule.PSFima(PSF, pixel_scale, size=image_PSF_size, pixelPSF=psf_pixel)
            psf_ima.write(psf_ima_file_tmp)
            logger.debug(f'PSF image saved as {psf_ima_file_tmp}')

    ## if all exist, quit
    if (not False in outpath_image_exist_list) and (outpath_PSF_exist):
        logger.info("All desired images exist, end the process.")
        return 1

    # +++ PSF
    PSF = psf_func(*psf_paras)

    # +++ background noise
    noise = NoiseModule.GaussianNoise(rms, rng_seed=int(rng_seed_band+120*gal_rotation_angle+94*id_exposure))

    # +++ PSF map
    if (not outpath_PSF_exist):
        mag_PSF_2 = 18. # for noise_flux = 2
        mag_PSF = mag_PSF_2 - 2.5*np.log10(rms/2.)
        image_PSF = PSFModule.PSFmap(PSF, pixel_scale, mag_PSF, 
            N_PSF=N_PSF, sep_PSF=sep_PSF, rng_seed=rng_seed_band,
            pixelPSF=psf_pixel)

        ## noise background
        image_PSF.addNoise(noise)

        ## save
        image_PSF.write(outpath_PSF_name)
        logger.debug(f"PSF map saved as {outpath_PSF_name}.")

        ## mark success to the header
        with fits.open(outpath_PSF_name, mode='update') as hdul:
            head_tmp = hdul[0].header
            ## update info
            head_tmp['flag_sim'] = 1

    # +++ sky image
    if (False in outpath_image_exist_list):

        # a list of canvas based on galaxy sky positions
        RA_gals = gals_info_band[0]['RA'].values
        DEC_gals = gals_info_band[0]['DEC'].values
        if gals_info_band[1] is not None:
            RA_gals = np.hstack([RA_gals, gals_info_band[1]['RA'].values])
            DEC_gals = np.hstack([DEC_gals, gals_info_band[1]['DEC'].values])
        RA0 = (np.amax(RA_gals) + np.amin(RA_gals))/2.
        DEC0 = (np.amax(DEC_gals) + np.amin(DEC_gals))/2.
        canvases_list = KiDSModule.getKiDScanvases(RA0, DEC0, SimpleCam, id_exposure=id_exposure)
        del RA_gals, DEC_gals, RA0, DEC0

        # all desired images
        for i_ima, outpath_image_exist in enumerate(outpath_image_exist_list):

            if (not outpath_image_exist):
                outpath_image_name = outpath_image_name_list[i_ima]

                # chip id
                i_chip = int(re.search(r'chip_(\d+)', outpath_image_name).group(1)) - 1
                ## get the canvas accordingly
                canvas = canvases_list[i_chip]

                # galaxy image
                image_galaxies = ObjModule.GalaxiesImage(canvas, band, pixel_scale, PSF,
                                                gals_info_band[0], gal_rotation_angle=gal_rotation_angle, 
                                                g_cosmic=g_cosmic, gal_position_type=gal_position_type,
                                                g_const=g_const,
                                                pixelPSF=psf_pixel)
                if gals_info_band[1] is not None:
                    image_galaxies += ObjModule.GalaxiesImage_casual(canvas, band, pixel_scale, PSF,
                                        gals_info_band[1], gal_rotation_angle=gal_rotation_angle, 
                                        g_cosmic=g_cosmic, gal_position_type=gal_position_type,
                                        g_const=g_const, 
                                        pixelPSF=psf_pixel)

                ## add stars
                if (stars_info_band is not None):
                    image_stars = ObjModule.StarsImage(canvas, band, pixel_scale, PSF, stars_info_band,
                                                    pixelPSF=psf_pixel)
                    image_galaxies += image_stars
                    del image_stars

                ## add noise background
                image_galaxies.addNoise(noise)

                ## save the noisy image
                image_galaxies.write(outpath_image_name)
                logger.debug(f"KiDS-like chip saved as {outpath_image_name}")

                ## mark success to the header
                with fits.open(outpath_image_name, mode='update') as hdul:
                    head_tmp = hdul[0].header
                    ## update info
                    head_tmp['flag_sim'] = 1

    logger.info(f'Finished for tile {tile_label} band {band} expo {id_exposure} rot {gal_rotation_angle}...')
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
        rms, psf_info_chips,
        g_cosmic,
        gals_info_band, gal_rotation_angle,
        stars_info_band,
        outpath_PSF_basename, N_PSF, sep_PSF,
        save_image_PSF, image_PSF_size,
        save_image_noise,
        outpath_dir,
        id_exposure,
        gal_position_type,
        g_const,
        SimpleCam) = para_list

    assert not save_image_noise, 'KiDS_varChips does not support for save noise_image for now!'

    logger.info(f'Simulating KiDS exposure with varChips for tile {tile_label} band {band} expo {id_exposure} rot {gal_rotation_angle}...')

    # PSF profiles
    if psf_info_chips[0].lower() == 'moffat':

        # function
        psf_func = PSFModule.MoffatPSF
        psf_pixel = False

        # paras
        seeing_chips, beta_chips, psf_e_chips = psf_info_chips[1:]
        psf_paras_chips = []
        for i_chip in range(32):
            seeing = seeing_chips[i_chip]
            beta = beta_chips[i_chip]
            psf_e = [psf_e_chips[0][i_chip], psf_e_chips[1][i_chip]]
            # if psf e is zero, replace with None
            if (psf_e[0] == 0.) and (psf_e[1] == 0.):
                psf_e = None

            psf_paras_chips.append((seeing, beta, psf_e))
            del seeing, beta, psf_e
        del seeing_chips, beta_chips, psf_e_chips, psf_info_chips

    elif psf_info_chips[0].lower() == 'airy':

        # function
        psf_func = PSFModule.AiryPSF
        psf_pixel = False

        # paras
        lam_chips, diam_chips, obscuration_chips, psf_e_chips = psf_info_chips[1:]
        psf_paras_chips = []
        for i_chip in range(32):
            lam = lam_chips[i_chip]
            diam = diam_chips[i_chip]
            obscuration = obscuration_chips[i_chip]
            psf_e = [psf_e_chips[0][i_chip], psf_e_chips[1][i_chip]]
            # if psf e is zero, replace with None
            if (psf_e[0] == 0.) and (psf_e[1] == 0.):
                psf_e = None

            psf_paras_chips.append((lam, diam, obscuration, psf_e))
            del lam, diam, obscuration, psf_e
        del lam_chips, diam_chips, obscuration_chips, psf_e_chips, psf_info_chips

    elif psf_info_chips[0].lower() == 'pixelima':

        # function
        psf_func = PSFModule.loadPixelPSF
        psf_pixel = True

        # paras
        psf_fits_file_chips = psf_info_chips[1]
        psf_paras_chips = []
        for i_chip in range(32):
            psf_paras_chips.append((psf_fits_file_chips[i_chip], pixel_scale, (0.5, 0.5)))
        del psf_fits_file_chips, psf_info_chips

    # outpath
    outpath_image_name_list = [os.path.join(outpath_dir, f'chips_tile{tile_label}_band{band}_rot{gal_rotation_angle:.0f}', f'exp{id_exposure}chip_{i_chip+1}OFCS.fits')
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
                logger.debug(f"{outpath_image_name} already exist.")
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
        psf_ima_file_tmp = os.path.join(psf_dir_tmp, f'exp{id_exposure}chip.fits')
        if os.path.isfile(psf_ima_file_tmp):
            logger.info(f'PSF image already exist.')
        else:
            # initial hdul
            hdu_list = fits.HDUList(fits.PrimaryHDU())
            ## add a card for pixel scale
            hdu_list[0].header['GS_SCALE'] = (pixel_scale, 'GalSim image scale')
            # produce 32 psf images
            for i_chip in range(32):
                psf_paras = psf_paras_chips[i_chip]

                PSF = psf_func(*psf_paras)
                psf_ima = PSFModule.PSFima(PSF, pixel_scale, size=image_PSF_size, 
                                    pixelPSF=psf_pixel)

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
    noise = NoiseModule.GaussianNoise(rms, rng_seed=int(rng_seed_band+120*gal_rotation_angle+94*id_exposure))

    # +++ sky image
    if (False in outpath_image_exist_list):

        # a list of canvas based on galaxy sky positions
        RA_gals = gals_info_band[0]['RA'].values
        DEC_gals = gals_info_band[0]['DEC'].values
        if gals_info_band[1] is not None:
            RA_gals = np.hstack([RA_gals, gals_info_band[1]['RA'].values])
            DEC_gals = np.hstack([DEC_gals, gals_info_band[1]['DEC'].values])
        RA0 = (np.amax(RA_gals) + np.amin(RA_gals))/2.
        DEC0 = (np.amax(DEC_gals) + np.amin(DEC_gals))/2.
        canvases_list = KiDSModule.getKiDScanvases(RA0, DEC0, SimpleCam, id_exposure=id_exposure)
        del RA_gals, DEC_gals, RA0, DEC0

        # all desired images
        for i_ima, outpath_image_exist in enumerate(outpath_image_exist_list):

            if (not outpath_image_exist):
                outpath_image_name = outpath_image_name_list[i_ima]

                # chip id
                i_chip = int(re.search(r'chip_(\d+)', outpath_image_name).group(1)) - 1

                # get the canvas accordingly
                canvas = canvases_list[i_chip]

                # get PSF
                psf_paras = psf_paras_chips[i_chip]
                PSF = psf_func(*psf_paras)

                # galaxy image
                image_galaxies = ObjModule.GalaxiesImage(canvas, band, pixel_scale, PSF,
                                                gals_info_band[0], gal_rotation_angle=gal_rotation_angle, 
                                                g_cosmic=g_cosmic, gal_position_type=gal_position_type,
                                                g_const=g_const,
                                                pixelPSF=psf_pixel)
                if gals_info_band[1] is not None:
                    image_galaxies += ObjModule.GalaxiesImage_casual(canvas, band, pixel_scale, PSF,
                                        gals_info_band[1], gal_rotation_angle=gal_rotation_angle, 
                                        g_cosmic=g_cosmic, gal_position_type=gal_position_type,
                                        g_const=g_const, 
                                        pixelPSF=psf_pixel)

                ## add stars
                if (stars_info_band is not None):
                    image_stars = ObjModule.StarsImage(canvas, band, pixel_scale, PSF, stars_info_band,
                                                pixelPSF=psf_pixel)
                    image_galaxies += image_stars
                    del image_stars

                ## add noise background
                image_galaxies.addNoise(noise)

                ## save the noisy image
                image_galaxies.write(outpath_image_name)
                logger.debug(f"KiDS-like chip saved as {outpath_image_name}")

                ## mark success to the header
                with fits.open(outpath_image_name, mode='update') as hdul:
                    head_tmp = hdul[0].header
                    ## update info
                    head_tmp['flag_sim'] = 1

    logger.info(f'Finished for tile {tile_label} band {band} expo {id_exposure} rot {gal_rotation_angle}...')
    return 0
