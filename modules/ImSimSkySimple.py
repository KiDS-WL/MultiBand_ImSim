# -*- coding: utf-8 -*-
# @Author: lshuns
# @Date:   2021-07-22 13:34:12
# @Last Modified by:   lshuns
# @Last Modified time: 2021-07-22 13:47:24

### Everything about simple images
__all__ = ['_PSFNoisySkyImages_simple']

###### dependence:
######      ImSimPSF: everthing about PSF
######      ImSimObject: everything about celestial objects
######      ImSimNoiseBackground: everthing about background noise
import ImSimPSF as PSFModule
import ImSimObject as ObjModule
import ImSimNoiseBackground as NoiseModule

import os
import galsim
import logging

import numpy as np

from astropy.io import fits

logger = logging.getLogger(__name__)

def _PSFNoisySkyImages_simple(para_list):
    '''
    Sky image with Gaussian noise and PSF
        Simple image without any survey strategy
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
        gal_position_type) = para_list

    logger.info(f'Simulating simple image for tile {tile_label} band {band} rot {gal_rotation_angle}...')

    # warning
    if save_image_chips:
        logger.warning('Simple image type does not produce chips!')

    # outpath
    outpath_image_name = outpath_image_basename + f'_rot{gal_rotation_angle:.0f}.fits'
    ## psf map
    if (outpath_PSF_basename is not None):
        outpath_PSF_name = outpath_PSF_basename + f'_rot{gal_rotation_angle:.0f}.fits'

    # first check if already exist
    outpath_image_exist = False
    try:
        with fits.open(outpath_image_name) as hdul:
            head_tmp = hdul[0].header
        flag_sim = head_tmp['flag_sim']
        if flag_sim >= 1:
            outpath_image_exist = True
            logger.info(f"{outpath_image_name} already exist.")
    except FileNotFoundError:
        pass
    except (KeyError, OSError) as e:
        os.remove(outpath_image_name)
        logger.warning(f"{outpath_image_name} broken, delete and re-produce.")
    ## psf map
    if (outpath_PSF_basename is not None) :
        outpath_PSF_exist = False
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
            logger.warning(f"{outpath_PSF_name} broken, delete and re-produce.")
    else:
        outpath_PSF_exist = True

    ## psf image
    ### different rotation has same psf, so only make once
    if (save_image_PSF) and (gal_rotation_angle==0.):
        psf_dir_tmp = os.path.join(outpath_dir, f'psf_tile{tile_label}_band{band}')
        psf_ima_file_tmp = os.path.join(psf_dir_tmp, f'psf_ima.fits')
        if os.path.isfile(psf_ima_file_tmp):
            logger.info('PSF images already exist.')
        else:
            PSF = PSFModule.MoffatPSF(seeing, moffat_beta=beta, psf_e=psf_e)
            psf_ima = PSFModule.PSFima(PSF, pixel_scale, size=image_PSF_size)
            psf_ima.write(psf_ima_file_tmp)
            logger.info(f'PSF image saved as {psf_ima_file_tmp}')

    ## if all exist, quit
    if (outpath_image_exist) and (outpath_PSF_exist):
        logger.info("All desired images exist, end the process.")
        return 1

    # +++ PSF
    PSF = PSFModule.MoffatPSF(seeing, moffat_beta=beta, psf_e=psf_e)

    # +++ background noise
    noise = NoiseModule.GaussianNoise(rms, rng_seed=int(rng_seed_band+120*gal_rotation_angle))

    # +++ PSF map
    if (not outpath_PSF_exist):
        mag_PSF_2 = 18. # for noise_flux = 2
        mag_PSF = mag_PSF_2 - 2.5*np.log10(rms/2.)
        image_PSF = PSFModule.PSFmap(PSF, pixel_scale, mag_PSF, N_PSF=N_PSF, sep_PSF=sep_PSF, rng_seed=rng_seed_band)

        ## noise background
        image_PSF_rot.addNoise(noise)

        ## save
        image_PSF.write(outpath_PSF_name)
        logger.debug(f"PSF map saved as {outpath_PSF_name}")

        ## mark success to the header
        with fits.open(outpath_PSF_name, mode='update') as hdul:
            head_tmp = hdul[0].header
            ## update info
            head_tmp['flag_sim'] = 1

    # +++ sky image
    if (not outpath_image_exist):

        # bounds and wcs from galaxy sky positions
        RA_gals = gals_info_band['RA'] # degree
        DEC_gals = gals_info_band['DEC'] # degree
        RA_min = np.min(RA_gals)
        RA_max = np.max(RA_gals)
        DEC_min = np.min(DEC_gals)
        DEC_max = np.max(DEC_gals)
        canvas = ObjModule.SimpleCanvas(RA_min, RA_max, DEC_min, DEC_max, pixel_scale)

        # star image
        if (stars_info_band is not None):
            image_stars = ObjModule.StarsImage(canvas, band, pixel_scale, PSF, stars_info_band)
        else:
            image_stars = None

        # galaxy images
        image_galaxies = ObjModule.GalaxiesImage(canvas, band, pixel_scale, PSF,
                                        gals_info_band, gal_rotation_angle=gal_rotation_angle, g_cosmic=g_cosmic, gal_position_type=gal_position_type)

        ## add stars
        if (image_stars is not None):
            image_galaxies += image_stars

        ## add noise background
        image_galaxies.addNoise(noise)

        ## save the noisy image
        image_galaxies.write(outpath_image_name)
        logger.info(f"Simple noisy sky image saved as {outpath_image_name}")

        ## mark success to the header
        with fits.open(outpath_image_name, mode='update') as hdul:
            head_tmp = hdul[0].header
            ## update info
            head_tmp['flag_sim'] = 1

    logger.info(f'Finished for tile {tile_label} band {band} rot {gal_rotation_angle}...')
    return 0