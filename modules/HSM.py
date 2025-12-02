# -*- coding: utf-8 -*-
# @Author: lshuns
# @Date:   1969-12-31 16:00:00
# @Last Modified by:   lshuns
# @Last Modified time: 2025-12-02 15:32:08

### Galaxty shape measurement using HSM module from GalSim

import os
import logging

import galsim
import numpy as np
import pandas as pd
from astropy.io import fits
from multiprocessing import shared_memory
from concurrent.futures import ProcessPoolExecutor

logger = logging.getLogger(__name__)

# Function to model sigma_e as a function of SNR
def _sigma_e_from_SNR(SNR, a, p, b):
    return a / (SNR**p) + b

# >>>>>>>>>>>>>>>>>>>>>>>>>>>> AdaptiveMomShape <<<<<<<<<<<<<<<<<<<<<<<<<<<<
# Worker function of AdaptiveMomShape
def _worker_ams(args):
    (idx, 
     xcen, ycen, id_detec_row, 
     shm_name_img, img_shape, img_dtype,
     shm_name_seg, seg_shape, seg_dtype,
     postage_size, guess_sig, precision, round_moments,
     outpath_stamp) = args

    # re-attach shared arrays
    shm_img = shared_memory.SharedMemory(name=shm_name_img)
    image_data = np.ndarray(img_shape, dtype=np.dtype(img_dtype), buffer=shm_img.buf)
    shm_seg = shared_memory.SharedMemory(name=shm_name_seg)
    seg_data = np.ndarray(seg_shape, dtype=np.dtype(seg_dtype), buffer=shm_seg.buf)

    # extract postage stamp
    xcen_int = int(np.round(xcen))
    ycen_int = int(np.round(ycen))
    half_size = postage_size // 2
    ## check boundaries
    if (xcen_int - half_size < 0 or
        ycen_int - half_size < 0 or
        xcen_int + half_size >= img_shape[1] or
        ycen_int + half_size >= img_shape[0]):
        shm_img.close(); shm_seg.close()
        return (idx, -999, -999, -999, -999, -999)    
    stamp_image = image_data[ycen_int-half_size:ycen_int+half_size,
                             xcen_int-half_size:xcen_int+half_size]
    stamp_seg = seg_data[ycen_int-half_size:ycen_int+half_size,
                         xcen_int-half_size:xcen_int+half_size]
    stamp_clean = stamp_image * ((stamp_seg == id_detec_row) | (stamp_seg == 0))
    ## check shape
    if stamp_image.shape != (postage_size, postage_size):
        shm_img.close(); shm_seg.close()
        return (idx, -999., -999., -999., -999., -999.)

    # save postage stamp if required
    if outpath_stamp is not None:
        hdu = fits.PrimaryHDU(stamp_image)
        hdu.writeto(outpath_stamp, overwrite=True)

        hdu = fits.PrimaryHDU(stamp_seg)
        base = os.path.splitext(os.path.basename(outpath_stamp))[0]
        image_fname = base.replace('_image', '_seg') + '.fits'
        outpath_tmp = os.path.join(os.path.dirname(outpath_stamp), image_fname)
        hdu.writeto(outpath_tmp, overwrite=True)

        hdu = fits.PrimaryHDU(stamp_clean)
        base = os.path.splitext(os.path.basename(outpath_stamp))[0]
        image_fname = base.replace('_image', '_clean_image') + '.fits'
        outpath_tmp = os.path.join(os.path.dirname(outpath_stamp), image_fname)
        hdu.writeto(outpath_tmp, overwrite=True)
        logger.debug(f'Saved postage stamp for object ID {int(id_detec_row)}')

    # run FindAdaptiveMom
    try:
        object_moments = galsim.hsm.FindAdaptiveMom(
            galsim.Image(stamp_clean),
            weight=None,
            badpix=None,
            guess_sig=guess_sig,
            precision=precision,
            guess_centroid=None,
            strict=False,
            check=True,
            round_moments=round_moments,
            hsmparams=None,
            use_sky_coords=False
        )
        if object_moments.moments_status == 0:
            shm_img.close(); shm_seg.close()
            return (idx,
                    object_moments.observed_e1,
                    object_moments.observed_e2,
                    object_moments.moments_sigma,
                    object_moments.moments_amp,
                    object_moments.moments_n_iter)
    except Exception:
        pass

    # failed
    ## return sentinel values
    shm_img.close(); shm_seg.close()
    return (idx, -999.0, -999.0, -999.0, -999.0, -999.0)

# main function for running AdaptiveMomShape
def AdaptiveMomShape(outpath_feather,
                    inpath_detection, 
                    inpath_image, 
                    inpath_seg_map,
                    inpath_weight_map=None, 
                    sigma_fromSNR_amp=1.16,
                    sigma_fromSNR_index=0.62,
                    sigma_fromSNR_base=0.0,
                    sigma_intrinsic=0.27,
                    guess_sig=5.0,
                    precision=1e-06,
                    round_moments=False,
                    save_Nstamps=0,
                    random_seed=914,
                    postage_size=48,
                    max_cores=12):

    logger.info(f'Measure adaptive moments for {os.path.basename(outpath_feather)}...')
    if os.path.isfile(outpath_feather):
        logger.info(f'The final feather catalogue {outpath_feather} already exists.')
        logger.info(f'End the process.')
        return 1

    # >>>>>>>>>>>>> 1. load detection catalogue and images
    ## load catalogue
    detec_cata = pd.read_feather(inpath_detection)[['NUMBER', 'X_IMAGE', 'Y_IMAGE',
                                                    'FLUX_AUTO', 'FLUXERR_AUTO']].reset_index(drop=True)
    logger.debug(f'Detection catalogue loaded from {inpath_detection}')

    ## load image, weight map, seg map
    with fits.open(inpath_image) as hdul:
        image_data = hdul[0].data
    with fits.open(inpath_seg_map) as hdul:
        seg_data = hdul[0].data
    logger.debug(f'Image, seg map loaded from {inpath_image}, {inpath_seg_map}') 

    ## mask image with weight map
    if inpath_weight_map is not None:
        with fits.open(inpath_weight_map) as hdul:
            weight_data = hdul[0].data
        image_data = np.where(weight_data > 0, image_data, 0.0)
        del weight_data  # free memory
        logger.debug(f'Image masked with weight map from {inpath_weight_map}')

    ## shared memory for big arrays
    shm_img = shared_memory.SharedMemory(create=True, size=image_data.nbytes)
    shm_img_arr = np.ndarray(image_data.shape, 
                             dtype=image_data.dtype, 
                             buffer=shm_img.buf)
    shm_seg = shared_memory.SharedMemory(create=True, size=seg_data.nbytes)
    shm_seg_arr = np.ndarray(seg_data.shape, 
                             dtype=seg_data.dtype, 
                             buffer=shm_seg.buf)
    shm_img_arr[:] = image_data
    shm_seg_arr[:] = seg_data
    del image_data, seg_data

    # >>>>>>>>>>>>> 2. run FindAdaptiveMom for each object
    ## randomly pick objects to save postage stamps 
    if save_Nstamps > 0: 
        np.random.seed(random_seed + save_Nstamps) 
        saved_indices = np.random.choice(detec_cata.index, 
                                         size=save_Nstamps, 
                                         replace=False) 
        stamp_dir = os.path.join(os.path.dirname(outpath_feather), 
                                 'stamps_AMS') 
        os.makedirs(stamp_dir, exist_ok=True) 
        logger.info(f'Saving {save_Nstamps} postage stamps to {stamp_dir}')

    ## prepare arguments for workers
    job_args = []
    for idx, row in detec_cata.iterrows():
        id_detec_row = row['NUMBER']
        if save_Nstamps > 0 and idx in saved_indices:
            base = os.path.splitext(os.path.basename(inpath_image))[0]
            image_fname = f"{base}_image_NUMBER{int(id_detec_row)}.fits"
            outpath_stamp = os.path.join(stamp_dir, image_fname)
        else:
            outpath_stamp = None
        job_args.append((
            idx,
            row['X_IMAGE'], row['Y_IMAGE'], id_detec_row,
            shm_img.name, shm_img_arr.shape, str(shm_img_arr.dtype),
            shm_seg.name, shm_seg_arr.shape, str(shm_seg_arr.dtype),
            postage_size, guess_sig, precision, round_moments,
            outpath_stamp))

    ## other useful parameters from detection catalogue
    id_detec = detec_cata['NUMBER'].values
    #### calculate weight based on SNR
    snr_detec = detec_cata['FLUX_AUTO'].values / detec_cata['FLUXERR_AUTO'].values
    weights = 1.0 / (_sigma_e_from_SNR(snr_detec, 
                                       sigma_fromSNR_amp, 
                                       sigma_fromSNR_index, 
                                       sigma_fromSNR_base)**2 + sigma_intrinsic**2)
    del detec_cata, snr_detec  # free memory

    ## run in parallel
    results = []
    with ProcessPoolExecutor(max_workers=max_cores) as exe:
        for res in exe.map(_worker_ams, job_args):
            results.append(res)

    # free shared memory
    shm_img.close(); shm_img.unlink()
    shm_seg.close(); shm_seg.unlink()

    # build output catalogue
    shape_cata = pd.DataFrame({
        'id_detec': id_detec,
        'e1_distortion_AMS': -999., 'e2_distortion_AMS': -999.,
        'moments_sigma_AMS': -999.,
        'moments_amp_AMS': -999.,
        'moments_n_iter_AMS': -999.,
        'weight_AMS': weights
    })
    del id_detec, weights  # free memory
    ## results from workers
    for idx, e1, e2, sig, amp, niter in results:
        shape_cata.loc[idx, ['e1_distortion_AMS','e2_distortion_AMS',
                             'moments_sigma_AMS',
                             'moments_amp_AMS',
                             'moments_n_iter_AMS']] = \
            [e1, e2, sig, amp, niter]
    ## set zero weight for failed measurements
    shape_cata.loc[shape_cata['e1_distortion_AMS'] == -999., 'weight_AMS'] = 0.0

    # save
    tmp = outpath_feather + '_tmp'
    shape_cata.to_feather(tmp)
    del shape_cata
    os.rename(tmp, outpath_feather)
    logger.info(f'Final catalogue saved as {outpath_feather}')
    return 0

# >>>>>>>>>>>>>>>>>>>>>>>>>>>> AdaptiveMomShape <<<<<<<<<<<<<<<<<<<<<<<<<<<<
# Worker function of EstimateShear
def _worker_hsm(args):
    (idx, 
     xcen, ycen, id_detec_row, 
     shm_name_img, img_shape, img_dtype,
     shm_name_seg, seg_shape, seg_dtype,
     psf_img_galsim,
     postage_size, shear_est, 
     recompute_flux, guess_sig_gal, guess_sig_PSF,
     precision,
     outpath_stamp) = args

    # re-attach shared arrays
    shm_img = shared_memory.SharedMemory(name=shm_name_img)
    image_data = np.ndarray(img_shape, dtype=np.dtype(img_dtype), buffer=shm_img.buf)
    shm_seg = shared_memory.SharedMemory(name=shm_name_seg)
    seg_data = np.ndarray(seg_shape, dtype=np.dtype(seg_dtype), buffer=shm_seg.buf)

    # extract postage stamp
    xcen_int = int(np.round(xcen))
    ycen_int = int(np.round(ycen))
    half_size = postage_size // 2
    ## check boundaries
    if (xcen_int - half_size < 0 or
        ycen_int - half_size < 0 or
        xcen_int + half_size >= img_shape[1] or
        ycen_int + half_size >= img_shape[0]):
        shm_img.close(); shm_seg.close()
        return (idx, -999.0, -999.0, -999.0, -999.0, -999.0, -999.0, -999.0, -999.0)
    stamp_image = image_data[ycen_int-half_size:ycen_int+half_size,
                             xcen_int-half_size:xcen_int+half_size]
    stamp_seg = seg_data[ycen_int-half_size:ycen_int+half_size,
                         xcen_int-half_size:xcen_int+half_size]
    stamp_clean = stamp_image * ((stamp_seg == id_detec_row) | (stamp_seg == 0))
    ## check shape
    if stamp_image.shape != (postage_size, postage_size):
        shm_img.close(); shm_seg.close()
        return (idx, -999.0, -999.0, -999.0, -999.0, -999.0, -999.0, -999.0, -999.0)

    # save postage stamp if required
    if outpath_stamp is not None:
        hdu = fits.PrimaryHDU(stamp_image)
        hdu.writeto(outpath_stamp, overwrite=True)

        hdu = fits.PrimaryHDU(stamp_seg)
        base = os.path.splitext(os.path.basename(outpath_stamp))[0]
        image_fname = base.replace('_image', '_seg') + '.fits'
        outpath_tmp = os.path.join(os.path.dirname(outpath_stamp), image_fname)
        hdu.writeto(outpath_tmp, overwrite=True)

        hdu = fits.PrimaryHDU(stamp_clean)
        base = os.path.splitext(os.path.basename(outpath_stamp))[0]
        image_fname = base.replace('_image', '_clean_image') + '.fits'
        outpath_tmp = os.path.join(os.path.dirname(outpath_stamp), image_fname)
        hdu.writeto(outpath_tmp, overwrite=True)
        logger.debug(f'Saved postage stamp for object ID {int(id_detec_row)}')

    # run EstimateShear
    try:
        object_moments = galsim.hsm.EstimateShear(
            galsim.Image(stamp_clean), psf_img_galsim,
            weight=None,
            badpix=None,
            sky_var=0.0, 
            shear_est=shear_est, 
            recompute_flux=recompute_flux, 
            guess_sig_gal=guess_sig_gal, 
            guess_sig_PSF=guess_sig_PSF, 
            precision=precision, 
            guess_centroid=None, 
            strict=False, 
            check=True, 
            hsmparams=None)
        if object_moments.correction_status == 0:
            shm_img.close(); shm_seg.close()
            if object_moments.meas_type == 'e':
                return (idx,
                        object_moments.corrected_e1,
                        object_moments.corrected_e2,
                        object_moments.resolution_factor,
                        object_moments.psf_sigma,
                        object_moments.psf_shape.g1,
                        object_moments.psf_shape.g2,
                        object_moments.moments_sigma,
                        object_moments.moments_amp)
            else:
                return (idx,
                        object_moments.corrected_g1,
                        object_moments.corrected_g2,
                        object_moments.resolution_factor,
                        object_moments.psf_sigma,
                        object_moments.psf_shape.g1,
                        object_moments.psf_shape.g2,
                        object_moments.moments_sigma,
                        object_moments.moments_amp)

    except Exception:
        pass

    # failed
    ## return sentinel values
    shm_img.close(); shm_seg.close()
    return (idx, -999.0, -999.0, -999.0, -999.0, -999.0, -999.0, -999.0, -999.0)

# main function for running EstimateShear
## assuming PSF is constant across the image
def EstimateShear_samePSF(outpath_feather,
                    inpath_detection, 
                    inpath_image, 
                    inpath_seg_map,  
                    inpath_psf_image,
                    inpath_weight_map=None, 
                    sigma_fromSNR_amp=1.16,
                    sigma_fromSNR_index=0.62,
                    sigma_fromSNR_base=0.0,
                    sigma_intrinsic=0.27,
                    shear_est='REGAUSS', 
                    recompute_flux='NONE', 
                    guess_sig_gal=5.0, 
                    guess_sig_PSF=3.0, 
                    precision=1e-06, 
                    save_Nstamps=0,
                    random_seed=914,
                    postage_size=48,
                    max_cores=12):

    logger.info(f'Measure shear using {shear_est} for {os.path.basename(outpath_feather)}...')
    if os.path.isfile(outpath_feather):
        logger.info(f'The final feather catalogue {outpath_feather} already exists.')
        logger.info(f'End the process.')
        return 1

    # >>>>>>>>>>>>> 1. load detection catalogue and images
    ## load catalogue
    detec_cata = pd.read_feather(inpath_detection)[['NUMBER', 'X_IMAGE', 'Y_IMAGE',
                                                    'FLUX_AUTO', 'FLUXERR_AUTO']].reset_index(drop=True)
    logger.debug(f'Detection catalogue loaded from {inpath_detection}')

    ## load image, weight map, seg map
    with fits.open(inpath_image) as hdul:
        image_data = hdul[0].data
    with fits.open(inpath_seg_map) as hdul:
        seg_data = hdul[0].data
    logger.debug(f'Image, seg map loaded from {inpath_image}, {inpath_seg_map}') 

    ## mask image with weight map
    if inpath_weight_map is not None:
        with fits.open(inpath_weight_map) as hdul:
            weight_data = hdul[0].data
        image_data = np.where(weight_data > 0, image_data, 0.0)
        del weight_data  # free memory
        logger.debug(f'Image masked with weight map from {inpath_weight_map}')

    ## shared memory for big arrays
    shm_img = shared_memory.SharedMemory(create=True, size=image_data.nbytes)
    shm_img_arr = np.ndarray(image_data.shape, 
                             dtype=image_data.dtype, 
                             buffer=shm_img.buf)
    shm_seg = shared_memory.SharedMemory(create=True, size=seg_data.nbytes)
    shm_seg_arr = np.ndarray(seg_data.shape, 
                             dtype=seg_data.dtype, 
                             buffer=shm_seg.buf)
    shm_img_arr[:] = image_data
    shm_seg_arr[:] = seg_data
    del image_data, seg_data

    ## load PSF image and convert to GalSim image
    with fits.open(inpath_psf_image) as hdul:
        psf_img_galsim = galsim.Image(hdul[0].data)
    logger.debug(f'PSF image loaded from {inpath_psf_image}')

    # >>>>>>>>>>>>> 2. run EstimateShear for each object
    ## randomly pick objects to save postage stamps 
    if save_Nstamps > 0: 
        np.random.seed(random_seed + save_Nstamps) 
        saved_indices = np.random.choice(detec_cata.index, 
                                         size=save_Nstamps, 
                                         replace=False) 
        stamp_dir = os.path.join(os.path.dirname(outpath_feather), 
                                 'stamps_AMS') 
        os.makedirs(stamp_dir, exist_ok=True) 
        logger.info(f'Saving {save_Nstamps} postage stamps to {stamp_dir}')

    ## prepare arguments for workers
    job_args = []
    for idx, row in detec_cata.iterrows():
        id_detec_row = row['NUMBER']
        if save_Nstamps > 0 and idx in saved_indices:
            base = os.path.splitext(os.path.basename(inpath_image))[0]
            image_fname = f"{base}_image_NUMBER{int(id_detec_row)}.fits"
            outpath_stamp = os.path.join(stamp_dir, image_fname)
        else:
            outpath_stamp = None
        job_args.append((
            idx,
            row['X_IMAGE'], row['Y_IMAGE'], id_detec_row,
            shm_img.name, shm_img_arr.shape, str(shm_img_arr.dtype),
            shm_seg.name, shm_seg_arr.shape, str(shm_seg_arr.dtype),
            psf_img_galsim,
            postage_size, shear_est, 
            recompute_flux, guess_sig_gal, guess_sig_PSF,
            precision,
            outpath_stamp))

    ## other useful parameters from detection catalogue
    id_detec = detec_cata['NUMBER'].values
    #### calculate weight based on SNR
    snr_detec = detec_cata['FLUX_AUTO'].values / detec_cata['FLUXERR_AUTO'].values
    weights = 1.0 / (_sigma_e_from_SNR(snr_detec, 
                                       sigma_fromSNR_amp, 
                                       sigma_fromSNR_index, 
                                       sigma_fromSNR_base)**2 + sigma_intrinsic**2)
    del detec_cata, snr_detec  # free memory

    ## run in parallel
    results = []
    with ProcessPoolExecutor(max_workers=max_cores) as exe:
        for res in exe.map(_worker_hsm, job_args):
            results.append(res)

    # free shared memory
    shm_img.close(); shm_img.unlink()
    shm_seg.close(); shm_seg.unlink()

    # build output catalogue
    if shear_est in ['REGAUSS', 'LINEAR', 'BJ']:
        e_tag = 'distortion_HSM'
    else:  # KSB
        e_tag = 'shear_HSM'
    shape_cata = pd.DataFrame({
        'id_detec': id_detec,
        f'e1_{e_tag}': -999., f'e2_{e_tag}': -999.,
        'resolution_factor_HSM': -999.,
        'psf_sigma_HSM': -999.,
        'psf_e1_shear_HSM': -999., 'psf_e2_shear_HSM': -999.,
        'moments_sigma_HSM': -999.,
        'moments_amp_HSM': -999.,
        'weight_HSM': weights
    })
    del id_detec, weights  # free memory

    ## results from workers
    for idx, e1, e2, resolution, psf_sigma, psf_g1, psf_g2, sigma, amp in results:
        shape_cata.loc[idx, [
            f'e1_{e_tag}', f'e2_{e_tag}',
            'resolution_factor_HSM',
            'psf_sigma_HSM',
            'psf_e1_shear_HSM',
            'psf_e2_shear_HSM',
            'moments_sigma_HSM',
            'moments_amp_HSM']] = \
            [e1, e2, resolution, psf_sigma, psf_g1, psf_g2, sigma, amp]
    ## set zero weight for unphysical measurements
    shape_cata.loc[(shape_cata[f'e1_{e_tag}']<-1)|(shape_cata[f'e1_{e_tag}']>1)|(shape_cata[f'e2_{e_tag}']<-1)|(shape_cata[f'e2_{e_tag}']>1), 'weight_HSM'] = 0.0

    # save
    tmp = outpath_feather + '_tmp'
    shape_cata.to_feather(tmp)
    del shape_cata
    os.rename(tmp, outpath_feather)
    logger.info(f'Final catalogue saved as {outpath_feather}')
    return 0