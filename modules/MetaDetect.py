# -*- coding: utf-8 -*-
# @Author: lshuns
# @Date:   1969-12-31 16:00:00
# @Last Modified by:   lshuns
# @Last Modified time: 2026-01-20 16:23:31

## Shear measurement using metadetection

import os
import sys
import copy
import json
import logging

import ngmix
import galsim
import numpy as np
import pandas as pd
from astropy.io import fits
from multiprocessing import shared_memory
from concurrent.futures import ProcessPoolExecutor

## Downloaded external modules 
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "sxdes"))
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "metadetect"))
import sxdes
import metadetect

logger = logging.getLogger(__name__)
logger.info(f'+++ metadetect version: {metadetect.__version__}')
logger.info(f'+++ sxdes version: {sxdes.__version__}')
logger.info(f'+++ ngmix version: {ngmix.__version__}')

def _run_metadetect_cell(args):
    """
    Run metadetect for each cell
    """
    (xcen, ycen, full_wcs, psf_img, cfg, mdet_seed,
     shm_name_img, img_shape, img_dtype,
     shm_name_weight_img, img_weight_shape, img_weight_dtype,
     cell_size, central_size, pixel_scale,
     outpath_cell, outpath_cell_cata,
     shm_name_noise_img, img_noise_shape, img_noise_dtype
     ) = args

    ## Re-attach shared arrays
    shm_img = shared_memory.SharedMemory(name=shm_name_img)
    image_data = np.ndarray(img_shape, 
                            dtype=np.dtype(img_dtype), 
                            buffer=shm_img.buf)
    shm_weight_img = shared_memory.SharedMemory(name=shm_name_weight_img)
    weight_data = np.ndarray(img_weight_shape, 
                                dtype=np.dtype(img_weight_dtype), 
                                buffer=shm_weight_img.buf)
    if shm_name_noise_img is not None:
        shm_noise_img = shared_memory.SharedMemory(name=shm_name_noise_img)
        noise_data = np.ndarray(img_noise_shape, 
                                dtype=np.dtype(img_noise_dtype), 
                                buffer=shm_noise_img.buf)
    else:
        noise_data = None

    ## Skip incomplete cells
    half_size = cell_size // 2
    if (ycen-half_size < 0) or (ycen+half_size > image_data.shape[0]) or \
    (xcen-half_size < 0) or (xcen+half_size > image_data.shape[1]):
        logger.debug(f'Skipping incomplete cell at xcen={xcen}, ycen={ycen}')
        return pd.DataFrame()

    ## Extract cells
    cell_image = image_data[ycen-half_size:ycen+half_size,
                             xcen-half_size:xcen+half_size]
    cell_weight_image = weight_data[ycen-half_size:ycen+half_size,
                            xcen-half_size:xcen+half_size]
    if noise_data is not None:
        cell_noise_image = noise_data[ycen-half_size:ycen+half_size,
                                xcen-half_size:xcen+half_size]
    else:
        cell_noise_image = None

    ## The corresponding WCS for the cell
    xmin = xcen - half_size
    ymin = ycen - half_size
    cell_wcs = full_wcs.shiftOrigin(galsim.PositionI(xmin, ymin))
    ## The local Jacobian at the centre of the cell
    jac_gs = cell_wcs.jacobian(galsim.PositionD(half_size + 0.5, 
                                                half_size + 0.5))
    ## Convert to ngmix Jacobian
    jac = ngmix.Jacobian(
        row=half_size-0.5,
        col=half_size-0.5,
        dudrow=jac_gs.dudy,
        dudcol=jac_gs.dudx,
        dvdrow=jac_gs.dvdy,
        dvdcol=jac_gs.dvdx,
    )
    del jac_gs

    ## Save cell image if required
    if outpath_cell is not None:
        img = galsim.Image(
            cell_image,
            wcs=cell_wcs,
            xmin=1,
            ymin=1
        )
        img.write(outpath_cell, clobber=True)
        logger.debug(f'Saved cell image for xcen={xcen}, ycen={ycen}')
        del img
        ## Save weight image
        img = galsim.Image(
            cell_weight_image,
            wcs=cell_wcs,
            xmin=1,
            ymin=1
        )
        img.write(outpath_cell+'.weight.fits', clobber=True)
        logger.debug(f'Saved cell weight image for xcen={xcen}, ycen={ycen}')
        del img
        ## Save noise image
        if cell_noise_image is not None:
            img = galsim.Image(
                cell_noise_image,
                wcs=cell_wcs,
                xmin=1,
                ymin=1
            )
            img.write(outpath_cell+'.noise.fits', clobber=True)
            logger.debug(f'Saved cell noise image for xcen={xcen}, ycen={ycen}')
            del img

    ## Prepare ngmix mbobs
    psf_jac = ngmix.DiagonalJacobian(scale=pixel_scale, 
                                    row=(psf_img.shape[0]-1)/2, 
                                    col=(psf_img.shape[1]-1)/2)
    obs = ngmix.Observation(
        image=cell_image.copy(),
        weight=cell_weight_image.copy(),
        jacobian=jac,
        ormask=np.zeros_like(cell_image, dtype=np.int32),
        bmask=np.zeros_like(cell_image, dtype=np.int32),
        psf=ngmix.Observation(
            image=psf_img.copy(),
            jacobian=psf_jac,
        ),
    )
    if cell_noise_image is not None:
        obs.noise = cell_noise_image.copy()
    mbobs = ngmix.MultiBandObsList()
    obslist = ngmix.ObsList()
    obslist.append(obs)
    mbobs.append(obslist)

    ## Run metadetect
    try:
        res = metadetect.do_metadetect(
                copy.deepcopy(cfg),
                mbobs,
                np.random.RandomState(seed=mdet_seed))
    except Exception as e:
        logger.error(f"Metadetect failed for cell at ({xcen}, {ycen}): {e}")
        return pd.DataFrame()
    finally:
        # free shared memory
        shm_img.close()
        if weight_data is not None:
            shm_weight_img.close()
        if noise_data is not None:
            shm_noise_img.close()

    ## Save cell results if required
    if outpath_cell_cata is not None:
        np.savez_compressed(outpath_cell_cata, **res)

    ## Select and return results for objects in the centre region
    res_centre_list = []
    for key in res.keys():

        ## Check if cell results are None
        if res[key] is None:
            logger.debug(f"No objects detected in {key} for cell at ({xcen}, {ycen})")
            continue

        ## Use unsheared position to select objects within the centre region
        idx_centre = np.where(
            (res[key]['sx_row_noshear'] >= (cell_size - central_size) / 2) &
            (res[key]['sx_row_noshear'] < (cell_size + central_size) / 2) &
            (res[key]['sx_col_noshear'] >= (cell_size - central_size) / 2) &
            (res[key]['sx_col_noshear'] < (cell_size + central_size) / 2)
        )[0]
        logger.debug(f'For {key}, number of objects in centre region / all: {len(idx_centre)}/{len(res[key])}')

        ## Skip if no objects in central region
        if len(idx_centre) == 0:
            logger.debug(f"No objects in central region in {key} for cell at ({xcen}, {ycen})")
            continue
        
        ## Create a dictionary for this key
        res_key = {'shear_type': [key] * len(idx_centre)}  

        ## Unpack the results
        for col in res[key].dtype.names:
            if len(res[key][col].shape) == 1:
                res_key[f'{col}'] = res[key][col][idx_centre]
            elif len(res[key][col].shape) == 2:
                res_key[f'{col}_1'] = res[key][col][idx_centre, 0]
                res_key[f'{col}_2'] = res[key][col][idx_centre, 1]
            ## Save cov as scalar uncertainties
            elif len(res[key][col].shape) == 3:
                res_key[f'{col}_as_sigma'] = np.sqrt(
                    np.trace(res[key][col][idx_centre], 
                             axis1=1, axis2=2) /2. )
        
        ## Convert to DataFrame and add to list
        df_key = pd.DataFrame(res_key)
        res_centre_list.append(df_key)
        del res_key, df_key

    ## Concatenate all DataFrames along rows
    if res_centre_list:
        return pd.concat(res_centre_list, axis=0, ignore_index=True)
    else:
        logger.debug(f"No objects for cell at ({xcen}, {ycen})")
        return pd.DataFrame()

def MetaDetectShear(outpath_feather,
                    inpath_image, 
                    inpath_psf,
                    inpath_config,
                    pixel_scale,
                    inpath_weight_map=None, 
                    inpath_noise_map=None,
                    save_Ncells=0,
                    random_seed=914,
                    cell_size=250,
                    central_size=150,
                    max_cores=12):
    """
    Main function for running metadetect
    """

    logger.info(f'Run metadetection for {os.path.basename(outpath_feather)}...')
    if os.path.isfile(outpath_feather):
        logger.info(f'The final feather catalogue {outpath_feather} already exists.')
        logger.info(f'End the process.')
        return 1

    ## Some sanity checks
    assert central_size < cell_size, "central_size should be smaller than cell_size!"
    assert cell_size % 2 == 0, "cell_size should be an even number!"
    assert central_size % 2 == 0, "central_size should be an even number!"

    ## >>>>>>>>>>> 0. Load config
    with open(inpath_config, 'r') as json_file:
        cfg = json.load(json_file)

    ## Check necessary images for metacal noise fixing
    if cfg['metacal']['fixnoise']:
        if cfg['metacal']['use_noise_image']:   
            assert inpath_noise_map is not None, "inpath_noise_map should be provided for noise fixing with use_noise_image=True"
        else:
            assert inpath_weight_map is not None, "inpath_weight_map should be provided for noise fixing with use_noise_image=False"

    ## >>>>>>>>>>>>> 1. Load images
    ## PSF image
    with fits.open(inpath_psf) as hdul:
        psf_img = hdul[0].data
    logger.debug(f'PSF image loaded from {inpath_psf}') 

    ## Main image
    with fits.open(inpath_image) as hdul:
        image_data = hdul[0].data
        full_wcs = galsim.FitsWCS(header=hdul[0].header)
    logger.debug(f'Image loaded from {inpath_image}') 
    ## Ensure the image is large enough for at least one cell
    Nimg_y, Nimg_x = image_data.shape
    if Nimg_y < cell_size or Nimg_x < cell_size:
        raise ValueError(f"Image size (padded) {Nimg_y}x{Nimg_x} is smaller than "
                        f"cell_size {cell_size}!")

    ## Weight map if provided
    if inpath_weight_map is not None:
        with fits.open(inpath_weight_map) as hdul:
            weight_data = hdul[0].data
        logger.debug(f'Weight image loaded from {inpath_weight_map}')
    else:
        weight_data = np.ones_like(image_data)

    ## Noise map if provided
    if inpath_noise_map is not None:
        with fits.open(inpath_noise_map) as hdul:
            noise_data = hdul[0].data
        logger.debug(f'Noise image loaded from {inpath_noise_map}')

    ## Pad the image to fit cells
    pad = (cell_size - central_size) // 2
    image_data = np.pad(image_data, 
                        pad_width=pad, 
                        mode='constant', 
                        constant_values=0)
    ## Update the size
    Nimg_y, Nimg_x = image_data.shape
    ## Adjust WCS for padding
    full_wcs = full_wcs.shiftOrigin(galsim.PositionI(-pad, -pad))
    ## Shared memory for big arrays
    shm_img = shared_memory.SharedMemory(create=True, 
                                         size=image_data.nbytes)
    shm_img_arr = np.ndarray(image_data.shape, 
                             dtype=image_data.dtype, 
                             buffer=shm_img.buf)
    shm_img_arr[:] = image_data
    del image_data

    ## Pad the weight image
    weight_data = np.pad(weight_data, 
                        pad_width=pad, 
                        mode='constant', 
                        constant_values=0)
    ## Shared memory for big arrays
    shm_weight_img = shared_memory.SharedMemory(create=True, 
                                                size=weight_data.nbytes)
    shm_weight_img_arr = np.ndarray(weight_data.shape, 
                            dtype=weight_data.dtype, 
                            buffer=shm_weight_img.buf)
    shm_weight_img_arr[:] = weight_data
    del weight_data

    ## Pad the noise image
    if inpath_noise_map is not None:
        noise_data = np.pad(noise_data, 
                            pad_width=pad, 
                            mode='edge')
        ## Shared memory for big arrays
        shm_noise_img = shared_memory.SharedMemory(create=True, 
                                                   size=noise_data.nbytes)
        shm_noise_img_arr = np.ndarray(noise_data.shape, 
                                dtype=noise_data.dtype, 
                                buffer=shm_noise_img.buf)
        shm_noise_img_arr[:] = noise_data
        del noise_data
    else:
        shm_noise_img = None

    ## >>>>>>>>>>>>> 1.5 Pre-compile numba functions to avoid fork issues
    logger.info("Pre-compiling numba functions...")
    ## Create a dummy observation to trigger JIT compilation
    dummy_img = np.ones((10, 10), dtype=np.float64)
    dummy_jac = ngmix.DiagonalJacobian(scale=pixel_scale, row=4.5, col=4.5)
    dummy_obs = ngmix.Observation(
        image=dummy_img,
        weight=np.ones_like(dummy_img),
        jacobian=dummy_jac,
    )
    ## This triggers the fill_pixels compilation
    _ = dummy_obs.pixels
    logger.info("Numba compilation complete")
    del dummy_img, dummy_jac, dummy_obs

    ## >>>>>>>>>>>>> 2. Run metadetection for each cell
    ## Calculate how many cells needed to tile the PADDED image
    step_size = central_size
    Ny = (Nimg_y - cell_size) // step_size + 1
    Nx = (Nimg_x - cell_size) // step_size + 1
    Ncells = Ny * Nx
    logger.info(f'Image size (padded): {Nimg_y} x {Nimg_x}')
    logger.info(f'Cell size: {cell_size}, Central size: {central_size}, Step size: {step_size}')
    logger.info(f'Grid: {Ny} x {Nx} = {Ncells} cells')
    ## Randomly pick objects to save cells 
    if save_Ncells > 0: 
        np.random.seed(random_seed + save_Ncells) 
        saved_cell_indices = np.random.choice(Ncells, 
                                      size=min(save_Ncells, Ncells), 
                                      replace=False)
        cells_dir = os.path.join(os.path.dirname(outpath_feather), 
                                 'cells_metadetect') 
        os.makedirs(cells_dir, exist_ok=True) 
        logger.info(f'Save {save_Ncells} out of {Ncells} cells to {cells_dir}')
        logger.info(f'Chosen cell indices: {sorted(saved_cell_indices)}')

    ## Prepare arguments for workers
    job_args = []
    for iy in range(Ny):
        for ix in range(Nx):
            ycen = iy * step_size + cell_size // 2
            xcen = ix * step_size + cell_size // 2
            if (save_Ncells > 0) and (
                ((iy * Nx + ix) in saved_cell_indices) 
                or ((iy==0) and (ix==0))
                or ((iy==Ny-1) and (ix==Nx-1)) # always save the first and the last
                    ):     
                base = os.path.splitext(os.path.basename(inpath_image))[0]
                cell_fname = f"{base}_cellY{iy}_X{ix}.fits"
                outpath_cell = os.path.join(cells_dir, cell_fname)
                cell_cata_fname = f"{base}_cellY{iy}_X{ix}_metadetect.npz"
                outpath_cell_cata = os.path.join(cells_dir, cell_cata_fname)
            else:
                outpath_cell = None
                outpath_cell_cata = None
            job_arg = [xcen, ycen, full_wcs, psf_img, cfg, random_seed + iy * Nx + ix,
                       shm_img.name, shm_img_arr.shape, str(shm_img_arr.dtype),
                       shm_weight_img.name, shm_weight_img_arr.shape, str(shm_weight_img_arr.dtype),
                       cell_size, central_size, pixel_scale,
                       outpath_cell, outpath_cell_cata
                       ]
            if shm_noise_img is not None:
                job_arg += [shm_noise_img.name, 
                             shm_noise_img_arr.shape, 
                             str(shm_noise_img_arr.dtype)]
            else:
                job_arg += [None, None, None]
            job_args.append(tuple(job_arg))

    ## Run in parallel
    logger.debug(f'Total number of cells to process: {Ncells}')
    results = []
    try:
        with ProcessPoolExecutor(max_workers=max_cores) as exe:
            for res in exe.map(_run_metadetect_cell, job_args):
                results.append(res)
    finally:
        ## Free shared memory
        shm_img.close(); shm_img.unlink()
        if shm_weight_img is not None:
            shm_weight_img.close(); shm_weight_img.unlink()
        if shm_noise_img is not None:
            shm_noise_img.close(); shm_noise_img.unlink()

    ## Save the final catalogue
    results = [r for r in results if not r.empty]
    if not results:
        raise RuntimeError("No objects detected in any cell!")
    else:
        results = pd.concat(results, ignore_index=True)
    tmp = outpath_feather + '_tmp'
    results.to_feather(tmp)
    del results
    os.rename(tmp, outpath_feather)
    logger.info(f'Final catalogue saved as {outpath_feather}')
    return 0
