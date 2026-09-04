# -*- coding: utf-8 -*-
# @Author: lshuns
# @Date:   1969-12-31 16:00:00
# @Last Modified by:   lshuns
# @Last Modified time: 2026-08-31 13:09:44

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
from metadetect import shearpos

logger = logging.getLogger(__name__)

## versions are logged once, from the first MetaDetectShear call, because at
##    import time the logging is usually not configured yet
_VERSIONS_LOGGED = False

def _log_versions():
    global _VERSIONS_LOGGED
    if _VERSIONS_LOGGED:
        return
    logger.info(f'+++ metadetect version: {metadetect.__version__}')
    logger.info(f'+++ sxdes version: {sxdes.__version__}')
    logger.info(f'+++ ngmix version: {ngmix.__version__}')
    _VERSIONS_LOGGED = True

def _check_psf_centred(psf_img, pixel_scale, tol=1.e-2, allow_uncentred=False):
    """
    Check that the PSF stamp is centred where ngmix assumes it to be.

    ngmix's metacal turns the PSF stamp into a GalSim InterpolatedImage, which
    places the profile at the *true centre* of the stamp, i.e. the zero-based
    position ((nrow-1)/2, (ncol-1)/2). The centre stored in the PSF Jacobian is
    ignored (ngmix.Jacobian.get_galsim_wcs() returns a purely local
    galsim.JacobianWCS). metadetect's own examples follow that same convention,
    see shear_meas_test/test_shear_meas.py (psf_dim=53, psf_cen=(psf_dim-1)/2)
    and metadetect/tests/sim.py (cen = (psf_im.shape-1)/2).

    A stamp centred anywhere else displaces every measured position by the same
    amount, which is why modules/ImSimPSF.py writes two flavours of every stamp:
    the default one is shifted by half a pixel for lensfit, and the '_centred'
    one sits on the true centre. metadetect needs the latter, which is what
    psf_image in the [metadetect] config section selects by default. Handing it
    the shifted one costs 0.5 pixel per axis, i.e. 0.71 pixel of astrometry.

    The measured cost is astrometric. Feeding the half-pixel-shifted stamp to
    metadetect displaces every reported position by 0.70 pixel (0.14 arcsec at
    0.2 arcsec/pixel), against 0.002 pixel with a centred stamp, which is enough
    to wreck the cross-match against the input catalogue.

    No shear bias has been demonstrated. A small ring-test suggested a paired
    difference of -2.6e-3 +- 1.3e-3 in m, but that was only 2.1 sigma, and a much
    larger production run (4 shear pairs, 50 tiles, 30 million objects) put the
    paired difference at +0.0e-3 +- 2.3e-3 in m1 and -0.7e-3 +- 2.2e-3 in m2,
    i.e. consistent with zero and with the earlier hint. This is expected: an
    off-centre PSF acts as a near-pure translation, and a linear shear acts on
    an object's shape regardless of where it sits, the weight function follows
    the detected centroid, and the reconvolution kernel is unchanged. Treat any
    residual m effect as bounded by a few 1e-3 rather than as established.

    This still raises rather than warns, because the astrometric damage alone
    invalidates the cross-match, unless the caller wants the uncentred stamp on
    purpose.

    The stamp is never modified.

    Parameters
    ----------
    psf_img : array
        The PSF stamp.
    pixel_scale (arcsec) : float
        Pixel scale of the PSF stamp.
    tol (pixel) : float, optional
        Offset above which the stamp counts as uncentred.
    allow_uncentred : bool, optional (default: False)
        Warn instead of raising. Only for deliberately reproducing a run made
        with the other stamp, never for production measurements.

    Returns
    -------
    psf_img : numpy.ndarray
        The stamp, in native byte order and unchanged otherwise. It is passed to
        every worker, so its precision is left alone: ngmix converts to float64
        internally, and promoting a float32 FITS stamp here would only double
        what is pickled per cell.

    Raises
    ------
    ValueError
        If the stamp is off-centre by more than tol and allow_uncentred is False.
    """

    ## native byte order (astropy hands back big-endian from a FITS file),
    ##    without touching the precision
    arr = np.ascontiguousarray(psf_img)
    if not arr.dtype.isnative:
        arr = arr.astype(arr.dtype.newbyteorder('='))
    nrow, ncol = arr.shape

    ## Measure the centroid, preferring adaptive moments over the first moment,
    ##    which is sensitive to the truncated wings of the stamp
    try:
        mom = galsim.hsm.FindAdaptiveMom(galsim.Image(arr, scale=pixel_scale))
        ## GalSim image positions are one-based
        row_obs = mom.moments_centroid.y - 1.
        col_obs = mom.moments_centroid.x - 1.
    except Exception as e:
        logger.debug(f'FindAdaptiveMom failed on the PSF stamp ({e}), '
                     'falling back on the flux-weighted centroid')
        rows, cols = np.mgrid[0:nrow, 0:ncol]
        total = arr.sum(dtype=np.float64)
        row_obs = float((arr * rows).sum(dtype=np.float64) / total)
        col_obs = float((arr * cols).sum(dtype=np.float64) / total)

    drow = row_obs - (nrow - 1) / 2.
    dcol = col_obs - (ncol - 1) / 2.

    if (abs(drow) > tol) or (abs(dcol) > tol):
        msg = (f'The PSF stamp is not centred where ngmix assumes it to be: its '
               f'centroid is at ({drow:+.4f}, {dcol:+.4f}) pixel from the true centre '
               f'of the stamp. Every measured position would be displaced by '
               f'{np.hypot(drow, dcol):.2f} pixel '
               f'({np.hypot(drow, dcol)*pixel_scale:.3f} arcsec), and the multiplicative '
               f'shear bias would be off by a few 1e-3. '
               "Use the '_centred' PSF image, i.e. psf_image = centred in the "
               '[metadetect] section of the configuration file.')
        if not allow_uncentred:
            raise ValueError(msg)
        logger.warning(msg + ' Continuing anyway because an uncentred PSF was '
                             'explicitly asked for; do not use this for a '
                             'production shear measurement.')
    else:
        logger.debug(f'PSF stamp is centred (offset {drow:+.4f}, {dcol:+.4f} pixel)')

    return arr

def _cal_sky_cell(rows, cols, shear_str, jac, dims, cell_wcs, step=None):
    """
    Sky coordinates (RA, Dec) for objects detected within a cell.

    metadetect detects and measures objects on the metacalibration-sheared
    images, so for the sheared types ('1p', '1m', '2p', '2m') the detected
    positions are displaced by the artificial shear. That displacement is undone
    here first (following metadetect/shearpos.py), so that the returned sky
    positions are those of the undistorted sky and can be directly cross-matched
    against the input catalogue with modules/CrossMatch.py.

    Parameters
    ----------
    rows, cols : array-like
        Detected positions within the cell, as returned by metadetect in
        'sx_row' and 'sx_col', i.e. zero-based indices of the cell image.
    shear_str : str
        The metadetect shear type: 'noshear', '1p', '1m', '2p', '2m',
        '1p_psf', ... No correction is applied for the types listed in
        shearpos.SKIP_SHEARS ('noshear' and the '*_psf' ones), for which the
        positions are already undistorted.
    jac : ngmix.Jacobian
        The local Jacobian of the cell, describing the WCS used for unshearing.
    dims : (nrow, ncol)
        Shape of the cell image, used to locate the canonical image centre
        about which the artificial shear was applied.
    cell_wcs : galsim.BaseWCS
        The WCS of the cell, with its first pixel at GalSim position (1, 1).
    step : float, optional
        The metacal shear step. Defaults to metadetect's DEFAULT_STEP (0.01).

    Returns
    -------
    RA, DEC : numpy.ndarray
        Sky coordinates in degrees.
    """

    if step is None:
        step = shearpos.DEFAULT_STEP

    rows = np.atleast_1d(np.asarray(rows, dtype=float))
    cols = np.atleast_1d(np.asarray(cols, dtype=float))

    ## Undo the position shift caused by the artificial metacal shear
    ##    (a no-op for 'noshear' and the '*_psf' types)
    rows, cols = shearpos.unshear_positions(rows, cols, shear_str,
                                            jac=jac, dims=dims, step=step)

    ## GalSim images are one-based, the metadetect positions are zero-based
    RA, DEC = cell_wcs.toWorld(np.asarray(cols, dtype=float) + 1.,
                               np.asarray(rows, dtype=float) + 1.,
                               units=galsim.degrees)

    return RA, DEC

def _run_metadetect_cell(args):
    """
    Run metadetect for each cell
    """
    (xcen, ycen, full_wcs, psf_img, cfg, mdet_seed_seq,
     shm_name_img, img_shape, img_dtype,
     shm_name_weight_img, img_weight_shape, img_weight_dtype,
     cell_size, central_size, pixel_scale,
     outpath_cell, outpath_cell_cata,
     shm_name_noise_img, img_noise_shape, img_noise_dtype
     ) = args

    half_size = cell_size // 2
    xmin = xcen - half_size
    ymin = ycen - half_size

    ## Skip incomplete cells
    ##    the grid built by MetaDetectShear always fits, so this is only a guard
    if (ymin < 0) or (ymin + cell_size > img_shape[0]) or \
       (xmin < 0) or (xmin + cell_size > img_shape[1]):
        logger.warning(f'Skipping incomplete cell at xcen={xcen}, ycen={ycen}')
        return pd.DataFrame()

    ## Attach the shared arrays, copy the cell out, and release them again
    ##    NOTE: SharedMemory.close() unmaps the buffer *silently*, so every view
    ##          into it has to be gone before closing, otherwise a later access
    ##          segfaults instead of raising
    shm_img = shared_memory.SharedMemory(name=shm_name_img)
    shm_weight_img = shared_memory.SharedMemory(name=shm_name_weight_img)
    if shm_name_noise_img is not None:
        shm_noise_img = shared_memory.SharedMemory(name=shm_name_noise_img)
    else:
        shm_noise_img = None
    try:
        image_data = np.ndarray(img_shape,
                                dtype=np.dtype(img_dtype),
                                buffer=shm_img.buf)
        weight_data = np.ndarray(img_weight_shape,
                                 dtype=np.dtype(img_weight_dtype),
                                 buffer=shm_weight_img.buf)
        ## copies, in native byte order and the precision ngmix works in
        cell_image = np.ascontiguousarray(
            image_data[ymin:ymin+cell_size, xmin:xmin+cell_size], dtype=np.float64)
        cell_weight_image = np.ascontiguousarray(
            weight_data[ymin:ymin+cell_size, xmin:xmin+cell_size], dtype=np.float64)
        del image_data, weight_data

        if shm_noise_img is not None:
            noise_data = np.ndarray(img_noise_shape,
                                    dtype=np.dtype(img_noise_dtype),
                                    buffer=shm_noise_img.buf)
            cell_noise_image = np.ascontiguousarray(
                noise_data[ymin:ymin+cell_size, xmin:xmin+cell_size], dtype=np.float64)
            del noise_data
        else:
            cell_noise_image = None
    finally:
        shm_img.close()
        shm_weight_img.close()
        if shm_noise_img is not None:
            shm_noise_img.close()

    ## The corresponding WCS for the cell
    ## NOTE: wcs.shiftOrigin(o).toWorld(p) == wcs.toWorld(p - o), so mapping the
    ##        cell pixel (1, 1) onto the full-image pixel (xmin+1, ymin+1)
    ##        requires a shift of -xmin, -ymin
    cell_wcs = full_wcs.shiftOrigin(galsim.PositionI(-xmin, -ymin))
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
            cell_image.astype(np.float32),
            wcs=cell_wcs,
            xmin=1,
            ymin=1
        )
        img.write(outpath_cell, clobber=True)
        logger.debug(f'Saved cell image for xcen={xcen}, ycen={ycen}')
        del img
        ## Save weight image
        img = galsim.Image(
            cell_weight_image.astype(np.float32),
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
                cell_noise_image.astype(np.float32),
                wcs=cell_wcs,
                xmin=1,
                ymin=1
            )
            img.write(outpath_cell+'.noise.fits', clobber=True)
            logger.debug(f'Saved cell noise image for xcen={xcen}, ycen={ycen}')
            del img

    ## Prepare ngmix mbobs and run metadetect
    ##    everything that can legitimately fail on a bad cell (a fully masked
    ##    cell already raises in ngmix.Observation) is kept inside the try, so
    ##    that one bad cell costs a cell and not the whole image
    try:
        psf_jac = ngmix.DiagonalJacobian(scale=pixel_scale, 
                                        row=(psf_img.shape[0]-1)/2, 
                                        col=(psf_img.shape[1]-1)/2)
        obs = ngmix.Observation(
            image=cell_image,
            weight=cell_weight_image,
            jacobian=jac,
            ormask=np.zeros_like(cell_image, dtype=np.int32),
            bmask=np.zeros_like(cell_image, dtype=np.int32),
            psf=ngmix.Observation(
                image=psf_img.copy(),
                jacobian=psf_jac,
            ),
        )
        if cell_noise_image is not None:
            obs.noise = cell_noise_image
        mbobs = ngmix.MultiBandObsList()
        obslist = ngmix.ObsList()
        obslist.append(obs)
        mbobs.append(obslist)

        res = metadetect.do_metadetect(
                copy.deepcopy(cfg),
                mbobs,
                np.random.RandomState(np.random.MT19937(mdet_seed_seq)))
    except Exception as e:
        logger.error(f"Metadetect failed for cell at ({xcen}, {ycen}): {e}")
        return pd.DataFrame()

    ## metadetect returns None when nothing can be measured at all, e.g. the
    ##    cell is entirely masked or the metacal step itself failed
    if res is None:
        logger.warning(f"Metadetect returned no result for cell at ({xcen}, {ycen})")
        return pd.DataFrame()

    ## Save cell results if required
    ##    types without any detection are None, drop them so that the npz stays
    ##    readable without allow_pickle
    if outpath_cell_cata is not None:
        res_save = {key: val for key, val in res.items() if val is not None}
        if len(res_save) < len(res):
            logger.debug(f'No detections for {sorted(set(res) - set(res_save))} '
                         f'in cell at ({xcen}, {ycen}), not saved')
        np.savez_compressed(outpath_cell_cata, **res_save)
        del res_save

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

        ## Sky positions, with the shift from the artificial shear removed,
        ##    for cross-matching with the input catalogue (see CrossMatch.py)
        RA, DEC = _cal_sky_cell(
                    res[key]['sx_row'][idx_centre],
                    res[key]['sx_col'][idx_centre],
                    key, jac, (cell_size, cell_size), cell_wcs,
                    step=cfg.get('metacal', {}).get('step', shearpos.DEFAULT_STEP))
        res_key['X_WORLD'] = RA
        res_key['Y_WORLD'] = DEC
        del RA, DEC

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
                    max_cores=12,
                    allow_uncentred_psf=False):
    """
    Main function for running metadetect

    allow_uncentred_psf : bool, optional (default: False)
        Proceed, with a warning, when the PSF stamp is not centred where ngmix
        assumes it to be. Off by default because such a run is biased in both
        position and shear, see _check_psf_centred.
    """

    _log_versions()

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
    ## ngmix places the PSF at the true centre of its stamp, refuse to run if it
    ##    is not there (see _check_psf_centred)
    psf_img = _check_psf_centred(psf_img, pixel_scale,
                                 allow_uncentred=allow_uncentred_psf)

    ## Main image
    with fits.open(inpath_image) as hdul:
        image_data = hdul[0].data
        full_wcs = galsim.FitsWCS(header=hdul[0].header)
    logger.debug(f'Image loaded from {inpath_image}') 
    Nimg_y, Nimg_x = image_data.shape
    logger.info(f'Image size: {Nimg_y} x {Nimg_x}')
    if (Nimg_y < central_size) or (Nimg_x < central_size):
        logger.warning(f'Image size {Nimg_y}x{Nimg_x} is smaller than '
                       f'central_size {central_size}, it will be padded into a single cell')

    ## The PSF is described by a plain pixel scale, so it has to match the scale
    ##    of the image WCS, otherwise the PSF has the wrong size relative to the
    ##    galaxies and the shear calibration is wrong
    jac_check = full_wcs.jacobian(galsim.PositionD(Nimg_x/2., Nimg_y/2.))
    wcs_scale = np.sqrt(abs(jac_check.dudx*jac_check.dvdy
                            - jac_check.dudy*jac_check.dvdx))
    if not np.isclose(wcs_scale, pixel_scale, rtol=1.e-3):
        logger.warning(f'pixel_scale={pixel_scale} does not match the scale of the image WCS '
                       f'({wcs_scale:.6g} arcsec/pixel)! The PSF will have the wrong size '
                       'relative to the image.')
    del jac_check

    ## Noise map if provided
    if inpath_noise_map is not None:
        with fits.open(inpath_noise_map) as hdul:
            noise_data = hdul[0].data
        logger.debug(f'Noise image loaded from {inpath_noise_map}')

    ## Weight map
    ##    metadetect sets the detection threshold from the weight map
    ##    (detnoise = 1/sqrt(median(weight)), metadetect/detect.py), and ngmix
    ##    uses it as the inverse variance for s2n and all parameter errors.
    ##    A flat weight of one would silently scale the detection threshold and
    ##    every s2n/error by the true background rms, so it has to be a real
    ##    inverse variance.
    if inpath_weight_map is not None:
        with fits.open(inpath_weight_map) as hdul:
            weight_data = hdul[0].data
        logger.debug(f'Weight image loaded from {inpath_weight_map}')
    else:
        if inpath_noise_map is not None:
            noise_var = float(np.var(noise_data))
            logger.info('No weight map provided, using the inverse variance of the noise image '
                        f'(rms = {np.sqrt(noise_var):.6g})')
        else:
            ## last resort: a robust estimate straight off the image
            med_tmp = float(np.median(image_data))
            mad_tmp = float(np.median(np.abs(image_data - med_tmp)))
            noise_var = (1.4826 * mad_tmp)**2.
            logger.warning('Neither a weight map nor a noise map was provided! '
                           f'Estimating the background rms from the image itself ({np.sqrt(noise_var):.6g}). '
                           'Provide inpath_noise_map or inpath_weight_map for a reliable '
                           'detection threshold and reliable errors.')
            del med_tmp, mad_tmp
        if not (noise_var > 0):
            raise ValueError(f'Non-positive noise variance ({noise_var}), cannot build a weight map!')
        weight_data = np.full(image_data.shape, 1./noise_var, dtype=np.float32)

    ## >>>>>>>>>>>>> 1.5 Pad the images so that the cells tile them completely
    ##    the buffer around the central region of a cell
    pad = (cell_size - central_size) // 2
    ##    plus whatever is needed to make the image a whole number of central
    ##    regions, so that no strip of the image is left untiled
    Ny = int(np.ceil(Nimg_y / central_size))
    Nx = int(np.ceil(Nimg_x / central_size))
    pad_after_y = pad + Ny*central_size - Nimg_y
    pad_after_x = pad + Nx*central_size - Nimg_x
    pad_width = ((pad, pad_after_y), (pad, pad_after_x))

    image_data = np.pad(image_data, 
                        pad_width=pad_width, 
                        mode='constant', 
                        constant_values=0)
    ## Update the size
    Nimg_y, Nimg_x = image_data.shape
    ## Adjust WCS for padding
    ## the padding moves the original pixel (1, 1) to (1+pad, 1+pad),
    ##    and wcs.shiftOrigin(o).toWorld(p) == wcs.toWorld(p - o)
    full_wcs = full_wcs.shiftOrigin(galsim.PositionI(pad, pad))
    ## Shared memory for the big arrays
    ##    everything that touches them is wrapped, because a segment that is not
    ##    unlinked survives in /dev/shm until reboot: an image-sized leak on
    ##    every failed tile would fill it up
    shm_img = shm_weight_img = shm_noise_img = None
    shm_img_arr = shm_weight_img_arr = shm_noise_img_arr = None
    try:
        shm_img = shared_memory.SharedMemory(create=True, 
                                             size=image_data.nbytes)
        shm_img_arr = np.ndarray(image_data.shape, 
                                 dtype=image_data.dtype, 
                                 buffer=shm_img.buf)
        shm_img_arr[:] = image_data
        del image_data

        ## Pad the weight image
        ##    zero weight outside the real image, so that the padding is ignored
        weight_data = np.pad(weight_data, 
                            pad_width=pad_width, 
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
                                pad_width=pad_width, 
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

        ## >>>>>>>>>>>>> 1.6 Pre-compile numba functions to avoid fork issues
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
        step_size = central_size
        Ncells = Ny * Nx
        logger.info(f'Image size (padded): {Nimg_y} x {Nimg_x}')
        logger.info(f'Cell size: {cell_size}, Central size: {central_size}, Step size: {step_size}')
        logger.info(f'Grid: {Ny} x {Nx} = {Ncells} cells')

        ## Independent random streams for the cells
        ##    consecutive integer seeds would make the streams of neighbouring
        ##    images overlap, so spawn them from a SeedSequence instead
        seed_seq_cells, seed_seq_save = np.random.SeedSequence(random_seed).spawn(2)
        cell_seed_seqs = seed_seq_cells.spawn(Ncells)

        ## Randomly pick objects to save cells 
        if save_Ncells > 0: 
            rng_cells = np.random.RandomState(np.random.MT19937(seed_seq_save))
            saved_cell_indices = rng_cells.choice(Ncells,
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
                job_arg = [xcen, ycen, full_wcs, psf_img, cfg, cell_seed_seqs[iy * Nx + ix],
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
        with ProcessPoolExecutor(max_workers=max_cores) as exe:
            for res in exe.map(_run_metadetect_cell, job_args):
                results.append(res)
    finally:
        ## Free shared memory
        ##    drop every view first: close() unmaps the buffer silently, and a
        ##    later access to a dangling view segfaults instead of raising
        del shm_img_arr, shm_weight_img_arr, shm_noise_img_arr
        for shm_tmp in (shm_img, shm_weight_img, shm_noise_img):
            if shm_tmp is not None:
                shm_tmp.close()
                shm_tmp.unlink()

    ## Save the final catalogue
    results = [r for r in results if not r.empty]
    if not results:
        raise RuntimeError("No objects detected in any cell!")
    else:
        results = pd.concat(results, ignore_index=True)
    ## A unique id for every row, following the NUMBER convention of the
    ##    SExtractor detection catalogues. Note that metadetect reports one row
    ##    per (object, shear_type), so NUMBER identifies a row, not a sky object;
    ##    that is what the cross-match and task 7 key on.
    results.insert(0, 'NUMBER', np.arange(len(results), dtype=int))
    tmp = outpath_feather + '_tmp'
    results.to_feather(tmp)
    del results
    os.rename(tmp, outpath_feather)
    logger.info(f'Final catalogue saved as {outpath_feather}')
    return 0
