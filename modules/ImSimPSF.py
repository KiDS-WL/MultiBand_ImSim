# -*- coding: utf-8 -*-
# @Author: lshuns
# @Date:   2020-11-26 15:00:22
# @Last Modified by:   lshuns
# @Last Modified time: 2021-11-24 13:47:16

### Everything about PSF
__all__ = ['MoffatPSF', 'AiryPSF', 'loadPixelPSF', \
            'PSFima', 'PSFmap', 'PSFmap_MultiPSF', 'PSFmap_DiffMag', 'PSFmap_MultiPSF_DiffMag']

import galsim
import logging
import numpy as np

logger = logging.getLogger(__name__)

def MoffatPSF(seeing, moffat_beta, psf_e=None):
    """
    Generate PSF model as a Moffat profile.

    Parameters
    ----------
    seeing : float
        Full-width-half-max of the PSF.
    moffat_beta : float
        Beta parameter for Moffat profile.
    psf_e : list, optional (default: [0, 0])
        List of psf ellipticity [e1, e2] 
        NOTE: e = 1 - q, as from KiDS AW data.

    Returns
    -------
    psf : galsim PSF model.
    """

    psf = galsim.Moffat(beta=moffat_beta, fwhm=seeing, trunc=4.5*seeing)
    if psf_e:
        psf_e1 = psf_e[0]
        psf_e2 = psf_e[1]

        psf_e = np.sqrt(psf_e1**2+psf_e2**2)
        # g_i = e_i/(2-e)
        psf_g1, psf_g2 = psf_e1/(2-psf_e), psf_e2/(2-psf_e)

        psf = psf.shear(g1=psf_g1, g2=psf_g2)

    return psf

def AiryPSF(lam, diam, obscuration, psf_e=None):
    """
    Generate PSF model as an Airy profile.

    Parameters
    ----------
    lam (nanometre) : float
        Wavelength in unit of nanometre.
    diam (metre) : float
        Telescope diameter in unit of metre.
    obscuration : float
        The linear dimension of a central obscuration as a fraction of the pupil dimension.
    psf_e : list, optional (default: [0, 0])
        List of psf ellipticity [e1, e2] 
        NOTE: e = 1 - q, as from KiDS AW data.

    Returns
    -------
    psf : galsim PSF model.
    """

    psf = galsim.Airy(lam=lam, diam=diam, obscuration=obscuration, scale_unit=galsim.arcsec)
    if psf_e:
        psf_e1 = psf_e[0]
        psf_e2 = psf_e[1]

        psf_e = np.sqrt(psf_e1**2+psf_e2**2)
        # g_i = e_i/(2-e)
        psf_g1, psf_g2 = psf_e1/(2-psf_e), psf_e2/(2-psf_e)

        psf = psf.shear(g1=psf_g1, g2=psf_g2)

    return psf

def loadPixelPSF(inpath, pixel_scale, offset=(0.5, 0.5)):
    """
    Load pixelized PSF images as PSF model.

    Parameters
    ----------
    inpath : str
        FITS file includes psf image.
    pixel_scale : float
        pixel scale corresponding to the psf image.
    offset : (float, float), optional (default: [0.5, 0.5])
        offset relative to the center.
        Use the default (0.5, 0.5) value for PSF from lensfit

    Returns
    -------
    PixelPSF : galsim InterpolatedImage.
    """

    PixelPSF = galsim.InterpolatedImage(inpath, scale=pixel_scale, flux=1., offset=offset)

    return PixelPSF

def PSFima(PSF, pixel_scale, size=32, pixelPSF=False):
    """
    Draw a single PSF image from a PSF model.

    Parameters
    ----------
    PSF : galsim object:
        PSF model.
    pixel_scale (arcsec) : float
        Pixel size in unit of arcsec.
    size : int, optional (default: 32)
        The image size in unit of pixel
    pixelPSF : bool, optional (default: False)
        if the PSF provided already including pixel response

    Returns
    -------
    psf_image : galsim Images
        Image for the PSF.
    """

    psf_image = galsim.Image(size, size)
    PSF_lf = PSF.shift(0.5*pixel_scale, 0.5*pixel_scale)

    if pixelPSF:
        draw_method = 'no_pixel'
    else:
        draw_method = 'auto'
    psf_image = PSF_lf.drawImage(image=psf_image, scale=pixel_scale, method=draw_method)

    return psf_image

def PSFmap(PSF, pixel_scale, mag_input, mag_zero=30., N_PSF=100, sep_PSF=120, rng_seed=940120, pixelPSF=False):
    """ 
    Generate N_PSF PSF stars from a single PSF model,
        with uniform magnitude.

    Parameters
    ----------
    PSF : galsim object:
        PSF model.
    pixel_scale (arcsec) : float
        Pixel size in unit of arcsec.
    mag_input : float
        The input magnitude for PSF.
    mag_zero : float, optional (default: 30.)
        The zero point for magnitude.
    N_PSF : int, optional (default: 100)
        Number of PSF generated, better be the square of some integer.
    sep_PSF (pixel) : int, optional (default: 120)
        Seperation between PSFs.        
    rng_seed : int, optional (default: 940120) 
        Seed for random number generator.
    pixelPSF : bool, optional (default: False)
        if the PSF provided already including pixel response

    Returns
    -------
    psf_image : galsim Images
        Image for PSFs.
    """

    # magnitude to flux
    flux = 10**(-0.4*(mag_input-mag_zero))

    # position
    # separation in pixels
    separation = int(sep_PSF)
    ## number in each row
    Nrow = int(N_PSF**0.5)
    ## get grid center
    x = np.arange(separation, separation+Nrow*separation, separation, dtype='int')
    y = np.repeat(x, Nrow)
    x = np.tile(x, Nrow)
    ## check outliers
    Nrow = N_PSF - len(x)
    if Nrow > 0:
        x = np.concatenate([x, np.arange(separation, separation+Nrow*separation, separation)])
        y = np.concatenate([y, np.full(Nrow, y[-1]+separation)])
    elif Nrow < 0:
        x = x[:N_PSF]
        y = y[:N_PSF]
    ## make random shift
    np.random.seed(rng_seed)
    shift_lim = int(separation/5.)
    dx_dy = np.random.randint(low=-shift_lim, high=shift_lim, size=N_PSF)
    x += dx_dy
    y += dx_dy    

    # initiate a canvas
    psf_image = galsim.ImageF(int(np.amax(x)+separation), int(np.max(y)+separation), scale=pixel_scale)
    logger.debug(f"PSF image bounds {psf_image.bounds}")

    # draw PSF 
    star_like_PSF = PSF.withFlux(flux)
    ## draw stamp
    if pixelPSF:
        draw_method = 'no_pixel'
    else:
        draw_method = 'auto'
    stamp_PSF = star_like_PSF.drawImage(scale=pixel_scale, method=draw_method)

    # place PSF
    for i in range(N_PSF):

        # position the stamp
        stamp_PSF.setCenter(x[i], y[i])

        # stitch the postage stamp to the canvas
        overlap = stamp_PSF.bounds & psf_image.bounds
        psf_image[overlap] += stamp_PSF[overlap]

    return psf_image

def PSFmap_MultiPSF(PSF_list, pixel_scale, mag_input, mag_zero=30., sep_PSF=120, rng_seed=940120, pixelPSF=False):
    """ 
    Generate PSF stars from a list of PSF models,
        with uniform magnitude.

    Parameters
    ----------
    PSF_list : a list of galsim object
        PSF models.
    pixel_scale (arcsec) : float
        Pixel size in unit of arcsec.
    mag_input : float, optional (default: 20.0)
        The input magnitude for PSF.
    mag_zero : float, optional (default: 30.)
        The zero point for magnitude.
    sep_PSF (pixel) : int, optional (default: 120)
        Seperation between PSFs.        
    rng_seed : int, optional (default: 940120) 
        Seed for random number generator.
    pixelPSF : bool, optional (default: False)
        if the PSF provided already including pixel response

    Returns
    -------
    psf_image : galsim Images
        Image for PSFs.
    """

    # magnitude to flux
    flux = 10**(-0.4*(mag_input-mag_zero))

    # number of psf
    N_PSF = len(PSF_list)

    # position
    # separation in pixels
    separation = int(sep_PSF)
    ## number in each row
    Nrow = int(N_PSF**0.5)
    ## get grid center
    x = np.arange(separation, separation+Nrow*separation, separation, dtype='int')
    y = np.repeat(x, Nrow)
    x = np.tile(x, Nrow)
    ## check outliers
    Nrow = N_PSF - len(x)
    if Nrow > 0:
        x = np.concatenate([x, np.arange(separation, separation+Nrow*separation, separation)])
        y = np.concatenate([y, np.full(Nrow, y[-1]+separation)])
    elif Nrow < 0:
        x = x[:N_PSF]
        y = y[:N_PSF]
    ## make random shift
    np.random.seed(rng_seed)
    shift_lim = int(separation/5.)
    dx_dy = np.random.randint(low=-shift_lim, high=shift_lim, size=N_PSF)
    x += dx_dy
    y += dx_dy    

    # initiate a canvas
    psf_image = galsim.ImageF(int(np.amax(x)+separation), int(np.max(y)+separation), scale=pixel_scale)
    logger.debug(f"PSF image bounds {psf_image.bounds}")

    # loop over PSFs
    for i in range(N_PSF):
        PSF = PSF_list[i]

        # draw PSF 
        star_like_PSF = PSF.withFlux(flux)
        logger.debug(f'PSF flux {flux} added.')
        ## draw stamp
        if pixelPSF:
            draw_method = 'no_pixel'
        else:
            draw_method = 'auto'
        stamp_PSF = star_like_PSF.drawImage(scale=pixel_scale, method=draw_method)
        logger.debug(f'PSF drawed (pixel {pixel_scale}).')

        # position the stamp
        stamp_PSF.setCenter(x[i], y[i])

        # stitch the postage stamp to the canvas
        overlap = stamp_PSF.bounds & psf_image.bounds
        psf_image[overlap] += stamp_PSF[overlap]
        logger.debug(f'stamp PSF stitched to canvas.')

    return psf_image

def PSFmap_DiffMag(PSF, pixel_scale, mag_inputs, mag_zero=30., sep_type='random', rng_seed=940120, area=1., pixelPSF=False):
    """ 
    Generate PSF stars from a single PSF model,
        with variable magnitudes.

    Parameters
    ----------
    PSF : galsim object:
        PSF model.
    pixel_scale (arcsec) : float
        Pixel size in unit of arcsec.
    mag_inputs : list
        A list of input magnitudes.
    mag_zero : float, optional (default: 30.)
        The zero point for magnitude.
    sep_type : str, optional (default: 'random')
        Seperation type between PSFs.        
    rng_seed : int, optional (default: 940120) 
        Seed for random number generator.
    area (deg^2) : float, optional (default: 1)
        Total area for simulated image (used when sep_type='random')
    pixelPSF : bool, optional (default: False)
        if the PSF provided already including pixel response

    Returns
    -------
    psf_image : galsim Images
        Image for PSFs.
    """

    mag_inputs = np.array(mag_inputs)
    
    # number of stars
    N_PSF = len(mag_inputs)

    # magnitude to flux
    fluxs = 10**(-0.4*(mag_inputs-mag_zero))

    # position
    if sep_type == 'random':
        # 1degree range
        high_pixel = int(area**0.5*3600./pixel_scale)
        ## for x 
        np.random.seed(rng_seed)
        x = np.random.randint(low=0, high=high_pixel, size=N_PSF)
        ## for y
        np.random.seed(rng_seed+94)
        y = np.random.randint(low=0, high=high_pixel, size=N_PSF)

    # initiate a canvas
    psf_image = galsim.ImageF(int(high_pixel+10), int(high_pixel+10), scale=pixel_scale)

    # place PSF
    for i in range(N_PSF):
        flux = fluxs[i]

        # draw PSF 
        star_like_PSF = PSF.withFlux(flux)
        ## draw stamp
        if pixelPSF:
            draw_method = 'no_pixel'
        else:
            draw_method = 'auto'
        stamp_PSF = star_like_PSF.drawImage(scale=pixel_scale, method=draw_method)

        # position the stamp
        stamp_PSF.setCenter(x[i], y[i])

        # stitch the postage stamp to the canvas
        overlap = stamp_PSF.bounds & psf_image.bounds
        psf_image[overlap] += stamp_PSF[overlap]

    return psf_image

def PSFmap_MultiPSF_DiffMag(PSF_list, pixel_scale, mag_list, mag_zero=30., sep_PSF=120, rng_seed=940120, pixelPSF=False):
    """ 
    Generate PSF stars from a list of PSF models,
        with variable magnitudes.

    Parameters
    ----------
    PSF_list : a list of galsim object
        PSF models.
    pixel_scale (arcsec) : float
        Pixel size in unit of arcsec.
    mag_list : list
        A list of the input magnitude for PSF.
    mag_zero : float, optional (default: 30.)
        The zero point for magnitude.
    sep_PSF (pixel) : int, optional (default: 120)
        Seperation between PSFs.        
    rng_seed : int, optional (default: 940120) 
        Seed for random number generator.
    pixelPSF : bool, optional (default: False)
        if the PSF provided already including pixel response

    Returns
    -------
    psf_image : galsim Images
        Image for PSFs.
    """

    # magnitude to flux
    flux_list = 10**(-0.4*(mag_list-mag_zero))

    # number of psf
    N_PSF = len(PSF_list)

    # position
    # separation in pixels
    separation = int(sep_PSF)
    ## number in each row
    Nrow = int(N_PSF**0.5)
    ## get grid center
    x = np.arange(separation, separation+Nrow*separation, separation, dtype='int')
    y = np.repeat(x, Nrow)
    x = np.tile(x, Nrow)
    ## check outliers
    Nrow = N_PSF - len(x)
    if Nrow > 0:
        x = np.concatenate([x, np.arange(separation, separation+Nrow*separation, separation)])
        y = np.concatenate([y, np.full(Nrow, y[-1]+separation)])
    elif Nrow < 0:
        x = x[:N_PSF]
        y = y[:N_PSF]
    ## make random shift
    np.random.seed(rng_seed)
    shift_lim = int(separation/5.)
    dx_dy = np.random.randint(low=-shift_lim, high=shift_lim, size=N_PSF)
    x += dx_dy
    y += dx_dy    

    # initiate a canvas
    psf_image = galsim.ImageF(int(np.amax(x)+separation), int(np.max(y)+separation), scale=pixel_scale)
    logger.debug(f"PSF image bounds {psf_image.bounds}")

    # loop over PSFs
    for i in range(N_PSF):
        PSF = PSF_list[i]
        flux = flux_list[i]

        # draw PSF 
        star_like_PSF = PSF.withFlux(flux)
        logger.debug(f'PSF flux {flux} added.')
        ## draw stamp
        if pixelPSF:
            draw_method = 'no_pixel'
        else:
            draw_method = 'auto'
        stamp_PSF = star_like_PSF.drawImage(scale=pixel_scale, method=draw_method)
        logger.debug(f'PSF drawed (pixel {pixel_scale}).')

        # position the stamp
        stamp_PSF.setCenter(x[i], y[i])

        # stitch the postage stamp to the canvas
        overlap = stamp_PSF.bounds & psf_image.bounds
        psf_image[overlap] += stamp_PSF[overlap]
        logger.debug(f'stamp PSF stitched to canvas.')

    return psf_image