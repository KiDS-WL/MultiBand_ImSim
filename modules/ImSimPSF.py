# -*- coding: utf-8 -*-
# @Author: lshuns
# @Date:   2020-11-26 15:00:22
# @Last Modified by:   lshuns
# @Last Modified time: 2021-01-20 15:31:15

### Everything about PSF

import galsim
import logging
import numpy as np

logger = logging.getLogger(__name__)

def MoffatPSF(seeing, moffat_beta=3.5, psf_e=[]):
    """
    Generate PSF model.

    Parameters
    ----------
    seeing : float
        Full-width-half-max of the PSF.
    moffat_beta : float
        Beta parameter for Moffat profile.
    psf_e : list, optional (default: [0, 0])
        List of psf ellipticity [e1, e2] (e = 1 - q).

    Returns
    -------
    psf: galsim Moffat
        PSF model.
    """

    psf = galsim.Moffat(beta=moffat_beta, fwhm=seeing, trunc=4.5*seeing)
    if psf_e != []:
        psf_e1 = psf_e[0]
        psf_e2 = psf_e[1]

        psf_e = np.sqrt(psf_e1**2+psf_e2**2)
        # g_i = e_i/(2-e)
        psf_g1, psf_g2 = psf_e1/(2-psf_e), psf_e2/(2-psf_e)

        psf = psf.shear(g1=psf_g1, g2=psf_g2)

    return psf


def MoffatPSF_list(seeing_list, moffat_beta_list, psf_e_list=[]):
    """
    Generate PSF model with variable inputs.

    Parameters
    ----------
    seeing : float
        Full-width-half-max of the PSF.
    moffat_beta : float
        Beta parameter for Moffat profile.
    psf_e : list, optional (default: [0, 0])
        List of psf ellipticity [e1, e2] (e = 1 - q).

    Returns
    -------
    psf: galsim Moffat
        PSF model.
    """

    N_psf = len(seeing_list)

    psf_list = []
    for i_psf in range(N_psf):
        seeing = seeing_list[i_psf]
        moffat_beta = moffat_beta_list[i_psf]

        psf = galsim.Moffat(beta=moffat_beta, fwhm=seeing, trunc=4.5*seeing)

        if psf_e_list != []:
            psf_e = psf_e_list[i_psf]
            psf_e1 = psf_e[0]
            psf_e2 = psf_e[1]

            psf_e = np.sqrt(psf_e1**2+psf_e2**2)
            # g_i = e_i/(2-e)
            psf_g1, psf_g2 = psf_e1/(2-psf_e), psf_e2/(2-psf_e)

            psf=psf.shear(g1=psf_g1, g2=psf_g2)

        psf_list.append(psf)

    return psf_list

def PSFima(PSF, pixel_scale, size=32):
    """
    draw a single PSF image
    """

    psf_image = galsim.Image(size, size)
    PSF_lf = PSF.shift(0.5*pixel_scale, 0.5*pixel_scale)
    psf_image = PSF_lf.drawImage(image=psf_image, scale=pixel_scale)

    return psf_image


def PSFmap(PSF, pixel_scale, mag_input, mag_zero=30., N_PSF=100, sep_PSF=120, rng_seed=940120):
    """ 
    Generate PSF stars for the image
        with uniform magnitude

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

    Returns
    -------
    psf_image: galsim Images
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
    stamp_PSF = star_like_PSF.drawImage(scale=pixel_scale)

    # place PSF
    for i in range(N_PSF):

        # position the stamp
        stamp_PSF.setCenter(x[i], y[i])

        # stitch the postage stamp to the canvas
        overlap = stamp_PSF.bounds & psf_image.bounds
        psf_image[overlap] += stamp_PSF[overlap]

    return psf_image

def PSFmap_MultiPSF(PSF_list, pixel_scale, mag_input, mag_zero=30., sep_PSF=120, rng_seed=940120):
    """ 
    Generate PSF stars for the image
        with uniform magnitude

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

    Returns
    -------
    psf_image: galsim Images
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
        stamp_PSF = star_like_PSF.drawImage(scale=pixel_scale)
        logger.debug(f'PSF drawed (pixel {pixel_scale}).')

        # position the stamp
        stamp_PSF.setCenter(x[i], y[i])

        # stitch the postage stamp to the canvas
        overlap = stamp_PSF.bounds & psf_image.bounds
        psf_image[overlap] += stamp_PSF[overlap]
        logger.debug(f'stamp PSF stitched to canvas.')

    return psf_image


def PSFmap_DiffMag(PSF, pixel_scale, mag_inputs, mag_zero=30., sep_type='random', rng_seed=940120, area=1.):
    """ 
    Generate PSF stars for the image
        with variable magnitudes

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

    Returns
    -------
    psf_image: galsim Images
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
        stamp_PSF = star_like_PSF.drawImage(scale=pixel_scale)

        # position the stamp
        stamp_PSF.setCenter(x[i], y[i])

        # stitch the postage stamp to the canvas
        overlap = stamp_PSF.bounds & psf_image.bounds
        psf_image[overlap] += stamp_PSF[overlap]

    return psf_image


def PSFmap_MultiPSF_DiffMag(PSF_list, pixel_scale, mag_list, mag_zero=30., sep_PSF=120, rng_seed=940120):
    """ 
    Generate PSF stars for the image
        with variable PSF as well as variable magnitude

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

    Returns
    -------
    psf_image: galsim Images
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
        stamp_PSF = star_like_PSF.drawImage(scale=pixel_scale)
        logger.debug(f'PSF drawed (pixel {pixel_scale}).')

        # position the stamp
        stamp_PSF.setCenter(x[i], y[i])

        # stitch the postage stamp to the canvas
        overlap = stamp_PSF.bounds & psf_image.bounds
        psf_image[overlap] += stamp_PSF[overlap]
        logger.debug(f'stamp PSF stitched to canvas.')

    return psf_image