# -*- coding: utf-8 -*-
# @Author: lshuns
# @Date:   2020-11-26 16:03:10
# @Last Modified by:   lshuns
# @Last Modified time: 2023-10-10 17:08:37

### Everything about celestial objects

import math

import galsim
import logging
import numpy as np

logger = logging.getLogger(__name__)

# some sensible constraints
## size cut
### (back to bulge/disk profile)
RE_CUT = [0, 100.] # arcsec
## sersic cut for sersic profile
### (back to bulge/disk profile)
SERSIC_N_CUT = [0.2, 8.]
## axis raito limits
### (larger or smaller values are set to equal the limiting values)
Q_MIN, Q_MAX = 0.05, 1.0
## sersic index limiting set by GalSim
### (larger or smaller values are set to equal the limiting values)
SERSIC_N_MIN, SERSIC_N_MAX = 0.3, 6.2
## truncate the profile
### (for faster calculation)
TRUNC_FACTOR = 5

def SimpleCanvas(RA_min, RA_max, DEC_min, DEC_max, pixel_scale, edge_sep=18.):
    """
    Build a simple canvas
    """

    ## center used as reference point
    RA0 = (RA_min + RA_max) / 2.
    DEC0 = (DEC_min + DEC_max) / 2.

    # decide bounds
    xmax = (RA_max - RA_min) * 3600. / pixel_scale + 2.*edge_sep/pixel_scale # edge_sep in both sides to avoid edge effects
    ymax = (DEC_max - DEC_min) * 3600. / pixel_scale + 2.*edge_sep/pixel_scale
    bounds = galsim.BoundsI(xmin=0, xmax=math.ceil(xmax), ymin=0, ymax=math.ceil(ymax))

    # build the wcs
    ## Linear projection matrix
    dudx = pixel_scale/3600. # degree
    dudy = 0.
    dvdx = 0.
    dvdy = pixel_scale/3600. # degree
    ## Reference pixel in image
    origin_ima = galsim.PositionI(x=int(bounds.getXMax()/2.), y=int(bounds.getYMax()/2.))
    ## Reference point in wcs
    world_origin = galsim.CelestialCoord(ra=RA0*galsim.degrees, dec=DEC0*galsim.degrees)
    ##
    wcs_affine = galsim.AffineTransform(dudx, dudy, dvdx, dvdy, origin=origin_ima)
    wcs = galsim.TanWCS(wcs_affine, world_origin, units=galsim.degrees)

    canvas = galsim.ImageF(bounds=bounds, wcs=wcs)

    return canvas

def GalaxiesImage(canvas, band, pixel_scale, PSF,
                    gals_info, gal_rotation_angle=0., g_cosmic=[0, 0], gal_position_type=['true', 18.],
                    g_const=True, 
                    pixelPSF=False):
    """
    Generate pure sky image with only galaxies

    Parameters
    ----------
    canvas: galsim Image object
        the canvas to put galaxies
    band: str
        Photometric band simulated.
    pixel_scale (arcsec) : float
        Pixel size in unit of arcsec.
    PSF : galsim object
        PSF model.
    gals_info : DataFrame
        Galaxies information.
    gal_rotation_angle : float, optional (default: 0)
        Rotating galaxies to concel shape noise.
    g_cosmic : float, optional (default: [0,0])
        Cosmic shear.
    gal_position_type : [str, float], optional (default: ['true', 18])
        galaxies position type, mainly useful for grid case.
    g_const : bool, optional (default: True)
        useing constant shear or variable ones
    pixelPSF : bool, optional (default: False)
        if the PSF provided already including pixel response

    Returns
    -------
    full_image: galsim Images
        Image for galaxies.
    """

    # check if the canvas is clean
    if np.any(canvas.array):
        raise Exception('The provided canvas is unclean for GalaxiesImage!')

    # constrain gal stamp for grid mode
    if gal_position_type[0] == 'grid':
        bounds_stamp = galsim.BoundsI(xmin=1, xmax=math.floor(gal_position_type[1]/pixel_scale), ymin=1, ymax=math.floor(gal_position_type[1]/pixel_scale))
        logger.info(f'galaxy stamp is constrained to {bounds_stamp}')

    # copy the canvas
    full_image = canvas.copy()
    wcs = full_image.wcs
    bounds = full_image.bounds

    # galaxy sky positions
    RA_gals = np.array(gals_info['RA']) # degree
    DEC_gals = np.array(gals_info['DEC']) # degree

    # galaxy image positions
    x_gals, y_gals = wcs.toImage(RA_gals, DEC_gals, units=galsim.degrees)
    logger.debug(f'Total number of input galaxies: {len(x_gals)}')
    del RA_gals, DEC_gals, wcs
    ## 0.5 for offset (difference between GalSim and Sextractor)
    x_gals += 0.5
    y_gals += 0.5

    # ignore those out of the boundaries
    mask_tmp = (x_gals>=bounds.xmin) & (x_gals<=bounds.xmax) & (y_gals>=bounds.ymin) & (y_gals<=bounds.ymax)
    x_gals = x_gals[mask_tmp]
    y_gals = y_gals[mask_tmp]
    gals_info_selec = gals_info[mask_tmp].copy()
    gals_info_selec.reset_index(drop=True, inplace=True)
    logger.debug(f'Number of galaxies within the bounds: {len(x_gals)}')
    del mask_tmp, bounds, gals_info

    # get detailed positions for galsim
    ## to int
    ix_gals = np.int32(np.floor(x_gals + 0.5))
    iy_gals = np.int32(np.floor(y_gals + 0.5))
    ## offset
    dx_gals = x_gals - ix_gals
    dy_gals = y_gals - iy_gals
    del x_gals, y_gals

    # loop over galaxies
    for i_gal, gal_info in gals_info_selec.iterrows():

        # galaxy position
        ## int position
        ix_gal = ix_gals[i_gal]
        iy_gal = iy_gals[i_gal]
        ## offset
        dx = dx_gals[i_gal]
        dy = dy_gals[i_gal]
        offset_gal = galsim.PositionD(dx, dy)

        # flux
        flux_gal = gal_info[band]

        # position angle
        PA_gal = gal_info['position_angle']
        ## rotating position angle for shape noise cancellation
        PA_gal += gal_rotation_angle

        # +++ Simulation
        n_gal = gal_info['sersic_n']
        re_gal = gal_info['Re']
        if (n_gal >= SERSIC_N_CUT[0]) and (n_gal <= SERSIC_N_CUT[1]) and (re_gal >= RE_CUT[0]) and (re_gal <= RE_CUT[1]):
            ### sersic profile
            # allowed sersic index range
            if (n_gal < SERSIC_N_MIN) or (n_gal > SERSIC_N_MAX):
                # print(f'Warning...........n_gal {n_gal} outrange of ({SERSIC_N_MIN}, {SERSIC_N_MAX})!')
                n_gal = float(np.where(n_gal<SERSIC_N_MIN, SERSIC_N_MIN, SERSIC_N_MAX))
                # print(f'...........assign {n_gal} for now!')
            q_gal = gal_info['axis_ratio']
            if  (q_gal < Q_MIN) or (q_gal > Q_MAX):
                # print(f"Warning...........q_gal {q_gal} outrange of ({Q_MIN}, {Q_MAX})!")
                q_gal = float(np.where(q_gal<Q_MIN, Q_MIN, Q_MAX))
                # print(f'...........assign {q_gal} for now!')
            re_gal *= (q_gal)**0.5 # account for the ellipticity
            galaxy = galsim.Sersic(n=n_gal, half_light_radius=re_gal, flux=flux_gal, trunc=TRUNC_FACTOR*re_gal, flux_untruncated=True)
            # intrinsic ellipticity
            galaxy = galaxy.shear(q=q_gal, beta=PA_gal*galsim.degrees)
        else:
            ### bulge + disk
            # bulge
            bulge_fraction = gal_info['bulge_fraction']
            bulge_n = gal_info['bulge_n']
            bulge_q = gal_info['bulge_axis_ratio']
            if  (bulge_q < Q_MIN) or (bulge_q > Q_MAX):
                bulge_q = float(np.where(bulge_q<Q_MIN, Q_MIN, Q_MAX))
            bulge_Re = gal_info['bulge_Re'] * (bulge_q)**0.5 # account for the ellipticity
            if (abs(bulge_n-4.)<1e-2):
                bulge_gal = galsim.DeVaucouleurs(half_light_radius=bulge_Re, flux=1.0, trunc=TRUNC_FACTOR*bulge_Re, flux_untruncated=True)
            else:
                bulge_gal = galsim.Sersic(n=bulge_n, half_light_radius=bulge_Re, flux=1.0, trunc=TRUNC_FACTOR*bulge_Re, flux_untruncated=True)
            # intrinsic ellipticity
            bulge_gal = bulge_gal.shear(q=bulge_q, beta=PA_gal*galsim.degrees)

            # disk
            if bulge_fraction < 1:
                disk_q = gal_info['disk_axis_ratio']
                if  (disk_q < Q_MIN) or (disk_q > Q_MAX):
                    disk_q = float(np.where(disk_q<Q_MIN, Q_MIN, Q_MAX))
                disk_Re = gal_info['disk_Re'] * (disk_q)**0.5 # account for the ellipticity
                disk_gal = galsim.Exponential(half_light_radius=disk_Re, flux=1.0)
                # intrinsic ellipticity
                disk_gal = disk_gal.shear(q=disk_q, beta=PA_gal*galsim.degrees)

                galaxy = flux_gal * (bulge_fraction * bulge_gal + (1 - bulge_fraction) * disk_gal)

            else:
                galaxy = flux_gal * bulge_gal

        # cosmic shear
        if g_const:
            galaxy = galaxy.shear(g1=g_cosmic[0], g2=g_cosmic[1])
        else:
            galaxy = galaxy.shear(g1=gal_info['gamma1'], g2=gal_info['gamma2'])

        # convolve with the PSF
        galaxy = galsim.Convolve(galaxy, PSF)

        # draw image
        if pixelPSF:
            draw_method = 'no_pixel'
        else:
            draw_method = 'auto'
        if gal_position_type[0] == 'grid':
            canvas_stamp = galsim.ImageF(bounds=bounds_stamp, scale=pixel_scale)
            stamp_gal = galaxy.drawImage(image=canvas_stamp, offset=offset_gal, method=draw_method)
        else:
            stamp_gal = galaxy.drawImage(scale=pixel_scale, offset=offset_gal, method=draw_method)

        # locate the stamp
        stamp_gal.setCenter(ix_gal, iy_gal)

        # stitch the postage stamp to the canvas
        overlap = stamp_gal.bounds & full_image.bounds
        full_image[overlap] += stamp_gal[overlap]

    return full_image

def GalaxiesImage_casual(canvas, band, pixel_scale, PSF,
                    gals_info_casual, gal_rotation_angle=0., g_cosmic=[0, 0], gal_position_type=['true', 18.],
                    g_const=True, 
                    pixelPSF=False):
    """
    Generate galaxies with casual mode
    require gals_info_casual produced by LoadCata.GalInfo_adjust4casual

    Parameters
    ----------
    canvas: galsim Image object
        the canvas to put galaxies
    band: str
        Photometric band simulated.
    pixel_scale (arcsec) : float
        Pixel size in unit of arcsec.
    PSF : galsim object
        PSF model.
    gals_info_casual : DataFrame
        Galaxies information with casual mode support.
    gal_rotation_angle : float, optional (default: 0)
        Rotating galaxies to concel shape noise.
    g_cosmic : float, optional (default: [0,0])
        Cosmic shear.
    gal_position_type : [str, float], optional (default: ['true', 18])
        galaxies position type, mainly useful for grid case.
    g_const : bool, optional (default: True)
        useing constant shear or variable ones
    pixelPSF : bool, optional (default: False)
        if the PSF provided already including pixel response

    Returns
    -------
    full_image: galsim Images
        Image for galaxies.
    """

    # we will add something, so do a copy
    gals_info_casual = gals_info_casual.copy()

    # check if the canvas is clean
    if np.any(canvas.array):
        raise Exception('The provided canvas is unclean for GalaxiesImage_casual!')

    # constrain gal stamp for grid mode
    if gal_position_type[0] == 'grid':
        bounds_stamp = galsim.BoundsI(xmin=1, xmax=math.floor(gal_position_type[1]/pixel_scale), ymin=1, ymax=math.floor(gal_position_type[1]/pixel_scale))
        logger.info(f'galaxy stamp is constrained to {bounds_stamp}')

    # copy the canvas
    full_image = canvas.copy()
    wcs = full_image.wcs
    bounds = full_image.bounds

    # galaxy sky positions
    RA_gals = np.array(gals_info_casual['RA']) # degree
    DEC_gals = np.array(gals_info_casual['DEC']) # degree

    # galaxy image positions
    x_gals, y_gals = wcs.toImage(RA_gals, DEC_gals, units=galsim.degrees)
    logger.debug(f'Total number of input galaxies: {len(x_gals)}')
    del RA_gals, DEC_gals, wcs
    ## 0.5 for offset (difference between GalSim and Sextractor)
    x_gals += 0.5
    y_gals += 0.5
    ## save to dataframe
    gals_info_casual.loc[:, 'x'] = x_gals
    gals_info_casual.loc[:, 'y'] = y_gals

    # ignore those out of the boundaries
    mask_tmp = (x_gals>=bounds.xmin) & (x_gals<=bounds.xmax) & (y_gals>=bounds.ymin) & (y_gals<=bounds.ymax)
    gals_info_selec = gals_info_casual[mask_tmp].copy()
    gals_info_selec.reset_index(drop=True, inplace=True)
    logger.debug(f'Number of galaxies within the bounds: {len(x_gals)}')
    del x_gals, y_gals, mask_tmp, bounds, gals_info_casual

    ## positions to int
    #### NOTE: offset is ignored in casual mode
    gals_info_selec.loc[:, 'x'] = np.int32(np.floor(gals_info_selec['x'] + 0.5))
    gals_info_selec.loc[:, 'y'] = np.int32(np.floor(gals_info_selec['y'] + 0.5))

    # bin the cata based on the qbin id
    info_grouped = gals_info_selec.groupby(by='i_qbin')
    del gals_info_selec
    
    # loop over the bins and draw galaxies
    for name, group in info_grouped:

        # get the seed galaxies
        gals0 = group.drop_duplicates(subset='index_seedGal', keep='first', ignore_index=True)

        # get stamps from seed galaxies
        stamps0 = {}
        for i_gal, gal_info in gals0.iterrows():

            # id 
            index_seedGal = gal_info['index_seedGal']

            # flux
            flux_gal = gal_info[band]

            # position angle
            PA_gal = gal_info['position_angle']
            ## rotating position angle for shape noise cancellation
            PA_gal += gal_rotation_angle

            # +++ Simulation
            n_gal = gal_info['sersic_n']
            re_gal = gal_info['Re']
            if (n_gal >= SERSIC_N_CUT[0]) and (n_gal <= SERSIC_N_CUT[1]) and (re_gal >= RE_CUT[0]) and (re_gal <= RE_CUT[1]):
                ### sersic profile
                # allowed sersic index range
                if (n_gal < SERSIC_N_MIN) or (n_gal > SERSIC_N_MAX):
                    # print(f'Warning...........n_gal {n_gal} outrange of ({SERSIC_N_MIN}, {SERSIC_N_MAX})!')
                    n_gal = float(np.where(n_gal<SERSIC_N_MIN, SERSIC_N_MIN, SERSIC_N_MAX))
                    # print(f'...........assign {n_gal} for now!')
                q_gal = gal_info['axis_ratio']
                if  (q_gal < Q_MIN) or (q_gal > Q_MAX):
                    # print(f"Warning...........q_gal {q_gal} outrange of ({Q_MIN}, {Q_MAX})!")
                    q_gal = float(np.where(q_gal<Q_MIN, Q_MIN, Q_MAX))
                    # print(f'...........assign {q_gal} for now!')
                re_gal *= (q_gal)**0.5 # account for the ellipticity
                galaxy = galsim.Sersic(n=n_gal, half_light_radius=re_gal, flux=flux_gal, trunc=TRUNC_FACTOR*re_gal, flux_untruncated=True)
                # intrinsic ellipticity
                galaxy = galaxy.shear(q=q_gal, beta=PA_gal*galsim.degrees)
            else:
                ### bulge + disk
                # bulge
                bulge_fraction = gal_info['bulge_fraction']
                bulge_n = gal_info['bulge_n']
                bulge_q = gal_info['bulge_axis_ratio']
                if  (bulge_q < Q_MIN) or (bulge_q > Q_MAX):
                    bulge_q = float(np.where(bulge_q<Q_MIN, Q_MIN, Q_MAX))
                bulge_Re = gal_info['bulge_Re'] * (bulge_q)**0.5 # account for the ellipticity
                if (abs(bulge_n-4.)<1e-2):
                    bulge_gal = galsim.DeVaucouleurs(half_light_radius=bulge_Re, flux=1.0, trunc=TRUNC_FACTOR*bulge_Re, flux_untruncated=True)
                else:
                    bulge_gal = galsim.Sersic(n=bulge_n, half_light_radius=bulge_Re, flux=1.0, trunc=TRUNC_FACTOR*bulge_Re, flux_untruncated=True)
                # intrinsic ellipticity
                bulge_gal = bulge_gal.shear(q=bulge_q, beta=PA_gal*galsim.degrees)

                # disk
                if bulge_fraction < 1:
                    disk_q = gal_info['disk_axis_ratio']
                    if  (disk_q < Q_MIN) or (disk_q > Q_MAX):
                        disk_q = float(np.where(disk_q<Q_MIN, Q_MIN, Q_MAX))
                    disk_Re = gal_info['disk_Re'] * (disk_q)**0.5 # account for the ellipticity
                    disk_gal = galsim.Exponential(half_light_radius=disk_Re, flux=1.0)
                    # intrinsic ellipticity
                    disk_gal = disk_gal.shear(q=disk_q, beta=PA_gal*galsim.degrees)

                    galaxy = flux_gal * (bulge_fraction * bulge_gal + (1 - bulge_fraction) * disk_gal)

                else:
                    galaxy = flux_gal * bulge_gal

            # cosmic shear
            if g_const:
                galaxy = galaxy.shear(g1=g_cosmic[0], g2=g_cosmic[1])
            else:
                galaxy = galaxy.shear(g1=gal_info['gamma1'], g2=gal_info['gamma2'])

            # convolve with the PSF
            galaxy = galsim.Convolve(galaxy, PSF)

            # draw image
            if pixelPSF:
                draw_method = 'no_pixel'
            else:
                draw_method = 'auto'
            if gal_position_type[0] == 'grid':
                canvas_stamp = galsim.ImageF(bounds=bounds_stamp, scale=pixel_scale)
                stamp_gal = galaxy.drawImage(image=canvas_stamp, method=draw_method)
            else:
                stamp_gal = galaxy.drawImage(scale=pixel_scale, method=draw_method)

            # collect stamps
            stamps0[int(index_seedGal)] = stamp_gal

        del gals0

        # put galaxies
        for i_gal, gal_info in group.iterrows():

            # get position
            ix_gal = gal_info['x']
            iy_gal = gal_info['y']

            # get stamp
            stamp_gal = stamps0[int(gal_info['index_seedGal'])].copy()

            # locate the stamp
            stamp_gal.setCenter(ix_gal, iy_gal)

            # stitch the postage stamp to the canvas
            overlap = stamp_gal.bounds & full_image.bounds
            full_image[overlap] += stamp_gal[overlap]

        del stamps0, group

    return full_image

def StarsImage(canvas, band, pixel_scale, PSF,
                    stars_info, 
                    pixelPSF=False):
    """
    Generate pure sky image with only stars

    Parameters
    ----------
    canvas: galsim Image object
        the canvas to put galaxies
    band: str
        Photometric band simulated.
    pixel_scale (arcsec) : float
        Pixel size in unit of arcsec.
    PSF : galsim object
        PSF model.
    stars_info : DataFrame
        Stars information.
    pixelPSF : bool, optional (default: False)
        if the PSF provided already including pixel response

    Returns
    -------
    full_image: galsim Images
        Image for stars.
    """

    # check if the canvas is clean
    if np.any(canvas.array):
        raise Exception('The provided canvas is unclean for StarsImage!')

    # copy the canvas
    full_image = canvas.copy()
    wcs = full_image.wcs
    bounds = full_image.bounds

    # sky positions
    RA_stars = np.array(stars_info['RA']) # degree
    DEC_stars = np.array(stars_info['DEC']) # degree

    # image positions
    x_stars, y_stars = wcs.toImage(RA_stars, DEC_stars, units=galsim.degrees)
    logger.debug(f'Total number of input stars: {len(x_stars)}')
    del RA_stars, DEC_stars, wcs

    # ignore those out of the boundaries
    mask_tmp = (x_stars>=bounds.xmin) & (x_stars<=bounds.xmax) & (y_stars>=bounds.ymin) & (y_stars<=bounds.ymax)
    x_stars = x_stars[mask_tmp]
    y_stars = y_stars[mask_tmp]
    stars_info_selec = stars_info[mask_tmp].copy()
    stars_info_selec.reset_index(drop=True, inplace=True)
    logger.debug(f'Number of stars within the bounds: {len(x_stars)}')
    del mask_tmp, bounds, stars_info

    # loop over stars
    for i_star, star_info in stars_info_selec.iterrows():

        # position
        ## 0.5 for offset (difference between GalSim and Sextractor)
        x_star = x_stars[i_star] + 0.5
        y_star = y_stars[i_star] + 0.5
        ## to int
        ix_star = int(math.floor(x_star + 0.5))
        iy_star = int(math.floor(y_star + 0.5))
        ## offset
        dx = x_star - ix_star
        dy = y_star - iy_star
        offset_star = galsim.PositionD(dx, dy)

        # flux
        flux_star = star_info[band]

        # generate star
        star = PSF.withFlux(flux_star)
        # star = galsim.Convolve(galsim.DeltaFunction(flux=flux_star), PSF)

        ## draw stamp
        if pixelPSF:
            draw_method = 'no_pixel'
        else:
            draw_method = 'auto'
        stamp_star = star.drawImage(scale=pixel_scale, offset=offset_star, method=draw_method)

        # position the stamp
        stamp_star.setCenter(ix_star, iy_star)

        # stitch the postage stamp to the canvas
        overlap = stamp_star.bounds & full_image.bounds
        full_image[overlap] += stamp_star[overlap]

    return full_image
