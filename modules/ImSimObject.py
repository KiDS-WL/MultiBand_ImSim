# -*- coding: utf-8 -*-
# @Author: lshuns
# @Date:   2020-11-26 16:03:10
# @Last modified by:   lshuns
# @Last modified time: 2021-04-29, 14:06:10

### Everything about celestial objects

import galsim
import numpy as np

# some sensible constraints
## size cut
### (back to bulge/disk profile)
RE_CUT = [1e-2, 10.] # arcsec
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

def WCS(RA_min, RA_max, DEC_min, DEC_max, pixel_scale):
    """
    Build the WCS system for the coordinate transformation
    """

    ## center used as reference point
    RA0 = (RA_min + RA_max) / 2.
    DEC0 = (DEC_min + DEC_max) / 2.

    # decide bounds
    xmax = (RA_max - RA_min) * 3600. / pixel_scale + 36./pixel_scale # 18 arcsec in both sides to avoid edge effects
    ymax = (DEC_max - DEC_min) * 3600. / pixel_scale + 36./pixel_scale
    bounds = galsim.BoundsI(xmin=0, xmax=int(xmax), ymin=0, ymax=int(ymax))

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

    return bounds, wcs

def GalaxiesImage(bounds, wcs, band, pixel_scale, PSF,
                    gals_info, gal_rotation_angle=0., g_cosmic=[0, 0]):
    """
    Generate pure sky image with only galaxies

    Parameters
    ----------
    bounds: galsim object
        Bounds for the whole image.
    wcs: galsim object
        WCS for coordinate transformation.
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

    Returns
    -------
    full_image: galsim Images
        Image for galaxies.
    """

    # galaxy sky positions
    RA_gals = np.array(gals_info['RA']) # degree
    DEC_gals = np.array(gals_info['DEC']) # degree

    # galaxy image positions
    x_gals, y_gals = wcs.toImage(RA_gals, DEC_gals, units=galsim.degrees)

    # initialize the canvas
    full_image = galsim.ImageF(bounds=bounds, wcs=wcs)

    # loop over galaxies
    for i_gal, gal_info in gals_info.iterrows():

        # galaxy position
        ## 0.5 for offset (difference between GalSim and Sextractor)
        x_gal = x_gals[i_gal] + 0.5
        y_gal = y_gals[i_gal] + 0.5
        ## to int
        ix_gal = int(x_gal)
        iy_gal = int(y_gal)
        ## offset
        dx = x_gal - ix_gal
        dy = y_gal - iy_gal
        offset_gal = galsim.PositionD(dx, dy)

        # flux
        flux_gal = gal_info[band]

        # position angle
        PA_gal = gal_info['position_angles']
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
            q_gal = gal_info['axis_ratios']
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
            bulge_fraction = gal_info['bulge_fractions']
            bulge_n = gal_info['bulge_n']
            bulge_q = gal_info['bulge_axis_ratios']
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
                disk_q = gal_info['disk_axis_ratios']
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
        galaxy = galaxy.shear(g1=g_cosmic[0], g2=g_cosmic[1])

        # convolve with the PSF
        galaxy = galsim.Convolve(galaxy, PSF)

        # draw image
        stamp_gal = galaxy.drawImage(scale=pixel_scale, offset=offset_gal)

        # locate the stamp
        stamp_gal.setCenter(ix_gal, iy_gal)

        # stitch the postage stamp to the canvas
        overlap = stamp_gal.bounds & full_image.bounds
        full_image[overlap] += stamp_gal[overlap]

    return full_image

def StarsImage(bounds, wcs, band, pixel_scale, PSF,
                    stars_info):
    """
    Generate pure sky image with only stars

    Parameters
    ----------
    bounds: galsim object
        Bounds for the whole image.
    wcs: galsim object
        WCS for coordinate transformation.
    band: str
        Photometric band simulated.
    pixel_scale (arcsec) : float
        Pixel size in unit of arcsec.
    PSF : galsim object
        PSF model.
    stars_info : DataFrame
        Stars information.

    Returns
    -------
    full_image: galsim Images
        Image for stars.
    """

    # sky positions
    RA_stars = stars_info['RA'] # degree
    DEC_stars = stars_info['DEC'] # degree

    # image positions
    x_stars, y_stars = wcs.toImage(RA_stars, DEC_stars, units=galsim.degrees)

    # initialize the canvas
    full_image = galsim.ImageF(bounds=bounds, wcs=wcs)

    # loop over stars
    for i_star, star_info in stars_info.iterrows():

        # position
        ## 0.5 for offset (difference between GalSim and Sextractor)
        x_star = x_stars[i_star] + 0.5
        y_star = y_stars[i_star] + 0.5
        ## to int
        ix_star = int(x_star)
        iy_star = int(y_star)
        ## offset
        dx = x_star - ix_star
        dy = y_star - iy_star
        offset_star = galsim.PositionD(dx, dy)

        # flux
        flux_star = star_info[band]

        # draw star
        star = flux_star * PSF

        ## draw stamp
        stamp_star = star.drawImage(scale=pixel_scale, offset=offset_star)

        # position the stamp
        stamp_star.setCenter(ix_star, iy_star)

        # stitch the postage stamp to the canvas
        overlap = stamp_star.bounds & full_image.bounds
        full_image[overlap] += stamp_star[overlap]

    return full_image
