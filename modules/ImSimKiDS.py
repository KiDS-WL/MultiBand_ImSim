# -*- coding: utf-8 -*-
# @Author: lshuns
# @Date:   2020-11-30 16:54:44
# @Last modified by:   lshuns
# @Last modified time: 2021-06-09, 17:37:40

### Everything about KiDS observations (instrumental setup & data acquisition procedure)
###### reference: https://www.eso.org/sci/facilities/paranal/instruments/omegacam/doc/VST-MAN-OCM-23100-3110_p107_v2.pdf

import math
import galsim
import numpy as np

# some OmegaCAM values
Pixel_scale = 0.214/3600. # degree
Nchips_x = 8 # number of CCDs (x-axis)
Nchips_y = 4 # number of CCDs (y-axis)
Npix_chip_x = 2048 # pixels each CCD (x-axis)
Npix_chip_y = 4100 # pixels each CCD (y-axis)
GAP_x = 100 # pix # gap between the long sides of the CCDs
GAP_y_narrow = 55 # pix # central gap along the short sides
GAP_y_wide = 376 # pix # wide gap along short sides
Npix_x = Nchips_x*Npix_chip_x + 7*GAP_x  # total number of pixels each tile (x-axis)
Npix_y = Nchips_y*Npix_chip_y + GAP_y_narrow + 2*GAP_y_wide # total number of pixels each tile (y-axis)

# some observation values
# N_exposures_gri = 5 # number of exposures in gri-bands
# N_exposures_u = 4 # number of exposures in u-band
Dither_x_pix = 25./0.214 # pix # dither step along RA
Dither_y_pix = 85./0.214 # pix # dither step along dec
Dither_x_arcsec = 25. # arcsec # dither step along RA
Dither_y_arcsec = 85. # arcsec # dither step along dec

def getKiDScanvases(RA0, DEC0, id_exposure=0):
    """
    Build a bunch of canvases mimicking OmegaCAM chips

    Parameters
    ----------
    RA0: float
        The original RA for the reference pixel
    DEC0: float
        The original DEC for the reference pixel
    id_exposure: int, optional (default: 0)
        ID for exposure (start with 0)

    Returns
    -------
    canvases_list: a list of GalSim Image object
        a list of canvases mimicking OmegaCAM chips
    """

    # center value shifted due to dither step
    RA0 -= (id_exposure-2)*(Dither_x_arcsec/3600.)
    DEC0 -= (id_exposure-2)*(Dither_y_arcsec/3600.)

    # transfer to the center value for the first chip
    RA0 -= (Npix_x/2.-Npix_chip_x/2.)*Pixel_scale
    DEC0 -= (Npix_y/2.-Npix_chip_y/2.)*Pixel_scale

    # bounds for each CCD
    bounds = galsim.BoundsI(xmin=0, xmax=Npix_chip_x-1, ymin=0, ymax=Npix_chip_y-1)
    ## use the center pixel as the reference pixel in image
    origin_ima = galsim.PositionI(x=math.floor(Npix_chip_x/2.), y=math.floor(Npix_chip_y/2.))

    # Linear projection matrix for wcs
    dudx = Pixel_scale # degree
    dudy = 0.
    dvdx = 0.
    dvdy = Pixel_scale # degree
    wcs_affine = galsim.AffineTransform(dudx, dudy, dvdx, dvdy, origin=origin_ima)

    # loop over all chips
    canvases_list = []
    for i_chip_x in range(Nchips_x):
        for i_chip_y in range(Nchips_y):

            ## Reference point in wcs
            RA0_tmp = RA0 + i_chip_x*(Npix_chip_x+GAP_x)*Pixel_scale
            DEC0_tmp = DEC0 + i_chip_y*Npix_chip_y*Pixel_scale
            ### gaps in y axis
            if i_chip_y >= 1:
                DEC0_tmp += GAP_y_wide*Pixel_scale
            if i_chip_y >= 2:
                DEC0_tmp += GAP_y_narrow*Pixel_scale
            if i_chip_y >= 3:
                DEC0_tmp += GAP_y_wide*Pixel_scale
            world_origin = galsim.CelestialCoord(ra=RA0_tmp*galsim.degrees, dec=DEC0_tmp*galsim.degrees)

            ## build the wcs
            wcs = galsim.TanWCS(wcs_affine, world_origin, units=galsim.degrees)

            ## build the canvas
            canvases_list.append(galsim.ImageF(bounds=bounds, wcs=wcs))

    return canvases_list

def cutKiDStile(image_ori_noiseless, noise, id_exposure=0):
    """
    Cut OmegaCAM-like tile from simulated image with dither and gaps.

    Parameters
    ----------
    image_ori_noiseless: GalSim image object
        The original image without noise.
    noise: GalSim noise object
        Noise added to the tile image.
    id_exposure: int, optional (default: 0)
        ID for exposure (start with 0)

    Returns
    -------
    image_tile: GalSim image object
        OmegaCAM-like image with gaps.
    weights_tile:
        Associated weight image.
    """

    # bounds for the tile image
    center_ori = image_ori_noiseless.center
    x_shift = center_ori.x - (id_exposure-2)*Dither_x_pix
    y_shift = center_ori.y - (id_exposure-2)*Dither_y_pix
    x_min = x_shift - Npix_x/2.
    x_max = x_shift + Npix_x/2.
    y_min = y_shift - Npix_y/2.
    y_max = y_shift + Npix_y/2.
    bounds_tile = galsim.BoundsI(xmin=math.floor(x_min), xmax=math.floor(x_max)-1, ymin=math.floor(y_min), ymax=math.floor(y_max)-1)

    # initialize the canvas
    image_tile = galsim.ImageF(bounds=bounds_tile, wcs=image_ori_noiseless.wcs)
    ## associated weight image
    weights_tile = image_tile.copy()
    weights_tile.fill(1.)

    # stitch the original image to the canvas
    overlap = image_tile.bounds & image_ori_noiseless.bounds
    image_tile[overlap] += image_ori_noiseless[overlap]

    # add noise
    image_tile.addNoise(noise)

    # gaps
    ## between the long sides of the CCDs
    for i_gap in range(7):
        col_gap_start = (i_gap+1)*Npix_chip_x + i_gap*GAP_x
        image_tile.array[:, col_gap_start:(col_gap_start+GAP_x)] = 0.
        weights_tile.array[:, col_gap_start:(col_gap_start+GAP_x)] = 0.

    ## central gap along the short sides
    row_gap_start = 2*Npix_chip_y + GAP_y_wide
    image_tile.array[row_gap_start:(row_gap_start+GAP_y_narrow), :] = 0.
    weights_tile.array[row_gap_start:(row_gap_start+GAP_y_narrow), :] = 0.

    ## wide gap along short sides
    #### 1
    row_gap_start = Npix_chip_y
    image_tile.array[row_gap_start:(row_gap_start+GAP_y_wide), :] = 0.
    weights_tile.array[row_gap_start:(row_gap_start+GAP_y_wide), :] = 0.
    #### 2
    row_gap_start = 3*Npix_chip_y + GAP_y_wide + GAP_y_narrow
    image_tile.array[row_gap_start:(row_gap_start+GAP_y_wide), :] = 0.
    weights_tile.array[row_gap_start:(row_gap_start+GAP_y_wide), :] = 0.

    return image_tile, weights_tile

def cutKiDSchips(image_tile):
    """
    Cut OmegaCAM-like chips from OmegaCAM-like tile

    Parameters
    ----------
    image_tile: GalSim image object
        The tile image produced by cutKiDStile.

    Returns
    -------
    image_chips: a list of GalSim image objects
        OmegaCAM-like chips.
    """

    # start pixel
    bounds = image_tile.bounds
    x_min = bounds.getXMin()
    y_min = bounds.getYMin()

    # y start point for each row of chips
    y_start_list = [y_min,
                        y_min + Npix_chip_y + GAP_y_wide,
                        y_min + 2*Npix_chip_y + GAP_y_wide + GAP_y_narrow,
                        y_min + 3*Npix_chip_y + 2*GAP_y_wide + GAP_y_narrow]

    # chips
    image_chips = []
    for y_start in y_start_list:
        y_min_chip = y_start
        y_max_chip = y_start + Npix_chip_y
        for i_chip in range(Nchips_x):
            x_min_chip = x_min + i_chip*(Npix_chip_x+GAP_x)
            x_max_chip = x_min_chip + Npix_chip_x
            bounds_chip = galsim.BoundsI(xmin=math.floor(x_min_chip), xmax=math.floor(x_max_chip)-1, ymin=math.floor(y_min_chip), ymax=math.floor(y_max_chip)-1)
            image_chip = image_tile[bounds_chip].copy()
            image_chips.append(image_chip)

    return image_chips
