# -*- coding: utf-8 -*-
# @Author: lshuns
# @Date:   2020-11-30 16:54:44
# @Last modified by:   lshuns
# @Last modified time: 2021-04-30, 19:10:17

### Everything about KiDS observations (instrumental setup & data acquisition procedure)
###### reference: https://www.eso.org/sci/facilities/paranal/instruments/omegacam/doc/VST-MAN-OCM-23100-3110_p107_v2.pdf

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
Dither_x = 25./0.214 # pix # dither step along RA
Dither_y = 85./0.214 # pix # dither step along dec

def getKiDStile(image_ori, id_exposure=0):
    """
    Get OmegaCAM-like tile: cut (dither) + gaps
        The WCS will be handled by GalSim automatically.

    Parameters
    ----------
    image_ori: GalSim image object
        The original image.
    id_exposure: int, optional (default: 0)
        ID for exposure (start with 0)

    Returns
    -------
    image_tile: GalSim image object
        OmegaCAM-like image with gaps.
    weights_tile:
        Associated weight image.
    """

    # original bounds
    x_min0 = image_ori.bounds.xmin
    x_max0 = image_ori.bounds.xmax
    y_min0 = image_ori.bounds.ymin
    y_max0 = image_ori.bounds.ymax

    # bounds for the new image
    center_ori = image_ori.center
    x_shift = center_ori.x - (id_exposure-2)*Dither_x
    y_shift = center_ori.y - (id_exposure-2)*Dither_y
    x_min = int(x_shift - Npix_x/2.)
    x_max = int(x_shift + Npix_x/2.) - 1
    y_min = int(y_shift - Npix_y/2.)
    y_max = int(y_shift + Npix_y/2.) - 1
    bounds = galsim.BoundsI(xmin=max(x_min0, x_min), xmax=min(x_max0, x_max), ymin=max(y_min0, y_min), ymax=min(y_max0, y_max))

    # take account for the out-box effect
    dx_min = x_min - bounds.xmin
    dy_min = y_min - bounds.ymin

    # cut out the new image
    image_tile = image_ori[bounds].copy()
    ## associated weight image
    weights_tile = image_tile.copy()
    weights_tile.fill(1.)

    # gaps
    ## between the long sides of the CCDs
    for i_gap in range(7):
        col_gap_start = (i_gap+1)*Npix_chip_x + i_gap*GAP_x + dx_min
        col_gap_end = min(col_gap_start+GAP_x, image_tile.array.shape[1])
        if col_gap_start < image_tile.array.shape[1]:
            image_tile.array[:, col_gap_start:col_gap_end] = 0.
            weights_tile.array[:, col_gap_start:col_gap_end] = 0.

    ## central gap along the short sides
    row_gap_start = 2*Npix_chip_y + GAP_y_wide + dy_min
    row_gap_end = min(row_gap_start+GAP_y_narrow, image_tile.array.shape[0])
    if row_gap_start < image_tile.array.shape[0]:
        image_tile.array[row_gap_start:row_gap_end, :] = 0.
        weights_tile.array[row_gap_start:row_gap_end, :] = 0.

    ## wide gap along short sides
    #### 1
    row_gap_start = Npix_chip_y + dy_min
    row_gap_end = min(row_gap_start+GAP_y_wide, image_tile.array.shape[0])
    if row_gap_start < image_tile.array.shape[0]:
        image_tile.array[row_gap_start:row_gap_end, :] = 0.
        weights_tile.array[row_gap_start:row_gap_end, :] = 0.
    #### 2
    row_gap_start = 3*Npix_chip_y + GAP_y_wide + GAP_y_narrow + dy_min
    row_gap_end = min(row_gap_start+GAP_y_wide, image_tile.array.shape[0])
    if row_gap_start < image_tile.array.shape[0]:
        image_tile.array[row_gap_start:row_gap_end, :] = 0.
        weights_tile.array[row_gap_start:row_gap_end, :] = 0.

    return image_tile, weights_tile

def getKiDSchips(image_ori, id_exposure=0):
    """
    Get OmegaCAM-like chips: cut image to chips by applying dither and gaps
        can be the original simulated image or tile images
    Parameters
    ----------
    image_ori: GalSim image object
        The original image.
    id_exposure: int, optional (default: 0)
        ID for exposure (start with 0)

    Returns
    -------
    image_chips: a list of GalSim image objects
        OmegaCAM-like chips.
    """

    # original bounds
    x_min0 = image_ori.bounds.xmin
    x_max0 = image_ori.bounds.xmax
    y_min0 = image_ori.bounds.ymin
    y_max0 = image_ori.bounds.ymax

    # bounds due to the dither
    center_ori = image_ori.center
    x_shift = center_ori.x - (id_exposure-2)*Dither_x
    y_shift = center_ori.y - (id_exposure-2)*Dither_y
    x_min = int(x_shift - Npix_x/2.)
    y_min = int(y_shift - Npix_y/2.)

    # y start point for each row of chips
    y_start_list = [y_min,
                        y_min + Npix_chip_y + GAP_y_wide,
                        y_min + 2*Npix_chip_y + GAP_y_wide + GAP_y_narrow,
                        y_min + 3*Npix_chip_y + 2*GAP_y_wide + GAP_y_narrow]

    # chips
    image_chips = []
    for y_start in y_start_list:
        y_min_chip = int(y_start)
        y_max_chip = int(y_start + Npix_chip_y) - 1
        for i_chip in range(Nchips_x):
            x_min_chip = int(x_min + i_chip*(Npix_chip_x+GAP_x))
            x_max_chip = int(x_min_chip + Npix_chip_x) - 1

            bounds_chip = galsim.BoundsI(xmin=max(x_min0, x_min), xmax=min(x_max0, x_max), ymin=max(y_min0, y_min), ymax=min(y_max0, y_max))
            image_chip = image_ori[bounds_chip].copy()
            image_chips.append(image_chip)

    return image_chips
