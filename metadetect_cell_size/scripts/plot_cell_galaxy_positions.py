# -*- coding: utf-8 -*-

### Where do the grid galaxies sit within the metadetect cells?
###
###    Rebuilds the galaxy pixel positions of every tile from the ImSim layout
###    rules (the grid placement in ImSim.py, ImSimObject.SimpleCanvas, and for
###    the sheared tags ImSim._canvas_bounds_unsheared and ImSim._shear_positions),
###    cuts every image into cells the way MetaDetect.MetaDetectShear does, and
###    stacks the positions of the galaxies in all cells of all tiles into a 2D
###    histogram over the cell. The dashed square is the central region that is
###    kept, the cross the point metacal shears the cell about.
###
###    Nothing is read from the shape catalogues. The number of galaxies of each
###    tile, which sets the size of its grid, is taken from the ImSim input
###    catalogues, and the rebuilt layout is checked against those catalogues and
###    against the image headers.
###
###    Usage: python plot_cell_galaxy_positions.py

import os
import re
import sys
import glob
import math
import logging

import numpy as np
import pandas as pd
import galsim
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from scipy.spatial import cKDTree

import analyse_cell_size as A
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', 'modules'))
import ImSim
import ImSimObject

logging.getLogger('ImSim').setLevel(logging.WARNING)

## ++++++++++++++ I/O and general setups

## The images, ImSim layout of test_scripts/grid/config_meta_073.ini (collaborator's)
imsim_dir = '/sdf/data/kipac/u/jlkitt/imsims/grid/seeing_073'
GRID_SIZE = 9.        # arcsec between neighbouring galaxies
PIXEL_SCALE = 0.2     # arcsec per pixel
BAND = 'LSST_r'
## the galaxies are drawn at sheared positions (shear_positions = True)

script_dir = os.path.dirname(os.path.abspath(__file__))
outpath_png = os.path.join(script_dir, 'cell_galaxy_positions.png')
outpath_zoom_png = os.path.join(script_dir, 'cell_galaxy_positions_zoom.png')

## the whole cell is shown in bins of REBIN pixels, to see where galaxies
##    concentrate, and the centre in single pixels, to see the positions sampled
REBIN = 5
ZOOM_HALF = 30

## The metadetect set-ups (cell_size, central_size), as in analyse_cell_size.py
SETUPS = sorted({tuple(map(int, re.search(r'cell(\d+)_central(\d+)', label).groups()))
                 for label in A.setups})

## the columns of the figure: one per shear tag, the zero-shear one first
TAG_ZERO = [t for t, g in A.shear_inputs.items() if g == (0., 0.)]
COLUMNS = [(f'{tag}\ng = ({g[0]:+.2f}, {g[1]:+.2f})', [tag])
           for tag, g in sorted(A.shear_inputs.items(), key=lambda kv: kv[1] != (0., 0.))]

## galaxies per pixel relative to the mean, blue below and red above
CMAP_DIV = LinearSegmentedColormap.from_list('div', ['#2a78d6', '#f0efec', '#e34948'])

plt.rc('font', size=11)

## ++++++++++++++ Workhorse

def grid_layout(N_gal, g):
    """
    Galaxy positions (numpy row, col of the image array) and the image shape,
    as ImSim lays out a grid tile with N_gal galaxies under the shear g.
    """
    ## the grid placement of ImSim.py
    apart = GRID_SIZE / 3600.
    N_rows = math.ceil(N_gal**0.5)
    N_cols = math.ceil(N_gal / N_rows)
    X_gals = np.arange(0, apart + (N_cols - 1) * apart, apart)
    Y_gals = np.arange(0, apart + (N_rows - 1) * apart, apart)
    X_gals, Y_gals = np.meshgrid(X_gals, Y_gals)
    ## every grid point is kept here: ImSim leaves N_rows*N_cols - N_gal of them
    ##    (less than one row) empty at random, which does not change the pattern
    gals = pd.DataFrame({'RA': X_gals.ravel(), 'DEC': Y_gals.ravel()})
    gals_list = [gals, None]

    ## the canvas, sized from the unsheared extent, and the position shear
    canvas_bounds = ImSim._canvas_bounds_unsheared(gals_list, g)
    if canvas_bounds is not None:
        canvas = ImSimObject.SimpleCanvas(canvas_bounds['RA_min'], canvas_bounds['RA_max'],
                                          canvas_bounds['DEC_min'], canvas_bounds['DEC_max'],
                                          PIXEL_SCALE, edge_sep=canvas_bounds['edge_sep'])
    else:
        canvas = ImSimObject.SimpleCanvas(gals['RA'].min(), gals['RA'].max(),
                                          gals['DEC'].min(), gals['DEC'].max(), PIXEL_SCALE)
    ImSim._shear_positions(gals_list, g)

    ## ImSimObject.GalaxiesImage draws each galaxy at wcs.toImage(RA, DEC)
    x, y = canvas.wcs.toImage(gals['RA'].values, gals['DEC'].values, units=galsim.degrees)
    bounds = canvas.bounds
    shape = (bounds.ymax - bounds.ymin + 1, bounds.xmax - bounds.xmin + 1)
    return y - bounds.ymin, x - bounds.xmin, shape

def cell_positions(coord, N_pix, cell_size, central_size):
    """
    Along one axis: every cell containing each galaxy, as MetaDetectShear lays
    them out, and the galaxy position within it (zero-based cell pixels, as
    sx_row/sx_col).

    The image is padded by pad = (cell_size - central_size)//2 in front, and cell
    i covers the padded pixels [i*central_size, i*central_size + cell_size).
    """
    pad = (cell_size - central_size) // 2
    N_cells = int(np.ceil(N_pix / central_size))
    i_last = np.floor((coord + pad) / central_size).astype(int)
    out = []
    for k in range(int(np.ceil(cell_size / central_size)) + 1):
        i_cell = i_last - k
        pos = coord + pad - i_cell * central_size
        ok = (pos >= 0) & (pos < cell_size) & (i_cell >= 0) & (i_cell < N_cells)
        out.append((pos, ok))
    return out

def stack_one_layout(hist, row, col, shape, cell_size, central_size):
    """Adds the galaxies of one tile to the stacked histogram over the cell."""
    idx = []
    for pos_r, ok_r in cell_positions(row, shape[0], cell_size, central_size):
        for pos_c, ok_c in cell_positions(col, shape[1], cell_size, central_size):
            mask = ok_r & ok_c
            idx.append(np.floor(pos_r[mask]).astype(int) * cell_size + np.floor(pos_c[mask]).astype(int))
    ## one bincount per tile: a galaxy is in up to (cell/central)^2 cells
    hist += np.bincount(np.concatenate(idx), minlength=cell_size**2).reshape(cell_size, cell_size)

def central_slice(cell_size, central_size):
    lo = (cell_size - central_size) // 2
    return slice(lo, lo + central_size)

def check_layout(layouts, tile_labels):
    """The rebuilt layout against the image headers and the ImSim input catalogues."""
    print('>>> Checking the rebuilt layout against the ImSim outputs...')
    N_bad = 0
    for tag in A.shear_inputs:
        for i_tile, tile_label in enumerate(tile_labels):
            row, col, shape = layouts[(tag, tile_label)]
            inpath_image = os.path.join(imsim_dir, tag, 'images', 'original',
                                        f'tile{tile_label}_band{BAND}_rot0.fits')
            wcs = galsim.FitsWCS(inpath_image)
            header = galsim.FitsHeader(inpath_image)
            if shape != (header['NAXIS2'], header['NAXIS1']):
                N_bad += 1
                print(f'    {tag} tile{tile_label}: image shape {shape} rebuilt, '
                      f'{(header["NAXIS2"], header["NAXIS1"])} in the image')
            ## positions, for the first tile
            if i_tile == 0:
                gals = pd.read_feather(os.path.join(imsim_dir, tag, 'catalogues', 'input',
                                                    f'gals_info_tile{tile_label}.feather'),
                                       columns=['RA_input', 'DEC_input'])
                x, y = wcs.toImage(gals['RA_input'].values, gals['DEC_input'].values,
                                   units=galsim.degrees)
                ## FITS pixels are one-based
                dist, _ = cKDTree(np.column_stack([row, col])).query(np.column_stack([y - 1, x - 1]))
                print(f'    {tag} tile{tile_label}: {len(gals)} input galaxies, distance to the '
                      f'nearest rebuilt one: median {np.median(dist):.2e}, max {dist.max():.2e} pixels')
    print(f'    image shape differs for {N_bad} of {len(tile_labels) * len(A.shear_inputs)} images')

def rebin(hist, factor):
    n = hist.shape[0] // factor
    return hist[:n * factor, :n * factor].reshape(n, factor, n, factor).sum(axis=(1, 3))

def draw_central(ax, cell_size, central_size):
    """The central region (dashed) and the point metacal shears about (cross)."""
    lo, hi = (cell_size - central_size) / 2., (cell_size + central_size) / 2.
    ax.plot([lo, hi, hi, lo, lo], [lo, lo, hi, hi, lo], color='black', ls='--', lw=1)
    ax.plot((cell_size - 1) / 2., (cell_size - 1) / 2., marker='+', color='black', ms=8, mew=1)

def plot_cells(hists, N_tiles):
    """The whole cell, in bins of REBIN pixels."""
    fig, axs = plt.subplots(len(SETUPS), len(COLUMNS),
                            figsize=(3.4 * len(COLUMNS) + 1, 3.3 * len(SETUPS)),
                            layout='constrained', squeeze=False)
    for i_row, (cell_size, central_size) in enumerate(SETUPS):
        for i_col, (title, _) in enumerate(COLUMNS):
            hist = rebin(hists[(cell_size, central_size, i_col)], REBIN)
            ax = axs[i_row, i_col]
            im = ax.imshow(hist / hist.mean(), origin='lower', extent=(0, cell_size, 0, cell_size),
                           cmap=CMAP_DIV, vmin=0.8, vmax=1.2, interpolation='nearest')
            draw_central(ax, cell_size, central_size)
            ax.tick_params(labelsize=8)
            if i_row == 0:
                ax.set_title(title)
            if i_col == 0:
                ax.set_ylabel(f'cell {cell_size} / central {central_size}\nrow in the cell [pixels]',
                              fontsize=9)
            if i_row == len(SETUPS) - 1:
                ax.set_xlabel('column in the cell [pixels]', fontsize=9)
    fig.colorbar(im, ax=axs, location='bottom', shrink=0.6, aspect=40, extend='both',
                 label=f'galaxies per {REBIN}x{REBIN}-pixel bin / mean (all cells stacked)')
    fig.suptitle(f'grid galaxies within the cells, {N_tiles} tiles stacked')
    plt.savefig(outpath_png, dpi=200)
    plt.close()
    print(f'plot saved in {outpath_png}')

def plot_zoom(hists, N_tiles):
    """The centre of the cell in single pixels, zero-shear tag."""
    i_col = [i for i, (_, tags) in enumerate(COLUMNS) if tags == TAG_ZERO][0]
    n_cols = int(np.ceil(len(SETUPS) / 2))
    fig, axs = plt.subplots(2, n_cols, figsize=(3.0 * n_cols, 6.6), layout='constrained',
                            squeeze=False)
    for ax, (cell_size, central_size) in zip(axs.ravel(), SETUPS):
        hist = hists[(cell_size, central_size, i_col)]
        im = ax.imshow(hist / hist.mean(), origin='lower', extent=(0, cell_size, 0, cell_size),
                       cmap=CMAP_DIV, vmin=0, vmax=2, interpolation='nearest')
        draw_central(ax, cell_size, central_size)
        cen = (cell_size - 1) / 2.
        ax.set_xlim(cen - ZOOM_HALF, cen + ZOOM_HALF)
        ax.set_ylim(cen - ZOOM_HALF, cen + ZOOM_HALF)
        ax.set_title(f'cell {cell_size} / central {central_size}', fontsize=10)
        ax.tick_params(labelsize=8)
    for ax in axs.ravel()[len(SETUPS):]:
        ax.axis('off')
    fig.colorbar(im, ax=axs, location='bottom', shrink=0.6, aspect=40, extend='max',
                 label='galaxies per pixel / mean (all cells stacked)')
    fig.suptitle(f'zero-shear tag, centre of the cell in single pixels, {N_tiles} tiles stacked')
    plt.savefig(outpath_zoom_png, dpi=200)
    plt.close()
    print(f'plot saved in {outpath_zoom_png}')

if __name__ == '__main__':
    ## the tiles, and their number of galaxies, which sets the grid size
    input_dir = os.path.join(imsim_dir, TAG_ZERO[0], 'catalogues', 'input')
    inpaths = sorted(glob.glob(os.path.join(input_dir, 'gals_info_tile*.feather')))
    tile_labels = [re.search(r'gals_info_tile(.*)\.feather', os.path.basename(p)).group(1)
                   for p in inpaths]
    N_gal = {label: len(pd.read_feather(p, columns=['RA_input']))
             for label, p in zip(tile_labels, inpaths)}
    print(f'>>> {len(tile_labels)} tiles, {min(N_gal.values())} to {max(N_gal.values())} galaxies each')

    ## the layout of every tile and shear tag
    layouts = {(tag, tile_label): grid_layout(N_gal[tile_label], g)
               for tag, g in A.shear_inputs.items() for tile_label in tile_labels}
    check_layout(layouts, tile_labels)

    ## stack every set-up and column over all tiles
    hists = {(cell_size, central_size, i_col): np.zeros((cell_size, cell_size))
             for cell_size, central_size in SETUPS for i_col in range(len(COLUMNS))}
    N_central = {key: 0 for key in hists}
    N_gal_tot = {i_col: 0 for i_col in range(len(COLUMNS))}
    for i_col, (_, tags) in enumerate(COLUMNS):
        for tag in tags:
            for tile_label in tile_labels:
                row, col, shape = layouts[(tag, tile_label)]
                N_gal_tot[i_col] += len(row)
                for cell_size, central_size in SETUPS:
                    stack_one_layout(hists[(cell_size, central_size, i_col)],
                                     row, col, shape, cell_size, central_size)

    ## every galaxy has to be in exactly one central region, and how evenly the
    ##    galaxies fill it
    print('\n>>> setup, column: central-region entries / galaxies (must be 1), '
          'fraction of central-region pixels with no galaxy, '
          f'RMS of the density in {REBIN}x{REBIN}-pixel bins over the central region '
          '(and from Poisson noise alone)')
    for (cell_size, central_size, i_col), hist in hists.items():
        central = hist[central_slice(cell_size, central_size), central_slice(cell_size, central_size)]
        central_rebin = rebin(central, REBIN)
        rms = np.std(central_rebin / central_rebin.mean())
        rms_poisson = 1. / np.sqrt(central_rebin.mean())
        print(f'    cell {cell_size} / central {central_size:3d}, {"+".join(COLUMNS[i_col][1]):9s}: '
              f'{central.sum() / N_gal_tot[i_col]:.6f}, {np.mean(central == 0):.3f}, '
              f'{rms:.4f} ({rms_poisson:.4f})')

    plot_cells(hists, len(tile_labels))
    plot_zoom(hists, len(tile_labels))
