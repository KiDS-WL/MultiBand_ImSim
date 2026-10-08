# -*- coding: utf-8 -*-

### Do cell sizes 500 and 250 see the galaxies in the same places?
###
###    The same stacked map of galaxy positions within the cell as
###    plot_cell_galaxy_positions.py, for pairs of set-ups with the same central
###    size and cell sizes 500 and 250, side by side, one row per shear tag.
###
###    All panels share one pixel scale, centred on the central region, so that the
###    central regions (dashed) have the same size and the patterns can be compared
###    pixel by pixel. A 250-pixel cell only fills the inner half of its panel.
###    The colour is the density relative to its mean inside the central region.
###
###    The central regions tile the image with a step of central_size from pixel 0
###    whatever the cell size, so inside them the two cell sizes must see exactly
###    the same galaxies at the same positions; this is checked at the end.
###
###    Usage: python plot_cell_galaxy_positions_pairs.py

import os
import re
import glob

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

import analyse_cell_size as A
import plot_cell_galaxy_positions as P

## ++++++++++++++ I/O and general setups

outpath_png = os.path.join(P.script_dir, 'cell_galaxy_positions_pairs.png')

## (cell_size, central_size), side by side in this order
SETUPS = [(500, 50), (250, 50), (500, 100), (250, 100), (500, 150), (250, 150)]

## the rows: one per shear tag, the zero-shear one first
TAGS = [tag for tag, g in sorted(A.shear_inputs.items(), key=lambda kv: kv[1] != (0., 0.))]

## the half-width shown, in pixels from the centre of the central region
HALF = max(cell_size for cell_size, _ in SETUPS) / 2.

## ++++++++++++++ Workhorse

def plot(hists, N_tiles):
    fig, axs = plt.subplots(len(TAGS), len(SETUPS),
                            figsize=(2.9 * len(SETUPS) + 1, 3.0 * len(TAGS)),
                            layout='constrained', squeeze=False)
    for i_row, tag in enumerate(TAGS):
        g = A.shear_inputs[tag]
        for i_col, (cell_size, central_size) in enumerate(SETUPS):
            hist = hists[(tag, cell_size, central_size)]
            central = hist[P.central_slice(cell_size, central_size),
                           P.central_slice(cell_size, central_size)]
            hist = P.rebin(hist, P.REBIN)
            ax = axs[i_row, i_col]
            ## pixel coordinates relative to the centre of the central region
            half_cell = cell_size / 2.
            im = ax.imshow(hist / (central.mean() * P.REBIN**2), origin='lower',
                           extent=(-half_cell, half_cell, -half_cell, half_cell),
                           cmap=P.CMAP_DIV, vmin=0.8, vmax=1.2, interpolation='nearest')
            ## the central region, and the point metacal shears about
            half = central_size / 2.
            ax.plot([-half, half, half, -half, -half], [-half, -half, half, half, -half],
                    color='black', ls='--', lw=1)
            ax.plot(-0.5, -0.5, marker='+', color='black', ms=8, mew=1)
            ## the edge of the cell, when it is smaller than the panel
            if half_cell < HALF:
                ax.plot([-half_cell, half_cell, half_cell, -half_cell, -half_cell],
                        [-half_cell, -half_cell, half_cell, half_cell, -half_cell],
                        color='gray', ls='-', lw=0.8)
            ax.set_xlim(-HALF, HALF)
            ax.set_ylim(-HALF, HALF)
            ax.tick_params(labelsize=7)
            if i_row == 0:
                ax.set_title(f'cell {cell_size} / central {central_size}', fontsize=10)
            if i_col == 0:
                ax.set_ylabel(f'{tag}, g = ({g[0]:+.2f}, {g[1]:+.2f})\ndy [pixels]', fontsize=9)
            if i_row == len(TAGS) - 1:
                ax.set_xlabel('dx [pixels]', fontsize=9)
    fig.colorbar(im, ax=axs, location='bottom', shrink=0.5, aspect=40, extend='both',
                 label=f'galaxies per {P.REBIN}x{P.REBIN}-pixel bin / mean in the central region '
                       '(all cells stacked)')
    fig.suptitle(f'grid galaxies within the cells, cell size 500 next to 250, {N_tiles} tiles stacked')
    plt.savefig(outpath_png, dpi=200)
    plt.close()
    print(f'plot saved in {outpath_png}')

if __name__ == '__main__':
    ## the tiles, and their number of galaxies, which sets the grid size
    input_dir = os.path.join(P.imsim_dir, P.TAG_ZERO[0], 'catalogues', 'input')
    inpaths = sorted(glob.glob(os.path.join(input_dir, 'gals_info_tile*.feather')))
    tile_labels = [re.search(r'gals_info_tile(.*)\.feather', os.path.basename(p)).group(1)
                   for p in inpaths]
    N_gal = {label: len(pd.read_feather(p, columns=['RA_input']))
             for label, p in zip(tile_labels, inpaths)}

    hists = {(tag, cell_size, central_size): np.zeros((cell_size, cell_size))
             for tag in TAGS for cell_size, central_size in SETUPS}
    for tag in TAGS:
        for tile_label in tile_labels:
            row, col, shape = P.grid_layout(N_gal[tile_label], A.shear_inputs[tag])
            for cell_size, central_size in SETUPS:
                stack_key = (tag, cell_size, central_size)
                P.stack_one_layout(hists[stack_key], row, col, shape, cell_size, central_size)

    ## inside the central region the two cell sizes must agree pixel by pixel
    print('>>> central region, cell 500 against cell 250: max |difference| in galaxies '
          'per pixel / galaxies in the central region')
    for central_size in sorted({c for _, c in SETUPS}):
        cells = sorted({cell for cell, c in SETUPS if c == central_size})
        if len(cells) != 2:
            continue
        for tag in TAGS:
            centrals = [hists[(tag, cell, central_size)][P.central_slice(cell, central_size),
                                                          P.central_slice(cell, central_size)]
                        for cell in cells]
            print(f'    central {central_size:3d}, {tag}: '
                  f'{np.abs(centrals[0] - centrals[1]).max():.0f} / {centrals[0].sum():.0f}')

    plot(hists, len(tile_labels))
