# -*- coding: utf-8 -*-

### Does the separation of the grid galaxies matter?
###
###    Uses the existing catalogues, nothing is re-run. Every detection is matched
###    to its input galaxy, which sits at grid position (i, j) of the 9-arcsec grid,
###    and the shear is measured on sub-grids that keep every second galaxy in both
###    directions, i.e. galaxies 18 arcsec apart:
###       ij00  i even, j even       (the selection shown in the plot)
###       ij01, ij10, ij11           the other three such sub-grids
###    together with
###       all       every detection, the analysis of analyse_cell_size.py (a check)
###       matched   every detection matched to an input galaxy
###
###    NOTE: only the measured sample changes, not the images. Every selected galaxy
###          still has its 9-arcsec neighbours around it in detection, metacal and
###          measurement. So this tests whether the bias is carried by particular
###          galaxies of the grid, not whether physically sparser images (e.g. a
###          15-arcsec grid) behave differently.
###
###    The selection is by the identity of the galaxy, applied in the same way to
###    every metacal catalogue, so it does not respond to the shear. Detections are
###    matched by their sky position with the metacal shift removed (X_WORLD,
###    Y_WORLD) to the input positions of the same shear tag; a galaxy sits at the
###    same grid point in every tag. Weights, cuts, response and fit are those of
###    analyse_cell_size.py.
###
###    Output: results_grid_separation.csv and m_grid_separation.png
###    Usage:  python test_grid_separation.py              # read the catalogues, save, plot
###            python test_grid_separation.py --plot-only  # plot the saved results again

import os
import re
import sys
import glob
from functools import lru_cache
from multiprocessing import Pool

import numpy as np
import pandas as pd
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from scipy.spatial import cKDTree

import analyse_cell_size as A

## ++++++++++++++ I/O and general setups

## ImSim input catalogues of the images (the same for every metadetect set-up)
input_dir = '/sdf/data/kipac/u/jlkitt/imsims/grid/seeing_073/{tag}/catalogues/input/gals_info_tile{tile}.feather'
GRID_SIZE = 9. / 3600.     # degrees between neighbouring grid galaxies
## a detection belongs to the input galaxy within this distance (arcsec)
R_MATCH = 2.

SELECTIONS = ['all', 'matched', 'ij00', 'ij01', 'ij10', 'ij11']
## the selection compared with the full sample in the plot
SELECTION_PLOT = 'ij00'
TYPES = ('noshear', '1p', '1m', '2p', '2m')

script_dir = os.path.dirname(os.path.abspath(__file__))
outpath_csv = os.path.join(script_dir, 'results_grid_separation.csv')
outpath_png = os.path.join(script_dir, 'm_grid_separation.png')

## ++++++++++++++ Workhorse

def wrap_ra(ra):
    """RA in (-180, 180] degrees: sheared grid positions can fall just below 0."""
    return (np.asarray(ra) + 180.) % 360. - 180.

@lru_cache(maxsize=16)
def input_galaxies(tag, tile):
    """KD-tree of the input positions of one tile and shear tag, and the grid (i, j) of each."""
    gals = pd.read_feather(input_dir.format(tag=tag, tile=tile),
                           columns=['index_input', 'RA_input', 'DEC_input'])
    ## the grid position of every galaxy, from the unsheared (zero-shear) positions
    zero = [t for t, g in A.shear_inputs.items() if g == (0., 0.)][0]
    grid = pd.read_feather(input_dir.format(tag=zero, tile=tile),
                           columns=['index_input', 'RA_input', 'DEC_input'])
    grid_i = pd.Series(np.round(grid['RA_input'].values / GRID_SIZE).astype(int), index=grid['index_input'])
    grid_j = pd.Series(np.round(grid['DEC_input'].values / GRID_SIZE).astype(int), index=grid['index_input'])
    tree = cKDTree(np.column_stack([wrap_ra(gals['RA_input']), gals['DEC_input']]) * 3600.)
    return (tree,
            grid_i.loc[gals['index_input']].values,
            grid_j.loc[gals['index_input']].values)

def file_sums_selected(args):
    """The sums of A.file_sums for every selection, for one catalogue."""
    inpath, tag = args
    tile = re.search(r'tile(.*)_band', os.path.basename(inpath)).group(1)
    cata = pd.read_feather(inpath, columns=A.weight_columns + ['X_WORLD', 'Y_WORLD'])
    g1 = cata[f'{A.fit_model}_g_1'].values
    g2 = cata[f'{A.fit_model}_g_2'].values
    weight = A.shear_weight(cata)
    shear_type = cata['shear_type'].values

    ## the input galaxy of every detection, and its grid position
    tree, grid_i, grid_j = input_galaxies(tag, tile)
    dist, idx = tree.query(np.column_stack([wrap_ra(cata['X_WORLD']), cata['Y_WORLD']]) * 3600.)
    matched = dist < R_MATCH
    i_mod = np.where(matched, grid_i[idx] % 2, -1)
    j_mod = np.where(matched, grid_j[idx] % 2, -1)
    masks = {'all': np.ones(len(cata), dtype=bool), 'matched': matched}
    for i in (0, 1):
        for j in (0, 1):
            masks[f'ij{i}{j}'] = (i_mod == i) & (j_mod == j)

    out = {}
    for sel, mask_sel in masks.items():
        sums = {}
        for t in TYPES:
            mask = mask_sel & (weight > 0) & (shear_type == t)
            w = weight[mask]
            sums[t] = np.array([w.sum(), (w * g1[mask]).sum(), (w * g2[mask]).sum(), mask.sum()])
        out[sel] = sums
    return inpath, out

def compute():
    results = []
    N_proc = int(os.environ.get('SLURM_CPUS_PER_TASK', 4))
    with Pool(N_proc) as pool:
        for label, (main_dir_tmp, shape_folder) in A.setups.items():
            inpath_dict = {tag: sorted(glob.glob(os.path.join(main_dir_tmp, tag, 'catalogues',
                                                              shape_folder, '*.feather')))
                           for tag in A.shear_inputs}
            N_files = {tag: len(v) for tag, v in inpath_dict.items()}
            print(f'\n>>>>>>>>>>>>>> {label}: catalogues found {N_files}', flush=True)
            if min(N_files.values()) == 0:
                print('>>> Some shear tags have no catalogue, skipped!')
                continue
            jobs = [(p, tag) for tag, v in inpath_dict.items() for p in v]
            per_file = dict(pool.imap_unordered(file_sums_selected, jobs, chunksize=4))
            for sel in SELECTIONS:
                res = A.bias_from_sums(label, inpath_dict, {p: v[sel] for p, v in per_file.items()})
                res['selection'] = sel
                results.append(res)
    results = pd.DataFrame(results)
    results.to_csv(outpath_csv, index=False)
    print(f'saved to {outpath_csv}')
    return results

def check(results):
    """'all' must reproduce results_cell_size.csv; how many detections are matched."""
    summary = pd.read_csv(A.outpath_csv).set_index('setup')
    res_all = results[results['selection'] == 'all'].set_index('setup')
    cols = ['m1', 'm2', 'm1_R11', 'm2_R22']
    diff = (res_all[cols] - summary.loc[res_all.index, cols]).abs().max()
    print('\n>>> selection "all" against results_cell_size.csv, max |difference|:')
    print(diff.to_string(float_format='%.1e'))
    print('\n>>> objects with non-zero weight kept by each selection, relative to "all":')
    n = results.pivot_table(index='setup', columns='selection', values='N_obj_weighted')
    print((n[SELECTIONS].div(n['all'], axis=0)).to_string(float_format='%.3f'))

def plot(results):
    plt.rcParams["text.usetex"] = True
    mpl.rcParams['xtick.direction'] = 'in'
    mpl.rcParams['ytick.direction'] = 'in'
    mpl.rcParams['xtick.top'] = True
    mpl.rcParams['ytick.right'] = True
    plt.rc('font', size=14, family='serif')

    data = results[results['selection'].isin(['all', SELECTION_PLOT])].copy()
    cell_central = data['setup'].str.extract(r'cell(\d+)_central(\d+)').astype(int)
    data['cell_size'], data['central_size'] = cell_central[0], cell_central[1]
    setups = (data[['setup', 'cell_size', 'central_size', 'complete', 'N_files']]
              .drop_duplicates('setup').sort_values(['cell_size', 'central_size']).reset_index(drop=True))
    cell_sizes = sorted(setups['cell_size'].unique())
    i_group = setups['cell_size'].map({c: i for i, c in enumerate(cell_sizes)}).values
    x_cen = np.arange(len(setups)) + 0.6 * i_group
    x_of = dict(zip(setups['setup'], x_cen))

    fig, ax = plt.subplots(figsize=(8, 5))
    for x0, setup in zip(x_cen, setups['setup']):
        if setup == 'cell250_central150':
            ax.axvspan(x0 - 0.4, x0 + 0.4, color='#f0efec', zorder=0)

    ## m1 and m2 with their own response, all galaxies (filled) and every second one (open)
    series = [('m1_R11', '#2a78d6', 'o', r'$m_1$'), ('m2_R22', '#eb6834', 's', r'$m_2$')]
    offsets = {('m1_R11', 'all'): -0.27, ('m1_R11', SELECTION_PLOT): -0.09,
               ('m2_R22', 'all'): 0.09, ('m2_R22', SELECTION_PLOT): 0.27}
    for comp, color, marker, label in series:
        for sel, fill, sel_label in (('all', True, r'all galaxies (9$^{\prime\prime}$ apart)'),
                                     (SELECTION_PLOT, False, r'every second galaxy (18$^{\prime\prime}$ apart)')):
            d = data[data['selection'] == sel]
            ax.errorbar(d['setup'].map(x_of) + offsets[(comp, sel)], d[comp], yerr=d[f'{comp}_err'],
                        color=color, marker=marker, markersize=7, mfc=color if fill else 'white',
                        elinewidth=1, ls='none', label=f'{label}, {sel_label}')
    ax.axhline(y=0, color='gray', ls='--', lw=1)

    for i, cell_size in enumerate(cell_sizes):
        mask = (i_group == i)
        if i > 0:
            ax.axvline(x=(x_cen[mask][0] + x_cen[i_group == i - 1][-1]) / 2., color='black', lw=1)
        ax.text(x_cen[mask].mean(), 1.02, f'cell size = {cell_size}',
                transform=ax.get_xaxis_transform(), ha='center', va='bottom')

    xticklabels, notes = [], []
    N_full = setups.loc[setups['complete'], 'N_files'].max()
    for _, row in setups.iterrows():
        label = f"{row['central_size']}"
        if not row['complete']:
            label += r'$^*$'
            notes.append(f"cell {row['cell_size']} / central {row['central_size']} "
                         f"({row['N_files']} of {N_full} catalogues)")
        if row['setup'] == 'cell250_central150':
            label += '\n(baseline)'
        xticklabels.append(label)
    ax.set_xticks(x_cen)
    ax.set_xticklabels(xticklabels)
    ax.tick_params(axis='x', length=0)
    ax.set_xlim(x_cen[0] - 0.6, x_cen[-1] + 0.6)
    xlabel = r'central size [pixels]'
    if notes:
        xlabel += '\n' + r'{\footnotesize $^*$incomplete run: ' + ', '.join(notes) + '}'
    ax.set_xlabel(xlabel)
    ax.set_ylabel(r'$m$ (each component with its own response)')
    ax.legend(frameon=False, loc='lower right', handletextpad=0.2, fontsize=10)
    plt.tight_layout()
    plt.savefig(outpath_png, dpi=300)
    plt.close()
    print(f'plot saved in {outpath_png}')

if __name__ == '__main__':
    if '--plot-only' in sys.argv:
        results = pd.read_csv(outpath_csv)
    else:
        results = compute()
    check(results)
    with pd.option_context('display.width', 200):
        print(results.pivot_table(index='setup', columns='selection', values='m2_R22')[SELECTIONS]
              .to_string(float_format='%+.4f'))
    plot(results)
