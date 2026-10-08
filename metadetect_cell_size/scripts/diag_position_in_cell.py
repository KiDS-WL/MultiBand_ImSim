# -*- coding: utf-8 -*-

### Where in the cell does the shear bias come from?
###
###    Splits the measurement of every cell set-up by where each detection sits
###    in its cell, relative to the cell centre, which is also (to half a pixel)
###    the point metacal shears the cell about:
###      - R11, R22, m1, m2 in bins of distance from the cell centre
###      - maps over the cell of the weighted sums of e1, e2 and of the number of
###        detections, for every shear tag and metacal type
###
###    The positions are the trimming positions (sx_row_noshear, sx_col_noshear,
###    in pixels of the cell image). Every cut is applied to each metacal
###    catalogue on its own positions, so R includes the selection response of the
###    cut. Weights and cuts are those of analyse_cell_size.py.
###
###    The centre is that of the central region, cell_size/2 in the zero-based
###    pixel coordinates of sx_row/sx_col, so that the kept detections are exactly
###    those with |dx|, |dy| < central_size/2. Metacal shears about the true centre
###    of the cell image, (cell_size-1)/2, half a pixel away.
###
###    The fit uses the sheared tags only. With the shears placed symmetrically
###    about zero, the zero-shear tag adds nothing to the slope, while its per-bin
###    response is noisy (the unsheared grid puts galaxies at only a few positions
###    within the cells), which would inflate the errors. The bin 'all' (every bin
###    together) therefore reproduces m1, m2, m1_R11 and m2_R22 of
###    results_cell_size.csv, which is checked at the end.
###
###    Usage:
###       python diag_position_in_cell.py              # read the catalogues, save, plot
###       python diag_position_in_cell.py --plot-only  # plot the saved results again

import os
import re
import sys
import glob
from multiprocessing import Pool

import numpy as np
import pandas as pd
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from scipy.optimize import curve_fit

import analyse_cell_size as A

## ++++++++++++++ I/O and general setups

script_dir = os.path.dirname(os.path.abspath(__file__))
outpath_csv = os.path.join(script_dir, 'diag_position_in_cell.csv')
outpath_maps = os.path.join(script_dir, 'diag_position_in_cell_maps.npz')
outpath_radial_png = os.path.join(script_dir, 'diag_position_radial.png')
outpath_maps_png = os.path.join(script_dir, 'diag_position_maps_zero_shear.png')

## Bins in distance from the cell centre (pixels)
##    the distance is max(|dx|, |dy|), square like the central regions, so a
##    set-up with central size c fills the bins up to c/2
R_EDGES = np.array([0., 12.5, 25., 50., 100., 150., 200.])
## Maps over the cell (pixels from the cell centre)
MAP_STEP = 5.
MAP_EDGES = np.arange(-200., 200. + MAP_STEP, MAP_STEP)
## map columns: weighted sums for the objects with a non-zero weight, and the
##    number of all detections
MAP_COLUMNS = ['w', 'w_e1', 'w_e2', 'N_all']

TYPES = ('noshear', '1p', '1m', '2p', '2m')
TAGS = list(A.shear_inputs)
## catalogues per task
CHUNK = 20

## Plot styles: hue for the cell size, lightness for the central size
##    blue: from step 250 to 700 of the reference palette
##    orange: around the categorical orange, the lightest still >= 2:1 on white
##    as many steps as there are central sizes are interpolated between the ends
CELL_RAMPs = {250: ['#86b6ef', '#3987e5', '#1c5cab', '#0d366b'],
              500: ['#f29d7b', '#eb6834', '#994422'],
              260: ['#5fbf8f', '#1b9e77', '#0b5e46']}

def ramp_color(cell_size, i, n):
    """The i-th of n steps of the ramp of this cell size, light to dark."""
    cmap = LinearSegmentedColormap.from_list(f'ramp{cell_size}', CELL_RAMPs[cell_size])
    return cmap(i / max(n - 1, 1))
CELL_SYMBOLs = {250: 'o', 500: 's', 260: 'D'}
## diverging: blue - neutral gray - red
CMAP_DIV = LinearSegmentedColormap.from_list('div', ['#2a78d6', '#f0efec', '#e34948'])
## bins of the maps with fewer detections are left blank
MAP_MIN_N = 20

plt.rc('font', size=11)

## ++++++++++++++ Workhorse

def parse_setup(label):
    cell_size, central_size = map(int, re.search(r'cell(\d+)_central(\d+)', label).groups())
    return cell_size, central_size

def chunk_sums(args):
    """
    Sums for a chunk of catalogues of one set-up and shear tag.

    Returns, per catalogue, the radial sums (type, bin, [w, w*g1, w*g2, N, w*r]),
    and the maps summed over the chunk (type, y, x, MAP_COLUMNS).
    """
    inpaths, cell_size = args
    ## the centre of the central region (see the header)
    cen = cell_size / 2.
    n_r = len(R_EDGES) - 1
    n_map = len(MAP_EDGES) - 1

    radial = {}
    maps = np.zeros((len(TYPES), n_map, n_map, len(MAP_COLUMNS)))
    for inpath in inpaths:
        cata = pd.read_feather(inpath, columns=A.weight_columns + ['sx_row_noshear',
                                                                   'sx_col_noshear'])
        weight = A.shear_weight(cata)
        g1 = cata[f'{A.fit_model}_g_1'].values
        g2 = cata[f'{A.fit_model}_g_2'].values
        dx = cata['sx_col_noshear'].values - cen
        dy = cata['sx_row_noshear'].values - cen
        r = np.maximum(np.abs(dx), np.abs(dy))
        i_r = np.digitize(r, R_EDGES) - 1
        shear_type = cata['shear_type'].values

        radial_tmp = np.zeros((len(TYPES), n_r, 5))
        for i_t, t in enumerate(TYPES):
            mask_t = (shear_type == t)
            mask_w = mask_t & (weight > 0)

            mask = mask_w & (i_r < n_r)
            w = weight[mask]
            for i_col, val in enumerate([w, w * g1[mask], w * g2[mask],
                                         np.ones_like(w), w * r[mask]]):
                radial_tmp[i_t, :, i_col] = np.bincount(i_r[mask], weights=val, minlength=n_r)

            for i_col, (mask_tmp, val) in enumerate([(mask_w, weight),
                                                     (mask_w, weight * g1),
                                                     (mask_w, weight * g2),
                                                     (mask_t, None)]):
                maps[i_t, :, :, i_col] += np.histogram2d(
                    dy[mask_tmp], dx[mask_tmp], bins=[MAP_EDGES, MAP_EDGES],
                    weights=None if val is None else val[mask_tmp])[0]
        radial[inpath] = radial_tmp
    return radial, maps

def bins_one_setup(label, radial_setup):
    """R and m for every distance bin of one set-up, and for all bins together."""
    n_r = len(R_EDGES) - 1
    tags_sheared = [t for t in TAGS if A.shear_inputs[t] != (0., 0.)]
    w_tot = sum(v[0, :, 0].sum() for t in tags_sheared for v in radial_setup[t].values())

    rows = []
    for k in list(range(n_r)) + ['all']:
        ks = list(range(n_r)) if k == 'all' else [k]
        g_in = {1: [], 2: []}
        g_out = {name: {1: [], 2: []} for name in ('common', 'own')}
        R11_list, R22_list = [], []
        w_bin = wr_bin = 0.
        for tag in tags_sheared:
            files = radial_setup[tag]
            tot = sum(v[:, ks, :].sum(axis=1) for v in files.values())
            if np.any(tot[:, 0] <= 0):
                break
            ## Calculate Response
            mean_g = tot[:, 1:3] / tot[:, [0]]
            R11 = (mean_g[1, 0] - mean_g[2, 0]) / 0.02
            R22 = (mean_g[3, 1] - mean_g[4, 1]) / 0.02
            R = (R11 + R22) / 2
            R11_list.append(R11)
            R22_list.append(R22)
            w_bin += tot[0, 0]
            wr_bin += tot[0, 4]

            ## measured shear for each tile, both rotations together
            tile_sums = {}
            for inpath, v in files.items():
                tile_label = re.search(r'tile(.*)_rot', os.path.basename(inpath)).group(1)
                tile_sums[tile_label] = tile_sums.get(tile_label, 0.) + v[0, ks, :3].sum(axis=0)
            tile_sums = np.array([s for s in tile_sums.values() if s[0] > 0])
            e1 = tile_sums[:, 1] / tile_sums[:, 0]
            e2 = tile_sums[:, 2] / tile_sums[:, 0]
            g1_input, g2_input = A.shear_inputs[tag]
            g_in[1].append(np.full(len(e1), g1_input))
            g_in[2].append(np.full(len(e2), g2_input))
            g_out['common'][1].append(e1 / R)
            g_out['common'][2].append(e2 / R)
            g_out['own'][1].append(e1 / R11)
            g_out['own'][2].append(e2 / R22)
        else:
            row = {'setup': label,
                   'r_bin': 'all' if k == 'all' else f'{R_EDGES[k]:g}-{R_EDGES[k+1]:g}',
                   'r_lo': R_EDGES[ks[0]], 'r_hi': R_EDGES[ks[-1] + 1],
                   'r_mean': wr_bin / w_bin, 'wfrac': w_bin / w_tot,
                   ## averaged over the sheared tags
                   'R11': np.mean(R11_list), 'R22': np.mean(R22_list)}
            for name in ('common', 'own'):
                for comp in (1, 2):
                    popt, pcov = curve_fit(A.line_func, np.concatenate(g_in[comp]),
                                           np.concatenate(g_out[name][comp]))
                    key = f'm{comp}' if name == 'common' else f'm{comp}_R{comp}{comp}'
                    row[key] = popt[0]
                    row[f'{key}_err'] = np.sqrt(pcov[0, 0])
            rows.append(row)
    return rows

def compute():
    n_map = len(MAP_EDGES) - 1
    tasks, task_keys, labels = [], [], []
    for label, (main_dir_tmp, shape_folder) in A.setups.items():
        cell_size, _ = parse_setup(label)
        inpath_dict = {tag: sorted(glob.glob(os.path.join(main_dir_tmp, tag, 'catalogues',
                                                          shape_folder, '*.feather')))
                       for tag in TAGS}
        N_files = {tag: len(v) for tag, v in inpath_dict.items()}
        print(f'>>> {label}: catalogues found {N_files}')
        if min(N_files.values()) == 0:
            print('>>> Some shear tags have no catalogue, skipped!')
            continue
        if any(N != A.N_files_expected for N in N_files.values()):
            print(f'>>> WARNING: incomplete, expected {A.N_files_expected} catalogues per shear tag!')
        labels.append(label)
        for tag, inpaths in inpath_dict.items():
            for i in range(0, len(inpaths), CHUNK):
                tasks.append((inpaths[i:i + CHUNK], cell_size))
                task_keys.append((label, tag))

    radial = {label: {tag: {} for tag in TAGS} for label in labels}
    maps = {label: np.zeros((len(TAGS), len(TYPES), n_map, n_map, len(MAP_COLUMNS)))
            for label in labels}
    N_proc = int(os.environ.get('SLURM_CPUS_PER_TASK', 4))
    with Pool(N_proc) as pool:
        for i_task, ((label, tag), (radial_chunk, maps_chunk)) in enumerate(
                zip(task_keys, pool.imap(chunk_sums, tasks))):
            radial[label][tag].update(radial_chunk)
            maps[label][TAGS.index(tag)] += maps_chunk
            if (i_task + 1) % 50 == 0:
                print(f'>>> {i_task + 1}/{len(tasks)} chunks done', flush=True)

    results = pd.DataFrame([row for label in labels
                            for row in bins_one_setup(label, radial[label])])
    results.to_csv(outpath_csv, index=False)
    print(f'saved to {outpath_csv}')
    np.savez_compressed(outpath_maps, map_edges=MAP_EDGES, tags=np.array(TAGS),
                        types=np.array(TYPES), columns=np.array(MAP_COLUMNS),
                        **{f'maps_{label}': maps[label] for label in labels})
    print(f'saved to {outpath_maps}')
    return results

def check_against_summary(results):
    """The bin 'all' must reproduce results_cell_size.csv."""
    if not os.path.isfile(A.outpath_csv):
        return
    summary = pd.read_csv(A.outpath_csv).set_index('setup')
    all_bins = results[results['r_bin'] == 'all'].set_index('setup')
    cols = [c for c in ('m1', 'm2', 'm1_R11', 'm2_R22') if c in summary.columns]
    diff = (all_bins[cols] - summary.loc[all_bins.index, cols]).abs()
    print(f'\n>>> max |difference| from {os.path.basename(A.outpath_csv)} for the bin "all":')
    print(diff.max().to_string(float_format='%.1e'))

def plot_radial(results):
    data = results[results['r_bin'] != 'all']
    setups = sorted(data['setup'].unique(), key=parse_setup)
    fig, axs = plt.subplots(1, 2, figsize=(11, 4.5), sharex=True, sharey=True)
    for ax, comp, title in zip(axs, ['m1_R11', 'm2_R22'], [r'$m_1$ (with $R_{11}$)',
                                                           r'$m_2$ (with $R_{22}$)']):
        for i_setup, setup in enumerate(setups):
            cell_size, central_size = parse_setup(setup)
            same_cell = [s for s in setups if parse_setup(s)[0] == cell_size]
            color = ramp_color(cell_size, same_cell.index(setup), len(same_cell))
            data_tmp = data[data['setup'] == setup]
            ## small shifts, so that the error bars of the same bin do not overlap
            x_val = data_tmp['r_mean'].values * (1 + 0.02 * (i_setup - (len(setups) - 1) / 2))
            ax.errorbar(x_val, data_tmp[comp].values, yerr=data_tmp[f'{comp}_err'].values,
                        color=color, marker=CELL_SYMBOLs[cell_size], markersize=6,
                        mfc=color if cell_size == 250 else 'white', elinewidth=1, lw=1,
                        label=f'{cell_size} / {central_size}')
        ax.axhline(y=0, color='gray', ls='--', lw=1)
        ax.set_xscale('log')
        ax.xaxis.set_major_formatter(mpl.ticker.FormatStrFormatter('%g'))
        ax.set_xlabel('distance from the cell centre, max(|dx|, |dy|) [pixels]')
        ax.set_title(title)
    axs[0].set_ylabel(r'$m$')
    axs[1].legend(title='cell / central size', frameon=False, fontsize=9,
                  loc='center left', bbox_to_anchor=(1.01, 0.5))
    plt.tight_layout()
    plt.savefig(outpath_radial_png, dpi=200)
    plt.close()
    print(f'plot saved in {outpath_radial_png}')

def plot_maps(maps_file):
    """Zero-shear tag, noshear catalogues: mean e1, e2 and the detection density."""
    maps_all = np.load(maps_file)
    edges = maps_all['map_edges']
    tags = list(maps_all['tags'])
    i_zero = [i for i, t in enumerate(tags) if A.shear_inputs[t] == (0., 0.)][0]
    setups = sorted([k[len('maps_'):] for k in maps_all.files if k.startswith('maps_')],
                    key=parse_setup)

    ## mean e1, e2 and the density relative to its mean, blank where too few detections
    panels = []
    for setup in setups:
        _, central_size = parse_setup(setup)
        m = maps_all[f'maps_{setup}'][i_zero, TYPES.index('noshear')]
        with np.errstate(invalid='ignore', divide='ignore'):
            e1 = m[..., 1] / m[..., 0]
            e2 = m[..., 2] / m[..., 0]
        N = m[..., 3]
        cen = 0.5 * (edges[1:] + edges[:-1])
        inside = (np.abs(cen)[None, :] < central_size / 2.) & (np.abs(cen)[:, None] < central_size / 2.)
        blank = (N < MAP_MIN_N) | (~inside)
        density = N / np.mean(N[inside])
        for arr in (e1, e2, density):
            arr[blank] = np.nan
        panels.append((setup, central_size, e1, e2, density))

    ## shared colour limits per column
    lim_e = np.nanpercentile(np.abs(np.concatenate([np.ravel(p[2:4]) for p in panels])), 98)
    lim_d = np.nanpercentile(np.abs(np.concatenate([np.ravel(p[4]) for p in panels]) - 1), 98)

    fig, axs = plt.subplots(len(panels), 3, figsize=(9, 2.7 * len(panels)), layout='constrained')
    axs = np.atleast_2d(axs)
    for i_row, (setup, central_size, e1, e2, density) in enumerate(panels):
        cell_size, _ = parse_setup(setup)
        for i_col, (arr, vmin, vmax, title) in enumerate([
                (e1, -lim_e, lim_e, r'$\langle e_1 \rangle$'),
                (e2, -lim_e, lim_e, r'$\langle e_2 \rangle$'),
                (density, 1 - lim_d, 1 + lim_d, 'detections / mean')]):
            ax = axs[i_row, i_col]
            im = ax.pcolormesh(edges, edges, arr, cmap=CMAP_DIV, vmin=vmin, vmax=vmax)
            half = central_size / 2.
            ax.set_xlim(-half, half)
            ax.set_ylim(-half, half)
            ax.set_aspect('equal')
            ax.tick_params(labelsize=8)
            if i_row == 0:
                ax.set_title(title)
            if i_col == 0:
                ax.set_ylabel(f'cell {cell_size} / central {central_size}\ndy [pixels]', fontsize=9)
            if i_row == len(panels) - 1:
                ax.set_xlabel('dx [pixels]', fontsize=9)
            if i_row == len(panels) - 1:
                fig.colorbar(im, ax=axs[:, i_col], location='bottom', shrink=0.8, aspect=30)
    fig.suptitle('zero-shear tag, noshear catalogues, relative to the cell centre')
    plt.savefig(outpath_maps_png, dpi=200)
    plt.close()
    print(f'plot saved in {outpath_maps_png}')

if __name__ == '__main__':
    if '--plot-only' in sys.argv:
        results = pd.read_csv(outpath_csv)
    else:
        results = compute()
    check_against_summary(results)
    with pd.option_context('display.width', 250, 'display.max_columns', None):
        print(results[['setup', 'r_bin', 'r_mean', 'wfrac', 'R11', 'R22',
                       'm1_R11', 'm1_R11_err', 'm2_R22', 'm2_R22_err']].to_string(
                           index=False, float_format='%.4f'))
    plot_radial(results)
    plot_maps(outpath_maps)
