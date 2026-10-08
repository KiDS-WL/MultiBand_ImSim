# -*- coding: utf-8 -*-

### Plots of the test-bench layout studies (test_layout.py), from their summaries:
###    m_vs_grid_spacing.png       m1, m2 against the grid spacing (coarse and fine scans)
###    m_vs_grid_spacing_fine.png  the same, zoomed on the fine scan
###    m_vs_jitter.png             9-arcsec grid with jittered / unsheared positions
###    m_vs_cell_size.png          9-arcsec grid in cells of different sizes
###    m_map_cell_spacing.png      m against cell size and grid spacing (run_layout_tests.sh map)
###    m_map_cell_spacing_500.png  the same for cells of 500 to 512 px (run_layout_tests.sh map500)
###    m_vs_table_margin.png       m against P - 4N (image table - 4 x cell size), all tests
###    m_vs_forced_fft.png         250-px cells with GalSim's FFT sizes changed (run_fft_size_test.sh)
###    A figure is skipped when its summary does not exist yet.
###
###    Galaxies within 20 px of the cell centre, noise-free scenes, cuts of
###    analyse_cell_size.py, unweighted; m_i = S_i / R_ii - 1.
###
###    Usage: python plot_layout_scans.py

import os

import numpy as np
import pandas as pd
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt

script_dir = os.path.dirname(os.path.abspath(__file__))
PIXEL_SCALE = 0.2
SPACING_SIM = 9.                 # arcsec, the grid of the simulation
COMPS = [('m1', '#2a78d6', 'o', r'$m_1$'), ('m2', '#eb6834', 's', r'$m_2$')]
NOTE = 'noise-free test bench: galaxies within 20 px of the cell centre'

plt.rc('font', size=11)

def read(tag):
    path = os.path.join(script_dir, f'layout_response{tag}.csv')
    if not os.path.isfile(path):
        print(f'no {os.path.basename(path)}, skipped')
        return None
    out = pd.read_csv(path, dtype={'cell': str})
    return out[~out['cell'].str.contains('-')]

def save(fig, name):
    path = os.path.join(script_dir, name)
    fig.savefig(path, dpi=200)
    plt.close(fig)
    print(f'plot saved in {path}')

def grid_spacing(df):
    df = df[df['layout'].str.match(r'^grid[0-9.]+$')].copy()
    df['spacing'] = df['layout'].str[len('grid'):].astype(float)
    return df

def plot_spacing(coarse, fine, xlim=None, name='m_vs_grid_spacing.png'):
    """m1 (top) and m2 (bottom) against the grid spacing, cells of 250 (filled) and 500 px (open)."""
    fig, axs = plt.subplots(2, 1, figsize=(8, 6.4), sharex=True, layout='constrained')
    for ax, (comp, color, marker, label) in zip(axs, COMPS):
        for df, size, which in ((coarse, 7, 'scan'), (fine, 4.5, 'fine scan')):
            if df is None:
                continue
            for cell, fill in (('250', True), ('500', False)):
                d = df[df['cell'] == cell].sort_values('spacing')
                if xlim is not None:
                    d = d[(d['spacing'] >= xlim[0] - 1e-6) & (d['spacing'] <= xlim[1] + 1e-6)]
                ax.errorbar(d['spacing'], d[comp], yerr=d[f'{comp}_err'], color=color, marker=marker,
                            ms=size, mfc=color if fill else 'white', ls='none', elinewidth=1,
                            label=f'cell {cell} px ({which})')
        ax.axhline(0, color='gray', ls='--', lw=1)
        ax.axvspan(SPACING_SIM - 0.04, SPACING_SIM + 0.04, color='#e1e0d9', zorder=0)
        ax.set_ylabel(label)
    ## one legend above the panels, in black, so that it covers no points
    handles = [plt.Line2D([], [], color='k', marker='o', ms=size, mfc='k' if fill else 'white', ls='none',
                          label=f'cell {cell} px ({which})')
               for df, size, which in ((coarse, 7, 'scan'), (fine, 4.5, 'fine scan')) if df is not None
               for cell, fill in (('250', True), ('500', False))]
    fig.legend(handles=handles, loc='outside lower center', ncol=len(handles), frameon=False, fontsize=8)
    axs[-1].set_xlabel('grid spacing [arcsec]   (shaded: the 9-arcsec grid of the simulation)')
    sec = axs[0].secondary_xaxis('top', functions=(lambda x: x / PIXEL_SCALE, lambda p: p * PIXEL_SCALE))
    sec.set_xlabel('grid spacing [pixels]')
    if xlim is not None:
        axs[0].set_xlim(xlim[0] - 0.05, xlim[1] + 0.05)
    fig.suptitle(NOTE, fontsize=10)
    save(fig, name)

def plot_jitter(jit, first):
    """9-arcsec grid: jittered neighbours, unsheared positions, and random positions for reference."""
    cats = []
    for layout, label in [('grid9', 'grid\n(no jitter)')] + \
            [(f'jitter9_{a:g}', f'jitter\n±{a:g} px') for a in sorted(
                float(x[len('jitter9_'):]) for x in jit['layout'].unique() if x.startswith('jitter9_'))] + \
            [('grid9_unsheared', 'grid,\npositions\nnot sheared')]:
        if layout in set(jit['layout']):
            cats.append((jit, layout, label))
    if first is not None and 'random9' in set(first['layout']):
        cats.append((first, 'random9', 'random,\nsame density'))
    x = np.arange(len(cats))
    fig, axs = plt.subplots(2, 1, figsize=(8, 6.4), sharex=True, layout='constrained')
    for ax, (comp, color, marker, label) in zip(axs, COMPS):
        for cell, fill, dx in (('250', True, -0.08), ('500', False, 0.08)):
            vals = [df[(df['layout'] == lay) & (df['cell'] == cell)].iloc[0] for df, lay, _ in cats]
            ax.errorbar(x + dx, [v[comp] for v in vals], yerr=[v[f'{comp}_err'] for v in vals],
                        color=color, marker=marker, ms=7, mfc=color if fill else 'white', ls='none',
                        elinewidth=1, label=f'cell {cell} px')
        ax.axhline(0, color='gray', ls='--', lw=1)
        ax.set_ylabel(label)
        ax.legend(frameon=False, fontsize=9)
    axs[-1].set_xticks(x)
    axs[-1].set_xticklabels([lab for _, _, lab in cats], fontsize=9)
    fig.suptitle(NOTE + ', 9-arcsec grid (45 px)', fontsize=10)
    save(fig, 'm_vs_jitter.png')

def fft_info():
    """
    Per cell size N (check_fft_sizes.py, fft_sizes_all.csv): the cell's InterpolatedImage table P,
    the drawImage FFT sizes W of the metacal images (most frequent first), and the share of
    those draws for which P/W is not an integer. None if not measured.
    """
    path = os.path.join(script_dir, 'fft_sizes_all.csv')
    if not os.path.isfile(path):
        return None
    df = pd.read_csv(path)
    cols = [c for c in df.columns if c.startswith('fft_')]
    out = {}
    for N, g in df.groupby('N'):
        P, W = int(g['table'].max()), g[cols].values.ravel()
        vals, counts = np.unique(W, return_counts=True)
        out[N] = dict(P=P, W=[int(v) for v in vals[np.argsort(-counts)]], frac=float(np.mean(P % W != 0)))
    return pd.DataFrame(out).T

def plot_cell_size(cs):
    """9-arcsec grid in cells of different sizes, with GalSim's sizes inside metacal (check_fft_sizes.py)."""
    d = cs[cs['layout'] == 'grid9'].copy()
    d['cell_size'] = d['cell'].astype(int)
    d = d.sort_values('cell_size')
    labels, xlabel = [str(N) for N in d['cell_size']], 'cell size [pixels]'
    fft = fft_info()
    if fft is not None:
        labels = [f'{N}\n{"/".join(map(str, fft.loc[N, "W"]))}\n{fft.loc[N, "P"]}\n{fft.loc[N, "frac"]:.0%}'
                  for N in d['cell_size']]
        xlabel = ('cell size  /  drawImage FFT W (most frequent first)  /  image table P   [pixels]\n'
                  '/  share of metacal draws with P/W not an integer')
    x = np.arange(len(d))
    fig, ax = plt.subplots(figsize=(8, 4.9), layout='constrained')
    for (comp, color, marker, label), dx in zip(COMPS, (-0.1, 0.1)):
        ax.errorbar(x + dx, d[comp], yerr=d[f'{comp}_err'], color=color, marker=marker, ms=7,
                    ls='none', elinewidth=1, label=label)
    ax.axhline(0, color='gray', ls='--', lw=1)
    ax.set_xticks(x)
    ax.set_xticklabels(labels, fontsize=8)
    ax.set_xlabel(xlabel, fontsize=9)
    ax.set_ylabel(r'$m$')
    ax.legend(frameon=False, loc='lower right')
    ax.set_title(NOTE + ', 9-arcsec grid (45 px)', fontsize=10)
    save(fig, 'm_vs_cell_size.png')

def plot_map(mp, ref_cell=None, name='m_map_cell_spacing.png', box=(45, 250)):
    """
    m against cell size and grid spacing: heatmaps (top) and m against spacing per cell size
    (bottom), compared with ref_cell (default the largest cell); box marks one (spacing, cell).
    """
    d = grid_spacing(mp)
    d['spacing_px'] = (d['spacing'] / PIXEL_SCALE).round(2)
    d['cell_size'] = d['cell'].astype(int)
    ref_cell = d['cell_size'].max() if ref_cell is None else ref_cell
    ref = d[d['cell_size'] == ref_cell]
    d = d[d['cell_size'] != ref_cell]
    xs, ys = np.sort(d['spacing_px'].unique()), np.sort(d['cell_size'].unique())
    ylabels, ylabel = [str(N) for N in ys], 'cell size [pixels]'
    fft = fft_info()
    if fft is not None:
        ylabels = [f'{N} ({fft.loc[N, "frac"]:.0%})' if N in fft.index else str(N) for N in ys]
        ylabel += '  (draws with P/W not an integer)'
    vmax = np.abs(d[['m1', 'm2']].values).max()
    dx, dy = np.diff(xs).min(), np.diff(ys).min()
    xe = np.concatenate([xs - dx / 2, [xs[-1] + dx / 2]])
    ye = np.concatenate([ys - dy / 2, [ys[-1] + dy / 2]])
    cmap = plt.get_cmap('viridis', len(ys))
    fig, axs = plt.subplots(2, 2, figsize=(11, 8.6), layout='constrained', height_ratios=(1.15, 1))
    for i, (comp, _, marker, label) in enumerate(COMPS):
        ## heatmap, values in units of 1e-2
        Z = d.pivot(index='cell_size', columns='spacing_px', values=comp).loc[ys, xs].values
        ax = axs[0, i]
        pc = ax.pcolormesh(xe, ye, Z, cmap='RdBu_r', vmin=-vmax, vmax=vmax)
        for (iy, ix), z in np.ndenumerate(Z):
            ax.text(xs[ix], ys[iy], f'{100 * z:.1f}', ha='center', va='center', fontsize=7,
                    color='white' if abs(z) > 0.6 * vmax else 'k')
        title = f'{label} (numbers: units of 0.01'
        if box is not None:
            ax.add_patch(plt.Rectangle((box[0] - dx / 2, box[1] - dy / 2), dx, dy, fill=False, ec='k', lw=2))
            title += f'; box: the simulation, {box[0]:g} px in {box[1]}-px cells'
        ax.set_xticks(xs[::2])
        ax.set_yticks(ys)
        ax.set_yticklabels(ylabels, fontsize=9)
        ax.set_xlabel('grid spacing [pixels]')
        ax.set_ylabel(ylabel, fontsize=10)
        ax.set_title(title + ')', fontsize=9)
        sec = ax.secondary_xaxis('top', functions=(lambda p: p * PIXEL_SCALE, lambda x: x / PIXEL_SCALE))
        sec.set_xlabel('grid spacing [arcsec]', fontsize=9)
        ## curves
        ax = axs[1, i]
        for k, N in enumerate(ys):
            c = d[d['cell_size'] == N].sort_values('spacing_px')
            ax.errorbar(c['spacing_px'], c[comp], yerr=c[f'{comp}_err'], color=cmap(k), marker=marker, ms=4,
                        lw=2.5 if box is not None and N == box[1] else 1.2, elinewidth=1, label=f'cell {N} px')
        r = ref.sort_values('spacing_px')
        ax.errorbar(r['spacing_px'], r[comp], yerr=r[f'{comp}_err'], color='0.5', marker=marker, ms=4,
                    mfc='white', ls='--', lw=1, label=f'cell {ref_cell} px')
        ax.axhline(0, color='gray', ls=':', lw=1)
        ax.axvline(45, color='#e1e0d9', lw=6, zorder=0)
        ax.set_xlabel('grid spacing [pixels]')
        ax.set_ylabel(label)
        if i == 1:
            ax.legend(frameon=False, fontsize=8, loc='lower right', ncol=2)
    fig.colorbar(pc, ax=axs[0, :], shrink=0.9, label=r'$m$')
    fig.suptitle(NOTE, fontsize=10)
    save(fig, name)

def plot_table_margin(spacings_px=(43.5, 44., 45., 88., 90.), xmax=144):
    """
    m against P - 4N (P: the cell's InterpolatedImage table, good_fft_size(4N); N: the cell size),
    from every layout test with these spacings, coloured by P; one column per spacing.
    """
    fft = fft_info()
    frames = [grid_spacing(d) for d in (read(t) for t in ('', '_cell_scan', '_map', '_fft_margin', '_map500'))
              if d is not None]
    d = pd.concat(frames)
    d['N'] = d['cell'].astype(int)
    d['spacing_px'] = (d['spacing'] / PIXEL_SCALE).round(2)
    d = d[d['N'].isin(fft.index) & d['spacing_px'].isin(spacings_px)].drop_duplicates(['spacing_px', 'N'])
    d['P'] = fft.loc[d['N'], 'P'].values.astype(int)
    d['margin'] = d['P'] - 4 * d['N']
    colors = {1024: '#2a78d6', 1536: '#3aa655', 2048: '#eb6834'}
    fig, axs = plt.subplots(2, len(spacings_px), figsize=(3.2 * len(spacings_px), 6.4), sharex=True,
                            sharey='row', layout='constrained')
    for j, s in enumerate(spacings_px):
        for i, (comp, _, marker, label) in enumerate(COMPS):
            ax = axs[i, j]
            for P, c in d[(d['spacing_px'] == s) & (d['margin'] <= xmax)].groupby('P'):
                c = c.sort_values('margin')
                ax.errorbar(c['margin'], c[comp], yerr=c[f'{comp}_err'], color=colors.get(P, 'k'), marker=marker,
                            ms=5, lw=1, elinewidth=1, label=f'table {P} px')
                for _, r in c.iterrows():
                    ax.annotate(str(r['N']), (r['margin'], r[comp]), textcoords='offset points',
                                xytext=(3, 4 if P != 2048 else -9), fontsize=6, color=colors.get(P, 'k'))
            ax.axhline(0, color='gray', ls=':', lw=1)
            ax.axvspan(40, 48, color='#e1e0d9', zorder=0)
            if j == 0:
                ax.set_ylabel(label)
            if i == 0:
                ax.set_title(f'spacing {s:g} px ({s * PIXEL_SCALE:g}″)', fontsize=10)
            if i == 1:
                ax.set_xlabel('P − 4N [pixels]')
            if i == 0 and j == 0:
                ax.legend(frameon=False, fontsize=8, loc='lower right')
    axs[0, 0].set_xlim(-6, xmax)
    fig.suptitle(NOTE + '\nP: the cell\'s image table, good_fft_size(4N); N: cell size (labels); shaded: between '
                 f'the largest biased and the smallest clean margin; cells with P − 4N > {xmax} px are all clean',
                 fontsize=9)
    save(fig, 'm_vs_table_margin.png')

## forced-FFT variants of run_fft_size_test.sh: name, cell, forced W, pad_factor
FFT_VARIANTS = [('default', 250, None, None), ('W256', 250, 256, None), ('W320', 250, 320, None),
                ('W384', 250, 384, None), ('W448', 250, 448, None), ('W512', 250, 512, None),
                ('W640', 250, 640, None), ('W768', 250, 768, None), ('W1024', 250, 1024, None),
                ('W2048', 250, 2048, None), ('pad6', 250, None, 6.), ('pad8', 250, None, 8.),
                ('pad8_W512', 250, 512, 8.), ('pad6_W768', 250, 768, 6.),
                ('cell500_W768', 500, 768, None), ('cell500_W1024', 500, 1024, None)]

def plot_forced_fft(spacings_px=(43.5, 44., 45.)):
    """m of the forced-FFT variants, per spacing; filled: P/W an integer, open: not."""
    import galsim
    from fractions import Fraction
    fft = fft_info()
    rows = []
    for name, N, W, pad in FFT_VARIANTS:
        d = read(f'_fft_{name}')
        if d is None:
            continue
        P = galsim.Image.good_fft_size(int(np.ceil((pad or 4.) * N)))
        W_eff = W if W is not None else fft.loc[N, 'W'][0]
        d = grid_spacing(d)
        d['spacing_px'] = (d['spacing'] / PIXEL_SCALE).round(2)
        for _, r in d[d['cell'] == str(N)].iterrows():
            rows.append(dict(name=name, N=N, P=P, W=W_eff, W_forced=W is not None, ratio=Fraction(P, W_eff),
                             spacing_px=r['spacing_px'], **{k: r[k] for k in ('m1', 'm1_err', 'm2', 'm2_err')}))
    ## the default 500-px cells, from the FFT-margin test
    d = read('_fft_margin')
    if d is not None:
        d = grid_spacing(d)
        d['spacing_px'] = (d['spacing'] / PIXEL_SCALE).round(2)
        for _, r in d[d['cell'] == '500'].iterrows():
            rows.append(dict(name='cell500', N=500, P=2048, W=512, W_forced=False, ratio=Fraction(2048, 512),
                             spacing_px=r['spacing_px'], **{k: r[k] for k in ('m1', 'm1_err', 'm2', 'm2_err')}))
    d = pd.DataFrame(rows)
    names = list(dict.fromkeys(d['name']))
    info = d.drop_duplicates('name').set_index('name')
    labels = [f'{info.loc[n, "W"]}{"" if info.loc[n, "W_forced"] else "*"}\n{info.loc[n, "P"]}\n{info.loc[n, "ratio"]}'
              for n in names]
    groups = [(names.index('default'), names.index('W2048'), 'cell 250 px: drawImage FFT W changed'),
              (names.index('pad6'), names.index('pad6_W768'), 'cell 250 px: table P changed'),
              (names.index('cell500_W768'), len(names) - 1, 'cell 500 px')]
    spacing_colors = dict(zip(spacings_px, ('#1b9e77', '#d95f02', '#7570b3')))
    fig, axs = plt.subplots(2, 1, figsize=(12, 7.5), sharex=True, layout='constrained')
    for ax, (comp, _, marker, label) in zip(axs, COMPS):
        for k, sp in enumerate(spacings_px):
            for _, r in d[d['spacing_px'] == sp].iterrows():
                integer = r['ratio'].denominator == 1
                c = spacing_colors[sp]
                ax.errorbar(names.index(r['name']) + 0.18 * (k - 1), r[comp], yerr=r[f'{comp}_err'], color=c,
                            marker=marker, ms=6, mfc=c if integer else 'white', ls='none', elinewidth=1)
        ax.axhline(0, color='gray', ls=':', lw=1)
        ax.set_yscale('symlog', linthresh=0.003, linscale=0.6)
        ax.axhspan(-0.001, 0.001, color='#e1e0d9', zorder=0)
        ax.set_ylabel(label)
        for i0, i1, text in groups:
            if i0 > 0:
                ax.axvline(i0 - 0.5, color='0.6', lw=1)
            if ax is axs[0]:
                ax.text((i0 + i1) / 2, 1.02, text, transform=ax.get_xaxis_transform(), ha='center', fontsize=9)
    axs[-1].set_xticks(range(len(names)))
    axs[-1].set_xticklabels(labels, fontsize=8)
    axs[-1].set_xlabel('drawImage FFT W (* GalSim\'s choice for most draws)  /  image table P  /  P/W   [pixels]',
                       fontsize=9)
    handles = [plt.Line2D([], [], color=spacing_colors[sp], marker='o', ls='none', label=f'spacing {sp:g} px')
               for sp in spacings_px] + \
              [plt.Line2D([], [], color='k', marker='o', ls='none', label='P/W an integer'),
               plt.Line2D([], [], color='k', marker='o', mfc='white', ls='none', label='P/W not an integer')]
    fig.legend(handles=handles, loc='outside lower center', ncol=5, frameon=False, fontsize=9)
    fig.suptitle(NOTE + '; GalSim sizes inside metacal changed, scenes unchanged\n'
                 'y axis linear within ±0.003, logarithmic beyond; shaded: |m| < 0.001', fontsize=9)
    save(fig, 'm_vs_forced_fft.png')

if __name__ == '__main__':
    coarse, fine = read('_spacing_scan'), read('_fine_scan')
    coarse = grid_spacing(coarse) if coarse is not None else None
    fine = grid_spacing(fine) if fine is not None else None
    if coarse is not None:
        plot_spacing(coarse, fine)
    if fine is not None:
        plot_spacing(coarse, fine, xlim=(fine['spacing'].min(), fine['spacing'].max()),
                     name='m_vs_grid_spacing_fine.png')
    jit = read('_jitter')
    if jit is not None:
        plot_jitter(jit, read(''))
    cs = read('_cell_scan')
    if cs is not None:
        plot_cell_size(cs)
    mp = read('_map')
    if mp is not None:
        plot_map(mp)
    mp500 = read('_map500')
    if mp500 is not None:
        plot_map(mp500, ref_cell=498, name='m_map_cell_spacing_500.png', box=None)
    if os.path.isfile(os.path.join(script_dir, 'layout_response_fft_margin.csv')) and fft_info() is not None:
        plot_table_margin()
    if read('_fft_default') is not None:
        plot_forced_fft()
