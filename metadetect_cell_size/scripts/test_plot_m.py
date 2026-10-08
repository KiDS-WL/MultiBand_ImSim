# -*- coding: utf-8 -*-

### Residual multiplicative shear bias m for each metadetect cell set-up
###    reads the summary saved by analyse_cell_size.py

import os
import re

import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

plt.rcParams["text.usetex"] = True

# +++ general settings for plot
mpl.rcParams['xtick.direction'] = 'in'
mpl.rcParams['ytick.direction'] = 'in'
mpl.rcParams['xtick.top'] = True
mpl.rcParams['ytick.right'] = True
plt.rc('font', size=14, family='serif')

# +++ I/O
script_dir = os.path.dirname(os.path.abspath(__file__))
inpath = os.path.join(script_dir, 'results_cell_size.csv')
outpath = os.path.join(script_dir, 'm_cell_size.png')
# outpath = 'show'
outpath_buffer = os.path.join(script_dir, 'm_buffer_size.png')
# outpath_buffer = 'show'
outpath_nobj = os.path.join(script_dir, 'm_N_per_central.png')
# outpath_nobj = 'show'

## the collaborator's run, which the others are compared with
baseline = 'cell250_central150'

# +++ plot related
COLORs = ['#2a78d6', '#eb6834']
SYMBOLs = ['o', 's']
LABELs = [r'$m_1$', r'$m_2$']
MS = 7
ELW = 1
## m1 and m2 side by side within each set-up, each corrected by R = (R11+R22)/2
##    (filled, as in B) and by its own response R11 or R22 (open)
COMPs = [['m1', 'm1_R11'], ['m2', 'm2_R22']]
LABELs_R = [[r'$m_1$, $R = (R_{11}+R_{22})/2$', r'$m_1$, $R_{11}$'],
            [r'$m_2$, $R = (R_{11}+R_{22})/2$', r'$m_2$, $R_{22}$']]
X_OFFSETs = [[-0.27, -0.09], [0.09, 0.27]]
## extra space between the two cell sizes
GROUP_GAP = 0.6
XLABEL = r'central size [pixels]'
YLABEL = r'$m$'
## m against a property of the set-up: one style per cell size
##    the colour stays that of the component, as above
CELL_SYMBOLs = {250: 'o', 260: 'D', 500: 's'}
CELL_LSs = {250: '-', 260: '--', 500: ':'}
CELL_FILLs = {250: True, 260: True, 500: False}
## the same buffer size occurs in several cell sizes, so shift them apart slightly
CELL_DXs = {250: -3, 260: 0, 500: 3}
XLABEL_BUFFER = r'buffer size (cell $-$ central) [pixels]'
XLABEL_NOBJ = r'detections per central region'
## detections per pixel^2: every set-up measures the same images, so this is the
##    same for all of them, 10.74 per 150x150 central region (baseline,
##    zero-shear catalogues, where the grid fills the whole image)
DENSITY = 10.74 / 150**2

def save_or_show(outpath):
    plt.tight_layout()
    if outpath == 'show':
        plt.show()
    else:
        plt.savefig(outpath, dpi=300)
        print('plot saved in', outpath)
    plt.close()

# +++ data
data = pd.read_csv(inpath)
cell_central = data['setup'].str.extract(r'cell(\d+)_central(\d+)').astype(int)
data['cell_size'] = cell_central[0]
data['central_size'] = cell_central[1]
data['buffer_size'] = data['cell_size'] - data['central_size']
data['N_per_central'] = DENSITY * data['central_size']**2
data = data.sort_values(['cell_size', 'central_size']).reset_index(drop=True)

## x position: one slot per set-up, the cell sizes in separate groups
cell_sizes = sorted(data['cell_size'].unique())
i_group = data['cell_size'].map({c: i for i, c in enumerate(cell_sizes)}).values
x_cen = np.arange(len(data)) + GROUP_GAP * i_group

# +++ plot
fig, ax = plt.subplots(figsize=(8, 5))

## mark the baseline column
for x0, setup in zip(x_cen, data['setup']):
    if setup == baseline:
        ax.axvspan(x0 - 0.4, x0 + 0.4, color='#f0efec', zorder=0)

for i, comps in enumerate(COMPs):
    for j, comp in enumerate(comps):
        ax.errorbar(x_cen + X_OFFSETs[i][j], data[comp].values, yerr=data[f'{comp}_err'].values,
                    color=COLORs[i], marker=SYMBOLs[i], markersize=MS,
                    mfc=COLORs[i] if j == 0 else 'white', elinewidth=ELW,
                    ls='none', label=LABELs_R[i][j])

ax.axhline(y=0, color='gray', ls='--', lw=1)

## separate the cell sizes and name each group
for i, cell_size in enumerate(cell_sizes):
    mask = (i_group == i)
    if i > 0:
        ax.axvline(x=(x_cen[mask][0] + x_cen[~mask & (i_group == i-1)][-1]) / 2.,
                   color='black', ls='-', lw=1)
    ax.text(x_cen[mask].mean(), 1.02, f'cell size = {cell_size}',
            transform=ax.get_xaxis_transform(), ha='center', va='bottom')

## tick labels: the central size, with the baseline and incomplete runs flagged
##    an incomplete run gets an asterisk, explained under the axis
xticklabels = []
notes_incomplete = []
N_files_full = data.loc[data['complete'], 'N_files'].max()
for _, row in data.iterrows():
    label = f"{row['central_size']}"
    if not row['complete']:
        label += r'$^*$'
        notes_incomplete.append(f"cell {row['cell_size']} / central {row['central_size']} "
                                f"({row['N_files']} of {N_files_full} catalogues)")
    if row['setup'] == baseline:
        label += '\n(baseline)'
    xticklabels.append(label)
ax.set_xticks(x_cen)
ax.set_xticklabels(xticklabels)
ax.tick_params(axis='x', length=0)
ax.set_xlim(x_cen[0] - 0.6, x_cen[-1] + 0.6)

if notes_incomplete:
    ax.set_xlabel(XLABEL + '\n' + r'{\footnotesize $^*$incomplete run: ' + ', '.join(notes_incomplete) + '}')
else:
    ax.set_xlabel(XLABEL)
ax.set_ylabel(YLABEL)
ax.legend(frameon=False, loc='lower right', handletextpad=0.2, fontsize=12)

save_or_show(outpath)

# +++ plot: m against a property of the set-up, the cell sizes told apart
def plot_m_vs(xcol, xlabel, outpath, dxs=None, xlog=False):
    fig, axs = plt.subplots(1, 2, figsize=(10, 4.5), sharex=True, sharey=True)

    for i_comp, (ax, comp) in enumerate(zip(axs, ['m1', 'm2'])):
        for cell_size in cell_sizes:
            data_tmp = data[data['cell_size'] == cell_size].sort_values(xcol)
            x_val = data_tmp[xcol].values + (dxs[cell_size] if dxs else 0)
            ax.errorbar(x_val, data_tmp[comp].values, yerr=data_tmp[f'{comp}_err'].values,
                        color=COLORs[i_comp], marker=CELL_SYMBOLs[cell_size], markersize=MS,
                        mfc=COLORs[i_comp] if CELL_FILLs[cell_size] else 'white',
                        elinewidth=ELW, ls=CELL_LSs[cell_size], lw=1)
            ## ring the baseline
            mask_base = (data_tmp['setup'] == baseline).values
            ax.plot(x_val[mask_base], data_tmp[comp].values[mask_base], ls='none',
                    marker='o', markersize=2.4*MS, mfc='none', mec='gray', mew=1)

        ax.axhline(y=0, color='gray', ls='--', lw=1)
        ax.set_title(LABELs[i_comp])
        ax.set_xlabel(xlabel)

    x_all = data[xcol].values
    if xlog:
        axs[0].set_xscale('log')
        axs[0].set_xlim(x_all.min() / 1.5, x_all.max() * 1.5)
        axs[0].xaxis.set_major_formatter(mpl.ticker.FormatStrFormatter('%g'))
    else:
        x_unique = np.unique(x_all)
        axs[0].set_xticks(x_unique)
        x_pad = 0.5 * np.min(np.diff(x_unique))
        axs[0].set_xlim(x_unique[0] - x_pad, x_unique[-1] + x_pad)
    axs[0].set_ylabel(YLABEL)

    ## the cell sizes are told apart by the marker and line style, so the legend is in black
    handles = [Line2D([], [], color='black', marker=CELL_SYMBOLs[c], ls=CELL_LSs[c], lw=1,
                      markersize=MS, mfc='black' if CELL_FILLs[c] else 'white')
               for c in cell_sizes]
    handles.append(Line2D([], [], ls='none', marker='o', markersize=2.4*MS,
                          mfc='none', mec='gray', mew=1))
    axs[0].legend(handles, [f'cell size = {c}' for c in cell_sizes] + ['baseline'],
                  frameon=False, loc='lower left', handletextpad=0.4)

    save_or_show(outpath)

plot_m_vs('buffer_size', XLABEL_BUFFER, outpath_buffer, dxs=CELL_DXs)
plot_m_vs('N_per_central', XLABEL_NOBJ, outpath_nobj, xlog=True)

## the numbers behind the plots
print(data[['setup', 'buffer_size', 'N_per_central', 'complete',
            'm1', 'm1_err', 'm2', 'm2_err']].to_string(index=False, float_format='%.4f'))
