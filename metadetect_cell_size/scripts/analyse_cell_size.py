# -*- coding: utf-8 -*-

### Shear response and residual shear bias for each metadetect cell set-up
###
###    Same weights, cuts and fit as utils_metadetect/A_assign_weights.py and
###    utils_metadetect/B_calculate_R_and_bias_whole.py, but the weights are
###    computed in memory: the baseline (250/150) catalogues are the
###    collaborator's and read-only, so they cannot be rewritten in place the way
###    A_assign_weights.py does.
###
###    The weighted averages are accumulated as sums, file by file. That gives the
###    same numbers as B without holding every catalogue in memory at once.
###
###    m1, m2 use R = (R11+R22)/2 for both components, as B does. m1_R11 and
###    m2_R22 are the same, except that each component is corrected by its own
###    response: g1 by R11, g2 by R22.

import os
import re
import glob
from multiprocessing import Pool

import numpy as np
import pandas as pd
from scipy.optimize import curve_fit

## ++++++++++++++ I/O and general setups

## The baseline run, and where the new set-ups are written
baseline_dir = '/sdf/data/kipac/u/jlkitt/imsims/grid/seeing_073/'
main_dir = '/sdf/data/kipac/u/liss/ImSim/output/test_dev/jamar_seeing_073/'

## label: (main directory, shape folder)
setups = {
    'cell250_central150': (baseline_dir, 'shapes'),
    'cell250_central50': (main_dir, 'shapes_cell250_central50'),
    'cell250_central100': (main_dir, 'shapes_cell250_central100'),
    'cell250_central200': (main_dir, 'shapes_cell250_central200'),
    'cell500_central250': (main_dir, 'shapes_cell500_central250'),
    'cell500_central300': (main_dir, 'shapes_cell500_central300'),
    'cell500_central400': (main_dir, 'shapes_cell500_central400'),
    ## same central sizes as cell 250, to separate the cell size from the central size
    'cell500_central50': (main_dir, 'shapes_cell500_central50'),
    'cell500_central100': (main_dir, 'shapes_cell500_central100'),
    'cell500_central200': (main_dir, 'shapes_cell500_central200'),
    ## a cell size outside the window where metacal's GalSim FFT sizes are incommensurate
    ##    (P/W = 4 for every draw), with the baseline's central size: the same galaxies as 250/150
    'cell260_central150': (main_dir, 'shapes_cell260_central150'),
}

## Shear inputs in simulations
shear_inputs = {'m020m020': (-0.02, -0.02),
                'm020p020': (-0.02, 0.02),
                'p000p000': (0.0, 0.0),
                'p020m020': (0.02, -0.02),
                'p020p020': (0.02, 0.02)}

## Number of catalogues per shear tag: 100 tiles x 2 rotations
N_files_expected = 200

## What is the fitting model used in metadetect
fit_model = 'wmom'

## Cut info (as A_assign_weights.py)
snr_min = 12.5
resolution_min = 1.2

## The intrinsic ellipticity dispersion (as A_assign_weights.py)
sigma_SN = 0.07

## Where to save the summary
outpath_csv = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                           'results_cell_size.csv')

## ++++++++++++++ Workhorse

## Fit lines
def line_func(x, m, c):
    return (1+m) * x + c

## catalogue columns needed for the weights and the shear
weight_columns = ['shear_type',
                  f'{fit_model}_g_1',
                  f'{fit_model}_g_2',
                  f'{fit_model}_g_cov_as_sigma',
                  f'{fit_model}_s2n',
                  f'{fit_model}_T_ratio',
                  f'{fit_model}_flags']

def shear_weight(cata):
    """
    The shear weight of A_assign_weights.py (shear_weight), zero for objects
    failing the cuts or with a nan measurement.
    """
    g1 = cata[f'{fit_model}_g_1'].values
    g2 = cata[f'{fit_model}_g_2'].values

    ## Weight based on ellipticity noise
    weight_sigma_e = 1./(sigma_SN**2 + np.square(cata[f'{fit_model}_g_cov_as_sigma'].values))
    ## Weight based on ellipticity
    e_sq = np.square(g1) + np.square(g2)
    weight_e = np.square(1-e_sq) * np.exp(-1*e_sq/2./0.09)
    ## Combine them together
    weight = weight_sigma_e * weight_e

    ## S/N and resolution cut
    mask_cut = ((cata[f'{fit_model}_s2n'].values <= snr_min)
                | (cata[f'{fit_model}_T_ratio'].values <= resolution_min)
                | (cata[f'{fit_model}_flags'].values != 0))
    weight[mask_cut] = 0.
    ## Zero weight for nan measurements
    weight = np.nan_to_num(weight, nan=0.)
    return weight

def file_sums(inpath):
    """
    Weighted sums of one catalogue, for the response and the per-tile shear.

    Only objects with a non-zero weight are kept, as
    B_calculate_R_and_bias_whole.py does.
    """
    cata = pd.read_feather(inpath, columns=weight_columns)
    g1 = cata[f'{fit_model}_g_1'].values
    g2 = cata[f'{fit_model}_g_2'].values
    weight = shear_weight(cata)

    mask = weight > 0
    sums = {}
    for shear_type in ('noshear', '1p', '1m', '2p', '2m'):
        mask_tmp = mask & (cata['shear_type'].values == shear_type)
        w = weight[mask_tmp]
        sums[shear_type] = np.array([w.sum(),
                                     (w * g1[mask_tmp]).sum(),
                                     (w * g2[mask_tmp]).sum(),
                                     mask_tmp.sum()])
    return inpath, sums

def bias_one_setup(label, main_dir_tmp, shape_folder, pool):
    print(f"\n>>>>>>>>>>>>>> {label} ({os.path.join(main_dir_tmp, '*', 'catalogues', shape_folder)})")

    ## Find all the catalogues
    inpath_dict = {}
    for shear_tag in shear_inputs:
        inpath_dict[shear_tag] = sorted(glob.glob(os.path.join(main_dir_tmp,
                                                               shear_tag,
                                                               'catalogues',
                                                               shape_folder,
                                                               '*.feather')))
    N_files = {shear_tag: len(v) for shear_tag, v in inpath_dict.items()}
    print(f">>> Number of catalogues found: {N_files}")
    if min(N_files.values()) == 0:
        print(">>> Some shear tags have no catalogue, skipped!")
        return None
    complete = all(N == N_files_expected for N in N_files.values())
    if not complete:
        print(f">>> WARNING: incomplete, expected {N_files_expected} catalogues per shear tag!")

    ## Weighted sums for every catalogue
    all_inpaths = [p for v in inpath_dict.values() for p in v]
    sums_dict = dict(pool.imap_unordered(file_sums, all_inpaths, chunksize=4))

    return bias_from_sums(label, inpath_dict, sums_dict)

def bias_from_sums(label, inpath_dict, sums_dict):
    """
    R and m of one set-up from the weighted sums of its catalogues.

    inpath_dict: {shear tag: [catalogue paths]}
    sums_dict: {catalogue path: sums, as returned by file_sums}
    """
    N_files = {shear_tag: len(v) for shear_tag, v in inpath_dict.items()}
    complete = all(N == N_files_expected for N in N_files.values())

    g1_input_all = []
    g2_input_all = []
    g1_measured_all = []
    g2_measured_all = []
    g1_measured_R11_all = []
    g2_measured_R22_all = []
    R_list = []
    R11_list = []
    R22_list = []
    N_obj = 0
    for shear_tag, (g1_input, g2_input) in shear_inputs.items():
        sums_list = [sums_dict[p] for p in inpath_dict[shear_tag]]

        ## Calculate Response
        tot = {k: np.sum([s[k] for s in sums_list], axis=0)
               for k in ('noshear', '1p', '1m', '2p', '2m')}
        R11 = (tot['1p'][1]/tot['1p'][0] - tot['1m'][1]/tot['1m'][0]) / 0.02
        R22 = (tot['2p'][2]/tot['2p'][0] - tot['2m'][2]/tot['2m'][0]) / 0.02
        R = (R11+R22)/2
        R_list.append(R)
        R11_list.append(R11)
        R22_list.append(R22)
        N_obj += int(sum(tot[k][3] for k in tot))
        print(f'>>> {shear_tag}: R11, R22, R', R11, R22, R)

        ## Calculate measured shear for each tile, both rotations together
        tile_sums = {}
        for inpath, sums in zip(inpath_dict[shear_tag], sums_list):
            tile_label = re.search(r'tile(.*)_rot', os.path.basename(inpath)).group(1)
            tile_sums[tile_label] = tile_sums.get(tile_label, 0.) + sums['noshear'][:3]
        tile_sums = np.array(list(tile_sums.values()))

        ## Apply the shear response correction
        g1_out_arr = tile_sums[:, 1] / tile_sums[:, 0] / R
        g2_out_arr = tile_sums[:, 2] / tile_sums[:, 0] / R
        ## and with each component's own response
        g1_out_R11_arr = tile_sums[:, 1] / tile_sums[:, 0] / R11
        g2_out_R22_arr = tile_sums[:, 2] / tile_sums[:, 0] / R22

        ## Save results
        g1_input_all.append(np.ones_like(g1_out_arr) * g1_input)
        g2_input_all.append(np.ones_like(g2_out_arr) * g2_input)
        g1_measured_all.append(g1_out_arr)
        g2_measured_all.append(g2_out_arr)
        g1_measured_R11_all.append(g1_out_R11_arr)
        g2_measured_R22_all.append(g2_out_R22_arr)
    g1_input_all = np.concatenate(g1_input_all)
    g2_input_all = np.concatenate(g2_input_all)
    g1_measured_all = np.concatenate(g1_measured_all)
    g2_measured_all = np.concatenate(g2_measured_all)
    g1_measured_R11_all = np.concatenate(g1_measured_R11_all)
    g2_measured_R22_all = np.concatenate(g2_measured_R22_all)
    print(">>>> Total number of points for fitting", len(g1_measured_all), len(g2_measured_all))

    popt, pcov = curve_fit(line_func, g1_input_all, g1_measured_all)
    m1, c1 = popt
    m1_err, c1_err = np.sqrt(np.diag(pcov))

    popt, pcov = curve_fit(line_func, g2_input_all, g2_measured_all)
    m2, c2 = popt
    m2_err, c2_err = np.sqrt(np.diag(pcov))

    print(f"m1 = {m1:.4f} pm {m1_err:.4f}, c1 = {c1:.5f} pm {c1_err:.5f}")
    print(f"m2 = {m2:.4f} pm {m2_err:.4f}, c2 = {c2:.5f} pm {c2_err:.5f}")

    ## the same, with each component corrected by its own response
    popt, pcov = curve_fit(line_func, g1_input_all, g1_measured_R11_all)
    m1_R11 = popt[0]
    m1_R11_err = np.sqrt(pcov[0, 0])

    popt, pcov = curve_fit(line_func, g2_input_all, g2_measured_R22_all)
    m2_R22 = popt[0]
    m2_R22_err = np.sqrt(pcov[0, 0])

    print(f"with R11, R22 separately: m1 = {m1_R11:.4f} pm {m1_R11_err:.4f}, "
          f"m2 = {m2_R22:.4f} pm {m2_R22_err:.4f}")

    return {'setup': label, 'complete': complete, 'N_files': sum(N_files.values()),
            'N_obj_weighted': N_obj, 'R': np.mean(R_list),
            'm1': m1, 'm1_err': m1_err, 'm2': m2, 'm2_err': m2_err,
            'c1': c1, 'c1_err': c1_err, 'c2': c2, 'c2_err': c2_err,
            'R11': np.mean(R11_list), 'R22': np.mean(R22_list),
            'm1_R11': m1_R11, 'm1_R11_err': m1_R11_err,
            'm2_R22': m2_R22, 'm2_R22_err': m2_R22_err}

if __name__ == '__main__':
    N_proc = int(os.environ.get('SLURM_CPUS_PER_TASK', 4))
    results = []
    with Pool(N_proc) as pool:
        for label, (main_dir_tmp, shape_folder) in setups.items():
            res = bias_one_setup(label, main_dir_tmp, shape_folder, pool)
            if res is not None:
                results.append(res)

    results = pd.DataFrame(results)
    print('\n>>>>>>>>>>>>>> Summary')
    with pd.option_context('display.width', 200, 'display.max_columns', None):
        print(results[['setup', 'complete', 'R', 'm1', 'm1_err', 'm2', 'm2_err',
                       'c1', 'c1_err', 'c2', 'c2_err']].to_string(index=False, float_format='%.5f'))
        print(results[['setup', 'R11', 'R22', 'm1_R11', 'm1_R11_err',
                       'm2_R22', 'm2_R22_err']].to_string(index=False, float_format='%.5f'))
    results.to_csv(outpath_csv, index=False)
    print(f'saved to {outpath_csv}')
