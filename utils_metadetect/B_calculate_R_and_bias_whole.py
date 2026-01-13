# -*- coding: utf-8 -*-
# @Author: lshuns
# @Date:   1969-12-31 16:00:00
# @Last Modified by:   lshuns
# @Last Modified time: 2026-01-13 13:37:00

### Calculate the mean response and residual biases for the whole sample 

import os
import re
import glob

import numpy as np 
import pandas as pd
from scipy.optimize import curve_fit

## ++++++++++++++ I/O and general setups

## Where to find the simulations
main_dir = '/sdf/data/kipac/u/liss/ImSim/output/test_dev/LSST_r/'

## Shear inputs in simulations
shear_tags = ['m283m283', 'm283p283', 'p283m283', 'p283p283']
g1_input_list = [-0.0283, -0.0283, 0.0283, 0.0283]
g2_input_list = [-0.0283, 0.0283, -0.0283, 0.0283]

## Shear types used in metadetect
shear_types = ['noshear', '1p', '1m', '2p', '2m']

## What is the fitting model used in metadetect
fit_model = 'wmom'

## Which shear weight to use
#### None = No weighting
# which_weight = None
#### weight_sigma_e = weight based on sigma_e
# which_weight = 'weight_sigma_e'
#### weight_e = weight based on e
# which_weight = 'weight_e'
#### shear_weight = weight_sigma_e * weight_e
which_weight = 'shear_weight'

## ++++++++++++++ Workhorse

## Loop over simulations to get all catalogues
cata = []
for i_shear, shear_tag in enumerate(shear_tags):
    inpath_list = glob.glob(os.path.join(main_dir, 
                                         shear_tag, 
                                         'catalogues/shapes_metadetect', 
                                         '*.feather'))
    print(f">>> Number of catalogues found in {shear_tag}: {len(inpath_list)}")
    for inpath in inpath_list:
        cata_tmp = pd.read_feather(inpath)
        ## Drop zero-weight objects and useless columns for memory
        if which_weight is not None:
            cata_tmp = cata_tmp.loc[cata_tmp[which_weight]>0, 
                            ['shear_type', 
                            f'{fit_model}_s2n', 
                            f'{fit_model}_g_1',
                            f'{fit_model}_g_2', 
                            f'{fit_model}_T', 
                            f'{fit_model}_T_ratio',
                            f'{fit_model}_band_flux', 
                            which_weight]]
            ## Renaming for easy use
            cata_tmp = cata_tmp.rename(columns={which_weight: 
                                                'weight'})
        else:
            cata_tmp = cata_tmp[ 
                            ['shear_type', 
                            f'{fit_model}_s2n', 
                            f'{fit_model}_g_1',
                            f'{fit_model}_g_2', 
                            f'{fit_model}_T', 
                            f'{fit_model}_T_ratio',
                            f'{fit_model}_band_flux']]
            ## No weighting
            cata_tmp['weight'] = 1

        ## Add input info
        cata_tmp['g1_input'] = g1_input_list[i_shear]
        cata_tmp['g2_input'] = g2_input_list[i_shear]
        cata_tmp['run_tag'] = shear_tag        
        cata_tmp['tile_label'] = re.search(r'tile(.*)_rot', os.path.basename(inpath)).group(1)
        ## rotation info
        cata_tmp['gal_rot'] = float(re.search(r'_rot(\d+)', os.path.basename(inpath)).group(1))

        cata.append(cata_tmp)
        del cata_tmp
cata = pd.concat(cata, ignore_index=True)

## Calculat shear response and residual shear bias for the whole sample
g1_input_all = []
g2_input_all = []
g1_measured_all = []
g2_measured_all = []
for i_shear, shear_tag in enumerate(shear_tags):
    ## The input shear
    g1_input = g1_input_list[i_shear]
    g2_input = g2_input_list[i_shear]

    ## Select simulations
    cata_tmp = cata[(cata['g1_input']==g1_input)
                    &(cata['g2_input']==g2_input)
                    ].reset_index(drop=True)
    print(f">>> Number of objects for {shear_tag}: {len(cata_tmp)}")

    ## Calculate Response
    g1_1p = np.average(cata_tmp.loc[
        cata_tmp['shear_type']=='1p', 
        f'{fit_model}_g_1'],
        weights = cata_tmp.loc[
        cata_tmp['shear_type']=='1p', 
        'weight'])
    g1_1m = np.average(cata_tmp.loc[
        cata_tmp['shear_type']=='1m', 
        f'{fit_model}_g_1'],
        weights = cata_tmp.loc[
        cata_tmp['shear_type']=='1m', 
        'weight'])
    R11 = (g1_1p - g1_1m) / 0.02
    g2_2p = np.average(cata_tmp.loc[
        cata_tmp['shear_type']=='2p', 
        f'{fit_model}_g_2'],
        weights = cata_tmp.loc[
        cata_tmp['shear_type']=='2p', 
        'weight'])
    g2_2m = np.average(cata_tmp.loc[
        cata_tmp['shear_type']=='2m', 
        f'{fit_model}_g_2'],
        weights = cata_tmp.loc[
        cata_tmp['shear_type']=='2m', 
        'weight'])
    R22 = (g2_2p - g2_2m) / 0.02
    R = (R11+R22)/2
    print('>>> R11, R22, R', R11, R22, R)

    ## Calculate measured shear for each tile
    cata_tmp = cata_tmp.loc[cata_tmp['shear_type']=='noshear', 
                        ['tile_label', 
                        f'{fit_model}_g_1',
                        f'{fit_model}_g_2', 
                        'weight']
                        ].copy().reset_index(drop=True)
    cata_tmp[f'{fit_model}_g_1'] = cata_tmp[f'{fit_model}_g_1'] * cata_tmp['weight']
    cata_tmp[f'{fit_model}_g_2'] = cata_tmp[f'{fit_model}_g_2'] * cata_tmp['weight']
    cata_tmp = cata_tmp.groupby(by=['tile_label']).sum()
    g1_out_arr = cata_tmp[f'{fit_model}_g_1'].values / cata_tmp['weight'].values
    g2_out_arr = cata_tmp[f'{fit_model}_g_2'].values / cata_tmp['weight'].values
    del cata_tmp

    ## Apply the shear response correction
    g1_out_arr = g1_out_arr/R
    g2_out_arr = g2_out_arr/R

    ## Save results
    g1_input_all.append(np.ones_like(g1_out_arr) * g1_input)
    g2_input_all.append(np.ones_like(g2_out_arr) * g2_input)
    g1_measured_all.append(g1_out_arr)
    g2_measured_all.append(g2_out_arr)
g1_input_all = np.concatenate(g1_input_all)
g2_input_all = np.concatenate(g2_input_all)
g1_measured_all = np.concatenate(g1_measured_all)
g2_measured_all = np.concatenate(g2_measured_all)
print(">>>> Total number of points for fitting", 
      len(g1_input_all), len(g2_input_all),
      len(g1_measured_all), len(g2_measured_all)
      )

## Fit lines
def line_func(x, m, c):
    return (1+m) * x + c

popt, pcov = curve_fit(line_func, 
                       g1_input_all, 
                       g1_measured_all)
m1, c1 = popt
m1_err, c1_err = np.sqrt(np.diag(pcov))

popt, pcov = curve_fit(line_func, 
                       g2_input_all, 
                       g2_measured_all)
m2, c2 = popt
m2_err, c2_err = np.sqrt(np.diag(pcov))

print(f"m1 = {m1:.4f} pm {m1_err:.4f}, c1 = {c1:.4f} pm {c1_err:.4f}")
print(f"m2 = {m2:.4f} pm {m2_err:.4f}, c2 = {c2:.4f} pm {c2_err:.4f}")

# #### No weighting
# >>> Number of catalogues found in m283m283: 100
# >>> Number of catalogues found in m283p283: 100
# >>> Number of catalogues found in p283m283: 100
# >>> Number of catalogues found in p283p283: 100
# >>> Number of objects for m283m283: 13676038
# >>> R11, R22, R 0.2429048780218206 0.24541933880659178 0.24416210841420619
# >>> Number of objects for m283p283: 13677949
# >>> R11, R22, R 0.24265313323144824 0.245113996298797 0.24388356476512263
# >>> Number of objects for p283m283: 13675674
# >>> R11, R22, R 0.24253989696479214 0.24491446765753516 0.24372718231116364
# >>> Number of objects for p283p283: 13676790
# >>> R11, R22, R 0.24288922378193 0.24544376039636748 0.24416649208914876
# >>>> Total number of points for fitting 200 200 200 200
# m1 = -0.0015 pm 0.0030, c1 = -0.0003 pm 0.0001
# m2 = -0.0014 pm 0.0031, c2 = 0.0003 pm 0.0001

# ########## with shear weight
# >>> Number of catalogues found in m283m283: 100
# >>> Number of catalogues found in m283p283: 100
# >>> Number of catalogues found in p283m283: 100
# >>> Number of catalogues found in p283p283: 100
# >>> Number of objects for m283m283: 11931865
# >>> R11, R22, R 0.24700703793440265 0.24793567562786542 0.24747135678113402
# >>> Number of objects for m283p283: 11933233
# >>> R11, R22, R 0.2467654896390211 0.2482768075850065 0.2475211486120138
# >>> Number of objects for p283m283: 11931149
# >>> R11, R22, R 0.24702901602454302 0.24788770956919223 0.24745836279686761
# >>> Number of objects for p283p283: 11932616
# >>> R11, R22, R 0.24737965443167542 0.24834342449650937 0.2478615394640924
# >>>> Total number of points for fitting 200 200 200 200
# m1 = -0.0037 pm 0.0028, c1 = -0.0001 pm 0.0001
# m2 = -0.0043 pm 0.0028, c2 = 0.0002 pm 0.0001

# # ########## with weight_sigma_e
# >>> Number of catalogues found in m283m283: 100
# >>> Number of catalogues found in m283p283: 100
# >>> Number of catalogues found in p283m283: 100
# >>> Number of catalogues found in p283p283: 100
# >>> Number of objects for m283m283: 11931865
# >>> R11, R22, R 0.2757788850685054 0.27725671801545076 0.2765178015419781
# >>> Number of objects for m283p283: 11933233
# >>> R11, R22, R 0.27550996927270255 0.2773742633973511 0.2764421163350268
# >>> Number of objects for p283m283: 11931149
# >>> R11, R22, R 0.2756998793528352 0.27715186021190613 0.27642586978237066
# >>> Number of objects for p283p283: 11932616
# >>> R11, R22, R 0.2760936638963224 0.27739741137293555 0.27674553763462895
# >>>> Total number of points for fitting 200 200 200 200
# m1 = -0.0034 pm 0.0032, c1 = -0.0001 pm 0.0001
# m2 = -0.0042 pm 0.0030, c2 = 0.0003 pm 0.0001

# # # ########## with weight_e
# >>> Number of catalogues found in m283m283: 100
# >>> Number of catalogues found in m283p283: 100
# >>> Number of catalogues found in p283m283: 100
# >>> Number of catalogues found in p283p283: 100
# >>> Number of objects for m283m283: 11931865
# >>> R11, R22, R 0.2336772137522222 0.23542313011278643 0.23455017193250433
# >>> Number of objects for m283p283: 11933233
# >>> R11, R22, R 0.2329592587125685 0.23579887987999512 0.23437906929628183
# >>> Number of objects for p283m283: 11931149
# >>> R11, R22, R 0.23335889248125774 0.2354236006388438 0.23439124656005078
# >>> Number of objects for p283p283: 11932616
# >>> R11, R22, R 0.23403879152601526 0.23599435839800603 0.23501657496201064
# >>>> Total number of points for fitting 200 200 200 200
# m1 = -0.0039 pm 0.0029, c1 = -0.0001 pm 0.0001
# m2 = -0.0043 pm 0.0032, c2 = 0.0004 pm 0.0001