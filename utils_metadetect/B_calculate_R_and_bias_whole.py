# -*- coding: utf-8 -*-
# @Author: lshuns
# @Date:   1969-12-31 16:00:00
# @Last Modified by:   lshuns
# @Last Modified time: 2026-09-04 13:35:39

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

## Shape folder name
shape_folder_list = ['shapes_metadetect', 'shapes_metadetect_shifted']

## Shear inputs in simulations
shear_tags = ['m283m283', 'm283p283', 'p283m283', 'p283p283']
g1_input_list = [-0.0283, -0.0283, 0.0283, 0.0283]
g2_input_list = [-0.0283, 0.0283, -0.0283, 0.0283]

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

## Fit lines
def line_func(x, m, c):
    return (1+m) * x + c

## ++++++++++++++ Run the whole analysis once per shape folder
for i_shape, shape_folder in enumerate(shape_folder_list):
    print(f"\n>>>>>>>>>>>>>> {shape_folder} "
          f"({i_shape+1}/{len(shape_folder_list)}), weight = {which_weight}")
    ## Loop over simulations to get all catalogues
    cata = []
    for i_shear, shear_tag in enumerate(shear_tags):
        inpath_list = glob.glob(os.path.join(main_dir, 
                                             shear_tag, 
                                             f'catalogues/{shape_folder}', 
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
                ## Drop the objects metadetect could not measure
                ##    the weighted branch gets this for free, because
                ##    A_assign_weights.py gives a nan measurement zero weight and
                ##    the >0 cut above then removes it. Here there is no such cut,
                ##    and a single nan turns every np.average below into nan.
                cata_tmp = cata_tmp.loc[
                                np.isfinite(cata_tmp[f'{fit_model}_g_1'].values)
                                & np.isfinite(cata_tmp[f'{fit_model}_g_2'].values), 
                                ['shear_type', 
                                f'{fit_model}_s2n', 
                                f'{fit_model}_g_1',
                                f'{fit_model}_g_2', 
                                f'{fit_model}_T', 
                                f'{fit_model}_T_ratio',
                                f'{fit_model}_band_flux']].copy()
                ## No weighting
                cata_tmp['weight'] = 1

            ## Add input info
            cata_tmp['g1_input'] = g1_input_list[i_shear]
            cata_tmp['g2_input'] = g2_input_list[i_shear]
            cata_tmp['run_tag'] = shear_tag        
            cata_tmp['tile_label'] = re.search(r'tile(.*)_rot', os.path.basename(inpath)).group(1)

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
    ## the concatenated catalogue is large, do not keep it across folders
    del cata, g1_input_all, g2_input_all, g1_measured_all, g2_measured_all


# >>>>>>>>>>>>>> shapes_metadetect (1/2), weight = shear_weight
# >>> Number of catalogues found in m283m283: 100
# >>> Number of catalogues found in m283p283: 100
# >>> Number of catalogues found in p283m283: 100
# >>> Number of catalogues found in p283p283: 100
# >>> Number of objects for m283m283: 17975943
# >>> R11, R22, R 0.2306168878572005 0.22687094800815916 0.22874391793267984
# >>> Number of objects for m283p283: 17976077
# >>> R11, R22, R 0.22977144424778856 0.2273988790131148 0.22858516163045167
# >>> Number of objects for p283m283: 17975751
# >>> R11, R22, R 0.22944450008164194 0.22756506320100003 0.228504781641321
# >>> Number of objects for p283p283: 17975151
# >>> R11, R22, R 0.229372982267951 0.22685700704574294 0.22811499465684698
# >>>> Total number of points for fitting 200 200 200 200
# m1 = -0.0011 pm 0.0036, c1 = -0.0001 pm 0.0001
# m2 = -0.0020 pm 0.0043, c2 = 0.0004 pm 0.0001

# >>>>>>>>>>>>>> shapes_metadetect_shifted (2/2), weight = shear_weight
# >>> Number of catalogues found in m283m283: 100
# >>> Number of catalogues found in m283p283: 100
# >>> Number of catalogues found in p283m283: 100
# >>> Number of catalogues found in p283p283: 100
# >>> Number of objects for m283m283: 17976537
# >>> R11, R22, R 0.22937246414765394 0.22748795364461838 0.22843020889613616
# >>> Number of objects for m283p283: 17976735
# >>> R11, R22, R 0.22940504681790316 0.22768304056454236 0.22854404369122278
# >>> Number of objects for p283m283: 17975469
# >>> R11, R22, R 0.22878005518777492 0.2269085663838761 0.2278443107858255
# >>> Number of objects for p283p283: 17975583
# >>> R11, R22, R 0.2289008175151937 0.22744213786328468 0.22817147768923918
# >>>> Total number of points for fitting 200 200 200 200
# m1 = -0.0004 pm 0.0037, c1 = -0.0000 pm 0.0001
# m2 = -0.0010 pm 0.0042, c2 = 0.0003 pm 0.0001
# finished  : Wed Sep  2 18:16:02 PDT 2026
