# -*- coding: utf-8 -*-
# @Author: lshuns
# @Date:   1969-12-31 16:00:00
# @Last Modified by:   lshuns
# @Last Modified time: 2026-01-13 13:36:09

### Assign shear weights to the metadetect catalogues 
###### follow the approach of https://ui.adsabs.harvard.edu/abs/2023OJAp....6E..17S/abstract

import os
import glob

import numpy as np 
import pandas as pd

## ++++++++++++++ I/O and general setups

## Where to find the simulations
main_dir = '/sdf/data/kipac/u/liss/ImSim/output/test_dev/LSST_r/'

## Shear inputs in simulations
shear_tags = ['m283m283', 'm283p283', 'p283m283', 'p283p283']

## What is the fitting model used in metadetect
fit_model = 'wmom'

## Cut info
snr_min = 12.5
resolution_min = 1.2

## The intrinsic ellipticity dispersion
sigma_SN = 0.07

## ++++++++++++++ Workhorse

## Loop over catalogues
for i_shear, shear_tag in enumerate(shear_tags):
    inpath_list = glob.glob(os.path.join(main_dir, 
                                         shear_tag, 
                                         'catalogues/shapes_metadetect', 
                                         '*.feather'))
    print(f">>> Number of catalogues found in {shear_tag}: {len(inpath_list)}")

    for inpath in inpath_list:
        cata = pd.read_feather(inpath)

        ## Weight based on ellipticity noise
        cata['weight_sigma_e'] = 1./(sigma_SN**2 
                                    + np.square(
                                        cata[f'{fit_model}_g_cov_as_sigma'].values
                                        )
                                    )
        
        ## Weight based on ellipticity
        e_sq = np.square(cata[f'{fit_model}_g_1'].values) \
            + np.square(cata[f'{fit_model}_g_2'].values)
        cata['weight_e'] = np.square(1-e_sq) * np.exp(-1*e_sq/2./0.09)

        ## Combine them together
        cata['shear_weight'] = cata['weight_sigma_e'] * cata['weight_e']

        ## S/N and resolution cut
        mask_cut = ((cata[f'{fit_model}_s2n'] <= snr_min)
                 & (cata[f'{fit_model}_T_ratio'] <= resolution_min)
                    )
        cata.loc[mask_cut, 
                'weight_sigma_e'] = 0.
        cata.loc[mask_cut, 
                'weight_e'] = 0.
        cata.loc[mask_cut, 
                'shear_weight'] = 0.
        
        ## Zero weight for nan measurements
        cata[['weight_sigma_e', 
              'weight_e', 
              'shear_weight']] = cata[['weight_sigma_e', 
                                        'weight_e', 
                                        'shear_weight']].fillna(value=0.)

        ## Save back
        cata.to_feather(inpath)
        print(f"++++++ Finished for {os.path.basename(inpath)}, zero weight ratio: {np.sum(cata['shear_weight']==0)/len(cata)}")
        del cata