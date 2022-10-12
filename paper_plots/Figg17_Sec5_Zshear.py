# -*- coding: utf-8 -*-
# @Author: lshuns
# @Date:   2022-06-09 18:12:58
# @Last Modified by:   lshuns
# @Last Modified time: 2022-10-04 10:03:52

### plot the shear amplitude as a function of redshift

import os
import glob

import numpy as np
import pandas as pd

import plotting

# >>>>>>>>>>>>>>> I/O

# the catalogue with diferent components of shears
infile = '/disks/shear10/ssli/ImSim/input/SURFS_cata/varShear/skills_v07Ds_input_part0_blended4_magCut25_varShear.feather'
cata_final = pd.read_feather(infile)

# >>>>>>>>>>>>>>> calculate mean as a function of redshift

# redshift array
zbins_edge = np.arange(0.2, 2.6, 0.2)
zbins_center = (zbins_edge[1:]+zbins_edge[:-1])/2.

# shear array
shear_CS = np.zeros_like(zbins_center)
# shear_CS_err = [np.zeros_like(zbins_center), np.zeros_like(zbins_center)]
shear_ggl = np.zeros_like(zbins_center)
# shear_ggl_err = [np.zeros_like(zbins_center), np.zeros_like(zbins_center)]
shear_all = np.zeros_like(zbins_center)
# shear_all_err = [np.zeros_like(zbins_center), np.zeros_like(zbins_center)]
## bulid a list for loop
shear_list = [shear_CS, shear_ggl, shear_all]
# shear_err_list = [shear_CS_err, shear_ggl_err, shear_all_err]
shear_suffix_list = ['CS0', 'ggl', 'fromCS0']

# the mean
# nboot = 500
for i_zbin in range(len(zbins_center)):
    zmin = zbins_edge[i_zbin]
    zmax = zbins_edge[i_zbin + 1]
    mask_tmp = (cata_final['zobs'].values >= zmin) & (cata_final['zobs'].values < zmax)

    ## shears
    for i_shear, shear_res in enumerate(shear_list):
        # get the amplitude
        shear_suffix = shear_suffix_list[i_shear]
        shear_amp = np.hypot(cata_final.loc[mask_tmp, f'gamma1_{shear_suffix}'].values, cata_final.loc[mask_tmp, f'gamma2_{shear_suffix}'].values)

        # get the average and errors
        shear_mean = np.nanmean(shear_amp)
        shear_res[i_zbin] = shear_mean
        # ## get error from boots
        # Relications_index = np.array([np.random.randint(len(shear_amp), size=len(shear_amp)) for _ in range(nboot)])
        # mean_BS = np.nanmean(shear_amp[Relications_index], axis = 1)
        del shear_amp
    #     sigma_BS = mean_BS - shear_mean
    #     err_low = np.abs(np.percentile(sigma_BS, 16))
    #     err_high = np.abs(np.percentile(sigma_BS, 84))
    #     ## assign
    #     shear_err_list[i_shear][0][i_zbin] = err_low
    #     shear_err_list[i_shear][1][i_zbin] = err_high
    del mask_tmp

## final
outpath = 'show'
outpath = './plots/Zshear.pdf'

xvals = [zbins_center] * 3
yvals = shear_list
COLORs = ['magenta', 'orange', 'black']
LINEs = ['--', ':', '-']
POINTs = ['', '', '']
LABELs = ['CS', 'GGL', 'All']
LINEWs = [1.5, 1.5, 1.5]

XLABEL = 'Redshift'

YLABEL = r'Mean $|\boldsymbol{\gamma}|=\sqrt{\gamma_1^2+\gamma_2^2}$'

texPacks = [r"\usepackage{amsmath}"]

plotting.LinePlotFunc(outpath,
                xvals, yvals,
                COLORs, LABELs=LABELs, LINEs=LINEs, LINEWs=LINEWs, POINTs=POINTs, POINTSs=None, fillstyles=None,
                XRANGE=None, YRANGE=None,
                XLABEL=XLABEL, YLABEL=YLABEL, TITLE=None,
                xtick_min_label=True, xtick_spe=None, ytick_min_label=True, ytick_spe=None,
                vlines=None, vline_styles=None, vline_colors=None, vline_labels=None, vline_widths=None,
                hlines=None, hline_styles=None, hline_colors=None, hline_labels=None, hline_widths=None,
                xlog=False, invertX=False, ylog=False, invertY=False, 
                loc_legend='best', legend_frame=False,
                font_size=15, usetex=True,
                texPacks=texPacks)
