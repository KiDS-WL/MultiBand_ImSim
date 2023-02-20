# -*- coding: utf-8 -*-
# @Author: lshuns
# @Date:   2022-10-27 14:13:20
# @Last Modified by:   lshuns
# @Last Modified time: 2023-02-20 10:47:50

### correct alpha in measured e with method D
##### weight should be corrected already
##### final catalogue only contains objects with non-zero weight and within given Z_B_edges
##### if the redshift calibration flag is provided, only gold class will be saved
##### method D contains two steps:
############ 1. general fitting to remove the main trend in resolution and SNR
############ 2. tomographic bins + 20*20 SNR and resolution bins to remove residual

import os
import time
import argparse

import numpy as np
import pandas as pd 
import numpy.linalg as la
import statsmodels.api as sm 

from astropy.io import fits
from astropy.table import Table

# +++++++++++++++++++++++++++++ parser for command-line interfaces
parser = argparse.ArgumentParser(
    description=f"step2_methodD.py: correct alpha in e1,2 with method D.",
    formatter_class=argparse.RawTextHelpFormatter)
parser.add_argument(
    "--inpath", type=str,
    help="the in path for the catalogue.\n\
    supported formats: feather, fits, LDAC cat")
parser.add_argument(
    "--outDir", type=str, 
    help="directory for the final catalogue.\n\
    outpath name will be inpath_name + .A2.feather")
parser.add_argument(
    "--col_goldFlag", type=str,
    help="columns to the redshift gold class in the catalogue.")
parser.add_argument(
    "--col_weight", type=str,
    help="columns to the weight in the catalogue.")
parser.add_argument(
    "--col_snr", type=str,
    help="columns to the SNR in the catalogue.")
parser.add_argument(
    "--col_ZB", type=str,
    help="columns to the Z_B in the catalogue.")
parser.add_argument(
    "--cols_e12", type=str, nargs=2, 
    help="column names for e1_gal, e2_gal.")
parser.add_argument(
    "--cols_psf_e12", type=str, nargs=2, 
    help="column names for e1_psf, e2_psf.")
parser.add_argument(
    "--Z_B_edges", type=float, nargs='*', default=None,
    help="edges for tomographic binning.")

## arg parser
args = parser.parse_args()
inpath = args.inpath
outDir = args.outDir
outpath = os.path.join(outDir, inpath + '.A2.feather')

col_goldFlag = args.col_goldFlag
col_weight = args.col_weight
col_snr = args.col_snr
col_ZB = args.col_ZB
col_e1, col_e2 = args.cols_e12
col_psf_e1, col_psf_e2 = args.cols_psf_e12
Z_B_edges = np.array(args.Z_B_edges)
del args

# >>>>>>>>>>>>>>>>>>>>> workhorse

# ++++++ 0. load catalogue
file_type = inpath[-3:]
if file_type == 'her':
    cata = pd.read_feather(inpath)
elif file_type == 'cat':
    with fits.open(inpath) as hdul:
        cata = Table(hdul['OBJECTS'].data).to_pandas()
elif file_type == 'its':
    with fits.open(inpath) as hdul:
        cata = Table(hdul[1].data).to_pandas()
else:
    raise Exception(f'Not supported input file type! {inpath}')
print('number original', len(cata))
## all float32 to float64
cata[cata.select_dtypes(np.float32).columns] = cata.select_dtypes(np.float32).astype(np.float64)
## only gold class
if col_goldFlag is not None:
    cata = cata[(cata[col_goldFlag].values>0)]
    cata.reset_index(drop=True, inplace=True)
    print('number after gold selection', len(cata))
## only preserve weight > 0
cata = cata[(cata[col_weight].values>0)]
cata.reset_index(drop=True, inplace=True)
print('number after weight selection', len(cata))
## those out of ZB range are removed
cata = cata[(cata[col_ZB].values>Z_B_edges[0])&(cata[col_ZB].values<=Z_B_edges[-1])]
cata.reset_index(drop=True, inplace=True)
print('number after Z_B selection', len(cata))

# ++++++ 1. get alpha map 
start_time = time.time()

## number of bins
N_R = 20
N_SNR = 20

## start with R bin
cata.loc[:, 'bin_R'] = pd.qcut(cata['R'], N_R, 
                                    labels=False, retbins=False)

## then snr bin
cata.loc[:, 'bin_snr'] = -999
for ibin_R in range(N_R):

    # select catalogue
    mask_binR = (cata['bin_R'].values == ibin_R)

    # bin in R
    cata.loc[mask_binR, 'bin_snr'] = pd.qcut(cata.loc[mask_binR, col_snr].values, N_SNR, 
                                                labels=False, retbins=False)

## group based on binning
cata_grouped = cata.groupby(by=['bin_R', 'bin_snr'])

## get the 2D alpha map
cata_alpha = cata_grouped.sum()
cata_alpha = cata_alpha[[col_weight]].copy()
cata_alpha.rename(columns={col_weight: 'total_weights'}, inplace=True)
for name, group in cata_grouped:

    # out shear
    e1_out = np.array(group[col_e1])
    e2_out = np.array(group[col_e2])
    weight_out = np.array(group[col_weight])

    # out PSF 
    e1_psf = np.array(group[col_psf_e1])
    e2_psf = np.array(group[col_psf_e2])

    # save center
    cata_alpha.loc[name, 'R_mean_wei'] = np.average(group['R'].values, weights=weight_out)
    cata_alpha.loc[name, 'SNR_mean_wei'] = np.average(group[col_snr].values, weights=weight_out)
    del group

    # calculate alpha using least square
    ## e1
    mod_wls = sm.WLS(e1_out, sm.add_constant(e1_psf), weights=weight_out)
    res_wls = mod_wls.fit()
    cata_alpha.loc[name, 'alpha1'] = res_wls.params[1]
    cata_alpha.loc[name, 'alpha1_err'] = (res_wls.cov_params()[1, 1])**0.5
    del e1_out, e1_psf, mod_wls, res_wls
    ## e2
    mod_wls = sm.WLS(e2_out, sm.add_constant(e2_psf), weights=weight_out)
    res_wls = mod_wls.fit()
    cata_alpha.loc[name, 'alpha2'] = res_wls.params[1]
    cata_alpha.loc[name, 'alpha2_err'] = (res_wls.cov_params()[1, 1])**0.5
    del e2_out, e2_psf, weight_out, mod_wls, res_wls

del cata_grouped
print('alpha map produced in', time.time() - start_time, 's')

# ++++++ 2. step1: fitting in whole and removing general trend
start_time = time.time()

## training model
# e1
fitting_weight = 1./np.square(cata_alpha['alpha1_err'].values)
A = np.vstack([fitting_weight * 1, 
                fitting_weight * np.power(cata_alpha['SNR_mean_wei'].values, -2),
                fitting_weight * np.power(cata_alpha['SNR_mean_wei'].values, -3),
                fitting_weight * cata_alpha['R_mean_wei'].values,
                fitting_weight * cata_alpha['R_mean_wei'].values * np.power(cata_alpha['SNR_mean_wei'].values, -2)]).T
poly1 = la.lstsq(A, fitting_weight * cata_alpha['alpha1'].values, rcond=None)[0] # poly coefficients
del fitting_weight, A
# e2
fitting_weight = 1./np.square(cata_alpha['alpha2_err'].values)
A = np.vstack([fitting_weight * 1, 
                fitting_weight * np.power(cata_alpha['SNR_mean_wei'].values, -2),
                fitting_weight * np.power(cata_alpha['SNR_mean_wei'].values, -3),
                fitting_weight * cata_alpha['R_mean_wei'].values,
                fitting_weight * cata_alpha['R_mean_wei'].values * np.power(cata_alpha['SNR_mean_wei'].values, -2)]).T
poly2 = la.lstsq(A, fitting_weight * cata_alpha['alpha2'].values, rcond=None)[0] # poly coefficients
del fitting_weight, A, cata_alpha

## remove the general trend
# e1
alpha = poly1[0] \
        + poly1[1] * np.power(cata[col_snr].values, -2) \
        + poly1[2] * np.power(cata[col_snr].values, -3) \
        + poly1[3] * cata['R'].values \
        + poly1[4] * cata['R'].values * np.power(cata[col_snr].values, -2)
cata.loc[:, 'AlphaRecalD1_e1'] = cata[col_e1].values - alpha * cata[col_psf_e1].values
del alpha, poly1
# e2
alpha = poly2[0] \
        + poly2[1] * np.power(cata[col_snr].values, -2) \
        + poly2[2] * np.power(cata[col_snr].values, -3) \
        + poly2[3] * cata['R'].values \
        + poly2[4] * cata['R'].values * np.power(cata[col_snr].values, -2)
cata.loc[:, 'AlphaRecalD1_e2'] = cata[col_e2].values - alpha * cata[col_psf_e2].values
del alpha, poly2

## meaningful e
mask_tmp = (cata['AlphaRecalD1_e1'].values>-1) & (cata['AlphaRecalD1_e1'].values<1) \
           & (cata['AlphaRecalD1_e2'].values>-1) & (cata['AlphaRecalD1_e2'].values<1)
print('number with meaningful e after D1', np.sum(mask_tmp), 'fraction', np.sum(mask_tmp)/len(cata))
cata = cata[mask_tmp]
cata.reset_index(drop=True, inplace=True)
del mask_tmp

print('D1 finished in', time.time() - start_time, 's')

# ++++++ 3. step2: direct correction for residual alpha
start_time = time.time()

# start with ZB bins
cata.loc[:, 'bin_ZB'] = pd.cut(cata[col_ZB], Z_B_edges, 
                                    right=True, labels=False)

# then R and SNR
cata.loc[:, 'bin_R'] = -999
cata.loc[:, 'bin_snr'] = -999
for ibin_ZB in range(len(Z_B_edges)-1):

    # select catalogue
    mask_binZB = np.array(cata['bin_ZB']) == ibin_ZB

    # bin in R
    cata.loc[mask_binZB, 'bin_R'] = pd.qcut(cata.loc[mask_binZB, 'R'].values, N_R, 
                                    labels=False, retbins=False)

    # in each R bin, do SNR binning
    for ibin_R in range(N_R):

        # select catalogue
        mask_binR = (cata['bin_R'].values == ibin_R)

        # bin in R
        cata.loc[mask_binZB&mask_binR, 'bin_snr'] = pd.qcut(cata.loc[mask_binZB&mask_binR, col_snr].values, N_SNR, 
                                                    labels=False, retbins=False)
        del mask_binR
    del mask_binZB

# correct in each bin
cata_corr = []
for name, group in cata.groupby(by=['bin_ZB', 'bin_R', 'bin_snr']):

    # >>>>>>>>>>>>>>>> calculate alpha
    # unique index
    index_obj = group['AlphaRecal_index'].values
    # out shear
    e1_out = np.array(group['AlphaRecalD1_e1'])
    e2_out = np.array(group['AlphaRecalD1_e2'])
    weight_out = np.array(group[col_weight])
    # out PSF 
    e1_psf = np.array(group[col_psf_e1])
    e2_psf = np.array(group[col_psf_e2])
    del group
    # calculate alpha using least square
    ## e1
    mod_wls = sm.WLS(e1_out, sm.add_constant(e1_psf), weights=weight_out)
    res_wls = mod_wls.fit()
    alpha1 = res_wls.params[1]
    del res_wls, mod_wls
    ## e2
    mod_wls = sm.WLS(e2_out, sm.add_constant(e2_psf), weights=weight_out)
    res_wls = mod_wls.fit()
    alpha2 = res_wls.params[1]
    del weight_out, res_wls, mod_wls

    # >>>>>>>>>>>>>>>> correct 
    e1_corr = e1_out - alpha1 * e1_psf
    e2_corr = e2_out - alpha2 * e2_psf
    del e1_out, alpha1, e1_psf, e2_out, alpha2, e2_psf

    # >>>>>>>>>>>>>>>> save
    cata_tmp = pd.DataFrame({'AlphaRecal_index': index_obj, 
                            'AlphaRecalD2_e1': e1_corr,
                            'AlphaRecalD2_e2': e2_corr})
    del index_obj, e1_corr, e2_corr
    cata_corr.append(cata_tmp)
    del cata_tmp
cata_corr = pd.concat(cata_corr)

# meaningful e
mask_tmp = (cata_corr['AlphaRecalD2_e1']>-1) & (cata_corr['AlphaRecalD2_e1']<1) \
           & (cata_corr['AlphaRecalD2_e2']>-1) & (cata_corr['AlphaRecalD2_e2']<1)
print('number with meaningful e after D2', np.sum(mask_tmp), 'fraction', np.sum(mask_tmp)/len(cata_corr))
cata_corr = cata_corr[mask_tmp]
del mask_tmp

# merge
cata = cata.merge(cata_corr, on='AlphaRecal_index')
cata.reset_index(drop=True, inplace=True)
del cata_corr
print('D2 finished in', time.time() - start_time, 's')

# save
cata.to_feather(outpath)
print('number in final cata', len(cata))
print('final results saved to', outpath)