# -*- coding: utf-8 -*-
# @Author: lshuns
# @Date:   2022-05-31 15:06:40
# @Last Modified by:   lshuns
# @Last Modified time: 2022-09-26 17:26:26

### correct alpha in variance with method C
############ (20*20 resolution and SNR bins)
############ direct correction
############ assume scalelength and R cut already applied

import os
import argparse

import pandas as pd 
import numpy as np

import statsmodels.api as sm

# +++++++++++++++++++++++++++++ parser for command-line interfaces
parser = argparse.ArgumentParser(
    description=f"step1_methodC.py: correct alpha in variance with method C.",
    formatter_class=argparse.RawTextHelpFormatter)
parser.add_argument(
    "--inpath", type=str,
    help="the in path for the catalogue.\n\
    supported formats: feather")
parser.add_argument(
    "--outDir", type=str, 
    help="directory for the final catalogue.\n\
    outpath name will be inpath_name + _A1")
parser.add_argument(
    "--col_weight", type=str,
    help="columns to the weight in the catalogue.")
parser.add_argument(
    "--col_var", type=str,
    help="columns to the variance in the catalogue.")
parser.add_argument(
    "--col_snr", type=str,
    help="columns to the SNR in the catalogue.")
parser.add_argument(
    "--cols_e12", type=str, nargs=2, 
    help="column names for e1_gal, e2_gal.")
parser.add_argument(
    "--cols_e12_corr", type=str, nargs=2, default=None,
    help="column names for e1, e2, correction.")
parser.add_argument(
    "--cols_e12_raw", type=str, nargs=2, default=None,
    help="column names for e1, e2, raw measurement.\n\
    this or cols_e12_corr should be exclusive.")
parser.add_argument(
    "--cols_psf_e12", type=str, nargs=2, 
    help="column names for e1_psf, e2_psf.")

## arg parser
args = parser.parse_args()
inpath = args.inpath
outDir = args.outDir
outpath = os.path.join(outDir, os.path.basename(inpath).replace('.feather', '_A1.feather'))

col_weight = args.col_weight
col_var = args.col_var
col_snr = args.col_snr
col_e1, col_e2 = args.cols_e12

if args.cols_e12_corr is not None:
    col_e1_corr, col_e2_corr = args.cols_e12_corr
elif args.cols_e12_raw is not None:
    col_e1_raw, col_e2_raw = args.cols_e12_raw
else:
    raise Exception('cols_e12_corr or cols_e12_raw has to provide one!')

col_psf_e1, col_psf_e2 = args.cols_psf_e12
del args

# +++++++++++++++++++++++++++++ workhorse
# >>>>>>>>>>>>>>>>>>> values for variance to weight
#  set discretisation correction (see flensfit)
tinterval = 0.02
tintervalsq = tinterval*tinterval
# give prior weight quantities
priorsigma = 0.253091
# give prior weight quantities
priormoment = 2.*priorsigma*priorsigma
efunc_emax = 0.804
maxmoment = efunc_emax*efunc_emax/2.
maxweight = 2*(maxmoment-tintervalsq)/(tintervalsq*maxmoment + priormoment*(maxmoment-tintervalsq))
minweight = 0

# >>>>>>>>>>>>>>>>>>> load info
cata = pd.read_feather(inpath)
print('number total', len(cata))
## all float32 to float64
cata[cata.select_dtypes(np.float32).columns] = cata.select_dtypes(np.float32).astype(np.float64)

# get the raw e
if 'col_e1_corr' in locals():
    cata.loc[:, 'raw_e1'] = np.array(cata[col_e1]-cata[col_e1_corr])
    cata.loc[:, 'raw_e2'] = np.array(cata[col_e2]-cata[col_e2_corr])
else:
    cata.loc[:, 'raw_e1'] = np.array(cata[col_e1_raw])
    cata.loc[:, 'raw_e2'] = np.array(cata[col_e2_raw])

# a unique id for easy merge
cata.loc[:, 'AlphaRecal_index'] = np.arange(len(cata))

# >>>>>>>>>>>>>>>>>>>> binning

# quantile bin in SNR and R
## number of bins
N_R = 20
N_SNR = 20

# start with R bin
cata.loc[:, 'bin_R'] = pd.qcut(cata['R'].values, N_R, 
                                    labels=False, retbins=False)

# initialize
cata.loc[:, 'bin_snr'] = -999

# in each R bin, do SNR binning
for ibin_R in range(N_R):

    # select catalogue
    mask_binR = (cata['bin_R'].values == ibin_R)

    # bin in R
    cata.loc[mask_binR, 'bin_snr'] = pd.qcut(cata.loc[mask_binR, col_snr].values, N_SNR, 
                                                labels=False, retbins=False)

# group based on binning
cata_corr = []
for name, group in cata.groupby(by=['bin_R', 'bin_snr']):

    # >>>>>>>>>>>>>>>> calculate alpha
    # unique index
    index_obj = group['AlphaRecal_index'].values
    # variance
    Var = np.array(group[col_var])
    # out shear
    e1_out = np.array(group['raw_e1'])
    e2_out = np.array(group['raw_e2'])
    emod_out = np.hypot(e1_out, e2_out)
    # out PSF 
    e1_psf = np.array(group[col_psf_e1])
    e2_psf = np.array(group[col_psf_e2])
    del group
    # rotate to PSF frame
    mask_tmp = emod_out > 0
    e_psf_rotated = e1_out * e1_psf + e2_out * e2_psf
    e_psf_rotated[mask_tmp] /= emod_out[mask_tmp]
    del mask_tmp, e1_psf, e2_psf, e1_out, e2_out, emod_out
    # calculate alpha using least square
    mod_ols = sm.OLS(Var, sm.add_constant(e_psf_rotated))
    res_ols = mod_ols.fit()
    alpha = res_ols.params[1]
    # alpha_err = (res_ols.cov_params()[1, 1])**0.5

    # >>>>>>>>>>>>>>>> correct 
    ## get logvar correction
    Var_corr = Var - alpha * e_psf_rotated
    del e_psf_rotated, Var

    # >>>>>>>>>>>>>>>> save
    cata_tmp = pd.DataFrame({'AlphaRecal_index': index_obj, 
                            'AlphaRecalC_variance': Var_corr})
    del index_obj, Var_corr
    cata_corr.append(cata_tmp)
    del cata_tmp
cata_corr = pd.concat(cata_corr)

# transfer back to weights
## hard cut on min var
cata_corr.loc[:, 'AlphaRecalC_variance'] = np.maximum(cata_corr['AlphaRecalC_variance'].values, tintervalsq)
## now convert these new 2Dvariances into weights
cata_corr.loc[:, 'AlphaRecalC_weight'] = 2*(maxmoment - cata_corr['AlphaRecalC_variance'].values)\
                / (cata_corr['AlphaRecalC_variance'].values*maxmoment \
                    + priormoment*(maxmoment - cata_corr['AlphaRecalC_variance'].values))
## post process of weights
cata_corr.loc[(cata_corr['AlphaRecalC_weight']>maxweight), 'AlphaRecalC_weight'] = maxweight
cata_corr.loc[(cata_corr['AlphaRecalC_weight']<maxweight/1000.), 'AlphaRecalC_weight'] = 0.

# merge
cata = cata.merge(cata_corr, on='AlphaRecal_index', how='left')

# if raw weight is zero, keep it zero
cata.loc[cata[col_weight]<maxweight/1000., 'AlphaRecalC_weight'] = 0

# save
cata.to_feather(outpath)
print('final results saved to', outpath)
