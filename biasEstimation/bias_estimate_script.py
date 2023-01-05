# -*- coding: utf-8 -*-
# @Author: lshuns
# @Date:   2022-04-13 14:16:35
# @Last Modified by:   lshuns
# @Last Modified time: 2022-11-08 17:00:30

### script to run the bias calculation

######## Supported biases (bias_type):
########        1. alpha_LS: the PSF leakage term (alpha) from least squares
########        2. m_LS_pair: the shear bias from least squares with pair based noise cancellation
########        3. m_LS_tile: the shear bias from least squares with tile based noise cancellation

import argparse

import numpy as np
import pandas as pd
import statsmodels.api as sm 

from scipy import stats
from astropy.io import fits
from astropy.table import Table

# self package
from bias_estimate_func import alphaCalFunc_least_squares, mCalFunc_pair_based, mCalFunc_tile_based

# +++++++++++++++++++++++++++++ parser for command-line interfaces
parser = argparse.ArgumentParser(
    description=f"bias_estimate_script.py: estimate the bias from catalogues.",
    formatter_class=argparse.RawTextHelpFormatter)
parser.add_argument(
    "--bias_type", type=str, choices=['alpha_LS', 'm_LS_pair', 'm_LS_tile'],
    help="Select the bias type to be calculated:\n\
    alpha_LS: the PSF leakage term (alpha) from least squares\n\
    m_LS_pair: the shear bias from least squares with pair-based average\n\
    m_LS_tile: the shear bias from least squares with tile-based average\n\
    ")
parser.add_argument(
    "--in_file", type=str,
    help="the input catalogue.\n\
    supported formats: fits, csv, feather\n\
        guess from the file suffix")
parser.add_argument(
    "--out_path", type=str, default='./test.txt',
    help="where to save the final results. \n\
    Default ./test.txt.")
parser.add_argument(
    "--col_goldFlag", type=str, default=None,
    help="columns to the redshift gold class in the catalogue.")
parser.add_argument(
    "--cols_e12", type=str, nargs=2, default=['e1_LF_r', 'e2_LF_r'], 
    help="column names for e1_gal, e2_gal.")
parser.add_argument(
    "--e_type", type=str, choices=['input', 'measured'], 
    help="what is the e values: input OR measured.")
parser.add_argument(
    "--cols_g12", type=str, nargs=2, default=['g2_in', 'g1_in'], 
    help="column names for g1_in, g2_in.")
parser.add_argument(
    "--cols_PSFe12", type=str, nargs=2, default=['psf_e1_LF_r', 'psf_e2_LF_r'], 
    help="column names for e1_PSF, e2_PSF.")
parser.add_argument(
    "--col_weight", type=str, default=None,
    help="columns to the weight.\n\
    If not provided, equal weight for all galaxies.")
parser.add_argument(
    "--col_label", type=str, default=None,
    help="columns to the tile label.")
parser.add_argument(
    "--col_id_input", type=str, default=None,
    help="columns to the unique id of input galaxies.")
parser.add_argument(
    "--col_binning", type=str,
    help="columns to be used for binning.\n\
    If not provided, no binning will be applied.")
parser.add_argument(
    "--binning_edges", type=float, nargs='*', default=None,
    help="edges for binning, \n\
            should not set along with binning_min_max_Nbin.")
parser.add_argument(
    "--binning_min_max_Nbin", type=float, nargs=3, default=None,
    help="min max N_bins for binning.\n\
            should not set along with binning_edges.")
parser.add_argument(
    "--binning_type", type=str, choices=['linear', 'log', 'quantile'], 
    help=" Binning type: linear, log, quantile.\n\
            Only used if binning_min_max_Nbin is provided.")
parser.add_argument(
    "--psf_frame", action="store_true",
    help="rotate to the PSF frame.")

## arg parser
args = parser.parse_args()
bias_type = args.bias_type
in_file = args.in_file
out_path = args.out_path

col_goldFlag = args.col_goldFlag

col_e1, col_e2 = args.cols_e12
e_type = args.e_type

col_g1, col_g2 = args.cols_g12

col_e1_psf, col_e2_psf = args.cols_PSFe12

col_weight = args.col_weight

col_label = args.col_label
col_id_input = args.col_id_input

col_binning = args.col_binning
if col_binning is not None:
    print('binning on:', col_binning)
    if (args.binning_edges is not None) and (args.binning_min_max_Nbin is not None):
        raise Exception('cannot assign both binning_edges and binning_min_max_Nbin!')
    if args.binning_edges is not None:
        bin_edges = np.array(args.binning_edges)
        print('provided bins', bin_edges)
    elif args.binning_min_max_Nbin is not None:
        min_bin, max_bin, N_bin = args.binning_min_max_Nbin
        binning_type = args.binning_type
        if binning_type == 'linear':
            bin_edges = np.linspace(min_bin, max_bin, int(N_bin))
            ## make sure edge values are included
            bin_edges[0] = min_bin - 0.001
            bin_edges[-1] = max_bin + 0.001
            print('linear bins', bin_edges)
        elif binning_type == 'log':
            bin_edges = np.logspace(np.log10(min_bin), np.log10(max_bin), int(N_bin))
            ## make sure edge values are included
            bin_edges[0] = min_bin - 0.001
            bin_edges[-1] = max_bin + 0.001
            print('log bins', bin_edges)
    else:
        raise Exception('either binning_edges or binning_min_max_Nbin should be assigned for col_binning==True!')

psf_frame = args.psf_frame
del args

# +++++++++++++++++++++++++++++ workhorse

# load working function
if bias_type == 'alpha_LS':
    CalFunc = alphaCalFunc_least_squares
elif bias_type == 'm_LS_pair':
    if col_id_input is not None:
        CalFunc = mCalFunc_pair_based
        print('using pair based mCalFunc')
    else:
        raise Exception('m_LS_pair chosen, but no col_id_input provided!')
elif bias_type == 'm_LS_tile':
    if col_label is not None:
        CalFunc = mCalFunc_tile_based
        print('using tile based mCalFunc')
    else:
        raise Exception('m_LS_tile chosen, but no col_label provided!')

# load simulation catalogue
file_type = in_file[-3:]
if file_type == 'csv':
    cata_sim = pd.read_csv(in_file)
elif file_type == 'her':
    cata_sim = pd.read_feather(in_file)
elif file_type == 'its':
    with fits.open(in_file) as hdul:
        cata_sim = hdul[1].data
else:
    raise Exception(f'Not supported input file type! {in_file}')
print('Number of sources in the catalogue', len(cata_sim))
### select gold class
if col_goldFlag is not None:
    cata_sim = cata_sim[(cata_sim[col_goldFlag].values>0)]
    cata_sim.reset_index(drop=True, inplace=True)
    print('number after gold selection', len(cata_sim))
### select those within the edges
if ('bin_edges' in locals()):
    cata_sim = cata_sim[(cata_sim[col_binning]>bin_edges[0])&(cata_sim[col_binning]<=bin_edges[-1])]
    print('selected objects (within bin edge)', len(cata_sim))

# get shear values
if e_type == 'input':
    # perfect shear values from input e
    g = np.array(cata_sim[col_g1]) + 1j*np.array(cata_sim[col_g2])
    e_in_gal = np.array(cata_sim[col_e1]) + 1j*np.array(cata_sim[col_e2])
    e_true = (e_in_gal+g) / (1+np.conj(g)*e_in_gal)

    e1_out = (e_true.real).astype(float)
    e2_out = (e_true.imag).astype(float)

elif e_type == 'measured':
    # measured values
    e1_out = np.array(cata_sim[col_e1])
    e2_out = np.array(cata_sim[col_e2])

else:
    raise Exception(f'unsupported e_type {e_type}')

# save useable parameters
cata_used = pd.DataFrame({
    'e1_out': e1_out.astype(float),
    'e2_out': e2_out.astype(float)})
## input shear
try:
    cata_used.loc[:, 'g1_in'] = np.array(cata_sim[col_g1]).astype(float)
    cata_used.loc[:, 'g2_in'] = np.array(cata_sim[col_g2]).astype(float)
except KeyError:
    pass
## psf e
try:
    cata_used.loc[:, 'e1_psf'] = np.array(cata_sim[col_e1_psf]).astype(float)
    cata_used.loc[:, 'e2_psf'] = np.array(cata_sim[col_e2_psf]).astype(float)
except KeyError:
    pass
## binning values
if col_binning is not None:
    cata_used.loc[:, col_binning] = np.array(cata_sim[col_binning]).astype(float)
## weight for galaxies
if col_weight is not None:
    cata_used.loc[:, 'shape_weight'] = np.array(cata_sim[col_weight]).astype(float)
else:
    cata_used.loc[:, 'shape_weight'] = 1
## IDs
if col_id_input is not None:
    cata_used.loc[:, 'id_input'] = np.array(cata_sim[col_id_input]).astype(int)
## labels
if col_label is not None:
    cata_used.loc[:, 'tile_label'] = np.array(cata_sim[col_label]).astype(str)
## delete original catalogue
del cata_sim

# select non-zero weights
cata_used = cata_used[cata_used['shape_weight']>0]
cata_used.reset_index(drop=True, inplace=True)
print('selected objects (weight>0)', len(cata_used))

# output
f = open(out_path, 'w')

# binning
if col_binning is not None:
    ## bin on quantile
    if (not 'bin_edges' in locals()) and (binning_type == 'quantile'):
        quartiles = np.linspace(0, 1, int(N_bin+1))
        bin_edges = stats.mstats.mquantiles(cata_used[col_binning].values, quartiles)
        ## make sure edge values are included
        bin_edges[0] = min_bin - 0.001
        bin_edges[-1] = max_bin + 0.001
        print('quantile bins', bin_edges)

        cata_used = cata_used[(cata_used[col_binning]>bin_edges[0]) & (cata_used[col_binning]<=bin_edges[-1])]
        cata_used.reset_index(drop=True, inplace=True)
        print('selected objects (within bin edge)', len(cata_used))

    ## the whole sample only contains objects within the bin edge
    min_bin = bin_edges[0]
    max_bin = bin_edges[-1]

    # mean and median of the x axis
    mean_bin = np.average(cata_used[col_binning].values, weights=cata_used['shape_weight'].values)
    median_bin = np.median(cata_used[col_binning].values)

    ## the whole results
    res = CalFunc(cata_used, psf_frame=psf_frame)
    N_s = len(cata_used)
    Nwei = np.sum(cata_used['shape_weight'].values)

    # collect columns names and values
    cols = ','.join(list(res.keys()))
    vals = ','.join(["{0:0.4f}".format(val) for val in res.values()])
    if col_binning is not None:
        cols = cols + f',{col_binning}_min,{col_binning}_max,{col_binning}_mean,{col_binning}_median,Nobj,Nwei'
        vals = vals + f',{min_bin},{max_bin},{mean_bin},{median_bin},{N_s},{Nwei}'    
    print(cols, file=f)
    print(vals, file=f)

    for i_bin in range(len(bin_edges)-1):
        min_bin = bin_edges[i_bin]
        max_bin = bin_edges[i_bin+1]
        mask_bin = (cata_used[col_binning]>min_bin) & (cata_used[col_binning]<=max_bin)

        cata_selec = cata_used[mask_bin].copy()
        del mask_bin
        cata_selec.reset_index(drop=True, inplace=True)

        # mean and median of the x axis
        mean_bin = np.average(cata_selec[col_binning].values, weights=cata_selec['shape_weight'].values)
        median_bin = np.median(cata_selec[col_binning].values)

        res = CalFunc(cata_selec, psf_frame=psf_frame)
        N_s = len(cata_selec)
        Nwei = np.sum(cata_selec['shape_weight'].values)
        del cata_selec

        # print out values
        vals = ','.join(["{0:0.4f}".format(val) for val in res.values()]) \
        + f',{min_bin},{max_bin},{mean_bin},{median_bin},{N_s},{Nwei}'
        print(vals, file=f)

else:
    #only whole
    res = CalFunc(cata_used, psf_frame=psf_frame)
    N_s = len(cata_used)
    Nwei = np.sum(cata_used['shape_weight'].values)

    # collect columns names and values
    cols = ','.join(list(res.keys()))
    vals = ','.join(["{0:0.4f}".format(val) for val in res.values()])
    cols = cols + f',bin_min,bin_max,bin_mean,bin_median,Nobj,Nwei'
    vals = vals + f',-999,-999,-999,-999,{N_s},{Nwei}'    
    print(cols, file=f)
    print(vals, file=f)

# close what is opened
f.close()
print(f'results saved to {out_path}')
