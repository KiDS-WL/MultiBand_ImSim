# -*- coding: utf-8 -*-
# @Author: lshuns
# @Date:   2022-07-05 14:44:10
# @Last Modified by:   lshuns
# @Last Modified time: 2022-10-13 14:32:49

### script to run the bias calculation
##### simulation reweighted with the data

######## Supported biases (bias_type):
########        2. m_LS_DataRewei_2D: the shear bias from least squares with data reweighting in 2 parameters

import argparse

import numpy as np
import pandas as pd
import statsmodels.api as sm 

from scipy import stats
from astropy.io import fits
from astropy.table import Table

# self package
from bias_estimate_func import mCalFunc_DataRewei_2D

# +++++++++++++++++++++++++++++ parser for command-line interfaces
parser = argparse.ArgumentParser(
    description=f"bias_estimate_script_dataRewei.py: estimate the bias from simulations with data reweighting.",
    formatter_class=argparse.RawTextHelpFormatter)
parser.add_argument(
    "--bias_type", type=str, choices=['m_LS_DataRewei_2D'],
    help="Select the bias type to be calculated:\n\
    m_LS_DataRewei_2D: the shear bias from least squares with data reweighting in 2 parameters\n\
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
    "--cols_e12", type=str, nargs=2, default=['e1_LF_r', 'e2_LF_r'], 
    help="column names for e1_gal, e2_gal.")
parser.add_argument(
    "--cols_g12", type=str, nargs=2, default=['g2_in', 'g1_in'], 
    help="column names for g1_in, g2_in.")
parser.add_argument(
    "--col_weight_sim_data", type=str, nargs=2,
    help="columns to the weight.\n\
    order: (sim, data)")
parser.add_argument(
    "--col_label", type=str, default=None,
    help="columns to the tile label.")
parser.add_argument(
    "--col_id_input", type=str, default=None,
    help="columns to the unique id of input galaxies.")
parser.add_argument(
    "--col_binning_sim_data", type=str, nargs=2,
    help="columns to be used for binning.\n\
    If not provided, no binning will be applied.\n\
    order: (sim, data)")
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
parser.add_argument(
    "--in_file_data", type=str,
    help="the data catalogue.\n\
    supported formats: fits, csv, feather\n\
        guess from the file suffix")
parser.add_argument(
    "--bin1_info", type=str, nargs=3, 
    help="Information on the first parameter for reweighting.\n\
            order: (sim_col, data_col, Nbins).")
parser.add_argument(
    "--bin2_info", type=str, nargs=3, 
    help="Information on the second parameter for reweighting.\n\
            order: (sim_col, data_col, Nbins).")
parser.add_argument(
    "--fitting_method", type=str, choices=['tile_based', 'pair_based'], 
    help="Fitting method for the bias calculation: tile_based, pair_based.")
parser.add_argument(
    "--save_surface_prefix", type=str, default=None,
    help="Place to save the 2D shear bias map.")
parser.add_argument(
    "--save_bounds_prefix", type=str, default=None,
    help="Place to save the 2D map bounds.")

## arg parser
args = parser.parse_args()
bias_type = args.bias_type
in_file = args.in_file
out_path = args.out_path

col_e1, col_e2 = args.cols_e12

col_g1, col_g2 = args.cols_g12

col_weight_sim, col_weight_data = args.col_weight_sim_data

col_label = args.col_label
col_id_input = args.col_id_input

if args.col_binning_sim_data is not None:
    print('binning on:', args.col_binning_sim_data)

    col_binning_sim, col_binning_data = args.col_binning_sim_data

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

# for reweighting
if bias_type == 'm_LS_DataRewei_2D':
    in_file_data = args.in_file_data
    bin1_col_sim, bin1_col_data, bin1_Nbins = args.bin1_info
    bin2_col_sim, bin2_col_data, bin2_Nbins = args.bin2_info
    bin1_Nbins = int(bin1_Nbins)
    bin2_Nbins = int(bin2_Nbins)
    fitting_method = args.fitting_method
    save_surface_prefix = args.save_surface_prefix
    save_bounds_prefix = args.save_bounds_prefix
del args

# +++++++++++++++++++++++++++++ workhorse

# load working function
if bias_type == 'm_LS_DataRewei_2D':
    CalFunc = mCalFunc_DataRewei_2D
    print('using data reweighting in 2 parameters')

# check info
if fitting_method == 'tile_based':
    if col_label is not None:
        print('     with the tile based mCalFunc')
    else:
        raise Exception('tile_based chosen, but no col_label provided!')
elif fitting_method == 'pair_based':
    if col_id_input is not None:
        print('     with the pair based mCalFunc')
    else:
        raise Exception('pair_based chosen, but no col_id_input provided!')

# >>>>>>>>>>>>>>>>>> the data
# load data catalogue
file_type = in_file_data[-3:]
if file_type == 'csv':
    cata_data = pd.read_csv(in_file_data)
elif file_type == 'her':
    cata_data = pd.read_feather(in_file_data)
elif file_type == 'its':
    with fits.open(in_file_data) as hdul:
        cata_data = hdul[1].data
else:
    raise Exception(f'Not supported input file type! {in_file_data}')
print('Number of sources in the data', len(cata_data))
### select those within the edges
if ('bin_edges' in locals()):
    cata_data = cata_data[(cata_data[col_binning_data]>bin_edges[0])&(cata_data[col_binning_data]<=bin_edges[-1])]
    print('     within the edges', len(cata_data))
### select non-zero weights
cata_data = cata_data[cata_data[col_weight_data]>0]
print('     with weight > 0', len(cata_data))

# save useable parameters
if ('col_binning_data' in locals()):
    cata_data = pd.DataFrame({
                    'col_binning': np.array(cata_data[col_binning_data]).astype(float),
                    'bin1_col': np.array(cata_data[bin1_col_data]).astype(float),
                    'bin2_col': np.array(cata_data[bin2_col_data]).astype(float),
                    'shape_weight': np.array(cata_data[col_weight_data]).astype(float)
                    })
else:
    cata_data = pd.DataFrame({
                    'bin1_col': np.array(cata_data[bin1_col_data]).astype(float),
                    'bin2_col': np.array(cata_data[bin2_col_data]).astype(float),
                    'shape_weight': np.array(cata_data[col_weight_data]).astype(float)
                    })

# >>>>>>>>>>>>>>>>>> the simulations
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
print('Number of sources in the simulation', len(cata_sim))
### select those within the edges
if ('bin_edges' in locals()):
    cata_sim = cata_sim[(cata_sim[col_binning_sim]>bin_edges[0])&(cata_sim[col_binning_sim]<=bin_edges[-1])]
    print('     within the edges', len(cata_sim))
### select non-zero weights
cata_sim = cata_sim[cata_sim[col_weight_sim]>0]
print('     with weight > 0', len(cata_sim))

# save useable parameters
cata_used = pd.DataFrame({
        'bin1_col': np.array(cata_sim[bin1_col_sim]).astype(float),
        'bin2_col': np.array(cata_sim[bin2_col_sim]).astype(float),
        'e1_out': np.array(cata_sim[col_e1]).astype(float),
        'e2_out': np.array(cata_sim[col_e2]).astype(float),
        'g1_in': np.array(cata_sim[col_g1]).astype(float),
        'g2_in': np.array(cata_sim[col_g2]).astype(float),
        'shape_weight': np.array(cata_sim[col_weight_sim]).astype(float)})
## binning info
if ('col_binning_sim' in locals()):
    cata_used.loc[:, 'col_binning'] = np.array(cata_sim[col_binning_sim]).astype(float)
## IDs
if col_id_input is not None:
    cata_used.loc[:, 'id_input'] = np.array(cata_sim[col_id_input]).astype(int)
## labels
if col_label is not None:
    cata_used.loc[:, 'tile_label'] = np.array(cata_sim[col_label]).astype(str)
## delete original catalogue
del cata_sim

# >>>>>>>>> output
f = open(out_path, 'w')

# >>>>>>>>> binning
if ('col_binning_sim' in locals()):
    ## bin on quantile
    if (not 'bin_edges' in locals()) and (binning_type == 'quantile'):
        quartiles = np.linspace(0, 1, int(N_bin+1))
        bin_edges = stats.mstats.mquantiles(cata_used['col_binning'].values, quartiles)
        ## make sure edge values are included
        bin_edges[0] = min_bin - 0.001
        bin_edges[-1] = max_bin + 0.001
        print('quantile bins', bin_edges)

        ## the whole sample only contains objects within the bin edge
        #### simulation
        cata_used = cata_used[(cata_used['col_binning']>bin_edges[0]) & (cata_used['col_binning']<=bin_edges[-1])]
        cata_used.reset_index(drop=True, inplace=True)
        print('simulation within bin edge:', len(cata_used))
        #### data
        cata_data = cata_data[(cata_data['col_binning']>bin_edges[0]) & (cata_data['col_binning']<=bin_edges[-1])]
        cata_data.reset_index(drop=True, inplace=True)
        print('data within bin edge:', len(cata_data))

    min_bin = bin_edges[0]
    max_bin = bin_edges[-1]

    # mean and median of the x axis
    mean_bin = np.average(cata_used['col_binning'].values, weights=cata_used['shape_weight'].values)

    ## the whole results
    if save_surface_prefix is not None:
        save_surface_path = save_surface_prefix + '_whole.csv'
    else:
        save_surface_path = None
    if save_bounds_prefix is not None:
        save_bounds_prefix_tmp = save_bounds_prefix + '_whole'
    else:
        save_bounds_prefix_tmp = None
    res = CalFunc(cata_used, cata_data, 
                        bin1_col = 'bin1_col', bin1_Nbins = bin1_Nbins,
                        bin2_col = 'bin2_col', bin2_Nbins = bin2_Nbins, 
                        fitting_method = fitting_method,
                        psf_frame = psf_frame,
                        save_surface_path = save_surface_path,
                        save_bounds_prefix = save_bounds_prefix_tmp)
    N_s = len(cata_used)
    Nwei = np.sum(cata_used['shape_weight'].values)

    # collect columns names and values
    cols = ','.join(list(res.keys()))
    vals = ','.join(["{0:0.4f}".format(val) for val in res.values()])
    cols = cols + f',{col_binning_sim}_min,{col_binning_sim}_max,{col_binning_sim}_mean,Nobj,Nwei'
    vals = vals + f',{min_bin},{max_bin},{mean_bin},{N_s},{Nwei}'    
    print(cols, file=f)
    print(vals, file=f)

    for i_bin in range(len(bin_edges)-1):
        min_bin = bin_edges[i_bin]
        max_bin = bin_edges[i_bin+1]

        cata_selec = cata_used[(cata_used['col_binning']>min_bin) & (cata_used['col_binning']<=max_bin)].copy()
        cata_selec.reset_index(drop=True, inplace=True)

        cata_data_selec = cata_data[(cata_data['col_binning']>min_bin) & (cata_data['col_binning']<=max_bin)].copy()
        cata_data_selec.reset_index(drop=True, inplace=True)

        # mean and median of the x axis
        mean_bin = np.average(cata_selec['col_binning'].values, weights=cata_selec['shape_weight'].values)

        # shear bias
        if save_surface_prefix is not None:
            save_surface_path = save_surface_prefix + f'_bin{i_bin}.csv'
        else:
            save_surface_path = None
        if save_bounds_prefix is not None:
            save_bounds_prefix_tmp = save_bounds_prefix + f'_bin{i_bin}'
        else:
            save_bounds_prefix_tmp = None
        res = CalFunc(cata_selec, cata_data_selec, 
                            bin1_col = 'bin1_col', bin1_Nbins = bin1_Nbins,
                            bin2_col = 'bin2_col', bin2_Nbins = bin2_Nbins, 
                            fitting_method = fitting_method,
                            psf_frame = psf_frame,
                            save_surface_path = save_surface_path,
                            save_bounds_prefix = save_bounds_prefix_tmp)
        N_s = len(cata_selec)
        Nwei = np.sum(cata_selec['shape_weight'].values)
        del cata_selec, cata_data_selec

        # print out values
        vals = ','.join(["{0:0.4f}".format(val) for val in res.values()]) \
        + f',{min_bin},{max_bin},{mean_bin},{N_s},{Nwei}'
        print(vals, file=f)

else:
    #only whole
    if save_surface_prefix is not None:
        save_surface_path = save_surface_prefix + '_whole.csv'
    else:
        save_surface_path = None
    if save_bounds_prefix is not None:
        save_bounds_prefix_tmp = save_bounds_prefix + '_whole'
    else:
        save_bounds_prefix_tmp = None
    res = CalFunc(cata_used, cata_data, 
                        bin1_col = 'bin1_col', bin1_Nbins = bin1_Nbins,
                        bin2_col = 'bin2_col', bin2_Nbins = bin2_Nbins, 
                        fitting_method = fitting_method,
                        psf_frame = psf_frame,
                        save_surface_path = save_surface_path,
                        save_bounds_prefix = save_bounds_prefix_tmp)
    N_s = len(cata_used)
    Nwei = np.sum(cata_used['shape_weight'].values)

    # collect columns names and values
    cols = ','.join(list(res.keys()))
    vals = ','.join(["{0:0.4f}".format(val) for val in res.values()])
    cols = cols + f',bin_min,bin_max,bin_mean,Nobj,Nwei'
    vals = vals + f',-999,-999,-999,{N_s},{Nwei}'    
    print(cols, file=f)
    print(vals, file=f)

# close what is opened
f.close()
print(f'results saved to {out_path}')
