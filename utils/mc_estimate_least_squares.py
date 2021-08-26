# -*- coding: utf-8 -*-
# @Author: lshuns
# @Date:   2021-08-25 12:37:31
# @Last Modified by:   lshuns
# @Last Modified time: 2021-08-25 12:45:54

### a script to calculate the mc shear bias directly from simulated catalogues
###    using new algorithm based on least squares
###### dependence: mc_CalFunc.mcCalFunc_least_squares

########################
# usage: mc_estimate_least_squares.py [-h] [--in_file_sim IN_FILE_SIM]
#                                     [--out_path OUT_PATH]
#                                     [--cols_g12_e12 COLS_G12_E12 COLS_G12_E12 COLS_G12_E12 COLS_G12_E12]
#                                     [--e_type {input,measured}]
#                                     [--col_binning COL_BINNING]
#                                     [--binning_edges [BINNING_EDGES [BINNING_EDGES ...]]]
#                                     [--binning_min_max_Nbin BINNING_MIN_MAX_NBIN BINNING_MIN_MAX_NBIN BINNING_MIN_MAX_NBIN]
#                                     [--log_binning] [--col_weight COL_WEIGHT]
#                                     [--show_plot] [--nboot NBOOT]
# mc_estimate_simple.py: estimate mc shear bias directly from simulated catalogues.
# optional arguments:
#   -h, --help            show this help message and exit
#   --in_file_sim IN_FILE_SIM
#                         the input simulation catalogue.
#                             supported formats: fits, csv, feather
#                                 guess from the file suffix
#   --out_path OUT_PATH   where to save the final results. 
#                             Default ./test.txt.
#   --cols_g12_e12 COLS_G12_E12 COLS_G12_E12 COLS_G12_E12 COLS_G12_E12
#                         column names for g1_in, g2_in, e1, e2.
#   --e_type {input,measured}
#                         what is the e values: input OR measured.
#   --col_binning COL_BINNING
#                         columns to be used for binning.
#                             If not provided, no binning will be applied.
#   --binning_edges [BINNING_EDGES [BINNING_EDGES ...]]
#                         edges for binning, 
#                                     should not set along with binning_min_max_Nbin.
#   --binning_min_max_Nbin BINNING_MIN_MAX_NBIN BINNING_MIN_MAX_NBIN BINNING_MIN_MAX_NBIN
#                         min max N_bins for binning.
#                                     should not set along with binning_edges.
#   --log_binning         binning on log scale, in which case binning_min_max_Nbin are treated in log scale.
#   --col_weight COL_WEIGHT
#                         columns to the weight.
#                             If not provided, equal weight for all galaxies.
#   --show_plot           Show an example plot.
#   --nboot NBOOT         Number of boots used for bootstrapping.
########################

import argparse

import numpy as np
import pandas as pd

from astropy.io import fits
from astropy.table import Table

# plot related
import matplotlib as mpl
mpl.rcParams['xtick.direction'] = 'in'
mpl.rcParams['ytick.direction'] = 'in'
mpl.rcParams['xtick.top'] = True
mpl.rcParams['ytick.right'] = True
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator

# Self-defined package
from mc_CalFunc import mcCalFunc_least_squares

# +++++++++++++++++++++++++++++ parser for command-line interfaces
parser = argparse.ArgumentParser(
    description=f"mc_estimate_simple.py: estimate mc shear bias directly from simulated catalogues.",
    formatter_class=argparse.RawTextHelpFormatter)
parser.add_argument(
    "--in_file_sim", type=str,
    help="the input simulation catalogue.\n\
    supported formats: fits, csv, feather\n\
        guess from the file suffix")
parser.add_argument(
    "--out_path", type=str, default='./test.txt',
    help="where to save the final results. \n\
    Default ./test.txt.")
parser.add_argument(
    "--cols_g12_e12", type=str, nargs=4, default=['g1_in', 'g2_in', 'e1', 'e2'], 
    help="column names for g1_in, g2_in, e1, e2.")
parser.add_argument(
    "--e_type", type=str, choices=['input', 'measured'], 
    help="what is the e values: input OR measured.")
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
    "--log_binning", action="store_true",
    help="binning on log scale, in which case binning_min_max_Nbin are treated in log scale.")
parser.add_argument(
    "--col_weight", type=str,
    help="columns to the weight.\n\
    If not provided, equal weight for all galaxies.")
parser.add_argument(
    "--show_plot", action="store_true",
    help="Show an example plot.")
parser.add_argument(
    "--nboot", type=int, default=500,
    help="Number of boots used for bootstrapping.")

## arg parser
args = parser.parse_args()
in_file_sim = args.in_file_sim
out_path = args.out_path
col_g1, col_g2, col_e1, col_e2 = args.cols_g12_e12
e_type = args.e_type

col_binning = args.col_binning
log_binning = args.log_binning
if col_binning is not None:
    if (args.binning_edges is not None) and (args.binning_min_max_Nbin is not None):
        raise Exception('cannot assign both binning_edges and binning_min_max_Nbin!')
    if args.binning_edges is not None:
        bin_edges = np.array(args.binning_edges)
    elif args.binning_min_max_Nbin is not None:
        min_bin, max_bin, N_bin = args.binning_min_max_Nbin
        if log_binning:
            bin_edges = np.logspace(min_bin, max_bin, int(N_bin))
        else:
            bin_edges = np.linspace(min_bin, max_bin, int(N_bin))
    else:
        raise Exception('either binning_edges or binning_min_max_Nbin should be assigned for col_binning==True!')

    print('binning on:', col_binning)
    print('edges:', bin_edges)

col_weight = args.col_weight
show_plot = args.show_plot
nboot = args.nboot

# +++++++++++++++++++++++++++++ workhorse
# load simulation catalogue
file_type = in_file_sim[-3:]
if file_type == 'csv':
    cata_sim = pd.read_csv(in_file_sim)
elif file_type == 'her':
    cata_sim = pd.read_feather(in_file_sim)
elif file_type == 'its':
    with fits.open(in_file_sim) as hdul:
        cata_sim = Table(hdul[1].data).to_pandas()
else:
    raise Exception(f'Not supported input file type! {in_file_sim}')
print('Number of sources in simulated catalogue', len(cata_sim))

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
    'g1_in': np.array(cata_sim[col_g1]),
    'g2_in': np.array(cata_sim[col_g2]),
    'e1_out': e1_out,
    'e2_out': e2_out})
## binning values
if col_binning is not None:
    cata_used.loc[:, col_binning] = np.array(cata_sim[col_binning])
## weight for galaxies
if col_weight is not None:
    cata_used.loc[:, 'shape_weight'] = np.array(cata_sim[col_weight])
else:
    cata_used.loc[:, 'shape_weight'] = 1
## delete original catalogue
del cata_sim

# output
f = open(out_path, 'w')

# whole sample
res = mcCalFunc_least_squares(cata_used, nboot=nboot)
N_s = len(cata_used)

# collect columns names and values
cols = '# ' + '    '.join(list(res.keys()))
vals = '    '.join(["{0:0.4f}".format(val) for val in res.values()])
if col_binning is not None:
    cols = cols + f'    {col_binning}_min    {col_binning}_max    Number'
    vals = vals + f'    -999      -999      {N_s}'    
print(cols, file=f)
print(vals, file=f)

# binning
if col_binning is not None:
    for i_bin in range(len(bin_edges)-1):
        min_bin = bin_edges[i_bin]
        max_bin = bin_edges[i_bin+1]
        mask_bin = (cata_used[col_binning]>=min_bin) & (cata_used[col_binning]<max_bin)

        cata_selec = cata_used[mask_bin]
        cata_selec.reset_index(drop=True, inplace=True)

        res = mcCalFunc_least_squares(cata_selec, nboot=nboot)
        N_s = len(cata_selec)

        # print out values
        vals = '    '.join(["{0:0.4f}".format(val) for val in res.values()]) + f'    {min_bin}    {max_bin}      {N_s}'
        print(vals, file=f)

# close what is opened
print(f'results saved to {out_path}')
f.close()

# make example plot if required
if show_plot:

    data = np.loadtxt(out_path)
    ### first point is for the whole sample
    data = data[1:]

    # x position
    x_val = (data[:, -2]+data[:, -3])/2.

    ####### 0. source number
    # y position
    y_val = data[:, -1]
    # plotting
    fig, ax = plt.subplots()
    ax.errorbar(x_val, y_val, color='r', label=None, linestyle='', marker='x', markersize=3, capsize=3)
    ax.xaxis.set_minor_locator(AutoMinorLocator())
    ax.yaxis.set_minor_locator(AutoMinorLocator())
    plt.xlabel(col_binning)
    plt.ylabel('source number')
    if log_binning:
        plt.xscale('log')
    plt.show()
    plt.close()

    ####### 1. m bias
    # y position
    y_val = (data[:, 0] + data[:, 1])/2.
    # error
    yerr = 0.5*(data[:, 4] + data[:, 5])
    # plotting
    fig, ax = plt.subplots()
    ax.errorbar(x_val, y_val, yerr=np.vstack([yerr, yerr]), color='r', label=None, linestyle='', marker='x', markersize=3, capsize=3)
    ax.xaxis.set_minor_locator(AutoMinorLocator())
    ax.yaxis.set_minor_locator(AutoMinorLocator())
    plt.xlabel(col_binning)
    plt.ylabel('m')
    if log_binning:
        plt.xscale('log')
    plt.show()
    plt.close()

    ####### 2. c bias
    # y position
    y_val = (data[:, 2] + data[:, 3])/2.
    # error
    yerr = 0.5*(data[:, 6] + data[:, 7])
    # plotting
    fig, ax = plt.subplots()
    ax.errorbar(x_val, y_val, yerr=np.vstack([yerr, yerr]), color='b', label=None, linestyle='', marker='x', markersize=3, capsize=3)
    ax.xaxis.set_minor_locator(AutoMinorLocator())
    ax.yaxis.set_minor_locator(AutoMinorLocator())
    plt.xlabel(col_binning)
    plt.ylabel('c')
    if log_binning:
        plt.xscale('log')
    plt.show()
    plt.close()    