# -*- coding: utf-8 -*-
# @Author: lshuns
# @Date:   2021-06-30 14:11:00
# @Last Modified by:   lshuns
# @Last Modified time: 2021-09-09 14:17:35

### create KiDS-phase1-like catalogue from ImSim outputs
###### reference: https://github.com/KiDS-WL/THELI_catalogues/blob/master/scripts/create_catalogue_products_K1000_Phase1.sh
### main procedures:
###     1. non-lensfit selection 
###     2. binning galaxies based on measured PSFs (5*50 K1000 fiducial)
###     3. perform lensfit re-weighting in each PSF bin
###     4. lensfit selection and combine

########################
# usage: phase1_forImSim.py [-h] [--in_dir IN_DIR] [--outfile_path OUTFILE_PATH]
#                           [--tmp_dir TMP_DIR] [--Ncores NCORES]
#                           [--recalibrate_LFweights RECALIBRATE_LFWEIGHTS]
#                           [--wanted_cols WANTED_COLS] [--clean_up]
# phase1_forImSim.py: create KiDS-phase1-like catalogue from ImSim outputs.
# optional arguments:
#   -h, --help            show this help message and exit
#   --in_dir IN_DIR       the top directory containing all ImSim outputs.
#   --outfile_path OUTFILE_PATH
#                         the path for the final catalogue.
#   --tmp_dir TMP_DIR     tmp directory for the intermediate outcomes.
#   --Ncores NCORES       Number of cores used for parallel run.
#   --recalibrate_LFweights RECALIBRATE_LFWEIGHTS
#                         path to the reweighting code.
#   --wanted_cols WANTED_COLS
#                         path to the list of wanted columns.
#                             If not provided, save all columns.
#   --clean_up            Clean all intermediate outputs.
########################

import os
import re
import glob
import time
import shutil
import argparse
import subprocess

import numpy as np
import pandas as pd
import multiprocessing as mp

from astropy.io import fits
from astropy.table import Table

# +++++++++++++++++++++++++++++ parser for command-line interfaces
parser = argparse.ArgumentParser(
    description=f"phase1_forImSim.py: create KiDS-phase1-like catalogue from ImSim outputs.",
    formatter_class=argparse.RawTextHelpFormatter)
parser.add_argument(
    "--in_dir", type=str,
    help="the top directory containing all ImSim outputs.")
parser.add_argument(
    "--outfile_path", type=str,
    help="the path for the final catalogue.")
parser.add_argument(
    "--tmp_dir", type=str,
    help="tmp directory for the intermediate outcomes.")
parser.add_argument(
    "--Ncores", type=int, default=12,
    help="Number of cores used for parallel run.")
parser.add_argument(
    "--recalibrate_LFweights", type=str, default='./recalibrate_LFweights_rotated_to_psf_frame_forImSim.py',
    help="path to the reweighting code.")
parser.add_argument(
    "--wanted_cols", type=str,
    help="path to the list of wanted columns.\n\
    If not provided, save all columns.")
parser.add_argument(
    "--clean_up", action="store_true",
    help="Clean all intermediate outputs.")

## arg parser
args = parser.parse_args()
in_dir = args.in_dir
outfile_path = args.outfile_path
tmp_dir = args.tmp_dir
Ncores = args.Ncores
recalibrate_LFweights = args.recalibrate_LFweights
wanted_cols = args.wanted_cols
if wanted_cols is not None:
    with open(wanted_cols) as f:
        wanted_cols = f.readlines()
    wanted_cols = [x.strip() for x in wanted_cols if x[0]!='#'] 
clean_up = args.clean_up

## check existence
if not os.path.isdir(in_dir):
    raise Exception(f'{in_dir} does not exist!')
if not os.path.isfile(recalibrate_LFweights):
    raise Exception(f'{recalibrate_LFweights} does not exist!')

## clean tmp dir for each run
if os.path.isdir(tmp_dir):
    shutil.rmtree(tmp_dir)
os.mkdir(tmp_dir)
## where to save the binning catalogues
split_dir = os.path.join(tmp_dir, 'split_cata')
os.mkdir(split_dir)
## where to save the re-weighted catalogues
reweighting_dir = os.path.join(tmp_dir, 'split_recal')
os.mkdir(reweighting_dir)

## the fiducial binning in K1000
#### 5*50
psfsize_edges=(1.070443, 1.698965, 1.881802, 2.080181, 2.292835, 4.076628)
psfebin_edges=(0., 0.00222, 0.003158, 0.003895, 0.004528, 0.005104, 0.005644, 0.006161, 0.006649, 0.007131, \
    0.007603, 0.008053, 0.008507, 0.00895, 0.009396, 0.009833, 0.010277, 0.01072, 0.011168, 0.011625, \
    0.012088, 0.012554, 0.013027, 0.013515, 0.014013, 0.014534, 0.01506, 0.015598, 0.016147, 0.01672, \
    0.01731, 0.017928, 0.018572, 0.01925, 0.019965, 0.020718, 0.021513, 0.022366, 0.023311, 0.024331, \
    0.025454, 0.026701, 0.028086, 0.029606, 0.031364, 0.033408, 0.035847, 0.03904, 0.043501, 0.051788, 0.194039)

# +++++++++++++++++++++++++++++ workhorse
start_time0 = time.time()

# get subdirs
subdirs = [x for x in glob.glob(os.path.join(in_dir, '*')) if (os.path.isdir(x) and (not 'tmp' in os.path.basename(x)))]
print('Number of subdirs:', len(subdirs))

# tile based manipulation
print('>>> start tile-based manipulation (non-lensfit selection + PSF binning)...')
start_time = time.time()
for subdir in subdirs:

    run_tag = os.path.basename(subdir)

    # basic info for input shear values
    with open(os.path.join(subdir, 'basic_info.txt'), 'r') as opened_file:
        all_lines = opened_file.read().splitlines()
    useful_line = [line for line in all_lines if 'g_cosmic' in line][0]
    g1g2 = useful_line.split('=')[1]
    ## variable shears
    if 'variable' in g1g2:
        g_const = False
    ## constant shear
    else:
        g_const = True
        g1 = float(g1g2.split()[0])
        g2 = float(g1g2.split()[1])

    # main catalogues from simulation
    file_list = glob.glob(os.path.join(subdir, 'catalogues', '*_combined.*'))
    if not file_list:
        raise Exception(f'Cannot find any combined catalogues in tag {run_tag}!\n\
        make sure taskID=7 is performed in the main pipeline!')
    print(f'Number of files in {run_tag}: {len(file_list)}' )
    for file in file_list:
        file_type = file[-3:]
        if file_type == 'csv':
            tile_cata = pd.read_csv(file)
        elif file_type == 'her':
            tile_cata = pd.read_feather(file)
        elif file_type == 'its':
            with fits.open(file) as hdul:
                tile_cata = Table(hdul[1].data).to_pandas()
        else:
            raise Exception(f'Not supported input file type! {file}')

        # get the shear in case of variable shears
        if not g_const:
            g1 = np.array(tile_cata['gamma1_input'])
            g2 = np.array(tile_cata['gamma2_input'])

        # select wanted columns
        if wanted_cols is not None:
            tile_cata = tile_cata[wanted_cols]

        # some columns for tile info
        ## input cosmic shear
        tile_cata.loc[:, 'g1_in'] = g1
        tile_cata.loc[:, 'g2_in'] = g2
        ## tile noise info
        tile_label = re.search(r'tile(.*)_rot', file).group(1)
        tile_cata.loc[:, 'tile_label'] = tile_label
        ## rotation info
        gal_rot = float(re.search(r'_rot(\d+)', file).group(1))
        tile_cata.loc[:, 'gal_rot'] = gal_rot

        # some columns later for binning the galaxies
        tile_cata.loc[:, 'PSFsize'] = np.array((tile_cata['psf_Q11_LF_r']*tile_cata['psf_Q22_LF_r'] - tile_cata['psf_Q12_LF_r']**2.)**0.5)
        tile_cata.loc[:, 'PSFeabs'] = np.hypot(tile_cata['psf_e1_LF_r'].values, tile_cata['psf_e2_LF_r'].values)

        # ======= 1. non-lensfit selection
        #### a) input match selection
        mask_input = (tile_cata['id_input']>-999)

        #### b) remove unmeasured
        mask_psf = (tile_cata['psf_Q11_LF_r'] != 0.0) & (tile_cata['psf_Q22_LF_r'] != 0.0)

        #### c) remove asteroids
        try:
            gmr = np.array(tile_cata['MAG_GAAP_0p7_g']) - np.array(tile_cata['MAG_GAAP_0p7_r'])
            imr = np.array(tile_cata['MAG_GAAP_0p7_i']) - np.array(tile_cata['MAG_GAAP_0p7_r'])
            mask_ast = (gmr <= 1.5) | (imr <= 1.5)
        except KeyError:
            mask_ast = np.full(len(tile_cata), True)

        #### d) remove binaries
        try:
            mask_binary = (np.hypot(tile_cata['e1_LF_r'].values, tile_cata['e2_LF_r'].values) <= 0.8) | (tile_cata['scalelength_LF_r'] >= (0.5 * np.exp(0.65788*(24.2 - tile_cata['MAG_GAAP_0p7_r']))))
        except KeyError:
            mask_binary = (np.hypot(tile_cata['e1_LF_r'].values, tile_cata['e2_LF_r'].values) <= 0.8) | (tile_cata['scalelength_LF_r'] >= (0.5 * np.exp(0.65788*(24.2 - tile_cata['MAG_AUTO']))))

        # apply the mask
        tile_cata = tile_cata[mask_input & mask_psf & mask_ast & mask_binary]
        tile_cata.reset_index(drop=True, inplace=True)

        # ======= 2. binning galaxies 
        for i_psfsize in range(len(psfsize_edges)-1):
            psfsize_min = psfsize_edges[i_psfsize]
            psfsize_max = psfsize_edges[i_psfsize + 1]

            for i_psfe in range(len(psfebin_edges)-1):
                psfebin_min = psfebin_edges[i_psfe]
                psfebin_max = psfebin_edges[i_psfe + 1]

                # select galaxies
                mask_tmp = (tile_cata['PSFsize'] >= psfsize_min) & (tile_cata['PSFsize'] < psfsize_max) & (tile_cata['PSFeabs'] >= psfebin_min) & (tile_cata['PSFeabs'] < psfebin_max)
                if np.sum(mask_tmp)>0:
                    tile_cata_selec = tile_cata[mask_tmp]
                    tile_cata_selec.reset_index(drop=True, inplace=True)
                    outpath_tmp = os.path.join(split_dir, f'{run_tag}_tile{tile_label}_rot{gal_rot:.0f}_{i_psfsize}_{i_psfe}.feather')
                    tile_cata_selec.to_feather(outpath_tmp)
                    # print('binning galaxies saved to', outpath_tmp)
print('<<< tile-based manipulation finished in', (time.time()-start_time)/60., 'min')

# ======= 3. re-weighting for each bin
print('>>> start re-weighting for each bin...')
start_time = time.time()

# function for parallel run
def FuncRewei(i_psfsize, i_psfe, recalibrate_LFweights, file4reweighting, reweighting_dir):
    ## log for reweighting run
    log_file = os.path.join(reweighting_dir, f'log_{i_psfsize}_{i_psfe}.txt')
    outLog = open(log_file, "w")

    # reweighting 
    cmd = ['python', recalibrate_LFweights, file4reweighting, file4reweighting]
    proc = subprocess.run(cmd, stdout=outLog, stderr=outLog)
    outLog.close()    

# actual run
print(f'Number of processes for re-weighting: {Ncores}')
work_pool = mp.Pool(processes=Ncores)
proc_list = []
for i_psfsize in range(len(psfsize_edges)-1):
    for i_psfe in range(len(psfebin_edges)-1):

        # combine all the catalogues belong to the same bin
        file_list = glob.glob(os.path.join(split_dir, f'*_{i_psfsize}_{i_psfe}.feather'))
        n_files = len(file_list)
        if n_files == 0:
            continue
        print(f'Number of files in bin {i_psfsize}, {i_psfe}:', n_files)
        cata_bin = []
        for file in file_list:
            tile_cata = pd.read_feather(file)
            cata_bin.append(tile_cata)
        cata_bin = pd.concat(cata_bin)
        cata_bin.reset_index(drop=True, inplace=True)
        ## save for reweighting
        file4reweighting = os.path.join(reweighting_dir, f'combined_{i_psfsize}_{i_psfe}.feather')
        cata_bin.to_feather(file4reweighting)

        ### running
        proc = work_pool.apply_async(func=FuncRewei,
                                args=(i_psfsize, i_psfe, recalibrate_LFweights, file4reweighting, reweighting_dir))
        proc_list.append(proc)

work_pool.close()
work_pool.join()
### check for any errors during run
for proc in proc_list:
    proc.get()
print('<<< re-weighting finished in', (time.time()-start_time)/60., 'min')

# ======= 4. perform lensfit selection & combine all catalogues
print('>>> start lensfit selection for each bin...')
start_time = time.time()
cata_final = []
for i_psfsize in range(len(psfsize_edges)-1):
    for i_psfe in range(len(psfebin_edges)-1):
        file4reweighting = os.path.join(reweighting_dir, f'combined_{i_psfsize}_{i_psfe}.feather')
        if not os.path.isfile(file4reweighting):
            continue
        cata_bin = pd.read_feather(file4reweighting)

        # 0) check if success of reweighting 
        if not 'recal_weight_LF_r' in cata_bin.columns:
            print(f'WARNING: recalibrate_LFweights unsuccess for {file4reweighting}!')
            print(f'            Number of source: {len(cata_bin)}')
            os.rename(file4reweighting, file4reweighting.replace('.feather', '_failed.feather'))
            continue

        # a) fitclass cut
        mask_class = (cata_bin['class_LF_r']!=-1) & (cata_bin['class_LF_r']!=-10) & (cata_bin['class_LF_r']!=-4) & (cata_bin['class_LF_r']!=1) & \
                        (cata_bin['class_LF_r']!=2) & (cata_bin['class_LF_r']!=-7) & (cata_bin['class_LF_r']!=-3)

        # b) magnitude cut
        mask_mag = (cata_bin['MAG_AUTO']>20.0)

        # c) blending cut
        mask_blending = (cata_bin['contamination_radius_LF_r']>4.25)

        # d) weighting cut
        mask_weight = (cata_bin['recal_weight_LF_r']>0)

        # e) 9-band photometry cut
        cata_bin.loc[:, 'FLAG_GAAP_0p7_ugriZYJHKs'] = 0
        try:
            for band in ['u', 'g', 'r', 'i', 'Z', 'Y', 'J', 'H', 'Ks']:
                cata_bin.loc[:, 'FLAG_GAAP_0p7_ugriZYJHKs'] += np.array(cata_bin[f'FLAG_GAAP_0p7_{band}'])
        except KeyError:
            pass
        mask_9band = cata_bin['FLAG_GAAP_0p7_ugriZYJHKs'] == 0

        ## apply and save
        cata_bin = cata_bin[mask_class & mask_mag & mask_blending & mask_weight & mask_9band]
        cata_final.append(cata_bin)
print('<<< lensfit selection finished in', (time.time()-start_time)/60., 'min')

cata_final = pd.concat(cata_final)
cata_final.reset_index(drop=True, inplace=True)
print(f'Total number of source {len(cata_final)}')

## check existence
if os.path.isfile(outfile_path):
    os.remove(outfile_path)
## save
if outfile_path[-3:] == 'her':
    cata_final.to_feather(outfile_path)
elif outfile_path[-3:] == 'csv':
    cata_final.to_csv(outfile_path, index=False)
elif outfile_path[-3:] == 'its':
    Table.from_pandas(cata_final).write(outfile_path, format='fits')
else:
    raise Exception(f'Not supported output file type! {outfile_path}')
print('final catalogue saved as', outfile_path)

# clean intermedate outcomes if set
if clean_up:
    shutil.rmtree(tmp_dir)
    print('All the intermediate outputs have been cleaned.')

print('all finished in', (time.time()-start_time0)/60., 'min')