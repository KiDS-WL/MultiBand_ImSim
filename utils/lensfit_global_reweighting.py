# @Author: lshuns
# @Date:   2021-03-10, 14:30:50
# @Last modified by:   ssli
# @Last modified time: 2021-03-22, 9:57:26

### re-calibrate the lensfit weights
###     calibration is done for the combined catalogue from all shear inputs
###         to make sure the average shear is zero across the calibrated catalogue
###     for sake of memory issues, each pointing (different noise realization) is calibrated separately

import os
import re
import glob
import time
import argparse
import subprocess

import pandas as pd
import numpy as np

# +++++++++++++++++++++++++++++ parser for command-line interfaces
parser = argparse.ArgumentParser(
    description=f"lensfit_global_reweighting.py: globally re-calibrate lensfit weights for ellipticity dependency.",
    formatter_class=argparse.RawTextHelpFormatter)
parser.add_argument(
    "--main_dir", type=str,
### NOTE: different sub-directories have identical galaxies&noise but different shear inputs
    help="the top directory containing all simulation outputs.")
parser.add_argument(
    "--lensfit_dir", type=str,
    help="directory to the lensfit code.")
parser.add_argument(
    "--python2_env", type=str, default='python2',
    help="cmd for python2 running.")
parser.add_argument(
    "--clean_up", action="store_true",
    help="Clean all intermediate outputs.")

## arg parser
args = parser.parse_args()
main_dir = args.main_dir
lensfit_dir = args.lensfit_dir
python2_env = args.python2_env
clean_up = args.clean_up

# +++++++++++++++++++++++++++++ workhorse
start_time = time.time()
subdirs = [x for x in glob.glob(os.path.join(main_dir, '*')) if os.path.isdir(x)]
print('Number of subdirs:', len(subdirs))
print('     NOTE: this should match with the number of input shear realizations.')

# get the tile info
## reweighting is done separately for each tile
tmp_noise_rot_files = glob.glob(os.path.join(subdirs[0], 'catalogues', 'shapes', 'tile*'))
tile_info = list(set([re.search(r'tile(.*)_band', tmp).group(1) for tmp in tmp_noise_rot_files]))

for tile_name in tile_info:

    # ++++++++++++ 1. collect original lensfit outputs
    feather_files_ordered = []
    lensfit_out_dirs_ordered = []
    for subdir in subdirs:
        # feather files for one tile in one shear input but different rotations
        tmp_feather = glob.glob(os.path.join(subdir, 'catalogues', 'shapes', f'tile{tile_name}_*'))
        # filenames for the collected feather files
        tmp_outdir = [os.path.join(subdir, 'catalogues', 'shapes', 'tmp_lensfit', os.path.basename(os.path.splitext(tmp)[0])) for tmp in tmp_feather]

        feather_files_ordered.extend(tmp_feather)
        lensfit_out_dirs_ordered.extend(tmp_outdir)

    print(f'tile: {tile_name}')

    # ++++++++++++ 2. combine original lensfit outputs
    mastercat = []
    for lensfit_out_dir in lensfit_out_dirs_ordered:

        # original lensfit output
        lensfit_asc = os.path.join(lensfit_out_dir, 'output.fits.asc')
        with open(lensfit_asc,'r') as fp:
            cat = fp.readlines()
        ## header info
        if not 'header' in locals():
            header = []
            for line in cat:
                if line[0]=='#':
                    header.append(line)
                else:
                    break
            n_header = len(header)
        ## combine
        mastercat += cat[n_header:]
    # save
    outpath = os.path.join(main_dir, f'tile{tile_name}.fits.asc')
    with open(outpath, 'w') as fp:
        fp.writelines(header+mastercat)
    print(f'combined lensfit outputs saved as {outpath}')
    del header

    # ++++++++++++ 3. correct weights for ellipticity dependency
    weights_recal_path = lensfit_dir + '/utils/apply2dtheta_global.py'
    cmd = [python2_env, weights_recal_path, f'--input={outpath}']
    proc = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, cwd=main_dir)
    print(f'Lensfit reweighting results saved as {outpath}.scheme2b_corr')

    # ++++++++++++ 3. save back to feather files
    new_weights = np.loadtxt(f'{outpath}.scheme2b_corr')[:, 4]
    N_lines_record = 0
    for feather_file in feather_files_ordered:

        data_tmp = pd.read_feather(feather_file)

        # how many rows in this file
        N_lines = len(data_tmp)

        # select and save
        data_tmp.loc[:, 'weight_global_LF'] = new_weights[N_lines_record:(N_lines_record+N_lines)].astype(float)
        data_tmp.to_feather(feather_file)
        print(f'Lensfit weights corrected for {feather_file}')

        # iterating
        N_lines_record += N_lines

    # ++++++++++++ 4. save back to combined files (if combined catalogues exist)
    for subdir in subdirs:
        tmp_combined_files = glob.glob(os.path.join(subdir, 'catalogues', f'tile{tile_name}_*_combined.feather'))
        if tmp_combined_files:
            for tmp_combined_file in tmp_combined_files:
                # get rotation
                rot = re.search(r'_rot(\d+)', os.path.basename(tmp_combined_file)).group(1)

                # get data
                data_final = pd.read_feather(tmp_combined_file)

                # corresponding shape file
                tmp_shape = glob.glob(os.path.join(subdir, 'catalogues', 'shapes', f'tile{tile_name}_*_rot{rot}.feather'))[0]
                band = re.search(r'_band(.*)_rot', os.path.basename(tmp_shape)).group(1)
                tmp_info = pd.read_feather(tmp_shape)[['id_detec', 'weight_global_LF']]
                tmp_info = tmp_info.add_suffix(f'_{band}')
                data_final = data_final.merge(tmp_info, left_on='NUMBER', right_on=f'id_detec_{band}', how='left')
                data_final.drop(columns=[f'id_detec_{band}'], inplace=True)

                # save
                data_final.to_feather(tmp_combined_file)
                print(f'Lensfit weight_global_LF saved back to {tmp_combined_file}')

    # ++++++++++++ -1. clean up
    if clean_up:
        print('Clean up tmp files')
        for lensfit_out_dir in lensfit_out_dirs_ordered:
            shutil.rmtree(lensfit_out_dir)
        for tmp_file in glob.glob(f'{outpath}*'):
            os.remove(tmp_file)
print('All finished in', time.time()-start_time, 's')
