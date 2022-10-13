# -*- coding: utf-8 -*-
# @Author: lshuns
# @Date:   2022-07-03 17:04:38
# @Last Modified by:   lshuns
# @Last Modified time: 2022-07-03 17:08:30

### transfer simulation filters to observation filters
# target bands: MAG_AUTO
# reference: https://transformcalc.icrar.org/

import os 
import glob

import pandas as pd 
import numpy as np 


# +++++++++++++++++++++++++++++ general info

# main directory contains all the outputs
main_dir = f'/disks/shear16/ssli/ImSim/output/skills_v07D7'

# the original and new folder
in_folder = 'SExtractor'
in_photo_folder = 'photometry'
out_folder = 'SExtractor_KiDS'

# sub-dir from different input
split_ran_tags = ['part0', 'part1', 'part2', 'part3', 'part4', 'part5']

# shear and rotation info
unique_shear_tags = ['m283m283', 'm283p283', 'p283m283', 'p283p283']
unique_rots = [0., 90.]

# reference bands and info
ref_bands = ['r', 'g']
mod_file = './transfer_TAGr_REFrg_SDSS2VST.csv'

# +++++++++++++++++++++++++++++ workhorse

## subdir is split run
subdirs = []
for part_tag in split_ran_tags:
    for shear_tag in unique_shear_tags:
        subdirs.append(os.path.join(main_dir, part_tag, shear_tag))
print('Total number of subdirs:', len(subdirs))

## loop over subdirs
for subdir in subdirs:
    print('Running for', subdir)

    indir_comb = os.path.join(subdir, 'catalogues')

    indir = os.path.join(indir_comb, in_folder)
    outdir = os.path.join(indir_comb, out_folder)
    os.mkdir(outdir)

    indir_photo = os.path.join(indir_comb, in_photo_folder)

    ## find all files
    file_list = glob.glob(os.path.join(indir, '*.feather'))
    print('number of files found', len(file_list))

    ## loop through all files
    for infile in file_list:

        basefile = os.path.basename(infile)
        print('for', basefile)

        cata = pd.read_feather(infile)
        Ntot = len(cata)
        print('total number of obj', Ntot)

        # get redshift info from combined catalogue
        z_in = pd.read_feather(os.path.join(indir_comb, basefile.replace('_bandr_', '_').replace('.feather', '_combined.feather')))
        z_in = z_in['redshift_input'].values
        ## stars have z = 0 
        z_in[z_in < 0] = 0 

        # get reference bands from photometry catalogue
        cata_tmp = pd.read_feather(os.path.join(indir_photo, basefile.replace('_bandr_', '_')))
        # used columns
        mag_ref1 = cata_tmp[f'MAG_GAAP_{ref_bands[0]}'].values
        mag_ref2 = cata_tmp[f'MAG_GAAP_{ref_bands[1]}'].values
        del cata_tmp

        # get modification info
        cata_filters = pd.read_csv(mod_file, comment='#')
        z_model = cata_filters['redshift'].values
        alpha_model = cata_filters['alpha'].values
        beta_model = cata_filters['beta'].values
        betaSingle_model = cata_filters['betaSingle'].values
        del cata_filters

        mag_auto = cata['MAG_AUTO'].values
        mag_auto_new = np.copy(mag_auto)
        # print('>>> mag_auto', np.sort(mag_auto))

        # loop over and correct
        N2ref_mag_auto = 0
        Nsingle_mag_auto = 0
        for i_obj, mag_auto_i in enumerate(mag_auto):

            # find the nearest redshift
            id_selected = np.argmin(np.abs(z_model - z_in[i_obj]))

            # 1. correct with 2 reference
            if (mag_ref1[i_obj] > -99) & (mag_ref1[i_obj] < 99)\
                    & (mag_ref2[i_obj] > -99) & (mag_ref2[i_obj] < 99):
                mag_auto_new[i_obj] = mag_auto_i \
                                    + alpha_model[id_selected]*(mag_ref1[i_obj] - mag_ref2[i_obj]) \
                                    + beta_model[id_selected]
                N2ref_mag_auto += 1

            # 2. correct with itself
            else:
                mag_auto_new[i_obj] = mag_auto_i + betaSingle_model[id_selected]
                Nsingle_mag_auto += 1

        del mag_auto
        del z_model, alpha_model, beta_model, betaSingle_model
        del mag_ref1, mag_ref2
        del z_in
        print('modified with 2ref', N2ref_mag_auto, N2ref_mag_auto/Ntot)
        print('modified with itself', Nsingle_mag_auto, Nsingle_mag_auto/Ntot)

        # save back to the catalogue
        cata.loc[:, 'MAG_AUTO'] = mag_auto_new
        del mag_auto_new

        # save
        cata.to_feather(os.path.join(outdir, basefile))
        del cata

