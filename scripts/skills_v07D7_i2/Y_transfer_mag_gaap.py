# -*- coding: utf-8 -*-
# @Author: lshuns
# @Date:   2023-01-26 13:45:12
# @Last Modified by:   lshuns
# @Last Modified time: 2023-02-04 14:40:51

### transfer simulation filters to observation filters
# target bands: u, g, r, i, Z (GAaP photometry)
# reference: https://transformcalc.icrar.org/

import os 
import glob

import pandas as pd 
import numpy as np 

# +++++++++++++++++++++++++++++ general info

# main directory contains all the outputs
main_dir = f'/disks/shear16/ssli/ImSim/output/skills_v07D7'

# the original and new folder
in_folder = 'photometry_dr5_ten'
out_folder = 'photometry_dr5_ten_KiDS'

# sub-dir from different input
# split_ran_tags = ['part0', 'part1', 'part2', 'part3', 'part4', 'part5']
# split_ran_tags = ['part0', 'part2', 'part3', 'part4', 'part5']
split_ran_tags = ['part1']

# shear and rotation info
unique_shear_tags = ['m283m283', 'm283p283', 'p283m283', 'p283p283']

# modified bands and modification info
tag_bands = ['u', 'g', 'r', 'i', 'Z', 'i2']
ref_bands = [['u', 'r'], 
                ['g', 'r'],
                ['r', 'g'],
                ['i', 'r'],
                ['Z', 'r'],
                ['i2', 'r']]
mod_file_list = ['../../change_filters/transfer_TAGu_REFur_SDSS2VST.csv',
                '../../change_filters/transfer_TAGg_REFgr_SDSS2VST.csv',
                '../../change_filters/transfer_TAGr_REFrg_SDSS2VST.csv',
                '../../change_filters/transfer_TAGi_REFir_SDSS2VST.csv',
                '../../change_filters/transfer_TAGZ_REFzr_SDSS2VISTA.csv',
                '../../change_filters/transfer_TAGi_REFir_SDSS2VST.csv']

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
        z_in = pd.read_feather(os.path.join(indir_comb, basefile.replace('.feather', '_combined.feather')))
        z_in = z_in['redshift_input'].values
        ## stars have z = 0 
        z_in[z_in < 0] = 0 

        # loop through all bands and correct
        mag_new_list = []
        for i_tag, tag_band in enumerate(tag_bands):

            # get modification info
            cata_filters = pd.read_csv(mod_file_list[i_tag], comment='#')
            z_model = cata_filters['redshift'].values
            alpha_model = cata_filters['alpha'].values
            beta_model = cata_filters['beta'].values
            betaSingle_model = cata_filters['betaSingle'].values
            del cata_filters

            # used columns
            col_tag = f'MAG_GAAP_{tag_band}'
            col_ref1 = f'MAG_GAAP_{ref_bands[i_tag][0]}'
            col_ref2 = f'MAG_GAAP_{ref_bands[i_tag][1]}'
            print('tag:', col_tag, 'ref:', col_ref1, col_ref2)
            mag_tag = cata[col_tag].values
            mag_ref1 = cata[col_ref1].values
            mag_ref2 = cata[col_ref2].values
            # for new results
            mag_tag_new = np.copy(mag_tag)

            # loop over and correct
            N2ref = 0
            Nsingle = 0
            for i_obj, mag_tag_i in enumerate(mag_tag):

                # if not measured, no need for correction
                if ((mag_tag_i <= -99) | (mag_tag_i >= 99)):
                    continue

                # find the nearest redshift
                id_selected = np.argmin(np.abs(z_model - z_in[i_obj]))

                # 1. correct with 2 reference
                if (mag_ref1[i_obj] > -99) & (mag_ref1[i_obj] < 99)\
                        & (mag_ref2[i_obj] > -99) & (mag_ref2[i_obj] < 99):

                    mag_tag_new[i_obj] = mag_tag_i \
                                        + alpha_model[id_selected]*(mag_ref1[i_obj] - mag_ref2[i_obj]) \
                                        + beta_model[id_selected]

                    N2ref += 1

                # 2. correct with itself
                else:
                    mag_tag_new[i_obj] = mag_tag_i + betaSingle_model[id_selected]

                    Nsingle += 1

            print('modified with 2ref', N2ref, N2ref/Ntot)
            print('modified with itself', Nsingle, Nsingle/Ntot)

            # collect results
            mag_new_list.append(mag_tag_new)
            del z_model, alpha_model, beta_model, betaSingle_model
            del mag_tag, mag_ref1, mag_ref2, mag_tag_new

        del z_in

        # save back to the catalogue
        for i_tag, tag_band in enumerate(tag_bands):
            # used columns
            cata.loc[:, f'MAG_GAAP_{tag_band}'] = mag_new_list[i_tag]
        del mag_new_list

        # save
        cata.to_feather(os.path.join(outdir, basefile))
        del cata

# for tile337.6_-32.1_rot90.feather
# total number of obj 106742
# tag: MAG_GAAP_u ref: MAG_GAAP_u MAG_GAAP_r
# modified with 2ref 65842 0.6168331116149219
# modified with itself 40 0.0003747353431638905
# tag: MAG_GAAP_g ref: MAG_GAAP_g MAG_GAAP_r
# modified with 2ref 102722 0.962339098012029
# modified with itself 59 0.0005527346311667385
# tag: MAG_GAAP_r ref: MAG_GAAP_r MAG_GAAP_g
# modified with 2ref 102722 0.962339098012029
# modified with itself 3640 0.03410091622791404
# tag: MAG_GAAP_i ref: MAG_GAAP_i MAG_GAAP_r
# modified with 2ref 102629 0.9614678383391729
# modified with itself 80 0.000749470686327781
# tag: MAG_GAAP_Z ref: MAG_GAAP_Z MAG_GAAP_r
# modified with 2ref 105893 0.9920462423413464
# modified with itself 50 0.0004684191789548631
# tag: MAG_GAAP_i2 ref: MAG_GAAP_i2 MAG_GAAP_r
# modified with 2ref 83673 0.7838807592138053
# modified with itself 56 0.0005246294804294467
# Elapsed:1:21:18.19,User=5095.748,System=1116.758,CPU=127.3%.