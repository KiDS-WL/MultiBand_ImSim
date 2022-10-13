# -*- coding: utf-8 -*-
# @Author: lshuns
# @Date:   2022-07-03 16:49:54
# @Last Modified by:   lshuns
# @Last Modified time: 2022-07-03 17:01:23

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
in_folder = 'photometry'
out_folder = 'photometry_KiDS'

# sub-dir from different input
split_ran_tags = ['part0', 'part1', 'part2', 'part3', 'part4', 'part5']

# shear and rotation info
unique_shear_tags = ['m283m283', 'm283p283', 'p283m283', 'p283p283']
unique_rots = [0., 90.]

# modified bands and modification info
tag_bands = ['u', 'g', 'r', 'i', 'Z']
ref_bands = [['u', 'r'], 
                ['g', 'r'],
                ['r', 'g'],
                ['i', 'r'],
                ['Z', 'r']]
mod_file_list = ['./transfer_TAGu_REFur_SDSS2VST.csv',
                './transfer_TAGg_REFgr_SDSS2VST.csv',
                './transfer_TAGr_REFrg_SDSS2VST.csv',
                './transfer_TAGi_REFir_SDSS2VST.csv',
                './transfer_TAGZ_REFzr_SDSS2VISTA.csv']

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



# for tile137.0_-1.5_rot0.feather
# total number of obj 109479
# tag: MAG_GAAP_u ref: MAG_GAAP_u MAG_GAAP_r
# modified with 2ref 75074 0.6857388174901122
# modified with itself 33 0.00030142767106020334
# tag: MAG_GAAP_g ref: MAG_GAAP_g MAG_GAAP_r
# modified with 2ref 106683 0.9744608555065355
# modified with itself 44 0.0004019035614136044
# tag: MAG_GAAP_r ref: MAG_GAAP_r MAG_GAAP_g
# modified with 2ref 106683 0.9744608555065355
# modified with itself 2434 0.02223257428365257
# tag: MAG_GAAP_i ref: MAG_GAAP_i MAG_GAAP_r
# modified with 2ref 100468 0.9176919774568639
# modified with itself 44 0.0004019035614136044
# tag: MAG_GAAP_Z ref: MAG_GAAP_Z MAG_GAAP_r
# modified with 2ref 106882 0.976278555704747
# modified with itself 17 0.00015528092145525627
# for tile137.0_-1.5_rot90.feather
# total number of obj 109496
# tag: MAG_GAAP_u ref: MAG_GAAP_u MAG_GAAP_r
# modified with 2ref 75865 0.6928563600496822
# modified with itself 42 0.00038357565573171625
# tag: MAG_GAAP_g ref: MAG_GAAP_g MAG_GAAP_r
# modified with 2ref 106749 0.9749123255644042
# modified with itself 56 0.0005114342076422883
# tag: MAG_GAAP_r ref: MAG_GAAP_r MAG_GAAP_g
# modified with 2ref 106749 0.9749123255644042
# modified with itself 2397 0.021891210637831518
# tag: MAG_GAAP_i ref: MAG_GAAP_i MAG_GAAP_r
# modified with 2ref 100601 0.9187641557682472
# modified with itself 64 0.0005844962373054724
# tag: MAG_GAAP_Z ref: MAG_GAAP_Z MAG_GAAP_r
# modified with 2ref 106860 0.9759260612259809
# modified with itself 37 0.0003379118871922262
# for tile137.0_1.5_rot0.feather
# total number of obj 120997
# tag: MAG_GAAP_u ref: MAG_GAAP_u MAG_GAAP_r
# modified with 2ref 79632 0.6581320198021439
# modified with itself 39 0.0003223220410423399
# tag: MAG_GAAP_g ref: MAG_GAAP_g MAG_GAAP_r
# modified with 2ref 113520 0.9382050794647802
# modified with itself 55 0.0004545567245468896
# tag: MAG_GAAP_r ref: MAG_GAAP_r MAG_GAAP_g
# modified with 2ref 113520 0.9382050794647802
# modified with itself 7097 0.05865434680198683
# tag: MAG_GAAP_i ref: MAG_GAAP_i MAG_GAAP_r
# modified with 2ref 107859 0.8914187955073266
# modified with itself 64 0.0005289387340181988
# tag: MAG_GAAP_Z ref: MAG_GAAP_Z MAG_GAAP_r
# modified with 2ref 118622 0.9803714141672935
# modified with itself 37 0.0003057927056042712
# for tile137.0_1.5_rot90.feather
# total number of obj 121218
# tag: MAG_GAAP_u ref: MAG_GAAP_u MAG_GAAP_r
# modified with 2ref 81019 0.6683743338448085
# modified with itself 27 0.00022273919714893827
# tag: MAG_GAAP_g ref: MAG_GAAP_g MAG_GAAP_r
# modified with 2ref 113363 0.9351993928294478
# modified with itself 37 0.0003052351960929895
# tag: MAG_GAAP_r ref: MAG_GAAP_r MAG_GAAP_g
# modified with 2ref 113363 0.9351993928294478
# modified with itself 7527 0.062094738405187345
# tag: MAG_GAAP_i ref: MAG_GAAP_i MAG_GAAP_r
# modified with 2ref 108151 0.8922024781798082
# modified with itself 51 0.00042072959461466117
# tag: MAG_GAAP_Z ref: MAG_GAAP_Z MAG_GAAP_r
# modified with 2ref 118967 0.9814301506376941
# modified with itself 18 0.00014849279809929218
# Elapsed:3:24.35,User=202.940,System=63.613,CPU=130.4%.
