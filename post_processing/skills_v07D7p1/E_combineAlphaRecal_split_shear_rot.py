# -*- coding: utf-8 -*-
# @Author: lshuns
# @Date:   2022-05-16 22:25:08
# @Last Modified by:   lshuns
# @Last Modified time: 2023-03-08 13:18:38

### assign the AlphaRecal results to the whole catalogue and then split it into unique inputs

import os
import re

import numpy as np
import pandas as pd

from astropy.io import fits
from astropy.table import Table

# +++++++++++++++++++++++++++++ general info

# the whole catalogue 
# inpath_whole = '/disks/shear16/ssli/ImSim/output/skills_v07D7/skills_v07D7p1_LF_321_kidsPhotometry_everything_col_flag.feather'
inpath_whole = '/disks/shear16/ssli/ImSim/output/skills_v07D7/skills_v07D7p1_LF_321_sharkPhotometry_everything_col_flag.feather'

# the AlphaRecal results
# inpath_alphaRecal = '/disks/shear16/ssli/ImSim/output/skills_v07D7/skills_v07D7p1_LF_321_kidsPhotometry_shear_noSG_noWeiCut_newCut_A1.feather'
inpath_alphaRecal = '/disks/shear16/ssli/ImSim/output/skills_v07D7/skills_v07D7p1_LF_321_sharkPhotometry_shear_noSG_noWeiCut_newCut_A1.feather'

save_form = 'fits'

# shear and rotation info
unique_shear_tags = ['m283m283', 'm283p283', 'p283m283', 'p283p283']
unique_rots = [0., 90.]

# +++++++++++++++++++++++++++++ workhorse

# load the catalogues
cata_final = pd.read_feather(inpath_whole)
print('number in ori', len(cata_final))
cata_alpha = pd.read_feather(inpath_alphaRecal)[['NUMBER', 'run_tag', 'gal_rot', 'tile_label', 'AlphaRecalC_variance', 'AlphaRecalC_weight']]
print('number in alphaRecal', len(cata_alpha))

# merge to get the alphaRecal info
cata_final = cata_final.merge(cata_alpha, on=['NUMBER', 'run_tag', 'gal_rot', 'tile_label'], how='left')
del cata_alpha
print('number after merge', len(cata_final))
print('number with meaningful AlphaRecal', cata_final['AlphaRecalC_weight'].notna().sum())

# dummy values for nan
cata_final.fillna(-999, inplace=True)

# select shear
Ntot = 0
for run_tag in unique_shear_tags:

    mask_shear = cata_final['run_tag'].values == run_tag

    # select rotations
    for rot in unique_rots:
        mask_rot = cata_final['gal_rot'].values == rot
        cata_selec = cata_final[mask_rot&mask_shear]

        cata_selec.reset_index(drop=True, inplace=True)
        print(f'number with shear {run_tag} rot {rot}:', len(cata_selec))
        Ntot += len(cata_selec)

        # save
        if save_form == 'feather':
            out_path = inpath_whole.replace('.feather', f'_shear_{run_tag}_rot_{rot:.0f}.feather')
            ## check existence
            if os.path.isfile(out_path):
                os.remove(out_path)
            cata_selec.to_feather(out_path)

        elif save_form == 'fits':
            out_path = inpath_whole.replace('.feather', f'_shear_{run_tag}_rot_{rot:.0f}.fits')
            ## check existence
            if os.path.isfile(out_path):
                os.remove(out_path)
            Table.from_pandas(cata_selec).write(out_path, format='fits')

        print('catalogue saved as', out_path)
        del cata_selec

print('number in split', Ntot)

# # ########### kids photometry
# number in ori 89303769
# number in alphaRecal 47869587
# number after merge 89303769
# number with meaningful AlphaRecal 47869587
# number with shear m283m283 rot 0.0: 11161407
# catalogue saved as /disks/shear16/ssli/ImSim/output/skills_v07D7/skills_v07D7p1_LF_321_kidsPhotometry_everything_col_flag_shear_m283m283_rot_0.fits
# number with shear m283m283 rot 90.0: 11160581
# catalogue saved as /disks/shear16/ssli/ImSim/output/skills_v07D7/skills_v07D7p1_LF_321_kidsPhotometry_everything_col_flag_shear_m283m283_rot_90.fits
# number with shear m283p283 rot 0.0: 11165963
# catalogue saved as /disks/shear16/ssli/ImSim/output/skills_v07D7/skills_v07D7p1_LF_321_kidsPhotometry_everything_col_flag_shear_m283p283_rot_0.fits
# number with shear m283p283 rot 90.0: 11163255
# catalogue saved as /disks/shear16/ssli/ImSim/output/skills_v07D7/skills_v07D7p1_LF_321_kidsPhotometry_everything_col_flag_shear_m283p283_rot_90.fits
# number with shear p283m283 rot 0.0: 11164742
# catalogue saved as /disks/shear16/ssli/ImSim/output/skills_v07D7/skills_v07D7p1_LF_321_kidsPhotometry_everything_col_flag_shear_p283m283_rot_0.fits
# number with shear p283m283 rot 90.0: 11163629
# catalogue saved as /disks/shear16/ssli/ImSim/output/skills_v07D7/skills_v07D7p1_LF_321_kidsPhotometry_everything_col_flag_shear_p283m283_rot_90.fits
# number with shear p283p283 rot 0.0: 11162941
# catalogue saved as /disks/shear16/ssli/ImSim/output/skills_v07D7/skills_v07D7p1_LF_321_kidsPhotometry_everything_col_flag_shear_p283p283_rot_0.fits
# number with shear p283p283 rot 90.0: 11161251
# catalogue saved as /disks/shear16/ssli/ImSim/output/skills_v07D7/skills_v07D7p1_LF_321_kidsPhotometry_everything_col_flag_shear_p283p283_rot_90.fits
# number in split 89303769
# Elapsed:2:19:30.38,User=4971.891,System=4544.471,CPU=113.6%.

# # # ########### shark photometry
# number in ori 89303769
# number in alphaRecal 47928768
# number after merge 89303769
# number with meaningful AlphaRecal 47928768
# number with shear m283m283 rot 0.0: 11161407
# catalogue saved as /disks/shear16/ssli/ImSim/output/skills_v07D7/skills_v07D7p1_LF_321_sharkPhotometry_everything_col_flag_shear_m283m283_rot_0.fits
# number with shear m283m283 rot 90.0: 11160581
# catalogue saved as /disks/shear16/ssli/ImSim/output/skills_v07D7/skills_v07D7p1_LF_321_sharkPhotometry_everything_col_flag_shear_m283m283_rot_90.fits
# number with shear m283p283 rot 0.0: 11165963
# catalogue saved as /disks/shear16/ssli/ImSim/output/skills_v07D7/skills_v07D7p1_LF_321_sharkPhotometry_everything_col_flag_shear_m283p283_rot_0.fits
# number with shear m283p283 rot 90.0: 11163255
# catalogue saved as /disks/shear16/ssli/ImSim/output/skills_v07D7/skills_v07D7p1_LF_321_sharkPhotometry_everything_col_flag_shear_m283p283_rot_90.fits
# number with shear p283m283 rot 0.0: 11164742
# catalogue saved as /disks/shear16/ssli/ImSim/output/skills_v07D7/skills_v07D7p1_LF_321_sharkPhotometry_everything_col_flag_shear_p283m283_rot_0.fits
# number with shear p283m283 rot 90.0: 11163629
# catalogue saved as /disks/shear16/ssli/ImSim/output/skills_v07D7/skills_v07D7p1_LF_321_sharkPhotometry_everything_col_flag_shear_p283m283_rot_90.fits
# number with shear p283p283 rot 0.0: 11162941
# catalogue saved as /disks/shear16/ssli/ImSim/output/skills_v07D7/skills_v07D7p1_LF_321_sharkPhotometry_everything_col_flag_shear_p283p283_rot_0.fits
# number with shear p283p283 rot 90.0: 11161251
# catalogue saved as /disks/shear16/ssli/ImSim/output/skills_v07D7/skills_v07D7p1_LF_321_sharkPhotometry_everything_col_flag_shear_p283p283_rot_90.fits
# number in split 89303769
# Elapsed:2:21:06.12,User=5257.129,System=7250.874,CPU=147.7%.
