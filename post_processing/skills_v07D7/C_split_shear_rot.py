# -*- coding: utf-8 -*-
# @Author: lshuns
# @Date:   2022-05-16 22:25:08
# @Last Modified by:   lshuns
# @Last Modified time: 2022-10-11 09:31:57

### split the catalogue to unique inputs

import os
import re

import numpy as np
import pandas as pd

from astropy.io import fits
from astropy.table import Table

# +++++++++++++++++++++++++++++ general info

# input 
inpath = '/disks/shear16/ssli/ImSim/output/skills_v07D7/skills_v07D7_LF_321_kidsPhotometry_everything_col_flag.feather'

save_form = 'fits'

# shear and rotation info
unique_shear_tags = ['m283m283', 'm283p283', 'p283m283', 'p283p283']
unique_rots = [0., 90.]

# +++++++++++++++++++++++++++++ workhorse

# load the catalogue
cata_final = pd.read_feather(inpath)
print('number in ori', len(cata_final))

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
            out_path = inpath.replace('.feather', f'_shear_{run_tag}_rot_{rot:.0f}.feather')
            ## check existence
            if os.path.isfile(out_path):
                os.remove(out_path)
            cata_selec.to_feather(out_path)

        elif save_form == 'fits':
            out_path = inpath.replace('.feather', f'_shear_{run_tag}_rot_{rot:.0f}.fits')
            ## check existence
            if os.path.isfile(out_path):
                os.remove(out_path)
            Table.from_pandas(cata_selec).write(out_path, format='fits')

        print('catalogue saved as', out_path)
        del cata_selec

print('number in split', Ntot)

########### kids photometry
# number in ori 89303975
# number with shear m283m283 rot 0.0: 11161396
# catalogue saved as /disks/shear16/ssli/ImSim/output/skills_v07D7/skills_v07D7_LF_321_kidsPhotometry_everything_col_flag_shear_m283m283_rot_0.fits
# number with shear m283m283 rot 90.0: 11160604
# Xcatalogue saved as /disks/shear16/ssli/ImSim/output/skills_v07D7/skills_v07D7_LF_321_kidsPhotometry_everything_col_flag_shear_m283m283_rot_90.fits
# number with shear m283p283 rot 0.0: 11166328
# catalogue saved as /disks/shear16/ssli/ImSim/output/skills_v07D7/skills_v07D7_LF_321_kidsPhotometry_everything_col_flag_shear_m283p283_rot_0.fits
# number with shear m283p283 rot 90.0: 11163075
# catalogue saved as /disks/shear16/ssli/ImSim/output/skills_v07D7/skills_v07D7_LF_321_kidsPhotometry_everything_col_flag_shear_m283p283_rot_90.fits
# number with shear p283m283 rot 0.0: 11164651
# catalogue saved as /disks/shear16/ssli/ImSim/output/skills_v07D7/skills_v07D7_LF_321_kidsPhotometry_everything_col_flag_shear_p283m283_rot_0.fits
# number with shear p283m283 rot 90.0: 11163644
# catalogue saved as /disks/shear16/ssli/ImSim/output/skills_v07D7/skills_v07D7_LF_321_kidsPhotometry_everything_col_flag_shear_p283m283_rot_90.fits
# number with shear p283p283 rot 0.0: 11163017
# catalogue saved as /disks/shear16/ssli/ImSim/output/skills_v07D7/skills_v07D7_LF_321_kidsPhotometry_everything_col_flag_shear_p283p283_rot_0.fits
# number with shear p283p283 rot 90.0: 11161260
# catalogue saved as /disks/shear16/ssli/ImSim/output/skills_v07D7/skills_v07D7_LF_321_kidsPhotometry_everything_col_flag_shear_p283p283_rot_90.fits
# number in split 89303975
# Elapsed:1:24:15.34,User=4014.117,System=4423.598,CPU=166.9%.