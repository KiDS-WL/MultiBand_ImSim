# -*- coding: utf-8 -*-
# @Author: lshuns
# @Date:   2022-09-28 16:15:41
# @Last Modified by:   lshuns
# @Last Modified time: 2022-09-28 16:57:04

### combine the extra running results back to the original result

import pandas as pd 

#### var
# all files to be combined
inpath_list = ['/disks/shear16/ssli/ImSim/output/skills_v07D7v/skills_v07D7v_LF_321_shear_noSG_noWeiCut_newCut_var.feather',
        '/disks/shear16/ssli/ImSim/output/skills_v07D7v_extra/skills_v07D7v_LF_321_shear_noSG_noWeiCut_newCut_var.feather']
# where to save
outpath = '/disks/shear16/ssli/ImSim/output/skills_v07D7v/skills_v07D7v_LF_321_shear_noSG_noWeiCut_newCut_var_withextra.feather'

# #### _var_withNeigIn
# # all files to be combined
# inpath_list = ['/disks/shear16/ssli/ImSim/output/skills_v07D7v/skills_v07D7v_LF_321_shear_noSG_noWeiCut_newCut_var_withNeigIn.feather',
#         '/disks/shear16/ssli/ImSim/output/skills_v07D7v_extra/skills_v07D7v_LF_321_shear_noSG_noWeiCut_newCut_var_withNeigIn.feather']
# # where to save
# outpath = '/disks/shear16/ssli/ImSim/output/skills_v07D7v/skills_v07D7v_LF_321_shear_noSG_noWeiCut_newCut_var_withNeigIn_withextra.feather'

# combine
cata_final = []
Ntot = 0
for inpath in inpath_list:
    cata = pd.read_feather(inpath)
    cata_final.append(cata)
    Ntmp = len(cata)
    del cata
    print('>>> sub number', Ntmp)
    Ntot += Ntmp
cata_final = pd.concat(cata_final, ignore_index=True)
print('+++ total', Ntot)

# save
print('cata final', len(cata_final))
cata_final.to_feather(outpath)
print('saved to', outpath)

# #### var
# >>> sub number 64116903
# >>> sub number 64110761
# +++ total 128227664
# cata final 128227664
# saved to /disks/shear16/ssli/ImSim/output/skills_v07D7v/skills_v07D7v_LF_321_shear_noSG_noWeiCut_newCut_var_withextra.feather
# Elapsed:7:41.70,User=852.579,System=2242.530,CPU=670.3%.

#### _var_withNeigIn
# >>> sub number 64116903
# >>> sub number 64110761
# +++ total 128227664
# cata final 128227664
# saved to /disks/shear16/ssli/ImSim/output/skills_v07D7v/skills_v07D7v_LF_321_shear_noSG_noWeiCut_newCut_var_withNeigIn_withextra.feather
# Elapsed:8:02.02,User=998.537,System=2996.564,CPU=828.8%.

