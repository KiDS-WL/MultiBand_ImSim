# -*- coding: utf-8 -*-
# @Author: lshuns
# @Date:   2023-05-02 14:01:45
# @Last Modified by:   lshuns
# @Last Modified time: 2023-05-03 16:20:41

### get the dm by applying constShear results first
#### do not bin in zTrue for constShear
#### in 2D bins (R and SNR) for each tomographic bin
###### algorithm: 
############### 1. get gOut(z) = gIn(z) * (1+m(z)) + c(z)
############### 2. get (1+m) = average (1+m(z))

import sys

import pandas as pd 
import numpy as np

import statsmodels.api as sm 

# >>>>>>>>>>>>>> I/O and setups

# the surface of doom
inpath_doom = './SurfaceOfDoom/m_surface_ZB0p1_SNR20_R20_skills_v07D7p1_LF_321_kidsPhotometry_gold.csv'

# the blending-only SKiLLS-gold catalogues
inpath_var = '/disks/shear16/ssli/ImSim/output/skills_v07D7v/skills_v07D7p1v_LF_321_shear_noSG_noWeiCut_newCut_var_withextra_goldSelected.feather'
inpath_const = '/disks/shear16/ssli/ImSim/output/skills_v07D7v/skills_v07D7p1v_LF_321_shear_noSG_noWeiCut_newCut_const_goldSelected.feather'

# bin width in true z
dz = 0.1
zmin = 0
zmax = 2.5
zTrue_range = np.linspace(zmin, zmax, int((zmax-zmin)/dz+1))
print('defined z range', zTrue_range)
NzTrueBins = len(zTrue_range)-1

# working on input or measured shear
e_type = 'measured'
## column name for the weights
col_weight = 'oldweight_LF_r'
col_e1 = 'e1_LF_r'
col_e2 = 'e2_LF_r'

# ### running info
# defined z range [0.  0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1.  1.1 1.2 1.3 1.4 1.5 1.6 1.7
#  1.8 1.9 2.  2.1 2.2 2.3 2.4 2.5]
# >>> cata loaded from /disks/shear16/ssli/ImSim/output/skills_v07D7v/skills_v07D7p1v_LF_321_shear_noSG_noWeiCut_newCut_var_withextra_goldSelected.feather
# number ori 75895923
# >>> cata loaded from /disks/shear16/ssli/ImSim/output/skills_v07D7v/skills_v07D7p1v_LF_321_shear_noSG_noWeiCut_newCut_const_goldSelected.feather
# number ori 9508161
# number within doom (var) 75830927
# number within doom (const) 9500010
# results saved back to ./SurfaceOfDoom/m_surface_ZB0p1_SNR20_R20_skills_v07D7p1_LF_321_kidsPhotometry_gold.csv
# Elapsed:33:21.08,User=4972.046,System=6362.446,CPU=566.4%.

# >>>>>>>>>>>>>>> load catalogues

# load the catalogues
cata_list = []
for inpath in [inpath_var, inpath_const]:
    cata = pd.read_feather(inpath)
    print('>>> cata loaded from', inpath)
    ## only weight larger than 0
    cata = cata[(cata[col_weight]>0)]
    cata.reset_index(drop=True, inplace=True)
    ## used columns
    # get shear values
    if e_type == 'input':
        # perfect shear values from input e
        g = np.array(cata['g1_in']) + 1j*np.array(cata['g2_in'])
        e_in_gal = np.array(cata[col_e1]) + 1j*np.array(cata[col_e2])
        e_true = (e_in_gal+g) / (1+np.conj(g)*e_in_gal)
        e1_out = (e_true.real).astype(float)
        e2_out = (e_true.imag).astype(float)
        del g, e_in_gal, e_true
        shape_weight = np.ones(len(e1_out))
    elif e_type == 'measured':
        # measured values
        e1_out = np.array(cata[col_e1])
        e2_out = np.array(cata[col_e2])
        shape_weight = cata[col_weight].values
    cata = pd.DataFrame({   'Z_B': cata['Z_B'].values,
                            'SNR': cata['SNR_LF_r'].values,
                            'R': cata['R'].values,
                            'id_input': cata['id_input'].values,
                            'shape_weight': shape_weight,
                            'e1_out': e1_out,
                            'e2_out': e2_out,
                            'g1_in': cata['g1_in'].values,
                            'g2_in': cata['g2_in'].values,
                            'redshift_input': cata['redshift_input'].values,
                            'binZB_id': -999,
                            'binSNR_id': -999,
                            'binR_id': -999
                            })
    del e1_out, e2_out, shape_weight
    print('number ori', len(cata))

    cata_list.append(cata)
    del cata
cata_var, cata_const = cata_list
del cata_list

# load surface of doom
mc_surface = pd.read_csv(inpath_doom)

# >>>>>>>>>>>>>> m calculate function
def mCalFunc_pair_based(cataSim):

    # build dataframe and select used columns
    cataSim = cataSim[['id_input', 'g1_in', 'g2_in',
                        'e1_out', 'e2_out', 'shape_weight']].copy()

    # sort to speed up
    cataSim.sort_values(by=['id_input', 'g1_in', 'g2_in'], inplace=True)

    # prepare the weighted mean
    cataSim.loc[:, 'e1_out'] *= cataSim['shape_weight'].values
    cataSim.loc[:, 'e2_out'] *= cataSim['shape_weight'].values

    # group based on tile and shear
    cataSim = cataSim.groupby(['id_input', 'g1_in', 'g2_in'], as_index=False).sum()
    if len(cataSim) < 3:
        raise Exception('less than 3 points for lsq, use pair_based!')
    # print('number of groups (points) for lsq', len(cataSim))
    ## last step of the weighted mean
    cataSim.loc[:, 'e1_out'] /= cataSim['shape_weight'].values
    cataSim.loc[:, 'e2_out'] /= cataSim['shape_weight'].values

    # get least square values
    ## e1
    mod_wls = sm.WLS(cataSim['e1_out'].values, \
                        sm.add_constant(cataSim['g1_in'].values), \
                        weights=cataSim['shape_weight'].values)
    res_wls = mod_wls.fit()
    m1 = res_wls.params[1] - 1
    c1 = res_wls.params[0]
    m1_err = (res_wls.cov_params()[1, 1])**0.5
    c1_err = (res_wls.cov_params()[0, 0])**0.5
    ## e2
    mod_wls = sm.WLS(cataSim['e2_out'].values, \
                        sm.add_constant(cataSim['g2_in'].values), \
                        weights=cataSim['shape_weight'].values)
    res_wls = mod_wls.fit()
    m2 = res_wls.params[1] - 1
    c2 = res_wls.params[0]
    m2_err = (res_wls.cov_params()[1, 1])**0.5
    c2_err = (res_wls.cov_params()[0, 0])**0.5

    # save
    res = {'m1': m1, 'm2': m2,
            'c1': c1, 'c2': c2,
            'm1_err': m1_err, 'm2_err': m2_err,
            'c1_err': c1_err, 'c2_err': c2_err
            }

    return res

# >>>>>>>>>>>>> calculate m

# bin galaxies
Z_B_edges = np.unique(mc_surface['binZB_min']).tolist() + [np.max(mc_surface['binZB_max'])]
cata_var.loc[:, 'binZB_id'] = pd.cut(cata_var['Z_B'].values, Z_B_edges, 
                                right=True, labels=False)
cata_const.loc[:, 'binZB_id'] = pd.cut(cata_const['Z_B'].values, Z_B_edges, 
                                right=True, labels=False)
for i_zbin in range(len(Z_B_edges)-1):
    SNR_edges = np.unique(mc_surface.loc[(mc_surface['binZB_id']==i_zbin), 
                    'binSNR_min']).tolist() \
                    + [np.max(mc_surface.loc[(mc_surface['binZB_id']==i_zbin), 
                        'binSNR_max'])]
    ## the varShear cata
    mask_binZB_var = cata_var['binZB_id'].values == i_zbin
    cata_var.loc[mask_binZB_var, 'binSNR_id'] = pd.cut(cata_var.loc[mask_binZB_var, 'SNR'].values, SNR_edges, 
                                right=True, labels=False)
    ## the constShear cata
    mask_binZB_const = cata_const['binZB_id'].values == i_zbin
    cata_const.loc[mask_binZB_const, 'binSNR_id'] = pd.cut(cata_const.loc[mask_binZB_const, 'SNR'].values, SNR_edges, 
                                right=True, labels=False)

    for i_SNR in range(len(SNR_edges)-1):
        R_edges = np.unique(mc_surface.loc[(mc_surface['binZB_id']==i_zbin)&(mc_surface['binSNR_id']==i_SNR), 
                    'binR_min']).tolist() \
                    + [np.max(mc_surface.loc[(mc_surface['binZB_id']==i_zbin)&(mc_surface['binSNR_id']==i_SNR), 
                        'binR_max'])]
        ## the varShear
        mask_binSNR_var = cata_var['binSNR_id'].values == i_SNR
        cata_var.loc[mask_binZB_var&mask_binSNR_var, 'binR_id'] = pd.cut(
                                    cata_var.loc[mask_binZB_var&mask_binSNR_var, 'R'].values, 
                                    R_edges, 
                                    right=True, labels=False)
        del mask_binSNR_var
        ## the constShear
        mask_binSNR_const = cata_const['binSNR_id'].values == i_SNR
        cata_const.loc[mask_binZB_const&mask_binSNR_const, 'binR_id'] = pd.cut(
                                    cata_const.loc[mask_binZB_const&mask_binSNR_const, 'R'].values, 
                                    R_edges, 
                                    right=True, labels=False)
        del mask_binSNR_const
    del mask_binZB_var, mask_binZB_const

# group
## drop -999 bins
cata_var = cata_var[(cata_var['binZB_id']>-999)&(cata_var['binSNR_id']>-999)&(cata_var['binR_id']>-999)].copy()
cata_var.reset_index(drop=True, inplace=True)
print('number within doom (var)', len(cata_var))
## sort to speed up
cata_var = cata_var.astype({'binZB_id': int, 'binSNR_id': int, 'binR_id': int})
cata_var.sort_values(by=['binZB_id', 'binSNR_id', 'binR_id'], inplace=True)
cata_var = cata_var.groupby(by=['binZB_id', 'binSNR_id', 'binR_id'])
## drop -999 bins
cata_const = cata_const[(cata_const['binZB_id']>-999)&(cata_const['binSNR_id']>-999)&(cata_const['binR_id']>-999)].copy()
cata_const.reset_index(drop=True, inplace=True)
print('number within doom (const)', len(cata_const))
## sort to speed up
cata_const = cata_const.astype({'binZB_id': int, 'binSNR_id': int, 'binR_id': int})
cata_const.sort_values(by=['binZB_id', 'binSNR_id', 'binR_id'], inplace=True)
cata_const = cata_const.groupby(by=['binZB_id', 'binSNR_id', 'binR_id'])

# loop over groups to get fraction
for name, group_var in cata_var:
    group_const = cata_const.get_group(name)
    binZB_id, binSNR_id, binR_id = name

    ## m from constShear
    res = mCalFunc_pair_based(group_const)
    del group_const

    ## dm for z_true bins 
    res_df = pd.DataFrame({'zTrueMean': np.zeros(NzTrueBins),
                            'm1': np.zeros(NzTrueBins),
                            'm2': np.zeros(NzTrueBins),
                            'm1_err': np.zeros(NzTrueBins),
                            'm2_err': np.zeros(NzTrueBins),
                            'TotalWei': np.zeros(NzTrueBins)})
    for i_zbin in range(NzTrueBins):
        zbin_min = zTrue_range[i_zbin]
        zbin_max = zTrue_range[i_zbin + 1]
        mask_tmp_var = (group_var['redshift_input']>zbin_min)&(group_var['redshift_input']<=zbin_max)
        if np.sum(mask_tmp_var) < 10:
            continue

        cata_var_tmp = group_var[mask_tmp_var].copy()
        del mask_tmp_var

        res_df.loc[i_zbin, 'TotalWei'] = np.sum(cata_var_tmp['shape_weight'].values)
        res_df.loc[i_zbin, 'zTrueMean'] = np.average(cata_var_tmp['redshift_input'].values, 
                                                weights=cata_var_tmp['shape_weight'].values)

        # apply the constant shear bias correction
        cata_var_tmp.loc[:, 'e1_out'] /= (1+res['m1']) 
        cata_var_tmp.loc[:, 'e2_out'] /= (1+res['m2']) 

        # calculate residual m
        res_residual = mCalFunc_pair_based(cata_var_tmp)
        del cata_var_tmp
        ## account for (1+m)
        res_residual['m1'] *= (1+res['m1']) 
        res_residual['m2'] *= (1+res['m2']) 
        res_residual['m1_err'] = (
                                    ((1+res['m1'])*res_residual['m1_err'])**2\
                                    + (res_residual['m1']*res['m1_err'])**2
                                    )**0.5
        res_residual['m2_err'] = (
                                    ((1+res['m2'])*res_residual['m2_err'])**2\
                                    + (res_residual['m2']*res['m2_err'])**2
                                    )**0.5

        # save
        res_df.loc[i_zbin, 'm1'] = res_residual['m1']
        res_df.loc[i_zbin, 'm2'] = res_residual['m2']
        res_df.loc[i_zbin, 'm1_err'] = res_residual['m1_err']
        res_df.loc[i_zbin, 'm2_err'] = res_residual['m2_err']
        del res_residual
    del group_var

    # get the average for the selected sample
    WeiNorm = np.sum(res_df['TotalWei'].values * dz)
    res_bin = {'m1': np.sum(res_df['m1'].values * res_df['TotalWei'].values * dz) / WeiNorm,
                        'm2': np.sum(res_df['m2'].values * res_df['TotalWei'].values * dz) / WeiNorm,
                        'm1_err': np.sqrt(
                                    np.sum(
                                        np.square(res_df['m1_err'].values * res_df['TotalWei'].values * dz)
                                    )
                                ) / WeiNorm, 
                        'm2_err': np.sqrt(
                                    np.sum(
                                        np.square(res_df['m2_err'].values * res_df['TotalWei'].values * dz)
                                    )
                                ) / WeiNorm, 
                        }
    del res_df

    # save the results
    mask_doom = (mc_surface['binZB_id']==binZB_id)\
            &(mc_surface['binSNR_id']==binSNR_id)\
            &(mc_surface['binR_id']==binR_id)
    mc_surface.loc[mask_doom, 'var_dm1'] = res_bin['m1']
    mc_surface.loc[mask_doom, 'var_dm2'] = res_bin['m2']
    mc_surface.loc[mask_doom, 'var_dm1_err'] = res_bin['m1_err']
    mc_surface.loc[mask_doom, 'var_dm2_err'] = res_bin['m2_err']
    del res_bin

# save back
# print('>>> sorted var_dm1', np.sort(mc_surface['var_dm1'].values))
# print('>>> sorted var_dm2', np.sort(mc_surface['var_dm2'].values))
# print('>>> sorted var_dm1_err', np.sort(mc_surface['var_dm1_err'].values))
# print('>>> sorted var_dm2_err', np.sort(mc_surface['var_dm2_err'].values))
mc_surface.to_csv(inpath_doom, index=False, float_format='%.6f')
print(f'results saved back to {inpath_doom}')