# -*- coding: utf-8 -*-
# @Author: lshuns
# @Date:   2022-09-04 16:13:56
# @Last Modified by:   lshuns
# @Last Modified time: 2022-09-30 13:18:53

### get the dm by applying constShear results first
#### do not bin in zTrue for constShear
#### in 2D bins (R and SNR)
###### algorithm: 
############### 1. get gOut(z) = gIn(z) * (1+m(z)) + c(z)
############### 2. get (1+m) = average (1+m(z))

import sys

import pandas as pd 
import numpy as np

import statsmodels.api as sm 

# >>>>>>>>>>>>>> I/O and setups

# the catalogues
# inpath_var = '/disks/shear16/ssli/ImSim/output/skills_v07D7v/skills_v07D7v_LF_321_shear_noSG_noWeiCut_newCut_var.feather'
inpath_var = '/disks/shear16/ssli/ImSim/output/skills_v07D7v/skills_v07D7v_LF_321_shear_noSG_noWeiCut_newCut_var_withextra.feather'
inpath_const = '/disks/shear16/ssli/ImSim/output/skills_v07D7v/skills_v07D7v_LF_321_shear_noSG_noWeiCut_newCut_const.feather'

# the binning bounds 
inpath_prefix_bounds = '../biasEstimation/outputs/m_weiRaw_DR4817Rewei_skills_v07D7_LF_321_kidsPhotometry_shear_noSG_noWeiCut_newCut'
bin1_col = 'SNR_LF_r'
bin2_col = 'R'

# output
outpath_prefix = './outputs/dm_forBlending'

# bin width in true z
dz = 0.1
zmin = 0
zmax = 2.5
zTrue_range = np.linspace(zmin, zmax, int((zmax-zmin)/dz+1))
print('defined z range', zTrue_range)
NzTrueBins = len(zTrue_range)-1
dz_label = str(dz).replace('.', 'p')

# the bins
col_binning = 'Z_B'
binning_edges = [0.1, 0.3, 0.5, 0.7, 0.9, 1.2, 2.0]

# working on input or measured shear
e_type = 'measured'

# column name for the weights
wei_col = 'oldweight_LF_r'

# >>>>>>>>>>>>>>> load catalogues

# load the catalogues
cata_list = []
for inpath in [inpath_var, inpath_const]:
    cata = pd.read_feather(inpath)
    print('>>> cata loaded from', inpath)
    ## only weight larger than 0
    cata = cata[(cata[wei_col]>0)]
    cata.reset_index(drop=True, inplace=True)
    ## used columns
    # get shear values
    if e_type == 'input':
        # perfect shear values from input e
        g = np.array(cata['g1_in']) + 1j*np.array(cata['g2_in'])
        e_in_gal = np.array(cata['e1_input']) + 1j*np.array(cata['e2_input'])
        e_true = (e_in_gal+g) / (1+np.conj(g)*e_in_gal)
        e1_out = (e_true.real).astype(float)
        e2_out = (e_true.imag).astype(float)
        del g, e_in_gal, e_true
        shape_weight = np.ones(len(e1_out))
    elif e_type == 'measured':
        # measured values
        e1_out = np.array(cata['e1_LF_r'])
        e2_out = np.array(cata['e2_LF_r'])
        shape_weight = cata[wei_col].values
    cata = pd.DataFrame({   bin1_col: cata[bin1_col].values,
                            bin2_col: cata[bin2_col].values,
                            'id_input': cata['id_input'].values,
                            'shape_weight': shape_weight,
                            'e1_out': e1_out,
                            'e2_out': e2_out,
                            'g1_in': cata['g1_in'].values,
                            'g2_in': cata['g2_in'].values,
                            'redshift_input': cata['redshift_input'].values,
                            col_binning: cata[col_binning].values
                            })
    del e1_out, e2_out, shape_weight

    cata_list.append(cata)
    del cata
cata_var, cata_const = cata_list
del cata_list

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

for i_bin in range(len(binning_edges)):

    # the whole results
    if i_bin == 0:
        outpath = outpath_prefix + f'_dz{dz_label}_whole.csv'
        bin_min = binning_edges[0]
        bin_max = binning_edges[-1]
        # 2D bins
        bin1_bounds = np.load(inpath_prefix_bounds + '_whole_bin1.npy')
        bin2_bounds = np.load(inpath_prefix_bounds + '_whole_bin2.npy')
    # others
    else:
        outpath = outpath_prefix + f'_dz{dz_label}_bin{i_bin-1}.csv'
        bin_min = binning_edges[i_bin - 1]
        bin_max = binning_edges[i_bin]
        # 2D bins
        bin1_bounds = np.load(inpath_prefix_bounds + f'_bin{i_bin-1}_bin1.npy')
        bin2_bounds = np.load(inpath_prefix_bounds + f'_bin{i_bin-1}_bin2.npy')

    bin1_Nbins = len(bin1_bounds) - 1
    bin2_Nbins = len(bin2_bounds[0]) - 1
    print('>> Number of 2D bins', bin1_Nbins, bin2_Nbins)

    ## select
    cata_var_selec = cata_var[(cata_var[col_binning]>bin_min)&(cata_var[col_binning]<=bin_max)].copy()
    cata_var_selec.reset_index(drop=True, inplace=True)
    cata_const_selec = cata_const[(cata_const[col_binning]>bin_min)&(cata_const[col_binning]<=bin_max)].copy()
    cata_const_selec.reset_index(drop=True, inplace=True)

    ## bin in the first col
    cata_var_selec.loc[:, 'bin1_id'] = pd.cut(cata_var_selec[bin1_col].values, bin1_bounds, 
                                right=True, labels=False)
    cata_const_selec.loc[:, 'bin1_id'] = pd.cut(cata_const_selec[bin1_col].values, bin1_bounds, 
                                right=True, labels=False)
    del bin1_bounds
    ## bin in the second col
    for ibin1, bin2_bound in enumerate(bin2_bounds):
        # mask in bin1
        mask_bin1_var = cata_var_selec['bin1_id'].values == ibin1
        mask_bin1_const = cata_const_selec['bin1_id'].values == ibin1
        # bin in bin2
        cata_var_selec.loc[mask_bin1_var, 'bin2_id'] = \
                            pd.cut(cata_var_selec.loc[mask_bin1_var, bin2_col].values, 
                                        bin2_bound, 
                                        right=True, labels=False)
        cata_const_selec.loc[mask_bin1_const, 'bin2_id'] = \
                            pd.cut(cata_const_selec.loc[mask_bin1_const, bin2_col].values,
                                        bin2_bound, 
                                        right=True, labels=False)
        del mask_bin1_var, mask_bin1_const

    # group
    N0 = len(cata_var_selec)
    cata_var_selec.dropna(inplace=True)
    print('var remaining after 2D binning', len(cata_var_selec)/N0)
    cata_var_selec = cata_var_selec.astype({'bin1_id': int, 'bin2_id': int})
    cata_var_selec = cata_var_selec.groupby(by=['bin1_id', 'bin2_id'])

    N0 = len(cata_const_selec)
    cata_const_selec.dropna(inplace=True)
    print('const remaining after 2D binning', len(cata_const_selec)/N0)
    cata_const_selec = cata_const_selec.astype({'bin1_id': int, 'bin2_id': int})
    cata_const_selec = cata_const_selec.groupby(by=['bin1_id', 'bin2_id'])

    # loop over bins and calculate m
    i_group = 0 
    mc_surface = pd.DataFrame(-999, 
                        index = np.arange(bin1_Nbins*bin2_Nbins), 
                        columns = ['bin1_id', 'bin2_id', 
                                    'bin1_mean', 'bin2_mean',
                                    'm1', 'm1_err', 'm2', 'm2_err',
                                    'c1', 'c1_err', 'c2', 'c2_err',
                                    'varWei'])
    for i in range(bin1_Nbins):
        for j in range(bin2_Nbins):

            # the name for group 
            name = (i, j)
            group_var = cata_var_selec.get_group(name)
            group_const = cata_const_selec.get_group(name)

            # bias results from constShear
            res = mCalFunc_pair_based(group_const)
            del group_const

            ## calculate m for z_true bins 
            res_df = pd.DataFrame({'zTrueMean': np.zeros(NzTrueBins),
                                    'm1': np.zeros(NzTrueBins),
                                    'm2': np.zeros(NzTrueBins),
                                    'c1': np.zeros(NzTrueBins),
                                    'c2': np.zeros(NzTrueBins),
                                    'm1_err': np.zeros(NzTrueBins),
                                    'm2_err': np.zeros(NzTrueBins),
                                    'c1_err': np.zeros(NzTrueBins),
                                    'c2_err': np.zeros(NzTrueBins),
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


                # apply the shear correction
                cata_var_tmp.loc[:, 'e1_out'] /= (1+res['m1']) 
                cata_var_tmp.loc[:, 'e2_out'] /= (1+res['m2']) 

                # calculate residual m
                res_residual = mCalFunc_pair_based(cata_var_tmp)
                del cata_var_tmp
                ## account for (1+m)
                res_residual['m1'] *= (1+res['m1']) 
                res_residual['m2'] *= (1+res['m2']) 
                res_residual['c1'] *= (1+res['m1']) 
                res_residual['c2'] *= (1+res['m2']) 
                res_residual['m1_err'] = (
                                            ((1+res['m1'])*res_residual['m1_err'])**2\
                                            + (res_residual['m1']*res['m1_err'])**2
                                            )**0.5
                res_residual['m2_err'] = (
                                            ((1+res['m2'])*res_residual['m2_err'])**2\
                                            + (res_residual['m2']*res['m2_err'])**2
                                            )**0.5
                res_residual['c1_err'] = (
                                            ((1+res['m1'])*res_residual['c1_err'])**2\
                                            + (res_residual['c1']*res['m1_err'])**2
                                            )**0.5
                res_residual['c2_err'] = (
                                            ((1+res['m2'])*res_residual['c2_err'])**2\
                                            + (res_residual['c2']*res['m2_err'])**2
                                            )**0.5

                # save
                res_df.loc[i_zbin, 'm1'] = res_residual['m1']
                res_df.loc[i_zbin, 'm2'] = res_residual['m2']
                res_df.loc[i_zbin, 'c1'] = res_residual['c1']
                res_df.loc[i_zbin, 'c2'] = res_residual['c2']
                res_df.loc[i_zbin, 'm1_err'] = res_residual['m1_err']
                res_df.loc[i_zbin, 'm2_err'] = res_residual['m2_err']
                res_df.loc[i_zbin, 'c1_err'] = res_residual['c1_err']
                res_df.loc[i_zbin, 'c2_err'] = res_residual['c2_err']
                del res_residual

            # get the average
            WeiNorm = np.sum(res_df['TotalWei'].values * dz)
            res_bin = {'m1': np.sum(res_df['m1'].values * res_df['TotalWei'].values * dz) / WeiNorm,
                                'm2': np.sum(res_df['m2'].values * res_df['TotalWei'].values * dz) / WeiNorm,
                                'c1': np.sum(res_df['c1'].values * res_df['TotalWei'].values * dz) / WeiNorm,
                                'c2': np.sum(res_df['c2'].values * res_df['TotalWei'].values * dz) / WeiNorm,
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
                                'c1_err': np.sqrt(
                                            np.sum(
                                                np.square(res_df['c1_err'].values * res_df['TotalWei'].values * dz)
                                            )
                                        ) / WeiNorm, 
                                'c2_err': np.sqrt(
                                            np.sum(
                                                np.square(res_df['c2_err'].values * res_df['TotalWei'].values * dz)
                                            )
                                        ) / WeiNorm
                                }
            del res_df

            # saving the surface
            mc_surface.loc[i_group, 'bin1_id'] = i
            mc_surface.loc[i_group, 'bin2_id'] = j
            mc_surface.loc[i_group, 'bin1_mean'] = np.average(group_var[bin1_col].values, weights=group_var['shape_weight'].values)
            mc_surface.loc[i_group, 'bin2_mean'] = np.average(group_var[bin2_col].values, weights=group_var['shape_weight'].values)
            mc_surface.loc[i_group, 'm1'] = res_bin['m1']
            mc_surface.loc[i_group, 'm1_err'] = res_bin['m1_err']
            mc_surface.loc[i_group, 'm2'] = res_bin['m2']
            mc_surface.loc[i_group, 'm2_err'] = res_bin['m2_err']
            mc_surface.loc[i_group, 'c1'] = res_bin['c1']
            mc_surface.loc[i_group, 'c1_err'] = res_bin['c1_err']
            mc_surface.loc[i_group, 'c2'] = res_bin['c2']
            mc_surface.loc[i_group, 'c2_err'] = res_bin['c2_err']
            mc_surface.loc[i_group, 'varWei'] = np.sum(group_var['shape_weight'].values)
            i_group += 1 
            del group_var

    # save
    mc_surface.to_csv(outpath, index=False)
    del mc_surface
    print('saved to', outpath)

# ######## out info
# defined z range [0.  0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1.  1.1 1.2 1.3 1.4 1.5 1.6 1.7
#  1.8 1.9 2.  2.1 2.2 2.3 2.4 2.5]
# >>> cata loaded from /disks/shear16/ssli/ImSim/output/skills_v07D7v/skills_v07D7v_LF_321_shear_noSG_noWeiCut_newCut_var_withextra.feather
# >>> cata loaded from /disks/shear16/ssli/ImSim/output/skills_v07D7v/skills_v07D7v_LF_321_shear_noSG_noWeiCut_newCut_const.feather
# >> Number of 2D bins 20 20
# var remaining after 2D binning 0.9999998699890628
# const remaining after 2D binning 0.999999628555743
# saved to ./outputs/dm_forBlending_dz0p1_whole.csv
# >> Number of 2D bins 20 20
# var remaining after 2D binning 0.9999938420509152
# const remaining after 2D binning 0.9999948690899613
# saved to ./outputs/dm_forBlending_dz0p1_bin0.csv
# >> Number of 2D bins 20 20
# var remaining after 2D binning 0.9999983438395481
# const remaining after 2D binning 0.999998300291568
# saved to ./outputs/dm_forBlending_dz0p1_bin1.csv
# >> Number of 2D bins 20 20
# var remaining after 2D binning 0.999997891913021
# const remaining after 2D binning 0.9999985704021201
# saved to ./outputs/dm_forBlending_dz0p1_bin2.csv
# >> Number of 2D bins 20 20
# var remaining after 2D binning 0.9999983101254003
# const remaining after 2D binning 0.999998069519685
# saved to ./outputs/dm_forBlending_dz0p1_bin3.csv
# >> Number of 2D bins 20 20
# var remaining after 2D binning 0.9999953989829452
# const remaining after 2D binning 0.9999948909289977
# saved to ./outputs/dm_forBlending_dz0p1_bin4.csv
# >> Number of 2D bins 20 20
# var remaining after 2D binning 0.9999964977328657
# const remaining after 2D binning 0.9999979999786665
# saved to ./outputs/dm_forBlending_dz0p1_bin5.csv
# Elapsed:1:01:29.16,User=8942.616,System=12848.242,CPU=590.6%.
