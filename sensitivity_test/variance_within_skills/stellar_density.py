# -*- coding: utf-8 -*-
# @Author: lshuns
# @Date:   2022-07-10 15:22:47
# @Last Modified by:   lshuns
# @Last Modified time: 2022-08-25 17:01:43

### check m difference between different stellar densities
########### three classes: low, median, high

import glob

import numpy as np
import pandas as pd
import statsmodels.api as sm 

from scipy import stats
from astropy.io import fits
from astropy.table import Table

# +++++++++++++++++++++++++++++ setups

# fiducial simulation
in_file_sim = '/disks/shear16/ssli/ImSim/output/skills_v07D7/skills_v07D7_LF_321_kidsPhotometry_shear_noSG_noWeiCut_newCut.feather'

# for saving results
outprefix = './outputs/dm_stellar_part'

######## general info
bin_edges = np.array([0.1, 0.3, 0.5, 0.7, 0.9, 1.2, 2.0])
col_binning_sim = 'Z_B'

col_label = 'tile_label'
col_e1, col_e2 = 'e1_LF_r', 'e2_LF_r'
col_g1, col_g2 = 'g1_in', 'g2_in'
col_weight_sim = 'oldweight_LF_r'

##################### functions

def mCalFunc_tile_based(cataSim, psf_frame=False):
    """
    NOTE: this can only be used when 
        1. the input shear are constant across the tile image
        2. number of tiles are sufficient large

    Calculate the residual shear bias for a given simulated catalogue
        it first calculates shear g_out for each tile (accounting for shape noise cancellation)
            then estimate bias using all tiles and shear inputs by requiring sum(g_out - (1+m)g_in) -> min

        Used columns and their names:
            tile_label: unique label for tile
            g1_in, g2_in: input shear
            e1_out, e2_out: measured ellipticity
            shape_weight: shape measurement weights
            if psf_frame:
                e1_psf, e2_psf: ellipticity of PSF
    """

    # out gal e
    e1_out = np.array(cataSim['e1_out'])
    e2_out = np.array(cataSim['e2_out'])

    # rotate to psf frame
    if psf_frame:
        print('Measured e will be rotated to PSF frame...')
        # out PSF 
        e1_psf = np.array(cataSim['e1_psf'])
        e2_psf = np.array(cataSim['e2_psf'])
        PSF_angle = np.arctan2(e2_psf, e1_psf)

        # rotate the ellipticities
        ct = np.cos(PSF_angle)
        st = np.sin(PSF_angle)
        e1_uc = e1_out*ct + e2_out*st
        e2_uc = -e1_out*st + e2_out*ct
        del e1_psf, e2_psf, PSF_angle, ct, st

        e1_out = e1_uc
        e2_out = e2_uc
        del e1_uc, e2_uc

    # build dataframe and select used columns
    cataSim = pd.DataFrame({'tile_label': np.array(cataSim['tile_label']),
                            'g1_in': np.array(cataSim['g1_in']),
                            'g2_in': np.array(cataSim['g2_in']),
                            'e1_out': e1_out,
                            'e2_out': e2_out,
                            'shape_weight': np.array(cataSim['shape_weight'])
                            })
    del e1_out, e2_out

    # sort to speed up
    cataSim.sort_values(by=['tile_label', 'g1_in', 'g2_in'], inplace=True)

    # prepare the weighted mean
    cataSim.loc[:, 'e1_out'] *= cataSim['shape_weight'].values
    cataSim.loc[:, 'e2_out'] *= cataSim['shape_weight'].values

    # group based on tile and shear
    cataSim = cataSim.groupby(['tile_label', 'g1_in', 'g2_in'], as_index=False).sum()
    if len(cataSim) < 10:
        raise Exception('less than 10 points for lsq, use pair_based!')
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
    del cataSim
    res = {'m1': m1, 'm2': m2,
            'c1': c1, 'c2': c2,
            'm1_err': m1_err, 'm2_err': m2_err,
            'c1_err': c1_err, 'c2_err': c2_err,
            }

    return res

# +++++++++++++++++++++++++++++ workhorse

# >>>> the whole
cata_used = pd.read_feather(in_file_sim)
print('Number of sources in whole', len(cata_used))
### select those within the edges
cata_used = cata_used[(cata_used[col_binning_sim]>bin_edges[0])&(cata_used[col_binning_sim]<=bin_edges[-1])]
print('     within the edges', len(cata_used))
### select non-zero weights
cata_used = cata_used[cata_used[col_weight_sim]>0]
print('     with weight > 0', len(cata_used))
# save useable parameters
cata_used = pd.DataFrame({
        'col_binning': np.array(cata_used[col_binning_sim]).astype(float),
        'tile_label': np.array(cata_used[col_label]).astype(str),
        'e1_out': np.array(cata_used[col_e1]).astype(float),
        'e2_out': np.array(cata_used[col_e2]).astype(float),
        'g1_in': np.array(cata_used[col_g1]).astype(float),
        'g2_in': np.array(cata_used[col_g2]).astype(float),
        'shape_weight': np.array(cata_used[col_weight_sim]).astype(float),
        'perfect_flag_star': np.array(cata_used['perfect_flag_star'].astype(int))})

# >>>>>>>>> split the whole catalogue based on the stellar density
## get number of stars in each tile
star_info = cata_used[['tile_label', 'perfect_flag_star']].groupby('tile_label', as_index=False).sum()
## split based on number of stars
star_info.loc[:, 'bin_star'] = pd.qcut(star_info['perfect_flag_star'].values, 3, 
                                    labels=False, retbins=False)
print('>>> mean densities:\n  ', star_info.groupby('bin_star').mean())
tile_labels_list = [star_info.loc[star_info['bin_star']==i_binstar, 'tile_label'].values for i_binstar in range(3)]
del star_info

# >>>>>>>>> output
f_list = [open(outprefix+f'_low.csv', 'w'), open(outprefix+f'_median.csv', 'w'), open(outprefix+f'_high.csv', 'w')]

# >>>>>>>>> m results

## the whole results
#### m from the whole sample
res = mCalFunc_tile_based(cata_used, psf_frame=False)

#### binning info
min_bin = bin_edges[0]
max_bin = bin_edges[-1]
#### loop over parts to get difference
for i_part, tile_labels in enumerate(tile_labels_list):
    # select objects in these part
    cata_used_tmp = cata_used[np.isin(cata_used['tile_label'].values, tile_labels)].copy()

    # basic info
    N_s = len(cata_used_tmp)
    mean_bin = np.average(cata_used_tmp['col_binning'].values, weights=cata_used_tmp['shape_weight'].values)
    median_bin = np.median(cata_used_tmp['col_binning'].values)

    # apply the shear correction
    cata_used_tmp.loc[:, 'e1_out'] /= (1+res['m1']) 
    cata_used_tmp.loc[:, 'e2_out'] /= (1+res['m2']) 

    # calculate m
    res_residual = mCalFunc_tile_based(cata_used_tmp, psf_frame=False)
    del cata_used_tmp
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

    # collect columns names and values
    cols = ','.join(list(res_residual.keys()))
    vals = ','.join(["{0:0.4f}".format(val) for val in res_residual.values()])
    cols = cols + f',{col_binning_sim}_min,{col_binning_sim}_max,{col_binning_sim}_mean,{col_binning_sim}_median,Nobj'
    vals = vals + f',{min_bin},{max_bin},{mean_bin},{median_bin},{N_s}'    
    print(cols, file=f_list[i_part])
    print(vals, file=f_list[i_part])

## the binning results
for i_bin in range(len(bin_edges)-1):
    min_bin = bin_edges[i_bin]
    max_bin = bin_edges[i_bin+1]

    # select simulations
    cata_selec = cata_used[(cata_used['col_binning']>min_bin) & (cata_used['col_binning']<=max_bin)]
    cata_selec.reset_index(drop=True, inplace=True)

    #### m from the whole sample
    res = mCalFunc_tile_based(cata_selec, psf_frame=False)

    #### loop over parts to get difference
    for i_part, tile_labels in enumerate(tile_labels_list):
        # select objects in these part
        cata_used_tmp = cata_selec[np.isin(cata_selec['tile_label'].values, tile_labels)].copy()

        # basic info
        N_s = len(cata_used_tmp)
        mean_bin = np.average(cata_used_tmp['col_binning'].values, weights=cata_used_tmp['shape_weight'].values)
        median_bin = np.median(cata_used_tmp['col_binning'].values)

        # apply the shear correction
        cata_used_tmp.loc[:, 'e1_out'] /= (1+res['m1']) 
        cata_used_tmp.loc[:, 'e2_out'] /= (1+res['m2']) 

        # calculate m
        res_residual = mCalFunc_tile_based(cata_used_tmp, psf_frame=False)
        del cata_used_tmp
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

        # print out values
        vals = ','.join(["{0:0.4f}".format(val) for val in res_residual.values()]) \
        + f',{min_bin},{max_bin},{mean_bin},{median_bin},{N_s}'
        print(vals, file=f_list[i_part])

    del cata_selec

# close what is opened
for f in f_list:
    f.close()
print(f'results saved as {outprefix}')

# # ############# outinfo
# Number of sources in whole 47870631
#      within the edges 47499698
#      with weight > 0 38907653
# >>> mean densities:
#              perfect_flag_star
# bin_star                   
# 0                959.000000
# 1               1848.444444
# 2               3360.500000
# results saved as ./outputs/dm_stellar_part
# Elapsed:15:06.80,User=693.780,System=1090.825,CPU=196.8%.