# -*- coding: utf-8 -*-
# @Author: lshuns
# @Date:   2023-04-10 17:26:01
# @Last Modified by:   lshuns
# @Last Modified time: 2023-04-11 17:01:56

### calculate m difference caused by the PSF modelling 
######## as a function of redshift bins
######### reweight PSFsize to match PSF distribution

import numpy as np
import pandas as pd
import statsmodels.api as sm 
import scipy.spatial as sst

############### I/O

## the fiducial catalogue
in_file_fiducial = '/disks/shear16/ssli/ImSim/output/skills_v07D7/skills_v07D7p1_LF_321_kidsPhotometry_shear_noSG_noWeiCut_newCut.feather'
## the PSF modelling cata
in_file_test = '/disks/shear16/ssli/ImSim/output/skills_v07D7p1_PSFmodelling/skills_v07D7p1_LF_321_combined_PSFmodelling41_321_noSG_noWeiCut_newCut.feather'
## where to save results
outpath = './results/dm_PSFmodelling_41_rewei.csv'

# general info
bin_edges = np.array([0.1, 0.3, 0.5, 0.7, 0.9, 1.2, 2.0])
col_binning_sim = 'Z_B'

############### function
def dmCalFunc_tile_based(cataSim0, cataSim1, nboot=500):
    """
    assume cataSim0 and cataSim0 having same tiles
    dm = cataSim0 - cataSim1
    """

    # >>>> average the first catalogue
    # build dataframe and select used columns
    cataSim0 = cataSim0[['tile_label', 'g1_in', 'g2_in',
                        'e1_out', 'e2_out', 'shape_weight']].copy()
    # sort to speed up
    cataSim0.sort_values(by=['tile_label', 'g1_in', 'g2_in'], inplace=True)
    # prepare the weighted mean
    cataSim0.loc[:, 'e1_out'] *= cataSim0['shape_weight'].values
    cataSim0.loc[:, 'e2_out'] *= cataSim0['shape_weight'].values
    # group based on input id and shear
    cataSim0 = cataSim0.groupby(['tile_label', 'g1_in', 'g2_in'], as_index=False).sum()
    Ntiles = len(cataSim0)
    print('>>> number of points for lsq', Ntiles)
    ## last step of the weighted mean
    cataSim0.loc[:, 'e1_out'] /= cataSim0['shape_weight'].values
    cataSim0.loc[:, 'e2_out'] /= cataSim0['shape_weight'].values

    # >>>> average the second catalogue
    # build dataframe and select used columns
    cataSim1 = cataSim1[['tile_label', 'g1_in', 'g2_in',
                        'e1_out', 'e2_out', 'shape_weight']].copy()
    # sort to speed up
    cataSim1.sort_values(by=['tile_label', 'g1_in', 'g2_in'], inplace=True)
    # prepare the weighted mean
    cataSim1.loc[:, 'e1_out'] *= cataSim1['shape_weight'].values
    cataSim1.loc[:, 'e2_out'] *= cataSim1['shape_weight'].values
    # group based on input id and shear
    cataSim1 = cataSim1.groupby(['tile_label', 'g1_in', 'g2_in'], as_index=False).sum()
    if len(cataSim1) != Ntiles:
        raise Exception(f'number of tiles different between cataSim0 and cataSim1')
    ## last step of the weighted mean
    cataSim1.loc[:, 'e1_out'] /= cataSim1['shape_weight'].values
    cataSim1.loc[:, 'e2_out'] /= cataSim1['shape_weight'].values

    # >>>> get the difference
    ## check the tile and input shears are the same
    if (np.sum(cataSim0['tile_label'].values == cataSim1['tile_label'].values)!=Ntiles) \
        | (np.sum(cataSim0['g1_in'].values == cataSim1['g1_in'].values)!=Ntiles) \
        | (np.sum(cataSim0['g2_in'].values == cataSim1['g2_in'].values)!=Ntiles):
        raise Exception('different tiles or input shears for the two cata')
    ## shear difference
    cataSim = pd.DataFrame({
        'g1_in': cataSim0['g1_in'].values, 
        'g2_in': cataSim0['g2_in'].values,
        'e1_out': cataSim0['e1_out'].values - cataSim1['e1_out'].values, 
        'e2_out': cataSim0['e2_out'].values - cataSim1['e2_out'].values, 
        'shape_weight': (cataSim0['shape_weight'].values + cataSim1['shape_weight'].values)/2.
        })
    del cataSim0, cataSim1

    # >>> get least square values
    ## e1
    mod_wls = sm.WLS(cataSim['e1_out'].values, \
                        sm.add_constant(cataSim['g1_in'].values), \
                        weights=cataSim['shape_weight'].values)
    res_wls = mod_wls.fit()
    m1 = res_wls.params[1]
    c1 = res_wls.params[0]
    m1_err = (res_wls.cov_params()[1, 1])**0.5
    c1_err = (res_wls.cov_params()[0, 0])**0.5
    ## e2
    mod_wls = sm.WLS(cataSim['e2_out'].values, \
                        sm.add_constant(cataSim['g2_in'].values), \
                        weights=cataSim['shape_weight'].values)
    res_wls = mod_wls.fit()
    m2 = res_wls.params[1]
    c2 = res_wls.params[0]
    m2_err = (res_wls.cov_params()[1, 1])**0.5
    c2_err = (res_wls.cov_params()[0, 0])**0.5

    # get error from boots
    if nboot is not None:
        ## sample with weights
        df_BS_list = [cataSim.sample(n=Ntiles, replace=True, weights=cataSim['shape_weight'].values, ignore_index=True) 
                            for _ in range(nboot)]
        del cataSim
        m1_BS = np.zeros(nboot)
        m2_BS = np.zeros(nboot)
        for i_BS, df_BS in enumerate(df_BS_list):

            # get least square values
            ## e1
            mod_wls = sm.WLS(df_BS['e1_out'].values, \
                                sm.add_constant(df_BS['g1_in'].values), \
                                weights=df_BS['shape_weight'].values)
            res_wls = mod_wls.fit()
            m1_BS[i_BS] = res_wls.params[1]
            ## e2
            mod_wls = sm.WLS(df_BS['e2_out'].values, \
                                sm.add_constant(df_BS['g2_in'].values), \
                                weights=df_BS['shape_weight'].values)
            res_wls = mod_wls.fit()
            m2_BS[i_BS] = res_wls.params[1]
            del df_BS
        del df_BS_list

        ## difference 
        sigma_m1_BS = m1_BS - m1
        sigma_m2_BS = m2_BS - m2
        del m1_BS, m2_BS
        ## choose the 1 sigma range as error
        m1_err_BS = ( np.abs(np.percentile(sigma_m1_BS, 16))
                        + np.abs(np.percentile(sigma_m1_BS, 84))
                        )/2.
        m2_err_BS = ( np.abs(np.percentile(sigma_m2_BS, 16))
                        + np.abs(np.percentile(sigma_m2_BS, 84))
                        )/2.
    else:
        m1_err_BS = 0
        m2_err_BS = 0

    # save
    res = {'m1': m1, 'm2': m2,
            'c1': c1, 'c2': c2,
            'm1_err': m1_err, 'm2_err': m2_err,
            'c1_err': c1_err, 'c2_err': c2_err,
            'm1_err_BS': m1_err_BS, 'm2_err_BS': m2_err_BS,
            }

    return res

############### workhorse

# >>>>>>>>> test catalogue
cata_test = pd.read_feather(in_file_test)
print('number test ori', len(cata_test))
## used columns
cata_test = cata_test[['id_input', 'PSFsize', col_binning_sim,
                            'oldweight_LF_r', 'e1_LF_r', 'e2_LF_r', 
                            'g1_in', 'g2_in', 'tile_label', 'gal_rot']]
## weight cut
cata_test = cata_test[cata_test['oldweight_LF_r'].values>0]
cata_test.reset_index(drop=True, inplace=True)
print('number test weight>0', len(cata_test))
## rename
cata_test.rename(columns={'oldweight_LF_r': 'shape_weight',
                    'e1_LF_r': 'e1_out', 'e2_LF_r': 'e2_out'}, inplace=True)
## get the used tiles
tile_labels = np.unique(cata_test['tile_label'].values)
print('number of tiles used', len(tile_labels))

# >>>>>>>>> fiducial catalogue
cata0 = pd.read_feather(in_file_fiducial)
print('number fiducial ori', len(cata0))
## weight cut
cata0 = cata0[cata0['oldweight_LF_r'].values>0]
cata0.reset_index(drop=True, inplace=True)
print('number fiducial weight>0', len(cata0))
## select within range
cata_fiducial = cata0[np.isin(cata0['tile_label'].values, tile_labels)].copy()
cata_fiducial.reset_index(drop=True, inplace=True)
print('number fiducial selected', len(cata_fiducial))
## used columns for m cal
cata_fiducial = cata_fiducial[['id_input', 'PSFsize', col_binning_sim,
                            'oldweight_LF_r', 'e1_LF_r', 'e2_LF_r', 
                            'g1_in', 'g2_in', 'tile_label', 'gal_rot']]
### rename
cata_fiducial.rename(columns={'oldweight_LF_r': 'shape_weight',
                    'e1_LF_r': 'e1_out', 'e2_LF_r': 'e2_out'}, inplace=True)
## used columns for reweighting
cata0 = cata0[['oldweight_LF_r', 'PSFsize', 'tile_label']].copy()

# >>>>>>>>> reweighting
# prepare the weighted mean
cata0.loc[:, 'PSFsize'] *= cata0['oldweight_LF_r'].values
cata0 = cata0.groupby(['tile_label'], as_index=False).sum()
cata0.loc[:, 'PSFsize'] /= cata0['oldweight_LF_r'].values
## get weighted mean for each tile in fiducial
cata_fiducial_wei = cata_fiducial[['shape_weight', 'PSFsize', 'tile_label']].copy()
cata_fiducial_wei.loc[:, 'PSFsize'] *= cata_fiducial_wei['shape_weight'].values
cata_fiducial_wei = cata_fiducial_wei.groupby(['tile_label']).sum()
cata_fiducial_wei.loc[:, 'PSFsize'] /= cata_fiducial_wei['shape_weight'].values
## find nearest
### build KDTree
kdt = sst.cKDTree(np.array([cata_fiducial_wei['PSFsize'].values]).T, leafsize=100)
### query
size_ori = np.array([cata0['PSFsize'].values]).T
dist, ind = kdt.query(size_ori, k=1, distance_upper_bound=np.inf)
del size_ori, dist
cata0.loc[:, 'tile_label_nearest'] = cata_fiducial_wei.index.values[ind].flatten()
del ind, cata_fiducial_wei
## get number of boost 
cata0 = cata0.groupby(['tile_label_nearest']).sum()
cata0.loc[:, 'oldweight_LF_r'] /= np.sum(cata0['oldweight_LF_r'].values)
print('>>> check boost factor', np.sum(cata0['oldweight_LF_r'].values), cata0['oldweight_LF_r'].values)
## apply the boost
Ntotboost = 0
for tile_label in tile_labels:
    ## boost from the base
    totWei0 = np.sum(cata0.loc[tile_label, 'oldweight_LF_r'])
    Ntotboost += totWei0
    print('tile_label', tile_label)
    print('>>> boost', totWei0)
    ## apply boost
    cata_test.loc[(cata_test['tile_label']==tile_label), 'shape_weight'] *= totWei0
    cata_fiducial.loc[(cata_fiducial['tile_label']==tile_label), 'shape_weight'] *= totWei0
del cata0
print('>>>>>>> total boost', Ntotboost)

# >>>>>>>>> output
f = open(outpath, 'w')

# >>>>>>>>> m results
## ++++ the whole results
### only objects within the whole range
min_bin = bin_edges[0]
max_bin = bin_edges[-1]
# select cata
cata_selec = cata_fiducial[(cata_fiducial[col_binning_sim]>min_bin) & (cata_fiducial[col_binning_sim]<=max_bin)].copy()
cata_selec.reset_index(drop=True, inplace=True)
cata_tmp = cata_test[(cata_test[col_binning_sim]>min_bin) & (cata_test[col_binning_sim]<=max_bin)].copy()
cata_tmp.reset_index(drop=True, inplace=True)
# basic info
N_s = len(cata_selec)
# results
res_residual = dmCalFunc_tile_based(cata_tmp, cata_selec, nboot=500)
del cata_selec, cata_tmp
# collect columns names and values
cols = ','.join(list(res_residual.keys()))
vals = ','.join(["{0:0.4f}".format(val) for val in res_residual.values()])
cols = cols + f',{col_binning_sim}_min,{col_binning_sim}_max,Nobj'
vals = vals + f',{min_bin},{max_bin},{N_s}'    
print(cols, file=f)
print(vals, file=f)
## +++ the binning results
for i_bin in range(len(bin_edges)-1):
    min_bin = bin_edges[i_bin]
    max_bin = bin_edges[i_bin+1]
    # select cata
    cata_selec = cata_fiducial[(cata_fiducial[col_binning_sim]>min_bin) & (cata_fiducial[col_binning_sim]<=max_bin)].copy()
    cata_selec.reset_index(drop=True, inplace=True)
    cata_tmp = cata_test[(cata_test[col_binning_sim]>min_bin) & (cata_test[col_binning_sim]<=max_bin)].copy()
    cata_tmp.reset_index(drop=True, inplace=True)
    # basic info
    N_s = len(cata_selec)
    # results
    res_residual = dmCalFunc_tile_based(cata_tmp, cata_selec, nboot=500)
    del cata_selec, cata_tmp
    # collect columns names and values
    vals = ','.join(["{0:0.4f}".format(val) for val in res_residual.values()]) \
     + f',{min_bin},{max_bin},{N_s}'
    print(vals, file=f)

# close what is opened
f.close()
print(f'results saved as {outpath}')

# ...
# results saved as ./results/dm_PSFmodelling_41_rewei.csv
# Elapsed:6:27.95,User=350.545,System=1710.330,CPU=531.2%.
