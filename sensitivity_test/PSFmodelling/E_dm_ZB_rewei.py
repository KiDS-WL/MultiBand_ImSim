# -*- coding: utf-8 -*-
# @Author: lshuns
# @Date:   2022-09-06 11:23:50
# @Last Modified by:   lshuns
# @Last Modified time: 2022-09-06 12:57:50

### calculate m difference between the test and fiducial
######## as a function of redshift bins
######### reweight PSFsize to match PSF distribution

import numpy as np
import pandas as pd
import statsmodels.api as sm 
import scipy.spatial as sst

############### I/O

# the fiducial catalogue
in_file_fiducial = '/disks/shear16/ssli/ImSim/output/skills_v07D7/skills_v07D7_LF_321_kidsPhotometry_shear_noSG_noWeiCut_newCut.feather'

# the PSF modelling cata
in_file_test = '/disks/shear16/ssli/ImSim/output/skills_v07D7_PSFmodelling/skills_v07D7_LF_321_combined_PSFmodelling41_321_noSG_noWeiCut_newCut.feather'

# the noise cata
in_file_noise = './tile_forTEST.csv'

# where to save results
outpath = './outputs/dm_ZB_test_fiducial_rewei_41.csv'

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

# >>>>>>>>> get the used tiles
tile_labels = pd.read_csv(in_file_noise)['label'].values
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
res_residual = dmCalFunc_tile_based(cata_test, cata_fiducial, nboot=500)
# basic info
N_s = len(cata_fiducial)
mean_bin = np.average(cata_fiducial[col_binning_sim].values, weights=cata_fiducial['shape_weight'].values)
# collect columns names and values
cols = ','.join(list(res_residual.keys()))
vals = ','.join(["{0:0.4f}".format(val) for val in res_residual.values()])
cols = cols + f',{col_binning_sim}_mean,Nobj'
vals = vals + f',{mean_bin},{N_s}'    
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
    mean_bin = np.average(cata_selec[col_binning_sim].values, weights=cata_selec['shape_weight'].values)

    # results
    res_residual = dmCalFunc_tile_based(cata_tmp, cata_selec, nboot=500)
    del cata_selec, cata_tmp

    # collect columns names and values
    vals = ','.join(["{0:0.4f}".format(val) for val in res_residual.values()]) \
     + f',{mean_bin},{N_s}'
    print(vals, file=f)

# close what is opened
f.close()
print(f'results saved as {outpath}')

# number of tiles used 30
# number fiducial ori 47870631
# number fiducial weight>0 39251930
# number fiducial selected 11454738
# number test ori 13784445
# number test weight>0 11396540
# >>> check boost factor 0.9999999999999999 [0.02601409 0.05262484 0.00868354 0.03314007 0.03482664 0.03094597
#  0.02953765 0.01875087 0.02798586 0.02361299 0.02303046 0.03698733
#  0.02908177 0.04410776 0.05484469 0.04508637 0.03056774 0.06033966
#  0.05515039 0.00881325 0.0252803  0.03529113 0.00511391 0.0405402
#  0.03863198 0.07627712 0.01136096 0.05182678 0.01954755 0.02199813]
# tile_label 136.0_1.5
# >>> boost 0.029537646995670137
# tile_label 180.0_-1.5
# >>> boost 0.04410775829082398
# tile_label 181.0_1.5
# >>> boost 0.04508637456026358
# tile_label 229.0_-1.5
# >>> boost 0.03056773708068469
# tile_label 8.2_-32.1
# >>> boost 0.01954754825161216
# tile_label 43.5_-32.1
# >>> boost 0.04054020410757336
# tile_label 333.8_-29.2
# >>> boost 0.03529113013588991
# tile_label 231.0_1.5
# >>> boost 0.05515038624641257
# tile_label 12.7_-30.2
# >>> boost 0.05262483677953524
# tile_label 134.0_1.5
# >>> boost 0.03314007233804659
# tile_label 9.5_-33.1
# >>> boost 0.021998134806630598
# tile_label 45.9_-32.1
# >>> boost 0.011360962585575984
# tile_label 180.5_2.5
# >>> boost 0.05484468630744402
# tile_label 136.0_0.5
# >>> boost 0.030945974933178858
# tile_label 229.0_0.5
# >>> boost 0.0603396603195526
# tile_label 45.4_-31.2
# >>> boost 0.03863198285721079
# tile_label 178.0_-1.5
# >>> boost 0.023030464820089025
# tile_label 10.6_-32.1
# >>> boost 0.026014088754807152
# tile_label 133.0_1.5
# >>> boost 0.008683539461412726
# tile_label 135.0_-1.5
# >>> boost 0.0348266432916619
# tile_label 45.6_-29.2
# >>> boost 0.07627711505078295
# tile_label 137.0_-0.5
# >>> boost 0.018750871022242252
# tile_label 232.0_-0.5
# >>> boost 0.00881325349494809
# tile_label 178.0_-0.5
# >>> boost 0.023612988699200678
# tile_label 46.3_-33.1
# >>> boost 0.05182677554703093
# tile_label 332.9_-32.1
# >>> boost 0.025280300092059074
# tile_label 337.0_-30.2
# >>> boost 0.005113908275371237
# tile_label 179.0_-1.5
# >>> boost 0.03698732611625274
# tile_label 179.0_1.5
# >>> boost 0.02908177103996911
# tile_label 137.0_-1.5
# >>> boost 0.027985857738066933
# >>>>>>> total boost 1.0
# >>> number of points for lsq 120
# >>> number of points for lsq 120
# >>> number of points for lsq 120
# >>> number of points for lsq 120
# >>> number of points for lsq 120
# >>> number of points for lsq 120
# >>> number of points for lsq 120
# results saved as ./outputs/dm_ZB_test_fiducial_rewei_41.csv
# Elapsed:5:43.44,User=308.888,System=1536.605,CPU=537.3%.