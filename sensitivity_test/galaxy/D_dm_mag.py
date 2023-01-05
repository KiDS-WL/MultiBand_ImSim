# -*- coding: utf-8 -*-
# @Author: lshuns
# @Date:   2022-09-06 12:33:44
# @Last Modified by:   lshuns
# @Last Modified time: 2022-09-06 14:29:20

### compare m between the test and fiducial
######## as a function of magnitude

import numpy as np
import pandas as pd
import statsmodels.api as sm 

from scipy import stats
from astropy.io import fits
from astropy.table import Table

############### I/O

# the test tags
test_tags = ['sizeU5', 'sizeU10', 'sizeU15',
                'qU5', 'qU10', 'qU15',
                'nU10', 'nU20', 'nU30']

# the fiducial catalogue
in_file_fiducial_list = [f'/disks/shear16/ssli/ImSim/output/skills_v07D7_{test_tag}/skills_v07D7_LF_321_part0_noSG_noWeiCut_newCut_fiducial.feather'
                            for test_tag in test_tags]

# the test catalogue
in_file_test_list = [f'/disks/shear16/ssli/ImSim/output/skills_v07D7_{test_tag}/skills_v07D7_LF_321_part0_noSG_noWeiCut_newCut_test.feather'
                            for test_tag in test_tags]

# where to save results
outpath_list = [f'./outputs/dm_{test_tag}_magauto.csv' for test_tag in test_tags]

# binning info
bin_col = 'MAG_AUTO'
bin_edges = np.arange(20, 26, 1)

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

for i_tag, in_file_fiducial in enumerate(in_file_fiducial_list):

    in_file_test = in_file_test_list[i_tag]
    outpath = outpath_list[i_tag]
    print('>>>> working on', test_tags[i_tag])

    # >>>>>>>>> fiducial catalogue
    cata_fiducial = pd.read_feather(in_file_fiducial)
    print('number fiducial ori', len(cata_fiducial))
    ## used columns
    cata_fiducial = cata_fiducial[[bin_col,
                                'oldweight_LF_r', 'e1_LF_r', 'e2_LF_r', 
                                'g1_in', 'g2_in', 'tile_label', 'gal_rot']]
    ## weight cut
    cata_fiducial = cata_fiducial[cata_fiducial['oldweight_LF_r'].values>0]
    cata_fiducial.reset_index(drop=True, inplace=True)
    print('number fiducial weight>0', len(cata_fiducial))
    ## rename
    cata_fiducial.rename(columns={'oldweight_LF_r': 'shape_weight',
                        'e1_LF_r': 'e1_out', 'e2_LF_r': 'e2_out'}, inplace=True)

    # >>>>>>>>> test catalogue
    cata_test = pd.read_feather(in_file_test)
    print('number test ori', len(cata_test))
    ## used columns
    cata_test = cata_test[[bin_col,
                                'oldweight_LF_r', 'e1_LF_r', 'e2_LF_r', 
                                'g1_in', 'g2_in', 'tile_label', 'gal_rot']]
    ## weight cut
    cata_test = cata_test[cata_test['oldweight_LF_r'].values>0]
    cata_test.reset_index(drop=True, inplace=True)
    print('number test weight>0', len(cata_test))
    ## rename
    cata_test.rename(columns={'oldweight_LF_r': 'shape_weight',
                        'e1_LF_r': 'e1_out', 'e2_LF_r': 'e2_out'}, inplace=True)

    # >>>>>>>>> output
    f = open(outpath, 'w')

    # >>>>>>>>> m results
    ## ++++ the whole results
    res_residual = dmCalFunc_tile_based(cata_test, cata_fiducial, nboot=500)
    # basic info
    N_s = len(cata_fiducial)
    mean_bin = np.average(cata_fiducial[bin_col].values, weights=cata_fiducial['shape_weight'].values)
    # collect columns names and values
    cols = ','.join(list(res_residual.keys()))
    vals = ','.join(["{0:0.4f}".format(val) for val in res_residual.values()])
    cols = cols + f',{bin_col}_mean,Nobj'
    vals = vals + f',{mean_bin},{N_s}'    
    print(cols, file=f)
    print(vals, file=f)

    ## +++ the binning results
    for i_bin in range(len(bin_edges)-1):
        min_bin = bin_edges[i_bin]
        max_bin = bin_edges[i_bin+1]

        # select cata
        cata_selec = cata_fiducial[(cata_fiducial[bin_col]>min_bin) & (cata_fiducial[bin_col]<=max_bin)].copy()
        cata_selec.reset_index(drop=True, inplace=True)
        cata_tmp = cata_test[(cata_test[bin_col]>min_bin) & (cata_test[bin_col]<=max_bin)].copy()
        cata_tmp.reset_index(drop=True, inplace=True)

        # basic info
        N_s = len(cata_selec)
        mean_bin = np.average(cata_selec[bin_col].values, weights=cata_selec['shape_weight'].values)

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


# >>>> working on sizeU5
# number fiducial ori 4188088
# number fiducial weight>0 3327371
# number test ori 4177685
# number test weight>0 3332170
# >>> number of points for lsq 40
# >>> number of points for lsq 40
# >>> number of points for lsq 40
# >>> number of points for lsq 40
# >>> number of points for lsq 40
# >>> number of points for lsq 40
# results saved as ./outputs/dm_sizeU5_magauto.csv
# >>>> working on sizeU10
# number fiducial ori 4188088
# number fiducial weight>0 3327371
# number test ori 4159632
# number test weight>0 3330305
# >>> number of points for lsq 40
# >>> number of points for lsq 40
# >>> number of points for lsq 40
# >>> number of points for lsq 40
# >>> number of points for lsq 40
# >>> number of points for lsq 40
# results saved as ./outputs/dm_sizeU10_magauto.csv
# >>>> working on sizeU15
# number fiducial ori 4188088
# number fiducial weight>0 3327371
# number test ori 4134490
# number test weight>0 3321653
# >>> number of points for lsq 40
# >>> number of points for lsq 40
# >>> number of points for lsq 40
# >>> number of points for lsq 40
# >>> number of points for lsq 40
# >>> number of points for lsq 40
# results saved as ./outputs/dm_sizeU15_magauto.csv
# >>>> working on qU5
# number fiducial ori 4188088
# number fiducial weight>0 3327371
# number test ori 4189906
# number test weight>0 3326854
# >>> number of points for lsq 40
# >>> number of points for lsq 40
# >>> number of points for lsq 40
# >>> number of points for lsq 40
# >>> number of points for lsq 40
# >>> number of points for lsq 40
# results saved as ./outputs/dm_qU5_magauto.csv
# >>>> working on qU10
# number fiducial ori 4188088
# number fiducial weight>0 3327371
# number test ori 4190250
# number test weight>0 3325811
# >>> number of points for lsq 40
# >>> number of points for lsq 40
# >>> number of points for lsq 40
# >>> number of points for lsq 40
# >>> number of points for lsq 40
# >>> number of points for lsq 40
# results saved as ./outputs/dm_qU10_magauto.csv
# >>>> working on qU15
# number fiducial ori 4188088
# number fiducial weight>0 3327371
# number test ori 4189703
# number test weight>0 3324558
# >>> number of points for lsq 40
# >>> number of points for lsq 40
# >>> number of points for lsq 40
# >>> number of points for lsq 40
# >>> number of points for lsq 40
# >>> number of points for lsq 40
# results saved as ./outputs/dm_qU15_magauto.csv
# >>>> working on nU10
# number fiducial ori 4188088
# number fiducial weight>0 3327371
# number test ori 4180728
# number test weight>0 3315180
# >>> number of points for lsq 40
# >>> number of points for lsq 40
# >>> number of points for lsq 40
# >>> number of points for lsq 40
# >>> number of points for lsq 40
# >>> number of points for lsq 40
# results saved as ./outputs/dm_nU10_magauto.csv
# >>>> working on nU20
# number fiducial ori 4188088
# number fiducial weight>0 3327371
# number test ori 4170830
# number test weight>0 3301496
# >>> number of points for lsq 40
# >>> number of points for lsq 40
# >>> number of points for lsq 40
# >>> number of points for lsq 40
# >>> number of points for lsq 40
# >>> number of points for lsq 40
# results saved as ./outputs/dm_nU20_magauto.csv
# >>>> working on nU30
# number fiducial ori 4188088
# number fiducial weight>0 3327371
# number test ori 4160472
# number test weight>0 3287526
# >>> number of points for lsq 40
# >>> number of points for lsq 40
# >>> number of points for lsq 40
# >>> number of points for lsq 40
# >>> number of points for lsq 40
# >>> number of points for lsq 40
# results saved as ./outputs/dm_nU30_magauto.csv
# Elapsed:4:16.77,User=382.089,System=1501.684,CPU=733.6%.