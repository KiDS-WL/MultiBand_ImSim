# -*- coding: utf-8 -*-
# @Author: lshuns
# @Date:   2023-04-03 09:27:37
# @Last Modified by:   lshuns
# @Last Modified time: 2024-02-26 16:53:22

### calculate the dm between variations and fiducial
###### in redshift bins
###### reweighting with data

import os
import numpy as np
import pandas as pd
import statsmodels.api as sm 

from scipy import stats
from astropy.io import fits
from astropy.table import Table

import plotting

############### I/O

# ### weight using the gold sample
# inpath_data = '/disks/shear10/ssli/K1000CS/LF321_Outputs/K1000_NS_V1.0.0A_ugriZYJHKs_photoz_SG_mask_LF_glab_321_v2_A1_goldclasses.cat.A2.feather'
# wei_col_data = 'AlphaRecalC_weight'
# save_suffix = 'gold_reweighted'

### weight using the whole sample
inpath_data = '/disks/shear10/ssli/KiDS/K1000_LF_321_mosaic/KIDS_1006tiles_r_SDSS.V1.0.0A_ugriZYJHKs_photoz_SG_mask_LF_glab_321_v2_selec_noWeiCut.feather'
wei_col_data = 'weight'
save_suffix = 'nogold_reweighted'

# the test tags
test_labels = ['skills_v07D7_nD', 'skills_v07D7_nU']
# test_labels = ['skills_v07D7_qD', 'skills_v07D7_qU']
# test_labels = ['skills_v07D7_sizeD', 'skills_v07D7_sizeU']

# the test catalogues
inpath_test_list = [f'/disks/shear16/ssli/ImSim/output/{test_label}/{test_label}_LF_321_part0_everything_col_flag.feather'
                            for test_label in test_labels]
# inpath_test_list = [f'/disks/shear16/ssli/ImSim/output/{test_label}/{test_label}_LF_321_part0_no_134.0_0.5_everything_col_flag.feather'
#                             for test_label in test_labels]

# the fiducial catalogue
inpath_fiducial = '/disks/shear16/ssli/ImSim/output/skills_v07D7/skills_v07D7p1_LF_321_part0_everything_col_flag.feather'
# inpath_fiducial = '/disks/shear16/ssli/ImSim/output/skills_v07D7/skills_v07D7p1_LF_321_part0_no_134.0_0.5_everything_col_flag.feather'

# where to save results
outDir = './results'

# binning info
bin_col = 'Z_B'
# bin_edges = [0.1, 0.3, 0.5, 0.7, 0.9, 1.2]
bin_edges = [0.1, 0.3, 0.5, 0.7, 0.9, 1.2, 2.0]

############### function
def _WgQuantile1DFunc(values, weights, Nbin):
    """
    Calculate the weighted quantile by given bin numbers
        designed for 1D numpy array.
    """

    # Define the quantile points based on the number of bins
    pq = np.linspace(0, 1, Nbin+1)

    # Sort the data
    ind_sorted = np.argsort(values)
    v_sorted = values[ind_sorted]
    wg_sorted = weights[ind_sorted]
    del ind_sorted, values, weights

    # Compute the quantiles
    Pn = (np.cumsum(wg_sorted) - 0.5*wg_sorted)/np.sum(wg_sorted)

    # interp the quantiles
    res = np.interp(pq, Pn, v_sorted)

    # include the edge points
    if v_sorted[0] > 0:
        res[0] = 0.99*v_sorted[0]
    else:
        res[0] = 1.01*v_sorted[0]
    if v_sorted[-1] > 0:
        res[-1] = 1.01*v_sorted[-1]
    else:
        res[-1] = 0.99*v_sorted[-1]

    return res

def _WgBin2DFunc(v1, v2, wgs, Nbin1, Nbin2, right=True):
    """
    Calculate the weighted quantile by given bin numbers
        designed for 2D numpy array
    """

    # Calculate quantiles for v1
    q1 = _WgQuantile1DFunc(v1, wgs, Nbin1)

    #Compute quantiles for v2 in each v1 bin
    q2s = np.zeros((Nbin1, Nbin2+1))
    for i in range(Nbin1):

        if right:
            mask = (v1>q1[i]) & (v1<=q1[i+1])
        else:
            mask = (v1>=q1[i]) & (v1<q1[i+1])

        q2s[i] = _WgQuantile1DFunc(v2[mask], wgs[mask], Nbin2)

    return q1, q2s

def Reweight(cataSim, cataData, 
                bin1_col, bin2_col, bin1_Nbins, bin2_Nbins):

    cataSim = cataSim.copy()
    cataData = cataData.copy()

    # 1. binning
    ## each bin contains equal sum of weights for simulation
    bin1_bounds, bin2_bounds = _WgBin2DFunc(cataSim[bin1_col].values, 
                                            cataSim[bin2_col].values, 
                                            cataSim['shape_weight'].values, 
                                            bin1_Nbins, bin2_Nbins,
                                            right=True)
    ## bin in the first col
    cataSim.loc[:, 'bin1_id'] = pd.cut(cataSim[bin1_col].values, bin1_bounds, 
                                right=True, labels=False)
    cataData.loc[:, 'bin1_id'] = pd.cut(cataData[bin1_col].values, bin1_bounds, 
                                right=True, labels=False)
    del bin1_bounds
    ## bin in the second col
    for ibin1, bin2_bound in enumerate(bin2_bounds):
        # mask in bin1
        mask_bin1_sim = cataSim['bin1_id'].values == ibin1
        mask_bin1_real = cataData['bin1_id'].values == ibin1
        # bin in bin2
        cataSim.loc[mask_bin1_sim, 'bin2_id'] = \
                            pd.cut(cataSim.loc[mask_bin1_sim, bin2_col].values, 
                                        bin2_bound, 
                                        right=True, labels=False)
        cataData.loc[mask_bin1_real, 'bin2_id'] = \
                            pd.cut(cataData.loc[mask_bin1_real, bin2_col].values,
                                        bin2_bound, 
                                        right=True, labels=False)
        del mask_bin1_real, mask_bin1_sim
    ## drop those outside the bins
    if cataSim.isnull().values.any():
        raise Exception('cataSim has nan, check what happened!')
    cataData = cataData[(cataData['bin1_id']>=0)&(cataData['bin1_id']<bin1_Nbins)&
                            (cataData['bin2_id']>=0)&(cataData['bin2_id']<bin2_Nbins)]
    cataSim = cataSim.astype({'bin1_id': int, 'bin2_id': int})
    cataData = cataData.astype({'bin1_id': int, 'bin2_id': int})

    # # >>> test plot: para for data
    # outpath = 'show'
    # paras = [np.array(cataData[bin1_col])]
    # wgs = [np.array(cataData['shape_weight'])]
    # COLORs = ['k']
    # LABELs = ['data']
    # LINEs = ['-']
    # nbins = 60
    # XRANGE = (np.min(cataData[bin1_col]), np.max(cataData[bin1_col]))
    # XLABEL = bin1_col
    # YLABEL = 'weighted DENSITY'
    # DENSITY = True
    # xlog = True

    # 2. get scaling values
    wei_sum = np.sum(cataSim['shape_weight'].values)
    cataSim_tmp = cataSim[['bin1_id', 'bin2_id', 'shape_weight']].sort_values(by=['bin1_id', 'bin2_id'])
    weight_sim = cataSim_tmp.groupby(by=['bin1_id', 'bin2_id']).sum()/wei_sum
    del cataSim_tmp
    # print('>>>> weight_sim', weight_sim)

    wei_sum = np.sum(cataData['shape_weight'].values)
    cataData = cataData[['bin1_id', 'bin2_id', 'shape_weight']].sort_values(by=['bin1_id', 'bin2_id'])
    weight_data = cataData.groupby(by=['bin1_id', 'bin2_id']).sum()/wei_sum
    # del cataData
    # print('>>>> weight_data', weight_data)

    scaling_df = weight_data / weight_sim
    del weight_data, weight_sim
    # print('>>>> scaling_df', scaling_df)
    # print('>>>>>>>>> mean', np.mean(scaling_df))

    # # >>> test plot: para before reweighting
    # outpath = 'show'
    # paras.append(np.array(cataSim[bin1_col]))
    # wgs.append(np.array(cataSim['shape_weight']))
    # COLORs.append('b')
    # LABELs.append('before')
    # LINEs.append('-')

    # 3. reweighting 
    for i in range(bin1_Nbins):
        for j in range(bin2_Nbins):
            cataSim.loc[(cataSim['bin1_id'].values == i) & (cataSim['bin2_id'].values == j), 'shape_weight'] *= scaling_df.loc[(i, j), 'shape_weight']
    del scaling_df

    # # >>> test plot: make the plot
    # paras.append(np.array(cataSim[bin1_col]))
    # wgs.append(np.array(cataSim['shape_weight']))
    # COLORs.append('orange')
    # LABELs.append('after')
    # LINEs.append('--')
    # plotting.HistPlotFunc(outpath,
    #             paras, wgs, COLORs, LABELs,
    #             nbins, XRANGE, YRANGE=None,
    #             XLABEL=XLABEL, YLABEL=YLABEL,
    #             DENSITY=DENSITY, HISTTYPE='step', STACKED=False,
    #             TITLE=None, xtick_min_label=True, ytick_min_label=True,
    #             xtick_spe=None, ytick_spe=None,
    #             vlines=None, vline_styles=None, vline_colors=None, vline_labels=None, vline_widths=None,
    #             hlines=None, hline_styles=None, hline_colors=None, hline_labels=None, hline_widths=None,
    #             xlog=xlog, ylog=False,
    #             loc_legend='best', 
    #             font_size=16, usetex=False,
    #             cumulative=False,
    #             LINEs=LINEs)

    return cataSim

def dmCalFunc_tile_based_dataRewei(cataSim0, cataSim1, cataData, 
                    bin1_col = 'SNR', bin1_Nbins = 20,
                    bin2_col = 'R', bin2_Nbins = 20, 
                    nboot=500):
    """
    assume cataSim0 and cataSim1 having same tiles
    dm = cataSim0 - cataSim1 (tile-by-tile)
    reweight using cataData
    """

    # >>>> reweight the simulations
    # build dataframe and select used columns
    cataSim0 = pd.DataFrame({bin1_col: np.array(cataSim0[bin1_col]),
                            bin2_col: np.array(cataSim0[bin2_col]),
                            'bin1_id': -999,
                            'bin2_id': -999,
                            'g1_in': np.array(cataSim0['g1_in']),
                            'g2_in': np.array(cataSim0['g2_in']),
                            'e1_out': np.array(cataSim0['e1_out']),
                            'e2_out': np.array(cataSim0['e2_out']),
                            'shape_weight': np.array(cataSim0['shape_weight']),
                            'tile_label': np.array(cataSim0['tile_label'])
                            })
    cataSim1 = pd.DataFrame({bin1_col: np.array(cataSim1[bin1_col]),
                            bin2_col: np.array(cataSim1[bin2_col]),
                            'bin1_id': -999,
                            'bin2_id': -999,
                            'g1_in': np.array(cataSim1['g1_in']),
                            'g2_in': np.array(cataSim1['g2_in']),
                            'e1_out': np.array(cataSim1['e1_out']),
                            'e2_out': np.array(cataSim1['e2_out']),
                            'shape_weight': np.array(cataSim1['shape_weight']),
                            'tile_label': np.array(cataSim1['tile_label'])
                            })
    cataData = pd.DataFrame({bin1_col: np.array(cataData[bin1_col]),
                            bin2_col: np.array(cataData[bin2_col]),
                            'bin1_id': -999,
                            'bin2_id': -999,
                            'shape_weight': np.array(cataData['shape_weight'])
                            })
    # reweighting
    print(">>>> Reweight cataSim0")
    cataSim0 = Reweight(cataSim0, cataData, 
                bin1_col, bin2_col, bin1_Nbins, bin2_Nbins)
    print(">>>> Reweight cataSim1")
    cataSim1 = Reweight(cataSim1, cataData, 
                bin1_col, bin2_col, bin1_Nbins, bin2_Nbins)
    del cataData

    # >>>> average the first catalogue
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
    if not cataSim0[['tile_label', 'g1_in', 'g2_in']].equals(cataSim1[['tile_label', 'g1_in', 'g2_in']]):
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

# the data catalogue
cata_data = pd.read_feather(inpath_data)
print('number data ori', len(cata_data))
## weight cut
cata_data = cata_data[cata_data[wei_col_data].values>0]
cata_data.reset_index(drop=True, inplace=True)
print('number data weight>0', len(cata_data))
## used columns
cata_data = cata_data[[bin_col,
                        wei_col_data, 'model_SNratio', 'R']]
## rename
cata_data.rename(columns={wei_col_data: 'shape_weight',
                    'model_SNratio': 'SNR'}, inplace=True)

# the fiducial catalogue
cata_fiducial = pd.read_feather(inpath_fiducial)
print('number fiducial ori', len(cata_fiducial))
## LF shear selections
cata_fiducial = cata_fiducial[(cata_fiducial['flag_gaap']==0)&(cata_fiducial['flag_asteroid']==0)&(cata_fiducial['flag_binary']==0)&(cata_fiducial['flag_LF_noWeiCut']==0)]
cata_fiducial.reset_index(drop=True, inplace=True)
print('number after selection', len(cata_fiducial))
## weight cut
cata_fiducial = cata_fiducial[cata_fiducial['oldweight_LF_r'].values>0]
cata_fiducial.reset_index(drop=True, inplace=True)
print('number weight>0', len(cata_fiducial))
## used columns
cata_fiducial = cata_fiducial[[bin_col,
                            'oldweight_LF_r', 'e1_LF_r', 'e2_LF_r', 
                            'g1_in', 'g2_in', 'tile_label', 'gal_rot',
                            'SNR_LF_r', 'R']]
## rename
cata_fiducial.rename(columns={'oldweight_LF_r': 'shape_weight',
                    'e1_LF_r': 'e1_out', 'e2_LF_r': 'e2_out',
                    'SNR_LF_r': 'SNR'}, inplace=True)

for i_cata, inpath in enumerate(inpath_test_list):

    test_label = test_labels[i_cata]
    print('>>>> working on', test_label)

    # load the catalogue
    cata_test = pd.read_feather(inpath)
    print('number test ori', len(cata_test))
    ## LF shear selections
    cata_test = cata_test[(cata_test['flag_gaap']==0)&(cata_test['flag_asteroid']==0)&(cata_test['flag_binary']==0)&(cata_test['flag_LF_noWeiCut']==0)]
    cata_test.reset_index(drop=True, inplace=True)
    print('number after selection', len(cata_test))
    ## weight cut
    cata_test = cata_test[cata_test['oldweight_LF_r'].values>0]
    cata_test.reset_index(drop=True, inplace=True)
    print('number weight>0', len(cata_test))
    ## used columns
    cata_test = cata_test[[bin_col,
                                'oldweight_LF_r', 'e1_LF_r', 'e2_LF_r', 
                                'g1_in', 'g2_in', 'tile_label', 'gal_rot',
                                'SNR_LF_r', 'R']]
    ## rename
    cata_test.rename(columns={'oldweight_LF_r': 'shape_weight',
                        'e1_LF_r': 'e1_out', 'e2_LF_r': 'e2_out',
                        'SNR_LF_r': 'SNR'}, inplace=True)

    # >>>>>>>>> output
    outpath = os.path.join(outDir, f'dm_ZBbins_{test_label}_{save_suffix}.csv')
    f = open(outpath, 'w')

    # >>>>>>>>> m results
    ## ++++ the whole results
    res_residual = dmCalFunc_tile_based_dataRewei(cata_test, cata_fiducial, cata_data, 
                                bin1_col = 'SNR', bin1_Nbins = 20,
                                bin2_col = 'R', bin2_Nbins = 20, 
                                nboot=500)
    # basic info
    N_s = len(cata_fiducial)
    # collect columns names and values
    cols = ','.join(list(res_residual.keys()))
    vals = ','.join(["{0:0.4f}".format(val) for val in res_residual.values()])
    cols = cols + f',{bin_col}_min,{bin_col}_max,Nobj'
    vals = vals + f',-999,-999,{N_s}'    
    print(cols, file=f)
    print(vals, file=f)

    ## +++ the binning results
    for i_bin in range(len(bin_edges)-1):
        min_bin = bin_edges[i_bin]
        max_bin = bin_edges[i_bin+1]

        print('++++ for ZB bin', i_bin+1)

        # select cata
        cata_selec = cata_fiducial[(cata_fiducial[bin_col]>min_bin) & (cata_fiducial[bin_col]<=max_bin)].copy()
        cata_selec.reset_index(drop=True, inplace=True)
        cata_tmp = cata_test[(cata_test[bin_col]>min_bin) & (cata_test[bin_col]<=max_bin)].copy()
        cata_tmp.reset_index(drop=True, inplace=True)
        cata_data_selec = cata_data[(cata_data[bin_col]>min_bin) & (cata_data[bin_col]<=max_bin)].copy()
        cata_data_selec.reset_index(drop=True, inplace=True)

        # basic info
        N_s = len(cata_selec)

        # results
        res_residual = dmCalFunc_tile_based_dataRewei(cata_tmp, cata_selec, cata_data_selec, 
                                bin1_col = 'SNR', bin1_Nbins = 20,
                                bin2_col = 'R', bin2_Nbins = 20, 
                                nboot=500)
        del cata_selec, cata_tmp, cata_data_selec

        # collect columns names and values
        vals = ','.join(["{0:0.4f}".format(val) for val in res_residual.values()]) \
         + f',{min_bin},{max_bin},{N_s}'
        print(vals, file=f)

    # close what is opened
    f.close()
    print(f'results saved as {outpath}')
