# -*- coding: utf-8 -*-
# @Author: lshuns
# @Date:   2023-05-03 13:33:54
# @Last Modified by:   lshuns
# @Last Modified time: 2023-05-03 13:43:44

# the function used to calculate m for any given sample by interpreting the surface of doom

import numpy as np 
import pandas as pd 

def mCalFunc_from_doom(cata, m_doom, 
                        col_ZB, col_SNR, col_R, col_weight,
                        col_m1, col_m2):
    """
    Calculate the shear bias for a given catalogue
        It first splits galaxies into bins based on their redshift, SNR, R
            then assign each galaxy a m using the surface of doom
        The final results are the weighted average of individual m

    Parameters:
    -----------
    cata : catalogue of the data for which the mean m is calculated

    m_doom : catalogue of the surface of doom with shear bias and binning info

    col_ZB : column name for the ZB redshift

    col_SNR : column name for the signal-to-noise ratio

    col_R : column name for the resolution factor

    col_weight: column name for the measurement weight

    col_m1 : column name for m1 in the surface of doom

    col_m2 : column name for m2 in the surface of doom
    """

    # used columns from data cata
    cata = pd.DataFrame({'Z_B': np.array(cata[col_ZB]).astype(float),
                        'SNR': np.array(cata[col_SNR]).astype(float),
                        'R': np.array(cata[col_R]).astype(float),
                        'weight': np.array(cata[col_weight]).astype(float),
                        'binZB_id': -999, 'binSNR_id': -999, 'binR_id': -999})

    # used columns from surface of doom
    m_doom = pd.DataFrame({'binZB_id': np.array(m_doom['binZB_id']).astype(int),
                        'binSNR_id': np.array(m_doom['binSNR_id']).astype(int),
                        'binR_id': np.array(m_doom['binR_id']).astype(int),
                        'binZB_min': np.array(m_doom['binZB_min']).astype(float),
                        'binZB_max': np.array(m_doom['binZB_max']).astype(float),
                        'binSNR_min': np.array(m_doom['binSNR_min']).astype(float),
                        'binSNR_max': np.array(m_doom['binSNR_max']).astype(float),
                        'binR_min': np.array(m_doom['binR_min']).astype(float),
                        'binR_max': np.array(m_doom['binR_max']).astype(float),
                        'm1': np.array(m_doom[col_m1]).astype(float),
                        'm1_err': np.array(m_doom[f'{col_m1}_err']).astype(float),
                        'm2': np.array(m_doom[col_m2]).astype(float),
                        'm2_err': np.array(m_doom[f'{col_m2}_err']).astype(float)
                        })

    # bin galaxies
    Z_B_edges = np.unique(m_doom['binZB_min']).tolist() + [np.max(m_doom['binZB_max'])]
    cata.loc[:, 'binZB_id'] = pd.cut(cata['Z_B'].values, Z_B_edges, 
                                right=True, labels=False)
    for i_zbin in range(len(Z_B_edges)-1):

        SNR_edges = np.unique(m_doom.loc[(m_doom['binZB_id']==i_zbin), 
                        'binSNR_min']).tolist() \
                        + [np.max(m_doom.loc[(m_doom['binZB_id']==i_zbin), 
                            'binSNR_max'])]

        mask_binZB = cata['binZB_id'].values == i_zbin
        cata.loc[mask_binZB, 'binSNR_id'] = pd.cut(cata.loc[mask_binZB, 'SNR'].values, SNR_edges, 
                                    right=True, labels=False)

        for i_SNR in range(len(SNR_edges)-1):

            R_edges = np.unique(m_doom.loc[(m_doom['binZB_id']==i_zbin)&(m_doom['binSNR_id']==i_SNR), 
                        'binR_min']).tolist() \
                        + [np.max(m_doom.loc[(m_doom['binZB_id']==i_zbin)&(m_doom['binSNR_id']==i_SNR), 
                            'binR_max'])]

            mask_binSNR = cata['binSNR_id'].values == i_SNR
            cata.loc[mask_binZB&mask_binSNR, 'binR_id'] = pd.cut(
                                        cata.loc[mask_binZB&mask_binSNR, 'R'].values, 
                                        R_edges, 
                                        right=True, labels=False)
            del mask_binSNR
        del mask_binZB
    # print(">>>> sorted bin ZB", np.sort(cata['binZB_id'].values))
    # print(">>>> sorted bin SNR", np.sort(cata['binSNR_id'].values))
    # print(">>>> sorted bin R", np.sort(cata['binR_id'].values))

    # group
    ## sort to speed up
    cata = cata.astype({'binZB_id': int, 'binSNR_id': int, 'binR_id': int})
    cata.sort_values(by=['binZB_id', 'binSNR_id', 'binR_id'], inplace=True)
    cata = cata.groupby(by=['binZB_id', 'binSNR_id', 'binR_id'])

    # loop over groups to get mean m
    m1_final = 0 
    m1_err_final = 0
    m2_final = 0 
    m2_err_final = 0
    wgRealSum = 0 
    for name, group in cata:
        binZB_id, binSNR_id, binR_id = name

        # total weights in each bin
        wgRealBin = np.sum(group['weight'].values)
        wgRealSum += wgRealBin

        # m from surface of doom
        mask_doom = (m_doom['binZB_id']==binZB_id)\
                &(m_doom['binSNR_id']==binSNR_id)\
                &(m_doom['binR_id']==binR_id)
        m1_final += (wgRealBin * m_doom.loc[mask_doom, 'm1'].values[0])
        m2_final += (wgRealBin * m_doom.loc[mask_doom, 'm2'].values[0])
        m1_err_final += (wgRealBin * m_doom.loc[mask_doom, 'm1_err'].values[0])**2
        m2_err_final += (wgRealBin * m_doom.loc[mask_doom, 'm2_err'].values[0])**2

    # take the mean
    m1_final = m1_final / wgRealSum
    m2_final = m2_final / wgRealSum
    m1_err_final = m1_err_final**0.5 / wgRealSum
    m2_err_final = m2_err_final**0.5 / wgRealSum
    return (m1_final, m2_final, m1_err_final, m2_err_final)
