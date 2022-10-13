# -*- coding: utf-8 -*-
# @Author: lshuns
# @Date:   2022-04-13 14:05:23
# @Last Modified by:   lshuns
# @Last Modified time: 2022-10-13 14:31:33

### functions for the bias calculation

__all__ = ['alphaCalFunc_least_squares', 
                'mCalFunc_pair_based', 'mCalFunc_tile_based', 
                'mCalFunc_DataRewei_2D']

import numpy as np
import pandas as pd
import statsmodels.api as sm 

def alphaCalFunc_least_squares(cataSim, psf_frame=False):
    """
    Calculate the alpha term for a given simulated catalogue
        by requiring sum(e_out - alpha*e_psf) -> min

        Used columns and their names:
            e1_out, e2_out: measured ellipticity
            e1_psf, e2_psf: measured PSF ellipticity
            shape_weight: shape measurement weights
    """

    # out shear
    e1_out = np.array(cataSim['e1_out'])
    e2_out = np.array(cataSim['e2_out'])
    wgSim = np.array(cataSim['shape_weight'])

    # out PSF 
    e1_psf = np.array(cataSim['e1_psf'])
    e2_psf = np.array(cataSim['e2_psf'])

    # rotate to psf frame
    if psf_frame:
        print('Measured e will be rotated to PSF frame...')
        PSF_angle = np.arctan2(e2_psf, e1_psf)

        # rotate the ellipticities
        ct = np.cos(PSF_angle)
        st = np.sin(PSF_angle)
        e1_uc = e1_out*ct + e2_out*st
        e2_uc = -e1_out*st + e2_out*ct

        e1_out = e1_uc
        e2_out = e2_uc
        del PSF_angle, ct, st, e1_uc, e2_uc

        # e1 psf is amplitude
        e1_psf = np.hypot(e1_psf, e2_psf)
        # e2 psf is nothing
        e2_psf = np.zeros(len(e1_psf))
    del cataSim

    # get least square values
    ## e1
    mod_wls = sm.WLS(e1_out, sm.add_constant(e1_psf), weights=wgSim)
    res_wls = mod_wls.fit()
    alpha1 = res_wls.params[1]
    c1 = res_wls.params[0]
    alpha1_err = (res_wls.cov_params()[1, 1])**0.5
    c1_err = (res_wls.cov_params()[0, 0])**0.5
    ## e2
    mod_wls = sm.WLS(e2_out, sm.add_constant(e2_psf), weights=wgSim)
    res_wls = mod_wls.fit()
    alpha2 = res_wls.params[1]
    c2 = res_wls.params[0]
    alpha2_err = (res_wls.cov_params()[1, 1])**0.5
    c2_err = (res_wls.cov_params()[0, 0])**0.5

    # save
    res = {'alpha1': alpha1, 'alpha2': alpha2,
            'c1': c1, 'c2': c2,
            'alpha1_err': alpha1_err, 'alpha2_err': alpha2_err,
            'c1_err': c1_err, 'c2_err': c2_err,
            }

    return res

def mCalFunc_pair_based(cataSim, psf_frame=False):
    """
    Calculate the residual shear bias for a given simulated catalogue
        it first calculates average e for each pair of objects with different rot but same g_in (accounting for shape noise cancellation)
            then estimate by bias requiring sum(g_out - (1+m)g_in) -> min

        Used columns and their names:
            id_input: identify unique input galaxies
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
    cataSim = pd.DataFrame({'id_input': np.array(cataSim['id_input']),
                            'g1_in': np.array(cataSim['g1_in']),
                            'g2_in': np.array(cataSim['g2_in']),
                            'e1_out': e1_out,
                            'e2_out': e2_out,
                            'shape_weight': np.array(cataSim['shape_weight'])
                            })
    del e1_out, e2_out

    # sort to speed up
    cataSim.sort_values(by=['id_input', 'g1_in', 'g2_in'], inplace=True)

    # prepare the weighted mean
    cataSim.loc[:, 'e1_out'] *= cataSim['shape_weight'].values
    cataSim.loc[:, 'e2_out'] *= cataSim['shape_weight'].values

    # group based on input id and shear
    cataSim = cataSim.groupby(['id_input', 'g1_in', 'g2_in'], as_index=False).sum()
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

##################### following functions are for data reweighting 
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

def mCalFunc_DataRewei_2D(cataSim_tmp, cataReal, 
            bin1_col = 'SNR', bin1_Nbins = 20,
            bin2_col = 'R', bin2_Nbins = 20, 
            fitting_method='tile_based',
            psf_frame = False,
            save_surface_path=None,
            save_bounds_prefix=None):
    """
    Calculating the shear bias with data reweighting.
        Reweighting performed in 2 parameters.
        Used columns and their names:
            Simulation:
                two reweighting columns + columns used by the chosen fitting_method
            Data:
                two reweighting columns + weight
    """

    # >>> 0. build dataframe and select used columns
    ## the simulation
    cataSim = pd.DataFrame({bin1_col: np.array(cataSim_tmp[bin1_col]),
                            bin2_col: np.array(cataSim_tmp[bin2_col]),
                            'bin1_id': -999,
                            'bin2_id': -999,
                            'g1_in': np.array(cataSim_tmp['g1_in']),
                            'g2_in': np.array(cataSim_tmp['g2_in']),
                            'e1_out': np.array(cataSim_tmp['e1_out']),
                            'e2_out': np.array(cataSim_tmp['e2_out']),
                            'shape_weight': np.array(cataSim_tmp['shape_weight'])
                            })
    if fitting_method == 'tile_based':
        cataSim.loc[:, 'tile_label'] = np.array(cataSim_tmp['tile_label'])
    elif fitting_method == 'pair_based':
        cataSim.loc[:, 'id_input'] = np.array(cataSim_tmp['id_input'])
    else:
        raise Exception(f'Unsupported fitting_method {fitting_method}!')
    del cataSim_tmp
    ## the data
    cataReal = pd.DataFrame({bin1_col: np.array(cataReal[bin1_col]),
                            bin2_col: np.array(cataReal[bin2_col]),
                            'bin1_id': -999,
                            'bin2_id': -999,
                            'shape_weight': np.array(cataReal['shape_weight'])
                            })

    # >>> 1. binning
    ## each bin contains equal sum of weights for simulation
    bin1_bounds, bin2_bounds = _WgBin2DFunc(cataSim[bin1_col].values, 
                                            cataSim[bin2_col].values, 
                                            cataSim['shape_weight'].values, 
                                            bin1_Nbins, bin2_Nbins,
                                            right=True)
    if save_bounds_prefix is not None:
        outpath_tmp = save_bounds_prefix + '_bin1.npy'
        np.save(outpath_tmp, bin1_bounds)
        print('bin1_bounds saved to', outpath_tmp)
        outpath_tmp = save_bounds_prefix + '_bin2.npy'
        np.save(outpath_tmp, bin2_bounds)
        print('bin1_bounds saved to', outpath_tmp)
        del outpath_tmp

    ## bin in the first col
    cataSim.loc[:, 'bin1_id'] = pd.cut(cataSim[bin1_col].values, bin1_bounds, 
                                right=True, labels=False)
    cataReal.loc[:, 'bin1_id'] = pd.cut(cataReal[bin1_col].values, bin1_bounds, 
                                right=True, labels=False)
    del bin1_bounds

    ## bin in the second col
    for ibin1, bin2_bound in enumerate(bin2_bounds):
        # mask in bin1
        mask_bin1_sim = cataSim['bin1_id'].values == ibin1
        mask_bin1_real = cataReal['bin1_id'].values == ibin1
        # bin in bin2
        cataSim.loc[mask_bin1_sim, 'bin2_id'] = \
                            pd.cut(cataSim.loc[mask_bin1_sim, bin2_col].values, 
                                        bin2_bound, 
                                        right=True, labels=False)
        cataReal.loc[mask_bin1_real, 'bin2_id'] = \
                            pd.cut(cataReal.loc[mask_bin1_real, bin2_col].values,
                                        bin2_bound, 
                                        right=True, labels=False)
        del mask_bin1_real, mask_bin1_sim

    # >>> 2. calculate mc in each bin
    ## group
    cataSim.fillna(-999, inplace=True)
    cataSim = cataSim.astype({'bin1_id': int, 'bin2_id': int})
    cataSim = cataSim.groupby(by=['bin1_id', 'bin2_id'])

    cataReal.fillna(-999, inplace=True)
    cataReal = cataReal.astype({'bin1_id': int, 'bin2_id': int})
    ### for data, we only need the weights
    cataReal = cataReal[['bin1_id', 'bin2_id', 'shape_weight']].groupby(by=['bin1_id', 'bin2_id'])

    ## loop over bins and calculate m
    wgRealSum = 0 
    res = {'m1': 0, 'm2': 0,
        'c1': 0, 'c2': 0,
        'm1_err': 0, 'm2_err': 0,
        'c1_err': 0, 'c2_err': 0}
    if save_surface_path is not None:
        i_group = 0 
        mc_surface = pd.DataFrame(-999, 
                            index = np.arange(bin1_Nbins*bin2_Nbins), 
                            columns = ['bin1_id', 'bin2_id', 
                                        'bin1_mean', 'bin2_mean',
                                        'm1', 'm1_err', 'm2', 'm2_err',
                                        'c1', 'c1_err', 'c2', 'c2_err',
                                        'dataWei', 'simWei'])
    for i in range(bin1_Nbins):
        for j in range(bin2_Nbins):

            # the name for group 
            name = (i, j)
            group = cataSim.get_group(name)

            # weights from data for re-weighting
            wgRealBin = np.sum(cataReal.get_group(name)['shape_weight'].values)
            wgRealSum += wgRealBin

            # mc fitting
            if fitting_method == 'tile_based':
                res_bin = mCalFunc_tile_based(group, psf_frame=psf_frame)
            elif fitting_method == 'pair_based':
                res_bin = mCalFunc_pair_based(group, psf_frame=psf_frame)
            else:
                raise Exception(f'Unsupported fitting_method {fitting_method}!')
            for key, value in res_bin.items():
                if 'err' in key:
                    res[key] += (wgRealBin * value)**2.
                else:
                    res[key] += (wgRealBin * value)

            # saving the surface
            if save_surface_path is not None:
                mc_surface.loc[i_group, 'bin1_id'] = i
                mc_surface.loc[i_group, 'bin2_id'] = j
                mc_surface.loc[i_group, 'bin1_mean'] = np.average(group[bin1_col].values, weights=group['shape_weight'].values)
                mc_surface.loc[i_group, 'bin2_mean'] = np.average(group[bin2_col].values, weights=group['shape_weight'].values)
                mc_surface.loc[i_group, 'm1'] = res_bin['m1']
                mc_surface.loc[i_group, 'm1_err'] = res_bin['m1_err']
                mc_surface.loc[i_group, 'm2'] = res_bin['m2']
                mc_surface.loc[i_group, 'm2_err'] = res_bin['m2_err']
                mc_surface.loc[i_group, 'c1'] = res_bin['c1']
                mc_surface.loc[i_group, 'c1_err'] = res_bin['c1_err']
                mc_surface.loc[i_group, 'c2'] = res_bin['c2']
                mc_surface.loc[i_group, 'c2_err'] = res_bin['c2_err']
                mc_surface.loc[i_group, 'dataWei'] = wgRealBin
                mc_surface.loc[i_group, 'simWei'] = np.sum(group['shape_weight'].values)
                i_group += 1 
            del group

    del cataSim, cataReal

    ## save surface
    if save_surface_path is not None:
        mc_surface.to_csv(save_surface_path, index=False)
        print('the bias surface saved as', save_surface_path)

    ## take the average
    for key, value in res.items():
        if 'err' in key:
            res[key] = value**0.5 / wgRealSum
        else:
            res[key] = value / wgRealSum

    return res