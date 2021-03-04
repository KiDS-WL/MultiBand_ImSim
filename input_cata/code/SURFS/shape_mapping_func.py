# @Author: lshuns
# @Date:   2021-02-19, 14:04:29
# @Last modified by:   ssli
# @Last modified time: 2021-02-25, 18:46:29

### module containing functions to do the COSMOS-morphology mapping

import time
import logging
import random
random.seed(94120)

import numpy as np
import pandas as pd
from scipy import stats
import pyvinecopulib as pv
import matplotlib.pyplot as plt

from astropy.io import fits
from astropy.table import Table
from astropy.cosmology import Planck15 as astropy_Planck15

## logging
logging.basicConfig(format='%(name)s : %(levelname)s - %(message)s', level=logging.INFO)
logger = logging.getLogger(__name__)

__all__ = ["TwoBinningMappingFunc", "ZBinningMappingFunc", "SimpleMappingFunc"]


def _histedges_equalN(x, nbin):
    """
    estimate hist edges for even binning
    """
    npt = len(x)
    return np.interp(np.linspace(0, npt, nbin + 1),
                     np.arange(npt),
                     np.sort(x))

def TwoBinningMappingFunc(simu_ori_cata, cosmos_cata, N_zbins=20, dmagbins = 0.5, Niter = 10000):
    """
    Mapping with redshift and magnitude binning first
    """

    # ++++++++++++++++++++++ redshift binning
    ### adaptive binning
    _, edges_z_bins, _ = plt.hist(cosmos_cata['Z_optimal'], _histedges_equalN(cosmos_cata['Z_optimal'], N_zbins))
    edges_z_bins[0] = np.min(cosmos_cata['Z_optimal'])
    edges_z_bins[-1] = np.max(cosmos_cata['Z_optimal']) + 0.01 # for edge value
    # ### log binning
    # edges_z_bins = np.exp(np.linspace(np.log(0.001), np.log(z_max), N_zbins))
    # edges_z_bins[0] = z_min
    # edges_z_bins[-1] = z_max
    logger.debug(f'edges_z_bins {edges_z_bins}')

    # ++++++++++++++++++++++ magnitude binning
    time_start = time.time()

    edges_mag_bins = []
    for i_zbin in range(N_zbins):
        zmin = edges_z_bins[i_zbin]
        zmax = edges_z_bins[i_zbin+1]
        logger.debug(f'redshift bin: [{zmin}, {zmax})')
        cosmos_cata_selec0 = cosmos_cata[(cosmos_cata['Z_optimal']>=zmin) & (cosmos_cata['Z_optimal']<zmax)]
        logger.debug(f'Number of data {len(cosmos_cata_selec0)}')
        ## outliers are not taken into account when calculate the correlation
        mag_min = np.percentile(cosmos_cata_selec0['m_r'].values, 0.1)
        mag_max = np.percentile(cosmos_cata_selec0['m_r'].values,  99.9)
        ## but the final bounds still include all the objects
        abs_mag_r_min = np.min(cosmos_cata_selec0['m_r'].values)
        abs_mag_r_max = np.max(cosmos_cata_selec0['m_r'].values) + 0.01 # for edge value

        # iterating started
        a = mag_min
        b = a
        edges_mag_bins_tmp = [a]
        N_permagbin = []
        for Nrun in range(Niter):

            # selecting and calculation
            N_permagbin_tmp = 0
            ## avoid too small bins
            while N_permagbin_tmp <= 10:
                # break out
                if b >= mag_max:
                    break
                b += dmagbins
                cosmos_cata_selec = cosmos_cata_selec0[(cosmos_cata_selec0['m_r']>=a) & (cosmos_cata_selec0['m_r']<b)]
                N_permagbin_tmp = len(cosmos_cata_selec)

            rho_s_n, pval_s_n = stats.spearmanr(cosmos_cata_selec['m_r'].values, cosmos_cata_selec['N_GALFIT_HI'].values)
            rho_s_Re, pval_s_Re = stats.spearmanr(cosmos_cata_selec['m_r'].values, cosmos_cata_selec['Re_physical'].values)
            rho_s_e, pval_s_e = stats.spearmanr(cosmos_cata_selec['m_r'].values, cosmos_cata_selec['BA_GALFIT_HI'].values)
            ##
            pval_s = np.min([pval_s_n, pval_s_Re, pval_s_e])

            # keep cutting toward left until it meets requirement
            for Nrun_left in range(Niter):
                if (pval_s >= 0.05) or (N_permagbin_tmp<=50):
                    edges_mag_bins_tmp.append(b)
                    N_permagbin.append(N_permagbin_tmp)
                    break

                # split to half
                xmid = (b+a)/2.
                # selecting and calculation
                cosmos_cata_selec_tmp = cosmos_cata_selec[(cosmos_cata_selec['m_r']>=a) & (cosmos_cata_selec['m_r']<xmid)]
                N_permagbin_tmp = len(cosmos_cata_selec_tmp)
                rho_s_n, pval_s_n = stats.spearmanr(cosmos_cata_selec_tmp['m_r'].values, cosmos_cata_selec_tmp['N_GALFIT_HI'].values)
                rho_s_Re, pval_s_Re = stats.spearmanr(cosmos_cata_selec_tmp['m_r'].values, cosmos_cata_selec_tmp['Re_physical'].values)
                rho_s_e, pval_s_e = stats.spearmanr(cosmos_cata_selec_tmp['m_r'].values, cosmos_cata_selec_tmp['BA_GALFIT_HI'].values)
                ##
                pval_s = np.min([pval_s_n, pval_s_Re, pval_s_e])

                # iterate
                b = xmid

            # break out
            if b >= mag_max:
                ## using original boundaries
                edges_mag_bins_tmp[-1] = abs_mag_r_max
                edges_mag_bins_tmp[0] = abs_mag_r_min
                break

            # iterate
            a = b

        ## combine for each redshift bin
        edges_mag_bins.append(edges_mag_bins_tmp)
        logger.debug(f'edges_mag_bins {edges_mag_bins}')
        logger.debug(f'Number per bin {N_permagbin}')
    logger.debug(f'mag binning edges finding finished in {time.time()-time_start}')

    # ++++++++++++++++++++ Sampling & mapping
    time_start = time.time()
    for i_zbin in range(N_zbins):
        zmin = edges_z_bins[i_zbin]
        zmax = edges_z_bins[i_zbin+1]
        mask_z_cosmos = (cosmos_cata['Z_optimal']>=zmin) & (cosmos_cata['Z_optimal']<zmax)
        mask_z_simu = (simu_ori_cata['redshift']>=zmin) & (simu_ori_cata['redshift']<zmax)

        edges_mag_bins_tmp = edges_mag_bins[i_zbin]
        for i_magbin in range(len(edges_mag_bins_tmp)-1):
            mag_min = edges_mag_bins_tmp[i_magbin]
            mag_max = edges_mag_bins_tmp[i_magbin+1]
            mask_mag_cosmos = (cosmos_cata['m_r']>=mag_min) & (cosmos_cata['m_r']<mag_max)
            mask_mag_simu = (simu_ori_cata['mag_r_abs']>mag_min) & (simu_ori_cata['mag_r_abs']<mag_max)

            ## check if simulation has sources in this bin
            if np.sum(mask_z_simu & mask_mag_simu) == 0:
                logger.debug(f'bin z[{zmin}, {zmax}), mag[{mag_min}, {mag_max}) is empty')
                logger.debug('move on')
                continue

            cosmos_cata_selec = cosmos_cata[mask_z_cosmos & mask_mag_cosmos]
            simu_cata_selec = simu_ori_cata[mask_z_simu & mask_mag_simu]

            logger.debug(f'redshift bin: [{zmin}, {zmax}), mag_r bin [{mag_min}, {mag_max})')
            SimpleMappingFunc(simu_cata_selec, cosmos_cata_selec, simu_ori_cata, mask_z_simu, mask_mag_simu)

    return simu_ori_cata

def ZBinningMappingFunc(simu_ori_cata, cosmos_cata, N_zbins=20):
    """
    Mapping with redshift binning first
    """

    # ++++++++++++++++++++++ redshift binning
    ### adaptive binning
    _, edges_z_bins, _ = plt.hist(cosmos_cata['Z_optimal'], _histedges_equalN(cosmos_cata['Z_optimal'], N_zbins))
    edges_z_bins[0] = np.min(cosmos_cata['Z_optimal'])
    edges_z_bins[-1] = np.max(cosmos_cata['Z_optimal']) + 0.01 # for edge value
    # ### log binning
    # edges_z_bins = np.exp(np.linspace(np.log(0.001), np.log(z_max), N_zbins))
    # edges_z_bins[0] = z_min
    # edges_z_bins[-1] = z_max
    logger.debug(f'edges_z_bins {edges_z_bins}')

    # ++++++++++++++++++++ Sampling & mapping
    time_start = time.time()
    for i_zbin in range(N_zbins):
        zmin = edges_z_bins[i_zbin]
        zmax = edges_z_bins[i_zbin+1]
        mask_z_cosmos = (cosmos_cata['Z_optimal']>=zmin) & (cosmos_cata['Z_optimal']<zmax)
        mask_z_simu = (simu_ori_cata['redshift']>=zmin) & (simu_ori_cata['redshift']<zmax)

        ## check if simulation has sources in this bin
        if np.sum(mask_z_simu) == 0:
            logger.debug(f'bin z[{zmin}, {zmax}) is empty')
            logger.debug('move on')
            continue

        cosmos_cata_selec = cosmos_cata[mask_z_cosmos]
        simu_cata_selec = simu_ori_cata[mask_z_simu]

        SimpleMappingFunc(simu_cata_selec, cosmos_cata_selec, simu_ori_cata, mask_z_simu, mask_mag_simu=None)

    return simu_ori_cata

def SimpleMappingFunc(simu_cata_selec, cosmos_cata_selec, simu_ori_cata=None, mask_z_simu=None, mask_mag_simu=None):
    """
    Mapping without redshift and magnitude binning first
    """

    time_start = time.time()

    N_cosmos = len(cosmos_cata_selec)
    logger.debug(f'Number for data {N_cosmos}')
    N_simu = len(simu_cata_selec)
    logger.debug(f'Number for simu {N_simu}')

    # parameters
    ## cosmos
    Re_cosmos = cosmos_cata_selec['Re_physical'].values
    sersic_n_cosmos = cosmos_cata_selec['N_GALFIT_HI'].values
    BA_cosmos = cosmos_cata_selec['BA_GALFIT_HI'].values
    mag_r_cosmos = cosmos_cata_selec['r_mag_auto'].values
    mag_i_cosmos = cosmos_cata_selec['ip_mag_auto'].values
    ## simu
    size_simu = simu_cata_selec['size_physical'].values
    bf_simu = simu_cata_selec['bf'].values
    mag_r_simu = simu_cata_selec['mag_r'].values
    mag_i_simu = simu_cata_selec['mag_i'].values

    # ++++++++++++++ learning with vine-copulas
    # combine parameters for sampling
    X = np.empty((N_cosmos, 5))
    X[:, 0] = Re_cosmos
    X[:, 1] = sersic_n_cosmos
    X[:, 2] = BA_cosmos
    X[:, 3] = mag_r_cosmos
    X[:, 4] = mag_i_cosmos

    # Transform copula data using the empirical distribution
    u = pv.to_pseudo_obs(X)

    # create pv object directly from data
    cop = pv.Vinecop(data=u)

    # Sample from the copula
    u_sample = cop.simulate(N_simu, seeds=[1, 2, 3, 4, 5])

    # Transform back simulations to the original scale
    X_sample = np.asarray([np.quantile(X[:, i], u_sample[:, i]) for i in range(5)])
    ## transform back to parameters
    Re_sample = X_sample[0]
    sersic_n_sample = X_sample[1]
    BA_sample = X_sample[2]
    r_sample = X_sample[3]
    i_sample = X_sample[4]

    # ++++++++++++ mapping
    ## ordering according to size
    index_order_sample = np.argsort(Re_sample)
    Re_sample_order = Re_sample[index_order_sample]
    sersic_n_sample_order = sersic_n_sample[index_order_sample]
    BA_sample_order = BA_sample[index_order_sample]
    r_sample_order = r_sample[index_order_sample]
    i_sample_order = i_sample[index_order_sample]

    ## desired parameters
    ### 1. bf
    Re_simu_selec0_bf = np.empty(N_simu, dtype=float)
    BA_simu_selec0_bf = np.empty(N_simu, dtype=float)
    sersic_n_simu_selec0_bf = np.empty(N_simu, dtype=float)
    simu_combine_bf = [np.empty(N_simu, dtype=float), np.empty(N_simu, dtype=float), np.empty(N_simu, dtype=float)]

    ### 2. colour
    Re_simu_selec0_c = np.empty(N_simu, dtype=float)
    BA_simu_selec0_c = np.empty(N_simu, dtype=float)
    sersic_n_simu_selec0_c = np.empty(N_simu, dtype=float)
    simu_combine_c = [np.empty(N_simu, dtype=float), np.empty(N_simu, dtype=float), np.empty(N_simu, dtype=float)]

    ### 3. random
    Re_simu_selec0_r = np.empty(N_simu, dtype=float)
    BA_simu_selec0_r = np.empty(N_simu, dtype=float)
    sersic_n_simu_selec0_r = np.empty(N_simu, dtype=float)
    simu_combine_r = [np.empty(N_simu, dtype=float), np.empty(N_simu, dtype=float), np.empty(N_simu, dtype=float)]

    ## size binning
    Nsize = int(N_simu/500. + 1) # ~500 sources in each bin
    _, edges_size_bins, _ = plt.hist(size_simu, _histedges_equalN(size_simu, Nsize))
    edges_size_bins[-1] += 0.01 # for edge value
    Ncum = 0
    for i_sizebin in range(Nsize):
        size_min = edges_size_bins[i_sizebin]
        size_max = edges_size_bins[i_sizebin+1]
        mask_size = (size_simu >= size_min) & (size_simu < size_max)

        ## simulation
        bf_simu_selec = bf_simu[mask_size]
        r_simu_selec = mag_r_simu[mask_size]
        i_simu_selec = mag_i_simu[mask_size]
        Ns = len(bf_simu_selec)
        logger.debug(f'simu number in size ({size_min}, {size_max}]: {Ns}')

        ## cosmos
        Nend = Ncum + Ns
        Re_selec = Re_sample_order[Ncum:Nend]
        n_selec = sersic_n_sample_order[Ncum:Nend]
        BA_selec = BA_sample_order[Ncum:Nend]
        r_selec = r_sample_order[Ncum:Nend]
        i_selec = i_sample_order[Ncum:Nend]
        cosmos_combine_tmp = [Re_selec, BA_selec, n_selec]
        ### iterating
        Ncum += Ns

        ############################# 1. bf vs. sersic_n
        #############################    (coupled with colour ranking for bf==1)
        ## ordering
        index_order_bf = np.argsort(bf_simu_selec)
        index_order_n = np.argsort(n_selec)
        Nbf1 = np.sum(bf_simu_selec<1)
        index_ri_simu = np.argsort(r_simu_selec[index_order_bf[Nbf1:]] - i_simu_selec[index_order_bf[Nbf1:]])
        index_ri_cosmos = np.argsort(r_selec[index_order_n[Nbf1:]] - i_selec[index_order_n[Nbf1:]])

        ## assignment
        ## tmp variables for ordering assignment
        tmp1 = np.empty(Ns, dtype=float)
        tmp2 = np.empty(len(index_ri_simu), dtype=float)
        for i_refer, tmp_simu_refer in enumerate(simu_combine_bf):

            ### a: bf<1 using bf ranking
            tmp1[index_order_bf[:Nbf1]] = cosmos_combine_tmp[i_refer][index_order_n[:Nbf1]]

            ### b: bf==1 using colour ranking
            tmp2[index_ri_simu] = cosmos_combine_tmp[i_refer][index_order_n[Nbf1:]][index_ri_cosmos]
            tmp1[index_order_bf[Nbf1:]] = tmp2

            tmp_simu_refer[mask_size] = tmp1

        ############################# 2. r-i colour
        ## ordering
        index_order_ri_simu = np.argsort(r_simu_selec - i_simu_selec)
        index_order_ri_cosmos = np.argsort(r_selec - i_selec)

        ## assignment
        tmp1 = np.empty(Ns, dtype=float)
        for i_refer, tmp_simu_refer in enumerate(simu_combine_c):

            tmp1[index_order_ri_simu] = cosmos_combine_tmp[i_refer][index_order_ri_cosmos]

            tmp_simu_refer[mask_size] = tmp1

        ############################# 3. random
        index_random = random.sample(range(len(BA_selec)), len(BA_selec))
        for i_refer, tmp_simu_refer in enumerate(simu_combine_r):
            tmp_simu_refer[mask_size] = cosmos_combine_tmp[i_refer][index_random]

    # ++++++++++++ results
    cols_basename = ['Re_physical', 'BA', 'sersic_n']
    methods_name = ['bf', 'c', 'r']
    methods_result = [simu_combine_bf, simu_combine_c, simu_combine_r]

    # save back to some original DataFrame
    if simu_ori_cata is not None:
        if (mask_z_simu is not None) and (mask_mag_simu is not None):
            mask_tmp = mask_z_simu & mask_mag_simu
        elif (mask_z_simu is not None) and (mask_mag_simu is None):
            mask_tmp = mask_z_simu
        elif (mask_z_simu is None) and (mask_mag_simu is not None):
            mask_tmp = mask_mag_simu
        else:
            raise Exception('mask_z_simu and mask_mag_simu cannot both be None!')

        for i_method, method_name in enumerate(methods_name):
            for i_col, col_basename in enumerate(cols_basename):
                simu_ori_cata.loc[mask_tmp, f'{col_basename}_{method_name}'] = methods_result[i_method][i_col]
                # arcsec size
                if col_basename=='Re_physical':
                    simu_ori_cata.loc[mask_tmp, f'Re_arcsec_{method_name}'] = methods_result[i_method][i_col]/simu_cata_selec['Dist_A'].values/1e3/np.pi*180.*3600.   # kpc2arcsec

    # save back to the selected DataFrame
    else:
        for i_method, method_name in enumerate(methods_name):
            for i_col, col_basename in enumerate(cols_basename):
                simu_cata_selec[f'{col_basename}_{method_name}'] = methods_result[i_method][i_col]
                # arcsec size
                if col_basename=='Re_physical':
                    simu_cata_selec[f'Re_arcsec_{method_name}'] = methods_result[i_method][i_col]/simu_cata_selec['Dist_A'].values/1e3/np.pi*180.*3600.   # kpc2arcsec

    logger.debug(f'Sampling and mapping finished {time.time()-time_start}')

    return simu_cata_selec
