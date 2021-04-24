# @Author: lshuns
# @Date:   2021-04-20, 12:25:33
# @Last modified by:   lshuns
# @Last modified time: 2021-04-22, 11:29:06

import time
import logging

import numpy as np
import pyvinecopulib as pv
import matplotlib.pyplot as plt

## logging
logging.basicConfig(format='%(name)s : %(levelname)s - %(message)s', level=logging.INFO)
logger = logging.getLogger(__name__)

__all__ = ["TwoBinningMappingFunc", "ZBinningMappingFunc", "SimpleMappingFunc_mag"]

def _histedges_equalN(x, nbin):
    """
    estimate hist edges for even binning
    """
    npt = len(x)
    return np.interp(np.linspace(0, npt, nbin + 1),
                     np.arange(npt),
                     np.sort(x))

def TwoBinningMappingFunc(simu_ori_cata, cosmos_cata, N_zbins=30, N_magbins=40):
    """
    Mapping with redshift binning and magnitude binning first

    Parameter used (names of the columns):
        learning sample (data): z, binning_mag, order_col, size, N_GALFIT_H, BA_GALFIT_HI
        target sample (simu): z, binning_mag, order_col
    """

    # ++++++++++++++++++++++ redshift binning
    ### adaptive binning
    _, edges_z_bins, _ = plt.hist(cosmos_cata['z'], _histedges_equalN(cosmos_cata['z'], N_zbins))
    edges_z_bins[0] = np.min(cosmos_cata['z'])
    edges_z_bins[-1] = np.max(cosmos_cata['z']) + 0.01 # for edge value
    logger.debug(f'edges_z_bins {edges_z_bins}')

    # loop over all bins
    time_start = time.time()
    for i_zbin in range(N_zbins):

        zmin = edges_z_bins[i_zbin]
        zmax = edges_z_bins[i_zbin+1]
        logger.debug(f'redshift bin: [{zmin}, {zmax})')
        mask_z_cosmos = (cosmos_cata['z']>=zmin) & (cosmos_cata['z']<zmax)
        mask_z_simu = (simu_ori_cata['z']>=zmin) & (simu_ori_cata['z']<zmax)

        # ++++++++++++++++++++++ magnitude binning
        # cosmos data
        cosmos_cata_selec0 = cosmos_cata[mask_z_cosmos]
        logger.debug(f'Number of data {len(cosmos_cata_selec0)}')

        # simulation
        simu_cata_selec0 = simu_ori_cata[mask_z_simu]
        logger.debug(f'Number of simulation {len(simu_cata_selec0)}')

        # adaptive binning
        _, edges_mag_bins_tmp, _ = plt.hist(cosmos_cata_selec0['binning_mag'], _histedges_equalN(cosmos_cata_selec0['binning_mag'], N_magbins))
        ### boundaries include all simulation galaxies
        edges_mag_bins_tmp[0] = np.min(simu_cata_selec0['binning_mag'])
        edges_mag_bins_tmp[-1] = np.max(simu_cata_selec0['binning_mag']) + 0.01 # for edge value
        logger.debug(f'edges_mag_bins_tmp {edges_mag_bins_tmp}')
        del cosmos_cata_selec0
        del simu_cata_selec0

        # ++++++++++++++++++++ Sampling & mapping
        for i_magbin in range(N_magbins):
            mag_min = edges_mag_bins_tmp[i_magbin]
            mag_max = edges_mag_bins_tmp[i_magbin+1]
            mask_mag_cosmos = (cosmos_cata['binning_mag']>=mag_min) & (cosmos_cata['binning_mag']<mag_max)
            mask_mag_simu = (simu_ori_cata['binning_mag']>=mag_min) & (simu_ori_cata['binning_mag']<mag_max)

            ## check if simulation has sources in this bin
            if np.sum(mask_z_simu & mask_mag_simu) == 0:
                logger.debug(f'bin z[{zmin}, {zmax}), mag[{mag_min}, {mag_max}) is empty')
                logger.debug('move on')
                continue

            logger.debug(f'working on redshift bin: [{zmin}, {zmax}), mag bin [{mag_min}, {mag_max})...')
            cosmos_cata_selec = cosmos_cata[mask_z_cosmos & mask_mag_cosmos]
            simu_cata_selec = simu_ori_cata[mask_z_simu & mask_mag_simu]

            SimpleMappingFunc(simu_cata_selec, cosmos_cata_selec, simu_ori_cata, mask_z_simu, mask_mag_simu)

    return simu_ori_cata

def ZBinningMappingFunc(simu_ori_cata, cosmos_cata, N_zbins=20):
    """
    Mapping with redshift binning first

    Parameter used (names of the columns):
        learning sample (data): z, order_col, size, N_GALFIT_H, BA_GALFIT_HI
        target sample (simu): z, order_col
    """

    # ++++++++++++++++++++++ redshift binning
    ### adaptive binning
    _, edges_z_bins, _ = plt.hist(cosmos_cata['z'], _histedges_equalN(cosmos_cata['z'], N_zbins))
    edges_z_bins[0] = np.min(cosmos_cata['z'])
    edges_z_bins[-1] = np.max(cosmos_cata['z']) + 0.01 # for edge value
    logger.debug(f'edges_z_bins {edges_z_bins}')

    # ++++++++++++++++++++ Sampling & mapping
    time_start = time.time()
    for i_zbin in range(N_zbins):
        zmin = edges_z_bins[i_zbin]
        zmax = edges_z_bins[i_zbin+1]
        mask_z_cosmos = (cosmos_cata['z']>=zmin) & (cosmos_cata['z']<zmax)
        mask_z_simu = (simu_ori_cata['z']>=zmin) & (simu_ori_cata['z']<zmax)

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
    Simple mapping without any binning
        assignment based on the chosen parameter (order_col)

    Parameters used (names of the columns):
        learning sample (data): order_col, size, N_GALFIT_H, BA_GALFIT_HI
        target sample (simu): order_col
    """

    time_start = time.time()

    N_cosmos = len(cosmos_cata_selec)
    logger.debug(f'Number for data {N_cosmos}')
    N_simu = len(simu_cata_selec)
    logger.debug(f'Number for simu {N_simu}')

    # parameters
    ## cosmos ### for sampling
    X = np.empty((N_cosmos, 4))
    X[:, 0] = np.array(cosmos_cata_selec['order_col'])
    X[:, 1] = np.array(cosmos_cata_selec['size'])
    X[:, 2] = np.array(cosmos_cata_selec['N_GALFIT_HI'])
    X[:, 3] = np.array(cosmos_cata_selec['BA_GALFIT_HI'])
    ## simu ### target
    order_col_simu = np.array(simu_cata_selec['order_col'])

    # ++++++++++++++ learning with vine-copulas
    # Transform copula data using the empirical distribution
    u = pv.to_pseudo_obs(X)
    # create pv object directly from data
    cop = pv.Vinecop(data=u)
    # Sample from the copula
    u_sample = cop.simulate(N_simu, seeds=[1, 2, 3, 4])
    # Transform back simulations to the original scale
    X_sample = np.asarray([np.quantile(X[:, i], u_sample[:, i]) for i in range(4)])
    ## transform back to parameters
    order_col_sample = X_sample[0]
    Re_sample = X_sample[1]
    sersic_n_sample = X_sample[2]
    BA_sample = X_sample[3]

    # ++++++++++++ mapping
    ## ordering according to the chosen parameter
    index_order_sample = np.argsort(order_col_sample)
    index_order_simu = np.argsort(order_col_simu)

    ## learned sample
    sample_combine = [Re_sample[index_order_sample], sersic_n_sample[index_order_sample], BA_sample[index_order_sample]]

    ## target sample
    simu_combine = [np.empty(N_simu, dtype=float), np.empty(N_simu, dtype=float), np.empty(N_simu, dtype=float)]
    for i_col, simu_col in enumerate(simu_combine):
        ## assignment
        simu_col[index_order_simu] = sample_combine[i_col]

    # ++++++++++++ results
    cols_basename = ['size', 'sersic_n', 'BA']

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

        for i_col, col_basename in enumerate(cols_basename):
            simu_ori_cata.loc[mask_tmp, f'{col_basename}'] = simu_combine[i_col]

    # save back to the selected DataFrame
    else:
        for i_col, col_basename in enumerate(cols_basename):
            simu_cata_selec[f'{col_basename}'] = simu_combine[i_col]

    logger.debug(f'Sampling and mapping finished {time.time()-time_start}')

    return simu_cata_selec
