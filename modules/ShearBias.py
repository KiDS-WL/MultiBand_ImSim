# @Author: lshuns
# @Date:   2021-02-02, 19:21:29
# @Last modified by:   lshuns
# @Last modified time: 2021-03-03, 18:18:01

### Module to calculate the m bias

import random
import logging

import numpy as np
import pandas as pd

from scipy import optimize

logger = logging.getLogger(__name__)

def _mcFitFunc(x, m, c):
    """
    Shear bias function.
    """
    return (1.+m)*x+c

def _WgQuantile1DFunc(values, weights, pq):
    """
    Calculate the weighted quantile by given probabilities
        designed for 1D numpy array.
    """

    # Sort the data
    ind_sorted = np.argsort(values)
    v_sorted = values[ind_sorted]
    wg_sorted = weights[ind_sorted]

    # Compute the auxiliary arrays
    Sn = np.cumsum(wg_sorted)
    Pn = (Sn-0.5*wg_sorted)/np.sum(wg_sorted)

    # Get the quantiles
    res = np.interp(pq, Pn, v_sorted)

    return res

def _WgBin2DFunc(v1, v2, wgs, Nbin1, Nbin2):
    """
    Calculate the weighted quantile by given bin numbers
        designed for 2D numpy array
    """

    # Define the probabilities for the quantiles based on the number of bins
    pq1 = np.linspace(0,1.0,Nbin1+1)
    pq2 = np.linspace(0,1.0,Nbin2+1)

    # Calculate quantiles for v1
    q1 = WgQuantile1DFunc(v1, wgs, pq1)

    #Compute quantiles for v2 in each v1 bin
    q2s=[]
    for i in range(len(q1)-1):
        mask = (v1>=q1[i]) & (v1<q1[i+1])
        q2 = WgQuantile1DFunc(v2[mask], wgs[mask], pq2)
        q2s.append(q2)

    return q1, np.array(q2s)

def mcCalFunc_simple(dataSim, nboot=0, rng_seed_boot=201294):
    """
    Calculate the residual shear bias for a given set of simulated catalogue
        Used columns and their names:
            g1_in, g2_in: input shear
            e1_out, e2_out: measured ellipticity
            shape_weight: shape measurement weights
    """

    # input shear
    g1_in = np.array(dataSim['g1_in'])
    g2_in = np.array(dataSim['g2_in'])
    ## unique pairs
    g_in_unique = np.empty([len(g1_in), 2])
    g_in_unique[:, 0] = g1_in
    g_in_unique[:, 1] = g2_in
    g_in_unique = np.unique(g_in_unique, axis=0)
    N_g_in_unique = len(g_in_unique)
    ### check sufficiency
    if(N_g_in_unique<3):
        raise Exception(f'Cannot do the mc fitting, less than 3 pairs of input shear!')

    # measured ellipticity
    e1Sim = np.array(dataSim['e1_out'])
    e2Sim = np.array(dataSim['e2_out'])

    # weights
    wgSim = np.array(dataSim['shape_weight'])

    # get measured shear by averaging over the measured ellipticity
    g_outANDsigma = np.empty([N_g_in_unique, 3])
    for i_g, g_in_tmp in enumerate(g_in_unique):
        maskShear = (g1_in == g_in_tmp[0]) & (g2_in == g_in_tmp[1])
        g_outANDsigma[i_g, 0] = np.average(e1Sim[maskShear], weights=wgSim[maskShear])
        g_outANDsigma[i_g, 1] = np.average(e2Sim[maskShear], weights=wgSim[maskShear])
        g_outANDsigma[i_g, 2] = 1./(np.sum(wgSim[maskShear]))**0.5

    # mc fitting with input and output shear
    m1c1, err1 = optimize.curve_fit(mcFitFunc, xdata=g_in_unique[:, 0], ydata=g_outANDsigma[:, 0], sigma=g_outANDsigma[:, 2])
    m1 = m1c1[0]
    m1_err = (err1[0, 0])**0.5
    c1 = m1c1[1]
    c1_err = (err1[1, 1])**0.5
    ##
    m2c2, err2 = optimize.curve_fit(mcFitFunc, xdata=g_in_unique[:, 1], ydata=g_outANDsigma[:, 1], sigma=g_outANDsigma[:, 2])
    m2 = m2c2[0]
    m2_err = (err2[0, 0])**0.5
    c2 = m2c2[1]
    c2_err =(err2[1, 1])**0.5

    # collect results
    res = {'m1': m1, 'm2': m2,
            'c1': c1, 'c2': c2,
            'm1_err': m1_err, 'm2_err': m2_err,
            'c1_err': c1_err, 'c2_err': c2_err}

    # Bootstrap for errors
    if nboot > 0:
        g_in_out_pairs = np.hstack([g_in_unique, g_outANDsigma])
        ## resampling
        sigma_m1_boots = np.zeros(nboot)
        sigma_m2_boots = np.zeros(nboot)
        sigma_c1_boots = np.zeros(nboot)
        sigma_c2_boots = np.zeros(nboot)
        for BS_index in range(nboot):
            # random seed
            random.seed(rng_seed_boot + 120*BS_index)
            # re-sampling with replacement
            g_in_out_BS = np.array(random.choices(g_in_out_pairs, k=N_g_in_unique))
            # fitting
            m1c1, err1 = optimize.curve_fit(mcFitFunc, xdata=g_in_out_BS[:, 0], ydata=g_in_out_BS[:, 2], sigma=g_in_out_BS[:, 4])
            sigma_m1_boots[BS_index] = m1c1[0] - m1
            sigma_c1_boots[BS_index] = m1c1[1] - c1
            m2c2, err2 = optimize.curve_fit(mcFitFunc, xdata=g_in_out_BS[:, 1], ydata=g_in_out_BS[:, 3], sigma=g_in_out_BS[:, 4])
            sigma_m2_boots[BS_index] = m2c2[0] - m2
            sigma_c2_boots[BS_index] = m2c2[1] - c2
        # choose the 1 sigma range
        name_boots = ['m1', 'm2', 'c1', 'c2']
        for i_name, sigma_boots in enumerate([sigma_m1_boots, sigma_m2_boots, sigma_c1_boots, sigma_c2_boots]):
            name = name_boots['i_name']
            res[f'{name}_err_BS_16'] = np.abs(np.percentile(sigma_boots, 16))
            res[f'{name}_err_BS_84'] = np.abs(np.percentile(sigma_boots, 84))

    return res

def mcCalFunc_reweight_R_SNR(dataSim, dataReal, Nbin_SNR=20, Nbin_R=20, nboot=0, rng_seed_boot=201294):
    """
    Calculating the residual shear bias with data reweighting.
        The two bins performed in resolution parameter R & SNR.
    """

    # parameters for binning
    ## simulation
    ### circularised galaxy size
    eSim = np.hypot(dataSim['e1_out'], dataSim['e2_out'])
    size_abSim = np.array(dataSim['size_out']) * np.sqrt((1.-eSim)/(1.+eSim))
    ### PSF size
    try:
        size_psfSim = np.sqrt(dataSim['psf_Q11'] * dataSim['psf_Q22'] - dataSim['psf_Q12']**2)
    except KeyError:
        size_psfSim = np.array(dataSim['psf_size'])
    ### resolution parameter
    RSim = size_psfSim/(size_abSim**2+size_psfSim)
    ### SNR
    snrSim = np.array(dataSim['shape_snr'])
    ### weights
    wgSim = np.array(dataSim['shape_weight'])

    ## Data
    ### circularised galaxy size
    eReal = np.hypot(dataReal['bias_corrected_e1'], dataReal['bias_corrected_e2'])
    size_abReal = np.array(dataReal['bias_corrected_scalelength_pixels']) * np.sqrt((1.-eReal)/(1.+eReal))
    ### PSF size
    size_psfReal = np.sqrt(dataReal['PSF_Q11'] * dataReal['PSF_Q22'] - dataReal['PSF_Q12']**2)
    ### resolution parameter
    RReal = size_psfReal/(size_abReal**2+size_psfReal)
    ### SNR
    snrReal = np.array(dataReal['model_SNratio'])
    ### weights
    wgReal = np.array(dataReal['recal_weight'])

    # 2D binning such that each bin contains equal sum of weights for simulation
    bin_SNR_bounds, bin_R_bounds = WgBin2DFunc(snrSim, RSim, wgSim, Nbin_SNR, Nbin_R)

    # calculate mc in each bin
    wgRealSum = 0
    res = {}
    for i in range(Nbin_SNR):
        mask1Sim = (snrSim>=bin_SNR_bounds[i]) & (snrSim<bin_SNR_bounds[i+1])
        mask1Real = (snrReal>=bin_SNR_bounds[i]) & (snrReal<bin_SNR_bounds[i+1])

        for j in range(Nbin_R):
            i_bin = i*Nbin_SNR + j
            mask2Sim = (RSim>=bin_R_bounds[i][j]) & (RSim<bin_R_bounds[i][j+1])
            mask2Real = (RReal>=bin_R_bounds[i][j]) & (RReal<bin_R_bounds[i][j+1])
            maskSim = mask1Sim & mask2Sim
            maskReal = mask1Real & mask2Real

            # weights from data for re-weighting
            wgRealBin = np.sum(wgReal[maskReal])
            wgRealSum += wgRealBin

            # mc fitting
            res_bin = mcCalFunc_simple(dataSim[maskSim], nboot=nboot, rng_seed_boot=rng_seed_boot)
            for key, value in res_bin.items():
                if 'err' in key:
                    try:
                        res[key] += (wgRealBin * value)**2.
                    except KeyError:
                        res[key] = (wgRealBin * value)**2.
                else:
                    try:
                        res[key] += (wgRealBin * value)
                    except KeyError:
                        res[key] = wgRealBin * value

    # take the average
    for key, value in res.items():
        if 'err' in key:
            res[key] = value**0.5 / wgRealSum
        else:
            res[key] = value / wgRealSum

    return res
