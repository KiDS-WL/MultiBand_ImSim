# -*- coding: utf-8 -*-
# @Author: lshuns
# @Date:   2021-07-14 17:04:09
# @Last Modified by:   lshuns
# @Last Modified time: 2021-08-26 13:37:41

### Module to calculate the shear bias with the linear regression
### Two main functions: 
###         mcCalFunc_simple: Calculate the residual shear bias for a given set of simulated catalogue
###         mcCalFunc_reweight_R_SNR_KiDS: KiDS-like data reweighting using R & SNR binning

import random

import numpy as np
import pandas as pd

from scipy import optimize

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
    q1 = _WgQuantile1DFunc(v1, wgs, pq1)

    #Compute quantiles for v2 in each v1 bin
    q2s=[]
    for i in range(len(q1)-1):
        mask = (v1>=q1[i]) & (v1<q1[i+1])
        q2 = _WgQuantile1DFunc(v2[mask], wgs[mask], pq2)
        q2s.append(q2)

    return q1, np.array(q2s)

def mcCalFunc_simple(dataSim, nboot=0, rng_seed_boot=201294):
    """
    Calculate the residual shear bias for a given set of simulated catalogue
        Used columns and their names:
            g1_in, g2_in: input shear
            e1_out, e2_out: measured ellipticity
            shape_weight: shape measurement weights
            tile_label: labels for grouping galaxies
    """

    # tile labels
    tile_labels = np.unique(dataSim['tile_label'])

    # input shear
    g1_in = np.array(dataSim['g1_in'])
    g2_in = np.array(dataSim['g2_in'])
    ## unique pairs
    g_in_unique = np.empty([len(g1_in), 2])
    g_in_unique[:, 0] = g1_in
    g_in_unique[:, 1] = g2_in
    g_in_unique = np.unique(g_in_unique, axis=0)

    # check sufficiency
    if(len(g_in_unique)*len(tile_labels)<3):
        raise Exception(f'Cannot do the linear regression, less than 3 points!')

    # get measured shear for each tile
    g1_out_list = []
    g2_out_list = []
    sigma_out_list = []
    g1_in_list = []
    g2_in_list = []
    N_points = 0
    for tile_label in tile_labels:
        dataSim_selec = dataSim[dataSim['tile_label']==tile_label]

        # input shear
        g1Sim = np.array(dataSim_selec['g1_in'])
        g2Sim = np.array(dataSim_selec['g2_in'])

        # measured ellipticity
        e1Sim = np.array(dataSim_selec['e1_out'])
        e2Sim = np.array(dataSim_selec['e2_out'])

        # weights
        wgSim = np.array(dataSim_selec['shape_weight'])

        # get measured shear by averaging over the measured ellipticity
        for g_in_tmp in g_in_unique:
            maskShear = (g1Sim == g_in_tmp[0]) & (g2Sim == g_in_tmp[1])
            if np.sum(wgSim[maskShear])>0:
                g1_out_list.append(np.average(e1Sim[maskShear], weights=wgSim[maskShear]))
                g2_out_list.append(np.average(e2Sim[maskShear], weights=wgSim[maskShear]))
                sigma_out_list.append(1./(np.sum(wgSim[maskShear]))**0.5)

                g1_in_list.append(g_in_tmp[0])
                g2_in_list.append(g_in_tmp[1])
                N_points += 1
    g1_out_list = np.array(g1_out_list)
    g2_out_list = np.array(g2_out_list)
    sigma_out_list = np.array(sigma_out_list)
    g1_in_list = np.array(g1_in_list)
    g2_in_list = np.array(g2_in_list)
    # print('Total number of points for mc fitting', N_points)

    # fitting input and output shear with a linear function
    m1c1, err1 = optimize.curve_fit(_mcFitFunc, xdata=g1_in_list, ydata=g1_out_list, sigma=sigma_out_list)
    m1 = m1c1[0]
    m1_err = (err1[0, 0])**0.5
    c1 = m1c1[1]
    c1_err = (err1[1, 1])**0.5
    ##
    m2c2, err2 = optimize.curve_fit(_mcFitFunc, xdata=g2_in_list, ydata=g2_out_list, sigma=sigma_out_list)
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
        g_in_out_points = np.stack((g1_in_list, g2_in_list, g1_out_list, g2_out_list, sigma_out_list), axis=-1)
        ## resampling
        sigma_m1_boots = np.zeros(nboot)
        sigma_m2_boots = np.zeros(nboot)
        sigma_c1_boots = np.zeros(nboot)
        sigma_c2_boots = np.zeros(nboot)
        for BS_index in range(nboot):
            # random seed
            random.seed(rng_seed_boot + 120*BS_index)
            # re-sampling with replacement
            g_in_out_BS = np.array(random.choices(g_in_out_points, k=N_points))
            # fitting
            m1c1, err1 = optimize.curve_fit(_mcFitFunc, xdata=g_in_out_BS[:, 0], ydata=g_in_out_BS[:, 2], sigma=g_in_out_BS[:, 4])
            sigma_m1_boots[BS_index] = m1c1[0] - m1
            sigma_c1_boots[BS_index] = m1c1[1] - c1
            m2c2, err2 = optimize.curve_fit(_mcFitFunc, xdata=g_in_out_BS[:, 1], ydata=g_in_out_BS[:, 3], sigma=g_in_out_BS[:, 4])
            sigma_m2_boots[BS_index] = m2c2[0] - m2
            sigma_c2_boots[BS_index] = m2c2[1] - c2
        # choose the 1 sigma range
        name_boots = ['m1', 'm2', 'c1', 'c2']
        for i_name, sigma_boots in enumerate([sigma_m1_boots, sigma_m2_boots, sigma_c1_boots, sigma_c2_boots]):
            name = name_boots[i_name]
            res[f'{name}_err_BS_16'] = np.abs(np.percentile(sigma_boots, 16))
            res[f'{name}_err_BS_84'] = np.abs(np.percentile(sigma_boots, 84))

    return res

def mcCalFunc_least_squares(dataSim, nboot=500):
    """
    Calculate the residual shear bias for a given set of simulated catalogue
        by requiring sum(g_out - (1+m)g_in) -> min

        It is slower than simple method, but gives consistent results.
        Better used for cases with random input shears.

        Used columns and their names:
            g1_in, g2_in: input shear
            e1_out, e2_out: measured ellipticity
            shape_weight: shape measurement weights
    """

    # for saving
    res = {}

    # input shear
    g1_in = np.array(dataSim['g1_in'])
    g2_in = np.array(dataSim['g2_in'])
    Ngal = len(g1_in)

    # out shear
    e1_out = np.array(dataSim['e1_out'])
    e2_out = np.array(dataSim['e2_out'])
    wgSim = np.array(dataSim['shape_weight'])
    totwgSim = np.sum(wgSim)

    # get least square values
    g1e1 = g1_in*e1_out*wgSim
    g1g1 = g1_in*g1_in*wgSim
    res['m1'] = np.sum(g1e1) / np.sum(g1g1) - 1
    g2e2 = g2_in*e2_out*wgSim
    g2g2 = g2_in*g2_in*wgSim
    res['m2'] = np.sum(g2e2) / np.sum(g2g2) - 1

    # get c
    res['c1'] = np.sum((e1_out - (1+res['m1']) * g1_in)*wgSim) / totwgSim
    res['c2'] = np.sum((e2_out - (1+res['m2']) * g2_in)*wgSim) / totwgSim

    # get error from boots
    Relications_index = np.array([np.random.randint(Ngal, size=Ngal) for _ in range(nboot)])

    m1_BS = np.sum(g1e1[Relications_index], axis = 1) / np.sum(g1g1[Relications_index], axis = 1) - 1
    m2_BS = np.sum(g2e2[Relications_index], axis = 1) / np.sum(g2g2[Relications_index], axis = 1) - 1

    totwgSim_BS = np.sum(wgSim[Relications_index], axis = 1)
    c1_BS = np.sum((e1_out[Relications_index] - (1+m1_BS.reshape(nboot, 1)) * g1_in[Relications_index])*wgSim[Relications_index], axis = 1) / totwgSim_BS
    c2_BS = np.sum((e2_out[Relications_index] - (1+m2_BS.reshape(nboot, 1)) * g2_in[Relications_index])*wgSim[Relications_index], axis = 1) / totwgSim_BS

    ## difference 
    sigma_m1_BS = m1_BS - res['m1'] 
    sigma_m2_BS = m2_BS - res['m2'] 
    sigma_c1_BS = c1_BS - res['c1'] 
    sigma_c2_BS = c2_BS - res['c2'] 

    ## choose the 1 sigma range
    name_BS = ['m1', 'm2', 'c1', 'c2']
    for i_name, sigma_BS in enumerate([sigma_m1_BS, sigma_m2_BS, sigma_c1_BS, sigma_c2_BS]):
        name = name_BS[i_name]
        res[f'{name}_err_BS_16'] = np.abs(np.percentile(sigma_BS, 16))
        res[f'{name}_err_BS_84'] = np.abs(np.percentile(sigma_BS, 84))

    return res

def mcCalFunc_reweight_R_SNR_KiDS(dataSim, dataReal, Nbin_SNR=20, Nbin_R=20, nboot=0, rng_seed_boot=201294, fitting_method='simple'):
    """
    Calculating the residual shear bias with KiDS-like data reweighting.
        The two bins performed in resolution parameter R & SNR.
        Used columns and their names:
            Simulation:
                g1_in, g2_in: input shear
                e1_out, e2_out: measured ellipticity
                r_ab_out: mean of the major and minor axis scalelengths
                shape_snr: shape measurement SNR
                shape_weight: shape measurement weights
                tile_label: labels for grouping galaxies
                psf_size: PSF size from measured quadrupole moments (psf_Q11 * psf_Q22 - psf_Q12**2)**0.5
            Data:
                e1_out, e2_out: measured ellipticity
                r_ab_out: mean of the major and minor axis scalelengths
                shape_snr: shape measurement SNR
                shape_weight: shape measurement weights
                psf_size: PSF size from measured quadrupole moments
    """

    # parameters for binning
    ## simulation
    ### circularised galaxy size
    eSim = np.hypot(dataSim['e1_out'], dataSim['e2_out'])
    r_abSim = np.array(dataSim['r_ab_out']) * np.sqrt((1.-eSim)/(1.+eSim))
    ### PSF size
    size_psfSim = np.array(dataSim['psf_size'])
    ### resolution parameter
    RSim = np.array(size_psfSim/(r_abSim**2+size_psfSim))
    ### SNR
    snrSim = np.array(dataSim['shape_snr'])
    ### weights
    wgSim = np.array(dataSim['shape_weight'])

    ## Data
    ### circularised galaxy size
    eReal = np.hypot(dataReal['e1_out'], dataReal['e2_out'])
    r_abReal = np.array(dataReal['r_ab_out']) * np.sqrt((1.-eReal)/(1.+eReal))
    ### PSF size
    size_psfReal = np.array(dataReal['psf_size'])
    ### resolution parameter
    RReal = np.array(size_psfReal/(r_abReal**2+size_psfReal))
    ### SNR
    snrReal = np.array(dataReal['shape_snr'])
    ### weights
    wgReal = np.array(dataReal['shape_weight'])

    # 2D binning such that each bin contains equal sum of weights for simulation
    bin_SNR_bounds, bin_R_bounds = _WgBin2DFunc(snrSim, RSim, wgSim, Nbin_SNR, Nbin_R)

    # calculate mc in each bin
    wgRealSum = 0
    res = {}
    for i in range(Nbin_SNR):
        mask1Sim = (snrSim>=bin_SNR_bounds[i]) & (snrSim<bin_SNR_bounds[i+1])
        mask1Real = (snrReal>=bin_SNR_bounds[i]) & (snrReal<bin_SNR_bounds[i+1])
        # print(f'In SNR bin: [{bin_SNR_bounds[i]:.2f}, {bin_SNR_bounds[i+1]:.2f})')

        for j in range(Nbin_R):
            mask2Sim = (RSim>=bin_R_bounds[i][j]) & (RSim<bin_R_bounds[i][j+1])
            mask2Real = (RReal>=bin_R_bounds[i][j]) & (RReal<bin_R_bounds[i][j+1])
            maskSim = mask1Sim & mask2Sim
            maskReal = mask1Real & mask2Real
            # print(f'In R bin: [{bin_R_bounds[i][j]:.2f}, {bin_R_bounds[i][j+1]:.2f})')
            # print('Number simulation, data', np.sum(maskSim), np.sum(maskReal))

            # weights from data for re-weighting
            wgRealBin = np.sum(wgReal[maskReal])
            wgRealSum += wgRealBin

            # mc fitting
            if fitting_method == 'simple':
                res_bin = mcCalFunc_simple(dataSim[maskSim], nboot=nboot, rng_seed_boot=rng_seed_boot)
            elif fitting_method == 'lsq':
                res_bin = mcCalFunc_least_squares(dataSim[maskSim], nboot=nboot)
            else:
                raise Exception(f'Unsupported fitting_method {fitting_method}!')
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


####### test
if __name__ == '__main__':
    import time
    ### how many bootstrap needed for stable errors
    infile = '/disks/shear16/ssli/ImSim/output/skills_v03s_part1_college_noise/Phase1.feather'
    dataSim = pd.read_feather(infile)

    ## used columns
    dataSim = pd.DataFrame({
        'g1_in': np.array(dataSim['g1_in']),
        'g2_in': np.array(dataSim['g2_in']),
        'e1_out': np.array(dataSim['e1_LF_r']),
        'e2_out': np.array(dataSim['e2_LF_r']),
        'shape_weight': np.array(dataSim['recal_weight_LF_r']),
        'tile_label': np.array(dataSim['tile_label'])})

    nboot_list = [100, 200, 300, 500, 1000]
    for nboot in nboot_list:
        time_start = time.time()
        print('nboot', nboot)
        res = mcCalFunc_simple(dataSim, nboot=nboot, rng_seed_boot=201294)
        print('running time', time.time()-time_start)
        print(res['m1_err'], res['m1_err_BS_16'], res['m1_err_BS_84'], res['m2_err'], res['m2_err_BS_16'], res['m2_err_BS_84'])

    # nboot 100
    # running time 8.808491945266724
    # 0.006147375969045828 0.005771725802050268 0.005623332583636224 0.004290754593977521 0.00402308073550633 0.004544156387100052
    # nboot 200
    # running time 9.705386400222778
    # 0.006147375969045828 0.005792218207473167 0.006321227094274578 0.004290754593977521 0.004399692004517568 0.004621210642324763
    # nboot 300
    # running time 9.26894235610962
    # 0.006147375969045828 0.005831290461518759 0.00620355103613899 0.004290754593977521 0.004399167733894536 0.00497967626855521
    # nboot 500
    # running time 9.705554962158203
    # 0.006147375969045828 0.005792218207473166 0.006070841648675191 0.004290754593977521 0.004288162007978422 0.004609461358955512
    # nboot 1000
    # running time 10.840693712234497
    # 0.006147375969045828 0.005981825033820744 0.0062529950123931186 0.004290754593977521 0.004146491134093678 0.004138712910925149


