# -*- coding: utf-8 -*-
# @Author: lshuns
# @Date:   2021-07-23 15:55:55
# @Last Modified by:   lshuns
# @Last Modified time: 2021-07-26 15:32:38

### transfer lensfit psf coefficients to Moffat profile
###### 1. build psf image from lensfit psf coefficients
###### 2. fit psf image with elliptical moffat profile

import os
import re
import time
import glob
import galsim

import numpy as np
import pandas as pd
import multiprocessing as mp
import matplotlib.pyplot as plt

from ctypes import c_float
from astropy.io import fits
from scipy.optimize import least_squares

def Coeff2Ima(x_inchip, y_inchip, chip_id, inpath_psfcoeffs, psf_ima_size=24):
    """
    transfer lensfit psf coefficients to psf image
    """

    # >>>>>>>>>>>>>>>>>> 1. extract useful info from the fits header
    with fits.open(inpath_psfcoeffs) as hdul:

        ## 1. general info from header
        head = hdul[0].header
        chorder = head['chorder']
        order = head['order']
        chipvar = head['chipvar']
        xsampl = head['xsampl']
        ysampl = head['ysampl']
        hxsize = head['hxsize']
        hysize = head['hysize']
        del head

        ## 2. coefficients array
        psfcoeffs = hdul[0].data
        y_npix = np.shape(psfcoeffs)[1]
        x_npix = np.shape(psfcoeffs)[2]

        ## 3. chip position array
        xchip = hdul[1].data['xchip']
        ychip = hdul[1].data['ychip']
        nchip = len(xchip)
    
    # >>>>>>>>>>>>>>>>>> 2. coefficients to image
    noiseless_psf_ima = np.zeros([y_npix, x_npix], dtype=c_float)
    ## chip position to tile position
    x_intile = (x_inchip + xchip[chip_id]*xsampl)/hxsize - 1
    y_intile = (y_inchip + ychip[chip_id]*ysampl)/hysize - 1
    ## add desired elements to build the image
    kd = 0 # index to the desired sets of coefficients
    nc = 0 # running index
    for i in range(order+1):
        for j in range(order-i+1):
            ## chip variance allowed by the lowest order coefficients
            if (chipvar == 1 and (i + j) <= chorder):
                kd = nc + chip_id
                nc += nchip
            else:
                kd = nc
                nc += 1
            noiseless_psf_ima += psfcoeffs[kd,:,:] * np.power(x_intile, i) * np.power(y_intile, j)
    ## shift the origin to the center
    noiseless_psf_ima = np.roll(np.roll(noiseless_psf_ima, y_npix//2, axis=0), x_npix//2 , axis=1)

    # >>>>>>>>>>>>>>>>> 3. crop the center region
    startx = x_npix//2 - (psf_ima_size//2)
    starty = y_npix//2 - (psf_ima_size//2)
    noiseless_psf_ima = noiseless_psf_ima[starty:starty+psf_ima_size, startx:startx+psf_ima_size]
        
    return noiseless_psf_ima

def EllipticalMoffat(fwhm_pixel, beta, q, pa, x_c, y_c, val_c, x_array, y_array):
    """
    elliptical moffat profile
    """

    # major axis
    fwhm_1 = (fwhm_pixel**2 * q)**0.5
    # minor axis
    fwhm_2 = (fwhm_pixel**2 / q)**0.5

    alpha_1 = 0.5*fwhm_1/np.sqrt(2.**(1./beta)-1.)
    alpha_2 = 0.5*fwhm_2/np.sqrt(2.**(1./beta)-1.)

    pa_rad = (180.-pa)/180.*np.pi # degree to radian ## 180.-pa to be consistent with galsim
    A = (np.cos(pa_rad)/alpha_1)**2. + (np.sin(pa_rad)/alpha_2)**2.
    B = (np.sin(pa_rad)/alpha_1)**2. + (np.cos(pa_rad)/alpha_2)**2.
    C = 2.0*np.sin(pa_rad)*np.cos(pa_rad)*(1./alpha_1**2. - 1./alpha_2**2.)

    # modelled psf image
    moffat_psf = val_c/((1.+ A*((x_array-x_c)**2) + B*((y_array-y_c)**2) + C*(x_array-x_c)*(y_array-y_c))**beta)
    
    return moffat_psf

def FitEllipticalMoffat(noiseless_psf_ima, pixel_scale=0.214):
    """
    Fitting a noise-free PSF image with a 2D elliptical moffat profile.
        assume noiseless_psf_ima is normalized to unity
    """

    # >>>>>>>>>>>>>>>> 1. center is fixed
    x_c = np.shape(noiseless_psf_ima)[0]//2
    y_c = np.shape(noiseless_psf_ima)[1]//2
    val_c = noiseless_psf_ima[x_c, y_c]

    # >>>>>>>>>>>>>>>> 2. initial guess for fitting values
    ## full width half max
    fwhm_pixel = np.sqrt(np.sum((noiseless_psf_ima>val_c/2.).flatten()))
    ## moffat beta
    beta = 3.
    ## axis ratio
    q = 0.98
    ## position angle
    pa = 90.
    ## combine
    val0 = (fwhm_pixel, beta, q, pa)

    # >>>>>>>>>>>>>>>> 3. bounds
    fwhm_bounds = (0.5*fwhm_pixel, 1.5*fwhm_pixel)
    beta_bounds = (1.5, 4)
    q_bounds = (0.5, 1.0)
    pa_bounds = (0, 180) # degrees
    bounds0 = np.transpose((fwhm_bounds, beta_bounds, q_bounds, pa_bounds))

    # >>>>>>>>>>>>>>>>> 4. residual function
    def err(val_fitting, x_c, y_c, val_c, data):
        # fitted values
        fwhm_pixel, beta, q, pa = val_fitting

        # model
        x_array, y_array = np.indices(data.shape)
        val_model = EllipticalMoffat(fwhm_pixel, beta, q, pa, x_c, y_c, val_c, x_array, y_array)

        # residual
        res = np.ravel(val_model - data)

        return res

    # >>>>>>>>>>>>>>>>> 5. least square fitting
    result = least_squares(err, val0, bounds=bounds0, args=(x_c, y_c, val_c, noiseless_psf_ima))
    ## solution to the fitted values
    fwhm_pixel, beta, q, pa = result.x

    return fwhm_pixel*pixel_scale, beta, q, pa

def GalSimCalib(fwhm_arcsec, beta, q, pa, pixel_scale=0.214, size=24):
    """
    first order calibration to the fitted parameter with galsim image
    """

    # psf image from galsim
    psf = galsim.Moffat(beta=beta, fwhm=fwhm_arcsec, trunc=4.5*fwhm_arcsec)
    psf = psf.shear(q=q, beta=pa*galsim.degrees)
    psf_image = galsim.Image(size, size)
    PSF_lf = psf.shift(0.5*pixel_scale, 0.5*pixel_scale)
    psf_image = PSF_lf.drawImage(image=psf_image, scale=pixel_scale).array

    # fitting on the galsim image
    fwhm_out, beta_out, q_out, pa_out = FitEllipticalMoffat(psf_image, pixel_scale=pixel_scale)

    # calibrated results
    fwhm_new = fwhm_arcsec - (fwhm_out - fwhm_arcsec)
    beta_new = beta - (beta_out - beta)
    q_new = q - (q_out - q)
    pa_new = pa - (pa_out - pa)

    return fwhm_new, beta_new, q_new, pa_new

def paraRun(run_id, outpath, psf_path_list):
    """
    function for parallel run
    """

    results_df = pd.DataFrame(columns=['label', 'file', 'starnum', 'expo_id', 'chip_id', 'InputSeeing_r', 'InputBeta_r', 'seeing_e1_r', 'seeing_e2_r'])
    i_row_res = 0

    # tmp output
    outpath_tmp = outpath.replace('.csv', f'tmp{run_id}.feather')

    for psf_path in psf_path_list:

        tile_label = os.path.basename(psf_path).replace('m', '-').replace('p', '.')[5:]
        
        psffile_list = np.sort(glob.glob(os.path.join(psf_path, 'r_SDSS/lensfit_V1.0.0A_svn_309c/psfs', 'OMEGA*')))
        print('tile', tile_label, 'exposures', len(psffile_list))
        for expo_id, inpath_psfcoeffs in enumerate(psffile_list):

            for chip_id in range(N_chips):

                # >>>>>>>>>>>>>>>>>> 1. check star sufficiency
                with fits.open(inpath_psfcoeffs) as hdul:
                    starnum = hdul[2].data[chip_id][1]
                # drop the chip if number of stars less than required
                if starnum < Nmin_stars:
                    print(f'WARNING: number of stars {starnum} less than required {Nmin_stars}, ignore the chip {chip_id} for {tile_label}!')
                    continue

                # >>>>>>>>>>>>>>>>>> 2. some general values
                results_df.loc[i_row_res, 'label'] = tile_label
                results_df.loc[i_row_res, 'file'] = os.path.basename(inpath_psfcoeffs)
                results_df.loc[i_row_res, 'starnum'] = starnum
                results_df.loc[i_row_res, 'expo_id'] = expo_id
                results_df.loc[i_row_res, 'chip_id'] = chip_id
            
                # >>>>>>>>>>>>>>>>>> 3. select positions and fit
                np.random.seed(np.array(re.findall(r"\d+", tile_label), dtype=np.int).sum()*10000 + chip_id*10 + 1)
                x_inchip_list = np.random.randint(int(0.1*x_inchip_max), int(0.9*x_inchip_max), N_pos)
                np.random.seed(np.array(re.findall(r"\d+", tile_label), dtype=np.int).sum()*10000 + chip_id*10 + 2)
                y_inchip_list = np.random.randint(int(0.1*y_inchip_max), int(0.9*y_inchip_max), N_pos)
                fwhm_final = 0
                beta_final = 0
                q_final = 0
                pa_final = 0
                for i_pos in range(N_pos):

                    # position
                    x_inchip = x_inchip_list[i_pos]
                    y_inchip = y_inchip_list[i_pos]

                    # build image                    
                    noiseless_psf_ima = Coeff2Ima(x_inchip, y_inchip, chip_id, inpath_psfcoeffs, psf_ima_size=psf_ima_size)

                    # fitting
                    fwhm, beta, q, pa = FitEllipticalMoffat(noiseless_psf_ima, pixel_scale=0.214)

                    fwhm_final += fwhm
                    beta_final += beta 
                    q_final += q
                    pa_final += pa

                # >>>>>>>>>>>>>>>>>> 4. average and calibrate for final values of the chip
                fwhm_final = fwhm_final/N_pos
                beta_final = beta_final/N_pos
                q_final = q_final/N_pos
                pa_final = pa_final/N_pos
                # calibrate
                fwhm_final, beta_final, q_final, pa_final = GalSimCalib(fwhm_final, beta_final, q_final, pa_final, pixel_scale=0.214, size=psf_ima_size)

                # >>>>>>>>>>>>>>>>>> 5. save to the form as used by ImSim
                results_df.loc[i_row_res, 'InputSeeing_r'] = fwhm_final
                results_df.loc[i_row_res, 'InputBeta_r'] = beta_final
                results_df.loc[i_row_res, 'seeing_e1_r'] = (1-q_final) * np.cos(2. * (pa_final/180.*np.pi))
                results_df.loc[i_row_res, 'seeing_e2_r'] = (1-q_final) * np.sin(2. * (pa_final/180.*np.pi))

                i_row_res += 1

    results_df.to_feather(outpath_tmp)

if __name__ == '__main__':

    ############# main run
    ### build a catalogue with moffat parameters for each chip each exposure

    N_proc = 64
    psf_ima_size = 12
    N_pos = 50
    Nmin_stars = 1
    N_chips = 32
    x_inchip_max, y_inchip_max = 2110, 4160

    outpath = '../kids_dr4_psf_moffat_fromcoeffs.csv'

    # >>>>>>>>> actual run
    start_time = time.time()
    psf_path_list = glob.glob('/disks/shear15/KiDS/KIDSCOLLAB_V1.0.0/KIDS_*')
    ## split according to running cores
    n = int(np.ceil(len(psf_path_list)/N_proc))
    psf_path_list_chunks = [psf_path_list[i:i+n] for i in range(0, len(psf_path_list), n)]
    ## parallel run
    p_list = []
    for run_id, psf_path_list in enumerate(psf_path_list_chunks):
        p = mp.Process(target=paraRun, args=(run_id, outpath, psf_path_list))
        p.start()
        p_list.append(p)
    for p in p_list:
        p.join()

    ## combine all outputs
    cata_final = []
    for run_id in range(len(psf_path_list_chunks)):
        infile = outpath.replace('.csv', f'tmp{run_id}.feather')
        cata = pd.read_feather(infile)
        cata_final.append(cata)

        ## remove tmp files
        os.remove(infile)

    cata_final = pd.concat(cata_final)
    cata_final.reset_index(drop=True, inplace=True)
    cata_final.to_csv(outpath, index=False)
    print('final results saved to', outpath)
    print('all finished in', time.time()-start_time, 's')
    # all finished in 4416.606633424759 s

    # ############# check FitEllipticalMoffat
    # ### use simulated psf image
    # import plotting

    # psf_ima_size = 12

    # def PSFima(seeing, moffat_beta, psf_q, psf_pa, pixel_scale=0.214, size=48):
    #     """
    #     draw a single PSF image
    #     """

    #     psf = galsim.Moffat(beta=moffat_beta, fwhm=seeing)
    #     if psf_q < 1:
    #         # psf = psf.shear(q=psf_q, beta=psf_pa*galsim.degrees)
    #         psf_e1 = (1-psf_q) * np.cos(2. * (psf_pa/180.*np.pi))
    #         psf_e2 = (1-psf_q) * np.sin(2. * (psf_pa/180.*np.pi))
    #         psf_e = np.sqrt(psf_e1**2+psf_e2**2)
    #         psf_g1, psf_g2 = psf_e1/(2-psf_e), psf_e2/(2-psf_e)
    #         psf = psf.shear(g1=psf_g1, g2=psf_g2)

    #     psf_image = galsim.Image(size, size)
    #     PSF_lf = psf.shift(0.5*pixel_scale, 0.5*pixel_scale)
    #     psf_image = PSF_lf.drawImage(image=psf_image, scale=pixel_scale).array

    #     ## normalized to unity
    #     psf_image = psf_image/np.sum(psf_image)

    #     return psf_image

    # ## 2. 100 random psf for systematics
    # N_test = 100
    # seeing_in_list = np.random.uniform(0.5, 1.0, N_test)
    # beta_in_list = np.random.uniform(2.0, 4.0, N_test)
    # q_in_list = np.random.uniform(0.85, 1.0, N_test)
    # pa_in_list = np.random.uniform(0., 180, N_test)
    # ## for output ori
    # seeing_ori_list = np.zeros(N_test)
    # beta_ori_list = np.zeros(N_test)
    # q_ori_list = np.zeros(N_test)
    # pa_ori_list = np.zeros(N_test)
    # ## for output calibrated
    # seeing_new_list = np.zeros(N_test)
    # beta_new_list = np.zeros(N_test)
    # q_new_list = np.zeros(N_test)
    # pa_new_list = np.zeros(N_test)

    # for i_psf in range(N_test):
    #     print('>>>>>>>>>>>>>>>>>>>')

    #     noiseless_psf_ima = PSFima(seeing_in_list[i_psf], beta_in_list[i_psf], q_in_list[i_psf], pa_in_list[i_psf], size=psf_ima_size)
    #     print(seeing_in_list[i_psf], beta_in_list[i_psf], q_in_list[i_psf], pa_in_list[i_psf])

    #     ## center
    #     x_c = np.shape(noiseless_psf_ima)[0]//2
    #     y_c = np.shape(noiseless_psf_ima)[1]//2
    #     ## center value
    #     val_c = noiseless_psf_ima[x_c, y_c]

    #     # >>>>>>>>>>> original fitting
    #     fwhm_arcsec, beta, q, pa = FitEllipticalMoffat(noiseless_psf_ima)
    #     print(fwhm_arcsec, beta, q, pa)
    #     seeing_ori_list[i_psf] = fwhm_arcsec
    #     beta_ori_list[i_psf] = beta
    #     q_ori_list[i_psf] = q
    #     pa_ori_list[i_psf] = pa

    #     # >>>>>>>>>>> calibrate
    #     fwhm_new, beta_new, q_new, pa_new = GalSimCalib(fwhm_arcsec, beta, q, pa, pixel_scale=0.214, size=psf_ima_size)
    #     print(fwhm_new, beta_new, q_new, pa_new)
    #     seeing_new_list[i_psf] = fwhm_new
    #     beta_new_list[i_psf] = beta_new
    #     q_new_list[i_psf] = q_new
    #     pa_new_list[i_psf] = pa_new

    # # build the difference
    # dseeing_ori = (seeing_in_list - seeing_ori_list)/seeing_in_list
    # dbeta_ori = (beta_in_list - beta_ori_list)/beta_in_list
    # dq_ori = (q_in_list - q_ori_list)/q_in_list
    # dpa_ori = (pa_in_list - pa_ori_list)/pa_in_list
    # ##
    # dseeing_new = (seeing_in_list - seeing_new_list)/seeing_in_list
    # dbeta_new = (beta_in_list - beta_new_list)/beta_in_list
    # dq_new = (q_in_list - q_new_list)/q_in_list
    # dpa_new = (pa_in_list - pa_new_list)/pa_in_list

    # # ### plots
    # ## Line
    # outpath_list = ['../plots/dseeing_seeing.png', '../plots/dbeta_beta.png', '../plots/dq_q.png', '../plots/dpa_pa.png', '../plots/dbeta_seeing.png']
    # xvals_list = [[seeing_in_list, seeing_in_list], [beta_in_list, beta_in_list], [q_in_list, q_in_list], [pa_in_list, pa_in_list], [seeing_in_list, seeing_in_list]]
    # yvals_list = [[dseeing_ori, dseeing_new], [dbeta_ori, dbeta_new], [dq_ori, dq_new], [dpa_ori, dpa_new], [dbeta_ori, dbeta_new]]
    # para_name_list = ['seeing', 'beta', 'q', 'position angle', 'dbeta vs. seeing']
    # yrange_list = [[-0.075, 0.01], [-0.1, 0.01], [-0.015, 0.003], [-0.01, 0.003], [-0.1, 0.01]]

    # for i_plot in range(len(outpath_list)):
    #     outpath = outpath_list[i_plot]
    #     xvals = xvals_list[i_plot]
    #     yvals = yvals_list[i_plot]

    #     COLORs = ['k', 'r']
    #     LABELs = ['original', 'calibrated']
    #     LINEs = ['', '']
    #     nbins = 20
    #     YRANGE = yrange_list[i_plot]
    #     XLABEL = 'input'
    #     YLABEL = '(input-output)/input'
    #     TITLE = para_name_list[i_plot]

    #     hlines = [0.0]

    #     plotting.LinePlotFunc(outpath,
    #             xvals, yvals,
    #             COLORs, LABELs=LABELs, LINEs=LINEs, LINEWs=None, POINTs=None, POINTSs=None, fillstyles=None,
    #             XRANGE=None, YRANGE=YRANGE,
    #             XLABEL=XLABEL, YLABEL=YLABEL, TITLE=TITLE,
    #             xtick_min_label=True, xtick_spe=None, ytick_min_label=True, ytick_spe=None,
    #             vlines=None, vline_styles=None, vline_colors=None, vline_labels=None, vline_widths=None,
    #             hlines=hlines, hline_styles=None, hline_colors=None, hline_labels=None, hline_widths=None,
    #             xlog=False, invertX=False, ylog=False, invertY=False, loc_legend='best')

