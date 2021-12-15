# -*- coding: utf-8 -*-
# @Author: lshuns
# @Date:   2021-11-23 14:21:21
# @Last Modified by:   lshuns
# @Last Modified time: 2021-11-27 11:13:09

### build psf image from lensfit psf coefficients

import os
import sys
import time
import glob

import numpy as np

from scipy import signal
from ctypes import c_float
from astropy.io import fits

############################ setup

Nmin_stars = 1
N_chips = 32
pixel_scale = 0.214
lensfit_V = 'lensfit_V1.0.0A_svn_309c'
band = 'r_SDSS'
save_band = 'r'
x_inchip, y_inchip = int(2110/2.), int(4160/2.)

# window function related
n_back=2
beta_kaiser=14

outdir_ori = f'/disks/shear10/ssli/ImSim/input/KiDS_PSF_ima/{save_band}_ori'
outdir_smooth = f'/disks/shear10/ssli/ImSim/input/KiDS_PSF_ima/{save_band}'

############################ the functions
def Coeff2Ima(x_inchip, y_inchip, chip_id, inpath_psfcoeffs):
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
        
    return noiseless_psf_ima

def WindowIma(ima, n_back=2, beta_kaiser=14):
    '''
    Apply window to the image to suppress noise at the edge
    '''

    # get the size of image
    ima_size = np.array(ima.shape)
    ima_cen = (ima_size/2).astype(int)

    # define boundary from negative pixel values
    i_neg, j_neg = np.where(ima < 0)
    ## distance to center
    ij_neg = np.zeros((len(i_neg), 2))
    ij_neg[:, 0] = i_neg
    ij_neg[:, 1] = j_neg
    dis_neg = ij_neg-ima_cen
    dis_neg = np.hypot(dis_neg[:, 0], dis_neg[:, 1])
    ## min dist as the boundary
    boundary = int(np.min(dis_neg) - n_back)
    ## preserve centre within boundary
    ima_good = ima[(ima_cen[0]-boundary):(ima_cen[0]+boundary), (ima_cen[1]-boundary):(ima_cen[1]+boundary)]
    ima_good_size = ima_good.shape
    print('good psf size', ima_good_size)

    # check if boundary is too close to the center
    ## by comparing the center value to the edge values
    Redge = np.max(np.concatenate((ima_good[0], ima_good[-1], ima_good[:, 0], ima_good[:, -1])))/ima[ima_cen[0], ima_cen[1]]
    print('edge value / centre value:', Redge)

    # get 2D kaiser window
    ## 0 axis
    len_bad = ima_size[0] - boundary*2
    kaiser_tmp = signal.windows.kaiser(len_bad, beta=beta_kaiser)
    window1D_i = np.concatenate((kaiser_tmp[:int(len_bad/2)], 
                    np.ones(ima_good_size[0]), 
                    kaiser_tmp[int(len_bad/2):]))
    ## 1 axis
    len_bad = ima_size[1] - boundary*2
    kaiser_tmp = signal.windows.kaiser(len_bad, beta=beta_kaiser)
    window1D_j = np.concatenate((kaiser_tmp[:int(len_bad/2)], 
                    np.ones(ima_good_size[1]), 
                    kaiser_tmp[int(len_bad/2):]))
    window2D = np.sqrt(np.outer(window1D_i, window1D_j))

    # apply window
    return ima*window2D, ima_good_size, Redge

############################ workhorse

# ## >>>>>>>>>>> from coeff
# start_time = time.time()
# psf_path_list = glob.glob('/disks/shear15/KiDS/KIDSCOLLAB_V1.0.0/KIDS_*')
# for psf_path in psf_path_list:

#     tile_label = os.path.basename(psf_path).replace('m', '-').replace('p', '.')[5:]

#     # to save original psf image
#     outdir_ori_tmp = os.path.join(outdir_ori, f'tile{tile_label}')
#     if not os.path.isdir(outdir_ori_tmp):
#         os.mkdir(outdir_ori_tmp)

#     # to save smoothed psf image
#     outdir_smooth_tmp = os.path.join(outdir_smooth, f'tile{tile_label}')
#     if not os.path.isdir(outdir_smooth_tmp):
#         os.mkdir(outdir_smooth_tmp)

#     psffile_list = np.sort(glob.glob(os.path.join(psf_path, f'{band}/{lensfit_V}/psfs', 'OMEGA*')))
#     print('tile', tile_label, 'exposures', len(psffile_list))
#     for expo_id, inpath_psfcoeffs in enumerate(psffile_list):

#         for chip_id in range(N_chips):

#             # >>>>>>>>>>>>>>>>>> 1. check star sufficiency
#             with fits.open(inpath_psfcoeffs) as hdul:
#                 starnum = hdul[2].data[chip_id][1]
#             # drop the chip if number of stars less than required
#             if starnum < Nmin_stars:
#                 print(f'WARNING: number of stars {starnum} less than required {Nmin_stars}, ignore the chip {chip_id} for {tile_label}!')
#                 continue

#             # >>>>>>>>>>>>>>>>>> 2. some values saved to header
#             hdr = fits.Header()
#             hdr['label'] = tile_label
#             hdr['file'] = os.path.basename(inpath_psfcoeffs)
#             hdr['lfV'] = lensfit_V
#             hdr['band'] = band
#             hdr['starnum'] = starnum
#             hdr['expoid'] = expo_id
#             hdr['IMAGEID'] = chip_id +1
#             hdr['pixscale'] = pixel_scale
#             hdr['xinchip'] = x_inchip
#             hdr['yinchip'] = y_inchip
#             ## for smoothed
#             hdr_s = hdr.copy()
#             hdr_s['nback'] = n_back
#             hdr_s['bkaiser'] = beta_kaiser
            
#             # >>>>>>>>>>>>>>>>>> 3. get the image for given position
#             noiseless_psf_ima = Coeff2Ima(x_inchip, y_inchip, chip_id, inpath_psfcoeffs)
#             ## smooth
#             noiseless_psf_ima_win = WindowIma(noiseless_psf_ima, n_back=n_back, beta_kaiser=beta_kaiser)

#             # >>>>>>>>>>>>>>>>>> 4. build fits and save
#             hdu = fits.PrimaryHDU(data=noiseless_psf_ima, header=hdr)
#             outpath = os.path.join(outdir_ori_tmp, f'psfIma_exp{expo_id}_chip{chip_id+1}.fits')
#             hdu.writeto(outpath)
#             ## smooth
#             hdu = fits.PrimaryHDU(data=noiseless_psf_ima_win, header=hdr_s)
#             outpath = os.path.join(outdir_smooth_tmp, f'psfIma_exp{expo_id}_chip{chip_id+1}.fits')
#             hdu.writeto(outpath)

# print('all finished in', time.time()-start_time, 's')
# # all finished in 2049.431708097458 s


## >>>>>>>>>>> from ima
start_time = time.time()
psf_tile_dirs = glob.glob(os.path.join(outdir_ori, 'tile*'))
N_badsmooth = 0
for psf_tile_dir in psf_tile_dirs:

    # to save smoothed psf image
    outdir_smooth_tmp = os.path.join(outdir_smooth, os.path.basename(psf_tile_dir))
    if not os.path.isdir(outdir_smooth_tmp):
        os.mkdir(outdir_smooth_tmp)

    psffile_list = glob.glob(os.path.join(psf_tile_dir, 'psfIma_*.fits'))
    print('tile', os.path.basename(psf_tile_dir), 'N_files', len(psffile_list))
    for psffile in psffile_list:

        # get image and hdr
        with fits.open(psffile) as hdul:
            hdr = hdul[0].header
            ima = hdul[0].data

        ## smooth
        noiseless_psf_ima_win, ima_good_size, Redge = WindowIma(ima, n_back=n_back, beta_kaiser=beta_kaiser)

        ## discard images with bad smoothing
        if (ima_good_size[0] < 10) and (Redge > 0.1) :
            N_badsmooth += 1
            print('bad smoothing for', os.path.basename(psffile))
            continue

        ## some head info
        hdr_s = hdr.copy()
        hdr_s['nback'] = n_back
        hdr_s['bkaiser'] = beta_kaiser
        hdr_s['goodNx'] = ima_good_size[0]
        hdr_s['goodNy'] = ima_good_size[1]
        hdr_s['Redge'] = Redge

        # >>>>>>>>>>>>>>>>>> 4. build fits and save
        ## smooth
        hdu = fits.PrimaryHDU(data=noiseless_psf_ima_win, header=hdr_s)
        outpath = os.path.join(outdir_smooth_tmp, os.path.basename(psffile))
        hdu.writeto(outpath)

print('Number of discard psfs', N_badsmooth)
print('all finished in', time.time()-start_time, 's')
# Number of discard psfs 84
# all finished in 1526.010176897049 s
