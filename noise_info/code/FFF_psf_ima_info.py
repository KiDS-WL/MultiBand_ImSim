# -*- coding: utf-8 -*-
# @Author: lshuns
# @Date:   2022-06-26 15:52:59
# @Last Modified by:   lshuns
# @Last Modified time: 2022-06-27 09:15:13

### extract information from PSF image

import os
import sys
import time
import glob

import numpy as np
import pandas as pd

from scipy import signal
from ctypes import c_float
from astropy.io import fits

import multiprocessing as mp
import matplotlib.pyplot as plt

### +++++++++++++++++++++++++++ I/O

indir = f'/disks/shear16/ssli/ImSim/input/KiDS_PSF_ima/r_tukey'
outpath = os.path.join(indir, 'PSFinfo_summarised.csv')

### +++++++++++++++++++++++++++ functions
def weighted_moments(image_2Darray):
    """
    Gaussian weighted moments as in LF
    """

    # weight function info
    ### sigma of Gaussian weight function
    sigmasq = (2.5)**2.
    ### number of moments iterations for centroid
    niter = 3
    ### distance limit for moments summations
    distsq_max = (image_2Darray.shape[1]/2)**2.

    # get positions
    xgrid, ygrid = np.meshgrid(np.arange(image_2Darray.shape[1]),
                                   np.arange(image_2Darray.shape[0]))

    # first guess of centre
    ## from no weight estimation
    M0 = np.sum(image_2Darray)
    Mx = np.sum(xgrid * image_2Darray) / M0
    My = np.sum(ygrid * image_2Darray) / M0

    # iterate with weight function
    for i_run in range(niter):
        # print('>>> i_run', i_run)

        # distance array
        x2cen = xgrid - Mx
        y2cen = ygrid - My
        x2censq = np.square(x2cen)
        y2censq = np.square(y2cen)
        distsq = x2censq + y2censq

        # weight array
        w = np.exp(-distsq/2./sigmasq)
        ## those out of the distance limit are not used
        w[distsq >= distsq_max] = 0

        # apply weight to the image
        image_2Darray_w = w * image_2Darray

        # get zero moment (sum)
        M0 = np.sum(image_2Darray_w)
        # print(">>> M0",  M0)
        # get first moment (centre)
        Mx = np.sum(xgrid * image_2Darray_w) / M0
        My = np.sum(ygrid * image_2Darray_w) / M0
        # print(">>> Mx, My",  Mx, My)
        # get second moment
        Mxx = np.sum(x2censq * image_2Darray_w) / M0
        Myy = np.sum(y2censq * image_2Darray_w) / M0
        Mxy = np.sum(x2cen * y2cen * image_2Darray_w) / M0
        # print(">>> Mxx, Myy, Mxy",  Mxx, Myy, Mxy)

    return dict(M0=M0, Mx=Mx, My=My, Mxx=Mxx, Myy=Myy, Mxy=Mxy)

def MomentShapes(image_2Darray, weighted=False):

    # get moments
    if weighted:
        moments = weighted_moments(image_2Darray)
    else:
        moments = unweighted_moments(image_2Darray)

    # get shape
    size = (moments['Mxx']*moments['Myy'] - moments['Mxy']**2.)**0.5
    # print('>>> size', size)

    # get e
    eDe = moments['Mxx'] + moments['Myy'] + 2*size
    e1 = (moments['Mxx'] - moments['Myy']) / eDe
    e2 = 2 * moments['Mxy'] / eDe

    return size, e1, e2

# >>>>>>>>>>>>>>>>>>>>>>>>>>> workhorse

# find all tiles
folder_list = glob.glob(os.path.join(indir, 'tile*'))
print('number of folders', len(folder_list))

# save the final results
dtypes = np.dtype(
    [
        ("file", str),
        ("FWHMarcsec", float),
        ("moment_size_pixel", float),
        ("moment_e1", float),
        ("moment_e2", float)
    ]
)
cata_final = pd.DataFrame(np.empty(len(folder_list)*5*32, dtype=dtypes))
print('expected number of files', len(cata_final))

# get info
i_index = 0
for i_indir, infolder in enumerate(folder_list):

    # find all chips
    chip_list = glob.glob(os.path.join(infolder, 'psfIma*'))
    for inpath_ima in chip_list:

        with fits.open(inpath_ima) as hdul:
            FWHM = hdul[0].header['FWHMarcsec']
            psf_ima = hdul[0].data.astype(float)

        ## get shapes
        size_ima_w, e1_ima_w, e2_ima_w = MomentShapes(psf_ima, weighted=True)

        # save
        cata_final.loc[i_index, 'file'] = f'{os.path.basename(infolder)}/{os.path.basename(inpath_ima)}'
        cata_final.loc[i_index, 'FWHMarcsec'] = FWHM
        cata_final.loc[i_index, 'moment_size_pixel'] = size_ima_w
        cata_final.loc[i_index, 'moment_e1'] = e1_ima_w
        cata_final.loc[i_index, 'moment_e2'] = e2_ima_w

        # print(cata_final.loc[i_index])

        i_index += 1

# only save those with information
cata_final = cata_final[~cata_final['file'].isnull()]
        
# save the final results
cata_final.to_csv(outpath, index=False)
print('saved to', outpath)