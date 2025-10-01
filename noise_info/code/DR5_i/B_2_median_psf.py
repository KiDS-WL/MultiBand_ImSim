# -*- coding: utf-8 -*-
# @Author: lshuns
# @Date:   2023-01-18 12:45:49
# @Last Modified by:   lshuns
# @Last Modified time: 2023-01-19 16:02:24

### get the median PSF value for each tile

import os
import re
import glob
import scipy
import shutil

import numpy as np 
import pandas as pd
from scipy.spatial import ConvexHull, convex_hull_plot_2d

# >>>>>> I/O
# the parent path where B_1_model_psf.py outputs are saved
outDir = '/disks/shear15/ssli/KiDS_noise_seeing/'

# >>>>>> workhorse
for i_label in ('i1', 'i2'):

    outDir_sub = os.path.join(outDir, f'psf2moffat_{i_label}')
    outpath = os.path.join(outDir_sub, f'moffat_{i_label}.csv')

    ## all the tiles with moffat measurements
    file_list = glob.glob(os.path.join(outDir_sub, 'KIDS_*.fits.psf.cat.moffat'))
    print('number of moffat files found', len(file_list))

    ## for saving final results
    data_final = pd.DataFrame(0, 
                        columns=["label", 'no', f"MeasSeeing_{i_label}", f"MeasBeta_{i_label}", f"seeing_e1_{i_label}"], 
                        index=np.arange(len(file_list)))
    data_final = data_final.astype({'label': str, 'no': 'int32'})

    ## loop over all files and estimate median using Convex Hull Peeling
    for i_file, file in enumerate(file_list):
        # the label
        data_final.loc[i_file, 'label'] = re.search(r'KIDS_(.*)_' + i_label, os.path.basename(file))[1]

        # +++ estimate median of Moffat parameters from Moffat fitting using Convex Hull Peeling
        # >>>> 1. read data
        try:
            cata = np.loadtxt(file)
        except ValueError:
            cata = []
            with open(file) as f:
                for line in f:
                    line = line.strip().split()
                    if len(line)==12:
                        cata.append(np.array(line).astype(float))
            cata = np.array(cata)
        # >>>> 2. clean outliers
        Nori = len(cata)
        ### only preserve 98% ## to remove outliers
        cata = cata[(cata[:, -1] > np.quantile(cata[:, -1], 0.01))&(cata[:, -1] < np.quantile(cata[:, -1], 0.99))\
                    &(cata[:, -2] > np.quantile(cata[:, -2], 0.01))&(cata[:, -2] < np.quantile(cata[:, -2], 0.99))]
        Nnew = len(cata)
        if Nori != Nnew:
            print(data_final.loc[i_file, 'label'], f'{Nori} >>> {Nnew}')
        data_final.loc[i_file, 'no'] = Nnew
        # >>>> 3. Convex Hull Peeling
        ## build points
        points = np.zeros((len(cata), 2))
        points[:, 0] = cata[:, -2]
        points[:, 1] = cata[:, -1]
        del cata
        ## peel hull
        while True:
            if len(points) <=3:
                break
            try:
                hull = ConvexHull(points)
            except scipy.spatial.qhull.QhullError:
                break
            # peel
            if len(points) == len(hull.vertices):
                break
            points = np.delete(points, hull.vertices, axis=0)
        ## the median
        if len(points) >20:
            print(points)
            raise Exception(f'too many points left {len(points)}')
        print('>>> number after peeling', len(points))
        data_final.loc[i_file, f'MeasSeeing_{i_label}'] = np.median(points[:, 0])
        data_final.loc[i_file, f'MeasBeta_{i_label}'] = np.median(points[:, 1])
        del points

        # +++ estimate median e from star catalogue
        ## read data
        cata = np.loadtxt(file.replace('.psf.cat.moffat', '.psf.cat'))
        ## the e defined as 1-b/a
        e_list = 1 - cata[:, 6]/cata[:, 5]
        del cata
        ## save the mean
        e_median =  np.median(e_list)
        del e_list
        data_final.loc[i_file, f'seeing_e1_{i_label}'] = e_median

        print('+++++ median seeing, beta, e', 
            data_final.loc[i_file, f'MeasSeeing_{i_label}'], 
            data_final.loc[i_file, f'MeasBeta_{i_label}'],
            data_final.loc[i_file, f'seeing_e1_{i_label}'])

    data_final.to_csv(outpath, index=False)
    del data_final
    print('saved to', outpath)


# number of moffat files found 1347
# saved to /disks/shear15/ssli/KiDS_noise_seeing/psf2moffat_i2/moffat_i2.csv
# Elapsed:9:03.61,User=478.364,System=14.169,CPU=90.6%.
