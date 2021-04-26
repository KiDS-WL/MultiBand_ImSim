# -*- coding: utf-8 -*-
# @Author: lshuns
# @Date:   2020-09-24 17:42:47
# @Last modified by:   ssli
# @Last modified time: 2021-04-26, 16:50:37

### cross-match catalogues based on object positions using KDTree

import os
import sys
import logging
import subprocess

import numpy as np
import pandas as pd
import scipy.spatial as sst

logger = logging.getLogger(__name__)

def KDTreeFunc(X1, X2, max_distance=np.inf, unique=True, k=1, leafsize=100):
    """
    Cross-match the values between X1 and X2 using a KD Tree.

    Parameters
    ----------
    X1 : array_like
        first dataset, shape(N1, D)
    X2 : array_like
        second dataset, shape(N2, D)
    max_distance : float (optional)
        maximum radius of search.  If no point is within the given radius,
        then inf will be returned.
    unique: bool (optional: True)
        return unique X2 index
    k : list of integer or integer
        The list of k-th nearest neighbors to return. If k is an integer it is treated as a list of [1, … k] (range(1, k+1)). Note that the counting starts from 1.

    Returns
    -------
    dist, ind: ndarrays
        The distance and index of the closest point in X2 to each point in X1
        Both arrays are length N1.
        Locations with no match are indicated by
        dist[i] = inf, ind[i] = N2
    """

    X1 = np.asarray(X1, dtype=float)
    X2 = np.asarray(X2, dtype=float)

    N1, D = X1.shape
    N2, D2 = X2.shape

    if D != D2:
        raise ValueError('Arrays must have the same second dimension')

    # build KDTree
    kdt = sst.cKDTree(X2, leafsize=leafsize)

    # query
    if unique:
        ## return one-one match
        ### original returns
        dist, ind = kdt.query(X1, k=k, distance_upper_bound=max_distance)
        if k==1:
            ### find duplicated matched
            u_ind, c_ind = np.unique(ind, return_counts=True)
            ind_duplicated = u_ind[c_ind>1]
            ### pick nearest
            for ind_tmp in ind_duplicated:
                ### mis-matched are passed
                if ind_tmp != len(X2):
                    ### which ind is duplicated
                    flag_duplicated = (ind == ind_tmp)
                    ### distance for duplicated ind
                    dist_tmp = dist[flag_duplicated]
                    dist_min = np.min(dist_tmp)
                    flag_wrong = (dist>dist_min)
                    ### only min dist is remained
                    ind[flag_duplicated & flag_wrong] = len(X2)
                    dist[flag_duplicated & flag_wrong] = np.inf
        else:
            ### return one list of entries by iterating all neighbors
            dist_final = np.full(len(dist), np.inf)
            ind_final = np.full(len(ind), len(X2))
            for i_k in range(k):
                ### the layer of neighbors
                dist_i = dist[:,i_k]
                ind_i = ind[:, i_k]
                ### copy to final list of entry (if it is empty)
                flag_empty = (dist_final==np.inf)
                dist_final = np.where(flag_empty, dist_i, dist_final)
                ind_final = np.where(flag_empty, ind_i, ind_final)

                ### find duplicated matched
                u_ind, c_ind = np.unique(ind_final, return_counts=True)
                ind_duplicated = u_ind[c_ind>1]
                ### pick nearest
                for ind_tmp in ind_duplicated:
                    ### mis-matched are passed
                    if ind_tmp != len(X2):
                        ### which ind is duplicated
                        flag_duplicated = (ind_final == ind_tmp)
                        ### distance for duplicated ind
                        dist_tmp = dist_final[flag_duplicated]
                        dist_min = np.min(dist_tmp)
                        flag_wrong = (dist_final>dist_min)
                        ### only min dist is remained
                        ind_final[flag_duplicated & flag_wrong] = len(X2)
                        dist_final[flag_duplicated & flag_wrong] = np.inf
            ind = ind_final
            dist = dist_final

    else:
        ## simply return query results
        dist, ind = kdt.query(X1, k=k, distance_upper_bound=max_distance)

    return dist, ind

def run_position2id(input_cata, detec_cata, id_list, position_list, mag_list,
                    outDir=None, basename=None, save_matched=False, save_false=False, save_missed=False,
                    dmag_max=0.5, r_max=0.5/3600., k=4, running_info=True):
    """
    Run cross-match based on positions.
        output files only include the unique id from both catalogues (and distance).
    """

    if running_info:
        logger.info("Running position-based cross match with KD Tree...")

    # id
    id_input = np.array(input_cata[id_list[0]], dtype=int)
    id_detec = np.array(detec_cata[id_list[1]], dtype=int)

    # magnitude
    mag_input = np.array(input_cata[mag_list[0]], dtype=float)
    mag_detec = np.array(detec_cata[mag_list[1]], dtype=float)

    # location
    xy_input = np.empty((len(id_input), 2), dtype=float)
    xy_input[:, 0] = np.array(input_cata[position_list[0][0]], dtype=float)
    xy_input[:, 1] = np.array(input_cata[position_list[0][1]], dtype=float)
    ##
    xy_detec = np.empty((len(id_detec), 2), dtype=float)
    xy_detec[:, 0] = np.array(detec_cata[position_list[1][0]], dtype=float)
    xy_detec[:, 1] = np.array(detec_cata[position_list[1][1]], dtype=float)

    if running_info:
        logger.info('Number of input {:}'.format(len(id_input)))
        logger.info('Number of detections {:}'.format(len(id_detec)))

    # KDTree match
    dist, ind = KDTreeFunc(xy_detec, xy_input, max_distance=r_max, unique=True, k=k, leafsize=100)
    ##
    flag_matched = ind<len(id_input)
    flag_false = np.invert(flag_matched)
    ind_matched = ind[flag_matched]

    # matched catalogue
    tmp_matched_cata = pd.DataFrame({
                    'id_detec': id_detec[flag_matched],
                    'id_input': id_input[ind_matched],
                    'distance_CM': dist[flag_matched],
                    'dmag_CM': np.abs(mag_input[ind_matched]-mag_detec[flag_matched])
                    })
    ## magnitude cut
    matched_cata = tmp_matched_cata[(tmp_matched_cata['dmag_CM']<dmag_max)].copy()
    matched_cata.reset_index(drop=True, inplace=True)
    if running_info:
        logger.info('Number of matched {:}'.format(len(matched_cata)))
        logger.info('matched/detections {:}'.format(len(matched_cata)/len(id_detec)))

    # false detection
    id_false = id_detec[flag_false]
    false_cata = pd.DataFrame({
                    'id_detec': id_false,
                    'id_input': np.full(len(id_false), -999).astype(int),
                    'distance_CM': np.full(len(id_false), -999).astype(int),
                    'dmag_CM': np.full(len(id_false), -999.).astype(float)
                    })
    ## dmag does not meet requirement
    tmp_cata = tmp_matched_cata[(tmp_matched_cata['dmag_CM']>=dmag_max)].copy()
    false_cata = pd.concat([false_cata, tmp_cata])
    false_cata.reset_index(drop=True, inplace=True)
    if running_info:
        logger.info('Number of false detections {:}'.format(len(false_cata)))
        logger.info('false_detec/detections {:}'.format(len(false_cata)/len(id_detec)))

    # missed inputs
    mask_tmp = np.ones(len(id_input), np.bool)
    mask_tmp[ind_matched] = 0
    id_miss = id_input[mask_tmp]
    miss_cata = pd.DataFrame({
                    'id_detec': np.full(len(id_miss), -999).astype(int),
                    'id_input': id_miss,
                    'distance_CM': np.full(len(id_miss), -999).astype(float),
                    'dmag_CM': np.full(len(id_miss), -999).astype(float)
                    })
    ## dmag does not meet requirement
    tmp_cata = tmp_matched_cata[(tmp_matched_cata['dmag_CM']>=dmag_max)].copy()
    miss_cata = pd.concat([miss_cata, tmp_cata])
    miss_cata.reset_index(drop=True, inplace=True)
    if running_info:
        logger.info('Number of missed inputs {:}'.format(len(miss_cata)))
        logger.info('miss_input/input {:}'.format(len(miss_cata)/len(id_input)))

    # save results if specified
    if save_matched:
        outfile = os.path.join(outDir, basename+'_matched.feather')
        matched_cata.to_feather(outfile)
        logger.info('Matched detections saved to {:}'.format(outfile))
    if save_false:
        outDir_miss_false = os.path.join(outDir, 'miss_false')
        if not os.path.exists(outDir_miss_false):
            os.mkdir(outDir_miss_false)
        outfile = os.path.join(outDir_miss_false, basename+'_false.feather')
        false_cata.to_feather(outfile)
        logger.info('False detections saved to {:}'.format(outfile))
    if save_missed:
        outDir_miss_false = os.path.join(outDir, 'miss_false')
        if not os.path.exists(outDir_miss_false):
            os.mkdir(outDir_miss_false)
        outfile = os.path.join(outDir_miss_false, basename+'_miss.feather')
        miss_cata.to_feather(outfile)
        logger.info('Missed inputs saved to {:}'.format(outfile))

    if running_info:
        logger.info("Finised position-based cross match with KD Tree.")
    return matched_cata, false_cata, miss_cata
