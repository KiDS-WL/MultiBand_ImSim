# -*- coding: utf-8 -*-
# @Author: lshuns
# @Date:   2026-10-19

### cross-match catalogues based on object positions using KDTree, with optimal one-to-one assignment
###     same interface and outputs as CrossMatch.py, plus blend information
###     differences from CrossMatch.py:
###         1. all input-detection pairs within the max distance are candidates
###             (CrossMatch.py only keeps the k nearest inputs of each detection)
###         2. the candidates are made one-to-one by solving the assignment problem:
###             mag_closest=True: minimise the total magnitude difference (distance only breaks ties),
###                 where each source can also stay unmatched at a cost of 1 mag,
###                 so pairs with |dmag| >= 2 mag are never kept, whether they compete with others or not;
###                 pairs without a usable magnitude (non-finite or a sentinel like 99) are matched by
###                 position only: kept when alone, but losing to any pair with a usable magnitude
###             mag_closest=False: maximise the number of matches, then minimise the total distance
###             (CrossMatch.py gives each detection its nearest input, and uses the magnitude only
###             when several detections claim the same input, so it is never used for unresolved blends)
###         3. blend information is saved as additional columns (*_CM), from friends-of-friends groups
###             linking all sources (both catalogues) within the max distance of each other

import os
import logging

import numpy as np
import pandas as pd
import scipy.spatial as sst
from scipy.optimize import linear_sum_assignment
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components

import galsim

logger = logging.getLogger(__name__)

# cost (mag) for a source to stay unmatched, used when mag_closest=True
##    a pair is only kept if |dmag| < 2*UNMATCHED_COST_MAG
##    in the SKiLLS-like tests, no correct match had |dmag| > 1.93 (99.9% < 1.02),
##    while pairs above 2 mag were forced or wrong matches
UNMATCHED_COST_MAG = 1.0
# magnitudes at or above this are taken as sentinels for failed measurements (e.g. 99)
MAG_INVALID = 90.

def KDTreePairsFunc(X1, X2, max_distance, leafsize=100):
    """
    Find all the pairs between X1 and X2 within max_distance using a KD Tree.

    Parameters
    ----------
    X1 : array_like
        first dataset, shape(N1, D)
    X2 : array_like
        second dataset, shape(N2, D)
    max_distance : float
        maximum radius of search
    leafsize : int (optional: 100)
        leaf size of the KD Tree

    Returns
    -------
    ind1, ind2, dist: ndarrays
        index in X1, index in X2, and distance of each pair
    """

    X1 = np.asarray(X1, dtype=float)
    X2 = np.asarray(X2, dtype=float)

    N1, D = X1.shape
    N2, D2 = X2.shape

    if D != D2:
        raise ValueError('Arrays must have the same second dimension')

    # build KDTree
    kdt1 = sst.cKDTree(X1, leafsize=leafsize)
    kdt2 = sst.cKDTree(X2, leafsize=leafsize)

    # all pairs within the max distance
    pairs = kdt1.sparse_distance_matrix(kdt2, max_distance, output_type='ndarray')

    return pairs['i'].astype(int), pairs['j'].astype(int), pairs['v'].astype(float)

def OptimalPairFunc(id_input_pair, id_detec_pair, cost_pair, N_input, N_detec, unmatched_cost=None):
    """
    Make one-to-one matches from candidate pairs, by solving the assignment problem.
        unmatched_cost=None (default):
            The number of matches is maximised first, then the total cost of the matches is minimised.
            Picking the best pair first instead can leave a source unmatched:
                detection A takes the input of detection B because it is the closest one in magnitude,
                then B has no candidate left, although A could have been matched to its own input.
            But it can also force wrong matches when a source has no true counterpart
                (e.g. a spurious detection, a shredded piece, or an input blended away).
        unmatched_cost given:
            Each source can also stay unmatched, at this cost,
            and the total cost (matches + unmatched sources) is minimised.
            So a pair is only formed if it is cheaper than leaving both its sources unmatched,
            i.e. cost < 2*unmatched_cost, whether it competes with others or not.

    Parameters
    ----------
    id_input_pair, id_detec_pair : array_like
        row index of the input and detection source in each candidate pair
        each (input, detection) pair should appear only once
    cost_pair : array_like
        cost of each pair (>=0, smaller is better)
    N_input, N_detec : int
        number of sources in the input and detection catalogues
    unmatched_cost : float (optional: None)
        cost for a source to stay unmatched, see above

    Returns
    -------
    flag_good: ndarray[bool]
        selected pairs
    """

    N_pair = len(cost_pair)
    flag_good = np.zeros(N_pair, dtype=bool)
    if N_pair == 0:
        return flag_good

    # pairs with no competition: the input and the detection have no other candidate
    n_per_input = np.bincount(id_input_pair, minlength=N_input)
    n_per_detec = np.bincount(id_detec_pair, minlength=N_detec)
    flag_simple = (n_per_input[id_input_pair] == 1) & (n_per_detec[id_detec_pair] == 1)
    if unmatched_cost is None:
        flag_good[flag_simple] = True
    else:
        ## same rule as in the competing sets: only kept if cheaper than both unmatched
        flag_good[flag_simple] = (cost_pair[flag_simple] < 2. * unmatched_cost)
    idx_conf = np.flatnonzero(~flag_simple)
    if len(idx_conf) == 0:
        return flag_good

    # the others: split into connected sets of competing pairs
    _, inp_c = np.unique(id_input_pair[idx_conf], return_inverse=True)
    _, det_c = np.unique(id_detec_pair[idx_conf], return_inverse=True)
    n_i = np.amax(inp_c) + 1
    n_d = np.amax(det_c) + 1
    graph = coo_matrix((np.ones(len(idx_conf)), (inp_c, n_i + det_c)), shape=(n_i + n_d, n_i + n_d))
    _, comp_node = connected_components(graph, directed=False)
    comp_pair = comp_node[inp_c]
    order = np.argsort(comp_pair, kind='stable')
    bounds = np.flatnonzero(np.diff(comp_pair[order])) + 1

    # solve each set
    for sel in np.split(order, bounds):
        _, rows = np.unique(inp_c[sel], return_inverse=True)
        _, cols = np.unique(det_c[sel], return_inverse=True)
        cost = cost_pair[idx_conf[sel]]
        n_r = np.amax(rows) + 1
        n_c = np.amax(cols) + 1
        pair_mat = np.full((n_r, n_c), -1, dtype=np.int64)
        pair_mat[rows, cols] = idx_conf[sel]
        if unmatched_cost is None:
            ## non-candidates cost more than all candidates together,
            ##    so the number of real matches is maximised first
            cost_mat = np.full((n_r, n_c), np.sum(cost) + 1.)
            cost_mat[rows, cols] = cost
        else:
            ## augmented problem: each input (row) and detection (column) has its own "unmatched" slot
            ##      [ pairs                 | input i unmatched ]
            ##      [ detection j unmatched | 0                 ]
            ##    forbidden entries cost more than leaving everything unmatched
            forbid = (n_r + n_c) * unmatched_cost + np.sum(cost) + 1.
            cost_mat = np.full((n_r + n_c, n_c + n_r), forbid)
            cost_mat[rows, cols] = cost
            cost_mat[np.arange(n_r), n_c + np.arange(n_r)] = unmatched_cost
            cost_mat[n_r + np.arange(n_c), np.arange(n_c)] = unmatched_cost
            cost_mat[n_r:, n_c:] = 0.
        row_best, col_best = linear_sum_assignment(cost_mat)
        flag_pair = (row_best < n_r) & (col_best < n_c)
        chosen = pair_mat[row_best[flag_pair], col_best[flag_pair]]
        flag_good[chosen[chosen >= 0]] = True

    return flag_good

def GroupFunc(X_input, X_detec, max_distance, leafsize=100):
    """
    Friends-of-friends groups of all the sources (both catalogues),
        linking any two sources within max_distance.

    Parameters
    ----------
    X_input, X_detec : array_like
        positions of the input and detection catalogues, shape(N, D)
    max_distance : float
        linking length
    leafsize : int (optional: 100)
        leaf size of the KD Tree

    Returns
    -------
    group_input, group_detec: ndarray[int]
        group id of each source in the input and detection catalogues
    """

    N_input = len(X_input)
    X_all = np.concatenate([np.asarray(X_input, dtype=float), np.asarray(X_detec, dtype=float)])
    N_all = len(X_all)

    pairs = sst.cKDTree(X_all, leafsize=leafsize).query_pairs(max_distance, output_type='ndarray')
    graph = coo_matrix((np.ones(len(pairs)), (pairs[:, 0], pairs[:, 1])), shape=(N_all, N_all))
    _, group = connected_components(graph, directed=False)

    return group[:N_input], group[N_input:]

def run_position2id(input_cata, detec_cata, id_list, position_list, mag_list,
                    outDir=None, basename=None, save_matched=False, save_false=False, save_missed=False,
                    r_max=0.6/3600., k=4, mag_closest=True, running_info=True,
                    useTan=False, pixel_scale=None, r_max_pixel=3):
    """
    Run cross-match based on positions.
        output files only include the unique id from both catalogues (and distance),
        plus blend information.

    All input-detection pairs within the max distance are candidates, made one-to-one by
        solving the assignment problem (see OptimalPairFunc):
        mag_closest=True: the total magnitude difference is minimised, and each source can
            stay unmatched at a cost of UNMATCHED_COST_MAG, so a pair is only kept if
            |dmag| < 2*UNMATCHED_COST_MAG (= 2 mag), whether it competes with others or not.
            This avoids forcing wrong matches on sources without a true counterpart
            (e.g. a spurious detection, a shredded piece, or an input blended away).
            Pairs without a usable magnitude (non-finite, or >= MAG_INVALID like 99) are matched
            by position only: kept when alone, but losing to any pair with a usable magnitude.
        mag_closest=False: the number of matches is maximised, then the total distance minimised
            (positions alone cannot tell a source without counterpart).

    Two options:
        useTan=True (recommend): mapping sky position to image position, then matching in image
            require:
                r_max_pixel defined in pixel
                pixel_scale: the corresponding image pixel
        useTan=False: directly match using sky position
            require:
                r_max defined in deg

    mag_closest:
        True (default): choose the magnitude-closest one in case of multiple matching
        False: choose the distance-closest one in case of multiple matching

    k:
        not used (all the candidates within the max distance are kept),
        kept for compatibility with CrossMatch.run_position2id

    Blend information (columns added to the outputs)
    ----------
        group_id_CM: friends-of-friends group id, linking all sources (both catalogues)
            within the max distance of each other
        n_input_group_CM: number of input sources in the group
        n_detec_group_CM: number of detections in the group
        n_input_near_CM: number of input sources within the max distance of the detection
            (matched and false detections only)
        n_detec_near_CM: number of detections within the max distance of the input
            (matched and missed inputs only)
        blend_flag_CM: bit flag from the group
            0: isolated (at most one input and one detection in the group)
            1: more than one input in the group (e.g. unresolved blend)
            2: more than one detection in the group (e.g. shredded object)
            3: both
    """

    if running_info:
        logger.info("Running position-based cross match with KD Tree and optimal assignment...")

    # id
    id_input = np.array(input_cata[id_list[0]], dtype=int)
    id_detec = np.array(detec_cata[id_list[1]], dtype=int)
    N_input = len(id_input)
    N_detec = len(id_detec)

    # magnitude
    mag_input = np.array(input_cata[mag_list[0]], dtype=float)
    mag_detec = np.array(detec_cata[mag_list[1]], dtype=float)

    # location
    if useTan:
        logger.info('CrossMatch using TAN projected distance')
        r_max = r_max_pixel
        col_r = 'distance_pixel_CM'
        logger.info(f'   max distance {r_max} pixel')
        if pixel_scale is None:
            raise Exception('need pixel_scale for TAN projection')

        # the coordinates provided
        RA_input = np.array(input_cata[position_list[0][0]], dtype=float)
        DEC_input = np.array(input_cata[position_list[0][1]], dtype=float)
        ##
        RA_detec = np.array(detec_cata[position_list[1][0]], dtype=float)
        DEC_detec = np.array(detec_cata[position_list[1][1]], dtype=float)

        # build the wcs
        ## Linear projection matrix
        dudx = pixel_scale/3600. # degree
        dudy = 0.
        dvdx = 0.
        dvdy = pixel_scale/3600. # degree
        ## Reference pixel in image
        origin_ima = galsim.PositionI(x=0, y=0)
        ## Reference point in wcs
        world_origin = galsim.CelestialCoord(ra=np.amin(RA_input)*galsim.degrees,
                                                dec=np.amin(DEC_input)*galsim.degrees)
        ## the wcs
        wcs_affine = galsim.AffineTransform(dudx, dudy, dvdx, dvdy, origin=origin_ima)
        wcs = galsim.TanWCS(wcs_affine, world_origin, units=galsim.degrees)

        # transfer ra dec to x y
        xy_input = np.empty((N_input, 2), dtype=float)
        xy_input[:, 0], xy_input[:, 1] = wcs.toImage(RA_input, DEC_input, units=galsim.degrees)
        del RA_input, DEC_input
        xy_detec = np.empty((N_detec, 2), dtype=float)
        xy_detec[:, 0], xy_detec[:, 1] = wcs.toImage(RA_detec, DEC_detec, units=galsim.degrees)
        del RA_detec, DEC_detec

    else:
        logger.info('CrossMatch using provided distance')
        logger.info(f'   max distance {r_max*3600} arcsec')
        col_r = 'distance_deg_CM'
        xy_input = np.empty((N_input, 2), dtype=float)
        xy_input[:, 0] = np.array(input_cata[position_list[0][0]], dtype=float)
        xy_input[:, 1] = np.array(input_cata[position_list[0][1]], dtype=float)
        ##
        xy_detec = np.empty((N_detec, 2), dtype=float)
        xy_detec[:, 0] = np.array(detec_cata[position_list[1][0]], dtype=float)
        xy_detec[:, 1] = np.array(detec_cata[position_list[1][1]], dtype=float)

    if running_info:
        logger.info('Number of input {:}'.format(N_input))
        logger.info('Number of detections {:}'.format(N_detec))

    # KDTree candidates
    if (N_input > 0) and (N_detec > 0):
        i_input_pair, i_detec_pair, dist_pair = KDTreePairsFunc(xy_input, xy_detec, r_max, leafsize=100)
        group_input, group_detec = GroupFunc(xy_input, xy_detec, r_max, leafsize=100)
    else:
        i_input_pair = np.array([], dtype=int)
        i_detec_pair = np.array([], dtype=int)
        dist_pair = np.array([], dtype=float)
        group_input = np.arange(N_input)
        group_detec = np.arange(N_input, N_input + N_detec)
    del xy_input, xy_detec

    # number of close candidates
    n_input_near = np.bincount(i_detec_pair, minlength=N_detec)
    n_detec_near = np.bincount(i_input_pair, minlength=N_input)

    # group sizes
    N_group = np.amax(np.concatenate([group_input, group_detec, [-1]])) + 1
    n_input_group = np.bincount(group_input, minlength=N_group)
    n_detec_group = np.bincount(group_detec, minlength=N_group)

    # one-to-one
    if mag_closest:
        ## magnitude difference, plus a tiny distance term to break ties
        ##    (at most 1e-3 mag at the max distance)
        dmag_pair = np.abs(mag_input[i_input_pair] - mag_detec[i_detec_pair])
        cost_pair = dmag_pair + 1e-3 * dist_pair / r_max
        ## pairs without a usable magnitude (non-finite or a sentinel like 99) are matched
        ##    by position only: just under the gate, so they are kept when alone,
        ##    but lose to any pair with a usable magnitude
        flag_bad_mag = ~np.isfinite(dmag_pair) \
                        | (np.abs(mag_input[i_input_pair]) >= MAG_INVALID) \
                        | (np.abs(mag_detec[i_detec_pair]) >= MAG_INVALID)
        cost_pair[flag_bad_mag] = 2. * UNMATCHED_COST_MAG - 1e-3 + 1e-4 * dist_pair[flag_bad_mag] / r_max
        if running_info and np.any(flag_bad_mag):
            logger.info(f'   {np.sum(flag_bad_mag)} candidate pairs without usable magnitudes, matched by position only')
        del dmag_pair, flag_bad_mag
        unmatched_cost = UNMATCHED_COST_MAG
    else:
        ## positions alone cannot tell a source without counterpart, so no unmatched option:
        ##    the number of matches is maximised, then the total distance minimised
        cost_pair = dist_pair
        unmatched_cost = None
    flag_good = OptimalPairFunc(i_input_pair, i_detec_pair, cost_pair, N_input, N_detec,
                                unmatched_cost=unmatched_cost)
    i_input_matched = i_input_pair[flag_good]
    i_detec_matched = i_detec_pair[flag_good]
    dist_matched = dist_pair[flag_good]
    del i_input_pair, i_detec_pair, dist_pair, flag_good, cost_pair

    ## sort by detection order, as in CrossMatch.py
    order = np.argsort(i_detec_matched)
    i_input_matched = i_input_matched[order]
    i_detec_matched = i_detec_matched[order]
    dist_matched = dist_matched[order]
    del order

    # blend flags
    blend_flag_group = (n_input_group > 1).astype(int) + 2*(n_detec_group > 1).astype(int)

    # matched catalogue
    matched_cata = pd.DataFrame({
                    'id_detec': id_detec[i_detec_matched],
                    'id_input': id_input[i_input_matched],
                    col_r: dist_matched,
                    'dmag_CM': (mag_input[i_input_matched]-mag_detec[i_detec_matched]),
                    'group_id_CM': group_detec[i_detec_matched],
                    'n_input_group_CM': n_input_group[group_detec[i_detec_matched]],
                    'n_detec_group_CM': n_detec_group[group_detec[i_detec_matched]],
                    'n_input_near_CM': n_input_near[i_detec_matched],
                    'n_detec_near_CM': n_detec_near[i_input_matched],
                    'blend_flag_CM': blend_flag_group[group_detec[i_detec_matched]]
                    })
    if running_info:
        logger.info('Number of matched {:}'.format(len(matched_cata)))
        logger.info('matched/detections {:}'.format(len(matched_cata)/max(N_detec, 1)))
        logger.info('blended matched (blend_flag_CM>0) {:}'.format(np.sum(matched_cata['blend_flag_CM'].values>0)))

    # false detection
    mask_tmp = np.ones(N_detec, bool)
    mask_tmp[i_detec_matched] = 0
    false_cata = pd.DataFrame({
                    'id_detec': id_detec[mask_tmp],
                    'group_id_CM': group_detec[mask_tmp],
                    'n_input_group_CM': n_input_group[group_detec[mask_tmp]],
                    'n_detec_group_CM': n_detec_group[group_detec[mask_tmp]],
                    'n_input_near_CM': n_input_near[mask_tmp],
                    'blend_flag_CM': blend_flag_group[group_detec[mask_tmp]]
                    })
    if running_info:
        logger.info('Number of false detections {:}'.format(len(false_cata)))
        logger.info('false_detec/detections {:}'.format(len(false_cata)/max(N_detec, 1)))

    # missed inputs
    mask_tmp = np.ones(N_input, bool)
    mask_tmp[i_input_matched] = 0
    miss_cata = pd.DataFrame({
                    'id_input': id_input[mask_tmp],
                    'group_id_CM': group_input[mask_tmp],
                    'n_input_group_CM': n_input_group[group_input[mask_tmp]],
                    'n_detec_group_CM': n_detec_group[group_input[mask_tmp]],
                    'n_detec_near_CM': n_detec_near[mask_tmp],
                    'blend_flag_CM': blend_flag_group[group_input[mask_tmp]]
                    })
    if running_info:
        logger.info('Number of missed inputs {:}'.format(len(miss_cata)))
        logger.info('miss_input/input {:}'.format(len(miss_cata)/max(N_input, 1)))
    del mask_tmp

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
        logger.info("Finised position-based cross match with KD Tree and optimal assignment.")
    return matched_cata, false_cata, miss_cata
