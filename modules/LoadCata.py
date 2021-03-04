# -*- coding: utf-8 -*-
# @Author: lshuns
# @Date:   2021-01-07 16:24:17
# @Last modified by:   lshuns
# @Last modified time: 2021-03-02, 20:51:17

### Everything about input catalogue

import logging
import numpy as np
import pandas as pd
from astropy.io import fits

logger = logging.getLogger(__name__)

def GalInfo(cata_pathfile, bands,
            id_name, primary_mag_name, mag_name_list,
            RaDec_names,
            shape_names,
            rng_seed=940120, mag_min_cut=0):
    """
    Get galaxy information from the input catalogue.

    Parameters
    ----------
    cata_pathfile : str
        Path to the input catalogue.
    bands : str
        A list of bands being queried, should be consistent with mag_name_list.
    id_name : str
        Column name for galaxy unique IDs.
    primary_mag_name : str
        Name for the primary magnitude used for pre-selection and save.
    mag_name_list : list of str
        A list of names for the magnitudes wanted.
    RaDec_names : list of str [RA, dec]
        Column names for position information.
    shape_names : list of str [Re, sersic_n, axis_ratio, PA, bulge_fraction, bulge_size, disk_size]
        Not all required, for those missed, simply feed 'none'.
    info_outfile (optional) : str (default: None)
        Path to the output of the simulated info.
    rng_seed (optional) : int (default: 940120)
        Seed for random number generator.
    mag_min_cut (optional) : float (default: 0)
        Min value of the primary band (the brightest galaxies will be simulated).
    """
    logger.info('Collect galaxy info from input catalogue...')
    logger.info(f'Query bands {bands}, the detection band {primary_mag_name}')
    logger.info(f"Magnitude brightest cut in {mag_min_cut}")

    # load file
    file_type = cata_pathfile[-3:]
    if file_type == 'csv':
        cata = pd.read_csv(cata_pathfile)
    elif file_type == 'her':
        cata = pd.read_feather(cata_pathfile)
    elif file_type == 'its':
        with fits.open(cata_pathfile) as hdul:
            cata = hdul[1].data
    else:
        raise Exception(f'Not supported file type! {cata_pathfile}')
    logger.info(f'gal_cata: {cata_pathfile}')

    # pre-selection based on magnitude
    ### only magnitude > mag_min_cut is simulated
    mask_mag = cata[primary_mag_name] > mag_min_cut
    cata = cata[mask_mag]
    if isinstance(cata, pd.DataFrame):
        cata.reset_index(drop=True, inplace=True)

    # random seed
    np.random.seed(rng_seed)

    # galaxy id
    index = cata[id_name]
    Ngal = len(index)
    logger.debug(f'Number of sources {Ngal}')

    # shape info
    name_Re, name_n, name_ar, name_PA, name_bf, name_bs, name_bar, name_bn, name_ds, name_dar = shape_names

    ## sersic profile
    try:
        Re = cata[name_Re]
    except KeyError:
        Re = np.zeros(Ngal)
    try:
        sersic_n = cata[name_n]
    except KeyError:
        sersic_n = np.zeros(Ngal)
    try:
        axis_ratios = cata[name_ar] # b/a
    except KeyError:
        axis_ratios = np.zeros(Ngal)

    ## commom
    try:
        position_angles = cata[name_PA]
    except KeyError:
        position_angles = np.random.uniform(low=-180., high=180., size=Ngal)

    ## bulge
    try:
        bulge_fractions = cata[name_bf] # bulge to total
    except KeyError:
        bulge_fractions = np.zeros(Ngal)
    try:
        bulge_Re = cata[name_bs] # arcsec
    except KeyError:
        bulge_Re = np.zeros(Ngal)
    try:
        bulge_axis_ratios = cata[name_bar] # arcsec
    except KeyError:
        bulge_axis_ratios = np.zeros(Ngal)
    try:
        bulge_n = cata[name_bn] # sersic n
    except KeyError:
        bulge_n = np.full(Ngal, 4.)

    ## disk
    try:
        disk_Re = cata[name_ds] # arcsec
    except KeyError:
        disk_Re = np.zeros(Ngal)
    try:
        disk_axis_ratios = cata[name_dar] # arcsec
    except KeyError:
        disk_axis_ratios = np.zeros(Ngal)

    # magnitudes
    mag_list = [np.array(cata[mag_name]).astype(float) for mag_name in mag_name_list]

    # sky position
    X_gals = cata[RaDec_names[0]] # degree
    Y_gals = cata[RaDec_names[1]] # degree

    # collect all information to dataframe
    data = [np.array(index).astype(int), np.array(X_gals).astype(float), np.array(Y_gals).astype(float),
                np.array(sersic_n).astype(float), np.array(Re).astype(float), np.array(axis_ratios).astype(float), np.array(position_angles).astype(float),
                np.array(bulge_fractions).astype(float), np.array(bulge_Re).astype(float), np.array(bulge_axis_ratios).astype(float), np.array(bulge_n).astype(float),
                np.array(disk_Re).astype(float), np.array(disk_axis_ratios).astype(float)] + mag_list
    data = np.transpose(data)
    name = ['index','RA','DEC',
                'sersic_n','Re','axis_ratios','position_angles',
                'bulge_fractions','bulge_Re','bulge_axis_ratios','bulge_n',
                'disk_Re','disk_axis_ratios'] + bands
    gals_info = pd.DataFrame(data=data, columns=name)
    gals_info = gals_info.astype({'index': int})
    logger.debug('Desired info collected to DataFrame.')

    return gals_info

def StarInfo(cata_pathfile, bands,
            id_name, primary_mag_name, mag_name_list,
            RaDec_names=None, mag_min_cut=0):
    """
    Get star information from the input catalogue.

    Parameters
    ----------
    cata_pathfile : str
        Path to the input catalogue.
    bands : str
        A list of bands being queried, should be consistent with mag_name_list.
    id_name : str
        Column name for galaxy unique IDs.
    primary_mag_name : str
        Name for the primary magnitude used for pre-selection and save.
    mag_name_list : list of str
        A list of names for the magnitudes wanted.
    RaDec_names (optional) : list of str [RA, dec] (default: None)
        Column names for position information.
        Not required, if stars will be randomly placed.
    mag_min_cut (optional) : float (default: 0)
        Min value of the primary band (the brightest stars will be simulated).
    """
    logger.info('Collect star info from input catalogue...')
    logger.info(f'Query bands {bands}, the primary band {primary_mag_name}')
    logger.info(f"Magnitude brightest cut in {mag_min_cut}")

    # load file
    file_type = cata_pathfile[-3:]
    if file_type == 'csv':
        cata = pd.read_csv(cata_pathfile)
    elif file_type == 'her':
        cata = pd.read_feather(cata_pathfile)
    elif file_type == 'its':
        with fits.open(cata_pathfile) as hdul:
            cata = hdul[1].data
    else:
        raise Exception(f'Not supported file type! {cata_pathfile}')
    logger.info(f'star_cata: {cata_pathfile}')

    # pre-selection based on magnitude
    ### only magnitude > mag_min_cut is simulated
    mask_mag = cata[primary_mag_name] > mag_min_cut
    cata = cata[mask_mag]
    if isinstance(cata, pd.DataFrame):
        cata.reset_index(drop=True, inplace=True)

    # unique id
    index = cata[id_name]
    Nstar = len(index)
    logger.debug(f'Number of stars {Nstar}')

    # magnitudes
    mag_list = [np.array(cata[mag_name]).astype(float) for mag_name in mag_name_list]

    # sky position
    if RaDec_names is not None:
        ## use true position
        X_stars = cata[RaDec_names[0]] # degree
        Y_stars = cata[RaDec_names[0]] # degree
    else:
        ## will be set late
        X_stars = np.zeros(Nstar)
        Y_stars = np.zeros(Nstar)

    # collect all information to dataframe
    data = [np.array(index).astype(int), np.array(X_stars).astype(float), np.array(Y_stars).astype(float)] + mag_list
    data = np.transpose(data)
    name = ['index','RA','DEC'] + bands
    stars_info = pd.DataFrame(data=data, columns=name)
    stars_info = stars_info.astype({'index': int})
    logger.debug('Desired info collected to DataFrame.')

    return stars_info

def NoiseInfo(cata_pathfile, bands, noise_psf_basenames, multiple_exposures_list=[], N_exposures=0):
    """
    Get star information from the input catalogue.

    Parameters
    ----------
    cata_pathfile : str
        Path to the input catalogue.
    bands : str
        A list of bands being queried.
    noise_psf_basenames : list of str [label, rms, seeing, MoffatBeta, psf_e1, psf_e2]
        A list of base names for noise and psf info.
        Not all required, if not available, simply set 'none'
    multiple_exposures_list : a list of boolean, optional (default: [])
        Use multiple noise info for multiple exposures or not, a list for MultiBand.
    N_exposures : int, optional (default: 0)
        Number of exposures, only used when multiple_exposures==True
    """
    logger.info('Collect observation info from input catalogue...')
    logger.info(f'Query bands {bands}')

    # load file
    file_type = cata_pathfile[-3:]
    if file_type == 'csv':
        cata = pd.read_csv(cata_pathfile)
    elif file_type == 'her':
        cata = pd.read_feather(cata_pathfile)
    elif file_type == 'its':
        with fits.open(cata_pathfile) as hdul:
            cata = hdul[1].data
    else:
        raise Exception(f'Not supported file type! {cata_pathfile}')
    logger.info(f'noise_cata: {cata_pathfile}')

    # noise and seeing info
    name_label, name_rms, name_seeing, name_MoffatBeta, name_psf_e1, name_psf_e2 = noise_psf_basenames

    noise_info = pd.DataFrame({'label': cata[name_label]})
    for i_band, band in enumerate(bands):
        try:
            multiple_exposures = multiple_exposures_list[i_band]
        except IndexError:
            multiple_exposures = False

        if multiple_exposures:
            for i_expo in range(N_exposures):
                noise_info.loc[:, f'rms_{band}_expo{i_expo}'] = np.array(cata[f'{name_rms}_{band}_expo{i_expo}']).astype(float)
                noise_info.loc[:, f'seeing_{band}_expo{i_expo}'] = np.array(cata[f'{name_seeing}_{band}_expo{i_expo}']).astype(float)
                try:
                    noise_info.loc[:, f'beta_{band}_expo{i_expo}'] = np.array(cata[f'{name_MoffatBeta}_{band}_expo{i_expo}']).astype(float)
                except KeyError:
                    # (if not provided using 3.5)
                    noise_info.loc[:, f'beta_{band}_expo{i_expo}'] = np.full(len(cata), 3.5)
                try:
                    noise_info.loc[:, f'psf_e1_{band}_expo{i_expo}'] = np.array(cata[f'{name_psf_e1}_{band}_expo{i_expo}']).astype(float)
                except KeyError:
                    # (if not provided assuming 0)
                    noise_info.loc[:, f'psf_e1_{band}_expo{i_expo}'] = np.zeros(len(cata))
                try:
                    noise_info.loc[:, f'psf_e2_{band}_expo{i_expo}'] = np.array(cata[f'{name_psf_e2}_{band}_expo{i_expo}']).astype(float)
                except KeyError:
                    # (if not provided assuming 0)
                    noise_info.loc[:, f'psf_e2_{band}_expo{i_expo}'] = np.zeros(len(cata))
        else:
            noise_info.loc[:, f'rms_{band}'] = np.array(cata[f'{name_rms}_{band}']).astype(float)
            noise_info.loc[:, f'seeing_{band}'] = np.array(cata[f'{name_seeing}_{band}']).astype(float)
            try:
                noise_info.loc[:, f'beta_{band}'] = np.array(cata[f'{name_MoffatBeta}_{band}']).astype(float)
            except KeyError:
                # (if not provided using 3.5)
                noise_info.loc[:, f'beta_{band}'] = np.full(len(cata), 3.5)
            try:
                noise_info.loc[:, f'psf_e1_{band}'] = np.array(cata[f'{name_psf_e1}_{band}']).astype(float)
            except KeyError:
                # (if not provided assuming 0)
                noise_info.loc[:, f'psf_e1_{band}'] = np.zeros(len(cata))
            try:
                noise_info.loc[:, f'psf_e2_{band}'] = np.array(cata[f'{name_psf_e2}_{band}']).astype(float)
            except KeyError:
                # (if not provided assuming 0)
                noise_info.loc[:, f'psf_e2_{band}'] = np.zeros(len(cata))

    logger.debug('Desired info collected to DataFrame.')

    return noise_info
