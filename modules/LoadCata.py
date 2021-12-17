# -*- coding: utf-8 -*-
# @Author: lshuns
# @Date:   2021-01-07 16:24:17
# @Last modified by:   lshuns
# @Last modified time: 2021-06-10, 18:17:08

### Everything about input catalogue

import os
import logging

import numpy as np
import pandas as pd

from astropy.io import fits
from astropy.table import Table

# let us ignore pandas PerformanceWarning for now
## as it is not the bottleneck, and it is hard to improve now...
### if you are curious, it is about the `frame.insert`,
###     basically I call df.loc[:, new_col]=new_val too many times...
#### it is recommended to pd.concat in the end to join all columns
import warnings
warnings.simplefilter(action='ignore', category=pd.errors.PerformanceWarning)

logger = logging.getLogger(__name__)

def GalInfo(cata_pathfile, bands,
            id_name, primary_mag_name, mag_name_list,
            RaDec_names,
            shape_names,
            z_name,
            mag_cut=[], size_cut=[],
            g_columns=None):
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
    z_name : str
        column name for redshift for saving.
    info_outfile (optional) : str (default: None)
        Path to the output of the simulated info.
    mag_cut (optional) : a list of float (default: [])
        Min & Max value of the primary band (the brightest & faintest galaxies will be simulated).
    size_cut (optional) : a list of float (default: [])
        Min & Max value of the effective radius in arcsec (the smallest & largest galaxies will be simulated).
    g_columns (optional) : a list of str (default: None)
        column names to the input shear, only provided for variable shear setup

    Returns
    -------
    gals_info : DataFrame
        galaxies info
    """

    logger.info('Collect galaxy info from input catalogue...')
    logger.info(f'Query bands {mag_name_list}, the detection band {primary_mag_name}')
    if mag_cut:
        logger.info(f"Magnitude cut in: {mag_cut}")
    if size_cut:
        logger.info(f"size cut in: {size_cut}")

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
    ### only magnitude within mag_cut is simulated
    if mag_cut:
        mask_tmp = (cata[primary_mag_name] >= mag_cut[0]) & (cata[primary_mag_name] <= mag_cut[1])
        cata = cata[mask_tmp]
        del mask_tmp
        if isinstance(cata, pd.DataFrame):
            cata.reset_index(drop=True, inplace=True)

    # galaxy id
    index = cata[id_name]
    Ngal = len(index)
    logger.debug(f'Number of sources {Ngal}')

    # name for the shape info
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
    position_angles = cata[name_PA]

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
    if z_name is not None:
        z_gals = cata[z_name] # for saving only
    else:
        z_gals = np.full(Ngal, -999)

    # collect all information to dataframe
    data = [np.array(index).astype(int), np.array(X_gals).astype(float), np.array(Y_gals).astype(float), np.array(z_gals).astype(float),
                np.array(sersic_n).astype(float), np.array(Re).astype(float), np.array(axis_ratios).astype(float), np.array(position_angles).astype(float),
                np.array(bulge_fractions).astype(float), np.array(bulge_Re).astype(float), np.array(bulge_axis_ratios).astype(float), np.array(bulge_n).astype(float),
                np.array(disk_Re).astype(float), np.array(disk_axis_ratios).astype(float)] + mag_list
    data = np.transpose(data)
    name = ['index','RA','DEC','redshift',
                'sersic_n','Re','axis_ratio','position_angle',
                'bulge_fraction','bulge_Re','bulge_axis_ratio','bulge_n',
                'disk_Re','disk_axis_ratio'] + bands
    gals_info = pd.DataFrame(data=data, columns=name)
    gals_info = gals_info.astype({'index': int})

    # input shear if provided
    if g_columns is not None:
        gals_info.loc[:, 'gamma1'] = np.array(cata[g_columns[0]]).astype(float)
        gals_info.loc[:, 'gamma2'] = np.array(cata[g_columns[1]]).astype(float)

    # select based on size
    ### only size within size_cut is simulated
    if size_cut:
        mask_tmp = (gals_info['Re'] >= size_cut[0]) & (gals_info['Re'] <= size_cut[1]) &\
                    (gals_info['bulge_Re'] >= size_cut[0]) & (gals_info['bulge_Re'] <= size_cut[1]) &\
                    (gals_info['disk_Re'] >= size_cut[0]) & (gals_info['disk_Re'] <= size_cut[1])
        gals_info = gals_info[mask_tmp]
        del mask_tmp
        gals_info.reset_index(drop=True, inplace=True)

    logger.debug('Desired info collected to DataFrame.')

    return gals_info

def GalInfo_adjust4casual(gals_info, 
                    qbin_band='r', qbin_mag=25., 
                    Nqbins=20, Frac_careful=0.01, rng_seed=512):
    """
    Adjust gals_info for casual mode (sample galaxies in faint end)

    Parameters
    ----------
    gals_info : DataFrame
        original galaxies information.
    qbin_band : str, optional (default: 'r')
        band used for quantile binning
    qbin_mag : float, optional (default: 25.)
        up to which magnitude galaxies are sampled
    Nqbins : int, optional (default: 20)
        Number of quantile bins for packing galaxies
    Frac_careful : float, optional (default: 0.01)
        Fraction of carefully simulated galaxies in each quantile bin
    rng_seed : int, optional (default: 512)
        seed for random sampling

    Returns
    -------
    gals_info_careful : DataFrame
        galaxies need careful simulation
    gals_info_casual : DataFrame
        galaxies can be casually simulated
    """

    logger.info('Adjust galaxy info for casual mode...')
    logger.info(f'bin band {qbin_band}, cut mag {qbin_mag}')
    logger.info(f'number bin {Nqbins}, Fraction of seed galaxies {Frac_careful}')

    # split the galaxies to careful and casual
    mask_casual = gals_info[qbin_band].values > qbin_mag
    gals_info_casual = gals_info[mask_casual].copy()
    gals_info_casual.reset_index(drop=True, inplace=True)
    gals_info_careful = gals_info[~mask_casual].copy()
    gals_info_careful.reset_index(drop=True, inplace=True)
    del mask_casual, gals_info

    # bin the causual sample with Nqbins provided
    gals_info_casual.loc[:, 'i_qbin'] = pd.qcut(gals_info_casual[qbin_band], Nqbins, labels=False)
    info_grouped = gals_info_casual.groupby(by=['i_qbin'])
    del gals_info_casual
    
    # loop over the bins and draw galaxies
    sampled_groups = []
    for name, group in info_grouped:

        # randomly sample the seed galaxies
        Ntot_tmp = len(group)
        ## seed galaxies
        gals0 = group.sample(frac=Frac_careful, random_state=rng_seed+name, ignore_index=True)
        Ncareful = len(gals0)
        if Ncareful < 10:
            raise Exception(f'Number of seed galaxies {Ncareful} < 10, think twice about your Frac_careful!')
        ### positions and tag are not used
        gals0.drop(columns=['RA', 'DEC'], inplace=True)
        ### index are renamed
        gals0.rename(columns={'index': 'index_seedGal'}, inplace=True)
        ## randomly sample from seed galaxies
        gals1 = gals0.sample(n=Ntot_tmp - Ncareful, random_state=rng_seed*36+name, replace=True)

        # assign to the group
        group = group[['index', 'RA', 'DEC']]
        group.reset_index(drop=True, inplace=True)
        group = pd.concat([group, pd.concat([gals0, gals1], ignore_index=True)], axis=1)
        sampled_groups.append(group)
        del gals0, gals1, group

    del info_grouped

    # concat to get the casual cata
    gals_info_casual = pd.concat(sampled_groups)

    return gals_info_careful, gals_info_casual

def StarInfo(cata_pathfile, bands,
            id_name, primary_mag_name, mag_name_list,
            RaDec_names=None, mag_cut=[]):
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
    mag_cut (optional) : a list of float (default: [])
        Min & Max value of the primary band (the brightest & faintest stars will be simulated).

    Returns
    -------
    stars_info : DataFrame
        stars info
    """

    logger.info('Collect star info from input catalogue...')
    logger.info(f'Query bands {mag_name_list}, the primary band {primary_mag_name}')
    if mag_cut:
        logger.info(f"Magnitude cut in: {mag_cut}")

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
    ### only magnitude within mag_cut is simulated
    if mag_cut:
        mask_tmp = (cata[primary_mag_name] >= mag_cut[0]) & (cata[primary_mag_name] <= mag_cut[1])
        cata = cata[mask_tmp]
        del mask_tmp
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

def NoiseInfo(cata_pathfile, bands, 
                only_labels, psf_type_list=['moffat'],
                noise_psf_basenames=None,
                label_basename=None, noise_basenames=None, 
                psf_basenames_moffat=None, psf_basenames_airy=None, 
                id_basenames=None,
                multiple_exposures_list=None, N_exposures_list=None, 
                file4varChips=None, varChips_list=None, N_chips_list=None,
                psf_PixelIma_dir=None):
    """
    Get noise information from the noise catalogue.

    Parameters
    ----------
    cata_pathfile : str
        Path to the input catalogue.
    bands : str
        A list of bands being queried.
    only_labels : bool
        only return labels or not
    psf_type_list : list of str, ['moffat', 'airy', 'pixelima']
        Type of PSF profile.
    noise_psf_basenames (deprecated) : list of str [label, rms, seeing, MoffatBeta, psf_e1, psf_e2]
        A list of base names for noise and psf info.
        Not all required, if not available, simply set 'none'
    label_basename : str
        column name for the noise label.
    noise_basenames : a list of str 
        A list of base names for noise info.
        [rms, ]
    psf_basenames_moffat : a list of str
        A list of base names for Moffat psf info.
        [seeing, beta, e1, e2]
    psf_basenames_airy : a list of str
        A list of base names for Airy psf info.
        [lam, diam, obscuration, e1, e2]
    id_basenames : a list of str
        A list of column names for varChips.
        [chip_id, expo_id]
    multiple_exposures_list : a list of boolean
        Use multiple noise info for multiple exposures or not, a list for MultiBand.
    N_exposures_list : a list of int
        Number of exposures, only used when multiple_exposures==True or varChips==True 
    file4varChips : str
        Path to the noise file with psf for each chip, only used for varChips==True
    varChips_list : list of boolean
        Use different psfs for different chips or not, a list for MultiBand.
    N_chips_list : a list of int
        Number of chips in each exposure, only used when varChips==True 
    psf_PixelIma_dir : str
        directory to psf images.

    Returns
    -------
    noise_info : DataFrame
        noise info
    """

    logger.info('Collect observational info from input noise catalogue...')

    # due to compatibility
    if len(psf_type_list) != len(bands):
        psf_type_list *= len(bands)

    # load file
    file_type = cata_pathfile[-3:]
    if file_type == 'csv':
        cata = pd.read_csv(cata_pathfile)
    elif file_type == 'her':
        cata = pd.read_feather(cata_pathfile)
    elif file_type == 'its':
        with fits.open(cata_pathfile) as hdul:
            cata = Table(hdul[1].data).to_pandas()
    else:
        raise Exception(f'Not supported file type! {cata_pathfile}')
    logger.info(f'noise cata: {cata_pathfile}')
    ## for varChips
    if file4varChips:
        file_type = file4varChips[-3:]
        if file_type == 'csv':
            cata_varChips = pd.read_csv(file4varChips)
        elif file_type == 'her':
            cata_varChips = pd.read_feather(file4varChips)
        elif file_type == 'its':
            with fits.open(file4varChips) as hdul:
                cata_varChips = Table(hdul[1].data).to_pandas()
        else:
            raise Exception(f'Not supported file type! {file4varChips}')
        logger.info(f'cata for varChips: {file4varChips}')

    # if only labels are required
    if only_labels:
        if noise_psf_basenames:
            name_label = noise_psf_basenames[0]
        else:
            name_label = label_basename
        return cata[name_label].to_list()

    # noise and seeing info
    ### >>> for old version
    if noise_psf_basenames:
        for psf_type in psf_type_list:
            if psf_type.lower() != 'moffat':
                raise Exception('old noise_psf_basenames does not support other psf_type!')
        name_label, name_rms, name_seeing, name_MoffatBeta, name_moffat_e1, name_moffat_e2, name_chip_id, name_expo_id = noise_psf_basenames
    ### >>> new version
    else:
        name_label = label_basename
        name_rms = noise_basenames[0]
        try:
            name_chip_id, name_expo_id = id_basenames
        except TypeError:
            name_chip_id = 'none'
            name_expo_id = 'none'
        # PSF
        if 'moffat' in psf_type_list:
            try:
                name_seeing, name_MoffatBeta, name_moffat_e1, name_moffat_e2 = psf_basenames_moffat
            except ValueError:
                name_seeing, name_MoffatBeta = psf_basenames_moffat
                name_moffat_e1 = 'none'
                name_moffat_e2 = 'none'
        if 'airy' in psf_type_list:
            try:
                name_lam, name_diam, name_obscuration, name_airy_e1, name_airy_e2 = psf_basenames_airy
            except ValueError:
                name_lam, name_diam, name_obscuration = psf_basenames_airy
                name_airy_e1 = 'none'
                name_airy_e2 = 'none'

    # collect info
    noise_info = pd.DataFrame({'label': cata[name_label]})
    for i_band, band in enumerate(bands):

        # what is the psf type ?
        psf_type = psf_type_list[i_band]
        ## PSF images does not need PSF parameters
        if psf_type == 'pixelima':
            logger.info(f'band {band} uses PSF images.')
            noise_info.loc[:, f'PixelIma_dir_{band}'] = \
                [os.path.join(psf_PixelIma_dir, band, f'tile{tile}') for tile in cata[name_label].to_list()]
            need_psf_para = False
        ## in case different profile has different psf_e name
        elif psf_type == 'moffat':
            name_psf_e1 = name_moffat_e1
            name_psf_e2 = name_moffat_e2
            need_psf_para = True
        elif psf_type == 'airy':
            name_psf_e1 = name_airy_e1
            name_psf_e2 = name_airy_e2
            need_psf_para = True

        try:
            multiple_exposures = multiple_exposures_list[i_band]
        except (IndexError, TypeError) as error:
            multiple_exposures = False
        try:
            varChips = varChips_list[i_band]
        except (IndexError, TypeError) as error:
            varChips = False

        if varChips:
            logger.info(f'band {band} uses varChips.')
            if not file4varChips:
                raise Exception('noise catalogue for varChips is not provided!')

            # number of exposures
            N_exposures = N_exposures_list[i_band]
            logger.info(f'number of exposures {N_exposures}')

            # number of chips
            N_chips = N_chips_list[i_band]
            logger.info(f'number of chips {N_chips}')

            for i_expo in range(N_exposures):
                ## noise info still from the main noise catalogue
                try:
                    noise_info.loc[:, f'rms_{band}_expo{i_expo}'] = np.array(cata[f'{name_rms}_{band}_expo{i_expo}']).astype(float)
                except KeyError:
                    ## same value for all exposures
                    noise_info.loc[:, f'rms_{band}_expo{i_expo}'] = np.array(cata[f'{name_rms}_{band}']).astype(float)
                    if i_expo == 0:
                        logger.warning('Use same rms for all exposures in varChips mode!')
                ## psf info from a separate noise catalogue
                if need_psf_para:
                    for i_chip in range(N_chips):
                        for tile_label in noise_info['label'].values:
                            # select the desired row
                            mask_chip = (cata_varChips[name_label]==tile_label) & (cata_varChips[name_expo_id]==i_expo) & (cata_varChips[name_chip_id]==i_chip)
                            if np.sum(mask_chip) != 1:
                                raise Exception(f'Number of psf in tile {tile_label}, expo {i_expo}, chip {i_chip} is wrong with value {np.sum(mask_chip)}!')
                            # always assume cata contains all desired parameters
                            if psf_type.lower() == 'moffat':
                                noise_info.loc[noise_info['label']==tile_label, f'seeing_{band}_expo{i_expo}_chip{i_chip}'] = np.array(cata_varChips.loc[mask_chip, f'{name_seeing}_{band}']).astype(float)[0]
                                noise_info.loc[noise_info['label']==tile_label, f'beta_{band}_expo{i_expo}_chip{i_chip}'] = np.array(cata_varChips.loc[mask_chip, f'{name_MoffatBeta}_{band}']).astype(float)[0]
                            elif psf_type.lower() == 'airy':
                                noise_info.loc[noise_info['label']==tile_label, f'lam_{band}_expo{i_expo}_chip{i_chip}'] = np.array(cata_varChips.loc[mask_chip, f'{name_lam}_{band}']).astype(float)[0]
                                noise_info.loc[noise_info['label']==tile_label, f'diam_{band}_expo{i_expo}_chip{i_chip}'] = np.array(cata_varChips.loc[mask_chip, f'{name_diam}_{band}']).astype(float)[0]
                                noise_info.loc[noise_info['label']==tile_label, f'obscuration_{band}_expo{i_expo}_chip{i_chip}'] = np.array(cata_varChips.loc[mask_chip, f'{name_obscuration}_{band}']).astype(float)[0]

                            # psf e 
                            try:
                                noise_info.loc[noise_info['label']==tile_label, f'psf_e1_{band}_expo{i_expo}_chip{i_chip}'] = np.array(cata_varChips.loc[mask_chip, f'{name_psf_e1}_{band}']).astype(float)[0]
                                noise_info.loc[noise_info['label']==tile_label, f'psf_e2_{band}_expo{i_expo}_chip{i_chip}'] = np.array(cata_varChips.loc[mask_chip, f'{name_psf_e2}_{band}']).astype(float)[0]
                            except KeyError:
                                noise_info.loc[noise_info['label']==tile_label, f'psf_e1_{band}_expo{i_expo}_chip{i_chip}'] = 0
                                noise_info.loc[noise_info['label']==tile_label, f'psf_e2_{band}_expo{i_expo}_chip{i_chip}'] = 0

        elif multiple_exposures:
            logger.info(f'band {band} uses diffExpo.')

            # number of exposures
            N_exposures = N_exposures_list[i_band]
            logger.info(f'number of exposures {N_exposures}')

            for i_expo in range(N_exposures):

                # noise
                try:
                    noise_info.loc[:, f'rms_{band}_expo{i_expo}'] = np.array(cata[f'{name_rms}_{band}_expo{i_expo}']).astype(float)
                except KeyError:
                    ## same value for all exposures
                    noise_info.loc[:, f'rms_{band}_expo{i_expo}'] = np.array(cata[f'{name_rms}_{band}']).astype(float)
                    if i_expo == 0:
                        logger.warning('Use same rms for all exposures in diffExpo mode!')

                # PSF
                if need_psf_para:
                    if psf_type.lower() == 'moffat':
                        try:
                            noise_info.loc[:, f'seeing_{band}_expo{i_expo}'] = np.array(cata[f'{name_seeing}_{band}_expo{i_expo}']).astype(float)
                        except KeyError:
                            ## same value for all exposures
                            noise_info.loc[:, f'seeing_{band}_expo{i_expo}'] = np.array(cata[f'{name_seeing}_{band}']).astype(float)
                            if i_expo == 0:
                                logger.warning('Use same seeing for all exposures in diffExpo mode!')
                        try:
                            noise_info.loc[:, f'beta_{band}_expo{i_expo}'] = np.array(cata[f'{name_MoffatBeta}_{band}_expo{i_expo}']).astype(float)
                        except KeyError:
                            ## same value for all exposures
                            noise_info.loc[:, f'beta_{band}_expo{i_expo}'] = np.array(cata[f'{name_MoffatBeta}_{band}']).astype(float)
                            if i_expo == 0:
                                logger.warning('Use same beta for all exposures in diffExpo mode!')
                    elif psf_type.lower() == 'airy':
                        try:
                            noise_info.loc[:, f'lam_{band}_expo{i_expo}'] = np.array(cata[f'{name_lam}_{band}_expo{i_expo}']).astype(float)
                        except KeyError:
                            ## same value for all exposures
                            noise_info.loc[:, f'lam_{band}_expo{i_expo}'] = np.array(cata[f'{name_lam}_{band}']).astype(float)
                            if i_expo == 0:
                                logger.warning('Use same lam for all exposures in diffExpo mode!')
                        try:
                            noise_info.loc[:, f'diam_{band}_expo{i_expo}'] = np.array(cata[f'{name_diam}_{band}_expo{i_expo}']).astype(float)
                        except KeyError:
                            ## same value for all exposures
                            noise_info.loc[:, f'diam_{band}_expo{i_expo}'] = np.array(cata[f'{name_diam}_{band}']).astype(float)
                            if i_expo == 0:
                                logger.warning('Use same diam for all exposures in diffExpo mode!')
                        try:
                            noise_info.loc[:, f'obscuration_{band}_expo{i_expo}'] = np.array(cata[f'{name_obscuration}_{band}_expo{i_expo}']).astype(float)
                        except KeyError:
                            ## same value for all exposures
                            noise_info.loc[:, f'obscuration_{band}_expo{i_expo}'] = np.array(cata[f'{name_obscuration}_{band}']).astype(float)
                            if i_expo == 0:
                                logger.warning('Use same obscuration for all exposures in diffExpo mode!')

                    # psf e
                    try:
                        try:
                            noise_info.loc[:, f'psf_e1_{band}_expo{i_expo}'] = np.array(cata[f'{name_psf_e1}_{band}_expo{i_expo}']).astype(float)
                        except KeyError:
                            ## same value for all exposures
                            noise_info.loc[:, f'psf_e1_{band}_expo{i_expo}'] = np.array(cata[f'{name_psf_e1}_{band}']).astype(float)
                            if i_expo == 0:
                                logger.warning('Use same psf_e1 for all exposures in diffExpo mode!')
                    except KeyError:
                        # (if not provided assuming 0)
                        noise_info.loc[:, f'psf_e1_{band}_expo{i_expo}'] = np.zeros(len(cata))
                        if i_expo == 0:
                            logger.warning('psf_e1 not provided, set to 0!')
                    try:
                        try:
                            noise_info.loc[:, f'psf_e2_{band}_expo{i_expo}'] = np.array(cata[f'{name_psf_e2}_{band}_expo{i_expo}']).astype(float)
                        except KeyError:
                            ## same value for all exposures
                            noise_info.loc[:, f'psf_e2_{band}_expo{i_expo}'] = np.array(cata[f'{name_psf_e2}_{band}']).astype(float)
                            if i_expo == 0:
                                logger.warning('Use same psf_e2 for all exposures in diffExpo mode!')
                    except KeyError:
                        # (if not provided assuming 0)
                        noise_info.loc[:, f'psf_e2_{band}_expo{i_expo}'] = np.zeros(len(cata))
                        if i_expo == 0:
                            logger.warning('psf_e2 not provided, set to 0!')
        else:
            logger.info(f'band {band} uses simple noise.')

            # noise
            noise_info.loc[:, f'rms_{band}'] = np.array(cata[f'{name_rms}_{band}']).astype(float)

            # psf
            if need_psf_para:
                if psf_type.lower() == 'moffat':
                    noise_info.loc[:, f'seeing_{band}'] = np.array(cata[f'{name_seeing}_{band}']).astype(float)
                    noise_info.loc[:, f'beta_{band}'] = np.array(cata[f'{name_MoffatBeta}_{band}']).astype(float)

                elif psf_type.lower() == 'airy':
                    noise_info.loc[:, f'lam_{band}'] = np.array(cata[f'{name_lam}_{band}']).astype(float)
                    noise_info.loc[:, f'diam_{band}'] = np.array(cata[f'{name_diam}_{band}']).astype(float)
                    noise_info.loc[:, f'obscuration_{band}'] = np.array(cata[f'{name_obscuration}_{band}']).astype(float)

                # psf e
                try:
                    noise_info.loc[:, f'psf_e1_{band}'] = np.array(cata[f'{name_psf_e1}_{band}']).astype(float)
                except KeyError:
                    # (if not provided assuming 0)
                    noise_info.loc[:, f'psf_e1_{band}'] = np.zeros(len(cata))
                    logger.warning('psf_e1 not provided, set to 0!')
                try:
                    noise_info.loc[:, f'psf_e2_{band}'] = np.array(cata[f'{name_psf_e2}_{band}']).astype(float)
                except KeyError:
                    # (if not provided assuming 0)
                    noise_info.loc[:, f'psf_e2_{band}'] = np.zeros(len(cata))
                    logger.warning('psf_e2 not provided, set to 0!')

    logger.debug('Desired noise info collected to DataFrame.')

    return noise_info
