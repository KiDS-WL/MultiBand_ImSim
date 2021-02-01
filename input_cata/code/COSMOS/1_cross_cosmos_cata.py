# -*- coding: utf-8 -*-
# @Author: lshuns
# @Date:   2021-01-27 10:48:13
# @Last Modified by:   lshuns
# @Last Modified time: 2021-01-30 18:15:32

### cross-match different cosmos catalogues:
###     Morphology: ACS-GC https://cdsarc.unistra.fr/ftp/J/ApJS/200/9/fits/ACS-GC_readme_2012.txt
###     Photometry: COSMOS15 DR2.1 https://cosmohub.pic.es/catalogs/160
###         u, B, V, r, ip, zp, zpp, Y, J, H, Ks
###     Redshifts: Jan Luca's compilation + pau + COSMOS15

import numpy as np
import pandas as pd
from astropy.io import fits

import os
import sys
# Self-defined package
sys.path.insert(0, os.path.realpath('../../../modules')) 
from CrossMatch import KDTreeFunc

# ++++++++++++ Files 
cata_dir = '/disks/shear10/ssli/ImSim/input/COSMOS_cata/'

## inputs
shape_file = cata_dir + 'outside/cosmos_i_public_catalog_V1.0.fits.gz'
cosmos15_file = cata_dir + 'outside/cosmos_2015_dr_2_1_COSMOS_flagged.fits'
pau_file = cata_dir + 'outside/pau_cosmos_photoz_v0_4.csv'
Jan_file = cata_dir + 'outside/COSMOS_specz_PAUS_COSMOS2015.fits'

## output
outpath = cata_dir + 'cosmos_shape_z_uBVriZYJHKs.feather'

# ++++++++++++ shape cata
with fits.open(shape_file) as hdul:
    shape_cata = hdul[1].data
print('ACS-GC number', len(shape_cata))
### desired columns
shape_cata = pd.DataFrame({
    'OBJNO': shape_cata['OBJNO'].astype(int),
    'RA': shape_cata['RA'].astype(float),
    'DEC': shape_cata['DEC'].astype(float),
    'SPECZ': shape_cata['SPECZ'].astype(float),
    'Z_gc':  shape_cata['Z'].astype(float),
    'ZQUALITY': shape_cata['ZQUALITY'].astype(float),
    'FLAG_GALFIT_HI': shape_cata['FLAG_GALFIT_HI'].astype(int), 'CHI2NU_HI': shape_cata['CHI2NU_HI'].astype(float),
    'MAG_GALFIT_HI': shape_cata['MAG_GALFIT_HI'].astype(float), 'MAGERR_GALFIT_HI': shape_cata['MAGERR_GALFIT_HI'].astype(float),
    'RE_GALFIT_HI': shape_cata['RE_GALFIT_HI'].astype(float), 'N_GALFIT_HI': shape_cata['N_GALFIT_HI'].astype(float), 'BA_GALFIT_HI': shape_cata['BA_GALFIT_HI'].astype(float),
    'REERR_GALFIT_HI': shape_cata['REERR_GALFIT_HI'].astype(float), 'NERR_GALFIT_HI': shape_cata['NERR_GALFIT_HI'].astype(float), 'BAERR_GALFIT_HI': shape_cata['BAERR_GALFIT_HI'].astype(float)
    })

# ++++++++++++ COSMOS15 & PAU
## COSMOS15
with fits.open(cosmos15_file) as hdul:
    cosmo15_cata = hdul[1].data
print('COSMOS15 total number', len(cosmo15_cata))
### select galaxies
mask_gal = (cosmo15_cata['type']==0)
cosmo15_cata = cosmo15_cata[mask_gal]
print('COSMOS15 galaxies number', len(cosmo15_cata))
### desired columns
cosmo15_cata = pd.DataFrame({
    'number': cosmo15_cata['number'].astype(int),
    'alpha_j2000': cosmo15_cata['alpha_j2000'].astype(float),
    'delta_j2000': cosmo15_cata['delta_j2000'].astype(float),
    'u_mag_auto': cosmo15_cata['u_mag_auto'].astype(float),
    'u_magerr_auto': cosmo15_cata['u_magerr_auto'].astype(float),
    'b_mag_auto': cosmo15_cata['b_mag_auto'].astype(float),
    'b_magerr_auto': cosmo15_cata['b_magerr_auto'].astype(float),
    'v_mag_auto': cosmo15_cata['v_mag_auto'].astype(float),
    'v_magerr_auto': cosmo15_cata['v_magerr_auto'].astype(float),
    'r_mag_auto': cosmo15_cata['r_mag_auto'].astype(float),
    'r_magerr_auto': cosmo15_cata['r_magerr_auto'].astype(float),
    'ip_mag_auto': cosmo15_cata['ip_mag_auto'].astype(float),
    'ip_magerr_auto': cosmo15_cata['ip_magerr_auto'].astype(float),
    'zp_mag_auto': cosmo15_cata['zp_mag_auto'].astype(float),
    'zp_magerr_auto': cosmo15_cata['zp_magerr_auto'].astype(float),
    'zpp_mag_auto': cosmo15_cata['zpp_mag_auto'].astype(float),
    'zpp_magerr_auto': cosmo15_cata['zpp_magerr_auto'].astype(float),
    'y_mag_auto': cosmo15_cata['y_mag_auto'].astype(float),
    'y_magerr_auto': cosmo15_cata['y_magerr_auto'].astype(float),
    'j_mag_auto': cosmo15_cata['j_mag_auto'].astype(float),
    'j_magerr_auto': cosmo15_cata['j_magerr_auto'].astype(float),
    'h_mag_auto': cosmo15_cata['h_mag_auto'].astype(float),
    'h_magerr_auto': cosmo15_cata['h_magerr_auto'].astype(float),
    'ks_mag_auto': cosmo15_cata['ks_mag_auto'].astype(float),
    'ks_magerr_auto': cosmo15_cata['ks_magerr_auto'].astype(float),
    'photoz': cosmo15_cata['photoz'].astype(float)
    })

## PAU
z_pau_cata = pd.read_csv(pau_file, comment='#', na_values=r'\N')
print('PAU number', len(z_pau_cata))
### desired columns
z_pau_cata = pd.DataFrame({
    'id_laigle': z_pau_cata['id_laigle'].astype(int),
    'zspec_mean': z_pau_cata['zspec_mean'].astype(float),
    'zspec_std': z_pau_cata['zspec_std'].astype(float),
    'photoz': z_pau_cata['photoz'].astype(float)
    })

## join pau photo-z to cosmo15
z_photo_cata = cosmo15_cata.merge(z_pau_cata, how='left', left_on='number', right_on='id_laigle', suffixes=[None, '_pau'])
z_photo_cata.drop(columns=['id_laigle'], inplace=True)
z_photo_cata.reset_index(drop=True, inplace=True)
print('combined z_photo_cata', len(z_photo_cata))

# ++++++++++++ Jan Luca's redshift cata
with fits.open(Jan_file) as hdul:
    z_cata_Jan = hdul[1].data
print('Jan cata number', len(z_cata_Jan))
## desired columns
z_cata_Jan = pd.DataFrame({
    'ID_Jan': np.arange(len(z_cata_Jan)),
    'RA_Jan': z_cata_Jan['RA'].astype(float),
    'DEC_Jan': z_cata_Jan['DEC'].astype(float),
    'Zbest_Jan': z_cata_Jan['Zbest'].astype(float)
    })

# ++++++++++++ cross match
# shape_cata
xy_shape = np.empty((len(shape_cata), 2), dtype=float)
xy_shape[:, 0] = np.array(shape_cata['RA'], dtype=float)
xy_shape[:, 1] = np.array(shape_cata['DEC'], dtype=float)

# 1. z_photo_cata
xy_z = np.empty((len(z_photo_cata), 2), dtype=float)
xy_z[:, 0] = np.array(z_photo_cata['alpha_j2000'], dtype=float)
xy_z[:, 1] = np.array(z_photo_cata['delta_j2000'], dtype=float)
## KDTree match
dist, ind = KDTreeFunc(xy_shape, xy_z, max_distance=0.25/3600., unique=True, k=1, leafsize=100)
flag_matched = ind<len(z_photo_cata)
flag_miss = np.invert(flag_matched)
print('Number of matched (COSMOS15 x ACS-GC)', np.sum(flag_matched))
print('     matched/COSMOS15', np.sum(flag_matched)/len(z_photo_cata))
print('     matched/ACS-GC', np.sum(flag_matched)/len(shape_cata))
ind_matched = ind[flag_matched]
## save id
shape_cata.loc[flag_matched, 'number'] = np.array(z_photo_cata.loc[ind_matched, 'number'])

# 2. z_cata_Jan
xy_z = np.empty((len(z_cata_Jan), 2), dtype=float)
xy_z[:, 0] = np.array(z_cata_Jan['RA_Jan'], dtype=float)
xy_z[:, 1] = np.array(z_cata_Jan['DEC_Jan'], dtype=float)
## KDTree match
dist, ind = KDTreeFunc(xy_shape, xy_z, max_distance=0.25/3600., unique=True, k=4, leafsize=100)
flag_matched = ind<len(z_cata_Jan)
flag_miss = np.invert(flag_matched)
print('Number of matched (z_cata_Jan x ACS-GC)', np.sum(flag_matched))
print('     matched/z_cata_Jan', np.sum(flag_matched)/len(z_cata_Jan))
print('     matched/ACS-GC', np.sum(flag_matched)/len(shape_cata))
ind_matched = ind[flag_matched]
## save id
shape_cata.loc[flag_matched, 'ID_Jan'] = np.array(z_cata_Jan.loc[ind_matched, 'ID_Jan'])

# combine
combine_cata = shape_cata.merge(z_photo_cata, on='number', how='left')
combine_cata = combine_cata.merge(z_cata_Jan, on='ID_Jan', how='left')
## rename z index
combine_cata.rename(columns={'number': 'id_laigle'}, inplace=True)
## dummy consistent with ACS-GC catalogue
combine_cata.fillna(-999, inplace=True)

# ++++++++++++ redshift selection
## pre-selection
##### large std in zspec_mean
N_removed = np.sum(combine_cata['zspec_std']>combine_cata['zspec_mean']/3.)
N_tot = np.sum(combine_cata['zspec_mean']>-999)
combine_cata.loc[combine_cata['zspec_std']>combine_cata['zspec_mean']/3., 'zspec_mean'] = -999
print(f'removed PAU_spec number (percent): {N_removed} ({N_removed/N_tot})')
##### bad quality in SPECZ
good_quality = ((combine_cata['ZQUALITY']>3) & (combine_cata['ZQUALITY']<5)) | (combine_cata['ZQUALITY']==1.5) | (combine_cata['ZQUALITY']==2.4) | (combine_cata['ZQUALITY']==2.5) | (combine_cata['ZQUALITY']==9.3) | (combine_cata['ZQUALITY']==9.5) | ((combine_cata['ZQUALITY']>13) & (combine_cata['ZQUALITY']<15)) | ((combine_cata['ZQUALITY']>23) & (combine_cata['ZQUALITY']<25))
bad_quality = np.invert(good_quality)
N_removed = np.sum((bad_quality)&(combine_cata['SPECZ']>-999))
N_tot = np.sum(combine_cata['SPECZ']>-999)
combine_cata.loc[bad_quality, 'SPECZ'] = -999
print(f'removed zCOSMOS number (percent): {N_removed} ({N_removed/N_tot})')
##### not good photoz_pau
N_removed = np.sum((combine_cata['ip_mag_auto']>23)&(combine_cata['photoz_pau']>-999))
N_tot = np.sum(combine_cata['photoz_pau']>-999)
combine_cata.loc[combine_cata['ip_mag_auto']>23., 'photoz_pau'] = -999
print(f'removed photoz_pau number (percent): {N_removed} ({N_removed/N_tot})')

## priority: Zbest_Jan > zspec_mean (zspec from PAU) > SPECZ (zspec from zCOSMOS) > photoz_pau (zphoto from PAU) > photoz (zphoto from COSMOS15) > Z_gc (zphoto from Ilbert 2009)
z_candidate = ['Zbest_Jan', 'zspec_mean', 'SPECZ', 'photoz_pau', 'photoz', 'Z_gc']
z_ori = ['JanLuca_Zbest', 'PAU_spec', 'zCOSMOS', 'PAU_photo', 'COSMOS15', 'Ilbert2009']
per_tot = 1
for i_z in range(len(z_candidate)):
    z_selec = z_candidate[i_z]
    z_name = z_ori[i_z]
    if i_z==0:
        combine_cata.loc[:, 'Z_optimal'] = np.array(combine_cata[z_selec])
        combine_cata.loc[(combine_cata[z_selec]>-999), 'Z_ORIGIN'] = z_name
    else:
        combine_cata.loc[mask_remain, 'Z_optimal'] = np.array(combine_cata.loc[mask_remain, z_selec])
        combine_cata.loc[mask_remain & (combine_cata[z_selec]>-999), 'Z_ORIGIN'] = z_name

    mask_remain = (combine_cata['Z_optimal']==-999)
    print(z_name, 'percent filled', per_tot-np.sum(mask_remain)/len(combine_cata))
    per_tot = np.sum(mask_remain)/len(combine_cata)
print('Total percent filled', 1-per_tot)
print('columns', combine_cata.columns)

# ++++++++++++ save
combine_cata.to_feather(outpath)
print('Saved to', outpath)

######### info
# ACS-GC number 304688
# COSMOS15 total number 536077
# COSMOS15 galaxies number 518404
# PAU number 40672
# combined z_photo_cata 518404
# Jan cata number 527253
# Number of matched (COSMOS15 x ACS-GC) 205045
#      matched/COSMOS15 0.39553128448082964
#      matched/ACS-GC 0.6729671007719372
# Number of matched (z_cata_Jan x ACS-GC) 208459
#      matched/z_cata_Jan 0.39536806808116787
#      matched/ACS-GC 0.6841720054613244
# removed PAU_spec number (percent): 9 (0.0007685082401161301)
# removed zCOSMOS number (percent): 1557 (0.1554978527913712)
# removed photoz_pau number (percent): 1387 (0.03541156045751634)
# JanLuca_Zbest percent filled 0.6841720054613244
# PAU_spec percent filled 0.0015294333875964794
# zCOSMOS percent filled 0.0028914824344903334
# PAU_photo percent filled 0.002802867195294878
# COSMOS15 percent filled 0.0023729191829018337
# Ilbert2009 percent filled 0.23798114792837266
# Total percent filled 0.9317498555899806
# columns Index(['OBJNO', 'RA', 'DEC', 'SPECZ', 'Z_gc', 'ZQUALITY', 'FLAG_GALFIT_HI',
#        'CHI2NU_HI', 'MAG_GALFIT_HI', 'MAGERR_GALFIT_HI', 'RE_GALFIT_HI',
#        'N_GALFIT_HI', 'BA_GALFIT_HI', 'REERR_GALFIT_HI', 'NERR_GALFIT_HI',
#        'BAERR_GALFIT_HI', 'id_laigle', 'ID_Jan', 'alpha_j2000', 'delta_j2000',
#        'u_mag_auto', 'u_magerr_auto', 'b_mag_auto', 'b_magerr_auto',
#        'v_mag_auto', 'v_magerr_auto', 'r_mag_auto', 'r_magerr_auto',
#        'ip_mag_auto', 'ip_magerr_auto', 'zp_mag_auto', 'zp_magerr_auto',
#        'zpp_mag_auto', 'zpp_magerr_auto', 'y_mag_auto', 'y_magerr_auto',
#        'j_mag_auto', 'j_magerr_auto', 'h_mag_auto', 'h_magerr_auto',
#        'ks_mag_auto', 'ks_magerr_auto', 'photoz', 'zspec_mean', 'zspec_std',
#        'photoz_pau', 'RA_Jan', 'DEC_Jan', 'Zbest_Jan', 'Z_optimal',
#        'Z_ORIGIN'],
#       dtype='object')
