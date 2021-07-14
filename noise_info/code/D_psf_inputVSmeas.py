# @Author: ssli
# @Date:   2021-05-11, 14:50:53
# @Last modified by:   ssli
# @Last modified time: 2021-05-11, 15:16:35

import os
import glob
import time
import galsim
import subprocess

import numpy as np
import pandas as pd
import scipy.spatial as sst
import multiprocessing as mp

from scipy import stats
from astropy.io import fits

# +++++++++++++++++ I/O

# original noise files
noise_files = ['../test_noise_ugri.csv', '../test_noise_NIR.csv']
## corresponding bands
bands = [['u', 'g', 'r', 'i'], ['Z', 'Y', 'J', 'H', 'Ks']]
## pixel scales
pixel_scale_list = [[0.214, 0.214, 0.214, 0.214], [0.34, 0.34, 0.34, 0.34, 0.34]] # arcsec/pixel
pixel_scale_final_list = [[0.2, 0.2, 0.2, 0.2], [0.34, 0.34, 0.34, 0.34, 0.34]] # arcsec/pixel

# file to save difference info
dpsf_files = ['../test_psf_diff_inVSmeas_ugri.csv', '../test_psf_diff_inVSmeas_NIR.csv']

# code for psf measurement
psf2moffat_code = '/disks/shear10/ssli/KiDS_noise_seeing/code/psf2moffat'

# directory for running saving
outdir = '/disks/shear10/ssli/KiDS_noise_seeing/PSF_inputVSmeas/'
### where to save images
outdir_ima = os.path.join(outdir, 'psf_map')
if not os.path.isdir(outdir_ima):
    os.mkdir(outdir_ima)
### where to save catalogues
outdir_cata = os.path.join(outdir, 'psf_map_cata')
if not os.path.isdir(outdir_cata):
    os.mkdir(outdir_cata)

# bands need to be swarped
bands_swarp = ['u', 'g', 'r', 'i']

# ++++++++++++++++ functions
def MoffatPSF(seeing, moffat_beta=3.5, psf_e=[]):
    """
    Generate PSF model.

    Parameters
    ----------
    seeing : float
        Full-width-half-max of the PSF.
    moffat_beta : float
        Beta parameter for Moffat profile.
    psf_e : list, optional (default: [0, 0])
        List of psf ellipticity [e1, e2] (e = 1 - q).

    Returns
    -------
    psf: galsim Moffat
        PSF model.
    """

    psf = galsim.Moffat(beta=moffat_beta, fwhm=seeing, trunc=4.5*seeing)
    if psf_e != []:
        psf_e1 = psf_e[0]
        psf_e2 = psf_e[1]

        psf_e = np.sqrt(psf_e1**2+psf_e2**2)
        # g_i = e_i/(2-e)
        psf_g1, psf_g2 = psf_e1/(2-psf_e), psf_e2/(2-psf_e)

        psf = psf.shear(g1=psf_g1, g2=psf_g2)

    return psf

def PSFmap(PSF, pixel_scale, mag_input, mag_zero=30., N_PSF=100, sep_PSF=120, rng_seed=940120):
    """
    Generate PSF stars for the image
        with uniform magnitude

    Parameters
    ----------
    PSF : galsim object:
        PSF model.
    pixel_scale (arcsec) : float
        Pixel size in unit of arcsec.
    mag_input : float
        The input magnitude for PSF.
    mag_zero : float, optional (default: 30.)
        The zero point for magnitude.
    N_PSF : int, optional (default: 100)
        Number of PSF generated, better be the square of some integer.
    sep_PSF (pixel) : int, optional (default: 120)
        Seperation between PSFs.
    rng_seed : int, optional (default: 940120)
        Seed for random number generator.

    Returns
    -------
    psf_image: galsim Images
        Image for PSFs.
    """

    # magnitude to flux
    flux = 10**(-0.4*(mag_input-mag_zero))

    # position
    # separation in pixels
    separation = int(sep_PSF)
    ## number in each row
    Nrow = int(N_PSF**0.5)
    ## get grid center
    x = np.arange(separation, separation+Nrow*separation, separation, dtype='int')
    y = np.repeat(x, Nrow)
    x = np.tile(x, Nrow)
    ## check outliers
    Nrow = N_PSF - len(x)
    if Nrow > 0:
        x = np.concatenate([x, np.arange(separation, separation+Nrow*separation, separation)])
        y = np.concatenate([y, np.full(Nrow, y[-1]+separation)])
    elif Nrow < 0:
        x = x[:N_PSF]
        y = y[:N_PSF]
    ## make random shift
    np.random.seed(rng_seed)
    shift_lim = int(separation/5.)
    dx_dy = np.random.randint(low=-shift_lim, high=shift_lim, size=N_PSF)
    x += dx_dy
    y += dx_dy

    # initiate a canvas
    psf_image = galsim.ImageF(int(np.amax(x)+separation), int(np.max(y)+separation), scale=pixel_scale)
    # print(f"PSF image bounds {psf_image.bounds}")

    # draw PSF
    star_like_PSF = PSF.withFlux(flux)
    ## draw stamp
    stamp_PSF = star_like_PSF.drawImage(scale=pixel_scale)

    # place PSF
    for i in range(N_PSF):

        # position the stamp
        stamp_PSF.setCenter(x[i], y[i])

        # stitch the postage stamp to the canvas
        overlap = stamp_PSF.bounds & psf_image.bounds
        psf_image[overlap] += stamp_PSF[overlap]

    return psf_image

def GaussianNoise(noise_rms, rng_seed=940120):
    """
    Generate Gaussian background noise.

    Parameters
    ----------
    noise_rms : float
        Observational noise rms for the whole image.
    rng_seed : int, optional (default: 940120)
        Seed for random number generator.

    Returns
    -------
    noise: galsim object
        Background noise.
    """

    rng = galsim.BaseDeviate(rng_seed)
    noise = galsim.GaussianNoise(rng, sigma=noise_rms)

    return noise

def SwarpImage(image_in, swarp_config_file,
                    image_out,
                    only_resample=False,
                    running_log=True, log_dir=None,
                    swarp_path='swarp',
                    clean_up_level=0):
    """
    SWarp for coadding or resampling.
        only_resample: set to True
                        if only resampling but no coadding.
        Clean up level:
            0: none
            1: rm original images
    """

    # first check if already exist
    try:
        with fits.open(image_out) as hdul:
            head_tmp = hdul[0].header
        flag_sim = head_tmp['flag_sim']
        if flag_sim >= 2:
            print(f"{image_out} already exist.")
            return 1
    except FileNotFoundError:
        pass
    except (KeyError, OSError) as e:
        os.remove(image_out)
        print("Remove existing images.")

    # running info
    if running_log:
        basename = os.path.basename(image_out)
        outLog = open(os.path.join(log_dir, basename.replace('.fits', '.log')), "w")
        errLog = open(os.path.join(log_dir, basename.replace('.fits', '.err.log')), "w")
    else:
        outLog = subprocess.PIPE
        errLog = subprocess.STDOUT

    # resampling out directory
    RESAMPLE_DIR = os.path.dirname(image_out)

    # build command
    cmd = [swarp_path]
    if isinstance(image_in, str):
        cmd.append(image_in)
    else:
        cmd.extend(image_in) # in case image_in containing multiple inputs
    if only_resample:
        cmd.extend(['-c', swarp_config_file,
            '-RESAMPLE_DIR', RESAMPLE_DIR])
        # print('Running SWarp for only resampling...')
    else:
        cmd.extend(['-c', swarp_config_file,
            '-RESAMPLE_DIR', RESAMPLE_DIR,
            '-IMAGEOUT_NAME', image_out, '-WEIGHTOUT_NAME', image_out.replace('.fits', '.weight.fits')])
        # print('Running SWarp for coadding...')

    # run
    proc = subprocess.run(cmd, stdout=outLog, stderr=errLog)

    if running_log:
        outLog.close()
        errLog.close()

    # rename for only resampling
    if only_resample:
        basename = os.path.basename(image_in)
        os.rename(os.path.join(RESAMPLE_DIR, basename.replace('.fits', '.resamp.fits')), image_out)
        # print(f"Swarp resampled image saved as {image_out}")
    # else:
    #     print(f"Swarp coadded image saved as {image_out}")

    # mark success to the header
    with fits.open(image_out, mode='update') as hdul:
        head_tmp = hdul[0].header
        ## update info
        head_tmp['flag_sim'] = 2

    if clean_up_level:
        try:
            os.remove(image_in)
        except TypeError:
            for image_tmp in image_in:
                os.remove(image_tmp)
        # print('Original images are removed.')

    # print('SWarp finished.')
    return 0

# ++++++++++++++++ workhorse
start_time = time.time()
for i_file, noise_file in enumerate(noise_files):
    print('file', noise_file)
    noise_info = pd.read_csv(noise_file)

    dpsf_file = dpsf_files[i_file]

    bands_selec = bands[i_file]
    pixel_scale_list_selec = pixel_scale_list[i_file]
    pixel_scale_final_list_selec = pixel_scale_final_list[i_file]

    # for saving difference
    dpsf_info = pd.DataFrame({'label': np.array(noise_info['label'])})

    # loop over all tiles
    for i_tile, noise_row in noise_info.iterrows():

        tile = noise_row['label']
        print('tile', tile)

        # all bands
        for i_band, band in enumerate(bands_selec):

            print('band', band)

            pixel_scale = pixel_scale_list_selec[i_band]
            pixel_scale_final = pixel_scale_final_list_selec[i_band]

            # original noise info
            try:
                rms0 = noise_row[f'rmsTHELI_{band}']
            except KeyError:
                rms0 = noise_row[f'rms_{band}']
            seeing0 = noise_row[f'MeasSeeing_{band}']
            beta0 = noise_row[f'MeasBeta_{band}']
            print('seeing0', seeing0)
            print('beta0', beta0)

            # ++++++++++++++++ generate images
            rms = rms0
            seeing = seeing0
            beta = beta0
            # hard cut on the max beta value
            if beta > 4:
                beta = 4

            # magnitude for simulated PSF stars
            mag_psf = 20. - 2.5*np.log10(rms/2.)

            # random seed
            rng_seed = 940120 + i_band * 94 + i_tile * 120

            # simulation
            PSF = MoffatPSF(seeing, moffat_beta=beta, psf_e=[])
            psf_map = PSFmap(PSF, pixel_scale, mag_psf, mag_zero=30., N_PSF=50, sep_PSF=120, rng_seed=rng_seed)

            # noise background
            noise_tmp = GaussianNoise(rms, rng_seed=rng_seed)
            psf_map.addNoise(noise_tmp)

            # save
            outpath_PSF_name = os.path.join(outdir_ima, f'psf_tile{tile}_band{band}.fits')
            psf_map.write(outpath_PSF_name)

            # ++++++++++++++++ swarp
            if band in bands_swarp:
                outfile = outpath_PSF_name.replace('.fits', '.resamp.fits')
                SwarpImage(outpath_PSF_name, swarp_config_file='../../config/resampling_aw.swarp',
                                    image_out=outfile,
                                    only_resample=True,
                                    running_log=False, log_dir='./',
                                    swarp_path='swarp',
                                    clean_up_level=0)

                # ++++++++++++++++ clean original images and weight images for swarped results
                # original images
                tmp = outfile.replace('.resamp.fits', '.fits')
                os.remove(tmp)
                # weight images
                tmp = outfile.replace('.resamp.fits', '.resamp.weight.fits')
                os.remove(tmp)

                # ++++++++++++++++ rename for easy use
                os.rename(outfile, outpath_PSF_name)

            # ++++++++++++++++ measure moffat
            # catalogue name
            cata_file = os.path.join(outdir_cata, os.path.basename(outpath_PSF_name) + '.psf.cat')
            moffat_file = cata_file + '.moffat'

            # necessary sex config files
            sex_config = []
            code_dir = os.path.dirname(psf2moffat_code)
            for config_file in ['default.sex', 'default.param', 'default.conv']:
                sex_config.append(os.path.join(code_dir, f'{config_file}'))
            # run sex
            cmd = ['sex', outpath_PSF_name, '-c', sex_config[0], '-PARAMETERS_NAME', sex_config[1], '-FILTER_NAME', sex_config[2], '-DETECT_THRESH', str(10), '-CATALOG_NAME', cata_file]
            proc = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)

            # create link
            inimage = os.path.join(outdir_cata, 'inimage.fits')
            os.symlink(outpath_PSF_name, inimage)
            proc = subprocess.run(f'(echo -1; cat {cata_file}) | {psf2moffat_code} > {moffat_file}',
                    stdout=subprocess.PIPE, stderr=subprocess.STDOUT, shell=True, cwd=outdir_cata)

            # ++++++++++++++++ clean images for next run
            os.remove(inimage)
            os.remove(outpath_PSF_name)

            # ++++++++++++++++ calculate the input and measured difference
            data = np.loadtxt(moffat_file, comments='%')
            seeing1 = np.median(data[:, 10]*pixel_scale_final)
            beta1 = np.median(data[:, 11])
            print('(seeing1-seeing0)/seeing0', (seeing1-seeing0)/seeing0)
            print('(beta1-beta0)/beta0', (beta1-beta0)/beta0)
            dpsf_info.loc[noise_info['label']==tile, f'dseeing_ori_{band}'] = (seeing1-seeing0)/seeing0
            dpsf_info.loc[noise_info['label']==tile, f'dbeta_ori_{band}'] = (beta1-beta0)/beta0

            # ++++++++++++++++++ clean catalogue for next run
            os.remove(moffat_file)

            # ++++++++++++++++ update input seeing & beta, run again
            rms = rms0
            seeing = seeing0 - (seeing1-seeing0)
            beta = 1./(1./beta0 - (1./beta1-1./beta0))
            if beta > 4:
                beta = 4

            mag_psf = 20. - 2.5*np.log10(rms/2.)

            rng_seed = 940120 + i_band * 94 + i_tile * 120

            PSF = MoffatPSF(seeing, moffat_beta=beta, psf_e=[])
            psf_map = PSFmap(PSF, pixel_scale, mag_psf, mag_zero=30., N_PSF=50, sep_PSF=120, rng_seed=rng_seed)

            noise_tmp = GaussianNoise(rms, rng_seed=rng_seed)
            psf_map.addNoise(noise_tmp)

            outpath_PSF_name = os.path.join(outdir_ima, f'psf_tile{tile}_band{band}.fits')
            psf_map.write(outpath_PSF_name)

            # ++++++++++++++++ swarp
            if band in bands_swarp:
                outfile = outpath_PSF_name.replace('.fits', '.resamp.fits')
                SwarpImage(outpath_PSF_name, swarp_config_file='../../config/resampling_aw.swarp',
                                    image_out=outfile,
                                    only_resample=True,
                                    running_log=False, log_dir='./',
                                    swarp_path='swarp',
                                    clean_up_level=0)

                # ++++++++++++++++ clean original images and weight images for swarped results
                # original images
                tmp = outfile.replace('.resamp.fits', '.fits')
                os.remove(tmp)
                # weight images
                tmp = outfile.replace('.resamp.fits', '.resamp.weight.fits')
                os.remove(tmp)

                # ++++++++++++++++ rename for easy use
                os.rename(outfile, outpath_PSF_name)

            # ++++++++++++++++ measure moffat
            # catalogue name
            cata_file = os.path.join(outdir_cata, os.path.basename(outpath_PSF_name) + '.psf.cat')
            moffat_file = cata_file + '.moffat'

            # necessary sex config files
            sex_config = []
            code_dir = os.path.dirname(psf2moffat_code)
            for config_file in ['default.sex', 'default.param', 'default.conv']:
                sex_config.append(os.path.join(code_dir, f'{config_file}'))
            # run sex
            cmd = ['sex', outpath_PSF_name, '-c', sex_config[0], '-PARAMETERS_NAME', sex_config[1], '-FILTER_NAME', sex_config[2], '-DETECT_THRESH', str(10), '-CATALOG_NAME', cata_file]
            proc = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)

            # create link
            inimage = os.path.join(outdir_cata, 'inimage.fits')
            os.symlink(outpath_PSF_name, inimage)
            proc = subprocess.run(f'(echo -1; cat {cata_file}) | {psf2moffat_code} > {moffat_file}',
                    stdout=subprocess.PIPE, stderr=subprocess.STDOUT, shell=True, cwd=outdir_cata)

            # ++++++++++++++++ clean images
            os.remove(inimage)
            os.remove(outpath_PSF_name)

            # ++++++++++++++++ calculate the seeing difference
            # load output
            # print('moffat_file', moffat_file)
            data = np.loadtxt(moffat_file, comments='%')
            seeing2 = np.median(data[:, 10]*pixel_scale_final)
            beta2 = np.median(data[:, 11])
            print('(seeing2-seeing0)/seeing0', (seeing2-seeing0)/seeing0)
            print('(beta2-beta0)/beta0', (beta2-beta0)/beta0)
            dpsf_info.loc[noise_info['label']==tile, f'dseeing_corr_{band}'] = (seeing2-seeing0)/seeing0
            dpsf_info.loc[noise_info['label']==tile, f'dbeta_corr_{band}'] = (beta2-beta0)/beta0

            # ++++++++++++++++ save
            noise_info.loc[noise_info['label']==tile, f'InputSeeing_{band}'] = seeing
            noise_info.loc[noise_info['label']==tile, f'InputBeta_{band}'] = beta
            print('InputSeeing', seeing)
            print('InputBeta', beta)

    # save
    noise_info.to_csv(noise_file, index=False)
    print(f'noise info saved as {noise_file}')
    dpsf_info.to_csv(dpsf_file, index=False)
    print(f'psf difference saved as {dpsf_file}')
