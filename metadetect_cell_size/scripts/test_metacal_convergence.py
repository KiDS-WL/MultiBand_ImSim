# -*- coding: utf-8 -*-

### Does metacal give the same image whatever the size of the cell around it?
###
###    Cuts squares of 250, 500 and 1000 pixels around the same point of a
###    simulated image, runs ngmix metacal on each exactly as metadetect does
###    (psf='fitgauss', step=0.01, lanczos15 InterpolatedImage of the whole
###    square), and compares the 250 and 500 results with the 1000 one over the
###    same central 250x250 pixels. The 1000-pixel square is the reference: with
###    that much sky around them, its central pixels do not depend on where the
###    square ends. If the result does not change when the square gets bigger, it
###    has converged.
###
###    The metacal image is built in stages, each tested on its own:
###       reconv    the square convolved with the target PSF only (no deconvolution)
###       noshear   deconvolved by the PSF and reconvolved by the target PSF
###       1p, 2p    deconvolved, sheared by 0.01 and reconvolved
###
###    and for PSF stamps of several sizes. 16 and 32 are crops of the 48-pixel
###    stamp of the simulation, 96 is drawn with ImSim's own code
###    (ImSimPSF.MoffatPSF and PSFima). The Moffat is truncated at 4.5 x FWHM,
###    16.4 pixels, so 48 and 96 hold the same PSF and only the stamp size differs,
###    32 loses a sliver at its edge and 16 loses the wings.
###
###    Differences are the RMS over pixels, in units of the image noise, averaged
###    over several points of the image, in bins of max(|dx|, |dy|) from the centre.
###
###    Usage: python test_metacal_convergence.py

import os
import sys

import numpy as np
import galsim
import ngmix
from astropy.io import fits
from ngmix.metacal.metacal import MetacalFitGaussPSF

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', 'modules'))
import ImSimPSF

## ++++++++++++++ I/O and general setups

image_dir = '/sdf/data/kipac/u/liss/ImSim/output/test_dev/jamar_seeing_073/p020p020/images/original'
tile = 'tile0.0_0.0_bandLSST_r'
PIXEL_SCALE = 0.2
## the PSF of the simulation (noise_info: InputSeeing_LSST_r, InputBeta_LSST_r)
SEEING, MOFFAT_BETA = 0.73, 2.22406805360007

STEP = 0.01
CELL_SIZES = [250, 500]
REF_SIZE = 1000
PSF_SIZES = [16, 32, 48, 96]
## points of the image at which the squares are centred (cell centres for central size 100)
CENTRES = [(c * 100 + 50, d * 100 + 50) for c, d in [(20, 20), (20, 60), (40, 40), (60, 20), (60, 60), (30, 50)]]
R_EDGES = [0, 12.5, 25, 50, 75, 100, 115, 125]
OPERATIONS = ['reconv', 'noshear', '1p', '2p']

## ++++++++++++++ Workhorse

def psf_stamps():
    """PSF stamps of every size, centred on the true centre of the stamp."""
    saved = fits.getdata(os.path.join(image_dir, f'psf_{tile}', 'psf_ima_centred.fits')).astype(float)
    psf = ImSimPSF.MoffatPSF(SEEING, MOFFAT_BETA)
    drawn = ImSimPSF.PSFima(psf, PIXEL_SCALE, size=saved.shape[0], half_pixel_shift=False).array
    print(f'>>> ImSim redraws the saved {saved.shape[0]}-pixel PSF stamp to '
          f'{np.abs(drawn - saved).max():.1e} (peak {saved.max():.3e})')
    stamps = {}
    for size in PSF_SIZES:
        if size <= saved.shape[0]:
            o = (saved.shape[0] - size) // 2
            stamps[size] = saved[o:o + size, o:o + size]
        else:
            stamps[size] = ImSimPSF.PSFima(psf, PIXEL_SCALE, size=size, half_pixel_shift=False).array
    return stamps

def metacal_images(image, sigma, psf_stamp, cy, cx, N):
    """The metacal stages of an N-pixel square centred on (cy, cx), as metadetect builds them."""
    sl = (slice(cy - N // 2, cy + N // 2), slice(cx - N // 2, cx + N // 2))
    obs = ngmix.Observation(
        np.array(image[sl], dtype=float), weight=np.full((N, N), 1 / sigma**2),
        jacobian=ngmix.DiagonalJacobian(scale=PIXEL_SCALE, row=N / 2 - 0.5, col=N / 2 - 0.5),
        psf=ngmix.Observation(psf_stamp, jacobian=ngmix.DiagonalJacobian(
            scale=PIXEL_SCALE, row=(psf_stamp.shape[0] - 1) / 2, col=(psf_stamp.shape[1] - 1) / 2)))
    m = MetacalFitGaussPSF(obs, rng=np.random.RandomState(7))
    out = {}
    for name, shear in [('1p', ngmix.Shape(STEP, 0.)), ('2p', ngmix.Shape(0., STEP))]:
        _, psf_obj = m.get_target_psf(shear, 'gal_shear')
        out[name] = m.get_target_image(psf_obj, shear=shear).array
        if name == '1p':
            ## metacal takes noshear from the 1p call, with the same target PSF
            out['noshear'] = m.get_target_image(psf_obj, shear=None).array
            ## the square convolved with the target PSF, without deconvolving it first
            out['reconv'] = galsim.Convolve([m.image_int, psf_obj]).drawImage(
                nx=N, ny=N, wcs=m.image.wcs, dtype=np.float64, method='no_pixel').array
    return out

def inner(array, N):
    o = (array.shape[0] - N) // 2
    return array[o:o + N, o:o + N]

if __name__ == '__main__':
    image = fits.getdata(os.path.join(image_dir, f'{tile}_rot0.fits'), memmap=True)
    noise = fits.getdata(os.path.join(image_dir, f'noise_{tile}', 'noise_image.fits'), memmap=True)
    patch = np.asarray(noise[:2000, :2000])
    sigma = 1.4826 * np.median(np.abs(patch - np.median(patch)))
    stamps = psf_stamps()

    N_in = min(CELL_SIZES)
    yy, xx = np.mgrid[0:N_in, 0:N_in] - (N_in - 1) / 2.
    r = np.maximum(np.abs(xx), np.abs(yy))
    masks = [(r >= lo) & (r < hi) for lo, hi in zip(R_EDGES[:-1], R_EDGES[1:])]

    sum_d2 = {}
    for cy, cx in CENTRES:
        for psf_size, stamp in stamps.items():
            ref = metacal_images(image, sigma, stamp, cy, cx, REF_SIZE)
            for N in CELL_SIZES:
                out = metacal_images(image, sigma, stamp, cy, cx, N)
                for op in OPERATIONS:
                    d2 = (inner(out[op], N_in) - inner(ref[op], N_in))**2 / sigma**2
                    key = (psf_size, op, N)
                    sum_d2[key] = sum_d2.get(key, 0.) + np.array([d2[mk].mean() for mk in masks])

    print(f'\n>>> RMS of metacal(cell) - metacal({REF_SIZE}) over the central {N_in}x{N_in} pixels, '
          f'in noise sigma ({sigma:.3f}), {len(CENTRES)} points, by max(|dx|,|dy|) from the centre [pixels]')
    print('   PSF  operation  cell ' + ''.join(f'{lo:>6g}-{hi:<5g}' for lo, hi in zip(R_EDGES[:-1], R_EDGES[1:])))
    for (psf_size, op, N), v in sum_d2.items():
        print(f'   {psf_size:3d}  {op:9s} {N:5d} ' + ''.join(f'{x:11.2e} ' for x in np.sqrt(v / len(CENTRES))))
