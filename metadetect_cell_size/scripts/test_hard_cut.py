# -*- coding: utf-8 -*-

### Is the hard cut of the grid stamps what biases metacal in 250-pixel cells?
###
###    In grid mode ImSim draws every galaxy, already convolved with the PSF, into a
###    square stamp of 2*floor(grid_size/pixel_scale/2) pixels and drops the light
###    outside it: 44 px for the 9-arcsec grid. Such a cut image is no longer
###    "sky convolved with the PSF", which metacal's deconvolution assumes, and the
###    cut is an axis-aligned square that the true shear leaves alone but metacal's
###    shear distorts.
###
###    Synthetic, noise-free scenes are drawn with ImSim's own code
###    (ImSimObject.GalaxiesImage, ImSimPSF.MoffatPSF, ImSim._shear_positions):
###    galaxies of the input catalogue of the simulation on the 9-arcsec grid,
###    the simulation's PSF, positions sheared as well. Only the stamp changes:
###       44px    the 9-arcsec grid stamp, as in the simulation
###       74px    the stamp of a 15-arcsec grid, on the same 9-arcsec positions
###       uncut   a 256-px stamp, holding essentially all the light
###
###    Part a (convergence): metacal images of 250- and 500-px cells against a
###       1000-px one, as test_metacal_convergence.py, for each stamp.
###    Part b (response): a galaxy of the catalogue is put within 20 px of the
###       centre of a cell, with grid neighbours around it, and sheared by
###       +-0.02 in g1 or g2, with its 90-degree rotated twin. Metacal is run on the
###       250- and 500-px cells, and the galaxy measured with wmom (1.2 arcsec) on a
###       32-px stamp at its known position, as metadetect does after detection.
###       The response to the true shear S is compared with metacal's R:
###          m_i = S_i / R_ii - 1
###       with errors from resampling the galaxies (everything is paired).
###       Cuts and weights are those of analyse_cell_size.py (S/N > 12.5,
###       T/T_psf > 1.2, flags == 0, shear_weight), with S/N and the shape
###       covariance computed for the noise level of the simulation, although the
###       scenes themselves are noise-free.
###
###    Usage:
###       python test_hard_cut.py a        # convergence, a few minutes on one core
###       python test_hard_cut.py b        # response, ~2 core-hours (run_hard_cut.sh)
###       python test_hard_cut.py b --summary-only   # summarise the saved part b again

import os
import sys
import logging
from multiprocessing import Pool

import numpy as np
import pandas as pd
import galsim
import ngmix

script_dir = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(script_dir, '..', '..', 'modules'))
sys.path.insert(0, os.path.join(script_dir, '..', '..', 'modules', 'metadetect'))
import analyse_cell_size as A
import ImSim
import ImSimObject
import ImSimPSF
from metadetect import shearpos

for name in ('ImSim', 'ImSimObject', 'ImSimPSF'):
    logging.getLogger(name).setLevel(logging.WARNING)

## ++++++++++++++ I/O and general setups

## galaxies of the simulation (input catalogue of one tile, zero-shear tag)
input_path = ('/sdf/data/kipac/u/jlkitt/imsims/grid/seeing_073/p000p000/catalogues/input/'
              'gals_info_tile0.0_0.0.feather')
PIXEL_SCALE = 0.2
GRID_SIZE = 9.               # arcsec between neighbouring galaxies
MAG_ZERO = 30.
BAND = 'LSST_r'
SEEING, MOFFAT_BETA = 0.73, 2.22406805360007
## noise level of the simulation, only to express image differences in its units
SIGMA_NOISE = 0.311

## the stamp is set through the grid size ImSim derives it from: 2*floor(grid/0.2/2)
STAMPS = {'44px': 9., '74px': 15., 'uncut': 51.2}
## galaxies drawn this far beyond the largest cell, so stamps reaching in are included
MARGIN = 160
## the scene sits here on the sky, away from RA = 0
RA0, DEC0 = 1.0, 0.3

STEP = 0.01                  # metacal shear step
G_TRUE = 0.02
SHEARS = {'g1p': (G_TRUE, 0.), 'g1m': (-G_TRUE, 0.), 'g2p': (0., G_TRUE), 'g2m': (0., -G_TRUE)}
ROTATIONS = (0., 90.)
TYPES = ('noshear', '1p', '1m', '2p', '2m')
BOX = 32                     # metadetect stamp size (meds min = max = 32)
WMOM_FWHM = 1.2

## part a
CONV_SIZES, CONV_REF = (250, 500), 1000
CONV_SCENES = 4
R_EDGES = [0, 12.5, 25, 50, 75, 100, 115, 125]
## part b
## from the environment (comma-separated), so that worker processes see the same;
##    the largest one is the reference the others are compared with
CELL_SIZES = tuple(int(x) for x in os.environ.get('CELL_SIZES', '250,500').split(','))
N_TARGETS = int(os.environ.get('N_TARGETS', 240))
OFFSET_MAX = 20.             # target within this many pixels of the cell centre (per axis)

## ++++++++++++++ Optional changes to GalSim inside metacal (FFT-size tests)
##    From the environment, so that worker processes see the same. They apply only inside
##    metacal() (not to the drawing of the scenes) and to images of at least MCAL_MIN_SIZE
##    pixels (the cells, not the PSF stamps):
##       MCAL_FFT_SIZE    size of the FFT of drawImage (default: GalSim's choice, which is
##                        good_fft_size(2N) for these N x N cells)
##       MCAL_PAD_FACTOR  pad_factor of the cell's InterpolatedImage (GalSim's default 4,
##                        a table of good_fft_size(4N))
MCAL_MIN_SIZE = 200
MCAL_FFT_SIZE = int(os.environ['MCAL_FFT_SIZE']) if os.environ.get('MCAL_FFT_SIZE') else None
MCAL_PAD_FACTOR = float(os.environ['MCAL_PAD_FACTOR']) if os.environ.get('MCAL_PAD_FACTOR') else None
_IN_METACAL = False          # set by metacal()

if MCAL_FFT_SIZE is not None:
    _drawFFT_makeKImage = galsim.GSObject.drawFFT_makeKImage

    def _drawFFT_makeKImage_fixed(self, image):
        """GSObject.drawFFT_makeKImage (GalSim 2.8) with the FFT size set to MCAL_FFT_SIZE."""
        if not _IN_METACAL or max(image.bounds.numpyShape()) < MCAL_MIN_SIZE:
            return _drawFFT_makeKImage(self, image)
        N = MCAL_FFT_SIZE
        dk = 2. * np.pi / (N * image.scale)
        Nk = N if N * dk / 2 > self.maxk else int(np.ceil(self.maxk / dk)) * 2
        kimage = galsim.ImageCD(bounds=galsim.BoundsI(0, Nk // 2, -Nk // 2, Nk // 2), scale=dk)
        return kimage, N

    galsim.GSObject.drawFFT_makeKImage = _drawFFT_makeKImage_fixed

if MCAL_PAD_FACTOR is not None:
    import ngmix.metacal.metacal as _ngmix_metacal

    def _galsim_stuff_impl(img, wcs, xinterp):
        """ngmix.metacal.metacal._galsim_stuff_impl (ngmix 2.4) with pad_factor MCAL_PAD_FACTOR."""
        image = galsim.Image(img, wcs=wcs)
        pad = MCAL_PAD_FACTOR if max(np.shape(img)) >= MCAL_MIN_SIZE else 4.
        return image, galsim.InterpolatedImage(image, x_interpolant=xinterp, pad_factor=pad)

    _ngmix_metacal._galsim_stuff_impl = _galsim_stuff_impl

outpath_conv = os.path.join(script_dir, 'hard_cut_convergence.csv')
outpath_resp = os.path.join(script_dir, 'hard_cut_response.npz')
outpath_summary = os.path.join(script_dir, 'hard_cut_response.csv')

## ++++++++++++++ Workhorse

def load_galaxies():
    gals = pd.read_feather(input_path)
    gals = gals.rename(columns={c: c[:-len('_input')] for c in gals.columns if c.endswith('_input')})
    gals[BAND] = 10**(-0.4 * (gals[BAND] - MAG_ZERO))
    return gals.reset_index(drop=True)

GALS = load_galaxies()
PSF = ImSimPSF.MoffatPSF(SEEING, MOFFAT_BETA)
PSF_STAMP = ImSimPSF.PSFima(PSF, PIXEL_SCALE, size=48, half_pixel_shift=False).array

def draw_scene(rows, half_px, stamp, g, rot):
    """
    Galaxies on the 9-arcsec grid around the target (rows[0], at the centre),
    positions sheared by g about the target, drawn with ImSim into a canvas of
    +-half_px pixels around the target. Returns the image array and the target's
    position in it (numpy row, col).
    """
    apart = GRID_SIZE / 3600.
    K = int(np.ceil(half_px * PIXEL_SCALE / GRID_SIZE)) + 1
    ii, jj = np.meshgrid(np.arange(-K, K + 1), np.arange(-K, K + 1))
    ii, jj = ii.ravel(), jj.ravel()
    ## the target first
    order = np.argsort((ii != 0) | (jj != 0), kind='stable')
    gals = GALS.loc[rows[:len(order)]].reset_index(drop=True).copy()
    gals['RA'] = RA0 + ii[order] * apart
    gals['DEC'] = DEC0 + jj[order] * apart
    ImSim._shear_positions([gals, None], g)

    half_deg = half_px * PIXEL_SCALE / 3600.
    canvas = ImSimObject.SimpleCanvas(RA0 - half_deg, RA0 + half_deg, DEC0 - half_deg, DEC0 + half_deg,
                                      PIXEL_SCALE, edge_sep=0.)
    image = ImSimObject.GalaxiesImage(canvas, BAND, PIXEL_SCALE, PSF, gals,
                                      gal_rotation_angle=rot, g_cosmic=list(g),
                                      gal_position_type=['grid', STAMPS[stamp]])
    pos = canvas.wcs.toImage(galsim.CelestialCoord(RA0 * galsim.degrees, DEC0 * galsim.degrees))
    return image.array, (pos.y - image.bounds.ymin, pos.x - image.bounds.xmin)

def metacal(cell, seed):
    """Metacal images of a cell, as metadetect makes them (noise-free, so no noise fixing)."""
    N = cell.shape[0]
    obs = ngmix.Observation(
        np.array(cell, dtype=float), weight=np.full(cell.shape, 1. / SIGMA_NOISE**2),
        jacobian=ngmix.DiagonalJacobian(scale=PIXEL_SCALE, row=(N - 1) / 2, col=(N - 1) / 2),
        psf=ngmix.Observation(PSF_STAMP, jacobian=ngmix.DiagonalJacobian(
            scale=PIXEL_SCALE, row=(PSF_STAMP.shape[0] - 1) / 2, col=(PSF_STAMP.shape[1] - 1) / 2)))
    global _IN_METACAL
    _IN_METACAL = True
    try:
        out = ngmix.metacal.get_all_metacal(obs, psf='fitgauss', step=STEP, fixnoise=False,
                                            rng=np.random.RandomState(seed), types=list(TYPES))
    finally:
        _IN_METACAL = False
    return {t: (out[t].image, out[t].jacobian, out[t].psf.image) for t in TYPES}

def cut_cell(image, centre, N):
    """N x N cell whose true centre ((N-1)/2) is at the half-integer position centre."""
    y0, x0 = int(round(centre[0] - (N - 1) / 2)), int(round(centre[1] - (N - 1) / 2))
    return image[y0:y0 + N, x0:x0 + N], (y0, x0)

def wmom_T(psf_image):
    """wmom size T of the (metacal target) PSF, for T_ratio."""
    n = psf_image.shape[0]
    obs = ngmix.Observation(psf_image, jacobian=ngmix.DiagonalJacobian(
        scale=PIXEL_SCALE, row=(n - 1) / 2, col=(n - 1) / 2))
    return ngmix.gaussmom.GaussMom(fwhm=WMOM_FWHM).go(obs)['T']

def wmom(image, row, col, T_psf):
    """
    wmom of the object at (row, col) on a BOX-pixel stamp, cut as metadetect's MEDS
    does, and its shear weight as analyse_cell_size.shear_weight (zero if cut).
    Returns (e1, e2, weight).
    """
    r0, c0 = int(row) - BOX // 2 + 1, int(col) - BOX // 2 + 1
    stamp = image[r0:r0 + BOX, c0:c0 + BOX]
    obs = ngmix.Observation(stamp, weight=np.full(stamp.shape, 1. / SIGMA_NOISE**2),
                            jacobian=ngmix.DiagonalJacobian(scale=PIXEL_SCALE, row=row - r0, col=col - c0))
    res = ngmix.gaussmom.GaussMom(fwhm=WMOM_FWHM).go(obs)
    if res['flags'] != 0:
        return np.array([np.nan, np.nan, 0.])
    e = np.array(res['e'])
    cata = pd.DataFrame({f'{A.fit_model}_g_1': [e[0]], f'{A.fit_model}_g_2': [e[1]],
                         f'{A.fit_model}_g_cov_as_sigma': [np.sqrt(np.trace(res['e_cov']) / 2.)],
                         f'{A.fit_model}_s2n': [res['s2n']], f'{A.fit_model}_T_ratio': [res['T'] / T_psf],
                         f'{A.fit_model}_flags': [res['flags']]})
    return np.array([e[0], e[1], A.shear_weight(cata)[0]])

## ---- part a

def convergence():
    rng = np.random.default_rng(11)
    half_px = CONV_REF // 2 + MARGIN
    n_in = min(CONV_SIZES)
    yy, xx = np.mgrid[0:n_in, 0:n_in] - (n_in - 1) / 2.
    r = np.maximum(np.abs(xx), np.abs(yy))
    masks = [(r >= lo) & (r < hi) for lo, hi in zip(R_EDGES[:-1], R_EDGES[1:])]
    sum_d2 = {}
    for i_scene in range(CONV_SCENES):
        rows = rng.integers(0, len(GALS), size=(2 * (half_px // 45 + 3) + 1)**2)
        offset = rng.uniform(-OFFSET_MAX, OFFSET_MAX, size=2)
        for stamp in STAMPS:
            image, target = draw_scene(rows, half_px, stamp, (0., 0.), 0.)
            centre = np.floor(np.array(target) + offset) + 0.5
            mc = {N: metacal(cut_cell(image, centre, N)[0], seed=i_scene)
                  for N in CONV_SIZES + (CONV_REF,)}
            for N in CONV_SIZES:
                for t in ('noshear', '1p', '2p'):
                    o = (N - n_in) // 2
                    a = mc[N][t][0][o:o + n_in, o:o + n_in]
                    o_ref = (CONV_REF - n_in) // 2
                    b = mc[CONV_REF][t][0][o_ref:o_ref + n_in, o_ref:o_ref + n_in]
                    d2 = (a - b)**2 / SIGMA_NOISE**2
                    key = (stamp, t, N)
                    sum_d2[key] = sum_d2.get(key, 0.) + np.array([d2[mk].mean() for mk in masks])
        print(f'>>> scene {i_scene + 1}/{CONV_SCENES} done', flush=True)
    rows_out = []
    for (stamp, t, N), v in sum_d2.items():
        rms = np.sqrt(v / CONV_SCENES)
        rows_out.append(dict(stamp=stamp, type=t, cell=N,
                             **{f'{lo:g}-{hi:g}': x for (lo, hi), x in zip(zip(R_EDGES[:-1], R_EDGES[1:]), rms)}))
    out = pd.DataFrame(rows_out)
    out.to_csv(outpath_conv, index=False)
    print(f'\n>>> RMS of metacal(cell) - metacal({CONV_REF}) over the central {n_in}x{n_in} pixels, '
          f'in units of the simulation noise ({SIGMA_NOISE}), noise-free scenes, by max(|dx|,|dy|) [px]')
    with pd.option_context('display.width', 200):
        print(out.to_string(index=False, float_format='%.2e'))
    print(f'saved to {outpath_conv}')

## ---- part b

def one_target(i_target):
    """e of one galaxy in every stamp, cell size, true shear, rotation and metacal type."""
    rng = np.random.default_rng(1000 + i_target)
    half_px = max(CELL_SIZES) // 2 + MARGIN
    K = int(np.ceil(half_px * PIXEL_SCALE / GRID_SIZE)) + 1
    rows = rng.integers(0, len(GALS), size=(2 * K + 1)**2)
    offset = rng.uniform(-OFFSET_MAX, OFFSET_MAX, size=2)
    ## e[stamp, cell, shear, rotation, type, (e1, e2, weight)]
    e = np.full((len(STAMPS), len(CELL_SIZES), len(SHEARS), len(ROTATIONS), len(TYPES), 3), np.nan)
    for i_s, stamp in enumerate(STAMPS):
        for i_g, g in enumerate(SHEARS.values()):
            for i_r, rot in enumerate(ROTATIONS):
                image, target = draw_scene(rows, half_px, stamp, g, rot)
                centre = np.floor(np.array(target) + offset) + 0.5
                for i_n, N in enumerate(CELL_SIZES):
                    cell, (y0, x0) = cut_cell(image, centre, N)
                    mc = metacal(cell, seed=i_target)
                    row, col = target[0] - y0, target[1] - x0
                    for i_t, t in enumerate(TYPES):
                        im, jac, psf_im = mc[t]
                        rows_t, cols_t = shearpos.shear_positions(
                            np.atleast_1d(row), np.atleast_1d(col), t, jac, im.shape, step=STEP)
                        e[i_s, i_n, i_g, i_r, i_t] = wmom(im, float(rows_t[0]), float(cols_t[0]),
                                                          wmom_T(psf_im))
    return i_target, e, GALS.loc[rows[0], ['Re', 'sersic_n', BAND]].values.astype(float), offset

def response():
    N_proc = int(os.environ.get('SLURM_CPUS_PER_TASK', 4))
    results = {}
    with Pool(N_proc) as pool:
        for n_done, (i, e, props, offset) in enumerate(pool.imap_unordered(one_target, range(N_TARGETS)), 1):
            results[i] = (e, props, offset)
            if n_done % 20 == 0:
                print(f'>>> {n_done}/{N_TARGETS} galaxies done', flush=True)
    idx = sorted(results)
    np.savez_compressed(outpath_resp,
                        e=np.array([results[i][0] for i in idx]),
                        target_props=np.array([results[i][1] for i in idx]),
                        offsets=np.array([results[i][2] for i in idx]),
                        stamps=np.array(list(STAMPS)), cell_sizes=np.array(CELL_SIZES),
                        shears=np.array(list(SHEARS)), types=np.array(TYPES))
    print(f'saved to {outpath_resp}')

def summarise(n_boot=2000):
    """
    m per stamp and cell size, and paired differences between cell sizes and between
    stamps (same galaxies, scenes and resampling, so most of the noise cancels).

    Two ways of averaging over galaxies:
       weights  the shear_weight of analyse_cell_size.py. It depends on the measured
                |e|, so it differs between the +g and -g scenes of a galaxy and the
                pairing no longer cancels the intrinsic shapes: noisy.
       cuts     the same cuts, unweighted. A galaxy near a cut can pass in its +g scene
                and fail in its -g one, which again breaks the pairing: noisy.
       fixed    unweighted, only the galaxies passing the cuts in every measurement
                (all stamps, cells, shears, rotations and metacal types), so every
                average is over the same galaxies and the pairing cancels the
                intrinsic shapes. Measures the shape response without the selection
                response.
    """
    d = np.load(outpath_resp)
    e = d['e']                                   # [target, stamp, cell, shear, rot, type, (e1, e2, w)]
    stamps, cells, shears, types = list(d['stamps']), list(d['cell_sizes']), list(d['shears']), list(d['types'])
    iT = {t: i for i, t in enumerate(types)}
    iG = {s: i for i, s in enumerate(shears)}
    rng = np.random.default_rng(5)
    boot = rng.integers(0, e.shape[0], size=(n_boot, e.shape[0]))
    COMPS = {1: ('g1p', 'g1m', '1p', '1m'), 2: ('g2p', 'g2m', '2p', '2m')}

    def m_values(x, mode):
        """m[stamp, cell, comp-1] (and S, R) from the measurements x of a set of galaxies"""
        if mode == 'weights':
            w = x[..., 2]
        elif mode == 'cuts':
            w = (x[..., 2] > 0).astype(float)
        else:
            keep = np.all(x[..., 2] > 0, axis=tuple(range(1, x.ndim - 1)))
            w = np.broadcast_to(keep.reshape((-1,) + (1,) * (x.ndim - 2)), x.shape[:-1]).astype(float)
        out = np.zeros((3, len(stamps), len(cells), 2))      # S, R, m
        for i_s in range(len(stamps)):
            for i_n in range(len(cells)):
                for c, (gp, gm, tp, tm) in COMPS.items():
                    def E(g, t):
                        xx = x[:, i_s, i_n, iG[g], :, iT[t], c - 1].ravel()
                        ww = w[:, i_s, i_n, iG[g], :, iT[t]].ravel()
                        return np.sum(np.where(ww > 0, xx * ww, 0.)) / np.sum(ww)
                    S = (E(gp, 'noshear') - E(gm, 'noshear')) / (2 * G_TRUE)
                    R = np.mean([(E(g, tp) - E(g, tm)) / (2 * STEP) for g in (gp, gm)])
                    out[:, i_s, i_n, c - 1] = S, R, S / R - 1
        return out

    rows = []
    for mode in ('fixed', 'cuts', 'weights'):
        full = m_values(e, mode)
        bs = np.array([m_values(e[b], mode) for b in boot])  # [boot, 3, stamp, cell, comp]
        for i_s, stamp in enumerate(stamps):
            for i_n, N in enumerate(cells):
                row = dict(average=mode, stamp=stamp, cell=int(N))
                for c in (1, 2):
                    row.update({f'S{c}': full[0, i_s, i_n, c - 1], f'R{c}{c}': full[1, i_s, i_n, c - 1],
                                f'm{c}': full[2, i_s, i_n, c - 1], f'm{c}_err': bs[:, 2, i_s, i_n, c - 1].std()})
                rows.append(row)
            ## the same galaxies in both cell sizes
            for c in (1, 2):
                dm = full[2, i_s, 0, c - 1] - full[2, i_s, 1, c - 1]
                dm_err = (bs[:, 2, i_s, 0, c - 1] - bs[:, 2, i_s, 1, c - 1]).std()
                rows.append(dict(average=mode, stamp=stamp, cell=f'{cells[0]}-{cells[1]}',
                                 **{f'm{c}': dm, f'm{c}_err': dm_err}))
        ## the same galaxies with and without the stamp cut
        i_cut, i_unc = stamps.index('44px'), stamps.index('uncut')
        for i_n, N in enumerate(cells):
            for c in (1, 2):
                dm = full[2, i_cut, i_n, c - 1] - full[2, i_unc, i_n, c - 1]
                dm_err = (bs[:, 2, i_cut, i_n, c - 1] - bs[:, 2, i_unc, i_n, c - 1]).std()
                rows.append(dict(average=mode, stamp='44px-uncut', cell=int(N),
                                 **{f'm{c}': dm, f'm{c}_err': dm_err}))
    out = pd.DataFrame(rows)
    ## one row per (average, stamp, cell), m1 and m2 side by side
    out = out.groupby(['average', 'stamp', 'cell'], sort=False).first().reset_index()
    out.to_csv(outpath_summary, index=False)
    n_gal = int(np.round(np.mean(e[:, :, :, :, :, iT['noshear'], 2] > 0) * e.shape[0]))
    n_fixed = int(np.sum(np.all(e[..., 2] > 0, axis=tuple(range(1, e.ndim - 1)))))
    print(f'\n>>> {e.shape[0]} galaxies within {OFFSET_MAX:g} px of the cell centre (~{n_gal} pass the cuts, '
          f'{n_fixed} pass them in every measurement), '
          f'noise-free scenes, m_i = S_i / R_ii - 1 (S: response to the true shear, R: metacal response)')
    print('    rows "250-500": m(cell 250) - m(cell 500); rows "44px-uncut": m(44-px stamp) - m(uncut)')
    with pd.option_context('display.width', 200, 'display.max_columns', None):
        print(out.to_string(index=False, float_format='%.4f', na_rep=''))
    print(f'saved to {outpath_summary}')

if __name__ == '__main__':
    part = sys.argv[1] if len(sys.argv) > 1 else 'a'
    if part == 'a':
        convergence()
    elif part == 'b':
        if '--summary-only' not in sys.argv:
            response()
        summarise()
    else:
        raise SystemExit(f'unknown part {part!r}, use a or b')
