# -*- coding: utf-8 -*-

### Does the layout of the neighbours decide the cell-250 bias?
###
###    The test bench of test_hard_cut.py part b (a galaxy within 20 px of a cell
###    centre, true shear +-0.02 in g1 or g2, 90-degree rotated twin, metacal on
###    cells of the sizes in CELL_SIZES, wmom at the known position,
###    m_i = S_i / R_ii - 1), with the neighbours laid out in different ways around
###    the same target galaxies (same galaxy and offset from the cell centre in
###    every layout and run):
###       grid9            9-arcsec grid, as in the simulation (44 galaxies / arcmin^2)
###       grid<s>          grid of spacing s arcsec, e.g. grid15, grid8.7
###       jitter9_<a>      9-arcsec grid, every neighbour moved at random by up to
###                        +-a pixels in each direction, e.g. jitter9_2
###       grid9_unsheared  9-arcsec grid whose positions are not sheared (the
###                        galaxies themselves still are)
###       random9          random positions, density of the 9-arcsec grid
###       random_cat       random positions, density of the input catalogue
###                        (galaxies of a 0.5 deg^2 tile, ~20 / arcmin^2)
###       isolated         no neighbours
###
###    All galaxies are drawn with ImSim's code in 256-px stamps (no cut, which
###    test_hard_cut.py showed makes no difference). Positions are sheared with
###    ImSim._shear_positions and shifted back so that the target stays in place;
###    a uniform shift changes nothing, and it keeps the +-g scenes paired.
###    grid9 repeats the scenes of the 'uncut' stamp of test_hard_cut.py.
###
###    m uses the galaxies that pass the cuts of analyse_cell_size.py in every
###    measurement of a layout, unweighted, so that the pairing cancels the
###    intrinsic shapes (see test_hard_cut.summarise). Every cell size is also
###    compared with a reference one (default the largest), on the same galaxies.
###
###    Configuration, from the environment so that worker processes see the same:
###       LAYOUTS      layout names separated by commas (default: the first study)
###       LAYOUT_TAG   added to the output names
###       CELL_SIZES   cell sizes separated by commas (default 250,500)
###       N_TARGETS    number of target galaxies (default 480)
###       REF_CELL     the cell size the others are compared with (default the largest)
###
###    Usage:
###       sbatch run_layout.sh                        # the first study, ~100 core-s per galaxy
###       sbatch run_layout_tests.sh <test>           # spacing_scan, fine_scan, jitter, cell_scan, map
###       python test_layout.py --summary-only        # summarise the saved results again
###                                                   #    (same environment as the run)
###    Plots: plot_layout_scans.py

import os
import sys
from multiprocessing import Pool

import numpy as np
import pandas as pd
import galsim

import test_hard_cut as H
import ImSim
import ImSimObject
from metadetect import shearpos

## ++++++++++++++ I/O and general setups

## density of the input catalogue: galaxies of one tile over its 0.5 deg^2
DENSITY_CAT = len(H.GALS) / (0.5 * 3600.**2)          # per arcsec^2
## (kind, parameter, jitter in pixels, positions sheared)
DEFAULT_LAYOUTS = {'grid9': ('grid', 9., 0., True), 'grid15': ('grid', 15., 0., True),
                   'grid18': ('grid', 18., 0., True),
                   'random9': ('random', 1. / 9.**2, 0., True),
                   'random_cat': ('random', DENSITY_CAT, 0., True),
                   'isolated': ('none', 0., 0., True)}

def layout_spec(name):
    if name in DEFAULT_LAYOUTS:
        return DEFAULT_LAYOUTS[name]
    if name == 'grid9_unsheared':
        return ('grid', 9., 0., False)
    if name.startswith('jitter9_'):
        return ('grid', 9., float(name[len('jitter9_'):]), True)
    if name.startswith('grid'):
        return ('grid', float(name[len('grid'):]), 0., True)
    raise ValueError(f'unknown layout {name!r}')

LAYOUTS = {name: layout_spec(name) for name in
           os.environ.get('LAYOUTS', ','.join(DEFAULT_LAYOUTS)).split(',')}
LAYOUT_TAG = os.environ.get('LAYOUT_TAG', '')
STAMP = 'uncut'
N_TARGETS = int(os.environ.get('N_TARGETS', 480))
## the scene of the first study covered cells of up to 500 px; never smaller, so
##    that the scenes (and the target galaxies) stay the same
HALF_PX_0 = 500 // 2 + H.MARGIN

outpath_resp = os.path.join(H.script_dir, f'layout_response{LAYOUT_TAG}.npz')
outpath_summary = os.path.join(H.script_dir, f'layout_response{LAYOUT_TAG}.csv')

## ++++++++++++++ Workhorse

def target_setup(i_target):
    """Catalogue rows and offset from the cell centre, drawn exactly as test_hard_cut.one_target."""
    rng = np.random.default_rng(1000 + i_target)
    K = int(np.ceil(HALF_PX_0 * H.PIXEL_SCALE / H.GRID_SIZE)) + 1
    rows = rng.integers(0, len(H.GALS), size=(2 * K + 1)**2)
    offset = rng.uniform(-H.OFFSET_MAX, H.OFFSET_MAX, size=2)
    half_px = max(HALF_PX_0, max(H.CELL_SIZES) // 2 + H.MARGIN)
    return rows, offset, half_px

def layout_offsets(layout, i_target, half_px):
    """Offsets (arcsec) of the galaxies from the target, the target first."""
    kind, par, jitter, _ = LAYOUTS[layout]
    half_arcsec = half_px * H.PIXEL_SCALE
    if kind == 'grid':
        K = int(np.ceil(half_px * H.PIXEL_SCALE / par)) + 1
        ii, jj = np.meshgrid(np.arange(-K, K + 1), np.arange(-K, K + 1))
        ii, jj = ii.ravel(), jj.ravel()
        order = np.argsort((ii != 0) | (jj != 0), kind='stable')
        dx, dy = ii[order] * par, jj[order] * par
        if jitter > 0:
            rng = np.random.default_rng(7000 + 97 * i_target + int(round(10 * jitter)))
            move = rng.uniform(-jitter, jitter, size=(2, len(dx))) * H.PIXEL_SCALE
            move[:, 0] = 0.                        # the target stays where it is
            dx, dy = dx + move[0], dy + move[1]
        return dx, dy
    if kind == 'random':
        rng = np.random.default_rng(5000 + 97 * i_target + list(DEFAULT_LAYOUTS).index(layout))
        n = rng.poisson(par * (2 * half_arcsec)**2)
        dx, dy = rng.uniform(-half_arcsec, half_arcsec, size=(2, n))
        return np.concatenate([[0.], dx]), np.concatenate([[0.], dy])
    return np.array([0.]), np.array([0.])

def draw_layout(rows, dx, dy, half_px, g, rot, shear_positions=True):
    """The scene, drawn with ImSim; returns the image and the target position (numpy row, col)."""
    gals = H.GALS.loc[rows[:len(dx)]].reset_index(drop=True).copy()
    gals['RA'] = H.RA0 + dx / 3600.
    gals['DEC'] = H.DEC0 + dy / 3600.
    if shear_positions:
        ImSim._shear_positions([gals, None], g)
        ## shift back so that the target stays at (RA0, DEC0)
        gals['RA'] -= gals.loc[0, 'RA'] - H.RA0
        gals['DEC'] -= gals.loc[0, 'DEC'] - H.DEC0
    half_deg = half_px * H.PIXEL_SCALE / 3600.
    canvas = ImSimObject.SimpleCanvas(H.RA0 - half_deg, H.RA0 + half_deg, H.DEC0 - half_deg, H.DEC0 + half_deg,
                                      H.PIXEL_SCALE, edge_sep=0.)
    image = ImSimObject.GalaxiesImage(canvas, H.BAND, H.PIXEL_SCALE, H.PSF, gals,
                                      gal_rotation_angle=rot, g_cosmic=list(g),
                                      gal_position_type=['grid', H.STAMPS[STAMP]])
    pos = canvas.wcs.toImage(galsim.CelestialCoord(H.RA0 * galsim.degrees, H.DEC0 * galsim.degrees))
    return image.array, (pos.y - image.bounds.ymin, pos.x - image.bounds.xmin)

def one_target(i_target):
    """e of one galaxy in every layout, cell size, true shear, rotation and metacal type."""
    rows, offset, half_px = target_setup(i_target)
    ## e[layout, cell, shear, rotation, type, (e1, e2, weight)]
    e = np.full((len(LAYOUTS), len(H.CELL_SIZES), len(H.SHEARS), len(H.ROTATIONS), len(H.TYPES), 3), np.nan)
    for i_l, layout in enumerate(LAYOUTS):
        dx, dy = layout_offsets(layout, i_target, half_px)
        ## dense grids need more galaxies than drawn for the 9-arcsec one; append
        ##    more, so that the first ones (and the 9-arcsec scenes) stay the same
        if len(dx) > len(rows):
            extra = np.random.default_rng(9000 + i_target).integers(0, len(H.GALS), size=len(dx) - len(rows))
            rows = np.concatenate([rows, extra])
        for i_g, g in enumerate(H.SHEARS.values()):
            for i_r, rot in enumerate(H.ROTATIONS):
                image, target = draw_layout(rows, dx, dy, half_px, g, rot,
                                            shear_positions=LAYOUTS[layout][3])
                centre = np.floor(np.array(target) + offset) + 0.5
                for i_n, N in enumerate(H.CELL_SIZES):
                    cell, (y0, x0) = H.cut_cell(image, centre, N)
                    mc = H.metacal(cell, seed=i_target)
                    row, col = target[0] - y0, target[1] - x0
                    for i_t, t in enumerate(H.TYPES):
                        im, jac, psf_im = mc[t]
                        rows_t, cols_t = shearpos.shear_positions(
                            np.atleast_1d(row), np.atleast_1d(col), t, jac, im.shape, step=H.STEP)
                        e[i_l, i_n, i_g, i_r, i_t] = H.wmom(im, float(rows_t[0]), float(cols_t[0]),
                                                            H.wmom_T(psf_im))
    return i_target, e

def response():
    N_proc = int(os.environ.get('SLURM_CPUS_PER_TASK', 4))
    results = {}
    with Pool(N_proc) as pool:
        for n_done, (i, e) in enumerate(pool.imap_unordered(one_target, range(N_TARGETS)), 1):
            results[i] = e
            if n_done % 40 == 0:
                print(f'>>> {n_done}/{N_TARGETS} galaxies done', flush=True)
    idx = sorted(results)
    np.savez_compressed(outpath_resp, e=np.array([results[i] for i in idx]),
                        layouts=np.array(list(LAYOUTS)), cell_sizes=np.array(H.CELL_SIZES),
                        shears=np.array(list(H.SHEARS)), types=np.array(H.TYPES))
    print(f'saved to {outpath_resp}')

def summarise(n_boot=2000):
    d = np.load(outpath_resp)
    e = d['e']                                   # [target, layout, cell, shear, rot, type, (e1, e2, w)]
    layouts, cells, shears, types = list(d['layouts']), list(d['cell_sizes']), list(d['shears']), list(d['types'])
    iT = {t: i for i, t in enumerate(types)}
    iG = {s: i for i, s in enumerate(shears)}
    COMPS = {1: ('g1p', 'g1m', '1p', '1m'), 2: ('g2p', 'g2m', '2p', '2m')}
    i_ref = cells.index(int(os.environ['REF_CELL'])) if 'REF_CELL' in os.environ else int(np.argmax(cells))
    rng = np.random.default_rng(5)

    def SRm(x):
        """S, R, m [cell, comp-1] from the measurements x[target, cell, shear, rot, type, 3], unweighted"""
        out = np.zeros((3, len(cells), 2))
        for i_n in range(len(cells)):
            for c, (gp, gm, tp, tm) in COMPS.items():
                E = lambda g, t: x[:, i_n, iG[g], :, iT[t], c - 1].mean()
                S = (E(gp, 'noshear') - E(gm, 'noshear')) / (2 * H.G_TRUE)
                R = np.mean([(E(g, tp) - E(g, tm)) / (2 * H.STEP) for g in (gp, gm)])
                out[:, i_n, c - 1] = S, R, S / R - 1
        return out

    rows = []
    for i_l, layout in enumerate(layouts):
        x = e[:, i_l]
        ## the galaxies passing the cuts in every measurement of this layout
        keep = np.all(x[..., 2] > 0, axis=tuple(range(1, x.ndim - 1)))
        x = x[keep]
        full = SRm(x)
        boot = rng.integers(0, len(x), size=(n_boot, len(x)))
        bs = np.array([SRm(x[b]) for b in boot])
        for i_n, N in enumerate(cells):
            row = dict(layout=layout, cell=str(N), N_gal=int(keep.sum()))
            for c in (1, 2):
                row.update({f'S{c}': full[0, i_n, c - 1], f'R{c}{c}': full[1, i_n, c - 1],
                            f'm{c}': full[2, i_n, c - 1], f'm{c}_err': bs[:, 2, i_n, c - 1].std()})
            rows.append(row)
        ## every cell size against the reference one, on the same galaxies
        for i_n, N in enumerate(cells):
            if i_n == i_ref:
                continue
            row = dict(layout=layout, cell=f'{N}-{cells[i_ref]}', N_gal=int(keep.sum()))
            for c in (1, 2):
                row.update({f'm{c}': full[2, i_n, c - 1] - full[2, i_ref, c - 1],
                            f'm{c}_err': (bs[:, 2, i_n, c - 1] - bs[:, 2, i_ref, c - 1]).std()})
            rows.append(row)
    out = pd.DataFrame(rows)
    out.to_csv(outpath_summary, index=False)
    print(f'\n>>> galaxies within {H.OFFSET_MAX:g} px of the cell centre, noise-free scenes, '
          f'galaxies passing the cuts in every measurement of a layout, unweighted; '
          f'm_i = S_i / R_ii - 1; rows "<cell>-{cells[i_ref]}": m(cell) - m(cell {cells[i_ref]})')
    with pd.option_context('display.width', 200, 'display.max_columns', None):
        print(out.to_string(index=False, float_format='%.4f', na_rep=''))
    print(f'saved to {outpath_summary}')

if __name__ == '__main__':
    if '--summary-only' not in sys.argv:
        response()
    summarise()
