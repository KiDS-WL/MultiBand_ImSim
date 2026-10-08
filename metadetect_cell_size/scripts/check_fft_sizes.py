# -*- coding: utf-8 -*-

### GalSim sizes inside metacal, recorded during real metacal calls on test-bench scenes of
###    the 9-arcsec grid (6 targets x 4 true shears, all five metacal images), per cell size:
###       table         size P of the cell's zero-padded InterpolatedImage table, good_fft_size(4N)
###       R_stepk       radius (px) behind the InterpolatedImage's stepk (flux-containing radius)
###       fft_<t>       size W of the drawImage FFT of metacal type t. drawImage centres the image
###                     bounds on 0, so W = good_fft_size(max(stepk-based size, N)), the
###                     stepk-based size being ~N + 10 px; it can vary between scenes
###       stepk_size_<t> that stepk-based size
###    One row per cell size, target, true shear and metacal type. Honours MCAL_FFT_SIZE and
###    MCAL_PAD_FACTOR (test_hard_cut.py). Saves fft_sizes<FFT_TAG>.csv.
###    About 10 seconds per cell size on one core.
###
###    Usage: python check_fft_sizes.py                          # the cells of cell_scan
###           CELL_SIZES=244,246 FFT_TAG=_map python check_fft_sizes.py

import os
os.environ['LAYOUTS'] = 'grid9'
os.environ.setdefault('CELL_SIZES', '240,250,260,270,300,350,400,450,500')
FFT_TAG = os.environ.get('FFT_TAG', '_cell_scan')
import numpy as np, galsim, pandas as pd
import ngmix.metacal.metacal as M
import test_layout as L, test_hard_cut as H

## record every drawImage FFT and InterpolatedImage of a cell
record = {'fft': [], 'ii': []}
_makeK = galsim.GSObject.drawFFT_makeKImage
def _makeK_rec(self, image):
    kimage, N = _makeK(self, image)
    if max(image.bounds.numpyShape()) >= H.MCAL_MIN_SIZE:
        record['fft'].append((N, self.getGoodImageSize(image.scale)))
    return kimage, N
galsim.GSObject.drawFFT_makeKImage = _makeK_rec
_stuff = M._galsim_stuff_impl
def _stuff_rec(img, wcs, xinterp):
    image, ii = _stuff(img, wcs, xinterp)
    if max(np.shape(img)) >= H.MCAL_MIN_SIZE:
        record['ii'].append((ii._xim.array.shape[0], np.pi / ii.stepk / H.PIXEL_SCALE))
    return image, ii
M._galsim_stuff_impl = _stuff_rec
M.USE_GALSIM_CACHE = False

rows = []
for i_target in range(6):
    r, offset, half_px = L.target_setup(i_target)
    dx, dy = L.layout_offsets('grid9', i_target, half_px)
    for gname, g in H.SHEARS.items():
        image, target = L.draw_layout(r, dx, dy, half_px, g, 0.)
        centre = np.floor(np.array(target) + offset) + 0.5
        for N in H.CELL_SIZES:
            cell, _ = H.cut_cell(image, centre, N)
            record['fft'].clear(); record['ii'].clear()
            H.metacal(cell, seed=1)
            ## one cell InterpolatedImage; one drawImage per metacal type, in the order of H.TYPES
            assert len(record['ii']) == 1 and len(record['fft']) == len(H.TYPES), \
                (len(record['ii']), len(record['fft']))
            table, R = record['ii'][0]
            out = dict(target=i_target, shear=gname, N=N, table=table, R_stepk=R)
            for t, (W, n_stepk) in zip(H.TYPES, record['fft']):
                out[f'fft_{t}'], out[f'stepk_size_{t}'] = W, n_stepk
            rows.append(out)
df = pd.DataFrame(rows)
fft_cols = [f'fft_{t}' for t in H.TYPES]
summary = df.groupby('N').agg(table=('table', 'max'), R_stepk_min=('R_stepk', 'min'), R_stepk_max=('R_stepk', 'max'))
summary['fft_sizes'] = df.groupby('N')[fft_cols].apply(lambda x: sorted(set(x.values.ravel())))
summary['frac_W_max'] = df.groupby('N')[fft_cols].apply(lambda x: (x.values == x.values.max()).mean()).round(2)
summary['P/W'] = [sorted({round(p / w, 3) for w in ws}) for p, ws in zip(summary['table'], summary['fft_sizes'])]
pd.set_option('display.width', 200)
print(summary.to_string())
df.to_csv(os.path.join(H.script_dir, f'fft_sizes{FFT_TAG}.csv'), index=False)
