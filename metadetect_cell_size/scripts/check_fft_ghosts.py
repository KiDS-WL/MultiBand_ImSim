# -*- coding: utf-8 -*-

### Are there shifted copies ("ghosts") of the cell in metacal images whose drawImage FFT W
###    does not divide the cell's image table P?
###
###    GalSim interpolates the cell's InterpolatedImage in k-space from a table of P px
###    (P = good_fft_size(4N)); in real space that leaves faint copies of the cell P px away.
###    drawImage folds the result with period W, so a copy lands P mod W from the cell: onto the
###    cell itself when P/W is an integer, elsewhere otherwise. For a 250-px cell (P = 1024):
###    W = 384 puts copies at +-128 px, inside the cell; W = 768 at +-256 px, just outside it.
###
###    One test-bench scene (9-arcsec grid, target 0, true shear g1p, rotation 0), metacal of a
###    250-px cell with W forced to 384, 512 and 768 (scene unchanged). The differences from
###    W = 512 (P/W = 2, clean) are cross-correlated with the W = 512 image over all shifts.
###    Then, for targets 0-2, the difference W = 384 - W = 512 in the inner 120 x 120 px is
###    fitted with the W = 512 image itself and with copies of it shifted by +-128 px (in x, y
###    and both), and the inner rms of W = 768 - W = 512 is given (its copies, +-256 px away,
###    cannot overlap the 250-px cell). Saves metacal_fft_ghosts.png. About 3 minutes on one core.
###
###    Usage: python check_fft_ghosts.py

import os
os.environ['LAYOUTS'] = 'grid9'
os.environ['CELL_SIZES'] = '250'
import numpy as np
import galsim
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import test_layout as L
import test_hard_cut as H

N, W_REF, W_TEST = 250, 512, (384, 768)
TYPE = os.environ.get('GHOST_TYPE', 'noshear')

## drawImage FFT size of the metacal images of the cell, switched between calls
FORCE = {'W': None}
_makeK = galsim.GSObject.drawFFT_makeKImage
def _makeK_forced(self, image):
    if FORCE['W'] is None or max(image.bounds.numpyShape()) < H.MCAL_MIN_SIZE:
        return _makeK(self, image)
    n = FORCE['W']
    dk = 2. * np.pi / (n * image.scale)
    nk = n if n * dk / 2 > self.maxk else int(np.ceil(self.maxk / dk)) * 2
    return galsim.ImageCD(bounds=galsim.BoundsI(0, nk // 2, -nk // 2, nk // 2), scale=dk), n
galsim.GSObject.drawFFT_makeKImage = _makeK_forced

def metacal_images(i_target):
    """metacal TYPE image of the 250-px cell of one scene, for each forced W"""
    rows, offset, half_px = L.target_setup(i_target)
    dx, dy = L.layout_offsets('grid9', i_target, half_px)
    image, target = L.draw_layout(rows, dx, dy, half_px, H.SHEARS['g1p'], 0.)
    cell, _ = H.cut_cell(image, np.floor(np.array(target) + offset) + 0.5, N)
    out = {}
    for W in (W_REF,) + W_TEST:
        FORCE['W'] = W
        out[W] = H.metacal(cell, seed=1)[TYPE][0]
    FORCE['W'] = None
    return out

def shifted(a, dy, dx):
    """a shifted by (dy, dx), zero-filled"""
    out, n = np.zeros_like(a), a.shape[0]
    out[max(dy, 0):n + min(dy, 0), max(dx, 0):n + min(dx, 0)] = \
        a[max(-dy, 0):n + min(-dy, 0), max(-dx, 0):n + min(-dx, 0)]
    return out

mc = metacal_images(0)

def xcorr(a, b):
    """Normalised cross-correlation of a with b shifted by (dy, dx), all shifts, zero-padded."""
    n = a.shape[0]
    fa, fb = np.fft.rfft2(a, s=(2 * n, 2 * n)), np.fft.rfft2(b, s=(2 * n, 2 * n))
    c = np.fft.fftshift(np.fft.irfft2(fa * np.conj(fb), s=(2 * n, 2 * n)))
    return c / np.sqrt((a**2).sum() * (b**2).sum())     # shift (0, 0) at index (n, n)

ref = mc[W_REF]
fig, axs = plt.subplots(2, 3, figsize=(13, 8.4), layout='constrained')
im = axs[0, 0].imshow(np.arcsinh(ref / H.SIGMA_NOISE), origin='lower', cmap='gray_r')
axs[0, 0].set_title(f'metacal {TYPE}, W = {W_REF} (P/W = 2)\narcsinh(image / noise sigma)', fontsize=9)
lags = np.arange(-N, N)
for j, W in enumerate(W_TEST, 1):
    diff = mc[W] - ref
    v = np.abs(diff).max() / H.SIGMA_NOISE
    axs[0, j].imshow(diff / H.SIGMA_NOISE, origin='lower', cmap='RdBu_r', vmin=-v, vmax=v)
    axs[0, j].set_title(f'W = {W} (P/W = {1024 / W:.3g}) minus W = {W_REF}, in noise sigma\n'
                        f'rms {diff.std() / H.SIGMA_NOISE:.1e}, max {v:.1e}', fontsize=9)
    c = xcorr(diff, ref)
    axs[1, j].imshow(c, origin='lower', cmap='RdBu_r', vmin=-np.abs(c).max(), vmax=np.abs(c).max(),
                     extent=(lags[0] - 0.5, lags[-1] + 0.5, lags[0] - 0.5, lags[-1] + 0.5))
    shift = (1024 % W if 1024 % W <= W / 2 else 1024 % W - W)
    for s in (shift, -shift):
        for ax_line in (axs[1, j].axvline, axs[1, j].axhline):
            ax_line(s, color='k', lw=0.6, ls='--')
    axs[1, j].set_title(f'cross-correlation of the difference with the W = {W_REF} image\n'
                        f'dashed: shift P mod W = ±{abs(shift)} px', fontsize=9)
    axs[1, j].set_xlabel('shift in x [px]')
    axs[1, j].set_ylabel('shift in y [px]')
    ## strongest correlations
    flat = np.argsort(np.abs(c).ravel())[::-1][:6]
    peaks = [(int(lags[k // (2 * N)]), int(lags[k % (2 * N)]), float(c.ravel()[k])) for k in flat]
    print(f'W = {W}: rms diff {diff.std() / H.SIGMA_NOISE:.2e} sigma; strongest correlations (dy, dx, r):',
          [(p[0], p[1], round(p[2], 3)) for p in peaks])
## radial profile of the difference rms from the cell centre
yy, xx = np.indices(ref.shape) - (N - 1) / 2
r = np.maximum(np.abs(yy), np.abs(xx))
edges = np.arange(0, N // 2 + 1, 12.5)
for W in W_TEST:
    diff = (mc[W] - ref) / H.SIGMA_NOISE
    prof = [diff[(r >= a) & (r < b)].std() for a, b in zip(edges[:-1], edges[1:])]
    axs[1, 0].semilogy(0.5 * (edges[:-1] + edges[1:]), prof, marker='o', label=f'W = {W}')
axs[1, 0].set_xlabel('max(|dx|, |dy|) from the cell centre [px]')
axs[1, 0].set_ylabel(f'rms of (W − {W_REF}) [noise sigma]')
axs[1, 0].legend(frameon=False)
fig.suptitle(f'one 250-px cell of the 9-arcsec grid (P = 1024): metacal {TYPE} image with the drawImage FFT W forced',
             fontsize=10)
path = os.path.join(H.script_dir, 'metacal_fft_ghosts.png')
fig.savefig(path, dpi=150)
print(f'plot saved in {path}')

## the inner difference fitted with shifted copies of the image
inner = slice(N // 2 - 60, N // 2 + 60)
shifts = [(0, 0)] + [(a, b) for a in (-128, 0, 128) for b in (-128, 0, 128) if (a, b) != (0, 0)]
for i_target in range(3):
    m = mc if i_target == 0 else metacal_images(i_target)
    diff = (m[384] - m[W_REF])[inner, inner].ravel()
    X = np.array([shifted(m[W_REF], a, b)[inner, inner].ravel() for a, b in shifts]).T
    def explained(cols):
        coef = np.linalg.lstsq(X[:, cols], diff, rcond=None)[0]
        return 1 - np.sum((diff - X[:, cols] @ coef)**2) / np.sum(diff**2)
    d768 = (m[768] - m[W_REF])[inner, inner]
    print(f'target {i_target}: W = 384 inner rms {diff.std() / H.SIGMA_NOISE:.1e} sigma, variance explained by the '
          f'image itself {explained([0]):.2f}, by copies shifted by +-128 px {explained(list(range(1, 9))):.2f}; '
          f'W = 768 inner rms {d768.std() / H.SIGMA_NOISE:.1e} sigma')
