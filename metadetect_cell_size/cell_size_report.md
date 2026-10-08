# Why 250-pixel metadetect cells bias shear in the grid simulations

*October 2026. This folder is `MultiBand_ImSim/metadetect_cell_size/`: scripts in [`scripts/`](scripts/), key figures in [`plots/`](plots/). The full log of tests is [`cell_size_test_log.md`](cell_size_test_log.md), and section numbers (§4.x) refer to it. The outputs of the original runs (summary tables, job logs, other figures) are in `MultiBand_ImSim/test_scripts/metadetect_cell_size/`.*

## Summary

- **Problem.** In the seeing-0.73″ grid simulations (galaxies on a regular 9″ grid), metadetect with 250-px cells gives a multiplicative bias in g2 of m2 = −0.004 to −0.014, depending on the central size. With 500-px cells, |m| ≤ 0.0025.
- **Cause.** It is a numerical artefact of how metacal renders a cell with GalSim. It needs both a particular cell size and a regular galaxy lattice.
  - **The two FFT sizes.** Metacal turns the whole cell into a GalSim `InterpolatedImage`, whose Fourier-space table is P = good_fft_size(4N) px for an N-px cell. It then draws the sheared, reconvolved image back with an FFT of W ≈ good_fft_size(N + 10) px.
  - **The bad cell sizes.** For N just below a power of two, P/W is not an integer: for 246–256 px, P = 1024 and W = 384.
  - **What goes wrong.** The interpolation in the table then leaves faint copies of the cell's content that the drawing FFT folds back into the cell, shifted by P mod W (±128 px for 250-px cells).
- **Why only grid simulations.** On a lattice the copies of neighbours land at the same offset from every galaxy, so their effect on the shapes adds up coherently. The offset doesn't follow the true shear, but it does follow metacal's artificial shear, so the measured response no longer matches the true one. At realistic, random positions the copies land at random offsets and average out: m is consistent with zero on the same test bench.
- **Evidence.**
  - A noise-free test bench reproduces the bias (m2 = −0.014).
  - The bias appears for exactly the cell sizes where P/W is not an integer (246–256 and 502–512 px out of 22 sizes tested).
  - Forcing either FFT size inside metacal switches it on and off at fixed cell size.
  - The copies are seen directly in the metacal images.
- **Confirmation on the real images.** Re-running metadetect on the same images with 260-px cells removes the bias: m2 = +0.0014 ± 0.0007 instead of −0.0075 ([`m_cell_size.png`](plots/m_cell_size.png)). The central size stays 150, so these are exactly the same galaxies as the 250/150 baseline. A 260-px cell has P/W = 4, and the result matches the 500-px cells.
- **Recommendation.** For grid simulations, use a cell size outside the bad windows (not within ~10 px below 256 or 512; e.g. 240, or 260–480), or random galaxy positions. Alternatively, change metacal so that W divides P. For real data the artefact is not expected to matter.

## 1. The problem

The images are the collaborator's grid simulations (`/sdf/data/kipac/u/jlkitt/imsims/grid/seeing_073/`):
- **Sky:** 100 tiles of 0.5 deg², LSST r band, 0.2″ pixels, Moffat PSF with FWHM 0.73″.
- **Galaxies:** catalogue galaxies with 21 < r < 24.5, placed on a 9″ (45-px) grid.
- **Shear:** five shear tags (g1, g2 = ±0.02 or 0), with positions sheared too, and 0°/90° rotated pairs.

Metadetect 0.13.0 / ngmix 2.4.1 / GalSim 2.8.5 measure the shapes with wmom. R, m and c follow the collaborator's analysis (weights and cuts of `A_assign_weights`, fit over tiles).

Re-running metadetect on the same images with different `cell_size`/`central_size` (`analyse_cell_size.py`, figure **[`m_cell_size.png`](plots/m_cell_size.png)**) gives:

| cell / central (px) | m1, own R11 | m2, own R22 |
|---|---|---|
| 250 / 50 | −0.0009 | **−0.0140** |
| 250 / 100 | +0.0015 | **−0.0098** |
| 250 / 150 (collaborator's baseline) | +0.0016 | **−0.0075** |
| 250 / 200 | +0.0019 | **−0.0043** |
| 500 / 200, 250, 300, 400 | −0.0008 to +0.0014 | +0.0012 to +0.0025 |

Errors are ±0.0007.

- **The cell size matters, the central size doesn't.** At the same central size, 500-px cells are unbiased (§4.6).
- **Within 250-px cells, the bias comes from galaxies near the cell centre.** m2 ≈ −0.015 within ~25 px of the centre, falling to ~0 by ~100 px (figure **[`diag_position_radial.png`](plots/diag_position_radial.png)**, §4.3). Smaller central regions keep only those galaxies, hence the trend with central size. For a given central size, 250- and 500-px cells keep pixel-identical galaxies (§4.4), so this is not a selection effect.

## 2. What it is not

Each of these was tested and excluded (§5 of `cell_size_test_log.md`):
- noise, noise fixing, detection, trimming and selection (the bias is reproduced without any of them, §4.11)
- the hard cut of ImSim's 44-px grid stamps (§4.11)
- the PSF stamp size (§4.8)
- a frame or sign mismatch between position and shape shear (§4.5)
- the response convention alone (§4.1)
- the buffer size (§4.2)
- the galaxy density, and a subset of the galaxies (§4.10, §4.12)

## 3. Isolating it on a noise-free test bench

`test_hard_cut.py` and `test_layout.py` build a test bench that reproduces the effect quickly:
- **Scene:** synthetic, noise-free scenes drawn with ImSim's own code and the simulation's galaxies and PSF.
- **Shear and pairing:** a galaxy within 20 px of a cell centre, sheared by ±0.02 in g1 or g2, with a 90° rotated twin.
- **Measurement:** metacal exactly as metadetect calls it, then wmom at the known position.
- **Output:** m = S/R − 1, where S is the response to the true shear and R metacal's response.
- **Sample:** the same 317 galaxies in every comparison, so intrinsic shapes cancel exactly.

The bench gives m2 = −0.0140 ± 0.0009 for 250-px cells on the 9″ grid, the same as 250/50 in the real data, and m = +0.0003 for 500-px cells. With it:

- **The bias needs a regular lattice** (§4.12).
  - The same galaxies are unbiased on random positions at the grid's density (m2 = −0.0001 ± 0.0011), at the catalogue's density (−0.0010 ± 0.0007), and with no neighbours (+0.0003).
  - They are biased only on certain grid spacings, and the dependence on spacing is sharp and resonance-like (figures **[`m_vs_grid_spacing.png`](plots/m_vs_grid_spacing.png)** and **[`m_vs_grid_spacing_fine.png`](plots/m_vs_grid_spacing_fine.png)**, §4.13, §4.15). The 9″ grid sits on the flank of a feature that peaks at 8.7–8.8″, where m1 ≈ m2 ≈ −0.034.
- **The bias needs particular cell sizes** (figures **[`m_vs_cell_size.png`](plots/m_vs_cell_size.png)**, **[`m_map_cell_spacing.png`](plots/m_map_cell_spacing.png)** and **[`m_map_cell_spacing_500.png`](plots/m_map_cell_spacing_500.png)**, §4.16–4.19).
  - On the 9″ grid, cells of 240, 244, 260, 270, 300, 350, 400, 450, 498 and 500 px are clean to 1e-4.
  - Cells of 246–256 px and of 502–512 px are biased, increasingly so towards 256 and 512 px. 512-px cells are at least as bad as 256-px ones (m1 up to +0.12).
  - So it is not "small cells don't converge".
- **The true shear enters through the lattice** (figure **[`m_vs_jitter.png`](plots/m_vs_jitter.png)**, §4.14).
  - If the galaxies are sheared but the grid positions are not, the response to the true shear becomes exact, while metacal's R11 is 4% too high.
  - Jitter of ±0.5–2 px barely changes the bias; ±10 px nearly removes it.

## 4. The cause

### 4.1 Two FFT sizes inside metacal

For each N×N cell, ngmix metacal (as called by metadetect) does the following with GalSim (§4.7):

1. **Make the cell an `InterpolatedImage`** with GalSim's defaults: x-interpolant Lanczos-15, k-interpolant Quintic, `pad_factor = 4`. GalSim zero-pads the cell to **P = good_fft_size(4N)** pixels and stores its discrete Fourier transform as a table, to be interpolated in k.
2. **Divide by the PSF**, shear by ±0.01, and multiply by the round target PSF (in k-space).
3. **Draw the result back** with `drawImage(nx=N, ny=N, method='no_pixel')`.
   - GalSim samples k on a grid of spacing 2π/W, where **W = good_fft_size(max(2π/stepk, N))**.
   - 2π/stepk is ~N + 8–12 px, set by the cell's content.
   - good_fft_size rounds up to 2ⁿ or 3·2ⁿ⁻¹.

Usually W = P/4. But when N + ~10 px crosses a power of two while 4N does not, W jumps to 1.5× the next lower power of two and P/W becomes 8/3. The table below gives the sizes recorded inside real metacal calls (`check_fft_sizes.py`, `fft_sizes_all.csv`; 120 draws per cell size):

| cell N (px) | table P | drawing FFT W | draws with P/W = 8/3 | m2 at 9″ (bench) |
|---|---|---|---|---|
| 240, 244 | 1024 | 256 | 0 | +0.0004 |
| 246 | 1024 | 256 or 384 | 40% | +0.0028 (onset) |
| 248 / 250 / 256 | 1024 | 384 (mostly) | 77% / 93% / 100% | −0.006 / **−0.014** / −0.020 |
| 260–350 | 1536 | 384 | 0 | +0.0003 |
| 400–500 | 2048 | 512 | 0 | +0.0003 |
| 502–512 | 2048 | 512 or 768 | 27–93% | m1 up to +0.06 |

Over all 22 cell sizes, the clean ones draw every metacal image with P/W = 4, and the biased ones draw some of them with P/W = 8/3. The bias grows with that share. The figure **[`m_vs_table_margin.png`](plots/m_vs_table_margin.png)** shows the same switch against P − 4N: biased at ≤ 40 px, clean at ≥ 48 px.

### 4.2 Why a non-integer P/W corrupts the image: copies of the cell

- **Where the copies come from.** The table is a discrete Fourier transform, so it implicitly treats the zero-padded cell as repeating every P pixels. Evaluating it between table nodes with the Quintic interpolant does not fully suppress those repeats. In real space, the interpolated cell therefore carries faint copies of its own content P pixels away.
- **Where they land.** Drawing on a k-grid of spacing 2π/W makes the output periodic with period W, so a copy P pixels away lands **P mod W** pixels from where it started.
  - **W divides P:** the copies fold back exactly onto the cell. For the unsheared image every k-sample then falls on a table node, so the result is exact.
  - **P = 1024, W = 384:** they land ±128 px away (in x, in y and diagonally), inside a 250-px cell and on top of other galaxies.

Each galaxy thus receives a faint, shifted copy of whatever lies 128 px away from it. The error is small per pixel: 1–4 × 10⁻³ of the noise σ near the cell centre. But it is systematic.

### 4.3 Evidence for the mechanism

- **Forcing the FFT sizes switches the bias** (figure **[`m_vs_forced_fft.png`](plots/m_vs_forced_fft.png)**, §4.20). These runs change only metacal's GalSim sizes; the scenes are untouched.
  - **250-px cells become clean** when W is 256, 512 or 1024 (P/W = 4, 2, 1), or when the table is enlarged to 1536 px (P/W = 4).
  - **They stay biased** with W = 320, 384, 448 or 640, by up to 0.13.
  - **500-px cells, normally clean, become biased** when W is forced to 768 (P/W = 8/3), like the natural 508–512-px cells.
  - The bias is therefore set by the pair (P, W), not by the cell size.
- **Every biased variant has copies landing inside the cell; every clean one has none** (§4.21). That includes two non-integer ratios that are clean because their copies fall outside the cell: W = 768 (copies ±256 px away) and W = 2048 (1024 px away).
- **The copies are visible in the images** (figure **[`metacal_fft_ghosts.png`](plots/metacal_fft_ghosts.png)**, `check_fft_ghosts.py`, §4.21).
  - In the inner 120×120 px of a 250-px cell, the metacal image drawn with W = 384 differs from the one drawn with W = 512 by 1–4 × 10⁻³ σ.
  - **Copies of the image shifted by ±128 px explain 48–87% of that difference** (three scenes); the unshifted image explains 4%.
  - With W = 768, whose copies cannot reach the cell, the inner difference is 100–300× smaller.
- **It also explains the earlier convergence test** (§4.8): 250-px metacal images differed from a 1000-px reference by ~1e-3 σ, while 500-px images (P/W = 4) agreed to ~3e-5.

Two observations are not explained in detail:
- why the bias peaks for galaxies near the cell centre
- why it hits g2 rather than g1 at 9″ (at other spacings or cell sizes it is g1, or both)

Both depend on exactly where the copies land relative to the lattice, which was not modelled.

## 5. Why it is a problem for grid simulations only

The copies exist in every metacal image drawn with a non-integer P/W, whatever the scene. They bias the *mean* shear only if they affect all galaxies in the same way:

1. **On a lattice the contamination is coherent.** The copy offsets are fixed in the cell's frame (±128 px along x, y and the diagonals for 250-px cells). Three spacings of the 9″ grid are 135 px, so for every galaxy a copy of a lattice neighbour lands about 7 px away, in the same directions and at the same distance. Every galaxy's shape is nudged the same way.
2. **The contamination responds to shear differently from the galaxies.**
   - The true shear moves the lattice: the neighbour three spacings away shifts by ~3 px for g = 0.02. The copy offset, set by the FFT sizes, stays put, so the contamination changes with the true shear.
   - In metacal's artificially sheared images the copy offset is sheared too, and moves by ~P × 0.01 ≈ 10 px.
   - So the response metacal measures (R) does not match the response to the true shear (S), giving m = S/R − 1 ≠ 0.
   - This matches the test where the grid positions are not sheared: S becomes exact and only R is off ([`m_vs_jitter.png`](plots/m_vs_jitter.png)).
3. **The spacing dependence fits.** The bias is strongest at 8.7–8.8″, where three spacings are 130–132 px, just beyond the 128-px offset. It flips sign between 8.6″ and 8.5″ (129 → 127.5 px, i.e. as three spacings pass the offset), and m2 fades by 9.3″ (140 px) ([`m_vs_grid_spacing_fine.png`](plots/m_vs_grid_spacing_fine.png)). The 15″ and 18″ grids, where 128 px falls between lattice sites, are clean with 250-px cells.
4. **With realistic positions there is no coherence.** What lies 128 px from a galaxy is random: often nothing, sometimes another galaxy at a random relative position and orientation. The copies act like an extremely faint random blend whose effect on the mean shear averages to zero.
   - On the same bench, with the same 250-px cells and the same galaxies, random positions give m2 = −0.0001 ± 0.0011 at the grid's density and −0.0010 ± 0.0007 at the catalogue's density, against −0.014 on the 9″ grid (§4.12).
   - Disordering the grid by ±10 px largely removes the bias (§4.14).
5. **Per galaxy the error is tiny.** It is ~10⁻³ of the pixel noise, far below anything visible in a single measurement. Only on a lattice does it add up over millions of galaxies.

So the artefact is specific to regular grids, and to cell sizes with a non-integer P/W. Real data, and simulations with realistic positions, are not expected to suffer from it at the level tested (|m| ≲ 0.001).

## 6. Confirmation with the pipeline on the real images

**Prediction.** A 260-px cell has P = 1536 and W = 384 in every metacal draw (P/W = 4), so the copies fold back onto the cell. With a central size of 150 it keeps exactly the same galaxies as the 250/150 baseline (the central regions tile the image identically, §4.4). Only metacal's rendering of the cell changes, so the bias should disappear.

**Run.** `config_cell260_central150.ini`, through the same pipeline task (6_2) and seed as the baseline: all 5 shear tags × 100 tiles × 2 rotations, 1000 catalogues, 1.4–1.9 h per tag on 120 cores. It was analysed exactly like the other set-ups (`analyse_cell_size.py`, `results_cell_size.csv`), and the figure is **[`m_cell_size.png`](plots/m_cell_size.png)**, middle group.

| cell / central | m1, own R11 | m2, own R22 | R11 | R22 |
|---|---|---|---|---|
| 250 / 150 (baseline) | +0.0016 | **−0.0075** | 0.3371 | 0.3338 |
| **260 / 150** | **−0.0002** | **+0.0014** | 0.3349 | 0.3347 |
| 500 / 200–400 | −0.0008 to +0.0014 | +0.0012 to +0.0025 | 0.3342–0.3347 | 0.3345–0.3347 |

Errors are ±0.0007.

- **The bias is gone.** A 4% change in cell size moves m2 by +0.009 on the same galaxies and images, into agreement with the 500-px cells.
- **The R11/R22 asymmetry that 250-px cells show** (1% at 250/150; it had partly masked the effect in the averaged-R convention, §4.1) **disappears too.**
- The small positive m2 of +0.001 to +0.0025, common to the 260- and 500-px runs, is unrelated to this artefact. It is a separate open point (`cell_size_test_log.md` §6 D).

**The dip at the cell centre is gone too** (figure **[`diag_position_radial.png`](plots/diag_position_radial.png)**, green diamonds). m2 (own R22) by distance from the cell centre:

| distance from the centre (px) | 0–12.5 | 12.5–25 | 25–50 | 50–100 |
|---|---|---|---|---|
| 250 / 150 | −0.018 ± 0.004 | −0.021 ± 0.003 | −0.007 ± 0.001 | −0.005 ± 0.001 |
| 260 / 150 | −0.004 ± 0.004 | −0.005 ± 0.003 | +0.004 ± 0.001 | +0.001 ± 0.001 |
| 500 / 200 | +0.003 ± 0.006 | −0.005 ± 0.004 | −0.001 ± 0.002 | +0.004 ± 0.001 |

The 260-px profile follows the 500-px one within its errors, while every 250-px set-up shows the centre dip.

## 7. Recommendations

- **Grid simulations:** avoid cell sizes where GalSim ends up with P/W = 8/3, i.e. within ~10 px below 256 or 512 (246–256 and 502–512 px tested; presumably also below 1024).
  - 250 px is in the first window.
  - 500 px is clean, but sits only 2 px from the second.
  - Sizes such as 240 or 260–480 px have P/W = 4 for every draw.
  - Random galaxy positions also avoid the problem.
- **A fix in metacal (bench-tested):** make the drawing FFT divide the table, e.g. `pad_factor = 6` (P = 1536) or W = 512 for 250-px cells. ngmix exposes neither, so this needs a code change in ngmix's metacal or an upstream report to ngmix/GalSim.
- **Existing grid-simulation results with 250-px cells** carry this artefact: m2 ≈ −0.004 to −0.014 depending on the central size. They should not be read as a metadetect bias for real data.

## 8. Reproducing

| purpose | script | output |
|---|---|---|
| m per cell set-up on the real images | `submit_cell_size.sh`, `run_cell_size.sh`, `config_cell*_central*.ini`; `run_analyse.sh` (`analyse_cell_size.py`); `test_plot_m.py` | `results_cell_size.csv`, [`m_cell_size.png`](plots/m_cell_size.png) |
| m against position in the cell | `run_diag_position.sh` (`diag_position_in_cell.py`) | [`diag_position_radial.png`](plots/diag_position_radial.png) |
| noise-free test bench | `test_hard_cut.py` (bench; options `CELL_SIZES`, `MCAL_FFT_SIZE`, `MCAL_PAD_FACTOR`), `test_layout.py` (layouts) | `layout_response*.csv` |
| layout, spacing, jitter, cell-size and map tests | `run_layout.sh`, `run_layout_tests.sh <test>` | `m_vs_grid_spacing*.png`, [`m_vs_jitter.png`](plots/m_vs_jitter.png), [`m_vs_cell_size.png`](plots/m_vs_cell_size.png), `m_map_cell_spacing*.png` |
| GalSim sizes inside metacal | `check_fft_sizes.py` | `fft_sizes_all.csv`, [`m_vs_table_margin.png`](plots/m_vs_table_margin.png) |
| forced FFT sizes | `run_fft_size_test.sh` | `layout_response_fft_*.csv`, [`m_vs_forced_fft.png`](plots/m_vs_forced_fft.png) |
| copies in the images | `check_fft_ghosts.py` | [`metacal_fft_ghosts.png`](plots/metacal_fft_ghosts.png) |

All scripts are in [`scripts/`](scripts/), and all figures from the bench are made by `plot_layout_scans.py`. The scripts write their outputs next to themselves, and the sbatch wrappers write job logs to `scripts/logs/`. The outputs of the original runs are in `../test_scripts/metadetect_cell_size/`.
