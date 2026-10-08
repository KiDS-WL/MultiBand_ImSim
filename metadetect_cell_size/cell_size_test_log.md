# Metadetect cell-size test on grid simulations: test log (5–8 Oct 2026)

*This is the running log of the investigation. The scripts and configs have moved to [`scripts/`](scripts/) and the key figures to [`plots/`](plots/) of `MultiBand_ImSim/metadetect_cell_size/`. The write-up is [`cell_size_report.md`](cell_size_report.md).*

## TL;DR

> **Update 7 (2026-10-08, §4.22): confirmed on the real images. The write-up is `cell_size_report.md`.**
> - Metadetect with 260-px cells (central 150, the same galaxies as the 250/150 baseline; P/W = 4) gives m2 = +0.0014 ± 0.0007 instead of −0.0075, and m1 = −0.0002. That matches the 500-px cells, and the R11/R22 asymmetry is gone (`m_cell_size.png`).
>
> **Update 6 (2026-10-07, §4.19–4.21): mechanism found by forcing GalSim's FFT sizes inside metacal. Corrects Updates 4–5.**
> - **Correction.** The "drawImage FFT sizes" in 4.7, 4.16–4.18 and Updates 3–5 were mismeasured: they were half the cell's image table. The real FFT size W is good_fft_size(~N + 10): 256 px for 240–244-px cells, 384 for most draws of 248–256 and for 260–350, 512 for 400–500, and 768 for most draws of 504–512. The empirical rule of Update 5 still holds (as P − 4N ≤ 40 px), but its "FFT box edge" explanation is withdrawn.
> - **What sets the bias is the ratio P/W.** P = good_fft_size(4N) is the zero-padded table GalSim interpolates the cell from. Clean cells draw every metacal image with P/W = 4; biased cells (246–256, 502–512 px) draw 27–100% of them with P/W = 8/3, and the bias grows with that share (`m_vs_cell_size.png`, `m_map_cell_spacing*.png`).
> - **Forcing the sizes in metacal switches the bias on and off, with the scenes unchanged** (`m_vs_forced_fft.png`).
>   - 250-px cells become clean with W = 256, 512 or 1024, or with a 1536-px table.
>   - They stay biased with W = 320, 384, 448 or 640.
>   - 500-px cells become biased with W = 768.
> - **Mechanism** (`metacal_fft_ghosts.png`). Interpolating the table in k leaves faint copies of the cell P px away, and the drawImage FFT folds them back by P mod W. They land on the cell itself when P/W is an integer; for 250-px cells they land ±128 px away, on other galaxies. Those copies explain 48–87% of the metacal-image error in the inner cell. On a lattice they hit every galaxy the same way; at random positions they don't add up.
> - **Practical:**
>   - Avoid cell sizes within ~10 px below 256, 512, … (250 px is in that window).
>   - Alternatively, fix metacal by making W divide P (bench-tested: `pad_factor` 6, or W = 512, for 250-px cells).
>   - Real, non-lattice data is not expected to be affected.
>
> **Update 5 (2026-10-07, §4.18): the bias is set by the cell size relative to GalSim's FFT size, confirmed at a second FFT size** (now `m_vs_table_margin.png`).
> *Partly wrong, see Update 6: the W here is half the image table P, not the drawImage FFT, and the 'Origin' bullet is withdrawn. The results and the empirical rule (as P − 4N ≤ 40 px) stand.*
> - **The rule.** Metacal draws each N×N cell with a GalSim FFT of size W = good_fft_size(2N): 512 px for N ≤ 256, 768 for 257–384, 1024 for 385–512.
>   - Every cell with W − 2N ≥ 24 px is clean at every spacing tested: 10 sizes, including 240, 244, 498 and 500.
>   - Every cell with W − 2N ≤ 20 px is biased at some spacings: 246–256 and 502–512.
> - **At W = 1024:** 498 and 500 px are clean, while 502 → 512 px are increasingly biased. 512-px cells are at least as bad as 256-px ones (m1 = +0.062 at 9″), and even an 18″ grid is biased with them.
> - **Origin.** The cell occupies pixel coordinates 1…N inside GalSim's periodic FFT box [−W/2, W/2), leaving (W − 2N)/2 px before the wrap boundary. The bias appears when that gap is ≤ 10 px. The mechanism inside metacal is not yet understood. The padded image table has the same margin (doubled), so it is not yet excluded.
> - **Practical:**
>   - 250 px lies in a bad window (246–256).
>   - 500 px is clean but only 2 px from the next bad window (502–512).
>   - On lattices, avoid cell sizes within ~10 px below 256, 384, 512, 768, …. Random positions are clean even at 250 px (4.12).
>
> **Update 4 (2026-10-07, §4.17): map of m against cell size (244–256 px) and grid spacing (43–47 px)** (`m_map_cell_spacing.png`).
> *Partly wrong, see Update 6: the 'FFT sizes' quoted here are half the image table, not the drawImage FFT; GalSim's sizes are not the same across 244–256 px.*
> - **The bias switches on between 244 px and 246 px and grows steadily up to 256 px.** 244 px is clean at every spacing. At 256 px, m2 = −0.048 and m1 = −0.050 at 43.5–44 px. With 4.16 (260 px and up clean), the biased cells are 246–256 px, not 250 alone. GalSim's sizes are the same across 244–256 px.
> - **Empirical rule over all 15 cell sizes tested.** Call W the size of GalSim's drawImage FFT in metacal (512 px up to 256, 768 for 260–350, 1024 for 400–500) and N the cell size. Then:
>   - biased ⇔ W − 2N ≤ 20 px, growing as W − 2N → 0
>   - every clean cell has W − 2N ≥ 24 px
>   - it is an absolute margin, not a ratio: 500 px has the same W/2N as 250 px but is clean
> - **The spacing feature stays at 43.5–44 px for every biased cell,** so it isn't a cell/spacing commensurability.
> - **Decisive next test:** cells of 504–512 px (W − 2N = 16–0 against W = 1024) should be biased if the rule holds.
>
> **Update 3 (2026-10-07, §4.14–4.16): three follow-up tests on the noise-free test bench.**
> - **Cell size (4.16, `m_vs_cell_size.png`):** on the 9″ grid the bias appears **only at 250 px** (refined by 4.17: 246–256 px) among cells of 240, 250, 260, 270, 300, 350, 400, 450 and 500 px. 240 and 260 are clean to 1e-4. *(Corrected: 240 and 250 px share the table, 1024, but not the drawImage FFT, 256 vs 384; see Update 6.)* "Small cells don't converge" is therefore not the whole story: the bias needs a particular combination of cell size and grid spacing.
> - **Fine spacing scan (4.15, `m_vs_grid_spacing_fine.png`):** the 9″ grid sits on the flank of a feature that peaks at 8.7–8.8″ (43.5–44 px), where m1 ≈ m2 ≈ −0.034, 2.4× the 9″ value. m2 recovers to ~0 by 9.3″, and flips sign between 8.5″ and 8.6″. A 1% change in spacing changes m by ~0.01.
> - **Jitter and unsheared positions (4.14, `m_vs_jitter.png`):**
>   - Jitter of ±0.5–2 px barely changes the bias. ±5 px removes m2 but doubles m1 (−0.013). Only ±10 px is close to clean (+0.002 / +0.003 ± 0.001).
>   - With unsheared grid positions, the response to the true shear becomes exact, but metacal's R11 is 4% too high (m1 = −0.038).
>   - So the error depends steeply on the lattice geometry, and the true shear feeds it into the measurement by distorting the lattice.
> - **Practical:** for grid simulations, use random positions or 500-px cells. 500 px is clean at all 22 spacings tested (6–18″). Don't rely on small jitter, or on another small cell size that has only been tested at 9″.
>
> **Update 2 (2026-10-07, §4.13): spacing scan.** With 250-px cells, *any* regular grid can be biased, and the size and sign change sharply and non-monotonically with the spacing. m1 reaches +0.05 at 7.5″ and 11″; m2 is +0.013 at 6/8/8.5/12″ and −0.014/−0.015 at 9/13.5″. A few spacings are clean (7, 9.5, 10, 15″), so the 15″ grid being clean is luck, not "sparse enough". 500-px cells are unbiased at every spacing, and random positions are unbiased (§4.12). **Use 500-px cells or random positions for grid simulations** (small jitter doesn't help, §4.14).
>
> **Update (2026-10-07, §4.12):** the bias needs **both** 250-px cells **and** the regular 9″ (45-px) galaxy grid. On the noise-free test bench, the same galaxies show no bias (|m| ≲ 0.001) on 15″ or 18″ grids, at random positions with the same density as the 9″ grid, at random positions with the catalogue's density, or with no neighbours. In every layout, 500-px cells show no bias. So this is an artefact of the 9″ grid simulations with 250-px cells, not a problem expected for randomly placed galaxies. The bullets below predate this update; read "the cell size is the cause" as "the cell size, given the 9″ grid".

- **The large negative m2 comes from the metadetect cell size, not the central size, buffer size or galaxy layout.**
  - With 250-px cells: m2 = −0.004 to −0.014, correcting each component by its own response.
  - With 500-px cells at the **same** central sizes: m2 = +0.001 to +0.004, and m1 ≈ 0.
- **At cell 250 the bias comes from galaxies within ~25–50 px of the cell centre.** m2 there is −0.013 to −0.021, falling to ~0 by ~100 px. The trend with central size is just this profile averaged over the kept c×c region. No such profile exists at cell 500.
- **The root cause is localised to metacal's PSF deconvolution on a finite cell.**
  - Metacal images of 250-px cells are not converged: they differ from a 1000-px reference by 1–4×10⁻³ of the noise per pixel, most near the centre.
  - 500-px cells are converged (~3×10⁻⁵).
  - Reconvolution alone is exact (~10⁻¹⁰). The PSF stamp size (32–96 px) makes no difference.
- **Still open:**
  - why the deconvolution error peaks at the cell centre
  - whether its size fully accounts for m2 ≈ −0.015
  - the smallest cell size that converges
- **Practical:** use cells of ≥ 500 px until the convergence scan gives the minimum.

---

## 1. Setup

**Images.** The collaborator's grid simulations, read-only:
`/sdf/data/kipac/u/jlkitt/imsims/grid/seeing_073/<tag>/`

- **Config:** `/sdf/home/j/jlkitt/software/MultiBand_ImSim/test_scripts/grid/config_meta_073.ini`
- **Survey:** `simple_0.5sqdeg`, LSST_r, 0.2″/px, Moffat PSF with FWHM 0.73″ and β = 2.224, truncated at 4.5 × FWHM (16.4 px).
- **Size:** 100 tiles, galaxy rotations 0° and 90°, `shear_positions = True`.
- **Grid:** 9″ (45 px) grid, galaxies with 21 < mag < 24.5.
  - ImSim takes the galaxies of a 0.5 deg² patch of the input catalogue and packs them onto a square 9″ grid starting at RA = Dec = 0. That gives 31,644–43,303 galaxies per tile, 178–209 per side, 0.45–0.52° across.
  - Image sizes therefore vary per tile, from 8596 to 9217 px.
  - The sheared tags use a larger canvas so the sheared grid fits (8937 px vs 8596 px for tile 0). The corners of that canvas are empty.
- **Shear tags:**

  | tag | (g1, g2) |
  |---|---|
  | m020m020 | (−0.02, −0.02) |
  | m020p020 | (−0.02, +0.02) |
  | p000p000 | (0, 0) |
  | p020m020 | (+0.02, −0.02) |
  | p020p020 | (+0.02, +0.02) |

- **Baseline:** cell 250 / central 150, the collaborator's task 6_2 (`--rng_seed 940120`, `psf_image = centred`, `trim_position = noshear`). Catalogues are in `.../seeing_073/<tag>/catalogues/shapes`.

**Our outputs.** `/sdf/data/kipac/u/liss/ImSim/output/test_dev/jamar_seeing_073/<tag>/`

- `images` is a **symlink** to the collaborator's `images`.
- Catalogues go to `catalogues/shapes_cell<C>_central<c>`.
- `catalogues/input` is not linked, so `cross_match = False` in our configs. The m analysis doesn't need it.

**Code versions in the jobs.** metadetect 0.13.0 and sxdes 0.3.1 (vendored in `modules/`), ngmix 2.4.1, GalSim 2.8.5, sep 1.4.1.

**Noise file.** A byte-identical copy of the collaborator's is at `MultiBand_ImSim/noise_info/Euclid_Q1_median_LSST10yr_repeated100.csv`. It fixes the 100 tile labels.

**How m is computed** (`analyse_cell_size.py`)

- **Weights and cuts:** those of `utils_metadetect/A_assign_weights.py`: wmom, S/N > 12.5, T_ratio > 1.2, flags == 0, `shear_weight`. They're computed in memory, because the baseline catalogues are read-only.
- **Fit:** as `B_calculate_R_and_bias_whole.py`: R per tag, mean e per tile with both rotations summed, then a straight-line fit of measured against input shear over the tiles.
- **Two conventions:**
  - `m1`, `m2` use R = (R11 + R22)/2 for both components, as B does.
  - `m1_R11`, `m2_R22` use each component's own response.
- **The zero-shear tag** has no influence on the slope, because the shears are symmetric about 0. Fits with only the four sheared tags give identical m.

## 2. Code and scripts added

**Pipeline feature, in `MultiBand_ImSim/modules`.**

- **New option:** `[metadetect] PSF_image_file` (`RunConfigFile.py`, `Run.py`).
  - It's a path relative to each `psf_tile*_band*/` folder, e.g. `metadetect/psf_ima.fits`, and the file is used as it is.
  - `psf_image` still decides whether an off-centre stamp is allowed.
  - When the option is absent or `None`, behaviour is unchanged. It was committed as v1.3.0.

**Test scripts**

| file | purpose |
|---|---|
| `cell_size_report.md` | The final write-up: the cause, the evidence (key figures) and why only grid simulations are affected. |
| `config_cell<C>_central<c>.ini` | 10 set-ups: cell 250 with central 50/100/200; cell 500 with central 50/100/200/250/300/400; cell 260 with central 150 (4.22). Each differs only in `cell_size`/`central_size`; the output folder names follow from those values via interpolation. |
| `run_cell_size.sh` | sbatch job: one set-up × one tag, task 6_2 (torino, 120 CPUs, 4 GB/CPU). Finished tiles are skipped, so a crashed job can be resubmitted. |
| `submit_cell_size.sh` | Submits named set-ups, 5 tags each, with walltimes from measured runtimes. With no arguments it lists the set-ups. |
| `analyse_cell_size.py` + `run_analyse.sh` | R, m and c per set-up, written to `results_cell_size.csv` (both m conventions, R11, R22). |
| `test_plot_m.py` | Draws `m_cell_size.png` (m per set-up, both conventions, incomplete runs marked *), `m_buffer_size.png` and `m_N_per_central.png`. |
| `diag_position_in_cell.py` + `run_diag_position.sh` | m and R in bins of distance from the cell centre (`diag_position_in_cell.csv`), and per-tag × per-type maps over the cell (`diag_position_in_cell_maps.npz`). Plots `diag_position_radial.png` and `diag_position_maps_zero_shear.png`. `--plot-only` replots without reading the catalogues. |
| `plot_cell_galaxy_positions.py` | Rebuilds galaxy pixel positions from the ImSim layout rules, cuts cells as `MetaDetectShear` does, and stacks positions within the cell. Writes `cell_galaxy_positions.png` (set-up × tag) and `_zoom.png`. It is checked against the ImSim input catalogues (4e-11 px) and the image sizes (500 of 500 match). |
| `plot_cell_galaxy_positions_pairs.py` | 500/c next to 250/c on one pixel scale (`cell_galaxy_positions_pairs.png`), with a pixel-by-pixel check inside the central region. |
| `test_metacal_convergence.py` | Metacal images of 250/500 cells against a 1000-px reference, split into stages and repeated for PSF stamps of 16/32/48/96 px. Runs locally in under a minute. |
| `test_hard_cut.py` + `run_hard_cut.sh` | Noise-free test bench drawn with ImSim's code. Part a: metacal convergence for 44/74/256-px stamps. Part b: m of galaxies near the cell centre for each stamp and cell size (`hard_cut_response.csv`; use the "fixed" rows). Optional GalSim changes inside metacal for the FFT tests: `MCAL_FFT_SIZE` fixes the drawImage FFT, `MCAL_PAD_FACTOR` the table (only inside `metacal()`, only for images ≥ 200 px). |
| `test_layout.py` + `run_layout.sh` | The same test bench for different neighbour layouts: 9/15/18″ grids, random positions at two densities, no neighbours (`layout_response.csv`). Layouts (`LAYOUTS`: `grid<s>`, `jitter9_<a>`, `grid9_unsheared`, …), cell sizes (`CELL_SIZES`) and the output tag (`LAYOUT_TAG`) come from the environment. |
| `run_layout_tests.sh <test>` | sbatch wrapper for the follow-up layout tests: `spacing_scan`, `fine_scan`, `jitter`, `cell_scan`, `map`, `fft_margin`, `map500` (§4.13–4.19). `REF_CELL` sets the cell the others are compared with. Each writes `layout_response_<test>.npz/.csv` in about 5–15 minutes on 120 cores. `./run_layout_tests.sh <test> --summary-only` re-summarises. |
| `check_fft_sizes.py` | GalSim sizes inside metacal, recorded during real metacal calls: the image table P, the stepk radius, and the drawImage FFT W of every metacal image, per cell size (`CELL_SIZES`, `FFT_TAG`). Writes `fft_sizes_all.csv` (22 cell sizes). The earlier `fft_sizes_{cell_scan,map,fft_margin}.csv` held wrong W values and were deleted. About 10 s per cell size locally. |
| `check_fft_ghosts.py` | Image-level check of the mechanism (4.21): metacal images of a 250-px cell with W forced to 384/512/768, their differences, cross-correlations, and a fit with ±128-px shifted copies. Writes `metacal_fft_ghosts.png`. About 3 minutes locally. |
| `run_fft_size_test.sh` | sbatch job looping `test_layout.py` over variants with GalSim sizes inside metacal forced (`MCAL_FFT_SIZE`, `MCAL_PAD_FACTOR`; 4.20). Writes `layout_response_fft_<variant>.npz/.csv`. ~40 minutes on 120 cores. |
| `plot_layout_scans.py` | Figures from the layout tests: `m_vs_grid_spacing.png`, `m_vs_grid_spacing_fine.png`, `m_vs_jitter.png`, `m_vs_cell_size.png`, `m_map_cell_spacing.png`, `m_map_cell_spacing_500.png`, `m_vs_table_margin.png` (replaces `m_vs_fft_margin.png`), `m_vs_forced_fft.png`. Plain matplotlib, so it also runs on torino. |
| `test_grid_separation.py` + `run_grid_separation.sh` | m for sub-grids keeping every second galaxy in both directions (18″ apart), selected by matching detections to input galaxies and their grid position. Writes `results_grid_separation.csv` and `m_grid_separation.png`. The job's plot fails on torino (LaTeX fonts missing there); replot locally with `--plot-only`. |

**Runtime per shear tag at 120 cores.**

| set-up | hours per tag |
|---|---|
| 250/50 | 8.1–9.8 |
| 250/100 | 2.3 |
| 250/200 | 0.95 |
| 500/250 | 1.4 |
| 500/300 | 1.1 |
| 500/400 | 0.85 |

That's about 0.24 h + 2.6e-4 h per cell at cell 250, and about 0.49 h + 6.8e-4 h per cell at cell 500. Expected: 500/200 ~1.9 h, 500/100 ~6 h, 500/50 ~22 h.

## 3. Results (`results_cell_size.csv`)

| set-up | m1 (R avg) | m2 (R avg) | m1_R11 | m2_R22 | R11 | R22 |
|---|---|---|---|---|---|---|
| 250/50 | +0.0041 | −0.0189 | −0.0009 | **−0.0140** | 0.3385 | 0.3331 |
| 250/100 | +0.0051 | −0.0133 | +0.0015 | **−0.0098** | 0.3375 | 0.3333 |
| 250/150 (baseline) | +0.0036 | −0.0095 | +0.0016 | **−0.0075** | 0.3371 | 0.3338 |
| 250/200 | +0.0034 | −0.0058 | +0.0019 | **−0.0043** | 0.3361 | 0.3338 |
| 500/50 * | −0.0011 | +0.0051 | −0.0003 | +0.0043 ± 0.0018 | 0.3344 | 0.3350 |
| 500/100 * | −0.0007 | +0.0018 | −0.0009 | +0.0020 ± 0.0009 | 0.3346 | 0.3346 |
| 500/200 | +0.0005 | +0.0025 | +0.0008 | +0.0023 | 0.3345 | 0.3345 |
| 500/250 | +0.0007 | +0.0027 | +0.0009 | +0.0025 | 0.3344 | 0.3345 |
| 500/300 | −0.0006 | +0.0014 | −0.0008 | +0.0016 | 0.3347 | 0.3347 |
| 500/400 | +0.0003 | +0.0023 | +0.0014 | +0.0012 | 0.3342 | 0.3346 |

Errors are ±0.0007 unless given.

\* Incomplete. 500/50 has 187 of 1000 catalogues (~37 per tag); 500/100 has 698 (~139 per tag). All their jobs crashed at 21:20 on 2026-10-06 with "No space left on device" because `/sdf/data/kipac` was full. As of 2026-10-07 it has space again (60% used). Resume with:

```bash
./submit_cell_size.sh cell500_central100 cell500_central50   # finished tiles are skipped
sbatch run_analyse.sh; sbatch run_diag_position.sh; python test_plot_m.py
```

## 4. What was tested, in order

**4.1 R11 vs R22 (response convention)**

- At cell 250, R11 is 0.3–1.0% larger than R22. B's averaged R therefore pushes m1 up and m2 down.
- With each component's own R, m1 ≈ 0, but m2 stays strongly negative. So the convention explains only part of the m1/m2 split.
- At cell 250 the zero-shear tag's R11 is 2–3% higher than the sheared tags'. That doesn't affect m.

**4.2 m against buffer size and against detections per central region** (`m_buffer_size.png`, `m_N_per_central.png`)

- The buffer (cell − central) is not the driver. At the same buffer (100 or 200), cells 250 and 500 disagree strongly.
- Detections per central region = 10.74 × (c/150)². All set-ups measure the same images, so the density is the same; the count follows from it without re-running anything.

**4.3 Where in the cell the bias comes from** (`diag_position_in_cell.csv`)

m2_R22 in bins of distance from the cell centre, max(|dx|, |dy|):

| set-up | 0–12.5 px | 12.5–25 px | 25–50 px | 50–100 px |
|---|---|---|---|---|
| 250/50 | −0.0145 (15) | −0.0138 (9) | – | – |
| 250/100 | −0.0133 (32) | −0.0132 (19) | −0.0087 (9) | – |
| 250/150 | −0.0179 (43) | −0.0205 (26) | −0.0072 (14) | −0.0052 (10) |
| 250/200 | −0.0101 (61) | −0.0208 (36) | −0.0113 (17) | −0.0014 (9) |
| 500/50 * | +0.0069 (77) | +0.0035 (22) | – | – |
| 500/100 * | +0.0041 (41) | +0.0041 (23) | +0.0013 (11) | – |
| 500/200 | +0.0029 (61) | −0.0049 (36) | −0.0006 (17) | +0.0035 (9) |

Errors in brackets are in units of 10⁻⁴.

- **Cell 250:** the four central sizes share one radial profile: strongly negative within ~25 px of the centre, near zero by ~100 px. Averaging that profile over each set-up's kept region reproduces the overall m2 (for example, −0.0041 predicted vs −0.0043 measured for 250/200).
- **Cell 500:** there's no centre effect.
- **m1 and R11:** m1 also dips in the innermost 12.5 px at cell 250 (e.g. −0.0077 ± 0.0015 at 250/50). R11 rises towards the centre (0.3373–0.3382 vs 0.3346–0.3350 at 50–100 px), while R22 stays flat.
- **Trimming edges:** the bias is **not** concentrated at the edges of the central region. That is evidence against trimming or selection.
- **A separate effect:** the zero-shear maps show mean e2 alternating in sign between neighbouring grid-lattice sites (±0.006) at 250/150 and 500/300. So measured e2 depends on a galaxy's exact sub-lattice position. It appears at both cell sizes and averages out in the sheared tags, so it's not the cause.

**4.4 Galaxy layout within the cells** (`plot_cell_galaxy_positions*.py`)

- **The sample doesn't depend on the central size.** The central regions tile the image with step = central size from pixel 0, so every detection is kept exactly once. For the same central size, 500/c and 250/c keep **pixel-identical** galaxies inside the central region (checked: difference 0 out of 3.73 M).
- **Zero-shear tag:** galaxies sit on a lattice repeating every 5 px (gcd(45, c) = 5) or 15 px (c = 150, 300). That leaves 84% or 98% of the cell's pixels empty, and puts one lattice point ~0.5 px from the cell centre in every cell. This explains that tag's odd R11, but it carries no weight in m.
- **Sheared tags:**
  - Below 5 px they fill positions uniformly, at the Poisson level, and identically to each other.
  - At 5–30 px scales the **+g2** tags show diagonal moiré fringes of 5–9% RMS, while the −g2 tags are at ~1%. The cause is the sheared 45-px grid beating against the cell lattice.
- **Do the fringes cause the bias? No.**
  - Reweighting every sheared tag to the same in-cell distribution changes m2 by ≤1e-4.
  - Per tag, the layout shifts the measured g2 by ~1e-6, while ~2e-4 would be needed.
  - All four tags are equally low in |g2|: about 1% at 250/100, 1.4% at 250/50. The −g2 tags have flat layouts.
  - Fringe strength doesn't track m2: 250/50 has no fringes and the largest bias.

**4.5 Position shear vs shape shear**

- `ImSim._shear_positions` and the galaxy-profile shear give the same pixel-frame matrix for g1 and g2: no sign or frame mismatch.
- The up-to-86-px displacements at the tile edge are just the large-distance part of a uniform shear about the tile centre. They don't matter physically.
- Metacal shears each cell about its own centre. That has the same local distortion as the true shear and differs only by a translation of the whole cell.

**4.6 Cell size directly** (500/50, 500/100, 500/200)

- At the same central size, cell 500 gives m2 ≈ +0.002 to +0.004 against −0.004 to −0.014 at cell 250. That's 7–12σ, even with the partial runs.
- m2 at cell 500 doesn't depend on the central size, from 50 to 400.
- **Conclusion: the cell size is the cause.**

**4.7 Code-level trace of what metadetect does to a cell**

1. **Cut the cell** in `MetaDetect._run_metadetect_cell`: an N×N image, constant weight, the noise image at the same pixels, a Jacobian centred at N/2 − 0.5, and the 48×48 centred PSF stamp. The cell's random seed comes from its index, so 250/c and 500/c give the same central region the same seed.
2. **Fit the PSF.**
3. **Metacal** (`ngmix.metacal.get_all_metacal`, `psf='fitgauss'`, `fixnoise` with the noise image):
   - The noise copy is rotated 90° (only the array), put through metacal, rotated back and added. With round PSFs this equals adding the noise metacal'ed with the opposite shear, pixel by pixel.
   - For the data and the noise copy, the **whole cell array** becomes a GalSim `InterpolatedImage` (lanczos15, Quintic in k, zero-padded to `good_fft_size(4N)`). It is deconvolved by the PSF `InterpolatedImage`, sheared about the cell centre, convolved with the round target Gaussian, and redrawn by FFT with `drawImage(nx=N, ny=N, method='no_pixel')`.
4. **Detect** with sep on each whole metacal image (no background subtraction).
5. **Measure:** 32-px stamps around each detection, then wmom (1.2″ FWHM).
6. **Un-shear** positions about (N−1)/2.
7. **Trim** to the central region on the noshear-frame positions.

**Only step 3 depends on N** for a galaxy in a given central region. Measured on real cells, GalSim's settings all scale with N:

| GalSim setting | cell 250 | cell 500 |
|---|---|---|
| padded image table | 1024 | 2048 |
| radius used for stepk | 125 px | 250 px |
| `drawImage` FFT size (corrected in 4.20; first given as 512 / 1024) | 384 px (93% of draws; 256 otherwise) | 512 px |

**4.8 Metacal image convergence** (`test_metacal_convergence.py`)

RMS of metacal(cell) − metacal(1000) over the central 250×250 pixels, in noise σ (0.311), averaged over 6 points:

| PSF stamp | stage | cell | 0–12.5 | 12.5–25 | 25–50 | 50–100 | 115–125 (250-cell rim) |
|---|---|---|---|---|---|---|---|
| 48 | reconvolve only | 250 | 2e-10 | 2e-10 | 2e-10 | 2e-10 | 0.57 |
| 48 | noshear (deconvolve + reconvolve) | 250 | 1.7e-3 | 2.7e-3 | 2.2e-3 | ~1e-3 | 0.44 |
| 48 | 1p | 250 | **3.8e-3** | 1.7e-3 | 2.0e-3 | ~1e-3 | 1.2 |
| 48 | 2p | 250 | 1.8e-3 | **3.9e-3** | 1.3e-3 | ~0.9e-3 | 0.88 |
| 48 | all stages | 500 | ~3e-5 | ~3e-5 | ~3e-5 | 3–8e-5 | ≤4e-4 |
| 32 or 96 | — | — | same as 48 | | | | |
| 16 | noshear | 500 | 8e-4 (not converged) | | | | |

- **Reconvolution is exact; the error comes with the deconvolution.**
- **250-px images are not converged; 500-px images are.**
- **The pattern matches the bias.** 1p is most wrong within 12.5 px, where R11 is anomalous. 2p and noshear are most wrong at 12.5–25 px, where m2 is most negative.
- **The PSF stamp size isn't the cause.** 32, 48 and 96 px give identical results; 16 px cuts the PSF wings and makes everything worse.

**4.9 Range of the deconvolution** (one-off check, not saved)

- Metacal's deconvolve-then-reconvolve, applied to a point source, still responds at 1e-6 to 1e-5 of its peak 25–300 px away. Reconvolution alone falls to ~1e-10 by 25 px.
- So deconvolution is non-local, which is why the cell edge matters at all.
- But the tail is about as strong at 250 px as at 125 px. It doesn't, by itself, explain why 500 converges and 250 doesn't, or why the error peaks at the centre.

**4.10 Galaxy separation by sub-sampling** (`test_grid_separation.py`, 2026-10-07)

Prompted by a colleague's finding that a 15″ grid doesn't show the problem.

- **Method:**
  - Every detection, in all five metacal catalogues, is matched to its input galaxy (100% of weighted detections match within 2″).
  - Each galaxy's grid position (i, j) comes from the zero-shear input catalogue. A galaxy sits at the same grid point in every shear tag (checked to 4e-13″).
  - m is then measured on the four sub-grids (i mod 2, j mod 2). Each holds 25% of the galaxies, 18″ apart.
- **Check:** the "all" selection reproduces `results_cell_size.csv` to 1e-16.
- **Result: no change.** m2_R22 for the (even, even) sub-grid against all galaxies:

  | set-up | sub-grid (18″) | all (9″) |
  |---|---|---|
  | 250/50 | −0.0136 ± 0.0014 | −0.0140 |
  | 250/100 | −0.0100 | −0.0098 |
  | 250/150 | −0.0064 | −0.0075 |
  | 250/200 | −0.0038 | −0.0043 |
  | cell 500 | +0.001 to +0.005 | unchanged |

  All four sub-grids agree within their noise, and m1 is unchanged too.
- **Conclusion:** the bias isn't carried by a particular subset of the grid; every galaxy is biased equally.
- **Limitation:** this changes only the measured sample, not the images. Every selected galaxy still has its 9″ neighbours in detection, metacal and measurement. So it does **not** test whether physically sparser images behave differently, which is what the colleague's 15″ grid tests.
- **No 15″ grid images found.** The colleague's `bright` config has `grid_size = 15` but `position_type = true`, so the grid size is unused, and its output folder (`/sdf/data/kipac/u/jlkitt/imsims/bright/seeing_073`) no longer exists.

**4.11 Hard cut of the grid stamps** (`test_hard_cut.py`, `run_hard_cut.sh`, 2026-10-07)

- **Background:** in grid mode ImSim draws each PSF-convolved galaxy into a square stamp of 2·floor(grid/0.2/2) px and drops the light outside it. That's 44 px for the 9″ grid and 74 px for 15″.
  - At 44 px, 41% of galaxies lose more than 0.1% of their flux (10% lose more than 5%). At 74 px, 25% do.
  - The cut is applied *after* the PSF convolution, so the image is no longer exactly "sky ⊗ PSF". The Sérsic truncation (5 Re) and the PSF truncation (4.5 FWHM) happen *before* convolution, so they keep the image consistent.
- **Test:** synthetic, noise-free scenes drawn with ImSim's own code (`GalaxiesImage`, `MoffatPSF`, `_shear_positions`), using the simulation's galaxies on the 9″ grid. Only the stamp changes: 44 px, 74 px, or 256 px (uncut).
- **Part a, convergence:** the 250-px metacal error against a 1000-px reference is **identical to three digits for all three stamps** (1–6e-3 σ in the interior). 500-px cells converge to ~5e-6. The scenes do differ (the uncut ones carry 0.7% more flux, spread over 5% of the pixels), but the cut-off light is negligible in the error.
- **Part b, bias:**
  - **Set-up:** 480 galaxies placed within 20 px of a cell centre with their grid neighbours, true shear ±0.02 in g1 or g2, with 0°/90° rotation pairs. Metacal is run on the 250- and 500-px cells, and wmom is measured on a 32-px stamp at the known position.
  - **Sample:** the analysis's cuts, keeping only the 317 galaxies that pass in every measurement. Weighted or per-measurement cuts break the pairing and are too noisy (±0.04).
  - **Result:**

    | | cell 250 | cell 500 |
    |---|---|---|
    | m1 | −0.0068 ± 0.0011 | +0.0005 |
    | m2 | **−0.0141 ± 0.0009** | +0.0003 |
    | m2 at 250/50 in the real data | −0.0140 | |

  - **The stamp cut changes m by < 2e-4. The hard cut is not the cause.**
  - **Where the mismatch is:** at cell 250, S2 = 0.3588 against R22 = 0.3639, and S1 = 0.3665 against R11 = 0.3691. At cell 500, S = R to 1e-4.
  - **Both responses suffer in 250-px cells,** by different amounts. The response to the true shear in the noshear image drops by 1.9% (S2 0.3588 vs 0.3659), and metacal's R22 by only 0.5%.
- **What this isolates:** the bias is reproduced **without noise, noise fixing, detection, trimming, selection or the stamp cut**. So it lies purely in metacal's deconvolve/shear/reconvolve rendering of a finite 250-px cell, combined with the wmom measurement.
- **A fast test bench:** a full run takes ~4 minutes on 64 cores (`sbatch run_hard_cut.sh`, then `python test_hard_cut.py b --summary-only`).

**4.12 Layout of the neighbours** (`test_layout.py`, `run_layout.sh`, 2026-10-07)

- **Set-up:** the same test bench as 4.11 (uncut stamps), with the same 480 target galaxies and offsets in every layout. Only the neighbours change. Positions are sheared by ImSim and shifted so the target stays in place. `grid9` reproduces the scenes of 4.11 exactly.
- **Result**, for galaxies passing the cuts in every measurement of a layout (317–344), unweighted:

  | layout | m2, cell 250 | m2, cell 500 | m2 difference, 250 − 500 | m1 difference, 250 − 500 |
  |---|---|---|---|---|
  | grid 9″ (44 per arcmin²) | −0.0140 ± 0.0009 | +0.0003 | **−0.0144 ± 0.0009** | **−0.0071 ± 0.0011** |
  | grid 15″ (16 per arcmin²) | +0.0004 | +0.0003 | 0.0000 | 0.0000 |
  | grid 18″ (11 per arcmin²) | +0.0003 | +0.0003 | 0.0000 | 0.0001 |
  | random, 9″-grid density | −0.0001 ± 0.0011 | +0.0004 | −0.0004 ± 0.0011 | −0.0005 ± 0.0009 |
  | random, catalogue density (~20 per arcmin²) | −0.0010 ± 0.0007 | +0.0003 | −0.0014 ± 0.0007 | −0.0003 ± 0.0006 |
  | isolated | +0.0003 | +0.0003 | 0.0000 | 0.0000 |

- **Conclusion: the bias needs the regular 9″ (45-px) lattice inside a 250-px cell.** It isn't the number of galaxies (random positions at the same density are unbiased), and it isn't a property of a 250-px cell alone (isolated galaxies and the other layouts are unbiased). This agrees with the colleague's finding that a 15″ grid has no problem.
- **To watch:** random at the catalogue density gives −0.0014 ± 0.0007 (2σ).
- **Practical implication:** galaxies at random positions shouldn't suffer this bias with 250-px cells. For grid simulations, use spacing ≥ 15″, random positions, or 500-px cells.

**4.13 Grid-spacing scan** (`sbatch run_layout_tests.sh spacing_scan`, i.e. `test_layout.py` with `LAYOUTS=grid6,…,grid15`; outputs `layout_response_spacing_scan.csv`; figure `m_vs_grid_spacing.png` from `plot_layout_scans.py`, 2026-10-07)

- **Set-up:** the same test bench as 4.12, with 13 grid spacings and the same 480 targets (317 pass the cuts), on uncut stamps. The 9″ grid reproduces 4.12 exactly.
- **Result, cell 250:**

  | spacing (″ / px) | m1 | m2 |
  |---|---|---|
  | 6 / 30 | +0.0075 | +0.0124 |
  | 7 / 35 | +0.0027 | −0.0002 |
  | 7.5 / 37.5 | **+0.0466** | +0.0003 |
  | 8 / 40 | +0.0069 | +0.0119 |
  | 8.5 / 42.5 | +0.0250 | +0.0131 |
  | **9 / 45** | **−0.0067** | **−0.0140** |
  | 9.5 / 47.5 | +0.0030 | +0.0005 |
  | 10 / 50 | −0.0009 | +0.0004 |
  | 10.5 / 52.5 | +0.0069 | +0.0003 |
  | 11 / 55 | **+0.0559** | +0.0003 |
  | 12 / 60 | +0.0076 | +0.0131 |
  | 13.5 / 67.5 | −0.0074 | −0.0152 |
  | 15 / 75 | +0.0005 | +0.0004 |

  Errors are 0.0003–0.003.
- **Cell 500:** m is +0.0003 to +0.0006 at every spacing.
- **Reading:**
  - The behaviour looks like a resonance and is not monotonic in the spacing. Sampled spacings can't show how narrow the features are; 4.15 resolves the one around 9″.
  - Some spacings give identical results: 9″ and 13.5″ (45 and 67.5 px, ratio 3:2), and 6″, 8″ and 12″ (30, 40 and 60 px). The effect seems to depend on how the lattice's discrete spatial frequencies line up with something in metacal's processing of a 250-px cell (FFT wrap 512 px, padded table 1024 px).
  - Whole-pixel and half-pixel spacings both appear among the biased ones, so this is not simply pixel phase.
  - The mechanism is not identified.

**4.14 Jittered and unsheared 9″ grids** (`sbatch run_layout_tests.sh jitter`; outputs `layout_response_jitter.csv` and `m_vs_jitter.png`, 2026-10-07)

- **Set-up:** the same test bench and the same 480 targets as 4.12 (317 pass the cuts).
  - `jitter9_<a>`: every neighbour is moved by a uniform random offset of up to ±a px in x and y. The target isn't moved. The offsets are the same in all eight scenes of a target, so the ±g pairing is kept.
  - `grid9_unsheared`: the galaxies are sheared, but their positions are not, so the lattice stays square and axis-aligned in every scene.
  - `grid9` reproduces 4.12 exactly.
- **Result, cell 250.** Cell 500 gives m1 = +0.0003 to +0.0007 and m2 = +0.0002 to +0.0003 in every row.

  | layout | m1 | m2 | S1 | R11 | S2 | R22 |
  |---|---|---|---|---|---|---|
  | grid, no jitter | −0.0067 ± 0.0011 | −0.0140 ± 0.0009 | 0.3663 | 0.3688 | 0.3586 | 0.3637 |
  | jitter ±0.5 px | −0.0069 ± 0.0011 | −0.0139 ± 0.0010 | 0.3661 | 0.3687 | 0.3586 | 0.3637 |
  | jitter ±1 px | −0.0073 ± 0.0012 | −0.0134 ± 0.0010 | 0.3658 | 0.3685 | 0.3587 | 0.3636 |
  | jitter ±2 px | −0.0092 ± 0.0013 | −0.0125 ± 0.0010 | 0.3643 | 0.3677 | 0.3588 | 0.3634 |
  | jitter ±5 px | **−0.0130 ± 0.0019** | −0.0000 ± 0.0012 | 0.3594 | 0.3641 | 0.3622 | 0.3622 |
  | jitter ±10 px (±2″) | +0.0023 ± 0.0015 | +0.0028 ± 0.0011 | 0.3608 | 0.3600 | 0.3654 | 0.3644 |
  | grid, positions not sheared | **−0.0383 ± 0.0027** | +0.0035 ± 0.0004 | 0.3612 | **0.3756** | 0.3659 | 0.3646 |
  | random, same density (4.12) | −0.0002 ± 0.0009 | −0.0001 ± 0.0011 | 0.3196 | 0.3197 | 0.3313 | 0.3314 |
  | isolated (4.12), for reference | +0.0003 | +0.0003 | 0.3613 | 0.3611 | 0.3659 | 0.3658 |

- **Jitter:**
  - Up to ±2 px it barely matters: m2 goes from −0.0140 to −0.0125. So the bias doesn't hinge on exact or sub-pixel positions. It is a collective property of the lattice.
  - At ±5 px, m2 vanishes but m1 doubles (−0.013, 7σ). At ±10 px, m is +0.002 to +0.003 (1.5σ and 2.5σ), still not as clean as random positions.
  - The two components react differently, so this is not a smooth fading-out of a single effect. Note also that ±5 px leaves the lattice's fundamental (45-px) Fourier component at 92% of its strength, yet changes m completely. So the bias does not simply follow the strength of the 45-px periodicity.
  - **Small jitter is not a fix for grid simulations.**
- **Unsheared positions:**
  - At cell 250, S1 and S2 become equal to the isolated-galaxy values (0.3612 vs 0.3613; 0.3659 vs 0.3659). **The error in the response to the true shear (S2 = 0.3586 on the sheared grid) appears only because the true shear distorts the lattice.**
  - Metacal's R is still wrong, and differently: R11 is 4% too high (0.3756 vs 0.3611), R22 is 0.3% low.
  - So the 250-px error depends steeply on the lattice geometry. A ±0.02 shear of the lattice moves the galaxies at the cell edge by ~2.5 px, and R11 differs by ~2% between the unsheared and sheared lattices (0.3756 vs 0.3688).
  - **Consistent with the real data:** the zero-shear tag, whose lattice is also unsheared, has an R11 2–3% higher than the sheared tags at cell 250 (4.1).
  - In this scene the "right" answer is m = 0: neighbours 45 px away don't blend (cell 500 gives +0.0007 / +0.0002), so the deviation at cell 250 is the numerical artefact.

**4.15 Fine spacing scan around 9″** (`sbatch run_layout_tests.sh fine_scan`; outputs `layout_response_fine_scan.csv` and `m_vs_grid_spacing_fine.png`, 2026-10-07)

- **Set-up:** as 4.13, with spacings 8.5–9.5″ in 0.1″ (0.5-px) steps. 8.5, 9 and 9.5″ reproduce 4.13 exactly.
- **Result, cell 250.** Cell 500 gives m = +0.0003 to +0.0007 at every spacing.

  | spacing (″ / px) | m1 | m2 |
  |---|---|---|
  | 8.5 / 42.5 | +0.0250 ± 0.0023 | +0.0131 ± 0.0009 |
  | 8.6 / 43 | −0.0100 ± 0.0021 | −0.0116 ± 0.0012 |
  | 8.7 / 43.5 | **−0.0340 ± 0.0027** | **−0.0306 ± 0.0022** |
  | 8.8 / 44 | **−0.0354 ± 0.0026** | **−0.0338 ± 0.0023** |
  | 8.9 / 44.5 | −0.0212 ± 0.0017 | −0.0247 ± 0.0016 |
  | 9.0 / 45 | −0.0067 ± 0.0011 | −0.0140 ± 0.0010 |
  | 9.1 / 45.5 | −0.0017 ± 0.0008 | −0.0062 ± 0.0005 |
  | 9.2 / 46 | −0.0045 ± 0.0010 | −0.0021 ± 0.0002 |
  | 9.3 / 46.5 | −0.0048 ± 0.0011 | −0.0003 ± 0.0001 |
  | 9.4 / 47 | −0.0008 ± 0.0010 | +0.0004 ± 0.0001 |
  | 9.5 / 47.5 | +0.0030 ± 0.0009 | +0.0005 ± 0.0000 |

- **Reading:**
  - The simulation's 9″ grid sits on the high-spacing flank of a feature that peaks at 8.7–8.8″ (43.5–44 px). There m1 ≈ m2 ≈ −0.034, 2.4 times the 9″ value of m2.
  - On the high side, m2 changes smoothly and reaches ~0 by 9.3″ (46.5 px). The feature is about 1″ (5 px) wide, and the 0.1″ steps resolve it.
  - On the low side it is abrupt: m flips from −0.01 to +0.01–0.025 between 8.6″ and 8.5″ (a 0.5-px change of spacing).
  - m1 and m2 are similar near the peak but diverge on the flank: at 9″, m2 ≈ 2·m1; at 9.2–9.3″, m1 ≈ −0.005 while m2 ≈ 0.
  - A ±0.02 shear changes the local lattice spacing by about ±2% (±0.9 px along the principal axes). Between 44 and 46 px, m2 changes by 0.03. This is in line with 4.14: the steep dependence on lattice geometry is what feeds the true shear into the error.

**4.16 Cell-size scan for the 9″ grid** (`sbatch run_layout_tests.sh cell_scan`; outputs `layout_response_cell_scan.csv` and `m_vs_cell_size.png`; GalSim sizes from `check_fft_sizes.py`, now in `fft_sizes_all.csv`, 2026-10-07)

- **Set-up:** as 4.12, `grid9` only, with cells of 240, 250, 260, 270, 300, 350, 400, 450 and 500 px around the same centres. Cells 250 and 500 reproduce 4.12 exactly.
- **Result:**

  | cell (px) | m1 | m2 | drawImage FFT W (px; corrected in 4.20) | padded table P (px) | stepk radius (px) |
  |---|---|---|---|---|---|
  | 240 | +0.0004 | +0.0004 | 256 | 1024 | 118–121 |
  | **250** | **−0.0067 ± 0.0011** | **−0.0140 ± 0.0009** | 384 (93% of draws), 256 | 1024 | 122–126 |
  | 260 | +0.0003 | +0.0003 | 384 | 1536 | 129–131 |
  | 270 | +0.0003 | +0.0003 | 384 | 1536 | 131–136 |
  | 300 | +0.0003 | +0.0003 | 384 | 1536 | 147–151 |
  | 350 | +0.0004 | +0.0003 | 384 | 1536 | 172–176 |
  | 400 | +0.0003 | +0.0003 | 512 | 2048 | 198–201 |
  | 450 | +0.0004 | +0.0003 | 512 | 2048 | 219–226 |
  | 500 | +0.0004 | +0.0003 | 512 | 2048 | 243–251 |

  - All m errors are ≤ 1e-4 except at 250.
  - Every cell other than 250 agrees with 500 to ≤ 1e-4.
  - The GalSim sizes come from 6 targets × 4 shears per cell. The FFT size W was first given (wrongly) as 512/768/1024; it can also differ between the draws of one cell (4.20).
- **Reading:**
  - **On the 9″ grid only the 250-px cell is biased.** Cells 4% smaller or larger are clean to 1e-4.
  - ~~The FFT and table sizes don't set it on their own: 240 and 250 have the same drawImage FFT size.~~ *(Corrected in 4.20: 240 and 250 share the table P = 1024 but not the drawImage FFT, 256 vs 384, and P/W is what matters.)*
  - **No single ratio organises the results:**
    - cell/spacing: 260/45 = 5.78 is clean, while 250/43 = 5.81 gives m2 = −0.012 (4.15).
    - So the resonance involves the cell size, the lattice and at least one more scale.
  - This qualifies the convergence picture of 4.8, which compared only 250 and 500. Whether 240 and 260 converge at image level, or carry a similar image error that happens not to bias m, has not been checked.
  - **Practical:** 240 or 260 would remove the bias for the 9″ grid, but they haven't been tested at other spacings, and 4.13 shows the features are narrow and spread across spacings. 500 px is the only size tested clean across spacings.
- **For the real data:** the bench predicts m ≈ 0 for the 9″ grid images with 240- or 260-px cells. One pipeline run (e.g. 260/150 against the 250/150 baseline) would confirm that the bench and the pipeline agree on the cell-size dependence too.

**4.17 Map of m against cell size and grid spacing** (`sbatch run_layout_tests.sh map`; outputs `layout_response_map.csv` and `m_map_cell_spacing.png`; GalSim sizes now in `fft_sizes_all.csv`, 2026-10-07)

- **Set-up:** as 4.15, with cells of 244–256 px in 2-px steps, plus 500 px as the reference, against spacings of 43–47 px (8.6–9.4″) in 0.5-px steps. The 27 entries shared with 4.15 and 4.16 (cells 250 and 500) agree exactly. All cells of 244–256 px use the same padded table, P = 1024. *(Corrected in 4.20: they do **not** all use the same drawImage FFT: it is 256 for 244 px and 384 for a share of the draws growing from 40% at 246 px to 100% at 256 px.)*
- **m2, in units of 0.01.** Errors are 0.0001–0.0028; cell 500 is +0.0003 everywhere.

  | cell \ spacing (px) | 43 | 43.5 | 44 | 44.5 | 45 | 45.5 | 46 | 46.5 | 47 |
  |---|---|---|---|---|---|---|---|---|---|
  | 244 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 |
  | 246 | 0.0 | +0.3 | +0.4 | +0.3 | +0.3 | +0.2 | +0.1 | +0.1 | +0.1 |
  | 248 | −0.4 | −1.1 | −1.4 | −1.0 | −0.6 | −0.2 | −0.0 | +0.1 | +0.1 |
  | **250** | −1.2 | −3.1 | −3.4 | −2.5 | **−1.4** | −0.6 | −0.2 | −0.0 | 0.0 |
  | 252 | −1.5 | −3.6 | −3.9 | −2.9 | −1.6 | −0.7 | −0.2 | −0.0 | 0.0 |
  | 254 | −1.8 | −4.1 | −4.3 | −3.2 | −1.8 | −0.8 | −0.3 | −0.0 | 0.0 |
  | 256 | −2.0 | −4.5 | −4.8 | −3.5 | −2.0 | −0.9 | −0.3 | −0.1 | 0.0 |

- **m1, in units of 0.01.** Errors are 0.0000–0.0034; cell 500 is +0.0003 to +0.0007.

  | cell \ spacing (px) | 43 | 43.5 | 44 | 44.5 | 45 | 45.5 | 46 | 46.5 | 47 |
  |---|---|---|---|---|---|---|---|---|---|
  | 244 | +0.1 | +0.1 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | +0.1 |
  | 246 | −0.0 | −0.5 | −0.9 | −1.0 | −0.9 | −0.8 | −0.8 | −0.5 | −0.1 |
  | 248 | −0.4 | −2.2 | −2.7 | −2.1 | −1.2 | −0.8 | −0.9 | −0.6 | −0.1 |
  | **250** | −1.0 | −3.4 | −3.5 | −2.1 | **−0.7** | −0.2 | −0.4 | −0.5 | −0.1 |
  | 252 | −1.4 | −4.0 | −3.9 | −2.2 | −0.5 | −0.1 | −0.4 | −0.6 | −0.3 |
  | 254 | −1.7 | −4.5 | −4.3 | −2.3 | −0.5 | −0.0 | −0.5 | −0.8 | −0.4 |
  | 256 | −2.0 | −5.0 | −4.7 | −2.5 | −0.5 | 0.0 | −0.5 | −0.9 | −0.6 |

- **Reading:**
  - **The bias switches on between 244 and 246 px, at fixed GalSim sizes.** 244 px is clean at every spacing (|m| ≤ 6e-4, as clean as 500 px). At 246 px it starts: m1 reaches −0.010, and m2 is small and *positive* (up to +0.004). From 248 px on, both are negative.
  - **From 248 to 256 px it grows steadily.** m2 at 44 px goes −0.014, −0.034, −0.039, −0.044, −0.048. With 4.16, where 260 px and above are clean, the biased cells are 246–256 px.
  - **m2 keeps its shape in spacing.** For 250–256 px the curves are roughly the same shape, scaled up: m2(256)/m2(250) ≈ 1.4 at most spacings. m1 is less regular: at 45 px it is largest at 248 px (−0.012) and flat at about −0.005 for 250–256 px.
  - **The spacing feature does not move with the cell size.** m2 peaks at 44 px for every biased cell, and m1 at 43.5–44 px. If it were a commensurability of cell and spacing (fixed N/s), the peak would move by ~3% from 248 to 256 px, i.e. from 44 to ~45.4 px, which the 0.5-px steps would show. So the spacing dependence is set by something that doesn't change across the map, e.g. the fixed FFT grid (512 px).
- **An empirical rule, over all 15 cell sizes tested (4.16 + 4.17):** call W the drawImage FFT size and N the cell size.
  *(Corrected in 4.20: the 'W' in this table is half the image table, P/2, not the drawImage FFT. The rule itself holds, as P − 4N ≤ 40 px; the real FFT sizes and the explanation are in 4.20–4.21.)*

  | cell N (px) | 240 | 244 | 246 | 248 | 250 | 252 | 254 | 256 | 260 | 270 | 300 | 350 | 400 | 450 | 500 |
  |---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
  | FFT size W (px) | 512 | 512 | 512 | 512 | 512 | 512 | 512 | 512 | 768 | 768 | 768 | 768 | 1024 | 1024 | 1024 |
  | W − 2N (px) | 32 | 24 | **20** | **16** | **12** | **8** | **4** | **0** | 248 | 228 | 168 | 68 | 224 | 124 | 24 |
  | m2 at 45 px | 0.0004 | 0.0004 | +0.0028 | −0.0057 | −0.0140 | −0.0162 | −0.0180 | −0.0197 | 0.0003 | 0.0003 | 0.0003 | 0.0003 | 0.0003 | 0.0003 | 0.0003 |

  - The bias appears exactly when W − 2N ≤ 20 px, and grows as W − 2N → 0. Every clean cell has W − 2N ≥ 24 px. W is always half the padded table (good_fft_size(4N)), so the same holds for table − 4N ≤ 40 px.
  - It is an absolute margin in pixels, not a ratio. 500 px has the same W/2N as 250 px (1.024) and is clean, with the same 24-px margin as the clean 244-px cell.
  - **Caveat:** apart from 500 px, which sits right at the clean edge, all the biased and near-threshold cases are 244–256 px. So the rule isn't yet tested independently at another FFT size. *(Tested in 4.18 and 4.19: it holds.)*
  - The mechanism is still not identified. The rule points at how GalSim sizes the FFTs relative to the cell: a cell that nearly fills half the FFT (or a quarter of the padded table) is mis-rendered when it contains a lattice.
- **Decisive test:** cells of 504, 508 and 512 px (W = 1024, W − 2N = 16, 8, 0) on the 9″ grid. If the rule is right they are biased, with 512 px about as bad as 256 px. Done in 4.18: the rule holds.

**4.18 The FFT-margin rule at W = 1024** (`sbatch run_layout_tests.sh fft_margin`; outputs `layout_response_fft_margin.csv` and `m_vs_fft_margin.png`, now replaced by `m_vs_table_margin.png`; GalSim sizes now in `fft_sizes_all.csv`, 2026-10-07)

- **Set-up:** as 4.17.
  - Cells: 498, 500, 502, 504, 508 and 512 px, all with table P = 2048, so (P − 4N)/2 = 28, 24, 20, 16, 8, 0. Also 250 px. *(Corrected in 4.20: these were first described as 'W = 1024'. The real drawImage FFT is 512 for 498/500 px and 512 or 768 for 502–512 px.)*
  - Spacings: 43.5, 44 and 45 px, plus 88 and 90 px (17.6″, 18″). The last two have the same spacing-to-FFT ratio as 44/45 px at W = 512, in case the spacing feature scales with W.
  - 500 px is the reference (`REF_CELL=500`).
  - A 512-px cell needs a canvas 6 px wider on each side. On one target, I checked that this leaves the 250- and 500-px results bit-identical.
  - At 44 and 45 px one galaxy fewer passes the cuts in every measurement (316 instead of 317). So grid9/250 reads −0.0066 / −0.0140 here, against −0.0067 / −0.0140 before.
- **Result** (m1 / m2; cell 250 as in 4.15–4.17):

  | cell (W − 2N) | 43.5 px | 44 px | 45 px | 88 px | 90 px |
  |---|---|---|---|---|---|
  | 498 (28) | +0.0005 / +0.0004 | +0.0004 / +0.0003 | +0.0004 / +0.0003 | +0.0004 / +0.0004 | +0.0004 / +0.0003 |
  | 500 (24) | +0.0005 / +0.0004 | +0.0004 / +0.0003 | +0.0004 / +0.0003 | +0.0004 / +0.0004 | +0.0004 / +0.0003 |
  | 502 (20) | +0.0003 / +0.0003 | **+0.0061** / +0.0004 | **+0.0131** / +0.0003 | +0.0005 / +0.0003 | +0.0004 / +0.0004 |
  | 504 (16) | −0.0000 / +0.0002 | **+0.0152** / +0.0005 | **+0.0322** / +0.0003 | +0.0005 / +0.0003 | +0.0004 / +0.0004 |
  | 508 (8) | −0.0020 / **−0.0133** | **+0.0216** / **−0.0053** | **+0.0494** / −0.0000 | **+0.0073** / **−0.0093** | **+0.0060** / −0.0000 |
  | 512 (0) | **−0.0058** / **−0.0486** | **+0.0340** / **−0.0172** | **+0.0624** / −0.0006 | **+0.0185** / **−0.0198** | **+0.0117** / −0.0006 |
  | 250 (12, W = 512) | −0.0340 / −0.0306 | −0.0353 / −0.0337 | −0.0066 / −0.0140 | +0.0004 / +0.0003 | +0.0004 / +0.0003 |

  Bold: differs from cell 500 by more than 0.005 (all ≥ 5σ). Some tiny shifts, e.g. m2 = −0.0006 ± 0.0001 at 512/45 and 512/90 px, are also formally significant. Errors are ≤ 0.0001 for the clean entries and 0.0005–0.0031 for the biased ones.
- **Reading:**
  - **The rule holds at a second FFT size.** 498 and 500 px (W − 2N = 28, 24) are clean at all five spacings. From 502 px (20) the bias appears and grows steadily to 512 px (0). The switch-on is between 20 and 24 px, the same as at W = 512 (244 clean, 246 biased).
  - **Over all 21 cell sizes tested (4.16–4.18):**
    - every cell with W − 2N ≥ 24 px is clean at every spacing tested: 240, 244, 260, 270, 300, 350, 400, 450, 498, 500
    - every cell with W − 2N ≤ 20 px is biased at some of the spacings: 246–256 and 502–512
  - **Where it starts, and which component it hits, depends on the spacing and on W.**
    - At W = 1024 and 44–45 px, m1 is biased *positive* from 502 px on (up to +0.062 at 512/45 px), while m2 is barely affected at 45 px.
    - At 43.5 px it starts only at 508 px and is mainly m2 (−0.049 at 512).
    - At W = 512 the same spacings give mainly negative m1 and m2.
    - 88/90 px are biased at 508–512 px, but differently from 44/45 px at W = 512, so the spacing feature doesn't simply scale with W.
    - An 18″ grid (90 px) is therefore not safe with 508–512-px cells, although it is clean in 250-px cells.
  - **512-px cells are at least as biased as 256-px ones.** "Bigger cells converge" is the wrong picture: what matters is the cell size relative to GalSim's FFT size.
- **Where the margin comes from: withdrawn, see 4.20–4.21.** The reasoning below assumed ngmix's image keeps coordinates 1…N inside `drawFFT_makeKImage`. In fact `drawImage` first centres the bounds on 0, so W = good_fft_size(max(stepk-based size ≈ N + 10, N)), and the cell is centred in the FFT box. Kept for the record:
  - ngmix draws each metacal image into an N×N GalSim image with pixel coordinates 1…N, the profile centred on the image.
  - GalSim's FFT box is centred on coordinate 0 and must cover twice the largest coordinate: W = good_fft_size(max(stepk-based size, 2N)). In every case here 2N is the larger, so W = good_fft_size(2N).
  - The cell thus fills the positive quarter of the box and leaves (W − 2N)/2 px between its far edge and the box edge, where the periodic FFT wraps around.
  - **The bias appears when that gap is ≤ 10 px and is absent when it is ≥ 12 px.**
  - The InterpolatedImage's padded table (good_fft_size(4N) = 2W) has the same margin, doubled, so these tests can't tell which of the two matters.
  - How content within ~10 px of the wrap boundary biases galaxies near the cell centre is not understood. The wrapped part lands at negative coordinates, outside the drawn image. Lanczos15 ringing of the cell edge (15-px reach) and the reconvolution kernel are candidates for what crosses the boundary.
- **Predictions:** ~~bad windows at 374–384, 758–768 and 1014–1024 px~~, withdrawn. The P/W picture (4.20) predicts bad windows only within ~10 px below 256, 512, 1024, …; for N ≈ 375–384 it predicts P/W = 1536/512 = 3, an integer, so clean. Not tested.

**4.19 Map for cells of 500–512 px** (`sbatch run_layout_tests.sh map500`; outputs `layout_response_map500.csv` and `m_map_cell_spacing_500.png`, 2026-10-07)

- **Set-up:** as 4.17, with cells of 500–512 px in 2-px steps against spacings of 43–47 px, and 498 px as the reference (`REF_CELL=498`). The 18 entries shared with 4.18 agree exactly.
  - All these cells have the image table P = 2048.
  - The drawImage FFT W is 512 for 498/500 px. From 502 px a growing share of the draws uses 768 instead: 27, 57, 67, 76, 89 and 93% for 502–512 px (4.20).
- **m1, in units of 0.01.** Errors are up to 0.0078; 498 and 500 px are ≤ 0.0006 everywhere.

  | cell \ spacing (px) | 43 | 43.5 | 44 | 44.5 | 45 | 45.5 | 46 | 46.5 | 47 |
  |---|---|---|---|---|---|---|---|---|---|
  | 500 | +0.1 | +0.1 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 |
  | 502 | −0.0 | 0.0 | +0.6 | +1.2 | +1.3 | +0.4 | +1.1 | +3.1 | +2.7 |
  | 504 | −0.1 | −0.0 | +1.5 | +3.5 | +3.2 | +0.8 | +2.8 | +8.2 | +6.3 |
  | 506 | −0.2 | −0.0 | +1.8 | +4.6 | +4.1 | +0.8 | +3.3 | +9.9 | +7.9 |
  | 508 | +0.7 | −0.2 | +2.2 | +5.8 | +4.9 | +1.0 | +3.7 | +10.9 | +8.8 |
  | 510 | +1.2 | −0.6 | +2.9 | +7.3 | +6.0 | +1.2 | +4.1 | +11.6 | +9.4 |
  | 512 | +1.5 | −0.6 | +3.4 | +7.8 | +6.2 | +1.3 | +4.4 | +12.2 | +9.6 |

- **m2, in units of 0.01.** Errors are up to 0.0032. It is 0.0 everywhere for 500–506 px, apart from +0.1 at 506/44 px.

  | cell \ spacing (px) | 43 | 43.5 | 44 | 44.5 | 45 | 45.5–47 |
  |---|---|---|---|---|---|---|
  | 508 | −0.9 | −1.3 | −0.5 | −0.1 | −0.0 | 0.0 |
  | 510 | −3.3 | −4.3 | −1.5 | −0.3 | −0.0 | 0.0 |
  | 512 | −4.1 | −4.9 | −1.7 | −0.3 | −0.1 | 0.0 |

- **Reading:**
  - **500 px is clean everywhere.** From 502 px the bias grows with the share of draws that use W = 768. That is the same switch as in the 246–256-px window (4.17, 4.20).
  - **It is mostly m1 and positive,** up to +0.12 at 512/46.5 px. m2 is biased only at 43–44 px, from 508 px.
  - **The spacing pattern differs from the 250-px window.** m1 peaks at 44.5 and 46.5 px and dips at 45.5 px. 4.21 suggests why: the copies that cause the bias land ±256 px away here, against ±128 px for 250-px cells, so they meet the lattice differently.

**4.20 Real FFT sizes inside metacal, and forcing them** (`check_fft_sizes.py` → `fft_sizes_all.csv`; `run_fft_size_test.sh` → `layout_response_fft_<variant>.csv`; figures `m_vs_forced_fft.png`, `m_vs_cell_size.png`, `m_vs_table_margin.png`, 2026-10-07)

- **Correction.** The drawImage FFT sizes given in 4.7 and 4.16–4.18 were measured by calling `drawFFT_makeKImage` directly on an image with coordinates 1…N. Inside `drawImage`, GalSim first centres the bounds on 0. So the real size is W = good_fft_size(max(stepk-based size, N)), where the stepk-based size is ~N + 8–12 px, content-dependent.
  - The old numbers were half the image table, P/2.
  - The real sizes, now recorded during real metacal calls by hooking `drawFFT_makeKImage` and ngmix's `_galsim_stuff_impl`, 6 targets × 4 shears × 5 metacal types per cell:

  | cell N (px) | 240 | 244 | 246 | 248 | 250 | 252 | 254 | 256 | 260–350 | 400–500 | 502 | 504 | 506 | 508 | 510 | 512 |
  |---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
  | table P (px) | 1024 | 1024 | 1024 | 1024 | 1024 | 1024 | 1024 | 1024 | 1536 | 2048 | 2048 | 2048 | 2048 | 2048 | 2048 | 2048 |
  | drawImage W (px) | 256 | 256 | 256/384 | 256/384 | 256/384 | 256/384 | 256/384 | 384 | 384 | 512 | 512/768 | 512/768 | 512/768 | 512/768 | 512/768 | 512/768 |
  | draws with P/W = 8/3 | 0 | 0 | 40% | 77% | 93% | 96% | 99% | 100% | 0 | 0 | 27% | 57% | 67% | 76% | 89% | 93% |
  | biased? | no | no | onset | yes | yes | yes | yes | yes | no | no | yes | yes | yes | yes | yes | yes |

  - **Every clean cell draws all its metacal images with P/W = 4. Every biased cell draws a share of them with P/W = 8/3, and the bias grows with that share.**
  - Why it happens: P = good_fft_size(4N), while W = good_fft_size(N + ~10). For N within ~10 px below 256 (or 512), N + 10 exceeds 256 and W jumps to 384 while P stays at 1024. That is also why "P − 4N ≤ 40 px" (the old "W − 2N ≤ 20") matched the data.
- **Forced-FFT test (`run_fft_size_test.sh`).**
  - Set-up: the test bench at 43.5, 44 and 45 px, with GalSim's sizes inside metacal changed through `test_hard_cut.py`. `MCAL_FFT_SIZE` fixes W for every draw of the cell; `MCAL_PAD_FACTOR` changes the table.
  - Both changes act only inside `metacal()`, so the scenes are unchanged. An early version also changed ImSim's drawing of the 256-px galaxy stamps; that was caught and fixed.
  - Check: forcing GalSim's own choice reproduces the default bit for bit (250 px with W = 384 on target 1; 500 px with W = 512).
  - Results, m1 / m2:

  | variant | P | W | P/W | 43.5 px | 44 px | 45 px |
  |---|---|---|---|---|---|---|
  | 250, default | 1024 | 384 (93%) | 8/3 | −0.0340 / −0.0306 | −0.0354 / −0.0338 | −0.0067 / −0.0140 |
  | 250, W 256 | 1024 | 256 | 4 | +0.0005 / +0.0004 | +0.0004 / +0.0003 | +0.0004 / +0.0004 |
  | 250, W 320 | 1024 | 320 | 16/5 | −0.0021 / −0.0027 | −0.0008 / −0.0028 | +0.0026 / −0.0010 |
  | 250, W 384 | 1024 | 384 | 8/3 | −0.0384 / −0.0360 | −0.0368 / −0.0372 | −0.0049 / −0.0146 |
  | 250, W 448 | 1024 | 448 | 16/7 | +0.0075 / **+0.1208** | −0.0115 / **+0.1297** | **−0.0612** / **+0.0588** |
  | 250, W 512 | 1024 | 512 | 2 | +0.0006 / +0.0003 | +0.0005 / +0.0003 | +0.0005 / +0.0003 |
  | 250, W 640 | 1024 | 640 | 8/5 | −0.0014 / −0.0027 | −0.0010 / −0.0028 | +0.0014 / −0.0010 |
  | 250, W 768 | 1024 | 768 | 4/3 | +0.0004 / +0.0003 | +0.0003 / +0.0003 | +0.0003 / +0.0003 |
  | 250, W 1024 | 1024 | 1024 | 1 | +0.0006 / +0.0003 | +0.0005 / +0.0003 | +0.0005 / +0.0003 |
  | 250, W 2048 | 1024 | 2048 | 1/2 | +0.0005 / +0.0003 | +0.0004 / +0.0003 | +0.0004 / +0.0003 |
  | 250, pad 6 | 1536 | 384 | 4 | +0.0004 / +0.0003 | +0.0003 / +0.0003 | +0.0003 / +0.0003 |
  | 250, pad 8 | 2048 | 384 | 16/3 | +0.0013 / +0.0021 | +0.0009 / +0.0023 | −0.0005 / +0.0013 |
  | 250, pad 8 + W 512 | 2048 | 512 | 4 | +0.0004 / +0.0003 | +0.0003 / +0.0003 | +0.0003 / +0.0003 |
  | 250, pad 6 + W 768 | 1536 | 768 | 2 | +0.0004 / +0.0003 | +0.0003 / +0.0003 | +0.0003 / +0.0003 |
  | 500, default | 2048 | 512 | 4 | +0.0005 / +0.0004 | +0.0004 / +0.0003 | +0.0004 / +0.0003 |
  | 500, W 768 | 2048 | 768 | 8/3 | −0.0032 / **−0.0322** | **+0.0355** / −0.0104 | **+0.0655** / −0.0003 |
  | 500, W 1024 | 2048 | 1024 | 2 | +0.0006 / +0.0004 | +0.0005 / +0.0003 | +0.0005 / +0.0003 |

  Errors are ≤ 0.0002 for the entries at the +0.0003–0.0006 level and up to 0.009 for the large ones.

- **Reading:**
  - **Changing either size alone switches the bias.** With the table fixed, W = 256, 512 or 1024 makes 250-px cells clean. With W fixed at 384, a 1536-px table makes them clean. And forcing W = 768 makes the clean 500-px cell biased, like the natural 508–512-px cells.
  - **Every variant with P/W an integer is clean** (8 cases). So are W = 768 (P/W = 4/3) and W = 2048 (1/2). The other non-integer ratios are biased, by very different amounts: W = 448 up to 0.13, W = 384 up to 0.038, and W = 320, W = 640 and pad 8 around 0.001–0.003. 4.21 explains which.
  - Forcing W = 384 for every draw gives a slightly larger bias than the default, where 93% of draws use 384. That fits the bias growing with the share of 8/3 draws.
  - **The cell size itself, its distance from the FFT box edge (4.18), and the lattice-to-cell ratio are not what matters:** at fixed N = 250 the bias comes and goes with P and W.

**4.21 Mechanism: copies of the cell from the image table, folded into it by the drawImage FFT** (`check_fft_ghosts.py`, `metacal_fft_ghosts.png`, 2026-10-07)

- **The idea.** Metacal makes the cell an InterpolatedImage, zero-padded to a table of P px, and evaluates it in k-space with a Quintic k-interpolant.
  - In real space, interpolating the table in k leaves faint copies of the cell's content P px away.
  - Drawing with an FFT of W px samples k on a grid of spacing 2π/W, which folds real space with period W. So a copy lands P mod W from the cell.
  - When P/W is an integer, the copies fold back exactly onto the cell, and the image is restored.
  - Otherwise they land somewhere else, and they matter when that is inside the cell.
- **Copy shifts for the variants of 4.20** (a 250-px cell spans ±125 px; a 500-px cell ±250 px):

  | P, W | shift P mod W | inside the cell? | result |
  |---|---|---|---|
  | integer P/W | 0 (onto itself) | — | clean |
  | 1024, 384 | ±128 | yes | biased (0.03–0.04) |
  | 1024, 448 | ±128 | yes | biased (up to 0.13) |
  | 1024, 320 | ±64 | yes | biased (0.001–0.003) |
  | 1024, 640 | ±256 first-order; ±128 from copies 2P away | first-order no; second-order yes | biased (0.001–0.003) |
  | 1024, 768 | ±256 | no (cell ±125) | clean |
  | 1024, 2048 | 1024 | no | clean |
  | 2048, 384 (pad 8) | ±128 | yes, but the table is twice the size | biased (0.001–0.002) |
  | 2048, 768 (500-px cell) | ±256 | yes (cell ±250) | biased (up to 0.07) |

  Every biased case has copies landing inside the cell, and every clean one has none. How strong a copy is (first- vs second-order, table size) and where it falls relative to the lattice set the size of the bias. That last part is not modelled.
- **Image-level check** (`check_fft_ghosts.py`): metacal noshear images of 250-px cells (9″ grid, 3 scenes), drawn with W = 384, 512 and 768.
  - In the inner 120×120 px, W = 384 differs from W = 512 by 1.3–3.8e-3 σ. **Copies of the image shifted by ±128 px (in x, y and both) explain 48–87% of that difference**, while the unshifted image explains 4%.
  - With W = 768 the inner difference is 1–2e-5 σ, 100–300× smaller, because its copies (±256 px) cannot reach into the cell.
  - In the figure, the largest W = 384 difference near the centre lies within ~10 px of where a bright galaxy at the cell's right edge lands after a −128-px shift.
- **Why lattices, and why it biases shear (consistent with the data, not tested directly):**
  - On a lattice, a copy of a neighbour lands at the same offset from every galaxy, so the contamination is coherent. At random positions it isn't (4.12).
  - The copy offset P mod W stays fixed while the true shear moves the galaxies, which explains why the true-shear response is affected (4.11) and why unsheared positions fix it (4.14). In metacal's sheared images the offset is sheared as well, moving by ~P × 0.01 ≈ 10 px, which explains why R is so sensitive (4.14).
  - For 250-px cells the copies sit 128 px away. Three lattice spacings are 127.5–141 px for spacings of 42.5–47 px. The bias is largest at 43.5–44 px (3 spacings = 130.5–132 px), flips sign at 42.5 px (127.5 px, 4.13/4.15), and fades by 46.5–47 px (139.5–141 px). That fits copies of neighbours landing on or next to the measured galaxy.
  - Jitter of ±5–10 px (4.14) scatters where the copies land.
- **This also explains the earlier convergence test (4.8):** 250-px metacal images differ from a 1000-px reference by ~1e-3 σ near the centre, while 500-px images (P/W = 4) agree to ~3e-5.
- **Practical:**
  - The bad cell sizes are those where GalSim ends up with P/W = 8/3: within ~10 px below 256, 512 (tested), and presumably 1024. 250 px is in the first window.
  - All other sizes tested are clean, even on lattices.
  - On the bench, metacal for 250-px cells is fixed by making W divide P: `pad_factor` 6 (P = 1536) or W = 512. ngmix doesn't expose either, so this needs a code change. The alternative is a cell size outside the windows, e.g. 240 or 260–500.
  - For real (non-lattice) data the effect does not add up coherently (4.12).

**4.22 Pipeline confirmation: 260/150 on the real grid images** (`config_cell260_central150.ini`, `./submit_cell_size.sh cell260_central150`, `run_analyse.sh`; figure `m_cell_size.png`, 2026-10-08)

- **Why 260/150:**
  - A 260-px cell has P = 1536 and W = 384 for every draw (P/W = 4).
  - Central 150 keeps the same galaxies as the baseline 250/150 (4.4).
- **Run:** 1000 catalogues, 1.4–1.9 h per tag. The results of the existing set-ups are unchanged by the re-analysis.
- **Result:**

  | cell / central | m1 | m2 | m1, own R11 | m2, own R22 | R11 | R22 |
  |---|---|---|---|---|---|---|
  | 250 / 150 | +0.0036 | −0.0095 | +0.0016 | −0.0075 | 0.3371 | 0.3338 |
  | **260 / 150** | −0.0004 | +0.0016 | −0.0002 | **+0.0014** | 0.3349 | 0.3347 |
  | 500 / 200 | +0.0005 | +0.0025 | +0.0008 | +0.0023 | 0.3345 | 0.3345 |

  Errors are ±0.0007.
- **Position diagnostic** (`run_diag_position.sh`, `diag_position_radial.png`): m2 by distance from the cell centre (0–12.5, 12.5–25, 25–50, 50–100 px)
  - 250/150: −0.018, −0.021, −0.007, −0.005
  - 260/150: −0.004 ± 0.004, −0.005 ± 0.003, +0.004 ± 0.001, +0.001 ± 0.001
  - The centre dip is gone; 260/150 follows the 500-px profiles.
- **Conclusion:** the bias disappears on the real images when the cell size moves out of the bad window. That confirms the bench and the mechanism (4.20, 4.21).

## 5. Ruled out

- **Noise image correlated with the image noise.** They use different seeds (`ImSimSkySimple.py`: `noise` seed + 120 × rotation vs `noise_2` seed + 99). The 90° noise rotation creates no centre-specific effect.
- **Position-shear vs shape-shear frame or sign mismatch.** See 4.5.
- **Cell or PSF Jacobian off-centre.** Both are at the true centres.
- **Background subtraction pattern.** sxdes and metadetect don't subtract a background.
- **Count imbalance between metacal catalogues** (1p/1m, 2p/2m): ~1e-5.
- **Trimming as a shear-dependent selection**, the explanation behind the earlier `shear_positions` fix. The tiling partitions the image and the bias isn't at the edges.
- **Galaxy layout, the moiré fringes, and differences in how tags sample the cell.** See 4.4.
- **Galaxy density or the cell size on its own.** Random positions at the 9″-grid density, other grid spacings and isolated galaxies show no bias at cell 250 (4.12).
- **The hard cut of the grid stamps (44 px).** It makes no difference to convergence or to m (4.11).
- **Noise, noise fixing, detection, trimming and selection.** The bias is reproduced without any of them (4.11).
- **The bias being carried by a subset of grid galaxies.** Every-second-galaxy sub-grids give the same m (4.10). Physical separation in the images is **not** ruled out (§6 G).
- **The buffer size.**
- **The image table or the drawImage FFT on its own.** Changing either one can switch the bias on or off; what matters is the combination, through the copies P/W leaves inside the cell (4.20, 4.21).
- **The cell's distance from the edge of the drawImage FFT box** (the explanation first given in 4.18). The FFT box is centred on the cell, and forcing W changes the bias at fixed N (4.20).
- **Exact or sub-pixel galaxy positions.** Jitter of ±0.5–2 px barely changes m (4.14).
- **Small cells in general, or 'bigger is better'.** On the 9″ grid, 240, 244, 260–450, 498 and 500 px are clean, while 246–256 and 502–512 px are biased (4.16–4.19).
- **PSF stamp size / PSF table period.** See 4.8.
- **The R convention alone.** It explains part only.

## 6. Open questions and next steps

**A0. Why does a regular lattice in a 250-px cell break metacal?** Answered (4.20, 4.21). Use `test_layout.py` (layouts via `LAYOUTS`, ~10 minutes on 120 cores for 13 layouts).

- **Done:**
  - spacing scan (4.13)
  - fine scan around 9″ (4.15)
  - jitter and unsheared positions (4.14)
  - cell-size scan (4.16)
  - map of cell size × spacing (4.17)
  - FFT-margin rule at the larger table (4.18)
  - map for 500–512-px cells (4.19)
  - real FFT sizes, and forcing them in metacal (4.20)
  - copies of the cell at image level (4.21)
- **What we now know:**
  - the bias needs a lattice and one particular cell size
  - it varies steeply with the lattice geometry
  - the true shear feeds it in by distorting the lattice
  - **the cause is the ratio of metacal's image table P to its drawImage FFT W.** When P/W is not an integer, interpolation copies of the cell fold back into it at P mod W. For 250-px cells, W = 384 and P = 1024 put them ±128 px away, on other galaxies (4.20, 4.21)
  - on a lattice the copies of neighbours hit every galaxy alike; their offset doesn't follow the true shear but does follow metacal's shear (4.21)
- **Next:**
  1. **Confirm on the real images with the pipeline:**
     - run 250/150 with metacal patched so that W divides P (`pad_factor` 6, or W = 512), or simply 260/160 or 240/140
     - expect m2 like the 500-px runs (+0.001 to +0.004) instead of −0.0075
     - the cell-size route needs no code change
  2. **Report upstream** (ngmix / GalSim): for some cell sizes metacal's InterpolatedImage table and drawImage FFT are not commensurate, so interpolation copies land inside the image. A fix is to draw with an FFT size that divides the table, or to set `pad_factor` so the table is a multiple of it.
  3. **Optional checks of the copy picture:**
     - the copy amplitude against table size (pad 8 gives a much smaller bias than W = 384 with the default table)
     - the predicted P/W = 3 cells (N ≈ 375–384), which should be clean
     - whether noise images (fixnoise) pick up the same copies

**A. Why does deconvolution on a 250-px cell go wrong, and most at the centre?**

Use the `test_hard_cut.py` part b test bench (4.11): it reproduces m2 = −0.014 noise-free in ~4 minutes. Vary one thing at a time:

- **No neighbours:** an empty scene with only the target galaxy. Done: unbiased (4.12, `isolated`).
- **Cell size:** 300, 350, 400, 450 px. Done on the m level for the 9″ grid: all clean, as are 240–270 (4.16).
- **Distance from the centre:** the target's offset from the cell centre (saved per galaxy).
- **GalSim settings:** `folding_threshold`, `kvalue_accuracy`, `maxk_threshold`, `pad_factor`.
- **Interpolants:** lanczos15 against quintic, for both x and k.

Candidates:

1. The cell's sharp edge and zero padding, combined with the non-local deconvolution, with focusing at the one point equidistant from all four edges.
2. The FFT wrap of 2.05N: periodic copies of the cell at ±512 px interacting with the long-range response. *(The drawImage FFT is really 384 px for 250-px cells. Periodic copies did turn out to be the cause, in a different form: copies from the image table folded by the drawImage FFT, 4.21.)*
3. GalSim accuracy settings (`folding_threshold`, `maxk_threshold`, `kvalue_accuracy`) or `pad_factor`. ngmix uses GalSim's defaults.

Tests, each quick at image level with `test_metacal_convergence.py`:

- **(a) Cell-size scan:** set `CELL_SIZES = [250, 300, 350, 400, 450, 500]` to find where the images converge. That's the practical minimum cell size.
- **(b) Isolated galaxy:** one noise-free galaxy at the centre vs the full cell content. Is it the content, or the numerics?
- **(c) Padding:** pad the cell with real sky or noise (`InterpolatedImage noise_pad`) instead of zeros.
- **(d) GalSim settings:** tighter gsparams, by monkeypatching `ngmix.metacal.metacal._galsim_stuff_impl`.
- **(e) 2D map:** map of the 250 − 1000 difference, to see the spatial pattern (cross or ring).

**B. Does the image error explain m2 ≈ −0.015 quantitatively?** Measure wmom shapes and R on the 250 vs 1000 metacal images for galaxies near the centre. Alternatively, run `metadetect.do_metadetect` on the same cell centre at both sizes and compare per-object e2 and R22.

**C. Finish 500/50 and 500/100** (command in §3), then redo the analysis, diagnostic and plots.

**D. Small positive m2 ≈ +0.002 at cell 500.** Check its significance. The set-ups share images and noise, so their errors are correlated.

**E. Zero-shear ⟨e2⟩ lattice pattern** (4.3, last point). A separate sub-pixel effect; low priority.

**F. Recommendation for production.**
- **Grid simulations:** avoid the cell sizes where GalSim ends up with P/W = 8/3, within ~10 px below 256 and 512 (246–256, 502–512 tested), or patch metacal so that W divides P.
  - 500 px is clean but only 2 px from the 502–512 window. Sizes such as 240, or 260–480, sit well inside a clean range (P/W = 4).
  - Random positions are also clean at 250 px.
- **Real data:** random positions are clean at 250 px (4.12), so the real-sky case isn't expected to suffer this bias. Noise and blends could still pick up a small incoherent version; not tested.

**G. Does physical galaxy separation matter?** (colleague: a 15″ grid has no problem)

- Sub-sampling the 9″ grid changes nothing (4.10), so this needs images with sparser galaxies.
- **Cheapest test:** get the colleague's 15″ grid images (ask where they are) and run `test_metacal_convergence.py` on them, by changing `image_dir`. Do 250-px metacal images converge there?
- **Fuller test:** simulate a 15″ or 18″ grid and run 250/50 and 500/50.
- Physical separation fits the deconvolution picture: fewer neighbours, less light near the cell edges, and less content for the non-local deconvolution to couple.

## 7. Caveats for reproducing

- **Unsaved checks.** Several checks were one-off scripts in the session scratchpad and are **not saved**: the reweighting test, per-tag own vs common g2, the sub-5-px phase check, the GalSim set-up numbers in 4.7, and the point-source range in 4.9. The numbers above are their results. The drawImage FFT size in 4.7 came from the same flawed measurement as 4.16–4.18 and was wrong. The corrected sizes are recorded during real metacal calls by the saved `check_fft_sizes.py` (4.20).
- **Partial runs.** In `diag_position_in_cell.csv` and `results_cell_size.csv`, 500/50 and 500/100 are partial. They use fewer tiles, and those tiles aren't paired with the cell-250 tiles.
- **Diagnostic positions.** The diagnostic measures positions from the centre of the central region, cell/2, which is 0.5 px from metacal's shear origin, (cell−1)/2. Its "all" bin reproduces `results_cell_size.csv` to ≤3.5e-6.
