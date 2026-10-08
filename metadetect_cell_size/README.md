# Metadetect cell size and the shear bias in grid simulations

This folder records an investigation (5–8 October 2026) into why metadetect, run with 250-pixel cells on the MultiBand_ImSim grid simulations, gives a multiplicative shear bias of m2 ≈ −0.004 to −0.014, while 500-pixel cells give |m| ≤ 0.0025.

## Key finding

**The bias is a numerical artefact of how metacal renders a cell with GalSim. It needs both a particular cell size and galaxies on a regular grid.**

- **The two FFT sizes.** Metacal turns each N-pixel cell into a GalSim `InterpolatedImage`, whose Fourier table is P = good_fft_size(4N) pixels. It draws the sheared, reconvolved image back with an FFT of W ≈ good_fft_size(N + 10) pixels.
- **The bad cell sizes.** For cells within ~10 px below 256 or 512 px (246–256, 502–512 px; 250 included), P/W = 8/3 instead of an integer.
- **The copies.** The interpolation in the table then leaves faint copies of the cell's content. The drawing FFT folds them back into the cell, 128 px away for 250-px cells.
- **Why grids.** On a regular grid, a copy of a neighbour lands at the same offset from every galaxy. Its offset doesn't follow the true shear but does follow metacal's artificial shear, so metacal's response no longer matches the true one.
- **Why not real data.** At realistic, random positions the copies land at random offsets and average out. On a noise-free test bench, random positions give m2 = −0.0001 ± 0.0011 with the same 250-px cells, against −0.014 on the 9″ grid.

**Confirmed on the real grid images:** with 260-px cells (P/W = 4) and the same central size, so the same galaxies as the 250/150 baseline, m2 = +0.0014 ± 0.0007 instead of −0.0075. That matches the 500-px cells.

**Recommendation:**
- For grid simulations, use a cell size outside the bad windows (e.g. 240, or 260–480 px), or random galaxy positions.
- Or make metacal draw with an FFT size that divides its table. This fixed it on the bench, but needs a change in ngmix.
- Real data is not expected to be affected.

## Contents

| file | what it is |
|---|---|
| [`cell_size_report.md`](cell_size_report.md) | **Start here.** The cause, the evidence with the key figures, why only grid simulations are affected, the confirmation run, and recommendations. |
| [`cell_size_test_log.md`](cell_size_test_log.md) | The full running log of every test, in order (§4.1–4.22), including ruled-out hypotheses and corrections. |
| [`cell_size_prompt.md`](cell_size_prompt.md) | The requests that drove the session, summarised in order. |
| [`plots/`](plots/) | The key figures (below). |
| [`scripts/`](scripts/) | All scripts and configs of the investigation (below). |

### Key figures (`plots/`)

| figure | shows |
|---|---|
| [`m_cell_size.png`](plots/m_cell_size.png) | Real images: m per cell/central set-up. 250 px is biased; 260 and 500 px are not. |
| [`diag_position_radial.png`](plots/diag_position_radial.png) | Real images: m against distance from the cell centre. The dip at the centre appears only for 250-px cells. |
| [`m_vs_grid_spacing.png`](plots/m_vs_grid_spacing.png), [`m_vs_grid_spacing_fine.png`](plots/m_vs_grid_spacing_fine.png) | Test bench: m against grid spacing; sharp, resonance-like features at 250 px. |
| [`m_vs_jitter.png`](plots/m_vs_jitter.png) | Test bench: jittered grids, an unsheared grid, and random positions. |
| [`m_vs_cell_size.png`](plots/m_vs_cell_size.png) | Test bench: 9″ grid in cells of 240–500 px, with GalSim's FFT sizes. |
| [`m_map_cell_spacing.png`](plots/m_map_cell_spacing.png), [`m_map_cell_spacing_500.png`](plots/m_map_cell_spacing_500.png) | Test bench: m over cell size × spacing near 250 and 500 px. |
| [`m_vs_table_margin.png`](plots/m_vs_table_margin.png) | Test bench: m against P − 4N for all cell sizes tested. |
| [`m_vs_forced_fft.png`](plots/m_vs_forced_fft.png) | Test bench: GalSim's FFT sizes forced inside metacal switch the bias on and off at a fixed cell size. |
| [`metacal_fft_ghosts.png`](plots/metacal_fft_ghosts.png) | Metacal images: the copies of the cell, shifted by ±128 px. |

### Scripts (`scripts/`)

- **Pipeline runs on the real images:**
  - `submit_cell_size.sh`, `run_cell_size.sh` and the ten `config_cell*_central*.ini` set-ups (metadetect task 6_2 on the collaborator's images)
  - `run_analyse.sh` → `analyse_cell_size.py` (R, m, c), then `test_plot_m.py`
  - `run_diag_position.sh` → `diag_position_in_cell.py` (m against position in the cell)
- **Noise-free test bench:**
  - `test_hard_cut.py`, the bench itself. It draws scenes with ImSim's code and runs metacal as metadetect does; the options `MCAL_FFT_SIZE` and `MCAL_PAD_FACTOR` force GalSim's sizes inside metacal.
  - `test_layout.py`, the neighbour layouts. `run_layout.sh`, `run_layout_tests.sh <test>` and `run_fft_size_test.sh` are the sbatch wrappers.
  - `plot_layout_scans.py` makes the bench figures.
- **GalSim internals:** `check_fft_sizes.py` (sizes recorded during real metacal calls) and `check_fft_ghosts.py` (the copies at image level).
- **Supporting tests of ruled-out hypotheses:** `plot_cell_galaxy_positions*.py`, `test_grid_separation.py` + `run_grid_separation.sh`, `run_hard_cut.sh`, `test_metacal_convergence.py`.

The scripts write their outputs next to themselves, and the sbatch wrappers (torino partition) write job logs to `scripts/logs/`. The outputs of the original runs are kept in `../test_scripts/metadetect_cell_size/`: summary CSV/NPZ files, job logs and supporting figures. Those include `results_cell_size.csv` and `layout_response*.csv`, which the plotting scripts read.

**Software:** MultiBand_ImSim v1.3.0, metadetect 0.13.0 (vendored in `modules/`), ngmix 2.4.1, GalSim 2.8.5 and sep 1.4.1, run on SLAC S3DF.

## Acknowledgement

This investigation was carried out as an interactive session with **Claude Code**, Anthropic's AI coding assistant (model Claude Opus 5.5). The author directed the work, set the tests, and questioned the reasoning at each step; the prompts are summarised in [`cell_size_prompt.md`](cell_size_prompt.md). Under the author's direction, Claude Code wrote and ran the scripts, made the figures and drafted the documents in this folder.
