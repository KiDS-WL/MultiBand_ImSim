# Prompts of the metadetect cell-size session (5–8 Oct 2026)

These are the requests that drove the investigation summarised in [`cell_size_report.md`](cell_size_report.md), in order and paraphrased to their key points. Requests unrelated to the cell-size investigation (a pipeline PSF option, folder links, job-submission instructions) are left out. The full test log is [`cell_size_test_log.md`](cell_size_test_log.md) (originally `test_grid_1.md`), and section numbers (§4.x) refer to it.

## A. Set-up of the test

| # | Prompt (key points) | Outcome |
|---|---|---|
| 1 | Test how results change with metadetect `cell_size`/`central_size`. Draft the `.ini` files, run scripts and sbatch submission in `test_scripts/metadetect_cell_size`, for 250 with central 50/100/200 (150 already exists) and 500 with central 250/300/400. | Configs, `run_cell_size.sh`, `submit_cell_size.sh`. |

## B. First results

| # | Prompt (key points) | Outcome |
|---|---|---|
| 2 | Draft `test_plot_m.py`: m on the y axis, set-ups on the x axis. | `m_cell_size.png` |
| 3 | Also plot m against the buffer size (cell − central). | `m_buffer_size.png` |
| 4 | Plot m against the number of objects per central region. Why is it only m2? | `m_N_per_central.png` |
| 5 | You already have runs with different central sizes, which give the numbers. Why re-run anything? | No re-run; computed from the existing runs. |
| 6 | Compute m1 and m2 using R11 and R22 separately, keeping everything else the same; add them to `m_cell_size.png`. | Both conventions in `results_cell_size.csv`. |

## C. Looking for the cause: layout and selection

| # | Prompt (key points) | Outcome |
|---|---|---|
| 7 | R11 vs R22 can't explain it, especially the large negative m2 for small central sizes. List all possible causes. | List of candidate causes. |
| 8 | Write it as a script with an sbatch wrapper for torino. | `diag_position_in_cell.py` + `run_diag_position.sh` |
| 9 | For cause 1: why would the sample depend on the central size or step, given the steps are adjoining? | Shown: the tiling partitions the image, so the sample doesn't depend on it. |
| 10 | Visualise where galaxies fall within cells: rebuild positions from the ImSim layout (9″ grid, 0.5 deg²) and metadetect's cell cutting, stack all cells in a 2D histogram, and mark the central region with dashed lines. | `plot_cell_galaxy_positions.py`, `cell_galaxy_positions.png` |
| 11 | Analyse the `run_diag_position.sh` results. What are the "four sheared tags"? The zero-shear pattern clearly changes between cell/central combinations; why is that not an issue? | m against distance from the cell centre: the bias is at the centre (§4.3, §4.4). |
| 12 | In the stacked sheared positions there are diagonal red/blue line-ups (250/100, 250/150). Can they introduce m2? | Moiré fringes tested (§4.4). |
| 13 | Why would weighting remove a spatial pattern? Compare with how metadetect calibrates R: if it doesn't account for the pattern, it leaves m2. | Reweighting and R tests: the fringes are not the cause. |
| 14 | Is the ~86-px grid shift a concern? The position shift from shear should match the shape distortion, which is per galaxy and doesn't depend on tile position. | Position and shape shear consistent (§4.5). |

## D. The cell size is the cause; what in metadetect depends on it

| # | Prompt (key points) | Outcome |
|---|---|---|
| 15 | Rule out the cell size: run 500/200, 500/100 and 500/50, draft the configs and scripts. If they agree with 250, focus on one 500 set-up. | 500-px cells unbiased: the cell size is the cause (§4.6). |
| 16 | All runs finished; add the three new results to `m_cell_size.png`. | Updated figure. |
| 17 | Make `cell_galaxy_positions.png` side by side for 500/50, 250/50, 500/100, 250/100, 500/150, 250/150 to check the patterns are the same. | `cell_galaxy_positions_pairs.png`: identical galaxies. |
| 18 | How would the same galaxy at the centre of a 250-px and of a 500-px cell experience a different shear? | |
| 19 | Read the metadetect code carefully and summarise step by step what it does to a cell, at code level. Same images, same patterns, different results, so it must be how metadetect operates on a cell. | Code-level trace (§4.7). |
| 20 | I don't understand the table: what are you testing, and what do "different centres" and "converging" mean? | Clarified. |
| 21 | So is the difference introduced by PSF deconvolution/reconvolution, the only cell-size-dependent step? Test PSF stamp sizes (e.g. 32×32, 16×16) for convergence. | `test_metacal_convergence.py` (§4.8). |
| 22 | Recap everything tested so far into `cell_size_test_log.md`, so I can pick up the conclusions and remaining hypotheses in a new session. | `cell_size_test_log.md` |

## E. It needs a lattice

| # | Prompt (key points) | Outcome |
|---|---|---|
| 23 | A colleague finds a 15″ grid has no problem. Check whether galaxy separation matters by selecting every second galaxy (twice as far apart) from the existing products, without re-running; make a plot like `m_cell_size.png`. | `test_grid_separation.py`, `m_grid_separation.png`: no change (§4.10). |
| 24 | "Less light from neighbours" is wrong: grid stamps cut each galaxy's light, so none leaks to neighbours. What does that imply for deconvolution? Is the hard flux cut a problem? | Hard-cut hypothesis. |
| 25 | Let us test this. How long would it take? | Noise-free test bench `test_hard_cut.py`. |
| 26 | Yes. | Run: the hard cut is not the cause, and the bias is reproduced noise-free (§4.11). |
| 27 | That contradicts the problem appearing only for the 9″ grid. Part a implies any grid size, and part b random positions, would have the problem, right? | Corrected: only 9″ had been tested. |
| 28 | Yes, test these. | `test_layout.py`: needs a lattice; random positions are clean (§4.12). |
| 29 | Do a spacing scan. | Resonance-like dependence on spacing (§4.13). |

## F. Pinning it down

| # | Prompt (key points) | Outcome |
|---|---|---|
| 30 | Plot these results, then do all the next steps. For each test update `cell_size_test_log.md` and save plots to the folder, not just tables. | Jitter and unsheared grid, fine spacing scan, cell-size scan (§4.14–4.16); `run_layout_tests.sh`, `plot_layout_scans.py`. |
| 31 | Produce a 2D map of m against cell size (244–256 px) and spacing (43–47 px). | `m_map_cell_spacing.png` (§4.17) |
| 32 | Submit it (the 502–512-px test). | Bias also at 502–512 px (§4.18). |
| 33 | Make the same map for the 500-px tests, and run 250-px cells with different FFT sizes to pin it down. | `m_map_cell_spacing_500.png`, `m_vs_forced_fft.png`, `metacal_fft_ghosts.png`. Cause found: the ratio of metacal's image table to its drawImage FFT. The earlier FFT-size numbers were corrected (§4.19–4.21). |

## G. Confirmation and write-up

| # | Prompt (key points) | Outcome |
|---|---|---|
| 34 | Run the next step (the pipeline on the real images). If confirmed, write `cell_size_report.md`: the cause with explanations and evidence pointing to the key plots, and why only grid simulations are affected, not realistic positions. | 260/150 run: m2 = +0.0014 instead of −0.0075 (§4.22); `cell_size_report.md` |
| 35 | Write `cell_size_prompt.md` summarising all my prompts (key points, not exact wording). | This file. |

## Recurring instructions

These came up repeatedly and are worth keeping for future sessions:
- Don't re-run what existing products can answer (#5, #23).
- Use torino sbatch wrappers for cluster jobs (#8).
- Reason at code level, and don't accept loose explanations: several hypotheses were challenged and corrected (#9, #13, #14, #18, #24, #27).
- For every test, update `cell_size_test_log.md` and save figures to the folder, not only tables (#22, #30).
