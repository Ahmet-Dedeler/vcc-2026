# 357 — Jurkat shRNA truth proxy: KILLED at inventory (2 panel intersects, need 8)

Status line: **NEGATIVE** (kill criterion fired at step 1; steps 2/4/5 not run).

## Known-answer reproduction (my own code, read-only, outputs in my `work/`)
`known357.py` recomputes 343's headline decomposition from 323's saved files:
mean NMAE **0.950836** (record 0.950836, |d| = 0.000000), magnitude ratio **0.205485**
(record 0.205), delivered |l|-weighted sign accuracy **0.587968** (record 0.588),
gate mass in our called set **0.346388** (record 0.346), zero-mass **6.4e-06**
(record "0.0000"), all-signs-right counterfactual **0.803387** (record 0.803),
over **58** eligible targets / 10,779-gene universe / 85 targets — all match.
Verdict: **REPRODUCED** (`work/known_answer_343.json`).

## What was measured (inputs, estimator, controls)
- Full chunked scan of `data/external/cmap/siginfo_beta.txt` (1,201,944 rows),
  filtered to `cell_iname in {JURKAT, HELA}` × `pert_type in {trt_sh, trt_xpr}`
  (`inventory357b.py`; quality cols `tas`, `cc_q75`, `ss_ngene`, `is_hiq`, `qc_pass`).
- Cross-checked against `data/external/cmap/genetic_sigs.parquet` (177,665 rows):
  its `pert_type` values are `trt_xpr`/`trt_sh.cgs` and it contains **zero**
  JURKAT or HELA rows — it is a consensus subset, not the inventory source.
- Intersected `cmap_name` against `data/challenge/pert_counts.csv` (300 targets;
  single panel file shared by A/B/C per AGENTS.md §2) and the 85 H1 eval targets
  from `experiments/323-current-family-validation/work/reference.npz['targets']`.
- Step 3 (read-only): searched `experiments/319-cap-variant/work/build*/` and
  `experiments/312-normrel-arm/work/build*/` for any 300×18,533 design array;
  inspected `scripts/make_submission_pm_cd4.py` header/docstring without running it.

## Results
| quantity | value |
|---|---:|
| JURKAT `trt_sh` sigs | **256** (reproduces worker 356's count) |
| distinct JURKAT `trt_sh` genes | **118** (mandate's "one per gene" is wrong: 35 genes have >1 hairpin, max 5) |
| JURKAT `trt_xpr` | **0** |
| HELA `trt_sh` / `trt_xpr` | **0 / 0** (confirms 356: no HELA knockdown bank) |
| JURKAT `trt_sh` intersect, 300-target panel | **2: POLD4, VAPB** |
| JURKAT `trt_sh` intersect, H1-85 eval targets | **3: LAP3, MFSD10, SLC5A6** |
| JURKAT `trt_sh` quality (all 256) | tas median 0.230 (IQR 0.189–0.285); cc_q75 0.51; ss_ngene 102; hiq 153/256; qc_pass 220/256; all pert_time = 120 h |
| overlapping sigs (1 hairpin each) | POLD4 tas 0.147/ss_ngene 57; VAPB tas 0.292/126; LAP3 tas 0.219/68; MFSD10 tas 0.139 (not hiq); SLC5A6 tas 0.409 but **qc_pass = 0** |

Full 118-gene list + per-gene sig counts for overlaps: `work/overlap_lists.json`;
filtered siginfo rows: `work/jurkat_hela_siginfo.parquet`.

Step 3 outcome: **no per-target design matrix exists on disk.** The build dirs hold
only CD4 bank arrays (`cd4_compact/<state>/{delta,log_fc,zscore,...}.npy` +
`target_row_map.csv`) and the pruned whitener (`whitener_union_pruned_def.npz`);
the 319 composition's A-panel design rows (300×18,533) are composed in memory at
build time by `scripts/make_submission_pm_cd4.py` from donor tables + CD4 banks +
whitener + normrel-cap weights + context-A basal. What would be needed for step 4:
run that builder chain (forbidden to this worker by rule 1) or re-implement its
staging in worker-owned code from the same frozen inputs. Per the mandate I did
**not** restore the 12.6 GB h5ad and did **not** extract any level-5 gctx columns.

## Verdict against the PLAN.md kill criterion
Kill criterion: fewer than 8 panel-A intersects → proxy too thin; report and stop.
**2 < 8 → KILLED.** With one hairpin each and tas 0.15–0.29, the two overlapping
sigs cannot support a per-target sign-agreement instrument, let alone a 1,000-draw
bootstrap. Steps 2 (gctx extraction), 4 (agreement table) and 5 (H1 calibration)
were not run — correctly, since a 2-gene proxy cannot adjudicate the 343
sign deficit either way.

## Limits
- Inventory only; no LINCS level-5 values were read, so landmark coverage of
  POLD4/VAPB signatures is unmeasured (moot at n = 2).
- `cmap_name` taken as the hairpin target symbol (standard for LINCS `trt_sh`;
  spot-checked: all values are HGNC-like symbols).
- Panel is one shared 300-gene file; "panel A targets" = the same 300.

## Exact commands to reproduce
```
NUMBA_NUM_THREADS=4 OMP_NUM_THREADS=4 uv run python experiments/357-jurkat-shrna-truth-proxy/inventory357b.py
NUMBA_NUM_THREADS=4 OMP_NUM_THREADS=4 uv run python experiments/357-jurkat-shrna-truth-proxy/known357.py
```
Both read-only toward other dirs; all writes in `experiments/357-jurkat-shrna-truth-proxy/work/`.

## Peak memory and wall time
Inventory scan (~465 MB siginfo, chunked): < 1 GiB RSS, ~2 min. Known-answer:
negligible (< 0.5 GiB, < 1 min). No run exceeded 4 GiB; no scratch files retained
beyond the two parquet/json outputs (~100 KB total).
