# 372 — Jurkat same-line calibration artifact built (259 panel rows, held-out sign +0.042, all gates pass)

**Status: COMPLETE — ENGINEERING_ONLY (artifact built and verified; no emission, no score, no production change).**

2026-09-10 UTC (≈00:48–00:53). Recipient = Nadig Jurkat (`data/donors/nadig_jurkat`), donors = the same three production banks as 346. Same estimator as 346 (verbatim `core341`/`core334` imports, same λ rule, same poison + shuffled controls). Output `work/ctx_a_calibrated.npz` + JSON sidecar. No `scripts/`/`src/` edits, no builder, no submission, no reserve, no download.

## Known-answer reproduction (own code, before any new measurement)
Re-ran 346's three-donor axis construction independently (Pan H1 `measured_mask.npy` row-constancy checked in slabs, same exclusion rule):
- training = **1,965**, covered = **5,275**, test-with-consensus = **80/85** — exact match to 346 `work/fit.json` / RESULTS.md.
- uncovered = **AASS, CNN3, COL6A1, HMGCS1, MED23** — exact match to 346 §"What was measured".
- 346 design-sign value 0.600185097698462 was not re-emitted (emission not in scope); the axis identity above is the gate this artifact depends on. Match: **PASS, proceed.**

## What was measured (inputs, estimator, controls)
- **Recipient:** `data/donors/nadig_jurkat/` — `delta.npy` 2,270×8,881, `genes.txt` 8,881, `perts.txt` 2,270, `b_ctrl.npy`. Raw cells `data/external/nadig2025/jurkat_raw_singlecell.h5ad` **unused** (346 needs no raw-derived quantity; recorded as unused per assignment).
- **Donors (X):** `data/donors/k562_gwps` (9,787 perts × 8,248 genes), `data/donors/orion_hct116` (17,767 × 38,606), `data/donors/orion_hek293t` (17,977 × 38,606). Confirmed donor set = 346's three production donors (NOT the four-bank registration; same deviation rationale: HepG2 would leave most query rows unmodified).
- **Panel:** `data/challenge/pert_counts.csv` 300 symbols, `data/challenge/gene_names.csv` 18,533 canonical.
- **Gene axis (346 rule minus Pan mask; Jurkat has no `measured_mask.npy`):** (3 donors ∩ Jurkat ∩ canonical) − training symbols − 300 panel symbols.
- **Consensus:** `core334.rms_direction` unit-RMS directions; `baseline = mean(donor directions)`; query = baseline[test]. Same as 346.
- **Ridge:** `core341.fit(x, y)` verbatim, λ = n_rows = 3 × n_train, RMS-scaled features; x = 3 stacked leave-one-out donor consensuses on training targets, y = Jurkat response tiled 3×. Same as 346.
- **Loader:** 346 `fit346.py::bank` shape (exact-name row/column selection + 3×3 scalar-index fixture assertion). Same.
- **Blend:** 346 `effects346.blend` verbatim (β=0.5, fixed-norm on covered axis), with the design row replaced by the donor-consensus row (no A design exists upstream of the builder). `matrix` = blended rows; `prediction_raw`/`consensus`/`consensus_norm` also stored so the builder can reconstruct either.
- **Controls:** held-target poison replay (1e6 noise on query rows → byte-identical fits/predictions); shuffled-label map (`core334.derangement(n_train, SEED=3460910)`); within-bank held-out sign check (two half-splits, `core341.sign_endpoints` weighted sign, `same_line_full` vs `crossline_diag`, shuffled as null).

## Results tables (units stated)

### Counts (targets, axis)
| quantity | value |
|---|---:|
| Jurkat perts total | 2,270 |
| Jurkat ∩ 3-donor perts, minus 300 panel | **2,006** training targets (NEXT-METHOD-SCREEN quoted 1,986; delta +20 = different donor build/exclusion vintage, recorded not tuned) |
| Jurkat ∩ 300 panel | **0** (asserted; `overlap_j_panel = []`) |
| 4-bank ∩ canonical (pre-exclusion) | 6,783 genes |
| …of which training symbols on axis (dropped) | 1,613 |
| …of which panel symbols on axis (dropped) | **219 of 300** (listed in sidecar `panel_genes_dropped`) |
| **Covered axis final** | **4,951** genes (26.7% of canonical) |
| Panel targets with full 3-donor consensus (query, panel order) | **259 / 300** |
| Panel targets with no consensus (unmodified downstream; listed) | 41: ABCD1, ANKRD52, BRD1, C5orf22, CAPRIN2, CDK8, DLAT, EPHB2, GLYR1, HEATR6, … (full list in sidecar) |
| Covered genes containing a panel symbol | **0** (gate) |

**DEVIATION (assignment §4 "100 A targets"):** the A/B/C panel is 300 symbols shared across contexts (`data/challenge/manifest.json`: 300 per context); 259 have a buildable 3-donor consensus and 0 are in the Jurkat bank. There is no buildable "100". Shipped all **259** in `pert_counts.csv` order; targets array order verified equal to the panel-order-filtered list. Do not subset to 100 downstream without a new instruction.

### Fit (same λ rule as 346)
| quantity | value |
|---|---:|
| ridge rows (3 × 2,006) | 6,018 |
| ridge features (covered genes) | 4,951 |
| λ (core341: n_rows) | **6,018** |
| normal-equation max residual | 2.865e-11 |
| subset-linearity max difference (341's check) | 1.377e-14 (< 1e-12 gate) |
| poison replay (query rows → 1e6 noise) | **PASS, byte-identical** gains/scales/betas/predictions/cross_gain |
| off-axis / design handling | n/a (no A design in scope; blend norm identity below instead) |
| blend norm identity max \|‖matrix‖ − ‖consensus‖\| | 1.42e-14 |

### Within-bank held-out weighted sign (Jurkat truth, `core341.sign_endpoints`, diagnostic not a score)
| split (fit n / eval n) | same | diag | shuffled | gain same−diag | 95% paired-target CI | frac improved | top-64 gain | shuffled−diag |
|---|---|---:|---:|---:|---|---|---:|---:|
| even/odd (1,003 / 1,003) | 0.61650 | 0.57446 | 0.55148 | **+0.04204** | [+0.03891, +0.04523] | 0.804 | +0.06418 | −0.02298 |
| random half (1,003 / 1,003, seed 3720910) | 0.61702 | 0.57504 | 0.54767 | **+0.04198** | [+0.03879, +0.04528] | 0.809 | +0.06389 | −0.02736 |

Comparison: 341's Jurkat gain vs diagonal is **+0.03380** [+0.03014, +0.03762] on its 512-coordinate source frame. Ours is **+0.042** on the 4,951-gene production axis, same sign, larger magnitude, both replications agree to 6e-05, shuffled is negative both times. Passes the sanity check (recipient's own labels carry per-target direction information on held-out targets; controls behave as in 346/341).

### Sanity (abort gates, all pass)
| gate | value | verdict |
|---|---|---|
| any NaN/Inf in matrix/raw/consensus | none | PASS |
| cosine(matrix, consensus) per row | median **0.8822**, min 0.8297, max 0.9111, mean 0.8804 | moderately positive, not ~1, not ~0: PASS |
| cosine(raw prediction, consensus) per row | median **0.5565**, min 0.3769, max 0.6601 | moderately positive: PASS |
| covered fraction of canonical | 4,951/18,533 = 0.2671 | reported |
| genes carry panel symbols | 0 | PASS |
| Jurkat panel overlap | 0 | PASS (asserted in code) |
| targets in panel order | verified equal | PASS |

### Artifact (exact keys)
`work/ctx_a_calibrated.npz` (30,987,488 B, sha256 `d270ded1…080dc27` in `work/sha256sums.txt`):
`targets` (259,) panel-order symbols · `genes` (4,951,) official symbols, no panel genes · `matrix` (259×4,951) float64 blended β=0.5 · `prediction_raw` (259×4,951) · `consensus` (259×4,951) · `consensus_norm` (259,) · `beta` 0.5 · `mode` 'blend_rows' · `context` 'A' · `norm_rule` 'fixed_norm_on_covered'.
Sidecar `work/ctx_a_calibrated.json`: bank shas (delta/perts/genes/b_ctrl for all four banks), λ/counts above, poison result, held-out table, cosine summaries, coverage fraction, the 41 missing targets, the 219 dropped panel genes, and the three deviations. `beta` is the 346 post-hoc-adjacent operating point the assignment fixed (346: β=0.5 overall +0.0071 [−0.004, +0.016], β=0.25 +0.0065; neither resolved — the artifact takes no position on which blends better on A).

## Verdict against the PLAN.md gates
No idea-kill criterion was registered (engineering artifact, not a priced arm). All five build gates pass: (1) known-answer axes exact, (2) panel overlap 0 asserted, (3) poison byte-identical, (4) held-out gain positive with negative shuffled control on both splits, (5) sanity cosines/coverage/NaN/order all pass. **Ship the artifact as-is; do not emit, score, or build from it in this experiment.**

## Limits
- Same-line, held-target, **cross-study** logic (Nadig Jurkat → context A assumed Jurkat via 041/315, not certainty). Not unseen-line generalization, not an A/B/C forecast, no board conversion.
- 259/300 panel targets predictable; 41 missing have no 3-donor consensus and are absent from the artifact (listed). 219/300 panel genes were on the pre-exclusion axis and removed; the covered axis holds 26.7% of canonical genes — everything outside keeps whatever the builder puts there.
- 346's β=0.5/0.25 trade-off was priced on H1 (85 targets, emitted counts); no emission was run here, so no member claim transfers. 346's registered verdict (full replacement not worth building) is not overturned by this artifact.
- Ridge λ/features/architecture are 341's, untuned. Held-out is same-bank/same-study (optimistic vs cross-study). Training reuse/correlation confounds remain as in 341/346.
- Peak RSS **12.22 GiB exceeded the 12 GiB cap by 0.22 GiB (1.8%)** transiently during the joint fit+poison pass; system free stayed 76–77%, no production impact. Future fits should free `raw` before the poison replay or fit halves serially.

## Exact commands to reproduce
```
cd ~/Code/virtual-cell-challenge-2026
memory_pressure | tail -1   # require ≥50% free (observed 77%)
NUMBA_NUM_THREADS=4 OMP_NUM_THREADS=4 uv run python -c "import csv,pathlib,numpy as np; ..."  # known-answer axes check (see build372.py::build path; expect 1965/5275/80)
NUMBA_NUM_THREADS=4 OMP_NUM_THREADS=4 uv run python experiments/372-jurkat-calibration-artifact/build372.py  # ~76 s; writes work/ctx_a_calibrated.npz + .json + fitted-*.txt + predictions372.npz
sha256sum experiments/372-jurkat-calibration-artifact/work/ctx_a_calibrated.npz  # expect d270ded1…080dc27
uv run python -c "import numpy as np; z=np.load('experiments/372-jurkat-calibration-artifact/work/ctx_a_calibrated.npz'); print(z['matrix'].shape, np.isfinite(z['matrix']).all())"  # (259, 4951) True
```
`model-*.npz` (3×179 MB, shas in `work/sha256sums.txt`) were deleted after hashing; they regenerate deterministically via `build372.py`. Nothing outside `experiments/372-jurkat-calibration-artifact/` was written; no git, no builder, no submission, no reserve, no download.

## Peak memory and wall time
| step | wall | peak RSS | notes |
|---|---|---:|---|
| known-answer axes recompute | <10 s | <1 GiB | 1965/5275/80 exact |
| `build372.py` (assemble + fit + poison + 2 held-out fits + artifact) | **75.7 s** | **12.22 GiB** (`ru_maxrss`; 1.8% over the 12 GiB cap, disclosed above) | system free 76–77% throughout; `NUMBA_NUM_THREADS=4 OMP_NUM_THREADS=4` |
| verify + hash | <15 s | <1 GiB | artifact 31 MB, predictions 56 MB retained; models deleted |
