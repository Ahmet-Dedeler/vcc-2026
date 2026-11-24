# 346 — PLAN: does a same-line (H1) calibration map, fitted on the recipient's own public
# non-evaluation perturbations, move the six competition members on H1?

Written before any fit or emission. 2026-09-10 UTC. Assignment 346-same-line-calibration-h1.

## Hypothesis
341 showed that a full ridge map from donor consensus to a recipient's response, fitted on that
recipient's **own training-target** labels, raises |lfc|-weighted sign accuracy by +0.033 on the
six-CRISPRi-line source frame. 343 showed the H1 board gap is exactly that sign quantity, at the
design stage (0.600 weighted sign accuracy for the designed cap15 effect). H1: applying such a
map, fitted on the **public Pan H1 study** (GSE295214) restricted to non-evaluation, non-panel
targets, to the 85 H1 evaluation targets' donor consensus, and blending it into 323's cap15
design at fixed row norm, raises the six-member H1-local overall score.

This is **known-line, held-target, cross-study** development evidence. It is not unseen-line
generalization and not an A/B/C forecast. Its purpose is to price the information regime that
context A would need (Nadig Jurkat non-panel rows) before anyone builds for A.

## Inputs (exact paths)
- Recipient training labels: `experiments/115-gse295214-h1/work/compact-donor/`
  (source id `gse295214_h1`; the compact donor extraction of `data/external/gse295214/`
  `pseudo_bulk_perturbation_database.h5ad`). 2,354 rows x 18,533 canonical genes,
  `delta.npy` = log1p(TS50000 perturbed) - log1p(TS50000 control), `b_ctrl.npy`,
  `measured_mask.npy` (row-constant, 16,268 measured coordinates).
- Donor banks (inputs X): `data/donors/k562_gwps`, `data/donors/orion_hct116`,
  `data/donors/orion_hek293t`. See "Deviation from the registered donor set" below.
- Evaluation targets, reference DE, design effects, scorer:
  `experiments/323-current-family-validation/` — `work/reference.npz` (85 targets,
  18,080-gene scored axis, `keep` = 10,779-gene DE universe), `work/reference_de.parquet`,
  `work/effects.npz` (`cap15`, 85 x 18,533 on the canonical axis), `work/metrics-cap15.json`,
  `work/per-target-cap15.csv`, `work/h1-anchors.json`, `evaluate323.py`, `run323.py`,
  `report323.py`, `exact_de323.py`.
- Panel exclusions: `data/challenge/pert_counts.csv` (300 A/B/C targets).
- Estimator primitives reused verbatim: `experiments/341-known-line-calibration/core341.py`
  (`diagonal`, `fit`, `predict`, `sign_endpoints`, `interval`) and
  `experiments/334-text-conditioned-transfer/core334.py` (`rms_direction`, `features`,
  `source_blocks`, `diagonal`, `derangement`, `endpoints`). Loader pattern from `run334.load`.
- Emission/scoring template: `experiments/342-hamilton-current-family/run342.py`,
  `report342.py`.

## Deviation from the registered donor set (stated up front)
The registration names K562-GWPS, both Orion banks and Nadig HepG2 as donors, from
`docs/sprint-2026-09-09/NEXT-METHOD-SCREEN.md`. Recomputed here, that four-bank set does
reproduce the screen exactly: **1,178** shared non-panel non-evaluation targets, **5,567**
common canonical genes after excluding the training/evaluation target union, **30/85**
evaluation targets with a direct Pan row (`work/coverage.json`). But Nadig HepG2 has rows
for only **21** of the 85 evaluation targets, so only **20/85** evaluation targets would have a
four-donor consensus at test time; 65 of 85 design rows would be left unmodified and the arm
would be nearly a no-op. The four-donor construction is therefore **not usable at test time**.

Primary construction uses the **three production donors** (K562-GWPS, Orion HCT116, Orion
HEK293T) — exactly 323's `DONORS`, i.e. the same donor consensus the production design is
built from. Recomputed counts: **1,965** training targets, **5,275** covered genes (after
also intersecting with the Pan measured mask), **80/85** evaluation targets with a full
three-donor consensus, **4,836** of the 10,779 DE-universe genes covered, and the covered
axis carries a median **0.621** of each cap15 row's L2 norm. The 5 evaluation targets without
a row in all three donors keep their unmodified design row in every arm; they are named in
RESULTS.md.

## Estimator
Let `direction` be the unit-RMS raw direction (334's `rms_direction`) of each bank's delta row
restricted to the covered gene axis; lines = (k562_gwps, orion_hct116, orion_hek293t,
gse295214_h1) with outer = gse295214_h1. On the 1,965 training targets, 334's `source_blocks`
gives three inner recipients, each with the mean direction of the other two donors as X
(`features[...,0]`) and its own direction as Y. The test query is the three-donor mean
direction of the 80 covered evaluation targets. Arms:

1. `crossline_diag` — per-gene diagonal (`core334.diagonal`) fitted across the three inner
   donor recipients on the training targets; applied to the query. 341's cross-line baseline.
2. `same_line_full` — `core341.fit(x, y)` ridge, x = the three stacked leave-one-out donor
   consensuses on the training targets, y = the Pan H1 response tiled three times.
3. `same_line_diag` — the diagonal component of (2) alone.
4. `labels_shuffled_full` — (2) with the Pan H1 training labels permuted across targets by
   `core334.derangement`.
5. `source_mean` — the raw donor consensus (design-stage reporting only).

## Integration into the production design (exact statement)
`C` = the 5,275 covered gene indices on `effects.npz['genes']`. For evaluation target i with a
three-donor consensus, `d = cap15[i, C]`, `p` = that arm's prediction for target i on `C`,
`u_d = d/||d||_2`, `u_p = p/||p||_2`, and for a blend weight beta

    m       = (1-beta)*u_d + beta*u_p
    new[C]  = m * (||d||_2 / ||m||_2)
    new[j]  = cap15[i, j]   for every gene j not in C

so beta=0 reproduces cap15 exactly, beta=1 replaces the design direction on `C` by the
calibrated direction, and every arm has **exactly** the original row norm on `C` and off `C`.
Amplitude — a closed axis (317, 227) — is therefore not what changes; only direction within
the covered subspace. Targets without a three-donor consensus keep `cap15` unchanged.

## Six-member pricing
Arms emitted and scored with 323's emission constants and `exact_de323`/cell-eval2, run from
this directory (never writing into 323's `work/`), sequentially, one at a time:

| arm | effects |
|---|---|
| `cap15` | unmodified 323 cap15 — KNOWN ANSWER |
| `diag_cross_b1` | arm 1 at beta=1 — control for the integration itself |
| `same_ridge_b1` | arm 2 at beta=1 |
| `same_ridge_b05` | arm 2 at beta=0.5 |
| `shuffled_b1` | arm 4 at beta=1 |

Known-answer requirement (protocol rule 3): before any modified arm is interpreted, the
`cap15` arm must reproduce `experiments/323-current-family-validation/work/metrics-cap15.json`
exactly — all six raw members, the emitted-count digest, the panel-total digest, the MSE legs,
every array in `stats-cap15.npz`, and every row of `per-target-cap15.csv`. A second known
answer: my re-implementation of 343's `signs343.py` design-stage weighted sign accuracy must
return 0.600185097698462 for cap15 (343's `work/signs-cap15.json::design_w_all`). A third:
the four-donor metadata counts must equal 1,178 / 5,567 / 30 (done, `work/coverage.json`).

Reporting: paired 10,000-target bootstrap (seed 20260909) of the H1-local overall and each of
the six raw members, as in `report342.py`, with the MSE panel ratio and its claim cap
recomputed in every resample.

Design-stage (pre-emission) reporting for every arm: |lfc|-weighted sign accuracy of the
designed effect on the reference's significant gate, computed exactly as 343's `signs343.py`
(`design_w_all`), plus the same restricted to the covered axis.

## Controls and falsifiers
- **Shuffled-label control** (`shuffled_b1`): same architecture, Pan H1 training labels
  deranged across targets. Any gain it reproduces is architecture/integration, not information.
- **Integration control** (`diag_cross_b1`): the integration applied to a map that uses no H1
  label at all. Isolates "replacing the design direction with a donor-consensus direction"
  from "replacing it with an H1-calibrated direction".
- **Poison test**: the fit is re-run with the Pan H1 rows of the evaluation targets replaced by
  noise; every fitted parameter (`gain`, `scale`, `beta`) and every prediction must be
  byte-identical. 30 of the 85 evaluation targets do have direct Pan rows; they are excluded
  from fitting and the poison test proves it.
- **Subset-linearity identity**: the mean over the three leave-one-out subset predictions must
  equal the prediction on the three-donor mean to < 1e-12 (341's check).
- **Norm identity**: every modified row's norm on `C` must equal cap15's to < 1e-10, and every
  off-`C` coordinate must be bit-identical to cap15.
- **beta=0 identity**: the integration at beta=0 must return cap15 bit-identically.

## Kill criteria (the regime is NOT worth building for context A if any fires)
1. The `same_ridge_b1` arm's H1-local overall paired 95% interval does not exclude zero above
   **+0.002** (i.e. lower bound <= +0.002).
2. The shuffled-label control `shuffled_b1` reproduces **>= half** of `same_ridge_b1`'s overall
   gain over `cap15`.
3. The design-stage |lfc|-weighted sign accuracy gain of `same_ridge_b1` over `cap15` is
   **< +0.015**.

RESULTS.md will state which fired. A clean negative with these controls is a good result.

## Budget and safety
5 hours. Five emission arms at ~142 s and ~5 GiB each, run strictly sequentially, with
`memory_pressure | tail -1` >= 50% free checked before each (a 22 GiB production build,
pid 9220/10210, is running). `NUMBA_NUM_THREADS=4`, `OMP_NUM_THREADS=4`. Nothing is written
outside `experiments/346-same-line-calibration-h1/`. No git command, no submission, no builder,
no reserve read, no download.
