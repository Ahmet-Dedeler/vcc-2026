# 046 — Supervised target-context adaptation is a no-go

**Date:** 2026-08-22

## Verdict

Do not use the fitted Jurkat map. The initial nested holdout looked extremely
positive, but a decisive negative control showed that paired target-context
supervision was not responsible for the gain.

The real finding is narrower and still interesting: a target-effect-derived
low-rank subspace acts as a strong generic retrieval transform. That is a new
geometry lead, but it is not a validated context-adaptation method and cannot
justify conditioning context A on Jurkat.

## Protocol

- Outer targets: K562-gwps, Orion-HCT116, Orion-HEK293T.
- Input: equalized consensus of the other two donors.
- Fit panel: 600 common non-VCC perturbations.
- Outer panel: 259 shared VCC perturbations.
- Overlap: exactly zero.
- Score axis: 7,312 genes after excluding all 300 VCC target genes.
- Adaptation axis: 5,582 genes shared with the auxiliary screens, after also
  excluding every one of the 600 training target genes.
- Selection: five perturbation-held-out folds.
- Models: scalar shrink, diagonal shrink, rank-4/8/16/32 randomized-PCA ridge,
  and residual-program maps.
- No full eigendecomposition, no `cell-eval2`, no single-cell counts, and no
  challenge submission.

## Initial result

Both raw and gene-normalized source representations passed the predeclared
outer-transfer gate.

| source | baseline | target-CV selected | delta | donor-nested | delta |
|---|---:|---:|---:|---:|---:|
| raw | 0.71389 | 0.75592 | **+0.04203** | 0.75230 | **+0.03841** |
| gene-normalized | 0.71985 | 0.75938 | **+0.03952** | 0.75801 | **+0.03816** |

Raw, target-CV-selected outer deltas:

| target | baseline | adapted | delta |
|---|---:|---:|---:|
| K562-gwps | 0.74278 | 0.78459 | +0.04181 |
| Orion-HCT116 | 0.69552 | 0.75278 | +0.05726 |
| Orion-HEK293T | 0.70338 | 0.73040 | +0.02703 |

An independent split/randomized-basis seed reproduced the raw macro gain:
**+0.04242**, with the same selected ranks and every donor positive.

## The control that invalidated the interpretation

The fixed raw rank-32, ridge-0.3 map was rerun with the target perturbation rows
randomly permuted before fitting. The target effect *set* and its low-rank
subspace stayed the same, but the source-to-target pairing was destroyed.

| outer target | paired delta | permuted-label delta | best wrong-map delta |
|---|---:|---:|---:|
| K562-gwps | +0.04168 | **+0.04810** | +0.03006 |
| Orion-HCT116 | +0.05944 | **+0.07412** | +0.07352 |
| Orion-HEK293T | +0.02514 | **+0.04747** | +0.06998 |
| **macro** | **+0.04209** | **+0.05656** | **+0.05786** |

The correctly paired model lost to the permuted model in all three contexts.
Wrong-context maps also beat the paired map in two of three contexts. Therefore
the experiment did not learn target-specific perturbation adaptation.

The likely mechanism is covariance geometry: the PCA basis is learned from the
target donor's non-VCC effect distribution, and merely rotating/shrinking
queries in that subspace improves retrieval. Row pairing is unnecessary and
slightly harmful here.

## Jurkat branch

A generic-to-Nadig-Jurkat residual-program artifact was briefly generated
after the initial gate passed. It is now explicitly quarantined under
`artifacts/` and must not be applied to context A. Jurkat has **zero** target
overlap with the 300-gene VCC panel, so there was no direct panel leakage; the
failure is conceptual, not leakage.

## Resource use

- Main grid peak RSS: **1.09 GiB**
- Independent seed peak RSS: **0.75 GiB**
- Negative-control peak RSS: **0.65 GiB**
- Hard guard: 4 GiB main / 2 GiB controls
- No heavy process remains running

## Artifacts

- `run.py` — nested perturbation/donor evaluation
- `results.json` — raw and gene-normalized main grid
- `results_seed146.json` — independent raw replicate
- `results_raw_confirm.json` — seed-46 regeneration
- `audit_controls.py` / `results_controls.json` — permutation and wrong-map controls
- `artifacts/` — quarantined failed Jurkat branch, retained for reproducibility
