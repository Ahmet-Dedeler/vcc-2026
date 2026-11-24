# 370 — Gene-specific Jurkat sign prior hits 0.577 held-out on A's block vs 0.515 all-down; PASS, A-only

Status: **COMPLETE** (kill criterion passes; production prior file written; no board claim — pds price unmeasured here).

## Known-answer reproduction (own code, `run370.py`, before any new number)

| quantity | reproduced | recorded | match |
|---|---|---:|---|
| context_A shape | 18400 × 18533 | AGENTS.md §2 / challenge bundle | exact |
| context_A median library | 20,109 | AGENTS.md §2 (~20,109) | exact to the integer |
| gene order head | TSPAN6, TNMD, DPM1 | `gene_names.csv` + anndata `var_names` | exact |
| 323 cap15 nmae (mean of `per-target-cap15.csv`, 58 rows) | 0.95083646 | `metrics-cap15.json` 0.95083646 | diff < 1e-9 |
| 323 cap15 fid / reach means | 0.49755452 / 0.25434726 (85 / 82 rows, 3 reach NaNs) | 0.49755452 / 0.25434726 | exact to 8 dp |
| Jurkat own-target delta < 0 | 2045/2270 covered, 0.99707 negative | mandate (99%) | confirms sign convention |
| K562 / H1 own-target delta < 0 | 7617/9787, 0.99028 / 255/255, 1.0000 | — (new checks) | convention holds both lines |

## What was measured (inputs, estimator, controls)

- **Block for A:** mean CPM per gene from `context_A.h5ad` controls (h5py over
  `X/data|indices|indptr`, 2,000-row chunks); mapped by name to the Nadig Jurkat
  donor axis (`delta.npy` 2270×8881 mmap, `genes.txt`, `b_ctrl.npy`); block =
  top-2,200 by A-control CPM among mapped genes (cutoff 130.4 CPM). Per-target
  eval set = block minus that target's head (top-400 |delta| on the donor axis;
  overlap mean 91.6/target, range 39–203; 6 exact-zero bank entries dropped).
- **Coverage:** 8,284/18,533 challenge genes map to the donor; 597/8,881 donor
  genes absent from challenge; 2,195 of the challenge-wide top-2,200 CPM genes
  are in the mapped block.
- **Split-half prior:** 5 repeats, `default_rng(20260910)` sequential
  permutations, 1135 fit / 1135 eval targets. Fit per-gene `P(down)` = fraction
  with delta < 0, plus magnitude-weighted `Pw` (zero-denom → 0.5; Jurkat: none).
  Eval on held-out half block entries: (i) all-down; (ii) gene prior
  `P>0.5→down`; (iii) thresholded `|P−0.5|>τ→prior else all-down`, τ ∈
  {0.05,0.10,0.15,0.20} + entry/gene coverage; (iv) shuffled-gene control
  (sub-RNG 9000+rep). Unweighted accuracy (what fid reads) and |delta|-weighted
  accuracy (what nmae-on-gate reads), mean ± sample-sd over splits.
- **H1 analogue:** same protocol on `h1_2025` donor (255×18080; genes verified
  identical to 323 ref axis), block = top-2,200 `ref_lin`, 127/128 halves.
  ~370 zero-denom genes/half → Pw=0.5 (affects `priorW`-H1 rows only).
- **Cross-line:** full-data priors both lines; shared axis 7,352 genes; eval
  blocks = top-2,200 own-`b_ctrl` among shared, head-excluded, all targets.

## Results tables (units: fractions; ± = sample-sd over 5 target-splits)

### Jurkat A-block, held-out (2.39M eval entries/repeat)

| estimator | unweighted acc | \|\|delta\|\|-weighted acc | coverage |
|---|---:|---:|---:|
| all-down | 0.5147 ± 0.0013 | 0.5244 ± 0.0024 | — |
| **gene prior (full)** | **0.5766 ± 0.0013** | **0.6269 ± 0.0020** | 1.000 |
| magnitude-weighted prior | 0.5754 ± 0.0014 | 0.6269 ± 0.0021 | 1.000 |
| prior, τ=0.05 | 0.5675 ± 0.0015 | 0.6133 ± 0.0023 | 0.635 ± 0.006 |
| prior, τ=0.10 | 0.5455 ± 0.0014 | 0.5775 ± 0.0023 | 0.330 ± 0.007 |
| prior, τ=0.15 | 0.5264 ± 0.0012 | 0.5449 ± 0.0023 | 0.150 ± 0.006 |
| prior, τ=0.20 | 0.5178 ± 0.0012 | 0.5298 ± 0.0023 | 0.057 ± 0.005 |
| shuffled-gene control | 0.4986 ± 0.0010 | 0.4978 ± 0.0014 | — |

Full prior wins at every τ; thresholding strictly loses (weak genes carry real
signal). Magnitude weighting adds nothing (±0.001). Full-data `P(down)` on the
block: 45.3% of genes up-leaning, deciles
0.209/0.391/0.427/0.456/0.484/0.512/0.537/0.568/0.603/0.660/0.825.

### H1 analogue (0.275M eval entries/repeat)

| estimator | unweighted acc | \|\|delta\|\|-weighted acc |
|---|---:|---:|
| all-down | 0.5171 ± 0.0034 | 0.5404 ± 0.0083 |
| **gene prior (full)** | **0.5558 ± 0.0010** | **0.5909 ± 0.0039** |
| shuffled | 0.5019 ± 0.0013 | 0.5055 ± 0.0029 |

Method works on a second line with 9× fewer fit targets (+0.039 over all-down).
Transport check: H1 donor-bank all-down 0.5171 vs 361's *delivered* all-down
tilt 0.5178 — bank sign statistics transport to delivered signs near-exactly.

### Cross-line (full-data priors, all targets)

| direction | all-down | cross prior | cross prior (weighted) | P-vector Pearson / sign agreement |
|---|---:|---:|---:|---|
| Jurkat prior → K562 block (21.0M entries) | 0.5091 | 0.5164 (+0.007) | 0.5357 | r = 0.47 / 0.6405 |
| K562 prior → Jurkat block (4.8M entries) | 0.5374 | 0.5478 (+0.010) | 0.5807 | (same shared axis) |

Partly generic (r = 0.47), mostly line-specific: within-Jurkat gain +0.062 vs
cross gain +0.007–0.010. K562's prior on Jurkat (0.5478) loses 0.03 to Jurkat's
own (0.5766).

### Block-definition note (full-data, in-sample, comparison only)

A-CPM block vs donor-`b_ctrl` block (both among 8,284 mapped): overlap only
1,390/2,200. All-down: 0.5146 (A-CPM) vs 0.5334 (`b_ctrl`). The PI's "54%"
was donor-basal-defined; **on the actual A-control block the tilt is 51.5%** —
the prior (+0.062) beats both baselines regardless.

## Verdict against the PLAN.md kill criterion

**PASS.** Criterion: held-out block accuracy ≥ 0.56 unweighted vs all-down
~0.54. Measured **0.5766 ± 0.0013 vs 0.5147 ± 0.0013** — clears by 0.0166
(~13 split-sd), beats shuffled by 0.078. For H1 the method reaches 0.5558 vs
all-down ~0.517, i.e. ~4.8× 361's +0.008 fid-equivalent tilt gain in accuracy
units ((0.5558−0.5171)/(0.5178−0.5083) ≈ 4.1× the all-down mechanism).
**Do not convert these accuracies to score points**: they are donor-delta sign
accuracies, not delivered-lfc fid/nmae/reach deltas; the pds/mse price of
flipping ~45% of block signs is unpriced here (363-style 323 arm needed).
Recommendation to PI: the prior in `work/jurkat_block_prior.npz` (all 2,270
targets) is the strongest A-side sign instrument measured to date; price it as
a design-stage sign-set arm on A's block with pds/mse read through emission.
**A-only**: B/C have no donor bank in `data/donors/` (HeLa/CAL33 absent); the
cross-line result forbids borrowing K562's or H1's prior for A, let alone B/C.

## Limits

- One study per line (Nadig Jurkat); splits are within-study target splits, no
  cross-study interval. Donor protocol ≠ challenge protocol (transport evidence
  is the H1 tilt match, a mean — not gene-specific accuracy).
- Head proxy: donor-|delta| top-400 stands in for the emitter's design head;
  delivered call-set overlap unmeasured here.
- H1 `priorW` rows carry ~2% zero-denom genes defaulted to 0.5 (unweighted
  `prior` rows unaffected; Jurkat/K562/full-file priors have no NaNs).
- Weighted accuracy uses donor |delta|, not reference |lfc| — directionally
  relevant to nmae, not identical to it.
- First script version had a `prior_w` key collision (weighted-acc-of-P
  overwritten by unw-acc-of-Pw) and an H1 divide-by-zero warning; both fixed,
  full rerun reproduced the uncontaminated keys bit-identically on `prior`.

## Exact commands to reproduce

```
cd ~/Code/virtual-cell-challenge-2026
memory_pressure | tail -1   # proceed only if free >= 50%
uv run python experiments/370-jurkat-block-sign-prior/run370.py
```

Writes `work/metrics370.json`, `work/per-split370.csv`,
`work/jurkat_block_prior.npz` (+ intermediates `contextA_mean_cpm.npy`,
`contextA_lib.npy`, `rerun.log`, `rerun.time`). Reads only: challenge +
donor banks + 323 `reference.npz`/`metrics-cap15.json`/`per-target-cap15.csv`.

## Peak memory and wall time

`/usr/bin/time -l` full rerun: **peak RSS 4.67 GB** (cap 12 GiB), **wall 6.9 s
compute** (+ interpreter/uv startup), single process, no run exceeded 4 GiB at
plan time estimates and none exceeded the cap; memory free ≥ 71% throughout.
Deterministic (seed 20260910): rerun reproduced `prior` 0.5766 / H1 0.5558
exactly.
