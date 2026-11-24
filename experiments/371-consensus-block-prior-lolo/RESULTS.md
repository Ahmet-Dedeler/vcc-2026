# 371 — Consensus prior transfers at +0.012 median over all-down and loses to the best single bank on 5/7 held-outs; NEGATIVE, no B/C rider

Status: **NEGATIVE** (kill criterion fails; `work/consensus_prior.npz` written for
the PI regardless, per mandate; no board claim — donor-delta sign accuracies,
not delivered score deltas).

## Known-answer reproduction (own code, before any new number)

| quantity | reproduced | recorded | match |
|---|---|---|---|
| 323 cap15 nmae (mean of `per-target-cap15.csv`, 58 rows) | 0.95083646 | `metrics-cap15.json` 0.95083646 | diff < 1e-9 |
| 323 cap15 fid / reach means | 0.49755452 / 0.25434726 (85 / 82 rows) | same | diff < 1e-9 |
| context_A/B/C shape | 18400 × 18533 all three | AGENTS.md §2 | exact |
| context_A median library | 20109 | 370 / AGENTS.md §2 | exact to the integer |
| challenge gene order head | TSPAN6, TNMD, DPM1 | `gene_names.csv` | exact |
| Jurkat own-target delta < 0 | 2045/2270 covered, 0.99707 | 370's 0.99707 | within 0.002 |
| 370 A-block (top-2200 A-CPM among Jurkat-mapped) | size 2200, min CPM 130.4 | 370's 130.4 | exact to 0.1 |
| K562-gwps own-target frac<0 0.99028; H1 255/255, 1.0000 | — | 370's checks | consistent |

Full rerun under `/usr/bin/time -l` reproduces `metrics371.json` **md5-identical**
(`87508ac7…`), so every number below is deterministic.

## What was measured (inputs, estimator, controls)

- **Banks (7):** `k562_gwps` (9787×8248), `k562_essential` (2021×8563),
  `nadig_jurkat` (2270×8881), `h1_2025` (255×18080),
  `194-debias-orion_hct116` (17767×38606), `194-debias-orion_hek293t`
  (17977×38606), `nadig_hepg2` (2160×9623, same 5-file layout verified).
  All `delta.npy` mmap'd, streamed in 1000-row chunks; sign convention
  pert-minus-control verified per bank via own-target negativity (table below).
- **Shared axis:** 7-way intersection by gene name, sorted, **6346 genes**
  (6-way without hepg2 would be 6683 — hepg2 costs 337 genes / 5.0%).
  All 6346 are on the challenge axis (34.2% coverage).
- **LOLO:** held-out bank L defines block_L = top-2200 by its own `b_ctrl` among
  shared genes; per-target head = top-400 |delta| on the bank's FULL donor axis
  (370's method), intersected with the block; exact-zero deltas dropped
  (<0.04% of block entries everywhere — see zero-fraction table).
  Truth = sign of L's delta. Priors fitted WITHOUT L: (a) all-down;
  (b) each other bank's own P(down); (c) equal-weight consensus mean;
  (d) target-count-weighted consensus; (e) agree-only consensus (all 6 others
  same side of 0.5, else all-down); (f) shuffled-gene consensus (5 draws).
  Ties (P==0.5) → predict up; **zero ties on any block**.
- **Endpoints:** unweighted accuracy (fid reads this) and |delta|-weighted
  accuracy (nmae-on-gate proxy), per bank and pooled; 1,000-draw target
  bootstrap per bank (seeded, per-target sums precomputed) for the
  consensus-equal-minus-all-down difference.
- **A/B/C:** chunked-h5py mean CPM for all three contexts (same code as 370);
  consensus coverage of each context's challenge-wide top-2200 block; A sanity
  scores full-consensus (in-sample for Jurkat — flagged) and LOLO-consensus
  (clean) on 370's exact A-block with Jurkat deltas, vs 370's 0.5766.

## Results tables (units: fractions)

### LOLO main endpoints (21.0M–38.9M eval entries/bank except H1)

| held-out L | all-down unw / w | con_eq unw / w (gain) | con_w unw | agree unw (entry cov) | shuffled unw |
|---|---:|---:|---:|---:|---:|
| k562_gwps | 0.5077 / 0.5127 | **0.5198** / 0.5426 (**+0.0122**) | 0.5183 | 0.5109 (0.109) | 0.5010 |
| k562_essential | 0.5293 / 0.5400 | **0.5537** / 0.5921 (**+0.0244**) | 0.5461 | 0.5356 (0.116) | 0.5043 |
| nadig_jurkat | 0.5313 / 0.5522 | **0.5517** / 0.5892 (**+0.0204**) | 0.5353 | 0.5369 (0.121) | 0.5051 |
| h1_2025 | 0.5159 / 0.5368 | **0.5227** / 0.5475 (**+0.0067**) | 0.5199 | 0.5187 (0.143) | 0.5017 |
| orion_hct116 | 0.5051 / 0.5148 | **0.5061** / 0.5163 (**+0.0010**) | 0.5059 | 0.5059 (0.124) | 0.5004 |
| orion_hek293t | 0.4993 / 0.5117 | **0.4952** / 0.4965 (**−0.0041**) | 0.4980 | 0.4999 (0.233) | 0.5002 |
| nadig_hepg2 | 0.5429 / 0.5743 | **0.5591** / 0.6041 (**+0.0163**) | 0.5517 | 0.5518 (0.109) | 0.5010 |

Pooled (micro over 113.5M entries): all-down 0.5072 → con_eq **0.5109 (+0.0037)**;
weighted 0.5184 → 0.5241 (+0.0057). Macro mean gain +0.0110,
**median gain +0.0122**. Target-count weighting (d) never beats equal weighting
(c) except trivially on hct116 (0.5059 vs 0.5061 — still loses): the two Orion
banks dominate the weights (~70%) and drag the mean toward coin-flip.

### Bootstrap: consensus-minus-all-down (1,000 target resamples; difference units)

| held-out | d_unw mean [95% CI] | P(d>0) | d_w mean [95% CI] |
|---|---:|---:|---:|
| k562_gwps | +0.0122 [+0.0111,+0.0133] | 1.000 | +0.0299 [+0.0274,+0.0325] |
| k562_essential | +0.0244 [+0.0212,+0.0276] | 1.000 | +0.0521 [+0.0464,+0.0575] |
| nadig_jurkat | +0.0204 [+0.0180,+0.0229] | 1.000 | +0.0370 [+0.0334,+0.0410] |
| h1_2025 | +0.0067 [+0.0034,+0.0099] | 1.000 | +0.0107 [+0.0030,+0.0188] |
| orion_hct116 | +0.0010 [+0.0004,+0.0014] | 1.000 | +0.0015 [+0.0004,+0.0025] |
| orion_hek293t | **−0.0041 [−0.0048,−0.0034]** | 0.000 | **−0.0152 [−0.0167,−0.0136]** |
| nadig_hepg2 | +0.0163 [+0.0137,+0.0187] | 1.000 | +0.0299 [+0.0255,+0.0342] |

Every interval excludes zero — this is a resolved measurement, not noise. The
hek293t interval is resolved **negative**: the consensus rider would cost
accuracy on a hek293t-like line, 3.7× worse weighted than unweighted.

### Consensus vs best single bank (control (b); unweighted)

| held-out | con_eq | best single donor | delta (con − best) |
|---|---:|---:|---:|
| k562_gwps | 0.5198 | **k562_essential 0.5211** | −0.0013 |
| k562_essential | 0.5537 | **k562_gwps 0.5560** | −0.0023 |
| nadig_jurkat | 0.5517 | **k562_essential 0.5539** | −0.0022 |
| h1_2025 | **0.5227** | k562_essential 0.5209 | +0.0018 |
| orion_hct116 | 0.5061 | nadig_hepg2 0.5062 | −0.0001 |
| orion_hek293t | 0.4952 | **orion_hct116 0.5028** | −0.0076 |
| nadig_hepg2 | **0.5591** | k562_essential 0.5537 | +0.0054 |

Consensus beats the best single on **2/7** banks, ties 1, loses 4 — three of
them to a same-line-different-study donor (K562↔K562-essential both directions;
k562_essential is the best donor on 4/7 held-outs). Pooling adds nothing over
picking the K562-essential prior; per-bank oracle-best gains
(0.0134/0.0267/0.0226/0.0068/0.0011/0.0035/0.0162) have median **0.0134**,
still below the +0.02 gate. Debiased-Orion donors transfer worst everywhere:
hek293t's prior scores 0.4723 on Jurkat (6 points BELOW all-down) and below
all-down on 4/6 recipients; hct116's below all-down on 2/6.

### Bank QC (zero fractions exact; own-target sign check)

| bank | zero full-axis | zero shared-axis | P_down mean | own-target cov | own frac<0 |
|---|---:|---:|---:|---:|---:|
| k562_gwps | 0.0000 | 0.0000 | 0.505 | 7617/9787 | 0.9903 |
| k562_essential | 0.0000 | 0.0000 | 0.500 | 1835/2021 | 0.9967 |
| nadig_jurkat | 0.0000 | 0.0000 | 0.498 | 2045/2270 | 0.9971 |
| h1_2025 | 0.0548 | 0.0000 | 0.503 | 255/255 | 1.0000 |
| **orion_hct116 (debiased)** | **0.0873** | 0.0003 | 0.503 | 17767/17767 | **0.9100** |
| **orion_hek293t (debiased)** | **0.0487** | 0.0003 | 0.504 | 17977/17977 | **0.9098** |
| nadig_hepg2 | 0.0000 | 0.0000 | **0.534** | 1958/2160 | 0.9939 |

Mandate's Orion check: the debiased/thinned banks carry 4.9–8.7% exact zeros
full-axis (concentrated off the shared axis — only 0.03% on shared genes), and
their own-target knockdown negativity is only **0.91 vs 0.99–1.00** everywhere
else: the debias representation weakens even the one sign that should be
certain. hepg2 has the strongest down tilt of any bank (P_down 0.534; block
all-down 0.5429/0.5743).

### A/B/C coverage + A sanity (units: genes; fractions)

| context | top-2200 block in consensus axis | min block CPM |
|---|---:|---:|
| A (Jurkat) | 1920/2200 (87.3%) | 130.6 |
| B (HeLa) | 1739/2200 (79.0%) | 125.1 |
| C (CAL33) | 1686/2200 (76.6%) | 125.9 |

A-block (370's exact 2200, Jurkat deltas, 1924 covered): full-consensus
(in-sample) 0.5602 covered-acc (w 0.6014) vs same-line split-half **0.5766**
(loses 0.016 despite the in-sample advantage); LOLO-consensus (clean transfer)
**0.5456** (w 0.5775) — below 370's K562-single 0.5478 and 0.031 below
same-line. With uncovered genes as all-down fallback: 0.5388 (LOLO) /
0.5517 (full) vs all-down 0.5146. Median libs: A 20109, B 19946, C 20034.

## Verdict against the PLAN.md kill criterion

**FAIL — record the negative.** Criterion: consensus LOLO gain ≥ +0.02
unweighted on the median held-out bank with bootstrap CI above zero. Measured:
median gain **+0.0122** (1.6× short); only 2/7 banks individually clear +0.02
(k562_essential +0.0244, jurkat +0.0204 — the two whose recipients have the
highest all-down tilt and a same-line sibling in the pool). The CI-above-zero
half passes on 6/7 banks, so this is not underpowering — the effect is
resolved and too small. Secondary controls agree: consensus ≤ best-single on
5/7; agree-only restriction loses to full consensus everywhere (weak genes
carry signal — same lesson as 370's τ sweep); shuffled sits at 0.500–0.505
(gene-specificity is real, just small). **No B/C rider.** Side note for the PI:
if a transfer prior were ever forced, the instrument says k562_essential's
P(down), not the consensus — but its median-bank gain (+0.0134 oracle-best)
also fails the gate, so that is a description, not a proposal.

## Limits

- LOLO banks are whole studies/lines; splits are within-bank target bootstraps
  — no cross-study interval for the same line (except the K562 pair, which is
  the strongest result here).
- Donor-delta sign accuracies, not delivered fid/nmae/reach deltas; the
  pds/mse price of flipping ~45% of block signs is unpriced (needs a 363-style
  323 arm, deliberately out of scope).
- Weighted accuracy uses donor |delta|, not reference |lfc| — nmae-relevant,
  not identical.
- Heads are full-axis top-400 |delta| (370's method): head density in-block
  varies by bank axis size (mean overlap 34.8–86.8/target); a shared-axis head
  would change eval sets slightly, not the ranking (all estimators share it).
- B/C tilts are unknowable without perturbation data; the hek293t row shows a
  no-tilt line can invert the consensus.
- H1 (255 targets) is noisy as a recipient but its CI still excludes zero;
  nothing rests on it alone.

## Exact commands to reproduce

```
cd ~/Code/virtual-cell-challenge-2026
memory_pressure | tail -1   # proceed only if free >= 50%
uv run python experiments/371-consensus-block-prior-lolo/run371.py
```

Reads only: 7 donor banks + challenge h5ads + 323 `metrics-cap15.json` /
`per-target-cap15.csv`. Writes only `work/`: `consensus_prior.npz` (genes,
P_down equal-weight all-7, P_down_weighted, P_banks 7×6346, bank_names,
n_targets, n_banks=7, seed), `metrics371.json`, `per-bank371.csv`,
`context{ A,B,C}_mean_cpm.npy`, `run371*.log`. Rerun is md5-identical on
`metrics371.json`.

## Peak memory and wall time

`/usr/bin/time -l` full run: **peak RSS 4.79 GB** (cap 12 GiB),
**wall 23.8 s** (22.2 s user), single process, mmap + 1000-row chunked
streaming; memory free ≥ 73% throughout. Well inside the 2.5 h budget.
