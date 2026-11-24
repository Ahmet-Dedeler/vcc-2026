# 810 PROPOSAL: strength-gated unwhitening of the donor base (P1), and its agreement-gated form (P2)

Written 2026-10-16 ~08:15Z by the 810 worker (lane B1). No candidate has been priced. Every number below is from the
measurement step (RESULTS.md §1–§2: the incumbent's composition split into its parts on 22 folds, plus composition-level,
direction-only ablations). The pricing code `price810.py` (sha256 in RESULTS.md) exists but has not run.

## 1. Mechanism (falsifiable, member by member)
- **Where the misdirection is.** On the top z_src tertile of BULK-W, the incumbent's direction is limited by the WHITENED
  donor base. On those rows, the unwhitened, normrel'd donor mean `pre` (the input the whitener receives) points at the
  held line's truth about 1.7× better than the whitened base: median cos 0.0329 vs 0.0193.
  - BULK-ABC cancer, strong tier: 0.0459 vs 0.0334.
- **What it is not.** H1 and CD4 do not misdirect the strong tier: dropping either part LOSES on 18/18 BULK-W folds and
  3/3 BULK-ABC folds. 706's hypothesis ("H1/CD4 dilute the strong rows") is refuted.
- **Why whitening hurts here.** γ 0.65 shrinks the directions with high signal variance across the donor screens. A strong
  response that two lines share lies along exactly those directions. The whitening loss within the tier grows with source
  strength (Spearman +0.22) and donor agreement (+0.16); on BULK-ABC cancer, z_src +0.23 and agreement +0.14. On BULK-W
  the gain is carried by strong targets whose two sources agree: +0.0055 fold-X vs +0.0005 on the strong but disagreeing
  half, and ≈ 0 on the weak and mid tiers (`base → pre`'s predecessor, the pre-only swap).
- **Why gating can beat the global knife-edge (212).** The same unwhitening applied to every target loses on the board,
  because pds pays what mse gains:
  - 097 closed whitening relaxation;
  - 212 found cosine maximised at γ 0 and pds at γ 0.8, cancelling to 3.8 %;
  - 099's SELF control saw pds flip sign across surfaces.

  The two members weight targets differently. mse reads the aggregate cosine, which is ENERGY-weighted: the strong tier
  holds 59 % (A/B/C) to 73 % (BULK-W) of the realised energy a. pds is COUNT-weighted: each target counts 1/n.
  Unwhitening only the top tertile keeps 47 % (A/B/C) to 80 % (BULK-W) of the global unwhitening's composition-level
  cosine gain, while exposing a third of the targets to pds risk. Those targets already have the largest pds margins
  (literal pds 0.944 on A/B/C, 0.863 on BULK-W).
- **Members.**
  - mse improves through the aggregate cosine (rows keep their norms, so the energy A and the cliff do not move, unlike 170).
  - pds loses some specificity on the moved rows.
  - fid / jac / reach / nmae are pinned by the unchanged carriers. Ballast K moves only through the render; the engine
    reports dK, and a |dK| > 5 % sends the candidate to DE-H1.
- **The construction keeps H1 and CD4.** Only the base part is re-directed, then the row is rescaled to its own norm. That is
  the difference from 707, which blends the whole final row toward U and so also thins the H1 and CD4 parts. kolf2's sign
  separates the two (composition level): 707's U blend at β 0.25 is −3.9 %, the base→pre swap at β 1 is +3.1 %.

## 2. Class
**Transfer.** A better bulk direction for the strong targets from the same public donors and operator. Nothing is fitted
to A/B/C, and the gate (z_src tertile, donor agreement) is computed from the day's donor banks.

## 3. Candidates and try budget (3 configurations, one family "810-sgu")
| arm | gate | β (base → pre) |
|---|---|---|
| sgu05 (P1) | top z_src tertile (707's tiering: over targets with a donor row) | 0.5 |
| sgu1 (P1) | same | 1 |
| sga1 (P2) | top tertile AND donor agreement > the tertile's median | 1 |

- Comparators, not tries of this family: the incumbent (L folds pct+P2g; ABC folds cap6+P2g(abc), as ka804) and 707's row
  (rebuilt; KA against 804's arrays).
- Best-of-3 shift e(3) = 0.85 sd on the winner.
- Nothing else is tried. In particular β is not tuned beyond {0.5, 1}, and the gate is not re-cut.

## 4. Prediction before data (check800 aggregate, bulk overall = mse(r) + literal pds; pdsl beside)
Composition → realised calibration from 707's own two reads (realised r / composition r):
- BULK-W: 7.97 / 11.27 = 0.71;
- BULK-ABC: 1.98 / 2.45 = 0.81.

**Projected realised r (relative aggregate cosine):**

| arm | BULK-W | BULK-ABC cancer | kolf2 |
|---|---|---|---|
| sgu05 | +15 % | +4.3 % | +1.7 % |
| sgu1 | +24 % | +7.0 % | +2.5 % |
| sga1 | +21 % | +6.1 % | — |

**pds.** 707's pds/r exchange on BULK-ABC is −0.0027 raw per +1.98 %, i.e. −0.00049 overall per % against mse's +0.000715
per %. P1's claim is that its exchange is no worse than 707's (it moves the same rows, keeps H1/CD4, and carries more
cosine per unit of rotation).

| arm | BULK-ABC cancer (3 folds): point (range) | bin | BULK-W (18 folds) | bin |
|---|---|---|---|---|
| sgu05 | +0.0009 (−0.0005 … +0.0025) | WEAK (3 folds: CI ≈ ±0.0013 as 707's) | +0.008 (+0.004 … +0.013), ≥ 16/18 | SHIP |
| sgu1 | +0.0015 (−0.0015 … +0.0035) | WEAK | +0.013 (+0.006 … +0.020), ≥ 16/18 | SHIP |
| sga1 | +0.0015 (−0.0005 … +0.0030) | WEAK, at best SHIP | +0.011 (+0.005 … +0.017), ≥ 16/18 | SHIP |

- Paired against 707 on BULK-ABC, P1's best arm should beat 707 by ≥ +0.0005, and its pds/r exchange should be ≥ 707's.
- BULK-W over-reads (key 4.75×; this direction class has no calibration point yet). The BULK-W numbers rank transfer; they
  are not board sizes.
- **Board projection:** none until 707's board read (10-17/18) gives this class its first ρ (702 rule 4: the interval
  includes zero until then).

**Falsifiers (any one closes the family for October):**
1. On BULK-ABC cancer, the pds leg ≥ the mse leg for all three arms (net ≤ 0): the energy-vs-count asymmetry does not hold
   on Arc's panel.
2. The per-moved-target pds cost of sgu1 ≥ the panel-average pds cost of global unwhitening. Read from the per-target pds
   arrays: Δpds summed over the tier's targets ÷ n_tier, against 212's global η. That would mean gating does not concentrate
   cosine relative to pds.
3. The best arm does not beat 707 on BULK-ABC (paired point ≤ 0): the base-only construction adds nothing over 707.
4. |dK_rel| > 5 % on any fold: the bulk engine cannot price it, and the candidate goes to DE-H1 first.

## 5. Novelty against the record (FINDINGS §2, the γ closures)
- **What is closed:**
  - a GLOBAL γ/λ change (170: γ is a local optimum for mse at α re-solve; 212: "do not reopen", cosine and pds cancel to
    3.8 %; 193: operator geometry);
  - global whitening relaxation (097; 099 §5.3 SELF);
  - per-target arm routing by a learned gate (066: PM in/out by stump);
  - 706/707's U blend of the final row (a separate family with a board read pending).
- **What is different here:**
  1. The gate: a third of the targets, chosen by the donors' own strength (or agreement), not a global γ. The mechanism
     (energy-weighted mse vs count-weighted pds) is a reason the global closure need not cover a gated form; 212 measured
     only the global knife-edge.
  2. Only the BASE is unwhitened. The H1/CD4/partial residuals are kept, and each row keeps its norm (no A/B energy
     change, so no cliff: 170's mechanism does not apply).
  3. It is predicted by a measured, per-target mechanism (whitening loss rising with strength and agreement), not tuned.
- 183's η table: un-whitening moves sat on a "bad part" of the pds–cos frontier (η −0.068). Falsifier 2 tests that
  frontier on the gated form directly.

## 6. Cost
- One EC2 run: 22 folds × 6 arms through engine801 (≈ 3 min per fold, 6 GB peak; parts already on S3 under
  `s3://…/output/lab/810parts/`).
- r7i.8xlarge, ~30 min, ≈ $1.5; ≤ $5 with one re-run. Aggregation with check800 on the Mac (MBs). Lane spend so far $0.73.

## 7. Day-of path (D/E/F)
- **Inputs.** The builder already computes `pre` (the normrel'd donor mean) and the whitened base before adding the
  residual. On the day, a dump hook (690g's interface, as 707's) recomposes the parts from the day's bundle, the day's
  3-donor operator and the day's CD4 banks with lib801's calls (decomp810.compose's code: 20–40 s, ≤ 6 GB).
- **KA.** The parts must sum to the builder's dump to ≤ 1e-12 (this lane measured 3.6e-14 on 22 folds).
- **Gate.** z_src from the three banks (522's def_strength_key522.zsrc, as 707 / 705); agreement from the three donors'
  raw rows (mean pairwise cosine).
- **The swap itself.** sgu_rows in price810.py. Nothing is fitted to A/B/C or D/E/F truth.
- **Unapplied** until a dossier and the coordinator's decision.
