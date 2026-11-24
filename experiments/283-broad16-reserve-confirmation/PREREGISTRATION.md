# 283 — Preregistration: claim C2, F5 on the four RESERVED Broad-16 lines

Written **2026-09-06T08:46:12Z**, before any byte of the reserved lines' data was fetched, read or computed on.
Nothing in this file was informed by a reserved-line number. Every reference number below is quoted from 268's
working-line tables or recomputed here from 268's working-line output files (`work/prereg-calibration.txt`).

The reserve: four Broad-16 lines `JVE127`, `MG63`, `SLR23`, `UMRC3`, declared by
`experiments/268-broad16-context-residuals/work/reserve.json` (sha256
`1ce639b509d353708feded808914c63233bcee03513ace1c3f47cb93d071b9ce`) at 2026-09-05T16:20Z by a fixed name-sort rule,
before any response vector existed. Their single-cell archive members were never downloaded; their pseudobulk
columns were asserted unread by `load268.py`. `docs/RESEARCH-PROGRAMME-2026-09.md` rule 2: **one access, this claim
only, appended to `docs/CONFIRMATION-LEDGER.md`.**

## What is being confirmed

F5 (`docs/FINDINGS-2026-09.md`, from 268 on the **12 working** lines): a reproducible line-specific response component
exists (11.5 % of full-depth response energy pooled, 4× the OR-control null's 2.8 %, 79 % surviving the guide unit),
it is high-rank (27 components for half its reproducible energy), and it is **not** basal-state-predictable under
leave-one-line-out (the basal-conditioned model's predicted residual has cosine +0.018 with the held residual against
+0.001 for the wrong recipient, but adding it to the shared response buys nothing: Δcos −0.0071 / +0.0000, Δretrieval
+0.0003 / −0.0010).

## Data and pipeline (frozen)

- **Source**: Figshare article 33273600, `single_cell_data.zip` (`https://ndownloader.figshare.com/files/67643568`,
  3,171,208,158 B), members `<line>_*GE_filtered_feature_bc_matrix` for the four reserved lines, fetched by HTTP range
  request and CRC32/size-verified into `work/sc/` (never under `data/`). Cell assignment from the authors' already-present
  `experiments/268-broad16-context-residuals/work/proc/all_cells_table.csv`.
- **Build**: `work/build_halves2.py` copied byte-identical from 268, `EXCLUDE_ARM=1`, seed **268** — the authors'
  `pass_qc` cells with neither `arm_trunc` nor `arm_gain`, per-KO pseudobulk sums for FULL / random cell halves A,B /
  guides v1,v2, duplicate gene symbols summed. Parse known-answer per line: FULL must equal the on-disk
  `pseudobulk_sum_by_ko` column exactly.
- **Targets**: 268's frozen 82 usable targets (≥ 6 cells on every working line), further restricted per reserved line to
  those with ≥ 6 cells in that line. The count kept is reported per line.
- **Gene axis, per reserved line**: genes non-NaN in all 12 working lines and in that reserved line, the 100 Broad
  target genes removed, mean control CPM ≥ 5 **computed on the 12 working lines only** — the leave-one-line-out
  analogue of 268 §4.
- **Response**: `log1p(5e4·p/Σp)` minus the same-half control log-profile (the 10 OR knockouts summed). The
  **shared response** is the mean over the 12 working lines of the same half. All fits (gene filter, centring, basal
  PCA, ridge, residual programme basis) use working lines only; the reserved line enters only as the thing predicted.
- **Bootstrap**: paired **target** bootstrap, 4,000 draws, seed **283**; each target carries its reserved folds.
  95 % = the 2.5/97.5 percentiles.

## Known-answer rows (run and reported BEFORE the reserved data is used)

- **KA-1**: `build_halves2.py` re-run on one working line (`T3M5`, batch apollo34 — the batch of reserved `JVE127` and
  `MG63`) from the already-fetched h5 must reproduce 268's `work/halves/T3M5.npz` at **0.0e+00** on every array.
- **KA-2**: the analysis script run in leave-one-working-line-out mode over the 12 working lines must reproduce every
  deterministic column of 268's `work/lolo268_K2.csv` (`ceiling_cos_AB`, `base_cos_full`, `base_retrieval`, and the
  `a`/`b`/`c1`/`c2`/`nn`/`oracle`/`WRONG` arms) at **0.0e+00**. The `SHUFFLED_labels` arms consume a shared RNG stream
  across folds and are reproduced only if the whole 12-fold loop is replayed; that is what is run.

An abort on either KA means no reserved-line number is computed (`AGENTS.md` §known-answer rule).

## Claim C2 — the frozen sub-claims, statistics and pass rules

Let `r_A[t] = R_A[t, reserved] − mean_{12 working} R_A[t, ·]` and likewise for half B (each half's shared mean is built
from the working lines' **same** half, so the two residuals share no cell on either side — 268's `lolo268.py`
construction).

### C2a — the line-specific residual is reproducible on the reserved lines

- **Statistic**: per reserved line, the median over targets of `cos(r_A[t], r_B[t])`, with a paired-target bootstrap
  interval (4,000 draws, seed 283).
- **PASS rule (frozen, decides C2a)**: on **≥ 3 of the 4** reserved lines, median > **0.028** *and* the bootstrap 95 %
  lower bound > **0.028**. Otherwise FAIL.
- **Calibration of this rule, disclosed now** (`work/prereg-calibration.txt`, computed from 268's working-line
  `lolo268_K2.csv`, no reserved data): the identical rule applied to 268's **12 working** lines passes on **8 of 12**
  (`C4I` 0.0211, `KMRC20` 0.0315 lower 0.0182, `SNU8` 0.0231, `T3M5` 0.0366 lower 0.0199 fail). The working-line
  per-fold medians span **[0.0211, 0.0652]**, overall median 0.0455. So ≥ 3/4 is a *stricter* bar than the working
  lines themselves clear, and a FAIL of the letter of this rule would not by itself refute F5 — which is exactly why
  the two calibrated secondaries below are frozen at the same time.
- **Secondary S1 (frozen, reported whatever it is)**: whether each reserved line's median lies inside the working-line
  per-fold range **[0.0211, 0.0652]**, and whether the number of reserved lines passing the 0.028 rule is consistent
  with the working lines' 8/12 (binomial, 4 draws at p = 8/12: P(≥3) = 0.593, P(≥2) = 0.889, P(0) = 0.012).
- **Secondary S2 (frozen)**: per reserved line, the fraction of targets whose residual split-half cosine exceeds that
  **same line's own** permutation null (half-A residual of target *t* paired with half-B residual of a permuted target,
  same line) at its 95th percentile. Null expectation 0.05; 268 measured 0.262 over the 984 working (target, line)
  pairs, and 0.075 on the OR controls. The same statistic is computed per working line in the KA run for comparison.
- **Secondary S3 (frozen)**: per reserved line, the same median residual cosine computed on the 10 OR-control knockouts
  (each against the other nine), the within-line target-independent null. 268 working-line value: median 0.0183.
- **Secondary S4 (frozen)**: the same median residual cosine with **guide halves** (`v1` vs `v2` cells, each against its
  own guide's controls) instead of cell halves — F5's "79 % survives the guide unit" claim. 268 working-line value:
  cell-half cosine median 0.0461 vs guide-half 0.0355 (0.77×).

### C2b — the 12-line shared response transfers to the reserved lines

- **Statistic**: per reserved line, the median over targets of `cos(shared response, reserved line's full response)`,
  paired-target bootstrap interval.
- **The range, stated before looking**: 268's working-line leave-one-line-out per-fold medians of exactly this quantity
  are **[0.0400, 0.0937]** (`SNU8` low, `HKGZCC` high; overall median 0.0658; per-fold means 0.0741–0.1450).
- **PASS rule (frozen)**: on **≥ 3 of the 4** reserved lines the median lies inside **[0.0400, 0.0937]**. FAIL otherwise.
- Declared asymmetry: the reserved lines' shared response is a mean of **12** lines where the working folds used **11**,
  which if anything helps the reserved lines slightly; the direction of that bias is stated, not corrected.
- Secondary (frozen): the same median for the **retrieval** statistic 268 used (1 − rank of the held truth among the
  82 by cosine to the prediction, `lolo268.py:retrieval`); working-line per-fold means 0.700–0.820.

### C2c — basal state does not predict the reserved lines' residuals

Model (c), fitted on the **12 working lines only** exactly as `lolo268.py` fits it on 11: basal programmes = the first
`K = 2` principal components of the working lines' centred control log-profiles, standardised; (c1) per-gene ridge
(λ = 1) of each target's interaction residual on the basal scores; (c2) the same at the level of a `J = 5` residual
programme basis (SVD of the working lines' interaction residual). Baseline (a) = the shared response alone (`r̂ = 0`).
Controls: **WRONG recipient** (the same fitted model driven by each working line's basal scores, mean over the 12) and
**SHUFFLED labels** (basal scores permuted across working lines before fitting, 5 permutations).

- **Primary statistics**: pooled over the four reserved lines, per-target paired differences (c) − (a) in
  **Δcos(full response)** and **Δretrieval**, with paired-target bootstrap 95 % intervals.
- **PASS rule for "no better" (frozen)**: C2c **PASSES** — i.e. 268's negative is confirmed — if, for **both** c1 and c2,
  the 95 % interval of Δcos(full) **and** of Δretrieval either covers 0 or lies entirely below 0. C2c **FAILS** — the
  negative does not reproduce, basal state does predict usefully on the reserve — if for either arm either interval
  lies **entirely above 0**.
- **Secondary (frozen, the "detectable but not usable" half of F5)**: cos(predicted residual, reserved residual) for
  c1/c2 against their WRONG-recipient and SHUFFLED-label controls. 268's working-line values: c1 +0.0181
  [+0.0134, +0.0229] vs wrong +0.0012, shuffled −0.0073; c2 +0.0174 [+0.0110, +0.0240] vs wrong +0.0007, shuffled
  −0.0007. Expected direction on the reserve: positive and above both controls. This is reported whatever it is and is
  **not** part of any pass rule.

## Scope of the access

One access. The four lines' single-cell members are fetched once, built once, and used for exactly the statistics
above. No model, no build, no `vcc`, no submission, no further use of the reserved lines, no selection of any kind on
them. Every number produced is written to `RESULTS.md` and the ledger, pass or fail.
