# 286 — A genuinely different response estimator: shared response directions learned from CROSS-PRODUCTS of independent views (halves, guides, donors) with a full-rank residual path — tested on held views (CD4 donors and states), then on held contexts (PLOLO, three cancer folds), then on the direction-resolving Flex reference (273). Priced on the coherent-pessimistic rule. No build, no submission. ≤ 8 hours.

You are a research subagent. Work in `experiments/286-crossproduct-estimator/`. **Edit nothing outside it.** No build,
no `vcc`, no `cell-eval2`, nothing under `submissions/`, `data/harness/` or `data/`; no download. `uv run python`. One
fold per process, serial; peak RSS under ~10 GiB; `memory_pressure` ≥ 40 % free before each fold (hold otherwise — two
sibling experiments, 284 and 285, run in parallel; do not read or write their directories). **Time box: 8 hours.** If
Stage B cannot finish all three folds, report the folds done and stop. `docs/RESEARCH-PROGRAMME-2026-09.md` binds you
(freeze the estimator family and controls before any held number; name the unit; strong baselines; the already-failed
constructions are baselines here). Both reserves are SPENT: never read `experiments/267-*/work/reserve/` or the Broad-16
lines JVE127 / MG63 / SLR23 / UMRC3.

## Why (third external review, 2026-09-06; R271 §0, R272 §0.10 and §5)
271 tested SCALAR reliability weights and DIAGONAL per-gene shrinkage of the donor rows (0/10 cells). 272 re-weighted
the eigenvalues of the FITTED whitener's eigenbasis (0/11) and named the object it did not measure: "a noise correction
that re-estimates the eigenvectors … is a different object" (`experiments/272-noise-corrected-operator/RESULTS.md`
§5). 272's mechanism is the reason to expect that object to differ: on every operator on disk the fitted eigenvectors
are 77–79 % noise directions carrying 87–92 % of the whitened prediction energy, and `pds` NEEDS those directions
because the target-specific signal is spread over them. **A cross-product between two INDEPENDENT views of the same
perturbation has expectation Σ_signal — the sampling noise is independent across views and averages out — so the
directions of the shared response can be learned without the within-view noise inflating them, and each target's
response can then be estimated by a FULL-MATRIX combination (not a diagonal one), with a floor so that every direction
keeps a path (no forcing into a small set of programmes).** This could fail. It has not failed because diagonal
shrinkage failed.

## READ FIRST
`docs/READ-FIRST.md` §3 (pricing, uncertainty language), §4, §5; 271's `est271.py` (the notation, the plumbing:
`lib263.stage` → per-target combination BEFORE whitening, norm-matched `:nm` → `plolo.Operator.apply` → CD4 deflated
k = 5 re-solved at f = 0.34 → the H1 pair union at g = 0.28 → one add, one preserve → `lib263.evaluate`), `halfrows271.py`
(267's half-banks → production row units), `sweep271.py` (the template for a fold sweep), `price271.py` (261's pricing
recipe imported: `price261.py`; the BR249 retrodiction row), `analyse271.py`; 272's `op272.py` (`fit_operator` = 062's
whitener recipe; `DampedOperator`; the noise spectrum) and RESULTS §0.10, §5; 263's `lib263.py`; 267's half-banks
(`experiments/267-…/work/half/<line>/`: batch-parity halves A / B of EVERY bank perturbation — 9,787 / 17,767 / 17,977
rows — with their own controls) and RESULTS §1–§2; 162's `plolo.py` (`build_fold`, `Operator`, `solve_fraction`,
`deliver`); 265's `boot265.py` (paired target bootstrap); 273's `lib273.py` / `instrument273.py` (the union-OFF h1-test
instrument: base + CD4 re-extracted for the 85 targets; per-target `pds` and own cosine); 280's bank
(`experiments/280-…/work/bank280_X.npy`, `bank280_rows.npy`, `obs280.npz` keys `donor, cond, gtype, run, pert, guide,
n_cells, total, keep_*`; `extra280_X.npy` / `extra280_rows.npy` = 3,502 rows of 500 NON-panel Stim48hr targets) and
`q2_280.py` (leave-donors-out, the k = 2 baseline: Δ vs k = 1 +0.0063 full / +0.0043 matched; absolute k = 2 cosine
+0.0354 full, +0.0270 matched, Stim48hr).

## The frozen estimator family (write `est286.py` and `PREREGISTRATION.md`; hash both into `work/` before any held number)
For a bank with two independent views a, b of each row (each view centred on its own mean over the bank's rows;
units as production's rows), over n rows and g genes:
- **Ŝ** (signal covariance) = ½ (AᵀB + BᵀA) / n, projected to PSD (report the energy clipped); **N̂** (noise covariance)
  = (A − B)ᵀ(A − B) / (4n) for the mean of two halves (271's identity Var((A+B)/2) = Var(A−B)/4); **Ĉ** = ½ (AᵀA + BᵀB) / n.
  KA: ‖Ĉ − Ŝ − N̂‖_F / ‖Ĉ‖_F (state it; it is the moment identity up to the centring).
- **E0** the equal-weight mean of the donor rows = the incumbent (KA 0.0e+00).
- **E1** full-matrix Wiener on the mean row with a FLOOR: G = Ŝ (Ŝ + N̂ / J + ε I)⁻¹ symmetrised; in G's eigenbasis
  every gain is raised to ≥ φ (φ ∈ {0.10, 0.25}: the full-rank residual path — no direction is zeroed); x̂ = G x̄.
- **E2** the LOW-RANK version (Ŝ truncated to rank r ∈ {20, 100}, no floor) — the construction the review warned
  against, run as the comparison that shows whether the floor matters.
- **E3** the signal-covariance WHITENER: 062's operator recipe (`op272.fit_operator`) with the full-row covariance
  replaced by Ŝ + floor, used in the operator's place on the incumbent rows — the corrected-directions object.
- **E4** 271's diagonal per-gene Wiener (its (c), `wiener_own`) and **E5** 272's eigenbasis floor (c = 1) — the
  ALREADY-FAILED baselines every arm must beat; **E6** the properly cross-fitted operator of 272 (construction B) if its
  arrays are reusable.
- Every row arm is norm-matched to the mean's norm (263's `:nm`, shape-only, gate +0.0019) AND reported literally
  (amplitude-changing, gate +0.0089); price both.
- **Controls, frozen:** (Cs) SHUFFLED TARGET PAIRING — Ŝ built from view a of row t and view b of row π(t) for a fixed
  derangement π (expectation zero signal): an arm built on Ŝ_shuf must not gain; (Cw) WRONG-SOURCE Ŝ (another line's
  halves) — a gain there means the covariance is generic, not response-specific; (Cr) an energy-matched random PSD
  covariance; (Cm) a mirrored / reversed gain (G' = 2I − G) as 271/272 used.

## Stage A — held VIEWS (CD4, 280's bank; ≤ 1.5 h; the cleanest test on disk)
Views = the (donor, guide) pseudobulks of a target within one state (Stim48hr primary; Rest and Stim8hr secondary).
**Fit Ŝ, N̂ on the 500 NON-panel targets' views only** (`extra280_*`; independent views = the two guides of each
(target, donor), and separately the donor pairs) and report the in-sample fit on the panel as a labelled comparison.
Estimate each panel target's response from the views of TWO donors (E1 / E2 vs E0 = the mean of those two donors' rows)
and score against the held two donors' consensus in 280's residual cosine — 280's Q2 construction, all three disjoint
pairings, matched UMI budget and full — paired target bootstrap of Δcos vs E0; also the 273-style `pds` with all panel
targets as distractors and the held consensus as truth. Then **held STATE**: Ŝ fitted on Stim48hr views, the estimate
scored on the held donors' Rest views — does the gain survive a state change. KA: E0 reproduces 280's Q2 k = 2 row at
0.0e+00. Controls Cs / Cw (Ŝ from another state) / Cr on every cell.

## Stage B — held CONTEXTS (PLOLO, `k562_gwps` / `orion_hct116` / `orion_hek293t`; ≤ 4 h)
267's half-banks are the views (batch-parity halves of every bank row; `halfrows271.py` converts them to production row
units verbatim). Per donor line fit Ŝ_j, N̂_j on ALL its bank rows EXCEPT the panel targets (held targets never enter
the fit; state n). Ŝ is g × g on the official axis (≈ 17.8 k, 1.3 GB float32; eigendecomposition in minutes; 272 handled
this rank at 5–6 GiB) — one fold per process. E1 / E2 replace each donor's production un-whitened row and go through
271's plumbing unchanged; E3 replaces the operator for the incumbent rows. Read per-target `pds`, own cosine,
`cos_agg`, 262's groups (anti-aligned / marginal band / saturated), knee crossings up / down at 0.02, the `nmae` identity
leg (271's `identity271` construction), paired target bootstrap (`boot265.py`), the 3/3-fold rule, leave-top-3-out
retention; price with 261's recipe (pessimistic / point / optimistic, class factors per READ-FIRST §3). **KA rows:**
271's `ebrow` and `relw` through this experiment's plumbing reproduce 271's per-fold `pds` at 0.0e+00; the BR249
retrodiction +0.00472 at 0.0e+00 vs 257; E1 with N̂ = 0 and φ = 1 equals E0; the shuffled-pairing Ŝ_shuf's top
eigenvalue vs the real Ŝ's (report the ratio). **Gate (frozen):** pessimistic ≥ +0.0019 shape-only / +0.0089
amplitude-changing, `pds` positive on 3/3 folds with the pooled interval above 0, retention ≥ 50 %, every control below
half the arm's Δ`pds`. `kolf2` as a fourth fold if ≥ 1 hour remains (labelled).

## Stage C — the direction-resolving reference (273's instrument, union OFF; ≤ 1 h)
The same arms composed through `lib273` / `instrument273.py` (base + CD4 re-extracted, 85 targets, the H1 union never
loaded): per-target `pds` and own cosine vs BASE and vs the incumbent, with the bootstrap. If Stage B's `pds` reads as
a cosine gain here but not as `pds`, or vice versa, say so (273 §5's BLURS reading) — this stage says whether a PLOLO
gain is the truth's noise or the prediction's transfer.

**Production rendering (the review's third stage) is NOT this mandate**: no package, no build. A Stage B arm that
clears the gate goes back to the coordinator for a preregistered build mandate.

## Mechanism reading (measured, labelled)
Principal angles between Ŝ's leading directions and the fitted operator's; the share of the incumbent's whitened
prediction energy on directions where Ŝ's eigenvalue is below N̂'s; how much of each estimator's change lies in 262's
marginal band vs the anti-aligned group.

## Deliverable
`PREREGISTRATION.md` (hashed) and `est286.py` (sha256 recorded); `RESULTS.md`: the KA table; Stage A (Δcos vs E0 per
pairing, held state, controls, intervals); Stage B per arm × fold (`pds`, cosine, groups, knee, identity leg, price
pessimistic / point / optimistic, retention, controls) and the verdict per the gate; Stage C; the mechanism reading;
one line: **HELD VIEWS (CD4): Δcos … [lo, hi] vs the two-donor mean (shuffled …, wrong-source …); HELD STATE: …; HELD
CONTEXTS (PLOLO): best arm … Δ`pds` … [lo, hi], 3/3: YES / NO, priced … pessimistic vs gate … → CLEARS / DOES NOT
CLEAR; FLEX REFERENCE: …**. `work/` arrays per fold. No build, no submission. Then stop.
