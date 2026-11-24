# 299 — Competition pilot, tonight: does fitting the whitening operator on RELIABILITY-SELECTED training rows (the non-panel bank perturbations whose two independent halves agree) change the held-context transfer, at fixed row count, donor proportions, gene axis and coverage? Preregistered, priced on the coherent-pessimistic rule, with matched-random and reversed-selection controls, on the three PLOLO cancer folds and the Flex reference. ≤ 3 hours, HARD STOP at 21:30Z. No build, no submission.

You are a research subagent. Work in `experiments/299-reliable-training-rows/`. **Edit nothing outside it.** No build,
no `vcc`, no `cell-eval2`, nothing under `submissions/`, `data/harness/` or `data/`. `uv run python`; one fold per
process, serial; peak RSS under ~10 GiB; `memory_pressure` ≥ 40 % free before each fold. **Time box: 3 hours, and a
HARD STOP at 2026-09-06 21:30Z** — at 21:15Z write RESULTS.md with whatever folds are complete and stop; the
coordinator needs the reading before the 00:00Z slot reset. `docs/READ-FIRST.md` §3 binds you (one variable; price all
six members on the coherent-pessimistic rule; known-answer row first; the gate). Both reserves and R3 are SPENT and are
not needed here; PLOLO's fold caches and 267's half-banks are the only inputs.

## Why (third external review, 2026-09-06 20:30 Ankara; F3, R272 §0.10, R286 §6)
The incumbent's whitener is fitted on the donors' full rows at n/p ≈ 1.6; 272 and 286 showed its eigenbasis is ~78 %
sampling-noise directions that carry ~80–90 % of the whitened prediction energy, and that re-weighting that basis (272),
replacing its covariance by an unbiased signal covariance (286 E3) or shrinking the OUTPUT rows (271) all lose `pds`.
**What has not been tested is changing WHICH rows teach the operator**: fit it on the training perturbations whose two
independent halves agree (reliable rows carry a higher signal share), at the SAME row count, the same donor
proportions, the same axis, the same γ / λ, with every panel target still covered by the donor rows. 091 found that
training-set composition matters at equal row count (`experiments/091-basis-prune/RESULTS.md`); 193 closed γ / λ /
exponent / row COUNT (`experiments/193-whitener-shrinkage/RESULTS.md:265`), not row selection by reliability. The
mechanism evidence cuts both ways — a cleaner covariance moves toward 286's E3, which lost — so this is a hypothesis
with a measured reason on each side, not a predicted win. It is cheap: the reliabilities and the fold machinery exist.

## READ FIRST
`experiments/162-retrodicting-validator/plolo.py` (`build_fold`, `Operator`, `fold_inputs`), `experiments/272-noise-
corrected-operator/op272.py` (`fit_operator` = 062's recipe on a chosen row set; the noise spectrum), `experiments/286-
crossproduct-estimator/{stage286.py, est286.py}` (267's half-banks → production row units for ALL bank rows; the
per-row halves A / B), `experiments/271-donor-reliability-estimation/{halfrows271.py, est271.py}` (row_cos, Spearman–
Brown), `experiments/263-donor-selection-arbiters/lib263.py` (`stage`, `combine`, `evaluate`), `experiments/271-…/
price271.py` + `experiments/261-union-deflation/price261.py` (the six-member price; the BR249 retrodiction row),
`experiments/265-majority-permutation-test/boot265.py`, `experiments/273-flex-reference-instrument/{lib273.py,
instrument273.py}` (the union-OFF Flex-400 instrument; how an operator variant is read there), 226's and 193's results
(row-count and shrinkage closures: the arms to NOT repeat).

## Protocol (write `PREREGISTRATION.md`, hash it into `work/PREREGISTRATION.sha256`, before any fold number)
1. **Known answers:** the fold operator refit with `op272.fit_operator` on the incumbent's row set reproduces the cached
   operator's action on the donor rows (state the tolerance; 272's KA); `cd4svd5` through `lib263` reproduces 263's
   `pds` per fold at ≤ 1e-9; the BR249 retrodiction +0.00472 at 0.0e+00; per-row split-half reliability r_j(t) for
   every bank row from 267's halves reproduces 271's panel-row reliabilities on the panel (0.0e+00 on the shared rows).
2. **The selection (frozen):** for each donor line j in a fold, rank its NON-panel bank rows by split-half reliability
   r_j (cos(A, B) of the row's two batch-parity halves, Spearman–Brown to full depth); the operator is refit on the
   TOP-q fraction of each donor's non-panel rows at q ∈ {0.5, 0.25}, plus ALL panel-target rows (coverage preserved),
   with the row count per donor otherwise held at the incumbent's proportions (state the counts). Everything else —
   γ 0.65, λ 0, the axis, the preserve, the CD4 deflation k = 5 at f = 0.34, the H1 pair union at g = 0.28 — verbatim
   from `cd4svd5`. Arms: `sel:q0.5`, `sel:q0.25`; **controls** at the same q: `rand:q` (a random subset of the same
   size, 3 seeds), `rev:q` (the LEAST reliable rows), `mag:q` (the highest-MAGNITUDE rows — the confound 285 / 288
   name: reliability tracks response size), and `full` = the incumbent (KA). Read every arm on the three cancer
   folds: per-target `pds`, own cosine, `cos_agg`, delivered amplitude, 262's groups, the 0.02-knee crossings, the
   `nmae` identity leg (271's construction), paired target bootstrap (`boot265.py`), 3/3-fold rule, leave-top-3-out
   retention; price with 261's recipe (pessimistic / point / optimistic; operator changes are a RE-SHAPE of occupied
   rows — k_cos 3.76, k_pds 1.22 on the pessimistic side).
3. **Gate (frozen, READ-FIRST §3):** pessimistic ≥ +0.0019 if delivered amplitude is within 0.5 % (shape-only) else
   ≥ +0.0089; `pds` positive on 3/3 folds with the pooled interval above 0; retention ≥ 50 %; `sel` above `rand` AND
   above `mag` AND above `rev` by more than half its own Δ`pds`. State each condition's outcome per arm.
4. **The Flex reference (if ≥ 40 min remain before the hard stop):** the passing arm, or the best arm, read on 273's
   union-OFF instrument (85 targets) — Δ`pds` and Δ own cosine vs the incumbent's matched rebuilt mean; same sign as
   PLOLO or not.
5. **Mechanism row (cheap):** for the `sel` operator, the share of directions where an independent noise estimate
   exceeds the eigenvalue (272's statistic) and the whitened prediction energy on them, next to the incumbent's 77–79 %
   / 87–92 % — does selecting reliable rows actually move the spectrum.

## Deliverable
`PREREGISTRATION.md` (hashed), `RESULTS.md`: the KA table; the row-selection counts per donor; the arm × fold table
(`pds`, cosine, amplitude, groups, knee, identity leg, price pessimistic / point / optimistic, retention, the three
controls) and the gate verdict per arm; the Flex re-read; the mechanism row; `work/pertarget299-<fold>.npz`; one line:
**RELIABLE-ROW WHITENER: best arm … Δ`pds` … [lo, hi], 3/3: YES / NO, controls (rand / mag / rev) …, priced …
pessimistic vs gate … → CLEARS / DOES NOT CLEAR; FLEX: …**. No build, no submission; the coordinator decides. Then stop.
