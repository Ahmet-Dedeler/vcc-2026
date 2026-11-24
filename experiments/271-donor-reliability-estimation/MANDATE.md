# 271 — Response ESTIMATION before selection: the reliability of each DONOR's rows at its own depth, and a noise-aware donor estimator (reliability-weighted, empirical-Bayes shrinkage) priced on PLOLO. Stage 2 of workstream 267. ≤ 5 hours. Measurement only. No build, no submission.

You are a research subagent. Work in `experiments/271-donor-reliability-estimation/`. **Edit nothing outside it.** No
build, no `vcc`, no `cell-eval2`, nothing under `submissions/`, `data/harness/` or `data/` (half-banks go under `work/`).
`uv run python`. One raw-cell extraction or PLOLO fold per process, peak RSS under ~10 GB, nothing else heavy is
running. **Time box: write up after ~5 hours.** `docs/RESEARCH-PROGRAMME-2026-09.md` binds you (name the unit, refit
inside the split, strong baselines, paired bootstraps, coherent pessimistic pricing from READ-FIRST §3).

## Why (R267 §3, §6, §"What stage 2 would be"; R159 §1.3)
267 showed the donor-selection oracle's +0.046 `pds` was the held REFERENCE's sampling noise selected twice: at 70–95
cells per half the reference's split-half residual cosine has a median of only +0.03 … +0.05 (30 % of Orion targets
≤ 0), and 159 puts the resid-space reliability at 150 cells at 0.05–0.12 for the Orion lines. **Our DONORS are those
same lines at those same depths.** So the un-whitened rows we average (K562-GWPS ~180 cells/target, Orion HCT116 and
HEK293T ~150) are themselves low-reliability estimates of each line's response, and an equal-weight mean of noisy rows
is not the best estimate of the shared response. The question is estimation, not selection: **per target, shrink each
donor's row toward the cross-donor consensus by the row's OWN measured noise (empirical Bayes / Wiener), so the
prediction leans on the reliable rows where they are reliable and on the consensus where they are not — does that
raise the transfer cosine on the marginal band (where +0.01 cosine buys +0.095 `pds`, R262 §4.2) and price positive on
the coherent-pessimistic rule?** The cautionary board point: a13 ("donor-agreement shrinkage", 2026-08-29) lost
−0.0203 vs a12 (`docs/BOARD-LEDGER.md` row a13). That arm shrank by BETWEEN-donor agreement (a signal 262/263 showed is
weak); this one uses WITHIN-donor sampling reliability (an independent, measurable quantity). The mandate requires
a13's construction as a baseline that the new estimator must beat on PLOLO before anything is said.

## READ FIRST
`docs/RESEARCH-PROGRAMME-2026-09.md`; 267's `RESULTS.md` (§1 inventory, §3 reliability, §6) and its `truth267.py` /
`work/half/<line>/bank-*` (the split-half banks per line — REUSE them as the donor halves: a fold's donors are the
other two lines, whose half-banks 267 already built), `sweep267.py` (the scoring path on 263's construction);
`experiments/159-transfer-ceiling/RESULTS.md` §1 (reliability laws, the 150 → 400 extrapolation, the h1-test 0.48);
262 §4.2 (the conversion curve) and its per-target arrays; 257's `price257.py` / 261's `price261.py` (the price recipe);
the a13 experiment (find it by `docs/BOARD-LEDGER.md` and `docs/LOG.md` 2026-08-29; what exactly it shrank and how);
`experiments/162-retrodicting-validator/plolo.py` (`build_fold`, how donor rows enter: `E062.prepare_main_source`,
the operator, `build_baseline`); `scripts/make_submission_pm_cd4.py` around `donor_weights` / `build_v9_pre_whiten_effects`
(where a per-target, per-donor weight would enter production — name the one place; do not edit).

## Stage 2 (this mandate)
1. **Donor-row reliability, per line and per target** (the scientific table): from 267's half-banks, cos(row_A, row_B)
   per target for K562-GWPS, HCT116, HEK293T at their own depth, in the un-whitened residual space AND after the
   fold's operator (whitened); distribution, share ≤ 0, relation to cell count and to 262's groups; Spearman–Brown to
   the full donor depth. KA: the union of halves reproduces the fold cache's donor rows (`pre`) to the digit.
2. **Noise-aware estimators (frozen before measuring; write them in `est271.py` with their formulas):**
   (a) **reliability-weighted mean**: per target, weights ∝ r_j / (1 − r_j) with r_j the donor row's split-half
   reliability (Spearman–Brown-corrected to full depth), normalised over the donors carrying the target;
   (b) **empirical-Bayes row shrinkage**: each donor row shrunk toward the equal-weight consensus by its estimated noise
   share (1 − r_j), then averaged; (c) **per-gene Wiener shrinkage** of each donor row using the row's own between-half
   gene-wise variance (a diagonal filter); (d) a13's construction re-implemented as the BASELINE (between-donor
   agreement weighting); (e) equal weights = `cd4svd5` (KA at 0.0e+00). Everything downstream as `cd4svd5` (whitening,
   CD4 deflated k = 5 resolved, union at 0.28).
3. **Measurement, three cancer folds, per fold and pooled:** `cos_agg`, macro, `pds` (the scorer's per-target statistic),
   per-target Δ vs `cd4svd5`, paired target bootstrap (≥ 10,000), by 262 group and by the reference-reliability tercile
   (267's), the share of marginal-band targets crossing 0.02; delivered amplitude ratio; the `nmae` identity.
   **Controls:** shuffled reliabilities (rows permuted across targets within donor), reversed weights, and (d).
   **The reserved `orion_hek293t` half-B truths are NOT to be used** (267's reserve; the coordinator scores reserves).
4. **Price** each estimator on the coherent-pessimistic / point / optimistic scenarios (re-dose band k_pds 0.85 / 1.22,
   k_cos 1.36 / 3.76 as pessimistic for the sign of the leg; the class for a base re-estimation is unmeasured — say so).
   CLEARED = pessimistic ≥ +0.0019 over `cd4svd5` AND `pds` > 0 on 3/3 folds with a bootstrap interval excluding 0 AND
   leave-top-3-out retention ≥ 50 % AND beats (d) and the shuffled control by more than their bootstrap width.
5. **The 400-cell question (scientific, cheap):** with 159's laws and 267's numbers, state the expected reliability of the
   CHALLENGE references (400 cells, Flex chemistry: h1-test measured 0.48 at 400) versus the Orion-like 0.10–0.27, and
   what that implies for the ceiling of per-target cosine on A/B/C. Numbers, not prose.

## Deliverable
`RESULTS.md`: the donor-reliability table; the estimators with formulas; the measurement and price tables with the
controls and (d); the 400-cell statement; one line **CLEARS / DOES NOT CLEAR** with the retained numbers; if it clears,
`PREREGISTRATION-DRAFT.md` in 249's form (five falsifiers, three scenarios, the production variable it would need — a
per-target, per-donor weight table consumed at `build_v9_pre_whiten_effects`, named, not made). A clean negative —
"donor noise is real but its correction does not move the transfer cosine" — is a valid result that closes this
construction. No model name, no submission. Then stop.
