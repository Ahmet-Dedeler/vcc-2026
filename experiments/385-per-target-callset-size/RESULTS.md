# 385 — Per-target call-set size: predictable, but the rule loses to the constant. NEGATIVE

Status: **NEGATIVE** (rule priced end-to-end; kill criterion fails; constant K = 2,600 stands).

## Known-answer reproduction

Own cap15 emit+score through `run385.py` (380's call sequence: profile /
spread / rank-axis / `rc_hybrid_block` + `CallSet(400, 2200, 1.9)` / block_stats /
score), 138 s. Against `323/work/metrics-cap15.json`: **both digests exact,
mse_legs exact, all six raw members bit-identical (abs diff 0.0)** —
fid 0.49755452, jac 0.07475983, mse 0.97494624, nmae 0.95083646, pds 0.88963585,
reach 0.25434726 (`work/known-answer.json`, `metrics-cap15.json`).

## What was measured

1. **Predictability** (`predict385.py`, no emission): per-target n_real from
   `reference.npz` (q < 0.05, own target excluded; cross-checked exactly against
   `reference.json`: median 89, mean 765, 3 zeros) vs build-time observables:
   cap15 row L2 (full + universe), #{universe genes with |effect| > τ} for
   τ ∈ {0.05, 0.1, 0.2, 0.5}, pre-whitening row L2 (`stage.npz` uniform arm),
   DepMap |GeneEffect| of the target (Jurkat row + 17-model mean; HeLa absent
   from the bank, H1 verified absent from all 17 models). Stats: Spearman +
   LOO OLS of log(n_real+1).
2. **Pricing** (`k385.py` → `price385.py` × 3 → `report385.py`): LOO DepMap-mean
   predictor → n̂ → K_p = clip(round(c·n̂), 400, 4474), c ∈ {1, 2}, realised as
   `CallSet(min(400,K_p), K_p−min(400,K_p), 1.9)` per target (protected 400-|H|
   head inside; count set through the z_cap mechanism per the rc_emit
   docstring), vs constant K = 2,600 and a shuffled-K control (same multiset,
   deranged, seed 385). Paired 10,000-target bootstrap, seed 20260909, 323
   H1-local anchors for scaling. CSNK1E/PGAM1 absent from DepMap columns →
   median-imputed (disclosed). K_c1: 75/85 at floor 400; K_c2: 61/85 at floor.

## Results

Predictability (n = 85; LOO = leave-one-out OLS in log space):

| observable | Spearman (p) | LOO Pearson / Spearman | resid SD (const 2.63) |
|---|---|---|---|
| depmap_absMean17 | **+0.51 (8e-7)** | +0.45 / +0.49 | 2.31 |
| cap15_n_gt_0.2 | +0.49 (2e-6) | +0.00 / +0.26 | 2.77 |
| cap15_n_gt_0.5 | +0.41 (1e-4) | +0.33 / −0.11 | 2.48 |
| cap15_L2_universe | +0.19 (0.08) | +0.08 / +0.10 | 2.63 |
| prewhite_L2_full | +0.09 (0.41) | −0.06 / −0.03 | 2.67 |
| depmap_absJurkat | +0.19 (0.09) | +0.02 / +0.07 | 2.62 |

Gate cleared (0.51 ≥ 0.4) — but note the scale: resid SD 2.31 in log space is
~10× multiplicative error, i.e. ~20% of variance. 192's verdict (6% variance,
unbuildable) is softened, not overturned.

Pricing vs constant cap15 (raw Δ, 95% CI; H1-scaled in text):

| arm | fid | jac | reach | nmae/pds/mse | overall |
|---|---|---|---|---|---|
| Kc1 | −0.0959 [−0.131,−0.062] | −0.0482 [−0.064,−0.033] | +0.0278 [+0.009,+0.051] | 0 / 0 / 0 | **−0.0331** [−0.047,−0.019] |
| Kc2 | −0.0852 [−0.118,−0.053] | −0.0411 [−0.056,−0.027] | +0.0255 [+0.007,+0.049] | 0 | −0.0285 [−0.042,−0.015] |
| Kshuf | −0.0972 [−0.131,−0.063] | −0.0505 [−0.066,−0.035] | +0.0280 [+0.010,+0.051] | 0 | −0.0342 [−0.048,−0.020] |

Scaled: Kc1 fid −0.118, jac −0.113, reach +0.032. Direct paired Kc1−Kshuf:
fid +0.0015 [−0.0157,+0.0176], jac +0.0054 [−0.0060,+0.0176] scaled —
combined **+0.0070, 95% CI ≈ [−0.022,+0.035]**, i.e. indistinguishable from
zero and not ≥ +0.005 by any evidential standard.

Mechanism (means over 85): cap15 precision 0.513 / coverage 0.972 at n_pred
2600 → Kc1 precision 0.539 / coverage 0.773 at n_pred 503. The smaller head is
marginally purer but coverage collapses, and the heavy tail (mean n_real 765)
makes under-calling the tail cost more than over-calling the majority saves —
192 §3.2's coverage model, confirmed with a real predictor. reach gains
(+0.028) because the small head is a purer prefix; net is −0.033 overall.

## Verdict against the kill criterion

**KILL: the rule does not beat the shuffled-K control** (combined +0.007 ±
~0.03 scaled — the targeting adds nothing; the whole action is the size cut,
and the size cut itself loses −0.033 overall to the constant). Constant
K = 2,600 stands; per-target sizing is closed on H1 with a real (DepMap)
predictor, strengthening 192 §3.3. No production change.

## Limits

H1-panel (85 targets, one line) development evidence, not A/B/C forecasts; no
board conversion applied. DepMap predictor uses other lines' KO fitness (no H1
outcomes; H1 absent from the bank) — legitimate build-time prior, ~10× error.
Per-target counts realised exactly as specified; delivered n_pred ≈ K_p (floor
503–602 mean). pds/mse/nmae pinned by construction (verified 0.000000).

## Environment incident (disclosed, gated, no src/ edits)

Mid-series, a concurrent worker's uncommitted 319o edit added optional
`order_key` to `rc_hybrid_block` and `VCC_RC_ORDER` to the builder; my frozen
gate tripped before Kc2. Handling, all read-only wrt `src/`/`scripts/`:
`equiv385.py` proves inertness — (A) the five builder functions I import are
AST-identical HEAD vs worktree; (B) HEAD vs worktree `rc_hybrid_block`
(default path) emit **bit-identical blocks** on real target 0 with
CallSet(400,2200,1.9) (block sum 8,245,700); cap15/Kc1 ran pre-drift under the
frozen sha. `verify385()` pins both files' shas — any further concurrent edit
is a hard stop, and pins held for Kc2/Kshuf/report. One self-inflicted issue:
the gate's aliased module load wrote numba cache entries that break the shared
`src/` cache lookup (`ModuleNotFoundError` at first `wilcoxon_z` compile);
**PI: `src/vcc2026/__pycache__/rc_emit.wilcoxon_z-89.*` (Sep 11 ~20:00 UTC)
needs deleting** — I did not touch `src/` per the rules. My runs after the
gate used `work/numba-cache/` (recompile of unchanged jitted sources,
bit-identical). No `git`, no builder, no submission commands run.

## Exact commands to reproduce

```
uv run python run385.py cap15        # known answer, ~140 s
uv run python predict385.py          # predictability table, seconds
uv run python k385.py                # K rule tables
uv run python equiv385.py            # drift gate (needs the pins for below)
uv run python price385.py Kc1        # ~140 s each, ONE AT A TIME
uv run python price385.py Kc2
uv run python price385.py Kshuf
uv run python report385.py           # paired bootstrap vs cap15
uv run python contrast385.py         # direct Kc1-Kshuf + fid mechanism
```

Scripts: `run385.py` (harness), `predict385.py`, `k385.py`, `price385.py`,
`report385.py`, `contrast385.py`, `equiv385.py`; outputs in `work/`
(`metrics-{cap15,Kc1,Kc2,Kshuf}.json`, `per-target-*.csv`, `stats-*.npz`,
`results385.json`, `drift385.json`, `predict385.json`, `k385.json`).
4 harness runs used (cap limit 4). Peak RSS not instrumented; consistent with
323's ~5 GiB profile, threads 4; system free 86% at start, 36% at end (another
job started elsewhere — no further jobs launched).
