# 359 — Expression-gated propagation: gate adds nothing; smoothing keeps 340's tradeoff

**Status: COMPLETE — NEGATIVE (all three kill clauses fire).**

## Known-answer reproduction (own code, `ka359.py`)

| quantity | reproduced | recorded (341) | match |
|---|---|---|---|
| pooled \|response\|-weighted sign, cross-line diagonal | 0.58939776 | .589398 | \|dw\| = 2.4e-7 (rounding) |
| pooled retrieval, cross-line diagonal | 0.79461386 | .794614 | \|dr\| = 1.4e-7 (rounding) |
| max \|ours − 341 saved `crossline_diag`\|, 6 lines | 0.0 | — | **bit-exact** |
| `B.known()` historical checks (5 banks) | 0.0 ×5 | 0.0 | exact |

Gate: **PASS**. End-to-end replay in 341's split regime
(`P.order(targets,'334-split|')` halves, gain on train, predict test, stitch).
Basal medians over the 512-gene frame are comparable across lines
(1.04–1.31, IQR 0.75–0.87), so per-line-median thresholds are density-matched
by construction (all hard-gate node densities exactly 0.5).

## What was measured (inputs, estimator, controls)

Six-line source frame (1,024 targets × 512 coords, `selection.json` +
`B.load`; truth = unit-RMS `direction`). Base `s0` = donor-mean ×
donor-fitted per-gene gain, identical for all arms. Propagation
`P_L = D_L^{-1}(A ∘ G_L)`, `s(λ) = (1−λ)s0 + λ·s0·P_Lᵀ` (1–2 steps),
λ ∈ {0.25, 0.5, 1}, λ chosen by inner-donor LOLO on train targets only
(held line never touched) → **λ = 0.25 for all four configured selections**
(hard.s1 / soft.s1 / hard.s2 / ungated.s1).
Priors: STRING KNN (directed 4,080 edges / symmetrized 5,842; reused
340 topology bit-for-bit) and CollecTRI (42 G×G edges over 29/512 genes,
56 T×G edges over 19/1024 targets; smoothing sign +mor, seed sign −mor —
see PLAN.md for the convention derivation). Gates: own-hard (median),
own-soft (sigmoid, s = 0.5), ungated, shuffled-node random ×3 (exact same
density), wrong-recipient hard ×5 (best reported). Endpoints:
`C.sign_endpoints` weighted sign (primary) + top-64, `P.endpoints`
retrieval (guard), `C.interval` paired + crossed CIs. 72 arms + diagonal.

## Results (pooled over 6 lines × 1,024 targets; deltas vs diagonal)

| arm @ donor-λ=0.25 | pooled wsign | Δ vs diag | target 95% | crossed 95% | retrieval | top-64 |
|---|---|---|---|---|---|---|
| diagonal | 0.589398 | 0 | — | — | 0.79461 | 0.64252 |
| **str.dir.hard.s1 (PRIMARY)** | 0.590106 | **+0.000708** | [+0.00044,+0.00098] | [−0.00025,+0.00194] | 0.78814 | 0.64415 |
| str.dir.soft.s1 | 0.590975 | +0.001577 | [+0.00118,+0.00198] | [+0.00030,+0.00306] | 0.78815 | 0.64573 |
| str.dir.ungated.s1 | 0.590954 | +0.001557 | [+0.00118,+0.00194] | [+0.00051,+0.00265] | 0.78852 | 0.64582 |
| str.dir.rand{0,1,2}.s1 | 0.58960/0.58973/0.59002 | +0.00020/+0.00033/+0.00063 | all ≈[0,+0.0009] | — | 0.7914–0.7920 | — |
| str.dir individual wrong ×6 | — | **+0.00024…+0.00088** | each > 0, tiny | — | 0.782–0.796 | — |
| str.dir.WRONGBEST.s1 | 0.596954 | +0.007557 | [+0.00730,+0.00781] | stable | — | — |
| **RANDMAX(3) (max of 3 random gates)** | 0.596503 | **+0.007105** | [+0.00687,+0.00734] | [+0.00663,+0.00760] | — | — |
| ct.hard.s1 / ct.ungated.s1 | 0.589150/0.589285 | −0.000248/−0.000113 | — | — | 0.79437/0.79439 | — |
| seeded-ct (19 covered targets) | — | −0.000010 | — | — | 0.79459 | — |

Further rows: λ = 1 collapses everywhere (sign −0.019…−0.046,
retrieval −0.05…−0.16 — pure neighbor-average destroys target
specificity); λ = 0.5 ≈ 0 or negative; 2-step ≈ 1-step; symmetrized ≈
directed (hard +0.00078, ungated +0.00173). Full sweep + per-line vectors:
`work/endpoints359.npz`, `work/results359.json`.

Per-line, own-hard vs ungated Δ vs diag (wsign):

| line | diag | hard − diag | ungated − diag |
|---|---|---|---|
| k562_gwps | 0.588676 | +0.000021 | +0.001295 |
| orion_hct116 | 0.579402 | +0.000386 | +0.000376 |
| orion_hek293t | 0.563000 | −0.000764 | +0.000046 |
| rpe1 | 0.609090 | +0.003433 | +0.003480 |
| nadig_jurkat | 0.591160 | +0.000173 | +0.001658 |
| nadig_hepg2 | 0.605059 | +0.001000 | +0.002485 |

Ungated ≥ hard on 5/6 lines (hct116 ties at +0.00001).

## Verdict against the PLAN.md kill criterion — DEAD, 3/3 clauses

1. Own-gated − best-wrong-gated = 0.000708 − 0.007557 = **−0.00685 ≤ +0.005** → dead.
2. Gated (+0.000708) ≤ ungated (+0.001557) → dead (also per-line 5/6, and soft ≈ ungated).
3. Gated ≤ diagonal is false (+0.0007) but the alive bar (≥ +0.015 over diagonal) is missed **21×**.

The +0.0076 WRONGBEST is a **max-of-5 selection artifact, not biology**:
max-over-3-RANDOM-gates gives +0.0071 (same size, zero biology), and the
per-(line,target) winner is uniform over sources (994–1087/6144 each ≈ 1/6 —
no source is systematically best). Own-hard sits inside the wrong-gate
distribution (+0.0002…+0.0009). Recipient specificity = zero.

Two honest sub-findings, neither shippable: (a) ungated STRING smoothing at
λ = 0.25 carries a small real sign gain (+0.0016, both CIs exclude zero) —
1/10 of the alive bar; (b) it pays 340's price: retrieval −0.006 on every
smoothed arm (340's cosine-up/retrieval-down tradeoff, reproduced on the
sign endpoint). The CollecTRI mechanism check passed as pre-registered
(≈ diagonal: only 29/512 genes movable).

## Limits

Six-line source-frame unit-RMS directions; no count-stage/DE-member pricing,
no A/B/C forecast. Targets correlated across lines; study/line confounding
as in 341. τ = per-line median is arbitrary but density-matched, and the
permuted-node random control (exact same density) brackets it. Soft-gate
slope s = 0.5 fixed, not tuned. Target-basal missing for 109–138/1024
symbols on 4 banks (fallback gate = 1; affects the seeded exploratory only).
CollecTRI arms are sparsity-bounded near zero by construction (42 edges).
WRONGBEST/RANDMAX are max-statistics, unbuildable without labels — reported
only to defuse the +0.0076. No 323 run, no downloads, no reserve contact.

## Exact commands to reproduce

```bash
NUMBA_NUM_THREADS=4 OMP_NUM_THREADS=4 uv run python experiments/359-expression-gated-propagation/ka359.py   # gate: PASS
NUMBA_NUM_THREADS=4 OMP_NUM_THREADS=4 uv run python experiments/359-expression-gated-propagation/run359.py  # arms + work/results359.json
```

## Peak memory and wall time

`ka359.py`: ~40 s wall. `run359.py`: line loop 28→46 s cumulative
(λ-select + 6 lines), ≈2–3 min total incl. intervals/summary. Resident
budget ≈ 0.5 GB (direction 6×1024×512 f64 = 201 MB, P-cache ~80 × 512² ×
8 B ≈ 170 MB, endpoints small; Orion 2.7 GB delta.npy via mmap window
only) — far under the 12 GiB cap; `memory_pressure` free 73% before,
≥68% throughout; threads capped at 4.
