# 361 — The CPM block is a coin flip (0.508); all-down buys +0.0081 fid and nothing else

Status: **COMPLETE** (kill criterion passes on the letter for one design-time counterfactual;
read the fragility notes before running an arm — margin over the bar is 1.4e-4 raw).

## Known-answer reproduction (own code, `fid361.py`)

| quantity | reproduced | recorded | match |
|---|---|---:|---|
| fid mean (85 targets) | 0.49755452 | 0.49755452 (`per-target-cap15.csv`) | max per-target diff **0.0** |
| nmae mean (58 eligible) | 0.95083646 | 0.95083646 (`metrics-cap15.json`) | max per-target diff **4.4e-16** |
| reach mean (82 scored, 3 NaN) | 0.25434726 | 0.25434726 | max per-target diff **0.0**, NaN set identical (PCCB/TRAFD1/TSKU) |

Conventions (from `cell_eval2/metrics/direction.py` `_components` + BRIEF §4): scored set =
{q_pred<0.05}∩{ref lfc finite & nonzero}, own target excluded from n_pred/n_real/k;
F=k/max(n_pred,n_real); n_real=0 & n_pred>0 → k/n_pred. On this panel every ref lfc is
finite-nonzero (adjudicable frac 1.0000) and every delivered lfc commits (frac 1.0000), so the
edge conventions never fire — except the 3 empty-pool reach NaNs. Reach ordering
(sig-first, q, p, |lfc| desc, gene) re-implemented independently; matches to 0.0.

## What was measured

Exact arithmetic on saved arrays only (`stats-cap15.npz`, `reference.npz`,
`effects.npz::cap15` mapped by gene name; 2 universe genes absent from the effects axis
zeroed). No emission, no 323 run. Head = top-400 |design| per target; CPM block = top-2200
basal-CPM (`ref_lin`) among non-head; reconstruction covers our calls at median 367 / 2200 /
32 per target (block calls/target exactly 2200 every target — the whole block is called, so
"calls-only" and "all-block" scopes are identical; verified bit-identical nmae/reach).
Counterfactuals flip ONLY block-gene lfc signs (|lfc|, q, p, ordering unchanged):
a. all-down; b. LOO reference-common-vector sign (truth ceiling); c. LOO design-common-vector
sign (design-time); d. down on top-half block by CPM, design sign elsewhere; e. current.
Paired 10,000-target bootstrap (seed 20260909) of Δ vs current; H1-local overall from
`work/h1-anchors.json` (scaled=(raw−b)/(r−b), overall Δ = Σ/6 with mse/pds/jac at 0).

## Results — decomposition (per-target means unless noted)

| quantity | value |
|---|---:|
| n_pred median / n_real median / mean / k median | 2599 / 89 / 765.2 / 1324 |
| share of our calls that are reference-significant (pooled) | 0.1134 |
| k-accuracy pooled, all calls / on refsig calls / on non-sig calls | 0.5131 / 0.5675 / 0.5061 |
| k-accuracy per-target mean, on refsig / on non-sig | 0.6451 / 0.5076 (refsig mean is small-n inflated) |
| head400: pooled calls 31,197, acc **0.5418**, contribution to mean F 0.0739 | |
| **cpmblock2200: pooled calls 187,000 (84.6%), acc 0.5083, contrib 0.4175** | |
| other-called: 2,763 calls, acc 0.5092, contrib 0.0061 | |
| reference down-tilt: block calls pooled 0.5178 (per-target range 0.444–0.640, median 0.511) / all adjudicable 0.5034 / ref-significant 0.5251 | |

The hypothesis is confirmed descriptively — the block is a coin flip — but the tilt it would
exploit is only 51.78% (nothing like 155's 59.7% on its pool), and the head itself is only 0.542.

## Results — counterfactuals (raw means; Δ with paired 95% CI)

| arm | fid (Δ) | nmae (Δ, lower better) | reach (Δ) | H1-local overall Δ |
|---|---:|---:|---:|---:|
| e current | 0.497555 | 0.950836 | 0.254347 | — |
| a all-down | 0.505696 (**+0.00814** [0.00179, 0.01481]) | 0.955668 (+0.00483 [−0.00208, +0.01289]) | 0.270146 (+0.01580 [−0.02314, +0.05662]) | **+0.00335** |
| b LOO-ref (ceiling) | 0.534794 (+0.03724 [0.02896, 0.04559]) | 0.941704 (−0.00913 [−0.01443, −0.00358], resolved beneficial) | 0.264264 (+0.00992 [−0.03874, +0.05579]) | +0.01215 |
| c LOO-design | 0.492437 (−0.00512, resolved harm) | 0.957126 (+0.00629, resolved harm) | 0.207671 (−0.04668, resolved harm) | −0.01194 |
| d top-half-down | 0.503692 (+0.00614 [0.00235, 0.00991]) | 0.954889 (+0.00405 [−0.00083, +0.00998]) | 0.257193 (+0.00285 [−0.02890, +0.03394]) | +0.00064 |

Scaled points (a): fid +0.00997, nmae −0.00832, reach +0.01847. Flip rates on block calls:
a 0.495, b 0.503, c 0.495, d 0.248. LOO ties on block genes: ref 0, design 256 (kept current).

## Verdict against the PLAN.md kill criterion

**Not dead — (a) passes, by 1.4e-4 raw.** Design-time arms: (a) +0.008142 ≥ +0.008 with nmae
interval [−0.00208, +0.01289] not entirely harmful and reach [−0.02314, +0.05662] not entirely
harmful → passes as written; (c) resolved-harmful on all three; (d) +0.00614, below the bar.
Three facts the PI must weigh: (1) the gain IS the tilt, exactly — tilt-arithmetic predicts
(a)'s Δfid at +0.008142 vs measured +0.008142 (0.0e+00 error), so on A/B/C the gain is whatever
the block-call down-tilt there is, minus 0.5083; at tilt 0.50 it is zero; (2) nmae points the
wrong way (−0.0083 scaled) — all-down sacrifices the 48% of block calls currently right on
up-genes, and the nmae CI covers +0.0129 of harm; (3) the fid CI [0.0018, 0.0148] sits almost
entirely below the bar — only the point clears it. The ceiling (b) proves real sign information
exists in the shared cross-target vector (+0.037 fid with nmae *improved*), but it is
truth-derived; the design-time version (c) is harmful everywhere — the design's common vector
points the wrong way (consistent with 347: design/common agreement is coin-flip).

If the PI runs the 323 arm, the emitter change is: after the design row is fixed, set
`design[block] = −|design[block]|` on the 2,200-gene CPM block (non-head by |design| rank),
before profile construction. The partition is |·|-invariant so head/block membership cannot
shift under the intervention; call set, magnitudes, q/p all downstream-identical except signs.
**pds is unpriced by this experiment** (sign flips change the profile; needs emission + scorer).

## Limits

H1 development panel (85 targets, one held-out line); not an A/B/C forecast — the tilt is the
entire mechanism and H1's 51.78% need not hold on Jurkat/HeLa/CAL33. fid/nmae/reach recomputed
exactly but mse/pds/jac carried at 0 for the H1-local overall (sign flips move the profile, so
a real arm's overall differs). Bootstrap CIs are target-resampling intervals on H1, not
transfer intervals. nmae cohort 58, reach cohort 82 — per-member valid-cohort means (report342
convention).

## Exact commands to reproduce

```
cd ~/Code/virtual-cell-challenge-2026
memory_pressure | tail -1   # proceed only if free >= 50%
uv run python experiments/361-fidelity-anatomy/fid361.py  # ~2 s, writes work/
```

Outputs: `work/decomposition.json`, `work/counterfactuals.json`, `work/per-target-361.csv`
(85 rows: per-arm fid/nmae/reach + n_pred/n_real/k + per-target block tilt).

## Peak memory and wall time

Peak RSS ≈ 263 MB (all arrays ≤ 85×10,779 float64); wall ≈ 2 s for the full script
(known-answer + decomposition + 5 scenarios + bootstrap). Sensitivity re-check
(scope-equivalence, tilt arithmetic) similar. No process exceeded 4 GiB; no scratch files kept.
