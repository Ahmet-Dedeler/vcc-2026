# 406 — The reach front breaks inside the head, but neighbourhood order keys do not move it

Status: **COMPLETE** (NEGATIVE on the order lever; descriptive table delivered). PI-run 2026-09-16,
`run406.py` → `work/results406.json` + `work/descriptive.pkl`, full log `work/run406.log`.

## Known-answer reproduction (own code, before any new number)
- Delivered reach (E0 signs, delivered order): **0.25434726**, 82/85 targets scored — matches the
  assignment's 0.25435 and 388's `walk_reach` value to <1e-8.
- E0 called-block precision: **0.5083** — matches 390 exactly.
- Reach walk copied VERBATIM from `knn397.py` (per-target curve check `assert |R−R399|<1e-12` passed
  on all 82 scored targets).
- Shipped-sign construction (six-bank K=25 shrink 400, `knn397d.py` logic) cross-checks the 399
  prereg addendum exactly: block precision **0.5338**, reach(SG399, delivered order) **0.25190**.
- Chunked Orion voter (mmap + ≤500-row chunks) reproduces verbatim `votes()` **bit-exactly**
  (int64 counts, NaN-identical P) on all small banks fully and on a 200-row Orion subset.

## What was measured (inputs, estimator, controls)
H1 analytic on 323's saved arrays, no emission. Signs fixed at SG399 (kNN on called block,
delivered elsewhere) for the descriptive and every variant, so fid/nmae/jac are unchanged by
construction. |z| := |E1| (cap15 design magnitude — the 323 frame carries no z-statistic; the
head-defining quantity). conf := |Ps−0.5| from the SHIPPED shrunk Ps (NaN→0); agree := kNN sign
== delivered sign. Variants re-permute ONLY head genes among their delivered-order slots:
(a) delivered; (b) |E1|·(1+β·agree); (c) |E1|·(1+β·conf·agree−β·conf·disagree), β∈{0.5,1,2};
(d) demote (disagree & conf≥0.2) head genes to end of head; (e) deranged control — per-target
permutation of conf within head, re-run (cβ)+(d), plus deranged agree for (b,β=1).
Endpoint: reach Δ vs (a), paired 10,000-target bootstrap CI (seed 406, 82 scored targets).
Break rank b := k*+1 (None if purity never dips below 0.9).

## Results 1 — descriptive: where the walk breaks (SG399 signs, delivered order)
Mean k* = 13.3 (388: ≈14 ✓). The breaking gene is ALWAYS a sign error under the walk's signs —
theorem, not data: if C/k ≥ 0.9 then (C+1)/(k+1) ≥ 0.9+0.1/(k+1), so a correct gene can never end
a qualifying prefix (hence `right = 0.0` everywhere by construction).

| scope (n_real tercile, edges 13/536) | n | mean n_real | mean k* | mean reach | never-break | break-rank mean/med | breaking gene head/block/uncalled | \|E1\| at break mean | conf at break | agree at break | P10 / P20 / P50 |
|---|---|---|---|---|---|---|---|---|---|---|---|
| all (82) | 82 | 793 | 13.3 | 0.252 | 0.16 | 16.1 / 4 | 0.62 / 0.29 / 0.09 | 0.104 | 0.055 | 0.41 | 0.757 / 0.700 / 0.671 |
| terc0 weak | 28 | 6 | 3.0 | 0.666 | 0.46 | 3.5 / 3 | 0.00 / 0.60 / 0.40 | 0.038 | 0.037 | 0.40 | 0.650 / — / — |
| terc1 mid | 27 | 187 | 5.2 | 0.059 | 0.00 | 6.2 / 4 | 0.70 / 0.30 / 0.00 | 0.097 | 0.054 | 0.44 | 0.715 / 0.628 / 0.611 |
| terc2 strong | 27 | 2216 | 32.0 | 0.016 | 0.00 | 33.0 / 8 | 0.89 / 0.11 / 0.00 | 0.148 | 0.066 | 0.37 | 0.815 / 0.761 / 0.716 |

Reading: the front lives in the head for mid/strong targets (70–89% of breaks) but NOT for weak
targets (0% head; the walk breaks in block/uncalled genes at rank ~3). 59% of breaking genes are
kNN-disagreements — yet at low confidence (0.055). Shrink-400 head conf is nearly flat (median
0.026, p90 0.083; only 1.6% ≥ 0.2), so any conf-graded key is a weak lever by construction; the
demote mask touches median 3 head genes/target (max 20).

## Results 2 — order variants (reach, Δ vs delivered order, 95% CI)
| variant | reach | Δ | break-rank mean/med | never |
|---|---|---|---|---|
| (a) delivered | 0.25190 | — | 16.1 / 4 | 0.16 |
| (b) β=0.5 agree | 0.25325 | +0.00136 [−0.00174, +0.00464] | 17.8 / 5 | 0.16 |
| (b) β=1 | 0.25278 | +0.00089 [−0.00238, +0.00426] | 16.8 / 5 | 0.16 |
| (b) β=2 | 0.25091 | −0.00099 [−0.00609, +0.00333] | 16.6 / 5 | 0.16 |
| (c) β=0.5 conf | 0.25085 | −0.00104 [−0.00505, +0.00211] | 16.2 / 4 | 0.16 |
| (c) β=1 | 0.25266 | +0.00077 [−0.00109, +0.00287] | 16.7 / 5 | 0.16 |
| (c) β=2 | 0.25294 | +0.00105 [−0.00092, +0.00326] | 17.3 / 6 | 0.16 |
| (d) demote dis@0.2 | 0.25236 | +0.00046 [−0.00003, +0.00140] | 16.2 / 4 | 0.16 |
| (e) deranged (c) β=0.5/1/2 | — | −0.00070 / +0.00052 / +0.00122 (CIs span 0, widths 0.006–0.007) | — | — |
| (e) deranged (d) | — | +0.00034 [−0.00037, +0.00134] | 15.9 / 4 | 0.16 |
| (e) deranged agree (b,β=1) | — | −0.00000 [−0.00269, +0.00296] | 12.5 / 4 | 0.16 |

## Verdict against the PLAN.md kill criterion
No variant clears +0.005 (best +0.00136, CI crosses zero), and the deranged controls read the
same ±0.001 scale as the live variants — order-statistic noise, not signal. **The order stays**
(delivered). The descriptive table above is the deliverable: it prices the head at 62% of breaks
overall but 0% for weak targets, which bounds any future within-head order lever.

## Limits
- H1 development panel only (85 targets, 82 scored) — not an A/B/C forecast; H1 reach noise is
  ±0.03–0.04 on Δ (cf. 397/399 notes).
- Conf:scale coupling — shrink-400 flattens head conf, so (c)/(d) were dosed weak by construction.
  An unshrunk-P sensitivity was NOT run (397e already closed unshrunk head flips: 245 calls, noise).
- Slot-preserving within-head permutation only; cross-region promotion was closed by 397b.
- Sign errors outside the reference pool are invisible to the walk; `right=0.0` is a walk theorem.

## Exact commands to reproduce
```
memory_pressure | tail -1   # gate: free ≥ 50%
NUMBA_CACHE_DIR=experiments/406-reach-front/work/numba-cache NUMBA_NUM_THREADS=4 OMP_NUM_THREADS=4 \
  python3 experiments/406-reach-front/run406.py   # writes work/results406.json, work/descriptive.pkl
python3 -c "import pickle,numpy as np; ..."       # breaking-gene |E1| aggregates (see turn log)
```
Only reads outside `experiments/406-reach-front/` (knn397.py frame exec, 323 work arrays,
donor banks mmap'd). No writes elsewhere. Probe scratch lives in /tmp (probe406.py), not committed.

## Peak memory and wall time
`/usr/bin/time -l` re-run: **6.3 s wall, peak RSS 6.1 GB** — dominated by touched mmap pages while
streaming the two 2.6 GiB Orion banks for the family marginal; anonymous footprint ≈ 1–2 GiB
(GWPS 323 MB + centred/normalised copies). Under the 8 GiB cap; one job at a time observed.
