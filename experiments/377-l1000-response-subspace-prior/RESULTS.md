# 377 — Same-line L1000 response-subspace prior: wrong-line matches it, retrieval pays 2.4× — NEGATIVE

Status: **NEGATIVE** (all three PLAN kill criteria fire; no HeLa artifact produced).

## Known-answer reproduction (own code, before any new measurement)
Mean of `per-target-cap15.csv` grouped by metric vs `work/metrics-cap15.json` raw members:

| member | per-target mean | recorded | match |
|---|---|---:|---|
| fid | 0.4975545182292536 | 0.4975545182292535 | 1.1e-16 |
| jac | 0.0747598329813194 | 0.0747598329813194 | 1.4e-17 |
| pds | 0.8896358543417368 | 0.8896358543417369 | 1.1e-16 |
| nmae (58 eligible targets) | 0.9508364630325663 | 0.9508364630325662 | 1.1e-16 |
| reach (nanmean; 3/85 targets NaN) | 0.25434726355678233 | 0.25434726355678233 | exact |
| mse | — (absent from the csv by construction) | 0.9749462418923823 | n/a |

## What was measured
**Bank gap (deviation recorded in PLAN.md before running):** the mandate assumed held
level-5 *compound* values. On disk only `level5_beta_trt_sh.gctx` (238,351 shRNA
consensus rows × 12,328 genes) and `level5_beta_trt_xpr.gctx` (142,901 OE rows) exist;
no `trt_cp` values were ever downloaded (a cp beta file would be ~35 GB). Same-line
coverage in the shRNA bank: **HEPG2 18,727 rows, JURKAT 256, HEK293T 97, K562 0,
HELA 0, HCT116 0**. So the executable analogue is the **same-line shRNA-response
subspace** (closer to CRISPRi than compounds); K562/HELA legs are untestable and the
mandate's HeLa-artifact gate cannot fire.

**Design (units: per-target means over held CRISPRi targets, 623-gene axis):**
- Gene axis = 978 L1000 landmarks ∩ test-line axis ∩ all five donor axes → **623 genes**.
- Targets = test-line bank targets carried by ≥2 of the other five CRISPRi donors:
  **HepG2 2,156; Jurkat 2,262** (donor-count histogram ≈ all 4–5× covered).
- Base predictions from unit-RMS donor directions (`rms_direction`, core334 verbatim):
  `source_mean` and 334-style pooled transfer-gain `diagonal` (median gain 0.55).
  Truth = own line's unit-RMS response. Uncorrected baselines — HepG2: wsign **0.6264**,
  retrieval 0.7808; Jurkat: wsign **0.6015**, retrieval 0.8151
  (wsign = |truth|-weighted sign accuracy, primary; retrieval = 334 cosine-rank formula).
- Subspaces from gene-centred L1000 landmark matrix, k ∈ {10,20,50,100}. Variance
  explained — HepG2: 0.285/0.353/0.465/0.579; Jurkat (256 rows): 0.462/0.551/0.716/0.864.
- Corrections: (a) projection, (b) shrinkage γ ∈ {0.25,0.5}, (c) sign-only hybrid.
  Controls: wrong-line MCF7 subspace (full 28,596 + n-matched), gene-shuffled same-line
  (3 seeds), random orthonormal (3 seeds), held-target same-line (96.8–100% rows kept).
- Statistics: paired 10,000-target bootstrap of Δ vs uncorrected; pooled = inverse-variance.

## Results
Best weighted-sign arm per subspace (mean base; HepG2 uncorrected 0.6264, Jurkat 0.6015):

| line | subspace | best mode | wsign | Δ vs uncorr [95% CI] | retrieval Δ [95% CI] |
|---|---|---|---|---:|---|
| HepG2 | **same** | shrink k50 γ0.5 | 0.6379 | **+0.0115 [+0.0107,+0.0124]** | **−0.0274 [−0.0299,−0.0249]** |
| HepG2 | **wrong (MCF7)** | shrink k50 γ0.5 | **0.6398** | +0.0134 | −0.0298 (shortlist row) |
| HepG2 | held-target | shrink k50 γ0.5 | 0.6380 | +0.0116 | −0.0275 |
| HepG2 | shuffled (best of 3) | shrink | 0.6265 | +0.0001 | ≈0 |
| HepG2 | random (best of 3) | shrink | 0.6265 | +0.0001 | ≈0 |
| HepG2 | same, projection | proj k100 | 0.6228 | −0.0036 | −0.0739 |
| HepG2 | same, sign-only | k100 | 0.6347 | +0.0083 | −0.0118 |
| Jurkat | **same** | shrink k50 γ0.5 | 0.6061 | **+0.0045 [+0.0038,+0.0053]** | **−0.0226 [−0.0247,−0.0204]** |
| Jurkat | **wrong (MCF7)** | shrink k100 γ0.5 | 0.6057 | +0.0042 | −0.0210 |
| Jurkat | shuffled/random | shrink | 0.6025–0.6029 | ≈+0.001 | ≈−0.005 |
| Jurkat | same, projection | proj k100 | 0.5756 | −0.0259 | −0.0533 |

Pooled weighted-sign Δ (same-line best arms): **+0.0075 [0.0069, 0.0080]**.
Diagonal-base arms mirror mean-base to ≤0.0005 (gain is all-positive, signs identical).

**Mechanism (why wrong-line works):** HepG2 vs MCF7 (liver vs breast) L1000 structure is
shared far above chance — mean drug-effect vectors corr **0.524**, top-50 subspace
canonical cosines mean **0.574** (random-expect 0.283). Pure projection *hurts* signs on
both lines, so the subspace directions hold less sign information than the consensus;
only shrinkage (which retains 50–75% of p) gains — a generic pull toward dominant
shared programmes, not line-specific structure.

## Verdict against the PLAN.md kill criterion — DEAD on all three legs
1. Pooled best-configuration gain **+0.0075 < +0.01**. (HepG2 alone is +0.0115 nominally,
   but that is a best-of-16 in-sample maximum; the pooled gate is the registered one.)
2. **Wrong-line matches/exceeds same-line on both lines** (HepG2 +0.0134 vs +0.0115;
   Jurkat +0.0042 vs +0.0045) under identical selection — the gain is not line-specific.
3. **Retrieval falls 2.4× more than signs rise** (−0.0274 vs +0.0115 HepG2; −0.0226 vs
   +0.0045 Jurkat) — the 155/362/363 exchange. Even sign-only (best trade) loses
   retrieval (−0.0118) while gaining +0.0083.
No `work/hela_l1000_subspace.npz` produced (mandate's K562+HepG2 gate untestable by bank gap).

## Limits
- Landmark-only, 623-gene axis; endpoints are source-stage proxies, not six-member scores;
  retrieval is the 346 pds proxy, not board pds. H1-panel numbers are development evidence.
- shRNA subspace, not the mandated compound subspace (values not held; would need ~35 GB
  download). Inference that compounds fare worse is an argument, not a measurement.
- Jurkat same-line subspace from 256 rows (variance explained inflated; reported as secondary).
- Best-arm Δs are in-sample maxima over 16 (k, mode) configs; bootstrap CIs do not include
  selection. Verdict does not depend on it (wrong-line wins under identical selection).
- Held-target control is weak (96.8–100% rows retained) — consistent but unsurprising.

## Exact commands to reproduce
```
uv run python experiments/377-l1000-response-subspace-prior/run377.py base
uv run python experiments/377-l1000-response-subspace-prior/run377.py l1000   # 11 GB sequential read, ~5 s NVMe
uv run python experiments/377-l1000-response-subspace-prior/score377.py
```
Work products: `work/base377.npz` (predictions/truth/axis/targets), `work/l1000_{HEPG2,MCF7,JURKAT}.npz`
(X landmark modz, genes, shRNA target per row), `work/scores377.json` (all arms + bootstraps),
`work/vec_{hepg2,jurkat}.npz` (per-target vectors for pooling). Nothing outside this directory read-for-write; no commits/builders/submissions.

## Peak memory and wall time
Wall: base 3.3 s, l1000 5.2 s, score 20.2 s (+~20 min exploration/planning). Memory guard:
system free 75–78% throughout (`memory_pressure`). Peak RSS not instrumented; bounded by
arrays — largest live: MCF7 28,596×623 float64 (143 MB) + HEPG2 18,727×623 (93 MB) +
gctx chunk 2,000×12,328 float32 (99 MB) + 2,156² cosine blocks (~37 MB); estimated <2 GiB,
well under the 12 GiB cap. Threads capped at 4.
