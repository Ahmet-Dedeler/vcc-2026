# 388 — The decoupled DE channel's H1 prize is amplitude only (≈ −0.02 raw nmae)

Status: **COMPLETE** (E4 ruled out by construction — see §What was measured).
Verdict: **kill criterion FIRED on both clauses**: no estimator beats E0's |l|-weighted
sign accuracy by ≥ +0.03 (best is E1 at +0.012), and the best magnitude map over all
estimators gains 0.024 raw nmae (< 0.03). **The decoupled channel's H1 prize is amplitude
only (≈ −0.02 raw nmae).** The mechanism itself is feasible (gate §4), but no priced
estimator has better signs to put through it.

## Known-answer reproduction (rule 3) — all match

Own code (`ka388.py`), 343's conventions, 58 eligible targets:

| row | here | KA | Δ |
|---|---:|---:|---:|
| mean NMAE | 0.950836 | 0.950836 | 0 |
| \|l\|-weighted sign acc, delivered | 0.587968 | 0.588 | 3e-5 |
| \|l\|-weighted sign acc, designed | 0.600185 | 0.600 | 2e-4 |
| cf every-sign-correct | 0.803387 | 0.803 | 4e-4 |
| cf best scalar/target (343 grid to 6×) | 0.927608 | 0.928 | 4e-4 |
| cf both (343's exact LS formula) | 0.658430 | 0.658 | 4e-4 |
| O1 reach (verbatim `walk_reach`) | 0.25434726 | 0.25434726 | <1e-8 |

Note: KA "both 0.658" is `decompose343.py`'s exact `cf_signs_and_scalar` (least-squares
scalar + correct signs). Grid re-optimized scalar + correct signs gives 0.654063 (extra row
in `work/ka388.json`). parquet-vs-npz reference agreement exact (max|dq| = 0, max|dlfc| = 0).

## What was measured (inputs, estimator, controls)

- Estimators as signed 85×10779 fields: **E0** delivered lfc; **E1** designed cap15
  (2 universe genes absent from the effects axis → 0, signs343 convention); **E2**
  un-whitened consensus `h1test-feats.npz::cons` — provenance verified in source
  (`probe172.py: SOURCES = (k562_gwps, orion_hek293t)`, two donors; axes verified
  name-identical to the 323 panel); **E3** = E1 magnitudes + E2 signs where they disagree
  (29.0% of entries flipped, 173 analogue); **E4 N/A** — `donorstats.*.npz` carry only
  gene-level stats (`m/marg_abs/marg_sd/resp_frac/scale`, 8248 genes), no signed
  per-(target,gene) rows; rebuilding needs GB-scale `delta.npy` reprocessing, outside this
  assignment's analytic scope; **E5** deranged E2 (seed-388 permutation, no fixed points);
  **E6** random signs × |E0|.
- Gate: 343's (ref q<0.05, finite nonzero lfc, own target out, ≥10 → 58 targets).
- nmae maps: estimator's own magnitudes; estimator signs × E0 magnitudes; optimal maps
  given signs — one global scalar, one scalar per |est| decile (pooled fit), one scalar
  per target (343 grid to 6×) + wide-grid (to 24×) diagnostic.
- Reach: `walk_reach` copied VERBATIM from `reach348.py:125-145`; (i) delivered p-value
  order + estimator signs, (ii) |estimator| order + estimator signs. 82/85 targets scored.
- Bootstrap: paired-target, 10,000 draws, seed 20260913 (342 convention).

Falsifier checks: E5/E6 sign accuracy 0.479/0.498 ≈ 0.5 ✓ (gate code sound); E3 sign
accuracy 0.5852 lies between E1 (0.6002) and E2 (0.5785) ✓ (flip code sound); E6
nmae_own == nmae_e0mag (1.0087) ✓ (identity: E6 magnitudes ARE |E0|); E0 per-target-map
nmae 0.9276 == KA cf_scalar 0.927608 ✓ (same computation).

## Ranked table (mean over 58 eligible targets; reach over 82 scored)

| est | sign acc gate | sign acc call-set (w) | nmae own | nmae best map | fid proxy (unwtd prec on calls) | reach (deliv order) | k* | reach (|est| order) |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| E0 delivered | 0.5880 | 0.6377 | 0.9508 | **0.9269** | 0.5735 | 0.25435 | 13.8 | 0.26815 |
| E1 designed | 0.6002 | 0.6493 | 0.9701 | 0.9317 | 0.5868 | 0.24327 | 13.7 | 0.26837 |
| E3 flip | 0.5852 | 0.6554 | 0.9721 | 0.9387 | 0.5857 | 0.20914 | 12.2 | 0.21461 |
| E2 consensus | 0.5785 | 0.6524 | 0.9712 | 0.9382 | 0.5833 | 0.20900 | 12.0 | 0.20415 |
| E6 random | 0.4979 | 0.4869 | 1.0087 | 0.9922 | 0.4837 | 0.09063 | 1.2 | 0.07459 |
| E5 deranged | 0.4792 | 0.4825 | 1.0036 | 0.9818 | 0.4946 | 0.08344 | 1.0 | 0.06758 |

(E4 N/A — see above. Best map = min over global/decile/per-target/wide maps.)

## Deltas vs E0 (paired bootstrap, 10k draws, seed 20260913)

| est vs E0 | Δ sign acc (+ better) | Δ own-nmae, E0−X (+ better) | Δ reach, delivered order (+ better) |
|---|---:|---:|---:|
| E1 | +0.0122 (se 0.0070, P=0.952) | −0.0192 (se 0.0040, P=0.000; E1's raw amplitude is worse) | −0.01108 (se 0.00879, P=0.087) |
| E3 | −0.0028 (se 0.0118, P=0.413) | −0.0213 (se 0.0052, P=0.000) | −0.04521 (se 0.01843, P=0.002) |
| E2 | −0.0094 (se 0.0130, P=0.235) | −0.0203 (se 0.0062, P=0.000) | −0.04535 (se 0.01843, P=0.002) |
| E5 | −0.1087 (se 0.0226, P=0.000) | −0.0527 (se 0.0083, P=0.000) | −0.17091 (se 0.03348, P=0.000) |
| E6 | −0.0900 (se 0.0200, P=0.000) | −0.0579 (se 0.0095, P=0.000) | −0.16371 (se 0.03546, P=0.000) |

Supporting numbers: E1 signs × E0 magnitudes → nmae 0.9485 (gain of 0.0023 over E0's
0.9508 — the sign edge is real but worthless without amplitude); global-d
...[truncated 3013 chars]