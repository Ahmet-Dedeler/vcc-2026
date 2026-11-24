# 347 — Confident strata exist by accuracy, but amplitude cannot buy NMAE in any of them (overshoot binds)

Status: **NEGATIVE** (route closed). Descriptive analysis on saved cap15 arrays only; no 323 run, no builder.

## Known-answer reproduction (own code, `strata347.py`)
Mean over the same 58 eligible targets, 343's exact conventions:

| quantity | reproduced | 343 recorded | match |
|---|---|---:|---|
| mean NMAE | 0.9508364630325662 | 0.9508364630325662 (`metrics-cap15.json` raw) | exact, all 16 digits |
| delivered \|l\|-weighted sign acc, all / called / uncalled | 0.587968 / 0.637666 / 0.551763 | 0.587968 / 0.637666 / 0.551763 | 1e-12 |
| designed w-acc, all / called / uncalled | 0.600185 / 0.649299 / 0.558007 | 0.600185 / 0.649299 / 0.558007 | 1e-12 |
| gate mass in own significant set | 0.346388 | 0.346388 | 1e-12 |
| eligible targets / pooled gate pairs | 58 / 64,931 | 58 | yes |

(Pooled delivered magnitude ratio 0.206 vs 343's mean-over-targets 0.205: convention difference, expected.)

## What was measured
Pooled gate-pair table (58 targets × reference-significant genes, own target excluded) with ref lfc `l`,
delivered lfc `lh` (amplitude proxy, as mandated), delivered p/q, design effect `ef` (effects.npz::cap15,
2 universe genes absent from the effects axis zeroed), ref_lin basal CPM. Stratum stats pooled
(n, \|l\|-mass share, DESIGN-weighted/unweighted accuracy primary, delivered secondary, delivered magnitude
ratio); exact per-target NMAE recomputes (mean over targets) scaling delivered lfc in-stratum ×{2,3} with
rest ×{1,0.5}. (d) SKIPPED as provided: `stage.npz` holds only consensus arms
(uniform/plain15/plain30/cap15/parity15) + CD4 aggregates — per-source rows for k562_gwps /
orion_hct116 / orion_hek293t are never saved (run323.py loads donors, stages consensus, `del donors`);
`cd4_votes` ∈ {0,3} only (802,841 × 0, 772,464 × 3): a carried-by-all-or-none flag across the 3 CD4
conditions, not a graded sign-agreement count. No donor-count stratification is recoverable.
(c) common vector is leave-one-out (mean of the other 84 design rows); `common≈0` = |common| ≤ pooled
10th percentile (0.000527). (e) emitter-2600 set reconstructed per target from the stated rule
(rc_emit.py:839,916-919): head = top-400 universe genes by |design effect| (the |H| key), rest = top-2200
by ref_lin among non-head.

## Results — strata (pooled; full table `work/strata-cap15.csv`)
Design |effect| deciles (per-target rank shown; global deciles near-identical):

| stratum | n | mass share | design w-acc | deliv w-acc | mag ratio |
|---|---:|---:|---:|---:|---:|
| a_pt_dec0 (weakest) | 6,511 | 0.105 | 0.459 | 0.497 | 0.046 |
| a_pt_dec1–dec7 | ~6,490 ea | ~0.094–0.102 | 0.531–0.583 | 0.513–0.576 | 0.057–0.299 |
| a_pt_dec8 | 6,491 | 0.094 | 0.613 | 0.602 | 0.376 |
| **a_pt_dec9 (strongest)** | **6,519** | **0.117** | **0.741** | **0.721** | **0.507** |
| a_glob_dec9 | 6,494 | 0.123 | 0.747 | 0.724 | 0.516 |

| stratum | n | mass share | design w-acc | deliv w-acc | mag ratio |
|---|---:|---:|---:|---:|---:|
| b_called (q<.05) / uncalled | 24,997 / 39,934 | 0.253 / 0.747 | 0.625 / 0.550 | 0.626 / 0.540 | 0.346 / 0.159 |
| b_pdec9 (most signif. delivered) | 6,494 | 0.095 | 0.768 | 0.768 | 0.574 |
| c_agree / c_disagree / common≈0 | 29,339 / 29,098 / 6,494 | 0.446 / 0.443 / 0.111 | 0.574 / 0.574 / 0.528 | 0.561 / 0.566 / 0.548 | 0.211 / 0.208 / 0.181 |
| e_ribo / mito / histone | 0 / 107 / 163 | 0 / 0.003 / 0.004 | — / 0.582 / 0.356 | — / 0.582 / 0.512 | — / 0.188 / 0.212 |
| e_head400 / cpm-rest2200 / outside2600 | 15,703 / 39,149 / 10,079 | 0.236 / 0.515 / 0.249 | 0.658 / 0.536 / 0.551 | 0.635 / 0.535 / 0.548 | 0.375 / 0.162 / 0.139 |
| **f_top3_deciles ∩ agree** | **8,911** | **0.138** | **0.651** | **0.637** | **0.405** |
| f_top3 alone | 19,495 | 0.305 | 0.653 | 0.640 | 0.402 |

## Results — exact rescaling counterfactuals (mean NMAE over 58 targets; base 0.95084)

| stratum scaled | ×2 (rest ×1) | ×3 (rest ×1) | ×2 (rest ×0.5) |
|---|---:|---:|---:|
| a_pt_dec9 | 0.9620 (**+0.0112**) | 1.0024 | 0.9695 |
| f_top3_and_agree | 0.9581 (+0.0073) | 0.9769 | 0.9747 |
| e_head400 | 0.9653 (+0.0145) | 1.0318 | 0.9671 |
| b_called | 0.9673 (+0.0164) | 1.0231 | 0.9717 |

Every rescaling of every gate-passing stratum makes NMAE worse, including with rest-shrunk.

## Mechanism — why accuracy ≥ 0.65 is not enough (188's identity, measured not inferred)
NMAE improves under scale-up iff e = ΔM/ΔU > ½ with M = Σ min(|h|,|r|)·1[sign ok]. Pooled at ×2:

| stratum | delivered sign acc (mass-wtd) | overshoot share of stratum \|lh\|¹ | e = ΔM/ΔU |
|---|---:|---:|---:|
| a_pt_dec9 | 0.721 | 0.340 | **0.341** |
| f_top3_and_agree | 0.726 | 0.330 | 0.337 |
| e_head400 | 0.635 | 0.376 | 0.380 |
| b_called | 0.626 | 0.377 | 0.380 |

¹ wrong-sign \|lh\| mass + correct-sign mass already above \|r\|. The confident stratum is exactly where
delivered magnitudes are already largest (mag ratio 0.37–0.57), so the min(|h|,|r|) cap binds and a third
of the added U is unmatched by construction. q_w > ½ is necessary but not sufficient; e > ½ is the gate
and it fails everywhere (0.34–0.38). The top-decile accuracy is also broad-based, not a few-target
artefact: 44/58 targets have top-decile design w-acc ≥ 0.65 (median 0.736, q10 0.588; mass-share median
0.115, range 0.062–0.422) — a real stratum, just not a buyable one.

## Verdict against the PLAN.md kill criterion
Literal reading: the kill condition ("no stratum with ≥ 10% mass reaches weighted sign accuracy ≥ 0.65")
is NOT met — a_pt_dec9 (11.7%, 0.741), e_head400 (23.6%, 0.658) and f_top3_and_agree (13.8%, 0.651) pass it.
Functional reading (the route's purpose): CLOSED anyway — the exact counterfactuals show scaling any of
them buys nothing (best +0.0073 raw worse), with the measured mechanism above. **Confidence-weighted
amplitude is closed; do not run a 323 arm on it.** Side negatives worth the record: (c) design/common-vector
agreement is exactly coin-flip informative (agree 0.574 vs disagree 0.574 design w-acc at 45%/44% mass —
the shared programme carries zero reliability information); (e) ZERO ribosomal (RPL/RPS) genes sit on the
reference gate (mito 107, histone 163 pairs — no mass to work with).

## Limits
H1 development panel only (58 targets, one held-out line); not an A/B/C forecast. Counterfactuals scale
delivered lfc as the amplitude proxy — a 323 emission run would re-realise noise but the sign/overshoot
arithmeti
...[truncated 615 chars]