# 425 — CRISPRi cis-neighbours are real but already in our rows (closed, tiny); the nmae magnitude floor is structured by expression (→ 426)

Status: 2026-09-23 03:00–03:15Z, PI, analytic only (no emission, no builder, no `cell-eval2`). Inputs: 323's h1-test
reference DE (`reference_de.parquet`, the scorer's Wilcoxon on arithmetic-mean CPM), 410's delivered H1 statistics of the
shipped 395 construction (`410/work/h1/stats-base395.npz`, base nmae 0.95210 = 412's), 417a's bulk dump
(`417/work/dump_ka417a_g065.npz`), GENCODE v44 basic gene TSSs (`tss425.py` → `work/tss_v44.tsv`, 61,228 names).

## 1. Cis-neighbours (bidirectional promoters) — real, large, and already carried: CLOSED
dCas9-KRAB at a target's TSS represses genes whose TSS sits within ~1–5 kb (bidirectional partners). On h1-test
(`h1cis425.py`, 85 targets):

| TSS distance from the target | pairs | share significant | share of significant calls DOWN | mean lfc of significant |
|---|---:|---:|---:|---:|
| < 1 kb | 14 | 0.857 | 0.917 | −2.13 |
| 1–5 kb | 11 | 0.727 | 0.875 | −1.30 |
| 5–20 kb | 26 | 0.192 | 0.600 | −0.24 |
| > 20 kb, same chromosome | 48,255 | 0.06–0.17 | 0.54–0.61 | −0.05 … −0.25 |

The biology is clear (e.g. MORN2 → DHX57 −4.6, ANKZF1 → ATG9A −4.3, NDUFB9 → TATDN1 −4.1 log2). It does not pay:
- the < 1 kb bin carries only **0.4 %** of the nmae denominator on average (max 9 % for one target), because the gates are
  hundreds of genes of modest |lfc|;
- **ORACLE** (the cis pairs set to the reference's own lfc): nmae **−0.0022**; the best oracle-free rule
  (`min(ours, −3)` on pairs < 1 kb) **−0.0005** (`h1nmae425.py`). About +0.0001 overall;
- **417a already carries them** (`panelcis425.py`): of the panel's 63 targets with an expressed non-panel neighbour
  within 5 kb, most neighbours are the most negative gene of their bulk row (rank 0–3; e.g. TRIP4 → PCLAF, DCAF16 → NCAPG,
  TBCK → AIMP1 −3.0). The donors are CRISPRi at the same TSSs, so the effect arrives through the rows.

Do not re-propose "model cis repression explicitly" without a new mechanism.

## 2. The finding that matters: nmae's L1 optimum is conditioned on significance, so magnitudes should follow 1/sqrt(count)
nmae is scored only on the reference's significant genes. Conditional on significance, |lfc| is bounded below by the
minimum detectable effect at the gene's control expression, MDE ≈ z·SE ≈ 0.24/sqrt(mean count per cell) (Poisson,
400 vs 18,400 cells). For a gene with sign accuracy a, a predicted magnitude x < |l| changes the expected error by
−(2a − 1)·x: **adding magnitude up to about the MDE is a sure gain whenever a > 1/2**, with no overshoot risk.

On H1 (`mde425.py`, gate genes of the 58 gated targets, 64,931 pairs):

| control mean count/cell | median \|l\| | q10 \|l\| | our median \|l̂\| | our sign accuracy |
|---|---:|---:|---:|---:|
| 0.10–0.31 | 0.624 | 0.371 | 0.034 | 0.520 |
| 1.0–1.6 | 0.266 | 0.180 | 0.056 | 0.546 |
| 3.5–6.1 | 0.161 | 0.107 | 0.030 | 0.587 |
| 6.1–148 | 0.134 | 0.080 | 0.027 | 0.592 |

The q10 of |l| tracks 0.24/sqrt(count) closely. Our delivered magnitudes are ~0.17× the truth and far under the MDE.

| rule on our delivered lfc (signs ours) | Δ nmae (H1) |
|---|---:|
| global × 1.25 (412's best) | −0.0041 |
| constant add sign·0.06 | −0.0066 |
| **add sign·c/sqrt(count), c = 0.1** | **−0.0117** |
| add, c = 0.075 / 0.125 / 0.15 | −0.0103 / −0.0109 / −0.0075 |
| floor \|l̂\| ≥ 0.2/sqrt(count) | −0.0095 |
| exponent scan (add c/count^e, best c each) | e 0: −0.0066 · 0.25: −0.0101 · **0.5: −0.0117** · 0.75: −0.0100 |

- The exponent 1/2 is the Poisson MDE law, and it wins the scan.
- Split-half over targets (`mde425b.py`, 20 splits): c fits at 0.107 ± 0.014, and the held-out gain is −0.0100 ± 0.0050.
  36 of 58 targets improve.
- Split by our call: gate genes we call (38.5 % of the gate, sign accuracy 0.600) give −0.0053; gate genes we do not call
  (sign accuracy 0.536) give −0.0064.

**Why it is implementable as a pure-nmae arm.**
- The K-class emitter's arithmetic-mean steering (`rc_emit_k`, prior mode) sets each gene's arithmetic mean A through its
  class split at a fixed per-gene total (the bulk).
- It designs the Wilcoxon |z| (calls, jac, reach order) after the split, and keeps every sign.
- So an MDE add on the steered targets moves only nmae, like cqawesome's pure-nmae board steps.

Built and measured on the full six-member H1 harness in **426**.

Files: `tss425.py`, `h1cis425.py`, `panelcis425.py`, `h1nmae425.py`, `mde425.py`, `mde425b.py`; outputs in `work/`
(`tss_v44.tsv` is regenerable from the GTF and is gitignored-size-safe at 3 MB).
