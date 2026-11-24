# 175 — per-donor sign accuracy on H1: p53 hypothesis REFUTED; the Orion pair adds almost no sign information over K562 alone

Date: **2026-09-02**, coordinator-run, measurement only (`perdonor175.py`, ~1 min). h1-test reference
DE pool (`164/de164.a15.npz`), donor tables `data/donors/{k562_gwps, orion_hct116, orion_hek293t}`
mean-centred per gene as production does (`DonorHandle.rows`), no scaling, no whitening.

## Hypothesis tested
Top-tier sign errors are target-level (SBDS 0.41, SOX4 0.46, CYCS 0.47 on their top-2% genes).
SBDS knockdown acts through p53; K562 is TP53-null, HEK293T's p53 is inactivated by large-T,
HCT116 is TP53-WT, H1 hESC is TP53-WT. So the p53-competent donor should transfer better on
p53-mediated responses, and pathway-status matching would be a mechanism for donor weighting.

## Result: refuted, and something else fell out

**On 374 pool pairs over 31 canonical p53 target genes:** K562 0.526, **HCT116 0.485**, HEK293T
0.447, equal mean **0.444**. The p53-WT donor is *worse* there than the TP53-null one, and the
consensus is below chance. p53 status does not explain the target-level failures.

**Per-donor sign accuracy on the pool (all three donors carrying, n = 39,945):**

| source | sign acc | median `\|s\|` on pool |
|---|---:|---:|
| **K562 alone** | **0.5670** | 0.024 |
| HCT116 alone | 0.5375 | 0.041 |
| HEK293T alone | 0.5308 | 0.037 |
| **Orion pair only (0,1,1)** | **0.5425** | |
| equal mean (1,1,1) — production's weights | 0.5697 | |
| (2,1,1) | 0.5740 | |
| (3,1,1) | 0.5742 | |
| majority vote | 0.5633 | |
| unanimous pairs only (30% of pairs) | 0.6213 | |
| **per-target ORACLE best single donor** | **0.6032** | routing ceiling |

**The Orion pair carries two-thirds of the consensus weight and lifts sign accuracy over K562 alone
by +0.0027.** Its effects are 1.7× larger in magnitude and less accurate — bigger and noisier, so in
an unweighted mean they dominate. Up-weighting K562 to 2–3× buys +0.004; the per-target oracle
router would buy +0.027 (066 closed routing for `pds`; this is its ceiling on sign).

Per target, HCT116 wins 19 / loses 29 against the other two's mean; K562-alone vs equal-mean is
22/48 with median +0.0003 — the asymmetry is not concentrated in a few targets.

## Reading
This is the `pds`-vs-sign trade a third time. 019 established that equal donor weights are best for
`pds` under whitening (the "4:1" K562-heavy configuration lost). For per-gene sign accuracy the
K562-heavy weighting wins. Same profile, two objectives, opposite optima — exactly as with gamma
(097/170), sign repair (173) and whitening itself (172). Composable rider at most; not a lever.

Not pursued here: production applies a per-donor scale (`MultiDonorModel._scale`) that this
unscaled analysis omits; the (2,1,1) row approximates its effect. Whether Orion's larger effects
are noise or genuine line-specific biology is not separable on this surface.

## Addendum (19:45) — per-donor rank-1 deflation: the Orion pair shares a protocol axis

`deflate175.py`. Per-donor PCs fitted on 4,000 of each donor's own perturbations. **HEK293T's top
direction carries 19.7% of its variance (K562 5.6%, HCT116 8.6%), and the two Orion lines' top
directions have cosine 0.872 with each other while K562's is at 0.23–0.31 with either** — a shared
X-Atlas platform axis, the donor-table face of the 26–44% harness blob. Deflating exactly one
component per donor lifts sign accuracy on the H1 pool for every donor and the consensus (equal mean
0.5697 → **0.5770**, (2,1,1) 0.5740 → 0.5782); k ≥ 2 removes signal. Interior optimum at k = 1 with a
shared, K562-absent direction = one artefact on top of a full-rank signal (consistent with 004).
Handed to 178 as a fourth composable switch (pre-consensus, per-donor — not 156's post-prediction
deflation).

## Addendum 2 (23:20) — the long-gene axis is NOT a batch-control artefact; a within-batch rebuild would not remove it

Checked at the raw Orion parquet level (HEK293T, 223 batch files, `gene_target`/`guide_target`/
`sample` per cell). Batch control pseudobulks DO sit at very different points on the axis (projection
sd **4.57**, range −12.5 to +4.7) — but a perturbation's ~154 cells are spread over **~115 batches at
~1.3 cells per batch** (median 2 cells/target/batch, 1 guide/target/batch), so `build_orion_donor.py`'s
global-control delta averages the batch offsets down to ~0.43, well below the observed per-perturbation
axis sd of **3.20**. The artefact is therefore a genuine per-perturbation shift of the long-gene /
nuclear-fraction axis in the Orion cells (plausibly a toxicity or cell-state readout the X-Atlas
libraries express and Arc's H1 data does not), not a batch-control mismatch. **A within-batch rebuild
of the donor tables would not remove it; rank-1 deflation (above) and production's whitener do.**
