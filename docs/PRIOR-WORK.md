# Prior work: what won VCC 2025, and what the literature says works

> **Historical review compiled August 21, 2026.** Section 5's donor counts,
> board comparisons, and suggested priorities are not current. Consult
> `docs/READ-FIRST.md` and the September 9 scientific review before treating an
> old recommendation as untested work.

Compiled 2026-08-21 from primary sources. **Verified** = read from the paper/code/leaderboard.
**Unverified** = inference or second-hand.

## 1. VCC 2025 final leaderboard (100 held-out perturbations)

| Rk | Team | Overall | DES | PDS | MAE | Model |
|---|---|---|---|---|---|---|
| 1 | BM_xTVC (BioMap) | 0.3401 | 0.3488 | 0.8720 | 1.026 | xTVCM_v01 |
| 2 | XLearning Lab | 0.3279 | 0.3555 | 0.8505 | 1.069 | X |
| 3 | Outlier | 0.3259 | 0.3611 | 0.8444 | 4.219 | TransPert |
| 6 | SIMON | 0.3053 | 0.3547 | 0.8175 | 1.111 | **"Naive"** |
| 15 | **Mean Predictors** | 0.2354 | 0.3054 | 0.7412 | 6.723 | final-v1 |

**Derived 2025 baselines** (Arc never published them; obtained by exact-fit inversion of the
public scoring code, residual ≤3e-17 over 100 rows): DES 0.0951, PDS 0.5082, MAE 0.0250.
The baseline model is the perturbation-mean.

**Only 1 of 100 entries beat the MAE baseline.** Everyone else clipped to 0, so 2025 was a
two-metric contest. A team called "Naive" finished 6th; "Mean Predictors" finished 15th and 3rd
on the Generalist board. That is the context for every 2025 method claim.

### What the winners actually did
- **1st, BioMap**: scFoundation/xTrimoGene encoder + protein-LM embeddings + disentangled
  cross-attention decoder. Trained on **pseudobulk**. Explicitly added *training-set DEG
  frequency and mean expression as input features* because, in their words, "purely AI-based
  approaches did not consistently outperform statistical baselines." No code, no paper.
- **2nd, XLearning Lab**: plain **fully-connected net** on pseudobulk over
  [control expression ‖ ESM-2 embedding of the target ‖ UMI indicator], **predicting the delta**.
  Worked out that PDS carried ~2× the weight of DES and optimised for it.
- **3rd, Outlier (TransPert)**: **pure statistics on summary data** — pseudobulk ATE matrices,
  similarity-weighted aggregation over 10 external CRISPRi screens, low-rank combiner, then a
  **global scalar cross-validated to maximise PDS**. Their SVD finding: the dominant rank-1
  component of the effect-magnitude matrix correlates 0.57–0.82 with baseline control expression
  across 7 lines. Much of the apparent signal is gene abundance.
- **Generalist, Altos Labs → PRiMeFlow** (arXiv 2604.13986, Arc explicitly recommends reading it):
  conditional flow matching in gene space, 1D U-Net. **The winning move was not the architecture
  — it was locating the held-out test perturbations in other public cell lines and fine-tuning on
  them.** That is exactly the 2026 task. Also: asymmetric classifier-free guidance (CFG 20 for
  perturbed, 5 for synthetic controls) explicitly to amplify the delta for PDS, and generating
  their own synthetic control cells rather than submitting real ones.

## 2. Why 2026's metrics changed

arXiv 2511.16954 showed that on the 2025 metrics, multiplying a prediction by a scalar `c ≈ 20–30`
raised PDS from 0.52 to ~0.90 **with no change in direction**. Arc's answer was `pds_cosine`,
which is exactly scale-invariant. The largest 2025 exploit is dead.

**Systema** (Nature Biotech 2025) is the intellectual origin of the 2026 design: reference the
*perturbed* centroid rather than the control, and score centroid accuracy. That is precisely
`b` = context-mean, `r` = split-half replicate, `pds_cosine` = centroid accuracy.

## 3. The negative-result literature (this is the important part)

- **Ahlmann-Eltze, Huber & Anders, Nature Methods 22:1657 (2025)** — "Deep-learning-based gene
  perturbation effect prediction does not yet outperform simple linear baselines." Five
  foundation models and two other deep models; **none beat the baselines.** Their baseline is a
  reduced-rank ridge regression on pseudobulk deltas.
- **scPerturBench, Nature Methods 23:451 (2026)** — 27 methods, 29 datasets. On the
  **cellular-context generalisation** task (closest to VCC 2026), over the genuinely
  cross-cell-line datasets: best model **0.608** PCC-delta, **trainMean 0.558**, control-copy
  −0.007, ridge 0.194. **The entire measured margin of deep learning over "average what you
  already have" is ~0.05 PCC-delta.**
- On the unseen-perturbation task, **GEARS (0.369) loses to a linear model (0.422) and to the
  mean baseline (0.402); scFoundation scores 0.032.**
- **"Response Magnitude as a Dominant Signal"** (arXiv 2608.00152): magnitude-only predictors
  transfer *positively* across cell types; expression-only predictors transfer *negatively*.
  Cross-cell-type pretraining on Replogle followed by VCC fine-tuning made things **worse**
  (R² = −0.181). ⚠️ Different prediction target than VCC — read as directional evidence only.

**This corroborates our experiments 004 and 006 from the outside.** We found spectral truncation,
supervised reweighting and magnitude-discarding all lose to the raw donor delta. The literature
finds the same thing at the level of whole architectures.

## 4. Cross-line conservation — the numbers

**Nadig et al. 2025, Nat Genet 57:1228** (TRADE, noise-corrected), median transcriptome-wide
impact correlation over 1,660 essential-gene perturbations:

| pair | r | | pair | r |
|---|---|---|---|---|
| HepG2 ↔ RPE1 | 0.75 | | K562 ↔ HepG2 | 0.64 |
| K562 ↔ Jurkat | 0.74 | | Jurkat ↔ RPE1 | 0.56 |
| Jurkat ↔ HepG2 | 0.69 | | **K562 ↔ RPE1** | **0.40** |

High pairs are biologically similar (K562/Jurkat haematopoietic p53-mutant suspension;
HepG2/RPE1 epithelial p53-WT adherent). **44% of perturbations are cell-type-preferential.**

**The single most important number:** within K562, the median correlation of raw log-fold-change
point estimates between two replicate experiments is **0.16**, while the TRADE-corrected true
correlation is **0.90**. Sampling noise, not biology, destroys the observable correlation.

⚠️ **All of these are measured on common-essential genes, which the VCC 2026 panel excludes
entirely** (our own finding: 0/300 overlap with four independent essential panels). Treat them as
an upper bound of unknown tightness.

## 5. What this tells us to do

1. **More donor lines.** X-Atlas/Orion (HCT116 + HEK293T, genome-wide, 8M cells) covers 271/300
   targets *and* 25 of the 28 K562-gwps misses. Takes us from 1 usable donor line to 3.
2. **Per-perturbation magnitude calibration** — our experiment 001/002 thesis, confirmed
   un-exploited: all 65 teams on the live board have raw mse > 1.
3. **Emission model before architecture.** Four of six metrics are computed from the cells we
   emit. Sidechain (public code, rank 5) measured an even-spread null at −0.01 against −0.30 for
   a resampled-control null.
4. **Shrink before transferring.** Nadig's 0.16-vs-0.90 gap says raw source log2FCs are mostly
   noise. Note this is *not* the same as the low-rank truncation that failed in our exp 004 —
   empirical-Bayes / positive-part shrinkage preserves the full spectrum.
5. **Weight donors by basal-profile similarity to the target.** The only context signal we get.
6. **Do not build a large end-to-end neural model first.** Three benchmarks and two competition
   years say the margin over "average what you know" is ~0.05.

## Sources

PRiMeFlow [arXiv 2604.13986](https://arxiv.org/abs/2604.13986) · [github.com/altoslabs/primeflow](https://github.com/altoslabs/primeflow)
PDS scaling [arXiv 2511.16954](https://arxiv.org/abs/2511.16954)
Ahlmann-Eltze [Nature Methods](https://www.nature.com/articles/s41592-025-02772-6) · [code](https://github.com/const-ae/linear_perturbation_prediction-Paper)
Systema [Nature Biotech](https://www.nature.com/articles/s41587-025-02777-8) · [code](https://github.com/mlbio-epfl/systema)
scPerturBench [code](https://github.com/bm2-lab/scPerturBench)
Nadig 2025 [PMC13063516](https://pmc.ncbi.nlm.nih.gov/articles/PMC13063516/)
X-Atlas/Orion [HF](https://huggingface.co/datasets/Xaira-Therapeutics/X-Atlas-Orion) · [bioRxiv](https://www.biorxiv.org/content/10.1101/2025.06.11.659105v1)
Stack [github.com/ArcInstitute/stack](https://github.com/ArcInstitute/stack) · [weights](https://huggingface.co/arcinstitute/Stack-Large-Aligned)
Sidechain (live rank 5, public) [github.com/saberhq/sidechain](https://github.com/saberhq/sidechain)
Turbine AI writeup [blog](https://blog.turbine.ai/p/how-did-we-get-a-regression-model)
