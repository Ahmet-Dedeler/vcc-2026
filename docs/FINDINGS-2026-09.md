# Findings — perturbation-response transfer (v3.1, F16 re-graded by 320: NEW 0 / EXTENDED 1 / KNOWN 4 — a re-measurement, not a discovery; was v3.0, 2026-09-07; Round 6: F16 — the seed-match repression signal, construct as unit, C8 aborted and C9 PASS; 308 closes donor quality weighting; 306 finds no acquirable data)

**Current interpretation overrides the dated narratives below.** 304's audit
and `experiments/304-f15-audit/COORDINATOR-REVIEW-2026-09-07.md` withdraw
F15's biological orthogonality/latent-mirror bound, ORF depth-law extension,
isoform-as-noise-ceiling reading and target-specific geometry claim. 305's
`COORDINATOR-REVIEW-2026-09-07.md` records two nominal primary repression
signals but **zero full-family discoveries**. Neither 301 nor 305 excludes
off-targeting as a mechanism. The 19–27% library-specific failure summaries
are not a measured frequency of off-target responses in 2026 challenge
truths. No competition setting changes; R4 remains spent and C6 remains VOID.

Written by experiment 275 from the RESULTS.md of 260–273 and the two board points of 2026-09-05. Every number below was
re-read from the saved JSON/CSV it came from (`experiments/275-findings-writeup/work/verify275.log`, 447 numbers).
Key: R<n> = `experiments/<n>-*/RESULTS.md`, BR<n> = `experiments/<n>-*/BOARD-RESULT.md`. Intervals are 95 % paired
target bootstraps [lo, hi] unless stated. `pds` is the scorer's per-target rank statistic (§5.1); "own cosine" is the
cosine of a prediction to its own truth; `cos_agg` the energy-weighted aggregate cosine.

## 1. Setting

The task (AGENTS.md §1–2): predict the post-knockdown expression of 300 CRISPRi targets, 400 cells each, in three held
cancer lines (contexts A/B/C) from their unperturbed cells and public screens in other lines. The incumbent (`cd4svd5`,
BR249) is an equal-weight mean of three donor screens (K562-GWPS; X-Atlas/Orion HCT116 and HEK293T; 150–200 cells per
target), whitened by an operator fitted on those donors, plus a CD4 T-cell residual and an H1 hESC residual. Four
instruments: (i) the board, the real six-member scorer, deterministic to 16 significant figures, two entries per day —
a claim is *board-validated* only when a prediction written before the entry existed matched it, and a *fitted
calibration* is a retrodiction with factors fitted on the point itself; (ii) PLOLO (R162), leave-one-library-out over
the panel with the held line's own 150–200-cell pseudobulk as truth (K562, HCT116, HEK293T held; an iPSC line apart);
(iii) the h1-test instrument (R273), 85 targets × 400 cells of H1 hESC on 10x Flex, scored with the scorer's own
statistic (cell-eval2's per-target `pds` reproduced to 5.6e-17); (iv) Broad-16 (R268), 82 Cas9 knockouts × 12 cancer
lines as pseudobulk halves. Two tracks (`docs/RESEARCH-PROGRAMME-2026-09.md`): the competition track is judged by the
board alone; the scientific track by held experiments, split-half reproducibility and preregistered confirmations on
reserved data (`docs/CONFIRMATION-LEDGER.md`).

## 2. Findings

### F1. A 150–200-cell Perturb-seq response on this panel is mostly not reproducible; the challenge's own Flex-400 reference is

**Claim.** On the non-essential panel a 10x/Orion response at 150–200 cells reproduces its own other half at a median
residual cosine of 0.03–0.05, with a third of the Orion targets at or below zero; an H1 Flex response at 400 cells
reproduces at 0.24 with no target at or below zero. (The clause "and the difference is not depth" is withdrawn: F10's correction leaves the cell axis unmeasured until 285.)

**Provenance (287, novelty audit).** F1a is a RE-MEASUREMENT, not a new claim: X-Atlas/Orion's own preprint reports 48 / 40 / 54 % of on-target perturbations (HCT116 / HEK293T / K562) with a detectable single-cell change, Replogle 2022 31.1 % with a significant phenotype, and Wang et al. 2026 (bioRxiv 2026-08-11) 65 % of 7,170 perturbations across 29 datasets unreliable by split-half Pearson with Spearman–Brown — the same construction on a superset of our data. What is ours is the conditioning (residual cosine, the challenge's 300 non-essential targets, a stated cells-per-half). Orion's headline median Spearman 0.99 is NTC profile stability, not response reproducibility. F1c (the construct) is published by the CD4 source: `guide_correlation_all` medians +0.06 / +0.07 / +0.06 across all genes, attributed to differential knockdown efficiency between guides (E287 `work/checks287.log`).

| reference (source) | chemistry | split unit | cells / half | targets | median split-half residual cosine | p10 / p90 | share ≤ 0 | Spearman–Brown to full depth |
|---|---|---|---:|---:|---:|---|---:|---:|
| K562-GWPS (R267 §3) | 10x 3′ | gemgroup parity | 81 | 269 | +0.047 | +0.008 / +0.162 | 5.6 % | 0.09 |
| HCT116 (R267 §3) | Orion | batch parity | 70 | 290 | +0.028 | −0.049 / +0.119 | 32 % | 0.05 |
| HEK293T (R267 §3) | Orion | batch parity | 95 | 294 | +0.028 | −0.122 / +0.214 | 34 % | 0.05 |
| HEK293T, random cells (R267 §3 v) | Orion | cells | 95 | 294 | +0.071 | — | 15 % | — |
| H1 h1-test (R273 §2) | 10x Flex | lane parity | 200 | 85 | +0.243 | +0.067 / +0.748 | 0 % | 0.39 median / 0.50 mean |
| H1, depth-matched (R273 §2) | 10x Flex | random cells | 80 | 85 | +0.114 | — | 0 % | — |
| guide construct P1 vs P2 (R273 §2; R267 §3) | Flex; Orion / 10x | construct | — | 5; 8–10 | +0.014; +0.010 / −0.001 / +0.004 | — | 2/5; 50 % | — |
| donor rows K562 / HCT116 / HEK293T (R271 §1.2) | as above | batch parity | 69–96 | 260–289 | +0.048 / +0.027 / +0.028 | — | 6 / 32 / 34 % | 0.09 / 0.05 / 0.05 |

| reliability is a target property (R267 §3, R273 §0.5) | K562 | HCT116 | HEK293T | Flex |
|---|---:|---:|---:|---:|
| median reliability, 262's anti-aligned / marginal / saturated groups | 0.035 / 0.053 / 0.121 | 0.014 / 0.025 / 0.107 | 0.036 / 0.032 / 0.110 | low tercile ≤ 0.118 |
| Spearman(reliability, incumbent per-target `pds`) | +0.58 | +0.23 | +0.17 | +0.62 |
| Spearman(reliability, cells per target) | +0.16 | +0.09 | +0.12 | −0.05 |

**Limit.** Halves of one experiment (batches, cells, constructs), never studies; three cancer lines and one embryonic
line; 85 Flex targets; the split unit matters where batches carry structure (HEK293T 0.028 by batch, 0.071 by random
cells). The quoted "0.481 measured at 400 cells" (R159 §1.3; repeated in R271 §5 and READ-FIRST §6) is a linear-law
extrapolation (0.478 in 159's saved fit); the measured value is 0.330 at 200/200 cells.

**What it changes.** Per-target readings on PLOLO's references are noise-dominated below ≈ 0.05 cosine; only totals and
group aggregates carry signal. The Flex-400 reference resolves every target's direction, so a deficit read on it is the
prediction's (F4). The guide construct is a unit of variance the challenge truth (one construct per target) does not
average over — and F14 (285) measures it as the DOMINANT irreducible term on one matched design (ρ_guide 0.09–0.41 for CRISPRi, batch free, a different study ρ 0.79).

### F2. The donor-selection "oracle" is the reference's own noise, selected twice

**Claim.** Choosing per target the donor row closer to the held truth (+0.046 `pds` on PLOLO's deep truths, R263 §3;
+0.0459 [+0.0382, +0.0539] reproduced in R271 §3.1) keeps none of its advantage when the choice is made on one
independent half of the reference and scored on the other.

| arm (R267 §4; batch halves, three cancer folds, 259 two-donor targets per fold, symmetric A→B / B→A) | Δ`pds` vs the equal-weight mean |
|---|---|
| in-sample oracle (chosen and scored on the same half) | +0.056 [+0.050, +0.063] |
| cross-validated oracle (chosen on the other half) | −0.018 [−0.027, −0.009], negative on 3/3 folds both ways |
| random pick (3 seeds) | −0.022 [−0.028, −0.017] |
| cross-validated minus random | +0.004 [−0.005, +0.013], P(≤ 0) 0.17 |
| retained fraction, cross-validated / in-sample | −0.32 [−0.49, −0.16]; random-cell halves −0.25 [−0.41, −0.10] |
| abstaining oracle (best of donor 1 / donor 2 / mean), cross-validated | −0.010 [−0.018, −0.002] |
| reliable tercile (n 285): cross-validated minus random | +0.021 [+0.007, +0.035] (batch), +0.027 [+0.014, +0.040] (random); vs the mean −0.002 [−0.015, +0.011] |
| anti-aligned group (n 149, 48 % of the in-sample gain): cross-validated / random | +0.038 / +0.040 |
| agreement of the two halves' picks (K562 / HCT116 / HEK293T; chance 50 %, break-even 66–70 %, R263 §5) | 56 / 58 / 53 % |
| build-time arbiters, 78 cells (R263 §4); agreement with the oracle | all negative; 44–61 % |
| three-donor majority rule on the iPSC fold, vs the incumbent (R265 §3) | +0.0121 [−0.0035, +0.0278], P(≤ 0) 0.064; permutation p 0.011, z +2.3 on a null SD 0.0068 (R263's "8 σ"); proper random-drop control +0.0068 |
| the same rule on the three cancer folds (R264, R265 §4) | permutation p 0.68 / 0.12 / 0.50; vs the incumbent −0.0107 [−0.0187, −0.0030] |

**Limit.** Three cancer folds, one instrument, halves of one experiment; the HEK293T half-B truths are reserved and were
scored once, for this claim (`docs/CONFIRMATION-LEDGER.md`).

**What it changes.** Selection cannot precede estimation on this data — precisely: F2 rejects per-target selection of a donor row by its proximity to the truth; it does not reject reliability-based filtering of TRAINING perturbations, which Wang et al. 2026 report as a positive (287, T5). Labels derived from 150–200-cell truths
inherit the noise: a utility router trained on them sits at its shuffled controls on held context and family
(+0.0001 [−0.0024, +0.0028], R269 §7). The only information is in the reliable third, where averaging already ties.

### F3. The whitening operator has learnt its training rows' noise; retrieval needs those directions

**Claim.** In the fitted eigenbasis of every operator on disk, an independent estimate of the donors' sampling noise
exceeds the fitted eigenvalue on 77–79 % of the directions, those directions carry 87–92 % of the whitened prediction,
and removing or softening them loses `pds` on PLOLO; on the reliable reference the soft correction is null on `pds`
and positive on cosine.

| operator (basis rows) (R272 §1.2) | p | directions with noise > eigenvalue | whitened prediction energy there | Σ noise / Σ eigenvalue |
|---|---:|---:|---:|---:|
| held K562 (HCT116 + HEK293T rows) | 17,806 | 76.9 % | 87.5 % | 0.95 |
| held HCT116 (K562 + HEK293T) | 17,942 | 77.5 % | 87.7 % | 0.88 |
| held HEK293T (K562 + HCT116) | 17,942 | 78.3 % | 87.9 % | 0.90 |
| production (three donors) | 17,942 | 78.7 % | 91.5 % | 1.00 |

| split-half cosine of a donor's batch halves after whitening (R271 §1.2, R272 §2.2) | K562 fold | HCT116 fold | HEK293T fold |
|---|---|---|---|
| through the operator fitted on both halves (in-sample), donor 1 / donor 2 | −0.14 / −0.14 | −0.14 / −0.24 | −0.13 / −0.24 (97–100 % ≤ 0) |
| the held line's halves through the same operator (out-of-sample) | +0.074 | +0.041 | +0.029 |
| through an operator fitted on the other batch-half (cross-fitted) | +0.011 / +0.009 | +0.022 (9 % ≤ 0) / +0.007 (22 %) | not refitted |

| un-learning arm (Δ`pds` vs the incumbent) | PLOLO, 3 cancer folds (R272 §1.3, §2.3) | Flex-400, 85 targets (R273 §4.1–4.2) |
|---|---|---|
| Wiener damping, κ = 1 (keeps 7 % of whitened energy) | −0.034 [−0.044, −0.024]; beats a shuffled profile by +0.032 [+0.015, +0.048] | −0.013 [−0.035, +0.012]; beats shuffled by +0.050 [+0.007, +0.093]; `cos_agg` −7.7 % |
| hard cut (noise-dominated directions removed) | −0.069 [−0.085, −0.054] | — |
| noise floor c = 1 (keeps 68 % on production) | −0.0022 [−0.0040, −0.0004]; vs a mirrored profile +0.0003 [−0.0005, +0.0011] | −0.0001 [−0.0053, +0.0043]; own cosine +0.0024 [+0.0014, +0.0034], `cos_agg` +4.4 %; without the CD4 residual +0.0043 [−0.0010, +0.0099] |
| cross-fitted operators averaged (HCT116 / K562 fold) | +0.0001 [−0.0033, +0.0035] / +0.0003 [−0.0038, +0.0045] | — |
| proper cross-fit (each half through the operator that never saw it) | +0.0003 [−0.0050, +0.0056] at `cos_agg` −3.3 % / +0.0022 [−0.0034, +0.0078] at −3.7 % | — |
| donor rows re-weighted by measured reliability, BLUE weights (R271 §3) | −0.0037 [−0.0103, +0.0028] | — |
| EB row shrinkage toward the consensus (R271 §3–4, R273 §4.2) | +0.0014 [+0.0005, +0.0024]; over its shuffled control +0.0013 [−0.0001, +0.0027]; retention 31 %; priced +0.00050 pessimistic | 0.0000 [−0.0022, +0.0025] |

Provenance (287): the eigenvalue-bias half of F3 is textbook (Ledoit–Wolf 2004; Marchenko–Pastur; Gavish–Donoho 2014) and is cited as such; what is NEW is F3b — that the standard shrinkage-whitening prescription LOSES rank retrieval here, against the embedding-retrieval literature's recommendation, because the objective differs (300-way scale-invariant rank retrieval with the discriminating signal in the low-variance directions; 287 T1).

Why: the target-specific signal is 5–10 % of a donor row's energy (norm moments 0.09 / 0.07 / 0.07, R271 §1.2; the
half-difference identity 0.105 / 0.054 / 0.064, R272 §1.1) and it is spread over the low-variance directions the
operator amplifies; damping them by any profile removes the discrimination before the noise.

**Limit.** The noise covariance comes from ~270–290 panel rows per donor and enters only on the diagonal of the
eigenbasis; two folds were refitted; the Flex re-read has ±0.005 `pds` intervals on 85 embryonic-line targets.

**What it changes.** The "whiten by the signal covariance" construction is closed from the noise side; on a reliable
reference the soft floor is free on `pds` and worth +4.4 % of aggregate cosine, an effect PLOLO's own noise hides.

**286 (the review's "genuinely different estimator"; prereg and `est286.py` hashed 13:54:07Z; 13 KA rows incl. 263 / 257 / 271 / BR249 at ≤ 2e-11).** The signal covariance Ŝ and noise covariance N̂ of each donor line were estimated from CROSS-PRODUCTS of 267's independent batch halves over every bank row except the panel (9.5–17.7 k rows; the moment identity Ĉ = Ŝ + 2N̂ at 5e-17), then used (E1) as a full-matrix Wiener filter on the donor rows with a floor φ ∈ {0, 0.10, 0.25} in every direction, (E2) as a low-rank projection (r = 20, 100), (E3) as the whitener's covariance, with shuffled-pairing, random-direction, mirrored and wrong-source controls and 271 / 272's constructions as baselines. The object is real and different: Ŝ_shuf's leading eigenvalue is 7–19× smaller, and Ŝ's top-100 directions sit at principal cosines 0.52–0.70 to the fitted operator's. **It does not clear anywhere.** Held CD4 donors (280's bank, two donors → the other two): Δcos +0.0004 [−0.0016, +0.0025] at φ = 0.25, Δ`pds` +0.0040 [−0.0021, +0.0101] matched by the random-direction control (+0.0032); held STATE negative (Rest −0.0038 [−0.0058, −0.0019] cosine, −0.0089 `pds`). PLOLO, three cancer folds: best arm −0.0080 `pds` [−0.0163, +0.0001] (1/3 folds), priced −0.0207 pessimistic vs the +0.0089 gate; E2 catastrophic (−0.14 … −0.18); E3 −0.0071 … −0.031; 0/12 cells; every control below the arms. The Flex reference (273's instrument, 85 targets) reads the same sign on `pds` and own cosine (−0.0064 / −0.0013), so the loss is the prediction's transfer, not the truth's noise. **Mechanism, measured:** 78–80 % of the directions have λ⁺(Ŝ) below the independent noise level and carry 81–82 % of the incumbent's whitened prediction energy; Ŝ's top-100 directions carry 1 %; the directions where Ŝ's eigenvalue is NEGATIVE carry 62 %. Only 4–11 % of a half row's covariance is shared response and its top-100 directions hold 19–32 % of that — the shared response is itself high-rank. Every cross-product arm helps the anti-aligned 17.5 % (+0.04 … +0.11) and pays more in the marginal band (−0.02 … −0.22), the shuffled control included: the group signature is shrinkage toward a subspace, not the learnt directions. **So F3's mechanism survives re-estimating the eigenvectors from independent views — the open question 272 left.** Provenance (297): EXTENDED of F3b (which stays NEW); Jolliffe 1982 (low-variance components can matter in prediction) must be cited with it, and "no detectable signal by an unbiased estimator" is not "no signal" at n/p ≈ 1.6 (the estimator class is itself biased — the 2025 cvPCA revision). Not closed by 286: a target-specific use of Ŝ; a per-target gate that harvests the anti-aligned gain without paying in the band (needs a signal 263 / 264 did not find); the one cell where the learnt directions beat their controls without losing (CD4 cosine, +0.0004 vs controls at −0.03). **299 (the fourth review's pilot; prereg hashed 18:08:36Z; 10 KA rows incl. the BR249 retrodiction at 0.0e+00): the whitener fitted on RELIABILITY-SELECTED training rows DOES NOT CLEAR — the preregistered arm `sel:q0.5` loses −0.0675 `pds` [−0.0817, −0.0535] on 0/3 folds at −25 % `cos_agg` and −28 % delivered amplitude, priced −0.089 pessimistic vs +0.0089, and sits BELOW its matched random (−0.061), magnitude-selected (−0.055) and reversed (−0.059) controls; at q < 1 the arm is a row-count effect (14k rows on 17.9k columns → 3,781 floored eigenvalues amplified 7.9 × 10³) and on the one full-rank fold all four arms are at zero with `sel` lowest. The mechanism row says why: the reliability-selected operator puts 99.9 % of its whitened prediction energy where an independent noise estimate finds no signal (incumbent 86.8 %). POST-HOC, at the incumbent's exact row count and donor proportions, a five-point rank-preserving dose curve is monotone with the OPPOSITE sign — tilting the rows TOWARD reliable ones loses (−0.0064 → −0.0011), tilting AWAY gains (+0.00094 [+0.00017, +0.00173] 3/3, shape-only, mirror control passed, +1.9 % `cos_agg`; steeper +0.00139 but amplitude-changing) — priced +0.00140 pessimistic = 0.74 × the shape-only gate with retention 47.8 % (< 50 %): two narrow misses, and the ordering-noise SD (0.0045) is three times the effect; the same sign reproduces on 273's Flex instrument (+0.0014 `pds`, +0.0006 own cosine vs a matched refit). The whitener's `pds` comes from the SPREAD of the donor rows over many directions, not their reproducibility (272, 286, 299 agree). Nothing build-worthy; zero slots on 2026-09-06.** **300 (prereg hashed 20:11:40Z; KA: 299's `wgt:mild` / `wgt:mildrev` per-fold `pds` at 0.0e+00, 062's prune row set reproduced exactly, BR249 at 0.0e+00) then measured the reversed tilt as a MECHANISM at the incumbent's exact row count, with mirror and random-weight controls: the dose curve is monotone and 3/3 positive at every realisable dose (+0.00049 / +0.00094 / +0.00103 / +0.00115 / +0.00150 at d = 0.5 / 1.0 / 1.25 / 1.5 / 2.0; d > 2 gives negative weights) and does NOT saturate; the shape-only band edge sits at d ≈ 1.1–1.2, so the largest shape-only dose is d = 1.0 (+0.00140 pessimistic = 0.74 × the +0.0019 gate, retention 47.8 %) and d = 1.25 prices +0.00180 = 0.95 × the shape-only gate but delivers +0.53 % amplitude and is judged against +0.0089 (0.20 ×) — price and amplitude are driven by the same cosine gain, so NO dose clears, structurally. All three controls pass (mirror −0.00105; random-weight seeds +0.00025 / −0.00023, both losing `cos_agg`): the effect is real, its whole span (+0.0015) a third of the ordering-noise SD. Mechanism: tilting away from reproducible rows moves 272's statistic in the paying direction (whitened energy on noise-dominated directions 86.8 % → 88.1 %; signal share 0.153 → 0.133) in proportion to Δ`pds`; a DIRECT spread manipulation (up-weighting rows with low leverage on the top-50 eigenvectors) is 72 % collinear per target and 70–81 % at the row level yet delivers a quarter of the `pds` (paired −0.00072 [−0.00125, −0.00019]) — the operative part of the reliability axis lives outside the leading directions, where 272 located the `pds`. The prune interaction is unresolved (−0.00023 [−0.00089, +0.00043]); **062's essential-gene prune loses on PLOLO (−0.00169 [−0.00457, +0.00118]) and gains on Flex (+0.0021 [−0.0018, +0.0062]) — its value is NOT established on either instrument (299 §6's side observation corrected)**; the Flex reading of the tilt with an interval, +0.0014 [−0.0004, +0.0034], contains zero. The dose axis of this arm is CLOSED; nothing justifies a slot.** Memory: three processes exceeded the 10 GiB cap (23 GiB, fixed in-session to 4.7 GiB; the machine stayed ≥ 73 % free), reported by the agent.

### F4. Per-target `pds` is a threshold in the per-target transfer cosine at ≈ 0.02, on both references

| own-cosine bin | PLOLO, incumbent, n 853 (R262 §4.2): mean `pds` / share of the deficit | Flex-400, n 85 (R273 §3): mean `pds` |
|---|---|---:|
| (−0.02, −0.01] | 0.15 / 16 % | 0.16 |
| (−0.01, 0] | 0.37 / 34 % | 0.45 |
| (0, 0.01] | 0.66 / 28 % | 0.76 |
| (0.01, 0.02] | 0.86 / 12 % | 0.92 |
| (0.02, 0.03] | 0.97 / 2 % | 0.98 |
| (0.03, 0.05] | 0.99 / 0.5 % | 1.00 |
| > 0.05 | 1.00 / 0.1 % | 1.00 |
| slope inside (0, 0.05], `pds` per unit cosine | 9.5 (290 decoys) | 6.0 (85 decoys) |

| the incumbent's regime and deficit (R262 §3–5, R273 §3) | PLOLO (3 cancer folds) | Flex-400 |
|---|---|---|
| share of targets: anti-aligned (base·truth ≤ 0) / marginal (0, 0.05] / saturated | 17.5 % / 66 % / 16 % | 12.9 % / 48 % / 39 % (55 % at `pds` = 1) |
| share of the deficit by the same bins | 53.1 % / 46.8 % / 0.1 % | — |
| deficit by mechanism, K562 / HCT116 / HEK293T: base points away; residual breaks a positive base; render | 57.5 / 51.5 / 49.0 %; 11.7 / 16.1 / 18.1 %; 8.7 / 10.6 / 12.3 % | — |
| fold consistency: bottom-decile overlap (expected); Spearman of per-target `pds` across held lines | 9 / 12 / 12 (2.2–2.5); +0.76 / +0.70 / +0.74 | — |
| covariates of per-target `pds`: base·truth cosine; held-vs-donor truth agreement; donor-donor agreement; CD4 DE count; hubness; truth norm | ρ +0.875 (AUC 0.05); +0.63; +0.38; +0.38; −0.03; −0.03 | — |

**Limit.** The knee is largely derivable from the decoy-cosine distribution and moves with the decoy count (slopes 9.5 vs 6.0), so its content is the empirical decoy distribution, not the existence of a knee (287; metric-geometry literature: Liu et al. 2025 on PDS scale and distance sensitivity). It is a property of the rank statistic with 85–300 decoys, not of truth noise; what an observed 0.02
means in signal cosine differs by reference (≈ 0.07–0.09 on PLOLO's 0.05-reliable truths, ≈ 0.03 on Flex; R271 §5),
so PLOLO's marginal band and the board's are different target sets. Per-target signs of a re-shaping arm are
line-specific (Spearman of the deflation's per-target Δ across folds −0.10 … +0.03, R262 §6).

**What it changes.** An arm's `pds` is a count of targets it moves across ≈ 0.02; the Flex-400 regime is
prediction-limited (median own cosine 0.037), so the deficit there is the transfer, not the truth.

### F5. On Broad-16 a reproducible line-specific response component exists, is high-rank, and is not predictable from basal EXPRESSION by a linear model at 11 training lines

| component (R268 §2; 82 targets × 12 lines, cell halves, 12,199 genes) | pooled fraction of full-depth response energy | per-target median | OR-control null |
|---|---:|---:|---:|
| shared across lines | 0.091 | 0.031 | 0.019 |
| line-specific and reproducible | 0.115 (direction-changing 0.099, amplitude 0.019) | 0.080 | 0.028 |
| noise | 0.794 | 0.880 | 0.953 |

| structure of the line-specific part (R268 §2–4) | value |
|---|---|
| guide halves (v1 vs v2, 78 targets): shared / line-specific | 0.080 / 0.091 (79 % survives the construct) |
| per (target, line) residual split-half cosine; share of 984 pairs above the permutation p95 (controls) | median 0.046 (p95 0.096); 26 % (7.5 %) |
| within-batch vs between-batch residual similarity | +0.020 vs −0.106, p = 0.001 |
| components for half of the reproducible energy (shared response); above the p99 null of the first 60 | 27 (3); 58 |
| two largest components | 13.5 % and 9.4 %: E2F/G2M proliferation against p53/stress arrest |
| leave-one-line-out, ridge on 2 basal programmes: predicted-residual cosine (wrong recipient; shuffled labels; ceiling) | +0.018 [+0.014, +0.023] (+0.001 [+0.001, +0.002]; −0.007 [−0.009, −0.005]; 0.045 median) |
| the same, gain on the full response: cosine (per-gene / programme-level); within-line retrieval | −0.007 [−0.009, −0.005] / +0.000 [−0.002, +0.002]; +0.0003 [−0.010, +0.011] |
| oracle line offset (held data): cosine / retrieval | +0.169 / −0.019 |

**Limit.** 12 lines of 11 lineages, Cas9 knockouts, 82 targets with zero overlap with the challenge panel, 21–65 cells
per half, one experiment (batch and biology not separable); four lines reserved (JVE127, MG63, SLR23, UMRC3).

**What it changes.** Context dependence is real and as large as the shared response, but it is not a few programmes
and not linearly predictable from two basal expression programmes at this n. **Scope narrowed 2026-09-06 14:30Z (287, CONTRADICTED at the headline):** the source paper (Maffa et al., bioRxiv 2026-08-24) DOES predict line-specific responses from cell state — p53 activation on UBE3A / MDM2 only in the HPV-E6 lines SKGII and C4I, IER3IP1 dependence from an 18q / 11q co-deletion (P = 3.3e-6) — per gene, from genotype and viral status; our null is for a transcriptome-wide linear ridge on two basal expression PCs. Both can hold. The two largest components (E2F/G2M vs p53/stress) are the source paper's own headline programmes, not a discovery here. **289 (Opus; prereg hashed 14:20:21Z; 268's arm reproduced at 0.0e+00) then reran the LOLO with the paper's OWN state features** — p53 functional status (DepMap dependency rule and, separately, mutation annotation), HPV-E6 status, basal E2F/G2M and P53_PATHWAY programme scores, arm-level 18q21 / 11q23-24 copy number (DepMap 24Q2 segments) — on the 12 working lines, 82 targets, with the wrong-recipient and shuffled-label controls and depth as a covariate (rule 9). **Transcriptome-wide, the same negative:** predicted-residual cosine +0.0110 [+0.0060, +0.0162] per gene / +0.0091 programme-level, NOT separated from the wrong-recipient control (+0.0041 [−0.0004, +0.0088] / +0.0048 [−0.0026, +0.0123]) and a LOSS on the full response (−0.0588 [−0.0714, −0.0467] / −0.0236 [−0.0302, −0.0172]); adding the features to 268's two PCs makes it worse; copy number is inert or harmful; 0 of 25 preregistered arm-cells pass; the loss is largest in the DEEPEST lines (Spearman with UMIs per perturbation −0.72, P = 0.008), so it is not a shallow-data artefact. **But the paper's mechanism is in our residual object:** P1 REPRODUCES — after UBE3A / MDM2 knockout the P53_PATHWAY score of the line-specific residual is highest in the two HPV-E6 lines (SKGII +0.10, C4I +0.07; every other line negative; exact permutation over 55 line pairs P = 0.018; the E2F/G2M direction shows the OPPOSITE ordering, P = 1.0), driven by UBE3A (MDM2 alone unresolved at 26–114 cells); the LOLO finds it BLIND (UBE3A is the best-predicted of 82 targets, carried by HPV status on 12/12 folds, +0.29 on the held C4I fold); and per-target SIGN gating on the p53 axis beats chance and its shuffled control (55.9 % vs 50.6 %; paired +0.053 [+0.015, +0.089]; 17.1 % of targets ≥ 9/12 vs 7.8 % shuffled). P2 (IER3IP1 from the co-deletion) does NOT reproduce on 12 lines (3 co-deleted lines, the wrong direction; underpowered by construction). **287's C-1 is RECONCILED, not contradicted:** a state feature predicts which way a particular target bends in a particular line; transcriptome-wide, target-averaged additive prediction from 11 training lines buys nothing; 268's aggregate readout could never have seen the mechanism (UBE3A is 1 of 82 targets, the p53 pathway 176 of 12,199 genes). Note for any state-feature model: low-cardinality features pass the shuffled-label control trivially; the wrong-recipient control is the one that discriminates. HPV-E6 and cervical lineage are the same column on these lines.

### F6. On the non-essential panel K562's DE sets are an order of magnitude sparser than the board's, and matched controls do not change that

| the scorer's Wilcoxon-BH recipe on raw K562-GWPS cells, ≤ 400 per target, 18,400 controls (R260 §2c, R266 §3–4) | value |
|---|---|
| 272 panel targets: median significant genes; at or above `nmae`'s 10-gene gate; with none | 5; 85 (31 %); 58 |
| 50 best-covered candidates: median; ≥ 10; ≥ 89 | 6; 17; 4 |
| at 400 cells: MYRF / NABP2 / PINK1 | 7 / 9 / 6 |
| knockdown detected: target-gene LFC median; negative | −3.05; 229/229 |
| strong K562 knockdowns, same pipeline: RPS3 / HSPA5 / LMO2 | 1,205 / 425 / 342 |
| gemgroup-matched and gemgroup + depth-matched controls: medians; ≥ 10 of 50 | 6.0 / 5.7 / 6.3; 17 / 17 / 18 (260's 50 counts reproduced) |
| what the matching removed: gemgroup composition distance; control UMI CV | 0.39 → 0.005; 0.440 → 0.378 |
| comparators: h1-test reference sets median / `nmae` cohort; the challenge cohort | 89 genes / 58 of 85; 209–261 of 300 |

**Limit.** One line, 10x 3′ at 11.3k UMI per cell; the panel is non-essential in K562 by design. **Provenance (287): KNOWN** — Replogle 2022 reports 30.5 % of 9,608 K562 perturbations with > 10 DEGs; our 31 % is that number on the panel with the scorer's recipe and matched controls. Nadig 2025: DEG counts understate the true effect (a typical perturbation affects ≈ 45 genes), so sparse DE sets ≠ no response.

**What it changes.** No local surface reads the four DE members in the board's regime; the board is the only DE
instrument (F7).

### F7. Board calibrations: what the two 2026-09-05 points established

| quantity (captured entries in `data/board-entries/`) | predicted before the slot | realised | reading |
|---|---|---|---|
| `nmae` identity, six-source union at g 0.35 (BR237 §2) | −0.004783 overall | −0.00716 (raw 0.93389 → 0.96016) | sign right, magnitude a lower bound (×1.50): validated at scale on this class |
| direction legs of the same arm, overall minus `nmae` (BR237 §2) | +0.0022 … +0.0038 (243's class band) | +0.00352, of which `pds` +0.00248 and `mse` +0.00088 | in band; k_cos 0.962, new-source class |
| CD4 deflation k = 5, `pds` (BR249 §1–2) | +0.0032 / +0.0046 / +0.0091 at k_pds 0.7 / 1.0 / 2.0 | raw +0.01492 = 1.22× PLOLO's +0.01226 (+0.00545 overall) | mechanism transfers; k_pds 1.22 |
| same arm, `mse` (BR249 §2) | −0.0003 | −0.00080 (raw +0.00454) | cosine cost under-read 3.8×: k_cos 3.76 on `cos_agg` (2.15 macro); the falsifier fired |
| same arm, overall (BR249 §1) | +0.0026 … +0.0040 | +0.0042634 | above the band's top |
| class factors from three points (R258 §3): k_pds; k_cos; DE-three per unit board Δcos | — | 0.70 / 1.22 / 2.48–3.79; 0.94–1.00 / 1.36 / 3.76; 0.0425 / 0.185 vs the 0.24 recipe | a class property, not a constant |
| retrodictions of BR249 with factors fitted on it (R257, R271, R272; R258; R259) | — | +0.00472; +0.00459 … +0.00494; +0.00470 vs +0.00426 | consistency of the arithmetic, not validation |
| PLOLO's own calibration on six geometry arms (R162) | — | `pds` LOO RMS 0.001494 scaled, 6/6 signs, p 0.0028 (null 0.003481); `mse` RMS 0.112 vs null 0.044; `reach` 2/6 signs | a veto instrument on `pds` only |

R265 corrected the σ language: 0.001494 is a prediction-error scale, not a sampling σ, and R263's "8 σ" cell is 2.3 σ
on its permutation null.

**Limit.** One board point per arm class; k_pds spans 0.7–3.8 over three points; the ordering members carry sd 0.0045
overall between digest-identical arms (AGENTS.md §8.1b), so no ordering claim is readable at these sizes.

**F5 — confirmation on the reserved lines (283, claim C2, one access, reserve spent).** C2a (reproducibility of the
line-specific residual at 268's working-line size, median split-half cosine > 0.028 on ≥ 3 of 4 reserved lines): **FAIL by
the frozen rule** — 1 of 4 clears it (P(≤ 1 of 4 | working rate 8/12) = 0.11); the component is nevertheless detectable
on every reserved line against its own permutation null (+0.0347 [+0.0229, +0.0468], P < 1e-4) and the OR-control null
reads identically on reserve and working lines (0.0148 vs 0.0147). C2b (the 12-line shared response transfers to the
reserved lines within 268's LOLO range): **PASS**. C2c (basal state predicts the reserved residual no better than the
shared response alone): **PASS**. Read F5 as: a detectable, high-rank line-specific component, not predicted by two basal expression programmes in a linear model,
whose size at 268's threshold did not reproduce on held lines. The source paper's own UMI–power law (r = 0.85, elbow ≈ 5 M UMIs per perturbation) pre-explains the miss: the reserved quarter's median 17,428 UMI/cell vs 30,944 on the working lines — C2a's threshold should have been depth-conditioned before it was frozen (287; programme rule 9).

### F8 — A target's cross-line conservation is a reproducible property that predicts its transfer, robust on the split half (276; ledger claim C1, relabelled)

**Claim.** For each of the 300 panel targets, the mean pairwise cosine of its response direction across the cell lines and
datasets on disk (19 lines; frozen definition `E276/work/DEFINITION-276.md`) is a reproducible target property, and,
computed with the held line LEFT OUT, it predicts how well the target's response transfers to that line.

**Measurement.** Unit: targets (295 with a score), pairs of (dataset, line). Reproducibility: across independent cells of
four sources ρ +0.273 [+0.158, +0.388]; across two disjoint sets of lines ρ +0.435 [+0.340, +0.522] (50/50 random line
splits positive, median +0.394). Prediction under leave-the-held-line-out: Spearman with `cd4svd5`'s per-target `pds`
pooled over the three cancer folds ρ +0.387 [+0.326, +0.444], AUC 0.694 for 262's anti-aligned group; on h1-test's 85
Flex-400 targets ρ +0.239 [+0.02, +0.435], AUC 0.634 (R276 §4). **Split-half robustness** (preregistered claim C1, one access, `docs/CONFIRMATION-LEDGER.md`; relabelled 2026-09-06 13:40Z — the half-B cells were inside 276's full-truth evaluation and 267 had scored them, so this is not an untouched confirmation): against the equal-weight mean's per-target `pds` on the reserved HEK293T
half-B batch truths, n = 293, ρ **+0.361 [+0.257, +0.458]**, pass rule "interval above 0" met; the non-reserved half A
reads +0.334 and 276's own deep-truth reading +0.365.

**Magnitude confound (287 T3; measured by the coordinator in E288 `check288.py`).** Within tertiles of an independent response-magnitude proxy (the CD4 screen's DE-gene count per target), ρ(C_lolo, per-target `pds`) is at chance in the weakest tertile (+0.03 [−0.18, +0.24] / +0.18 [−0.02, +0.36] / +0.15 [−0.04, +0.34] on K562 / HCT116 / HEK293T), +0.21…+0.28 in the middle and +0.33…+0.47 in the strongest; the partial Spearman given magnitude is +0.26 / +0.25 / +0.23 (kolf2 +0.11) against the raw +0.41 / +0.39 / +0.37 (+0.28) — roughly 40 % of F8's association is response magnitude and the rest survives. By split-half reliability tertile there is no such gradient (low tertile +0.38 / +0.48 / +0.43), so the confound is magnitude, not reliability. The CD4 source's own table shows the same coupling (ρ(conservation, DE count) +0.51…+0.54). A per-target cross-line conservation score with random controls already exists in that table (K562 vs CD4, median +0.30; 287: EXTENDED) — new here is its reproducibility across disjoint line sets and its prediction with the held line left out.

**Limit.** The 17 extra lines beat the two-donor agreement the pipeline already computes on one of three readings only
(Δρ +0.066 [+0.004, +0.133] PLOLO frozen axis; +0.052 [−0.011, +0.117] wide axis; −0.028 [−0.166, +0.106] Flex-400), so on
this panel more lines at 150 cells add little predictive information over two; the reliability-corrected map (median
+0.228) says the ceiling is far above the raw map, i.e. depth, not breadth, is the binding axis. Cross-line agreement of
L1000 (CMap) signatures does not predict transfer (ρ +0.013 / +0.125) — closed at that construction.

**What it changes.** Conservation is the first target-level quantity on disk that (i) is reproducible, (ii) predicts
transfer without seeing the held line, and (iii) survived a split-half robustness check (C1, relabelled). **Off-map (291):** the map's own per-gene score predicts the GWCD4i authors' published CD4-vs-K562 conservation at ρ +0.274 [+0.210, +0.332] over 1,773 genes (Stim48hr; +0.110 after their own magnitude covariates) and Broad-16's 12-line Cas9 conservation at +0.210 [−0.004, +0.414] over 82 targets — as well as any feature model fitted on it. The map transfers to lines, datasets and a perturbation modality it never saw. **Claim C5 on reserve R3 (295; preregistered 16:09:10Z before the first byte of any reserve file's `X`; one access, all eight files now SPENT): on 174 genes in eight datasets no bank, map row, donor table or fold ever touched — five studies never opened, four of them Cas9 knockout, one a NEW cell type (primary CD8⁺ T cells) — the map's per-gene score predicts each gene's transfer into the reserve at ρ +0.502, file-block 95 % [+0.153, +0.566] (gene bootstrap [+0.371, +0.608]); 0/200 shuffles, z +6.9; partial given the file's reliability and magnitude +0.309 [+0.139, +0.396]; within-file ranks +0.351 [+0.077, +0.485] → C5a PASS.** Baselines: annotation +0.421 (paired +0.080 [−0.075, +0.177] — the map's margin over annotation is NOT separated), magnitude +0.358 (paired +0.148 [+0.038, +0.292], separated). CRISPRi-only +0.517, Cas9-only +0.214, same sign, neither separated on 3 / 5 blocks. Two honest limits: 8 file blocks with one file (Adamson, 78 of 174 genes) carrying the pool, the lower bound just clear of the rule's +0.15; the map's advantage is over magnitude and reliability, not over annotation. This is the programme's first confirmation on data no agent had inspected. Provenance (297): EXTENDED — PerturbMap's prospective reliability weighting is prior art for a transfer-predicting score; a per-gene score validated on never-opened studies under a preregistered rule is not. **Single-file jackknife (298, on 295's saved per-gene table; no reserve file reopened): C5a HOLDS without Adamson — ρ +0.268, gene bootstrap [+0.061, +0.440], file-block [+0.105, +0.429]; Adamson alone +0.518; every other leave-one-out stays at +0.496…+0.540. C5b is MARGINAL without Adamson — +0.195, file-block lower bound +0.006, gene bootstrap [−0.026, +0.380]; dropping Tian instead puts the pooled file-block bound at −0.005 — so C5b rests on the pooled eight-file set and is not robust to leaving out either CRISPRi file that carries it.** It does not yet buy `pds` (it
identifies the anti-aligned rows, which no residual repairs: 262/263), so its competition use, if any, is a rendering
decision on identifiable rows (277), not a dose or a selector.

### F9 — Cell-identity gene families have less conserved perturbation responses than core-machinery families; preregistered and confirmed on 17,106 held genes (279)

**Claim (preregistered, hashed 2026-09-06 06:25:54Z before any held number).** Perturbations of zinc-finger transcription
factors have lower cross-line conservation (276's frozen statistic) than other perturbations; pass rule: raw and
expression/magnitude-matched median differences negative with 95 % intervals below 0.

**Measurement.** Unit: genes; held set = every non-panel gene measured in ≥ 2 lines on disk (17,106; ZNF 696). Raw
Δmedian −0.01047 [−0.01269, −0.00844], AUC 0.564 [0.548, 0.581], one-sided Mann–Whitney p 4.1e-09; matched Δmedian
−0.00635 [−0.01172, −0.00079]. Size-matched random families: 0/100 draws reach the effect (null band [−0.0030,
+0.0045]). Inside one dataset (Orion HCT116 × HEK293T): −0.00901 [−0.01573, −0.00091]; on the balanced three-source
core (identical pairs for every gene): −0.01785 [−0.02197, −0.01231], AUC 0.604. Both arms of the frozen rule fire →
CONFIRMED.

**Scope.** The effect is not specific to zinc fingers: among 93 HGNC families with ≥ 30 held genes, ZNFs rank 13–21;
homeoboxes (AUC 0.603) and bHLH (0.574) show it too, and the family axis runs from core machinery (spliceosome,
proteasome, ribosome biogenesis: +0.07 … +0.11 above the rest) to cell-identity families (receptors, homeoboxes,
zinc fingers). It is state dependence rather than measurement: the deficit GROWS with the reliability tercile (AUC
0.517 → 0.551 → 0.614) and 74 % of it survives the reliability correction; the tercile contrasts are not
expression-matched (stated limit).

**Corroboration on an independent estimator (287, `work/checks287.py` CHECK 2).** Inside the CD4 source's own K562-vs-CD4 table (DESeq2 z-scores, one cross-line pair, 56 ZNFs) the deficit reads Δmedian −0.030 / −0.032 / −0.032, AUC 0.58–0.59, one-sided P 0.013–0.019, the bootstrap interval below 0 in Stim48hr only — corroboration, not confirmation. The direction is what the lineage-dependency and co-essentiality literature expects; the family-resolved, preregistered, magnitude- and reliability-controlled quantification is the new part (287: EXTENDED).

**What it changes.** A biological statement about which perturbation responses are conserved across cell states — the
organising question's first direct answer — and a 17,106-gene cross-line conservation map (`E279/work/znf279.csv`,
28 sources, 21 lines) as a reusable object. Competition: nothing to build (277 closed the class-conditional rendering).

### F10 — In primary CD4 T cells, a 10x Flex perturbation response reproduces across four human donors at +0.054; the Flex-400 regime of F1/F4 is a property of the H1 reference, not of chemistry or cell count (280; preregistered)

**Claim (preregistered, hashed 2026-09-06 06:57:35Z before any expression value was read).** Unit: biological donors (4)
× activation states (3) of the GWCD4i screen (10x GEM-X Flex, MIT licence; the CD4 residual's source), 297 panel
targets, ~250 cells per side of a two-donor split.

**Measurement.** Across-donor residual cosine (Stim48hr, 3 disjoint pairings): median **+0.054 [+0.043, +0.062]**, 28.7 %
of targets ≤ 0; within-donor at matched cells (guide split, 347 cells/side) +0.048 [−0.000, +0.067]. Donor as a unit:
Δ = +0.0033 [+0.0004, +0.0072] (passes the frozen rule; ≈ 1/30 of the cell-sampling term D − A = +0.099 [+0.092, +0.106];
not clearing zero in Rest / Stim8hr; equal to the guide-construct effect). Replicates vs cells: three donors beat one by
+0.0067 [+0.0040, +0.0108] cosine at matched UMI budget (all three states), +0.0129 at full cells — 52 % replicate
averaging, 48 % more cells. Depth curve — **retracted 2026-09-06 13:40Z**: it thinned COUNTS from pseudobulks (UMI depth), not independent cells, and its Spearman–Brown paired the donor-split reliability with the guide cohort's 347 cells (the donor sides carry ~250; `work/depth280.py:98`); no cell-budget ceiling is established. **285's replacement (F14):** Spearman–Brown in cells holds on four real cell ladders (measured − SB within ±0.02), cells and UMIs are one axis, the CD4 donor term is flat in the budget (ρ 0.02–0.07 at 1–4 × 10⁵ UMIs) — more measurement of one construct does not move across-donor agreement, and above 76–330 cells per side it stops moving anything but the construct's own response. Reproducibility
tracks response magnitude (Spearman with DE-gene count +0.38 / +0.49 / +0.59), not cells (−0.03 … +0.07, p ≥ 0.26).
Estimator-independent: the study's own donor-pair DESeq2 rows give +0.09 … +0.11 on their DE-eligible targets
(+0.14 … +0.15 by this statistic on the same targets); selection, not estimator, moves the number. F8's map holds on
this donor surface (ρ +0.155 / +0.295 / +0.295).

**Limit.** n = 4 people; one run-clean condition (Stim48hr); the within-donor comparator is a guide split; a
pseudobulk log1p(CPM) estimator (agrees with DESeq2 rows at +0.35–0.53); one cell type.

**What it changes.** F1's "Flex at 400 cells resolves every target" (measured on H1) does not generalise: same
chemistry, comparable cells, 4.5× lower across-donor agreement in a primary somatic cell type on this non-essential panel — and across-donor agreement is a different quantity from the challenge's within-reference split-cell repeatability. Any
statement of the form "the challenge's reference resolves directions" must name the cell type and the response
magnitude; the somatic A/B/C references are more plausibly in the +0.05 regime than in H1's. 273's re-reads (F4's
regime shares, the "prediction-limited" reading) therefore describe H1-like references only. *(A closing sentence — "the depth axis is closed on this disk: cells are not the binding constraint at any genome-scale budget; a resolving reference needs a different cell type" — was DELETED 2026-09-06 14:30Z: it had survived the retraction above and is contradicted by two published depth laws, Wang et al. 2026's 84 % → 10 % unreliable from < 50 to > 500 cells and Broad-16's UMI–power relation r = 0.85; 285 remeasures the cell axis.)* **Provenance (287): F10a is KNOWN** — the source's own `donor_correlation_all_mean` medians read +0.107 / +0.112 / +0.114 across all genes (+0.10 / +0.09 / +0.13 on the panel targets) and +0.70–0.76 on hit genes; F10a re-measures that published pair with a cosine estimator on the panel. **NEW is F10b**: the donor as a marginal replicate unit and the replicates-vs-cells exchange at matched budget — in tension with the 34-line iPSC map (Cell Genomics 2025), which finds heritable donor modulation of trans effects at 26 donors; n = 4 here.

### F11 — Within one individual the across-state loss of a CD4 response's agreement is not distinguishable from zero at ~107 cells per side (the "half" ratio downgraded), cross-line conservation does not predict it, and at small cell numbers the reproducible part of a pseudobulk is a donor × guide construct effect, not a gene effect (282; preregistered, the preregistered index retracted as an artefact)

**Claim (preregistered, hashed 2026-09-06 08:03:50Z).** Across-state conservation of a target's response (Rest / Stim8hr /
Stim48hr) in the same donors' CD4 T cells, against the across-donor ceiling within a state; frozen prediction: low
cross-line conservation ⇒ high within-individual state dependence (ρ < 0).

**Measurement.** Unit: (donor pair, state, target), 280's residual geometry, 5/5 KA rows (S2 reproduces 280's ceiling).
The preregistered index read −2.68 [−3.56, −2.07] — an artefact: same donor + same guide across states agree 5.6× better
than across donors, and the crossed control assigns that agreement to the (donor × guide) construct. Deconfounded (guide
broken across the comparison): **state-dependence index SD = 1 − S1/S2, median +0.49 [+0.33, +0.61]** — **downgraded 2026-09-06 13:40Z: the ratio conditions on the 49 targets with positive denominators (not 74; `work/probe3_282.py:70`), and the absolute across-state difference is +0.0080 [−0.0026, +0.0186], which includes zero; "state removes half" is not established** (~107 cells/side; depth-matched ratio +0.43 [+0.28, +0.53]; pooled ratio 0.385). Prediction: preregistered ρ(SD, C) = +0.262 [+0.120,
+0.399] — the wrong sign; deconfounded ρ −0.144 [−0.403, +0.126] — covers 0 → **NOT confirmed**. Conservation predicts the
across-donor ceiling (ρ +0.301 [+0.177, +0.419]) and not across-state agreement (ρ +0.084 [−0.045, +0.209]). Crossed
cells at ~66 cells/side, matched depth: same donor + same guide, state changed +0.0424; same donor, other guide for the
same gene, same state +0.0116; same guide, different donor, same state +0.0141; any two changed +0.008 … +0.011.

**Limit.** 74 targets with the full crossed design; ~66–107 cells per side; one cell type; pseudobulk log1p(CPM).

**What it changes.** Two things. (i) The SIZE of the within-individual state effect could not be measured at ~107 cells per side (absolute difference +0.0080 [−0.0026, +0.0186]); state dependence itself is the CD4 source paper's headline (287: KNOWN). What is new and preregistered is F11c: cross-line conservation does not predict across-state agreement (deconfounded ρ −0.144 [−0.403, +0.126] — "no relation detected", not "no relation"), while it does track the across-donor ceiling (+0.301; the authors' own columns corroborate at +0.11…+0.17). "Conserved across cell lines" and "stable across cell states within a person" are empirically different axes. (ii) At the cell numbers of a genome-scale
screen, the reproducible component of a pseudobulk is dominated by the construct (guide efficiency / off-target), not
by the gene: any reference built from one construct per target — the challenge's — carries construct-specific signal
that a second construct would not reproduce (F1's guide halves, +0.014, said the same from the other side). The
preregistered statistic is retracted; the deconfounded one is reported.

### F12 — Cross-line conservation is predictable from build-time gene features at the family level (ρ +0.81 across 51 families) and only partly at the gene level (ρ +0.32, no better than measured reliability) (281; preregistered)

**Claim (preregistered, hashed 2026-09-06 08:15:18Z).** From features known before any perturbation (DepMap essentiality
and gene effect, CCLE expression breadth and level, regulon size, TF class, complex membership, paralogue count, HGNC
family), predict 279's 17,106-gene conservation with whole HGNC families held out.

**Measurement.** Unit: genes, block-bootstrapped over held families (3,422 genes in 481 families). Ridge on 21 features:
held-family Spearman **+0.320 [+0.256, +0.379]**, R² **0.111 [0.074, 0.146]** — the frozen rule's YES (**caveat 2026-09-06 13:40Z: the blocks use each gene's FIRST HGNC group, so multi-family genes can share another family across the split; 284 reruns with all memberships excluded**). Against the
confounder: reliability alone predicts as well (features − reliability: −0.031 [−0.071, +0.012]); features add +0.056
[+0.040, +0.070] on top of reliability, shrinking to +0.012 [−0.003, +0.026] on the measurement-matched balanced core
(reliability ρ +0.520 vs features +0.382). Family level: out-of-fold predicted family medians reproduce 279's axis at
Spearman **+0.811 [+0.641, +0.916]** over 51 families; the zinc-finger deficit is 84 % predicted from biology alone
(observed −0.01785, predicted −0.01506) — **57–80 % under 284's stricter folds (below)**.

**284 (rerun with EVERY membership of a held family excluded from training; prereg hashed 2026-09-06 13:46:19Z; 281's primary cell reproduced at 0.0e+00 first).** The leak was real — 43.7 % of 281's test genes shared an HGNC family with a training gene, 784 training rows carried a held family as a non-first group — and it did NOT inflate the result: on the all-membership split (0 training rows carry a held family) ridge ρ **+0.334 [+0.216, +0.355]**, R² **0.112 [0.047, 0.125]** (component block bootstrap, degenerate because one gene–family component holds 63 % of the genes; 281's first-family bootstrap on the same fit [+0.271, +0.392]); features − reliability −0.035 [−0.082, +0.033]; on the measurement-matched core +0.413 vs reliability +0.529, features adding +0.015 [−0.009, +0.029]. Family axis under folds in which a test gene shares NO family with its training rows: **+0.805 [+0.679, +0.883]** over 78 families; the ZNF deficit **57 % predicted** (−0.01009 of −0.01785; 80 % on first-family folds; the strict fold also trains on 6,289 non-giant genes only, so the drop mixes the leak with a smaller training set). **The remainder test (new):** the 28 sources split into two halves sharing no source and no cell line; each half's gene-level remainder (its map minus the features' out-of-fold prediction) correlates across halves at **ρ +0.328 [+0.241, +0.373]** against a ceiling (raw maps) of +0.382 [+0.267, +0.410], a within-component permutation null of +0.044 ± 0.018 (0/200 draws), 20/20 admissible partitions positive (median +0.319); +0.277 [+0.175, +0.304] after measured reliability and response magnitude are also subtracted; +0.322 on the 1,816 genes OFF 279's axis, whose own knockdown cannot contribute. Frozen rule met → **the remainder is a gene property the features do not predict**, not measurement noise (297: NEW, the load-bearing new result of the round; the described tail is a description, not a tested contrast). Described tail: mRNA-surveillance / exon-junction machinery conserved MORE than essentiality predicts (SMG5, SMG7, UPF2, CASC3, ZC3H13, EIF3H), proteostasis and chromosome-segregation machinery LESS (CDC37, VCP, PSMB4, NCAPG, NCAPD2, HAUS6) — the remainder separates WITHIN the essential genes, where the essentiality × expression axis has no resolution.

**Limit.** The family-axis number is released at 284's +0.805 [+0.679, +0.883] (strict folds); the direction — functionally related genes behave alike across lines — is what co-essentiality maps predict, so it is a calibration of a known structure (287: EXTENDED). F12b's "the gene residual is mostly measurement" is WITHDRAWN: 284 shows the remainder reproduces across disjoint source sets (ρ +0.328 vs ceiling +0.382). Both halves share the same 28 banks' pipelines and the 2,739-gene axis (a shared processing artefact would appear in both); "gene property" is a statement about reproducibility, not mechanism, and gene-level technical covariates correlated with gene identity across studies are not all excluded. Gene-level R² 0.11; features and reliability are partly
collinear (expressed, essential genes are both better measured and more conserved); one map, one statistic.

**What it changes.** The core-machinery-to-cell-identity axis (F9) is a property of gene biology that can be stated
before any perturbation is measured; a gene's own conservation beyond its family is NOT predictable from 21 annotation features — but it is real and reproducible (284), worth about as much rank information as the whole feature model (+0.33 reproducible vs +0.33 predicted). The map is therefore the necessary object: for a model, the family axis is a prior annotation can carry, and the per-gene remainder is information only a measured map, or response-derived features, can carry. Regulon size and TF status still read at chance.


### F13 — A response's cross-line conservation is predictable from its own GEOMETRY measured in other lines: broad, spectrally central responses transfer, narrow idiosyncratic ones do not; programme identity adds little (290; preregistered, two of four predictions refuted with the sign inverted). The SIGN is confirmed on two surfaces outside the map; "geometry beats annotation" holds only on 290's remainder target (291, claims C3 / C4)

**Claim (preregistered, hashed 2026-09-06 14:23:29Z before any feature existed).** Response-derived features computed on
one half of the 28 sources (284's primary partition: 14 sources / 10 lines vs 14 sources / 11 lines, no shared source or
line) predict the other half's gene-level conservation REMAINDER (284's object: the map minus the 21 annotation features'
out-of-fold prediction); and — the frozen mechanism — responses engaging the recipient's proliferation / stress axis, and
generic responses, are the LESS conserved ones.

**Measurement.** Unit: 2,848 genes with ≥ 1 cross-line pair in both halves; grouped folds over family components; component
and first-family block bootstraps. Features per gene, aggregated over the half's sources: magnitude, breadth (a bank-level
DE-count proxy), genericness (cosine to the bank's mean row; energy in the bank's top-5 / top-20 PCs), signed projections
on nine hallmark sets and a proliferation summary, on-target knockdown where the gene is on the axis. **P1 PASS:** ridge
out-of-fold Spearman **+0.317 [+0.290, +0.391]** against the annotation model's **−0.057 [−0.138, −0.021]** on the same
target — and at the empirical CEILING, half 1's own measured remainder (+0.328; paired Δ +0.005 [−0.033, +0.101]).
**P2 PASS:** beats annotation by +0.395 [+0.344, +0.490]. Symmetric direction +0.341; declared replication partitions
+0.244 / +0.255 and +0.247 / +0.237; 40/40 partition × direction cells positive (median +0.248). Not magnitude (marginal
+0.045 [−0.060, +0.066]; partialled out +0.316), not reliability (+0.278 after `r_bar`), not coverage, not one dataset
(without the two Orion banks +0.308; single banks +0.15 … +0.41), not a leak (200 shuffles −0.007 ± 0.022, 0/200; the
same-half leakage reading would have been +0.488). **P3 FAIL, sign inverted:** proliferation engagement → MORE conserved
(ρ **+0.118 [+0.100, +0.233]**). **P4 FAIL, sign inverted:** genericness → MORE conserved (top-5 PC energy **+0.115
[+0.069, +0.188]**; top-20 +0.158). Attribution: genericness alone +0.235, breadth +0.183 permutation drop; the nine
hallmark sets replaced by random matched gene sets cost only +0.023 (model +0.294 ± 0.007), and a random 228-gene
direction already reads +0.079 of the proliferation summary's +0.118 — **the signal is spectral energy geometry, not
programme identity** — narrowed by 298: fitted ALONE the ten hallmark projections reach out-of-fold ρ +0.158 [+0.139, +0.300] against +0.099 ± 0.019 for random matched gene sets (0/20 draws reach it), so pathway identity carries real signal on its own and is REDUNDANT given geometry (spectral columns alone +0.242; full +0.317): "identity adds little once geometry is in the model", not "identity does not predict"; 297's contradiction with Wang et al. 2026's Hallmark contrast dissolves. **302 (F15, a fourth surface and a second modality): `gen_pc20` of the knockdown row predicts the MAGNITUDE of overlap between a TF's overexpression and knockdown responses in hESC (ρ(|S|, gen_pc20) +0.38 [+0.23, +0.52], partial +0.38 given reliability, magnitudes and cells) but not the sign (+0.06, interval across zero) — size, not sign.** On-target knockdown depth also goes with MORE conservation (ρ −0.196 on the 1,032 on-axis genes).
Systematic exception (the 20 most under-predicted genes): sequence-specific DNA binders — four C2H2 zinc fingers, TCF3,
CAMTA1, CREM, MTF2, TET1, TGIF2 — are narrow, non-generic and MORE conserved than geometry predicts. Correction carried
into F12: 281's `r_bar` is measured on four banks (K562-GWPS, HCT116, HEK293T, DLD-1), all in half 1, not on all 28
sources; no 284 verdict changes.

**291 (two surfaces OUTSIDE the 28-bank map; claims C3 / C4 written into `docs/CONFIRMATION-LEDGER.md` before any comparison; prereg hashed 15:13:32Z; 290's cell, 279's map, 287's random-control row and 268's decomposition reproduced first).** The 290 ridge fitted on the map's own conservation, then APPLIED: (C3) Broad-16, 82 Cas9 targets across 12 cancer lines disjoint from every bank — ρ +0.210 [−0.013, +0.426] (interval touches zero; 10/200 target shuffles reach it) → **FAIL on that condition, low-powered** (±0.22 on 82 targets); the single feature `gen_pc20` +0.295 [+0.070, +0.507] PASSES; annotation +0.300 [+0.083, +0.489]; the map's own score +0.210. (C4) the GWCD4i authors' published CD4-vs-K562 conservation, 1,773 genes, features from the 23 K562-free banks — ρ **+0.203 [+0.130, +0.252]** (0/200 shuffles, z +8.4), partial given the authors' own magnitude counts +0.061 [+0.030, +0.156], `gen_pc20` +0.186 [+0.119, +0.232] and `breadth` +0.128 [+0.016, +0.156] all PASS; **annotation +0.270 [+0.192, +0.311] is NOT beaten** (paired −0.064 [−0.124, +0.011]) → **FAIL on that condition**, identically in all three states. **Reading:** the inverted sign — spectrally central, broad responses are MORE conserved — is CONFIRMED off-map on both surfaces as a single feature and survives the authors' magnitude; the headline "geometry beats annotation by +0.39" was a property of 290's TARGET (the remainder after the annotation model, orthogonal to annotation by construction), not of conservation; off-map the published statistic is mostly magnitude (the condition's DE count +0.54) and geometry keeps about a quarter of its rank correlation after it; the DNA-binder exception does NOT reproduce on the CD4 surface (annotated TFs sit BELOW the model); and the map's own score transfers as well as any model fitted on it (+0.210 vs +0.210; +0.274 vs +0.203). Neither surface was untouched (268 / 283 / 289; 287) and neither is a reserve; both have now been compared with these features once.

**Limit.** One map, one statistic (279's residual cosine on the 2,739-gene common axis); both halves share the banks'
pipelines; genericness is defined against each bank's own row population; breadth is a relative quantile count, not a
DE test; 74 % of the unit is one family component (a family-blind split gains only 0.03); R² 0.069 on a rank statistic
whose own ceiling is +0.38.

**What it changes.** 279's axis read from the response side: core-machinery perturbations produce large, broad,
spectrally central responses that transfer across lines; cell-identity perturbations produce small, narrow, idiosyncratic
ones that do not. The DNA-binder tail (conserved beyond its geometry inside the map) did not reproduce off-map and is dropped from the claim. "A spectral feature set rather than a pathway-level one" is WITHDRAWN outside the map: on new surfaces annotation and geometry are two partly independent views (ρ(M_resp, y | M_ann) +0.16 / +0.11) and annotation is as good or better; the response geometry and the measured remainder are likewise two partly independent views inside the map (each adds ≈ +0.22 beyond the other; together +0.383). Competition: the
donor rows ARE on disk at build time, so the geometry is computable for every panel target — but the actions a
per-target conservation predictor could key were closed (277's class-conditional rendering; 262/263's anti-aligned rows
that no residual repairs), so nothing is proposed. **Provenance (297): the SIGN is KNOWN — Wang et al. 2026 report model performance ordered Shared > Specific > Unreliable in 39 of 41 method–benchmark combinations on 7,170 perturbations, and the CD4 source's own conservation column tracks its DE count at ρ +0.51…+0.54; 290's source-disjoint, annotation-orthogonal construction at its ceiling is EXTENDED; the identity-vs-geometry ablation is NEW (as an ablation).** Scientific: what is established after 291 is the SIGN (a broad, spectrally central response transfers; a narrow idiosyncratic one does not), measured inside the map and confirmed on two outside surfaces — **and on a third, untouched one: C5b PASS (295), `gen_pc20` ρ +0.371, file-block [+0.064, +0.415] on the 174 reserve genes (`breadth` +0.437 [+0.123, +0.521]); it does not beat annotation there either (paired −0.041 [−0.224, +0.094])** — and that the map's own per-gene score transfers off-map (F8).


### F14 — Cells and UMI depth are ONE removable axis (Spearman–Brown in total UMIs, measured on real cell ladders); the guide construct is the dominant irreducible term of a perturbation response — larger than a different study in the same line — the batch is free, a different line costs about half (285; preregistered)

**Claim (preregistered, hashed 2026-09-06 13:54:51Z with only file metadata read).** On one matched design — the same
targets, the same cells per side, the residual-cosine AND the scorer's own per-target statistic — measure which
disagreement between two measurements of one perturbation more measurement can remove, and which it cannot.

**Measurement.** Units: same-guide random cell halves at n ∈ {25, 50, 100, 200}; thinned counts at fixed n; two
constructs of the same gene (K562-GWPS **748** two-construct perturbations, Nadig Jurkat 124, HepG2 70, h1-test 5,
Broad-16 Cas9 97 targets × 12 lines); gemgroups / lanes / batches; a different STUDY on the same line (K562-GWPS vs
Replogle's K562 essential screen, 1,765 shared perturbations); different lines (Broad-16, 66 pairs); donors and states
(the CD4 bank at matched UMI budgets). KA: 267's medians at 0.0e+00, 273's +0.242920 recomputed at 0.0e+00, 159's
0.330039, 280's +0.0539 within 1.5e-3, the halves' union identity at 0. (i) **Cells and depth are one axis:** at equal
total UMIs, 200 cells at half depth vs 100 at full depth differ by +0.0078 [−0.0068, +0.0119] (h1-test), +0.0012 /
−0.0017 (K562); thinning ladders lie on the cell ladders. (ii) **Spearman–Brown in cells IS the depth law:** measured −
SB from n = 25 reads −0.006 [−0.030, +0.033] (h1-test, 25 → 200: predicted +0.2455, measured +0.2392), −0.006 (K562),
−0.004 (Jurkat), −0.018 (HepG2) — the replacement for 280's retracted curve. (iii) **The construct is the dominant
irreducible term:** at n = 25 two constructs of the same gene agree at **+0.0006 [−0.0002, +0.0014]** (K562, 748; 48 %
≤ 0) against **+0.0063 [+0.0055, +0.0071]** for two halves of ONE construct's cells (CELLS − GUIDE +0.0058 [+0.0045,
+0.0069]); in the scorer's statistic mean `pds` **0.633 → 0.509** (chance 0.500). Replicated: Jurkat +0.0345 [+0.0269,
+0.0452] (`pds` 0.800 → 0.560), HepG2 +0.0963 [+0.0674, +0.1447], Broad-16 +0.0022 [+0.0009, +0.0032] pooled (9/12
lines). Each unit's own ceiling ρ (its agreement over the same-construct cell split at the same n — what more measurement
cannot change): **guide construct** 0.09 [−0.04, +0.23] K562 · 0.17 [0.09, 0.29] Jurkat · 0.41 [0.20, 0.55] HepG2 ·
0.02–0.06 CD4 · 0.76–0.81 Broad-16 (Cas9 knockout: a modality hypothesis, not a controlled contrast); **batch / lane /
gemgroup** 0.87–1.12 on five datasets (free), one exception Orion HEK293T 0.42–0.47 and FLAT across a 5× budget ladder;
**a different study, same line** 0.79 [0.75, 0.82] (cheaper than a second construct; `pds` 0.843 → 0.825); **a different
line** 0.35–0.45; **a different activation state**, same donor and construct, 0.49–0.53; **a different donor** 0.02–0.07,
flat in the budget. Effective guides raise the construct ceiling (ρ 0.12 → 0.24 across on-target-knockdown terciles)
but do not abolish it, and knockdown strength predicts the cell-split agreement as strongly as the guide agreement
(magnitude, not off-target cleanliness). More cells of ONE construct stop helping at **76–330 cells per side** (HepG2 /
Jurkat / h1-test). Reproducibility tracks response magnitude at fixed n (h1-test terciles +0.044 / +0.140 / +0.506; the
K562 PANEL targets +0.0095 vs K562 essential targets +0.0605 in the same file at the same n). The bounded CD4
cell-level extraction was NOT admissible: the object is uncompressed CSR but its chunks are not laid out linearly
(≈ 27 h of range requests for the 62,842 planned rows); ≈ 0.78 GB of metadata transferred, no expression value read,
the plan stored (R285 §7).

**Provenance (297).** The construct as the dominant term is KNOWN (the CD4 source's own `guide_correlation_all` column, +0.06–0.07, is our CD4 ρ_guide; Replogle 2022's dual-sgRNA design and Mixscale exist because of it); the Spearman–Brown reliability law is KNOWN (Wang et al. 2026, ρ(N) = Nτ²/(Nτ² + 1) — ours verifies it on real cell ladders); the line / state / donor ladder is KNOWN within ~1.5× of published counterparts; cells-and-depth-as-one-axis is EXTENDED (for a pseudobulk statistic — the sequencing-budget literature's opposite headline is for per-cell objectives; T6); batch free is EXTENDED; NEW: a different study in the same line agrees better than a second construct in the same study (a reliability ratio, not Replogle's raw r 0.50; T8 — **settled by 298 in raw units: on the same 112 K562 perturbations at full depth, a different study agrees at raw cosine +0.168 [+0.120, +0.215] and a second construct at +0.011 [+0.004, +0.018] (difference +0.157 [+0.108, +0.206]); the pooled study raw cosine +0.215 reaches +0.411 on the strongest-responding third, which reconciles it with Replogle's ≈ 0.50 on a selected cohort, while the construct term is FLAT with magnitude, +0.006 in every tertile**). **Limit.** The guide unit is measured on the targets that carry two constructs (1–10 % of a library, not random, mostly
non-panel in K562); ρ is a ratio of medians with a small denominator on weak-response sources (the CELLS − GUIDE
differences are the tight statistic); Broad-16's cell split mixes constructs; CD4's cell reference is a count split
corrected by a transferred offset (+0.044); nothing is measured on A/B/C. Two jobs peaked at 11.8–12.0 GiB, serially,
with the machine ≥ 73 % free.

**What it changes.** (i) F10's depth question is answered: the challenge's 400 cells of one construct are a
well-measured quantity, and what they measure is that CONSTRUCT's response, not the gene's — above ~100–300 cells the
truth's reproducible remainder is construct-specific. (ii) Across-donor and across-state agreement are not depth
problems. (iii) For transfer, stated as a HYPOTHESIS with the numbers as its motivation: a predictor built from other
constructs (public screens with their own libraries) faces a construct ceiling at both ends — with ρ_guide ≈ 0.1–0.4
(CRISPRi) and ρ_line ≈ 0.35–0.45 a cross-construct, cross-line prediction can share only a small fraction of a truth's
reproducible variance, of the order of the incumbent's median own cosine on the Flex reference (0.037, F4). **Fact checked 16:15Z:** the 2026 panel (`data/challenge/pert_counts.csv`, 300 targets; `manifest.json`: 300 constructs, 400 cells per perturbation, six CRISPRi lines on 10x Flex, zero-shot — no perturbed cells from any 2026 context are released) shares only **25 of 300 targets** with the 2025 H1 library (train + validation + test, 300 targets), so the 2025 H1 source cannot carry the 2026 constructs' responses for 275 targets, and no source on disk carries them for any. **The construct ceiling therefore binds every zero-shot predictor in this challenge equally:** a truth's construct-specific component (1 − ρ_guide of its reproducible variance above the cell-sampling floor; ρ_guide 0.09–0.41 on public CRISPRi libraries, 0.06 on the 2025 library's 5 two-construct test targets at 400 cells, R273 §2) is unpredictable without that library's data. **293 measured the ceiling in the challenge's own regime (prereg hashed 16:02:51Z with only obs read; KA: 273's +0.242920 and its 5-target construct row at 0.0e+00; the 2025 files thinned to h1-test's depth reproduce 273's lane reading to 1.2 %):** on the 2025 library's 17 two-construct targets (H1, Flex; h1-test and the 2026 controls sit at ≈ 20 k UMI/cell, the 2025 files at 54 k — "400 cells" needs its depth stated), ρ_guide is **+0.075 [−0.010, +0.288]** at 200/200 — at K562-GWPS's level, not above it (the frozen prediction that a validated dual-guide library would read 0.2–0.5 is REFUTED) — flat across an 8× cell and 2.7× depth range; lanes free ((a) − (d) +0.011 [−0.003, +0.033]); at the median two constructs of one gene agree no better than constructs of DIFFERENT genes ((b) − (x) +0.007 [−0.033, +0.083]) while the scorer's rank statistic still separates them (`pds` 1.000 same construct → 0.677 other construct → 0.539 other gene). **The mechanism is knockdown failure, not construct-specific biology:** where both constructs reach |lfc| ≥ 1.0, ρ_guide is **+0.687 [+0.197, +1.067]** (5 targets); where the weaker does not, **−0.003 [−0.085, +0.181]** (11); the weaker construct's knockdown correlates with the per-target ratio at Spearman +0.63 (P = 0.009) — and two-construct targets are the library's WEAKEST (34 constructs vs 283 single-construct targets: 38 % vs 17 % below |lfc| 1.0, Mann–Whitney P = 0.0015), so 0.075 is a lower bound and +0.687 is the estimate for the ~83 % of targets whose one construct works. **A construct that misses its target still produces a REPRODUCIBLE transcriptome** (HIRA_P1 at lfc +0.03: split-half +0.268; MED15_P2 at −0.05: +0.047; their knocking-down partners +0.882 / +0.815 and pointing elsewhere): a property of the benchmark that caps every entry equally — **but the extrapolation to the 2025 panel's 18.7 % failure rate ("roughly one 2026 truth in five") is WITHDRAWN (297, T7): it rests on two constructs, against a characterised off-target mechanism reported at 0.8–4.7 % of guides (Hartman 2026); "at least some failing constructs" is what is established, and the 34-construct split-half-vs-matched-controls test with a seed scan would settle the rate.** 293's selection effect (two-construct targets are the library's weakest, so any library-wide two-guide statistic — the CD4 source's +0.06 included — is a lower bound) and its construct-level rather than cell-split ceiling are EXTENDED; knockdown failure as the mechanism is KNOWN (the CD4 source states it; Mixscale is built on it). **Ceilings (variance model g + u + e, SB to 400/400, challenge depth):** a 400-cell truth against another of the SAME construct +0.476 [+0.260, +0.717]; of a DIFFERENT construct +0.046 [−0.006, +0.121]; a perfect construct-level prediction cos +0.690 [+0.510, +0.847]; a perfect gene-level (zero-shot) prediction cos **+0.213 [+0.000, +0.348]** as measured, **≈ +0.57** knockdown-matched (bracket 0.31–0.69). **CORRECTION to the sentence above (16:45Z):** the construct ceiling caps every zero-shot entry but does NOT currently bind — even the pessimistic 0.21 is ≈ 5× the incumbent's median per-target transfer cosine on this reference (0.04, F4) — so the marginal band is the prediction's transfer, not the construct term. **296 (prereg hashed 16:41:42Z; KA: 285's K562 guide row and §6b terciles, 293's strata at 0.0e+00; the CD4 bank's on-target |lfc| vs the authors' own t-statistic ρ +0.64) — GUIDE EFFICACY EXPLAINS THE CONSTRUCT TERM: PARTLY.** Stratified by the weaker construct's knockdown at a library-relative cut: Nadig Jurkat ρ_guide **+0.458 [+0.245, +1.180]** where both constructs knock down vs **+0.096 [−0.175, +0.210]** where the weaker fails (passes at all three cuts and both rungs; knockdown predicts the construct agreement at +0.446 and the cell split at −0.022); K562-GWPS **+0.340 [+0.113, +0.548]** vs **+0.001 [−0.082, +0.105]** at n = 50 (but knockdown predicts the cell split as strongly as the guide split there — 285 §6b's magnitude caution reproduced); HepG2 +0.645 vs +0.307 (not separated); Broad-16 (Cas9) +0.951 vs +0.763 — no split, its constructs work regardless of the mRNA readout; **primary CD4 T cells +0.018 vs +0.014 — flat at every condition, budget and cut with 41–199 pairs per stratum, a powered null.** So "guide efficacy is the dominant irreducible term" is NOT a general claim: in immortalised-line CRISPRi libraries a validated construct pair agrees much better than an unvalidated one; in primary T cells it does not; Arc's 2025 library (+0.687 vs −0.003) is the favourable case. **Cross-gene null correction (new):** in all four CRISPRi libraries the between-construct agreement is NOT separated from a cross-GENE null at the median ((b) − (x) intervals through zero; only CD4's clears) — HepG2's 0.41 comes from a growth-arrest programme its essential targets share, and null-corrected the ordering is K562 ≈ HepG2 ≈ CD4 ≈ 0.03–0.06, Jurkat 0.08, Broad-16 0.79 (the one term that survives its own null; its `pds` 0.70 vs a null 0.50). F14's cross-library ρ_guide ordering therefore partly ranks shared biology, not construct-specific biology; K562's ALL row (a ratio of ≈ 0 to a small number) is not stably estimated and only its stratified contrast is quoted. **Failing constructs (|lfc| < 0.5): reproducible in all five libraries** (own split-half above zero: K562 +0.0069, Jurkat +0.035, HepG2 +0.138, Broad-16 +0.091, CD4 +0.104 — 293's two anecdotes generalised to 2,400+ constructs) but systematically LESS reproducible than working ones, no closer to the library mean, and failing constructs of different genes agree no better than working ones un-centred — **no shared vector signature**; what a failing construct's reproducible response is (an off-target of that spacer, a sub-threshold knockdown, a transduced-population effect) is not distinguished. **Failure rates:** the absolute |lfc| cut is not comparable across libraries (Arc's panels were selected for high basal expression; log1p compresses low-expressed genes: K562-GWPS reads 70–75 % below |lfc| 1.0 by that cut); the transferable numbers are Arc's validated panel 18.7 % and the CD4 authors' own significance call, 26.7 % of 73,765 guide × condition rows not a significant knockdown — so **roughly one target in five to one in four carries a truth whose construct did not knock the gene down**, and that construct's response is reproducible but is not the vector; in every CRISPRi library on disk a "construct" is a transcription-start-site (P1 vs P2) unit. **T7 measured (301): the 2025 library's protospacers are Replogle 2022's Table S1 (Arc's "Behind the Data" cites its protospacers and cloning; 328/348 guide ids resolve, all 34 two-construct arms; non-targeting ids match verbatim); with CAGE-TSS promoter windows (Horlbeck 2016's annotation, GRCh37; KA: 68/68 arm spacers recover their own site, after a GRCh38 gene-TSS window set failed that row at 22/68) a PAM-anchored seed scan finds a distal 12-nt seed + NGG inside a candidate gene's CRISPRi window for 4/13 failing arms vs 4/21 working arms (Fisher P = 0.68) and 4/31 non-targeting spacers — NO enrichment at any seed length, both window definitions; the one apparent hit (ZNF581_P2 → ZNF580) is a shared-locus artefact. Failing constructs' reproducible responses are not promoter seed off-targets; the two-construct term is knockdown failure. The per-construct expression test of the eight candidates needs 293's residual rows rebuilt (not saved).** **302 (F15): the Spearman–Brown depth law holds for ORF overexpression in hESCs (c ≈ 9,600 cells) and two isoform constructs of one TF agree at +0.021 — F14's law and construct term in a second modality.** Limits: n = 17 / 16 / 5; one line, one chemistry, one library; the (x) null and (a2) column were added during the run (labelled); the variance model assumes exchangeable constructs.


### F15 — No detected median reversal in the observed processed responses (302; C6 VOID; reviewed by 304)

**Claim (preregistered, hashed 2026-09-07 03:33:19Z and written into the ledger before the first byte of the reserve's
`X`).** On the 148 transcription factors both overexpressed in the JoungZhang 2023 hESC TF-ORF atlas (≥ 100 cells) and
knocked down in the map's hESC CRISPRi screen (`gse295214_h1`), the cosine S between the two residual responses on a
15,710-gene axis is negative at the median (C6a), cross-line conservation predicts −S (C6b), and the knockdown row's
`gen_pc20` predicts −S (C6c). One streaming pass over the atlas (77 s, 7.0 GiB); R4 is SPENT.

**Measurement.** KA rows: the control-code cell counts, the 249 / 233 / 148 overlap, the H1 bank reproducing 279's
`Ccross` on 26 pairs at 0.0e+00, and the positive control — overexpressing TF X raises X's own transcript in 95.1 % of
TFs and knocking it down lowers it in 100 % — all PASS; **KA-E FAILED**: the file's own `obs.ncounts` disagrees with
its `X` row sums by a batch-keyed 6–18 % (and its two count columns disagree with each other by 3.2 %; `X` holds
half-integers), a property of the released file that the frozen clause turns into VOID verdicts. Reported as
descriptive: median S **−0.0007 [−0.0047, +0.0016]** against a TF-shuffled null of +0.0008 ± 0.0017 (z −0.83); 16 TFs
below −0.05 vs 6 above +0.05 (the only sign of the predicted reversal; the pluripotency network NANOG / POU5F1 / MYC /
KLF4 / TCF3 / TEAD2 and core transcription machinery on the mirrored side, MYC's partner MAX on the aligned side);
C6b ρ −0.114 [−0.285, +0.067] (the wrong sign); C6c +0.058 [−0.111, +0.228].
The low ORF split-half value (+0.025) is a limitation, not a measured explanation
of the null or a valid latent-mirror bound: the KD reliability was borrowed from
another bank and mixed a half-sample with a full-sample convention.

**304 audit, independently reproduced by `scripts/audit_round5.py`.**
The c ≈ 9,600 claim used different TFs in depth bins, not within-TF learning
curves; the saved-object ladder does not support that prescription. Its
normalization/centering assumptions also need controls, so neither the
universal law nor its universal refutation is claimed. Isoforms encode
different proteins and their +0.021 agreement is not an exchangeable-construct
noise ceiling. The +0.3823 geometry association is largely present under
shuffled TF pairing: a separate 2,000-draw uniform null has mean +0.2894,
95% range [+0.1525, +0.4266], plus-one p=.09645. A target-specific
cross-modality geometry effect is not established.

**Limit.** Both atlas and KD bank are H1, but common assay state has not been
established; Joung's atlas used seven days in differentiation medium. The
observed near-zero median is not biological equivalence or proof of
orthogonality. C6's frozen validity failure remains; fixing interpretation
does not rehabilitate the claim or reopen R4.

**Next measurement.** A second OE library can test another observed
cross-study comparison. It cannot identify measurement reliability or isolate
modality unless its line, state, assay and replication support that inference.

**303 measured that second comparison (C7, 06:01Z September 7).** The complete
protocol/code was frozen in commit `45788c6` before the single R5 read.
C7a **PASS**: median matched residual cosine −.02633
[−.05203, −.00949], conditional target bootstrap; 64/100 negative,
41 below −.05 versus 13 above +.05, tail excess .28 [.14, .41].
Both directional matching-null p values are 1/2001, Holm=.001999 across
four cells. C7b/c/d **FAIL**: conservation→reversal rho −.12883,
geometry→reversal +.07606, geometry→absolute overlap +.08070,
all adjusted p=1. The uncentered median remains negative (−.02029),
a declared descriptive sensitivity. All validity gates pass; independent
saved-array reconstruction gives zero difference in the full cosine matrix.
R5 is spent. This establishes modest target-matched negative alignment
in this cross-study resource under target-label exchangeability, not a
universal inverse, a causal explanation for 302's null, or independent
study replication. The novelty of the quantitative result is not asserted.
See `experiments/303-second-overexpression-library/COORDINATOR-RESULT.md`.

### T7 update — sequence-linked repression observations, no full-family discovery (305)

Eleven primary construct–gene pairs were frozen from sequence; eight are
measurable. GRB10_P2 → SNX8 and DPF2_P1 → CEP120 show effects −0.4924 and
−0.6328, CPM ratios .4507 and .3446, with negative effects in both cell halves.
The other construct of each intended target does not show the same repression.
These numbers independently reproduce from saved integer sums. **Both have
full-11 Holm p=.10476; zero of 11 passes.** The empirical rank test's resolution
is insufficient, not evidence that repression is absent.

Two of four unique candidate-bearing NTC pairs also show nominal repression;
their twelve split-evaluations reuse identical control cells. Post-hoc
expression-matched and lane-matched diagnostics help prioritize follow-up
but do not rescue the primary criterion. Their Poisson-binomial significance
assumes an unestablished independent/selection-matched null. Off-target
mechanism, its frequency and the failing-vs-working contrast remain unresolved.
See 305's coordinator review for the proposed independent, construct-level test.

### F16 — In the 2025 VCC CRISPRi library, the genes a construct's protospacer seed matches in their promoter are repressed in that construct's own cells, and the effect follows the PROTOSPACER rather than the intended target (307 claim C8, 309 claim C9; construct as the experimental unit)

**PROVENANCE (320, the audit F16 lacked): NEW 0 / EXTENDED 1 / KNOWN 4 / CONTRADICTED 0 — F16 is a
RE-MEASUREMENT of a known phenomenon on a new sample.** The phenomenon (**KNOWN**: Rohatgi *et al.*, *Cell
Genomics* 4(11):100693, 2024 — screens + RNA-seq **+ dCas9 ChIP-seq occupancy**, 98–1,141 reproducible peaks per
guide, **in our constructs' own parent hCRISPRi-v2 library**, two years earlier; Hartman *et al.* 2026 across three
genome-wide Perturb-seq sets); the 12-nt-seed + NGG promoter rule, its window and its **measured low precision**
(**KNOWN**: Rohatgi 2024, Hartman 2026); the partner-construct specificity control (**KNOWN**: Jackson *et al.*,
*Nat Biotechnol* 21(6):635–637, **2003** — the founding observation of the field, and already used for exactly this
purpose on CRISPRi Perturb-seq); aggregate-finds-what-per-pair-cannot (**KNOWN**: Goeman & Bühlmann 2007). Only the
**design** is **EXTENDED** — the construct-level unit with a joint null over covariate-matched candidate *sets*,
new only in using an **empirically measured** 209,070-spacer off-target propensity as the matching covariate
(a recombination of Cazares 2024 and Davis 2023, not a new method). **Nothing biological survives as new.**
**Two further corrections the audit forces:** the headline p is a **gene-sampling** p, not a construct-level one,
which weakens the "construct as the experimental unit" framing; and the partner-construct control **must not be
claimed as ours**. The one item the audit could not find published is the **Horlbeck supplementary-table
Excel-escaping defect** — data hygiene, not science.

**Claim (preregistered twice, hashed 2026-09-07 16:00:01Z and 16:27:26Z, each written into the ledger before any
repaired primary value was computed).** Over **all 34** gene-targeting constructs of the 2025 library, the genes
selected by 301's seed rule (longest seed in {12, 11, 10, 9} giving a hit in the CRISPRi active window after a 10 kb
distal filter) are more repressed in the cells carrying that construct than sequence-, promoter- and
expression-matched genes the rule did not select. **ONE verdict cell**; the unit resampled is the construct (n = 34,
83 construct–gene pairs, no construct dropped for lacking a strong candidate).

**The null is the design.** Pseudo-candidate *sets* are resampled — not genes independently — matched jointly on the
measured probability that a real hCRISPRi-v2 protospacer hits each gene off-target under 301's exact rule (a
209,070-spacer rescan), promoter GC, TSS count and density, a mappability proxy, and basal control expression, with
every library-targeted gene, the construct's own target and its own spacers' genes excluded. Attainable floor
p = 1.0e-5, three orders below α — the defect that made 305's per-pair design unresolvable (its floor was
11/158 ≈ 0.0696 after Holm).

| | C8 (307) | C9 (309, calibration-repaired) |
|---|---|---|
| `S` per construct | **+0.786** [+0.387, +1.166] | **+0.71293** [+0.32419, +1.07351] |
| Δ log1p-CP50k | −0.0790 [−0.1366, −0.0310] | −0.078205 [−0.135970, −0.030017] |
| joint-null p | 1.0e-5 (floor) | **1.99998e-05**; max leave-one-construct-out 8.99991e-05 |
| known-answer rows | **2 of 7 FAILED** → `UNRESOLVABLE` | **9/9 PASS** → **PASS** |

**Why there are two rows.** C8 met both pass conditions but **aborted itself** on two calibration defects — a
construct-independent bias in the no-signal row (+0.0683 vs a ≤ 0.050 rule; +0.0717 ± 0.0022 at 5,000 derangements)
and one matching axis of nine (promoter-window count, SMD +0.2905 vs ≤ 0.250) — and the frozen abort clause was
honoured rather than loosened. C9 repaired exactly those two (scores centred on each construct's own
wrong-recipient distribution; promoter-window count added to the hard strata) and re-tested the identical
hypothesis at the identical α with the identical single verdict cell. **The bias repair works AGAINST the effect and
was declared to do so before running: it moved the estimate −9.5 %, from +0.78797 to +0.71293, and did not narrow
the null (sd 0.17316 vs 0.17173).** It still cleared.

**What else survived.** 26 of 34 constructs positive (exact binomial p = 0.0015). A structurally different
**real-spacer resampling null** — a pseudo-construct is a random pair of stringency-0 hCRISPRi-v2 protospacers
rescanned through 301's exact rule — gives +0.67242, p = 0.00149 at its own floor. **The specificity test is the
design's best idea:** the identical candidate genes read at wrong-recipient level in the cells of the **partner
construct that knocks down the SAME intended gene with DIFFERENT protospacers** — own +0.78797 vs partner +0.04296
[−0.29799, +0.38107], difference **+0.74500 [+0.22896, +1.25766]**, paired Wilcoxon p = 0.00459.

**Limit — and it is a hard one.** (i) **C8 and C9 are two tests of one hypothesis on one data set**, not a
replication; the repaired estimate was predictable in advance from C8's own post-hocs (P5 +0.71399), which is why
C9 buys "not an artefact of the two defects" and nothing more. (ii) **This does NOT establish dCas9-KRAB off-target
binding.** Expression is observational and nothing here is an occupancy assay. (iii) The 34 constructs differ in
on-target knockdown (−3.04 … +0.03), cell number (95 … 3,415), composition and handling; **no construct has a
biological replicate**, and all three splits share one control pool. (iv) The partner-construct test narrows the
explanation to **something carried by the protospacer sequence** but **cannot separate seed-mediated binding from a
construct artefact correlated with spacer identity**. (v) With the construct as the unit the interval is set by
n = 34: more cells do not narrow it, more constructs would.

**What it changes.** T7 moves from "not testable on disk" (298) and "no enrichment among failing constructs" (301,
which stands — 307's KA-2 reproduced its selection exactly: failing 4/13, working 4/21, NTC 4/31) to a measured
aggregate association with the construct as the unit. **301's negative and F16 are not in conflict:** 301 asked
whether candidates are enriched among *failing* constructs and found they are not; F16 asks whether candidates are
repressed *in aggregate across all constructs* and finds they are. The failing-vs-working contrast remains
unsupported. The decisive next experiment is not more analysis of these cells: it is a **sequence-mutant control**
(re-clone the candidate-bearing constructs with silent seed mutations that abolish the predicted match while
preserving the on-target site; ~24 constructs × ~500 cells, one 10x lane, order $3–6 k) or a **dCas9 occupancy
assay** ($5–10 k). Recorded as ledger claims C8 and C9.

**A defect found in a published resource, worth carrying.** 307's spacer calibration found an **Excel-escaping
artefact in Horlbeck et al.'s supplementary gene labels (`="SYMBOL"`)** that silently stopped an own-target
exclusion from firing and inflated the per-spacer L = 12 hit rate from **0.107 to 0.597**. A known-answer row caught
it before use.

## 3. What is and is not established

| finding | unit resampled | n | independent replication | confirmable on the reserved data (preregistered test) |
|---|---|---|---|---|
| F1 Orion/K562 reliability medians 0.03–0.05, a third ≤ 0 | batches or cells of one experiment | 269–294 targets × 3 lines | batch and random splits agree on 2/3 lines; the donor side (R271) repeats it | the HEK293T half-B reliability was the item-3 access; no further test needed |
| F1 Flex-400 reliability +0.243, 0/85 ≤ 0 | lanes / cells | 85 targets, 1 line | lane and random halves agree (Spearman +0.994); 159's rows reproduced at 0.0 | no reserved Flex data; the board's replicate anchor (r 0.927–0.984 on `pds`) is the only external check |
| F1 guide construct +0.014 vs +0.243 | construct | 5 (H1) + 8–10 (Orion, K562) | same sign in four sets, all tiny n | no |
| F2 the oracle disappears (retained −0.32) | targets; halves | 3 folds × 259 | both split kinds, both directions, τ0 gate; 265's permutation on the cancer folds | this was the single half-B scoring; a second access needs a new claim |
| F2 reliable-tercile information +0.021 / +0.027 | targets | 285 | replicated on the random split (p 0.001 after a 12-stratum correction) | yes, on the random-cell half B of HEK293T (the batch half was already read): "cross − random on T3 > 0" |
| F3 spectrum 77–79 % / 87–92 % | directions | 4 operators | the same on all four | not needed (a property of the fitted objects) |
| F3 anti-correlation is in-sample noise | donor halves | 2 folds | both donors, both folds | not applicable (the reserve is a truth, not a donor row) |
| F3 damping loses `pds` on PLOLO, null on Flex | targets | 3 folds; 85 | PLOLO 3/3; Flex one line | yes: "Wiener κ = 1 ≤ −0.01 and floor c = 1 within ±0.005 on HEK293T half B", one access, both arms |
| F4 knee at 0.02; deficit 53 % anti-aligned / 47 % marginal | targets | 853 + 85 | two references, two decoy counts; overlap 4–5× chance | yes, in the same access: mean `pds` ≥ 0.97 above own cosine 0.02 on half B |
| F5 line-specific component 0.115, 4× null | cells; guides | 82 targets × 12 lines | guide halves keep 79 % | yes: the four reserved lines' own cell halves, predicted line-specific fraction 0.04–0.14 and > 3× the OR null |
| F5 not basal-predictable (zero full-response gain) | lines | 12 folds | two model classes, two controls | yes: predict the reserved lines from the 12-line fit; test residual cosine > wrong recipient and full-response Δcos ≤ +0.005 |
| F6 K562 DE sparsity | cells / control draws | 272 targets; 9 control draws | 266 reproduced 260's counts; strong-target known answers | no reserved data; a second line would be the replication |
| F7 `nmae` identity ×1.50 | one board point | 1 | no (two earlier small-magnitude mis-signs) | only the board |
| F7 class factors k_pds / k_cos | one board point per class | 3 | no | only the board |
| 270's cross-study sign-stability feature (AUC 0.89 [0.77, 0.97], ρ +0.36 [+0.00, +0.63], 40 targets, 4 anti-aligned) | targets | 40 | no (its CMap analogue is at chance) | held targets rather than reserved data: preregister before computing it on the remaining panel |

## 4. Implications

**Competition track.** Every knob-level closure since BR249 (R257–R264, R271–R272) was read on truths whose per-target
reliability is 0.05–0.09; only fold totals and group aggregates were readable there, and a +0.002 `pds` arm sits inside
that noise. The board's regime is Flex at 400 cells, which resolves every target: the incumbent's median own cosine on
H1 is 0.037, 13 % of targets point away and 48 % sit on the knee, so the lever is the transfer cosine on the marginal
band. Two cosine effects are visible on the reliable reference and invisible on PLOLO: the soft noise floor (+4.4 %
`cos_agg`, `pds` null) and the un-deflated CD4 residual (+30 % `cos_agg` against +24 % deflated); neither is priced, and
the board class factor of an operator change is unmeasured. No selector, router or reliability weighting buys `pds`;
the +0.046 ceiling needs an arbiter that knows the held line. Pricing stays coherent-pessimistic because k_pds spans
0.7–3.8 and DE-three 0.04–0.19 across three points; the `nmae` identity is a validated lower bound on a loss, and no
local DE surface exists in the board's regime.

**Scientific track.** F1 is a reliability law for a non-essential panel: at 150–200 cells a third of Orion targets have
no measurable direction, and Flex at 400 cells resolves all 85 H1 targets by chemistry and line rather than depth
(+0.114 at matched depth); the guide construct is a unit of variance the challenge truth does not average over. F2 is a
methodological result: a per-target choice validated on the truth it was chosen on inflates by that truth's noise, and
the retained fraction under independent halves is the statistic to report. F3 shows an operator fitted at n/p ≈ 1.6
encodes its rows' noise realisation and that a rank-retrieval metric rewards exactly those directions. F4 makes the
metric's threshold explicit and splits the deficit into an unmeasurable-response tier and a marginal band. F5 puts
context dependence at the size of the shared response, high-rank and batch-entangled, with predictable direction but no
usable prediction at 11 lines. After the novelty audit (E287 `NOVELTY-AUDIT.md`: NEW 4, EXTENDED 11, KNOWN 7, CONTRADICTED 2): defensible as new — F3b (rank retrieval needs the noise-dominated directions of its own whitener, against the shrinkage-whitening prescription), F11c (conservation does not predict within-individual state dependence; preregistered negative), F10b (the donor as a marginal unit; replicates vs cells at matched budget), F9 as a preregistered quantification of a qualitatively known axis (hold F12a until 284); re-measurements credited to their sources — F1a, F1c, F6, F10a, F11a; not scientific claims — F7, C1, F11a's ratio, F4a's knee as a discovery. The structural point no single finding states: within-reference split-cell repeatability (+0.243), across-donor agreement (+0.054) and across-context transfer (≈ +0.03 on the panel) are three different quantities the literature conflates, and they are numerically far apart on the same cells.

## 5. Methods appendix

**5.1 The statistic.** For target t with prediction p_t (bulk log1p(50,000 · P / ΣP) minus the same profile of the real
control cells, the panel's target-gene columns removed; on PLOLO rendered through production's `profile_from_effect`
into the challenge control frames and read back as the scorer does), `pds_t` = [#{q ≠ t : cos(p_t, y_q) < cos(p_t, y_t)}
+ ½ #{ties}] / (n − 1) over the held truths y, averaged over frames on PLOLO; a fold's `pds` is the mean over its n
targets and one rank step is 1/(n − 1) per frame (R265 §1, R273 §1). Own cosine is cos(p_t, y_t) on the scored columns;
`cos_agg` is R243's energy-weighted aggregate cosine.

**5.2 Split units.** Batch or gemgroup parity (even vs odd) with controls split by the same rule; random cell halves
(one draw per cell, seed 20260905); guide constructs P1 vs P2 where two exist; Flex lanes from the barcode suffix (48).
Residual space: each side's profile minus the mean over that side's targets (the bank-wide mean on PLOLO's instrument);
Spearman–Brown 2r/(1 + r) to the union of the halves. Broad-16: the cells of each (line, KO) split at random, each half
scored against the same half of the ten OR-control KOs; noise is the energy the other half does not reproduce, halved to
convert half depth to full depth (checked at 0.4996).

**5.3 Intervals, permutations, multiplicity.** Paired target bootstrap, 20,000 draws, stratified by fold, percentile
[2.5, 97.5] and P(Δ ≤ 0); 268 uses 4,000 draws over targets carrying their 12 folds. Permutations: 265 draws 100,000 drop
assignments (uniform and frequency-preserving), p = (1 + #{null ≥ observed}) / (B + 1), on an additive proxy whose total
error is 0.04 null SD on the iPSC fold and 0.18–0.27 on the matched cancer folds (refused above 1 SD); 268 permutes batch
labels over lines and the rows of half B for component energies; 162 is exact over 720 relabelings. Multiplicity is
reported as cells screened (267: 12 strata × 2 splits; 263: ≈ 130 cells; 271: 10 + 6; 272: 11 + 9; 273: 26 arms).

**5.4 Pricing (competition track).** Six members per arm: `pds` = 2.2138 × k_pds × Δ`pds` / 6; `mse` through R197's
identity at the arm's own amplitude terms with Δcos_board = k_cos × relative `cos_agg` × 0.10547; DE-three at a rate per
unit Δcos; `nmae` by the R216/R188 identity / 0.60 / 6. Pessimistic scenario: per leg the band end that lowers the total
for that cell's sign — k_pds {0.70 or 0.85 on a gain, 1.22 or 2.48 on a loss}, k_cos {0.94 on a gain, 3.76 on a loss},
DE-three {0.0425, 0.24}, identity ×1.5 on a loss and ×0 on a gain; gates +0.0019 (shape-only, delivered amplitude within
0.5 %) and +0.0089 (amplitude-changing); a cleared cell also needs 3/3 folds, leave-top-3-out retention ≥ 50 % and every
control below half the arm. Three labels apart: regression check (0.0e+00 on cached rows), fitted calibration,
independent validation.

**5.5 Leakage rules.** The H1 union stays off on h1-test (`h1_2025` covers 85/85 of its targets); PLOLO operators are
fitted without the held line and every preprocessing step is refit inside the training split; a same-target lookup
across PLOLO folds leaks the held line's own response through the training folds' donor means and is excluded (R269 §5:
+0.033 with the leak, −0.014 [−0.023, −0.005] without); reserved data (the HEK293T half-B truths, four Broad-16 lines)
is scored once per preregistered claim and logged in `docs/CONFIRMATION-LEDGER.md`; a known-answer row precedes every
screen.
