# 390 — The welding distortion is real and worthless; the reference's down-tilt on our CPM block is worth +0.008 raw fid with exact delivery (H1, analytic)

Status: **COMPLETE** (PI-run, 2026-09-13); board arm 392, harness 391. Analytic on 323's saved arrays; no emission, no scorer run.

## Known answers (own code, `price390.py`)
fid = directional precision on called ∩ adjudicable ∩ not-own × min(1, n_pred/n_real): **0.4975545182292535** = 323's
recorded 0.49755452; nmae 0.9508364630 = recorded; reach (348 walk, delivered order) 0.2543472636 = recorded. Exact.

## 1. The welding distortion (mechanism confirmed, value null)
Delivered log2FC = prediction bulk LFC + d_g, d_g = log2(bulk_ctrl / arith_ctrl): per-entry corr **1.000**, slope
**1.00**. 9.3 % of delivered signs differ from the bulk LFC sign; 99.99 % of those have |bulk LFC| < |d_g|.
Control distortion |d_g| median / 90th pct: H1 0.010 / 0.029; A 0.019 / 0.054; B 0.019 / 0.061; C 0.019 / 0.053.

| estimator vs delivered E0 (paired bootstrap, 10k) | Δ fid precision | Δ nmae (+ better) | Δ reach |
|---|---|---|---|
| E7 bulk LFC (distortion removed) | +0.0012 [−0.0006, +0.0031] | +0.0008 [−0.0004, +0.0020] | −0.0147 [−0.038, +0.004] |
| E8 designed signs × |E0| | +0.0011 [−0.0007, +0.0030] | +0.0024 [+0.0011, +0.0039] | −0.0111 [−0.030, +0.005] |
The flips are coin-flip-sized effects; removing them buys nothing. 388's "fid proxy +0.013" used called ∩ reference gate,
not fid's called ∩ adjudicable denominator; on the real denominator it is +0.001.

## 2. Sign rules on the call set with exact delivery (`signrule390.py`), fid raw / nmae raw / reach raw deltas
Base rates: P(ref down) on our calls 0.516 (head 0.508, CPM block 0.518); delivered precision head 0.542, block 0.508.
| rule | Δ fid | Δ nmae (lower better) | Δ reach |
|---|---|---|---|
| R1 every call down | +0.0040 [−0.0035, +0.0117] | **+0.0190 worse** [+0.009, +0.030] | −0.020 [−0.068, +0.026] |
| **R2 head keeps sign, CPM block down** | **+0.0081 [+0.0015, +0.0150]** | +0.0053 [−0.0014, +0.0130] | +0.016 [−0.022, +0.056] |
| R3 cross-fitted (|E0| decile × sign), 20 splits | +0.0079 [+0.0030, +0.0129] | +0.0022 [−0.0013, +0.0056] | +0.019 [−0.013, +0.055] |
R2 reproduces 361's exact-arithmetic +0.0081. 363 lost it through the BULK (coherent down-shift → pds −0.073 and
CPM coupling scrambling untouched genes); an unwelded channel delivers arithmetic-mean signs with per-gene totals pinned.

## 3. Compositional collateral of steering the block (`collateral390.py`)
sum_g A_g ≡ 1e6 per cell, so the block's arithmetic-mass change must be absorbed by every other gene. Targets: every
block gene at ref_lin·2^−ε. The block's mass-weighted delivered LFC is ~0 (median +0.0004; 49.6 % of block genes up), so:
| ε (log2) | shift on all other genes (median [q10, q90]) | head calls flipped (median) | other calls flipped (median [q90]) |
|---|---|---|---|
| 0.002 | +0.0037 [0.0016, 0.0063] | 0 % | 0 % [3 %] |
| 0.01 | +0.0159 | 0 % | 2.7 % [7.6 %] |
| 0.05 | +0.0751 | 0 % | 9.4 % [18.5 %] |
**Use ε ≈ 0.002:** fid reads sign only, nmae prefers near-zero magnitude on coin-flip genes, collateral is negligible.

## 4. Board conversion (regression of score on raw over 920 entries, 2026-09-13)
overall points per raw unit: pds +0.369, mse −0.012 (clamped region excluded), nmae −0.255, **fid +0.561**, reach +0.188,
jac +0.454. R2 on H1 → +0.0045 overall from fid alone if it transfers 1:1.

## 5. Per-gene cross-line sign prior on the called CPM block (`consensus390.py`), exact delivery, H1 reference truth
371's `consensus_prior.npz` P_banks mapped by gene name (covers 83 % of called block genes; uncovered → down).
Rule: down if prior P_down > 0.5 else up; head (top-400 |design|) and non-block calls keep delivered signs.
Block precision (called ∩ adjudicable): delivered **0.5083**; all-down 0.5178; LOLO6 0.5226; K562-essential 0.5215; CLEAN4 0.5223.

| arm (magnitude) | Δ fid raw | Δ nmae raw (lower better) | Δ reach raw | overall est. (0.561 fid − 0.255 nmae + 0.188 reach) |
|---|---|---|---|---|
| all-down (|E0|) | +0.0081 [+0.0016, +0.0148] | +0.0048 [−0.0019, +0.0125] | +0.016 [−0.022, +0.056] | +0.0063 |
| all-down (ε = 0.002) | +0.0081 | **+0.0048 [+0.0016, +0.0081]** | +0.016 | +0.0063 |
| LOLO6 = 6 banks, h1_2025 out (|E0|) | +0.0119 [+0.0060, +0.0179] | +0.0017 [−0.0042, +0.0089] | +0.010 [−0.022, +0.042] | +0.0081 |
| K562-essential alone (|E0|) | +0.0111 [+0.0059, +0.0163] | +0.0019 [−0.0032, +0.0088] | +0.017 [−0.018, +0.054] | +0.0089 |
| **CLEAN4 = K562-GWPS, K562-ess, Jurkat, HepG2 (|E0|)** | **+0.0116 [+0.0057, +0.0177]** | +0.0017 [−0.0047, +0.0090] | +0.018 [−0.017, +0.056] | **+0.0095** |
ε magnitudes cost a resolved +0.0046 nmae with every prior (the new signs are right more often than not, so magnitude
pays): **ship |E0| magnitudes.** H1 is the line where 371 measured the prior weakest (+0.0067 over all-down on h1_2025 vs
Jurkat +0.020, HepG2 +0.016, K562 +0.012); the known downside is a HEK293T-like line (−0.004, debiased bank).
Transfer caveat: the H1 CLEAN4 prior includes no H1 data, but its banks are the same public Perturb-seq family the
reference was built from; A/B/C block coverage 87 / 79 / 77 % (371).

## 6. Variants in the sign-preserving mass-balanced form (`variants390.py`, H1 analytic, board-slope overall estimates)
| variant | Δ fid | Δ nmae (lower better) | Δ reach | overall est. |
|---|---|---|---|---|
| base: CPM block, CLEAN4 (= 392) | +0.0116 [+0.0059, +0.0177] | +0.0024 | +0.018 [−0.018, +0.055] | +0.0093 |
| head + block, CLEAN4 | +0.0088 | **+0.039 [+0.023, +0.058] harm** | **−0.096 [−0.168, −0.032] harm** | −0.023 |
| block, flip only |E0| ≤ target's block q50 | +0.0081 [+0.0052, +0.0111] | −0.0001 | **+0.025 [+0.0005, +0.055]** | +0.0092 |
| block, flip only |E0| ≤ q80 | +0.0115 [+0.0070, +0.0161] | +0.0011 | +0.023 [−0.007, +0.058] | **+0.0105** |
| block, LOLO6 | +0.0119 | +0.0024 | +0.010 | +0.0079 |
| block, K562-essential | +0.0111 | +0.0024 | +0.017 | +0.0088 |
| block, CLEAN4 where |P − 0.5| ≥ 0.02 | +0.0121 | +0.0016 | −0.002 | +0.0061 |
| block, CLEAN4 where |P − 0.5| ≥ 0.05 | +0.0118 | +0.0013 | +0.014 | +0.0089 |
Head signs must stay as designed (the head carries the design's information and the reach front). Keeping the strongest
block genes' delivered sign (q80) is the best variant by a margin inside the noise; candidate for the next dose/variant slot.

## 7. The same prior on genes we do NOT call (`variants390.py noncalled`), mass-balanced
| arm | Δ fid | Δ nmae (lower better) | Δ reach | overall est. |
|---|---|---|---|---|
| non-called universe genes, CLEAN4 | 0.0000 | −0.0011 [−0.0045, +0.0028] | −0.022 [−0.050, 0.000] | −0.0038 |
| non-called, CLEAN4 where |P − 0.5| ≥ 0.05 | 0.0000 | −0.0008 | −0.007 [−0.020, 0.000] | −0.0010 |
| block + non-called, CLEAN4 | +0.0095 | +0.0004 | −0.027 [−0.061, +0.001] | +0.0001 |
Closed: reach walks every reference-significant gene in our order, including the uncalled tail, and the prior's signs
are worse than the delivered ones there. The channel's use stays the called CPM block.

## 8. Reach order boosted by design–prior sign agreement, on 392's signs (`variants390.py order`)
Order = sig-first, then design |H| × (1 + β·[sign(design) == prior sign]); signs fixed at 392's, so fid / nmae / jac are
unchanged by construction. Reference: 392 signs with the |H| order read reach 0.27409 (delivered order 0.27218).
| β | reach | Δ vs |H| order |
|---|---|---|
| 0.25 | 0.27250 | −0.0016 [−0.0067, +0.0017] |
| 1 | 0.27065 | −0.0035 [−0.0094, +0.0013] |
| 4 | 0.26147 | **−0.0126 [−0.0203, −0.0060]** |
| ∞ (agreement first) | 0.25097 | **−0.0231 [−0.0382, −0.0112]** |
Shuffled agreement at β = 1: −0.0102. Closed: agreement is not a confidence signal for the reach front; |H| stays the order.

## Next
Price 371's per-gene consensus sign prior (LOLO, h1_2025 excluded) the same way — it beat all-down by +0.007…+0.024 on
held-out banks but was closed only for its pds price through the bulk, which this channel does not pay. Then the harness
arm on 389's emitter (six members end to end), once 389's S2 steering code is stable.
