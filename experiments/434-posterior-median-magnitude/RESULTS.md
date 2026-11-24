# 434 — A learned L1 magnitude (conditional median of the truth given oracle-free features): −0.032 nmae on H1, but only −0.007 with panel donors and panel-like training, worse than the plain MDE. Not built

Status: 2026-09-23 10:30Z, PI, analytic on H1 only (no emission, no builder). Scripts `feat434.py`, `pm434.py`,
`l1reg434.py`, `l1panel434.py`, `l1panel434d.py`; the cell table is `work/allcells434.npz`.

**Question.** nmae reads only the reference's significant genes. So the L1-optimal prediction of a cell is the
**median of l given the features and significance**. 426's MDE add is one fixed shape of that idea. Can a fitted
conditional median, using our delivered lfc/q, the donor deltas, control expression and target strength, beat it?

**Set-up.**
- Truth: h1-test's gates. There are 58 targets with at least 10 gate genes, 64,931 cells.
- Base: 426's `base_s0_bps` delivered lfc (the shipped 410b construction), nmae 0.95237.
- Features: |ours|, our −log10 q, control count/cell, donor count, donor sign agreement with our sign, mean/min/max
  donor delta oriented by our sign, mean |δ|, target row norm, target call count, donor row norm, gene variance across
  donor rows.
- Out-of-fold over targets (GroupKFold 5), gradient boosting.

## Results (Δ nmae vs base, H1)

| rule | all 58 targets |
|---|---:|
| plain MDE add (426, c 0.1 cap 0.1 min 0.3) | −0.0098 |
| sign-confidence classifier (AUC 0.589, calibrated) → posterior-median magnitude, with flips | −0.0188 |
| same, raise-only on our sign | −0.0201 |
| **direct L1 regression (quantile 0.5 GBM, weights 1/Σ_p\|l\|)** | **−0.0319** |
| of which: our signs only, no flips | −0.0186 |
| bound: our signs, oracle magnitudes | −0.1081 (nmae 0.844) |

- The regression flips our sign on 26 % of cells, always at small magnitude. Flipping at MDE-sized magnitudes costs
  +0.082.
- Feature subsets: our own statistics alone −0.0069; **donors alone −0.0286**. The gain is donor information.

**The transfer check kills most of it.** 430 found that H1's donor-sign advantage comes from essential-screen donors
(K562-essential, RPE1, Nadig Jurkat). They cover 20 of H1's 58 gated targets and **0 of the 300 A/B/C targets**.
Keeping only the panel-available donors (K562-GWPS, Orion HCT116, Orion HEK293T), on the 38 H1 targets no essential
screen covers:

| model | train on | Δ nmae, 38 panel-like targets |
|---|---|---:|
| plain MDE | — | **−0.0111** |
| full features | all 58 targets | −0.0167 |
| full features | the 38 panel-like targets only | −0.0066 |
| donors only | the 38 panel-like targets only | +0.0029 |

- Trained the way it would have to be for the panel (panel-like targets, panel donors), the learned median is worse
  than the MDE.
- The "train on all" row borrows the relationships from the essential-gene targets. Those are the strong, conserved
  responses the panel does not have.

**Reading.** With the donors the panel actually has, a fitted magnitude/sign model adds nothing over the MDE's
fixed 1/sqrt(count) law. This is the same information limit 430 found for signs. Not built.

## Side finding: where H1's nmae loss sits (base, per target)

| real gate size | targets | our nmae | \|l\|-weighted sign accuracy | nmae with our signs at oracle magnitude |
|---|---:|---:|---:|---:|
| 10–100 | 17 | 0.916 | 0.62 | 0.755 |
| 100–500 | 13 | 0.963 | 0.55 | 0.903 |
| 500–2,000 | 15 | 0.972 | 0.56 | 0.884 |
| ≥ 2,000 | 13 | 0.966 | 0.57 | 0.857 |

- The small-gate (weak) targets are where our signs are best: RFX5 0.87, EPHX1 0.90, CREM 0.90 |l|-weighted.
  Their gate genes are few and large (median |l| 0.3–1.4), and our delivered magnitudes there are 0.02–0.1.
- We call ~2,600 genes on every target whatever its real gate size (10 to 5,484). That is the call-set design (jac is
  the board's best), not a bug.
- The oracle per-target scale (427: −0.0229) lives mostly in those weak targets. What is missing is an oracle-free
  signal that a target is weak **and** that its top genes are its gate. The small-gate targets' high accuracy suggests
  looking at the top of our own per-target ranking. That is a per-gene question, not 427's per-target one.

## 434e: the top of our own ranking is already right and already large (`topk434.py`)
Per target, genes ranked by our q (then |lfc|). On H1's gates:

| our rank | gate cells | sign acc | \|l\|-weighted | share of gate \|l\| | median \|l\| | our median \|l̂\| |
|---|---:|---:|---:|---:|---:|---:|
| 0–9 | 174 | 0.920 | 0.971 | 0.010 | 0.99 | 0.52 |
| 10–49 | 378 | 0.812 | 0.890 | 0.010 | 0.38 | 0.30 |
| 50–199 | 1,434 | 0.718 | 0.792 | 0.021 | 0.24 | 0.20 |
| 200–999 | 7,413 | 0.611 | 0.661 | 0.079 | 0.17 | 0.08 |
| 1,000–2,999 | 17,756 | 0.577 | 0.585 | 0.171 | 0.17 | 0.02 |
| ≥ 3,000 (uncalled) | 37,776 | 0.533 | 0.537 | **0.709** | 0.35 | 0.05 |

- **71 % of the gate's |l| mass sits in genes we rank below 3,000**, i.e. genes we do not call. Our sign there is
  0.537. These are mostly low-expressed genes, whose significant |l| is large because the MDE is large.
- The head is already sized: boosting the top 5–100 genes to M = 0.15–0.6 is flat or worse than the MDE alone
  (−0.0092 … −0.0016 vs −0.0098).
- Direction balance is not the lever either. Per-target P(up) on the gate spans 0.06–0.89 and knowing it would give
  0.612 accuracy vs our 0.559, but nothing we have predicts it: Spearman +0.30 at best (our all-gene up fraction);
  the panel donors −0.13.

## 434f: what the rivals' nmae implies about their signs (H1 mapping)
- Take the true gate l, flip its sign at rate 1 − a (iid), and give every cell the L1-optimal magnitude for that a:
  the per-expression-bin quantile F⁻¹(1 − 1/(2a)).
- The resulting H1 nmae:

  | a | 0.55 | 0.56 | 0.58 | 0.60 | 0.62 | 0.65 | 0.70 | 0.75 |
  |---|---:|---:|---:|---:|---:|---:|---:|---:|
  | nmae | 0.954 | 0.950 | 0.929 | 0.900 | 0.882 | 0.844 | 0.788 | 0.722 |

- Our H1 nmae 0.952 sits at a = 0.56, our measured gate accuracy. The board (0.945) agrees.
- **The rivals' board nmae implies gate sign accuracy of ~0.62 (the 0.87–0.88 pack) to ~0.65–0.66 (cqawesome 0.830,
  Illumina 0.824, chill 0.835).**
- Their fid (0.52–0.55) is only a little above ours, but fid multiplies accuracy by the call yield, so it does not
  bound their accuracy.
- The gap is direction information on the uncalled, low-expressed gate genes. Magnitude rules cannot close it: the
  best one here buys −0.02 on H1 and less on panel-like conditions.

## 434g: an October-only lever. On essential-screen-covered targets, sign the MDE add by the donor mean (`ess434.py`)
- The essential screens (K562-essential, RPE1, Nadig Jurkat) cover 0/300 A/B/C targets but ~30 D/E/F targets each
  (DEF-RUNBOOK 5.-2). So this lever is an exact no-op on the board and only acts in October.
- On H1's 20 covered targets (base nmae 0.9725; our gate sign accuracy 0.548; the six-donor mean's 0.589):

| rule on the covered targets | Δ nmae there |
|---|---:|
| MDE add, our sign (production) | −0.0075 |
| **MDE add signed by the covering donors' mean** | **−0.0242** |
| MDE add signed by the essential screens' mean | −0.0272 |
| our lfc + 1× donor mean (log2) + MDE on the blend's sign | −0.0266 |
| our lfc + 2× / 3× donor mean | +0.043 / +0.156 (the donors' magnitudes must not be added) |

- **The increment over the production MDE is ~−0.017 per covered target.** On a 300-target panel with ~45 covered
  targets, that is about −0.0025 raw nmae, **+0.0007 overall**. It is small but free: no bulk change, no call change,
  and a no-op wherever the essential screens are silent.
- Implementation, if taken: a per-(target, gene) `mde_sign` vector, the covering donors' mean sign where any essential
  screen covers the target and our sign elsewhere. The emitter has a `mde_sign` spec path (426b), but production
  passes only `mde_add`. Queue it in the October knob list (DEF-RUNBOOK), priced on H1, and check it on 432's
  synthetic D/E/F bundle for coverage.
