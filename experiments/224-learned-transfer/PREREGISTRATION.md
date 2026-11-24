# 224 — PREREGISTRATION

**Written 2026-09-04, after Phase 0's known-answer rows passed and BEFORE any model was fitted.**
Commit this file before `ladder224.py` is run for the first time.

## The question

> Does ANY trained model beat the equal-weight donor mean on held-**LINE** transfer cosine, using
> the data we already hold — and if so, by how much, and does the advantage grow with the number
> of training lines?

## The substrate

Leave-one-**LINE**-out over the donor tables on disk. Two nested surfaces, both mean-removed per
line (the across-target mean *is* the metric's zero), panel target coordinates and intersection
target coordinates dropped from the gene axis, rows unit-normalised before averaging
(AGENTS.md §4.2).

| surface | tables | distinct LINES | targets | genes |
|---|---:|---:|---:|---:|
| **S5** | 6 | 5 (K562, HCT116, HEK293T, RPE1, JURKAT) | 1,592 | 5,002 |
| **S7** | 8 | 7 (+ HEPG2, H1) | 968 | 5,070 |
| **S8** | 9 | 8 (+ DLD1) | 219 | 5,681 |

`k562_gwps` and `k562_essential` are two libraries of ONE line and never straddle a split.
**S7 is the primary surface** because the decision variable is the slope in `k` and S7 carries the
most distinct lines at a workable target count.

## The baseline (rung 0)

For held line `A` and donor set `D` (|D| = k), `P_t = mean_{l in D} unit(C_l[t])`. Scored as
**transfer cosine, macro over targets**: `mean_t cos(P_t, C_A[t])`. `pds` (300-decoy) is reported
alongside as a free secondary and is **not** the gate.

## THE GATE — fixed now, not negotiable later

A model **clears** iff, on **held LINES**:

```
relative cosine gain  =  (cos_model - cos_donor_mean) / cos_donor_mean   >=  +2.9%
```

at the surface's full `k`, **and** the paired bootstrap over target clusters excludes 0.

`+2.9%` is AGENTS.md §8.4's corrected **shape-only** bar (`+2.5–2.9%`; the amplitude-floor bar is
`+11.6–13.6%`). The conservative end of the range is used. A learned transfer model that only
reshapes the direction is a shape-only arm, so `+2.9%` is the right bar for it.

**Held-perturbation gain is NOT the gate and never substitutes for it.** 195 cleared on unseen
perturbations (+0.00904) and failed at the line boundary (−0.00618); a single pooled number hides
exactly that wall. Both are reported separately, always.

## THE STOP RULE

The rungs are fitted cheapest-first and the ladder **stops at the first rung that fails to clear**
unless that rung shows a **positive slope in `k`**. The slope, not the level, is the decision
variable:

```
slope = d (relative cosine gain) / d k,  fitted over k = 1 .. k_max on the primary surface
```

- slope <= 0 and gate not cleared  ->  **STOP. Report the negative. Do not rent a GPU.**
- slope > 0 and gate not cleared   ->  continue up the ladder; a larger model may reach the bar.
- gate cleared                     ->  run every Phase-2 falsifier before any claim.

A bigger model does not fix a flat slope. 184's entire result was that the advantage grew with
lines on Tahoe and did not on CRISPRi.

## The rungs, in fitting order

| rung | model | inputs | why |
|---|---|---|---|
| **R0** | equal-weight donor mean | donor rows | the baseline |
| **R1** | ridge on target features -> response, **replacing** the donor mean | DepMap co-essentiality (1,178-d) + STRING (64-d), PCA-50, Ridge(10) | **132's exact family — the known-answer FAILURE row** |
| **R2** | the same ridge **stacked on** the donor mean (`a·P + b·f(x)`), and as a residual model | donor rows + features | what 132's own RESULTS.md named as the only legitimate reopening |
| **R3** | learned gene-space linear operator `P -> P·M`, ridge in a PCA basis, rank in {64, 256, 1024, full} | donor rows | the direct "trained model replaces the consensus"; subsumes the whitening family; a *regression* objective, where 195 used a *listwise retrieval* objective |
| **R4** | recipient-basal conditioning (184's shape, on CRISPRi): per-gene multiplicative gain `g(b_A)`, and basal-conditioned donor weights | donor rows + recipient basal | 184's operator, on the data it was absent in |
| **R5** | factorised bilinear / small MLP over (target embedding x recipient basal) -> correction to `P` | everything | the highest-capacity rung reachable without a GPU |
| **R6** | larger model on rented GPU | — | **only if R1–R5 show a positive slope in `k`** |

All rungs are fit by **inner leave-one-training-line-out** inside the donor pool: for each training
line `B` in `D`, the input is the donor mean over `D \ {B}` and the target is `C_B`. The held line
`A` contributes **no row to any fit**, at any stage.

## Phase 2 falsifiers — all mandatory before any positive is believed

1. **Wrong-recipient null.** Pair each prediction with a *different* held line's basal / identity.
   184's arm, control and null read +0.0665 / +0.0667 / +0.0663 — identical to three decimals, and
   that is what killed it. A gain that survives the wrong recipient is measuring amplitude, not
   transfer.
2. **Permutation over target labels.** Derange the target correspondence in the fit; the gain must
   collapse.
3. **Held-perturbation vs held-line, reported separately and never pooled.**
4. **Leak audit.** Every donor <-> surface overlap stated explicitly, including
   `data/donors/h1_2025` (85/85 of h1-test) and 222's shipped union.

## What is NOT in scope

No `vcc`, no submission package, nothing written outside `experiments/224-learned-transfer/`, no
`cell-eval2`. The endpoint is transfer cosine on donor tables, which is the direction channel's
own instrument (197's inversion), not the scorer.
