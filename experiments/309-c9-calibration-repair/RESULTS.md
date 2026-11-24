# 309 — RESULTS: Claim C9, the calibration repair of C8, re-tested on the same hypothesis at the same α

Preregistration `PREREGISTRATION.md`, sha256
`2fdf4531f3c319a052bd73eb3c76326954e2f75362c1a0540b3974288b8febca`, hashed **2026-09-07 16:27:26Z**
and appended to `docs/CONFIRMATION-LEDGER.md` as **Claim C9** before any repaired primary value was
computed. No reserve declared or spent; no h5ad opened; 305's saved integer sufficient sums and
301's/307's saved sequence products reused throughout. Peak RSS **1.25 GiB** (budget 4 GiB); total
compute 20 s.

---

## THE LINE

**C9 (aggregate candidate-repression enrichment, construct as unit, calibration-repaired): PASS**
— `S` = **+0.71293** standardized units per targeting construct, construct-level bootstrap 95 %
**[+0.32419, +1.07351]** (n = 34 constructs); joint-null **p = 1.99998e-05** (1 of 100,000
selection-matched joint draws reached the observed value; attainable floor 1.0e-5); maximum
leave-one-construct-out p = **8.99991e-05**; **nine KA rows 9/9**. **C8's estimate was
S = +0.786 [+0.387, +1.166]** and it stands, unretracted, beside this one.

`Δ` = **−0.078205 log1p-CP50k [−0.135970, −0.030017]** per construct. Both preregistered pass
conditions were met and no known-answer row failed, so §8's abort clause did not fire.

---

## Abstract

**What was repaired, and why that is legitimate.** C8 measured one aggregate statistic — are the
genes 301's seed-match rule selects for a construct more repressed in that construct's cells than
sequence-, promoter- and expression-matched genes the rule did not select? — with the **targeting
construct** as the unit, and got **S = +0.786 [+0.387, +1.166]**, joint-null p at its floor. It then
**correctly aborted** to `UNRESOLVABLE-AT-THIS-RESOLUTION` because 2 of its 7 preregistered
known-answer rows failed: **KA-4**, the no-signal row, measured a construct-independent bias of
**+0.0683** against a `≤ 0.050` rule (**+0.0717 ± 0.0022** at 5,000 derangements — 9 % of the
observed effect), and **KA-6**, matching balance, missed on one axis of nine (promoter-window count,
SMD **+0.2905** vs `≤ 0.250`). Those are defects of the null and the matching, not of the
hypothesis. C9 repairs exactly those two and re-tests **the identical hypothesis, at the identical
α = 0.05 one-sided, with the identical single verdict cell and the identical statistic family**.
No threshold was loosened, then or now.

**The repairs, and their direction.** (i) The bias repair centres every construct's score on **its
own wrong-recipient distribution** — the same candidate genes, in the same matched pools, read in
the cells of every other targeting construct of the same split (15 / 9 / 7 recipients). It is C8's
partner-construct test generalised from one wrong recipient to all of them. (ii) The matching repair
adds promoter-window count to the hard strata (10 basal deciles × 5 `pi_off` quintiles × 3 window
classes = 150 cells), plus one global exclusion needed to make a pool recipient-independent.
**Repair (i) works against the effect**, as the preregistration stated before running: it moved the
estimate from **+0.78797 (uncentred, repaired pools) to +0.71293**, i.e. **−9.5 %**. It did not
narrow the null (repaired null sd 0.17316 vs C8's 0.17173), so the p rose from the floor 1.0e-5 to
2.0e-5 — still 2,500× below α.

**The number, and what survived.** `S` = **+0.71293 [+0.32419, +1.07351]**, p = 1.99998e-05,
max leave-one-out p = 8.99991e-05, 26 of 34 constructs positive (exact binomial p = 0.0015). Two
structurally different nulls agree: the **real-spacer resampling null** (a pseudo-construct is a
random pair of stringency-0 hCRISPRi-v2 protospacers rescanned through 301's exact rule) gives
**+0.67242, p = 0.00149** at its floor; the **wrong-recipient permutation** gives p = 2.0e-4 at its
floor — though under this repair that permutation is the same object as KA-4R and is *not* an
independent check (§5.2). The **partner-construct specificity** test on the repaired pools is
unchanged from C8: own **+0.78797** vs partner **+0.04296 [−0.29799, +0.38107]**, difference
**+0.74500 [+0.22896, +1.25766]**, paired Wilcoxon p = 0.00459, n = 34.

**This is a SECOND TEST OF THE SAME HYPOTHESIS ON THE SAME DATA.** Same 34 constructs, same cells,
same sums, same candidate genes. It is **not** an independent replication and nothing here should be
read as one. Across C8 and C9 there are two verdict cells on one hypothesis and one data set;
Bonferroni/Holm over the two attempts requires the smaller p ≤ 0.025, which 2.0e-5 satisfies
arithmetically — but the arithmetic is not the interesting part. The real multiplicity is
**analytic**: C9's design was chosen after seeing C8's diagnostics, and the repaired point estimate
was **predictable in advance** from C8's own published post-hocs (P5: +0.71399; P2: +0.78813). The
preregistration says so explicitly. What C9 buys is that the estimate is **not an artefact of the
two calibration defects C8 found** — it does not make one experiment into two, does not narrow the
interval (the unit is still 34 constructs), and does not upgrade an association into a mechanism.

**What a positive does and does not establish — stated here, not only in the limitations.** A PASS
says only that genes selected by a 12-nt-seed + NGG promoter match are, in aggregate, more repressed
in the cells carrying that construct than matched genes the rule did not select, over and above the
level those same genes read in the other constructs of the split. **It does NOT establish
dCas9-KRAB off-target binding.** Expression here is observational; nothing measured is an occupancy
assay; the 34 constructs differ in on-target knockdown (`kd_lfc` −3.04 … +0.03), in cell number
(95 … 3,415), in population composition and in handling; there is no biological replicate of any
construct; and all three splits share one control pool (KA-3). The partner-construct test narrows
the explanation to **something carried by the protospacer sequence** — and it **cannot separate
seed-mediated binding from a construct artefact correlated with spacer identity.** A FAIL would have
meant "no aggregate expression corroboration under this test", never "off-targets do not occur".

---

## 1. The nine known-answer rows — 9/9

| id | row | result | rule | |
|---|---|---|---|---|
| KA-1 | 305's 23 published rows (pooled effect, both half effects, control CP50k) recomputed from `sums305_<split>.npz` | max \|diff\| **0.000e+00** | ≤ 1e-9 | **PASS** |
| KA-2(a) | 301's selection reproduced from its saved hits | failing **4/13**, working **4/21**, NTC **4/31**; **8** L = 12 targeting constructs; L\* = 12/11/10/9 for 8/16/7/3; **83** pairs over 34 constructs | exact | **PASS** |
| KA-2(b) | **the inherited pool is the CORRECTED one** (Horlbeck `="SYMBOL"` Excel escape) | **0** escaped labels of 209,070; **0** own-target-in-own-hit-set violations over all four L; 68/68 real spacers found verbatim; pool L = 12 rate **0.160879** and real-spacer rate **0.132353**, both at \|diff\| **0.0e+00** vs C8's saved values; **bug-equivalent rate 0.597** | see §9 of the prereg | **PASS** |
| KA-3 | the 31 NTC constructs' per-half sums across the three splits | identical names; max \|sum diff\| **0**; max \|cell-count diff\| **0**; **38,176** control cells | exact | **PASS** |
| KA-4 | **unchanged machinery** — C8's no-signal row on C8's pools, 5,000 within-split derangements | mean **+0.07168** (SE 0.00222), 95 % [−0.23670, +0.38396]; C8 published **+0.0717** | reproduction within ±0.010 | **PASS** |
| KA-5 | null calibration, 200 draws each p-valued against the rest | **repaired null** KS p **0.1794**, fraction ≤ 0.05 = **0.0350**; **C8's null** KS p **0.2456**, fraction **0.0600** | KS p ≥ 0.05 and fraction ∈ [0.02, 0.10], on both | **PASS** |
| KA-6 | **unchanged machinery** — the nine-axis SMD table on C8's pools | reproduces C8 to max \|diff\| **0.00e+00**, including `log1p_n_windows` **+0.2905** | reproduction within 0.002 | **PASS** |
| KA-7 | 293's on-target `kd_lfc` for all 34 constructs | max \|diff\| **2.220e-16** (worst `DAXX_P2`) | ≤ 1e-9 | **PASS** |
| **KA-4R** | **repaired** — the no-signal row on the repaired statistic, 5,000 within-split derangements | mean **+0.00290** (SE 0.00222), 95 % **[−0.30574, +0.31191]** | **C8's original rule**: \|mean\| ≤ 0.050 **and** the interval brackets 0 | **PASS** |
| **KA-6R** | **repaired** — the nine-axis SMD table on the repaired 150-cell pools | worst axis `map_mult` **+0.2460**; `log1p_n_windows` **+0.0000** | **C8's original rule**: \|SMD\| ≤ 0.250 on every axis | **PASS** |

**Two honesty notes on these rows.**

- **KA-4 (unchanged) still fails C8's original `≤ 0.05` threshold** — +0.0717 against ≤ 0.050. That
  is the point: the row governs in C9 by **reproduction**, because it is the defect being repaired,
  and the preregistration said so before it was run. C9 does not re-litigate the threshold C8
  aborted on.
- **KA-4R passes by construction in expectation.** The repaired statistic subtracts precisely the
  level this row measures, so its expectation is exactly zero. Its evidential content is (i) the
  correction is correctly implemented and (ii) Monte Carlo confirms the interval brackets zero
  (mean +0.00290 against an SE of 0.00222). The preregistration declared this cost **before** the
  repair variant was fixed, and it is the reason the structurally different checks — KA-5, the
  real-spacer null, the partner test — are kept and reported.
- **KA-6R passes but one axis is marginal.** `map_mult` (the promoter-corpus 12-mer multiplicity
  proxy) reads **+0.2460** against a `≤ 0.250` rule — inside the rule with 1.6 % to spare, and
  essentially unchanged from C8's +0.2417. The axis the repair targeted, `log1p_n_windows`, went
  from +0.2905 to **0.0000**, and `gc` (−0.1289 → −0.0251) and `log1p_n_sites` (+0.1141 → +0.0012)
  improved as side effects. Mappability remains the least well matched covariate in this design.

**Regression of C8's own primary.** C8's uncentred `S` rebuilt from C8's pools by code written here
reproduces **+0.7856650218824774 at \|diff\| 0.0e+00**, with p_joint 9.9999e-06 — C8's published
values exactly.

**KA-0 (reported, not a verdict row) — an audit finding in C8's write-up.** The 34 targeting
constructs split **train 16 / test 10 / validation 8** (13 failing + 21 working). C8's `RESULTS.md`
§1 and `PREREGISTRATION.md` §2 both state "train 21 / validation 7 / test 6". C8's own
`constructs307.csv` — the object every C8 computation actually read, including its within-split
derangements — carries 16/10/8. **The C8 numbers are unaffected** (the CSV is what its code used);
the sentence in its write-up is wrong. C9's within-split centring uses the measured 16/10/8, giving
15 / 9 / 7 wrong recipients per construct.

## 2. The two repairs, measured

### 2.1 The bias repair (KA-4), and why variant (a) was the right choice

The preregistration chose **variant (a)** — per-construct centring — over **variant (b)** — scalar
subtraction — and justified it from the diagnostic: C8's P3 put the per-construct spread of the
wrong-recipient level at [−0.295, +0.369] around a mean of +0.07, so a single scalar would leave
each construct's own null level inside its contribution and would not repair the construct-level
interval, which is the deliverable when the construct is the unit.

**Measured here, the justification holds and is stronger than the argument assumed:** over the 34
constructs the wrong-recipient level `b_c` has mean **+0.07504**, **sd 0.42538**, range
**[−0.8093, +1.2951]**. The per-construct spread is **5.7×** the mean it scatters around.

**The two variants have exactly the same point estimate**, by linearity:
`mean(z) − mean(b) = mean(z − b) = +0.71293`. They differ only in the interval and in the null. So
the choice between them **could not have been made for its effect on the answer** — it has none.
Variant (b)'s propagated interval is **[+0.18217, +1.23115]** (the uncentred bootstrap CI shifted by
`b̄` and widened by 1.96 × SE(`b̄`) = 0.14298); variant (a)'s is **[+0.32419, +1.07351]** — narrower,
because pairing each construct with its own wrong-recipient level removes shared variance.

### 2.2 The matching repair (KA-6), and the attainable resolution

| axis | C8 SMD | C9 SMD (repaired pools) |
|---|---:|---:|
| `pi_off_Lstar` | +0.0010 | −0.0061 |
| `basal_log1p_cp50k` | −0.0266 | −0.0292 |
| `gc` | −0.1289 | −0.0251 |
| `log1p_tss_density` | +0.1988 | +0.1602 |
| `map_mult` | +0.2417 | **+0.2460** |
| `lowcomp` | +0.2132 | +0.2146 |
| `mean_offset_100` | +0.0822 | +0.0974 |
| **`log1p_n_windows`** | **+0.2905** (the C8 failure) | **+0.0000** |
| `log1p_n_sites` | +0.1141 | +0.0012 |

**Attainable resolution at the finer stratification, reported rather than papered over.** Over the
34 constructs the repaired pools are: **83 pools, mean size 180.4, minimum 25, maximum 200**;
**3 of 83** needed the declared widening to adjacent basal deciles, and those three are named —
`DPF2_P2`/`ITIH3` (stratum cell held 18 eligible genes), `ZNF581_P2`/`RUNX2` (9),
`STAG2_P1`/`MAP3K7CL` (10). **No construct was dropped and none needed to be**: every one of the 34
retains a full candidate set at the 150-cell stratification. The cost of tripling the number of hard
strata is therefore a mean pool 10 % smaller than C8's 196.7, and three pools built from an adjacent
decile — not a loss of constructs.

**The one repair element whose direction favours the effect is worth 0.05 %.** Making the pools
recipient-independent required excluding, globally, every gene that is the `L*` candidate of any of
the 65 constructs (134 genes, ≈ 0.9 % of the eligible universe). Switching it off:
**S = +0.71259** versus the primary's **+0.71293** — a difference of **+0.00034**, 0.05 % of the
estimate. The preregistration flagged this element as non-conservative and bounded it in advance;
the bound is measured and it is negligible.

## 3. The primary

| quantity | C9 (repaired) | C8 (for comparison, stands) |
|---|---|---|
| `S` observed | **+0.71293** | +0.78567 |
| construct-level bootstrap 95 % CI (20,000 resamples, n = 34) | **[+0.32419, +1.07351]** | [+0.38669, +1.16604] |
| joint-null p (one-sided, repression) | **1.99998e-05** | 9.9999e-06 (the floor) |
| attainable floor | 1.0e-5 (B = 100,000) | 1.0e-5 |
| leave-one-construct-out max p | **8.99991e-05** (deleting `GRB10_P2`) | 1.99998e-05 (same construct) |
| null distribution of `S` | mean +0.00134, sd **0.17316**, 95 % [−0.34018, +0.33849] | mean +0.00005, sd 0.17173 |
| `Δ` | **−0.078205 log1p-CP50k [−0.135970, −0.030017]** | −0.079012 [−0.136647, −0.030962] |
| candidates / constructs | 83 genes / 34 constructs | identical |
| constructs with `z > 0` | **26/34**, exact binomial p **0.00147**, median `z^rep` +1.09591 | 24/34, p 0.0122 |

**Both pass conditions are satisfied and no known-answer row failed, so the verdict is PASS.**

**Decomposition of the two repairs (declared Sens-4).**

| machinery | `S` | p_joint |
|---|---:|---:|
| C8 as published (50-cell pools, uncentred) | +0.78567 | 9.9999e-06 |
| KA-6 repair only (150-cell pools, uncentred) | +0.78797 | 9.9999e-06 |
| KA-4 repair only (C8's 50-cell pools, centred) | +0.70998 | 1.99998e-05 |
| **both repairs — the C9 primary** | **+0.71293** | **1.99998e-05** |

The matching repair moves the estimate by **+0.0023** (+0.3 %); the bias repair by **−0.0750**
(−9.5 %), exactly the direction and roughly the magnitude the preregistration stated before either
was run. The whole of C9's movement relative to C8 is the bias repair, and it is a subtraction.

## 4. Declared secondaries (no verdict)

| # | secondary | C9 result | C8 |
|---|---|---|---|
| **Sec-1** | **real-spacer resampling null** — a pseudo-construct is a random pair of stringency-0 hCRISPRi-v2 protospacers rescanned through 301's exact rule, accepted only at matching `L*` and set size; scoring is matching-free (percentile in the whole eligible universe) | **all 34 constructs resolvable**, 672 joint draws (floor 1.49e-3). **Uncentred** `W` = **+0.74441** vs null −0.01084 ± 0.17886, **p = 0.00149**; **centred** `W` = **+0.67242** vs null −0.01335 ± 0.17509, **p = 0.00149**. Both at the null's floor | +0.74470, p 0.00149 |
| **Sec-2** | **wrong-recipient permutation** | p = **2.0e-4** at its floor (0 of 5,000 derangements reach +0.71293). **Under repair (a) this is the same object as KA-4R and is NOT an independent check** — the repair subtracts exactly the level the permutation measures. C8's p = 2.0e-4 stands as the reading on C8's machinery | 2.0e-4 |
| **Sec-3** | **partner-construct specificity** — each construct's candidate genes re-scored in the cells of its **partner construct, same intended target, different protospacers**; uncentred, on the repaired pools | n = **34** (every construct has a partner, and **all 34 partners are in the same split**, so each is one of that construct's wrong recipients). Own **+0.78797** vs partner **+0.04296 [−0.29799, +0.38107]**; own − partner **+0.74500 [+0.22896, +1.25766]**, paired Wilcoxon **p = 0.00459** | +0.786 vs +0.041, diff +0.745, p 0.00459 |
| **Sec-4** | **failing vs working — 301's original T7 question** | failing (13) `z^rep` **+0.8333**, working (21) **+0.6384**, difference **+0.1949 [−0.5584, +0.9434]**, Mann–Whitney **p = 0.7231**. **Still no contrast.** The sign is opposite to C8's (−0.077) and both intervals are far wider than the difference: this axis is unresolved at n = 34, not "positive" now | −0.077 [−0.869, +0.699], p 0.595 |
| **Sec-5** | **the 31 non-targeting constructs**, identical pipeline, **candidate-bearing observations selected by the same rule, NOT negative controls** | `S_NTC` = **+0.43522 [+0.03554, +0.82811]**, joint-null **p = 0.00666** — about **61 %** of the targeting constructs' aggregate and now nominally separated from its own null (C8's reading was +0.286, p 0.054). See §5.3 | +0.2857, p 0.0537 |
| **Sec-6** | **by seed level** | L\* = 12 (n 8) **+0.685**, p 0.0302 · L\* = 11 (n 16) **+0.963**, p 4.0e-5 · L\* = 10 (n 7) **+0.197**, p 0.303 · L\* = 9 (n 3) **+0.658**, p 0.129. **The effect is still not concentrated at L = 12** — the single most awkward fact for a simple seed-binding reading, unchanged by the repair | +0.633 / +1.014 / +0.496 / +0.652 |
| **Sec-7** | **measurement floor** | **33 of 83** candidates sit below 305's RAW CP50k ≥ 0.5 floor; `S` on above-floor candidates only = **+0.65229** — not a low-expression artefact | +0.6907 |
| **Sens-1** | **repair variant (b)**, scalar subtraction with a propagated interval | **identical point estimate +0.71293** (linearity); interval **[+0.18217, +1.23115]** vs variant (a)'s [+0.32419, +1.07351] | — |
| **Sens-2** | **centring over all 33 other targeting constructs** (cross-split) instead of within-split | `S` = **+0.78480 [+0.40999, +1.13698]**, p 9.9999e-06. Cross-split recipients read a *lower* wrong-recipient level, so this variant subtracts less. **The preregistered within-split choice is the more conservative one**, and it is the one that governs | — |
| **Sens-3** | **global `L*`-candidate exclusion off** | +0.71259 vs +0.71293 (difference +0.00034) | — |
| **Sens-4** | **each repair alone** | §3's table | — |

**Per-construct table.** `work/enrichment309.csv` carries, for all 65 constructs: unit, group,
intended target, split, `kd_lfc`, cells, `L*`, candidate count and names, mean matched-pool size,
pools widened, **number of wrong recipients**, `r_c`, **`z_c` uncentred**, **`b_c` (the
wrong-recipient level)**, **`z_c` repaired**, the three `Δ` columns, the per-construct null's mean /
sd / 2.5 % / 97.5 %, the per-construct p, the **partner-construct z**, the number of candidates
below the measurement floor and the minimum candidate CP50k.

Extremes (repaired): `GRB10_P2` **+2.790** (k = 2: `NEUROG3`, `SNX8`; `b_c` = −0.360; its partner
`GRB10_P1`, same intended target, different protospacers, reads **−1.920** on the same gene set),
`HIRA_P2` +2.389, `SMARCB1_P2` +1.754. Bottom: `TMSB4X_P2` −1.835, `SOX4_P1` −1.726, `MED15_P2`
−1.282. Six of 34 constructs clear p ≤ 0.05 individually — the claim is the **aggregate**, not any
single construct, which is why the leave-one-out condition exists and why the worst deletion
(`GRB10_P2`) still gives p = 9.0e-5.

## 5. What this does not establish, what changed, and what did not

### 5.1 Not established (unchanged from C8)

1. **Off-target binding.** No occupancy measurement exists here. The association is with the
   protospacer **sequence**, not with a demonstrated dCas9-KRAB binding event. The partner-construct
   test excludes intended-target downstream biology and generic candidate-gene properties; it
   **cannot** separate seed-mediated binding from a construct artefact that travels with spacer
   identity.
2. **A mechanism graded by seed length.** Sec-6 again puts the largest aggregate at L\* = 11, not
   12, with L\* = 10 unresolved. A strictly seed-length-graded model predicts the opposite ordering.
3. **A failing-construct explanation.** Sec-4 finds no failing-vs-working contrast, with the sign
   flipped relative to C8 and an interval seven times the difference. **T7 stays exactly where 301
   left it: off-targeting is neither proved as the failing-construct mechanism nor excluded.**
4. **Independent replication.** One experiment's cells, one control pool shared by all three splits
   (KA-3), no biological replicate of any construct. The construct-level bootstrap resamples
   constructs, not experiments. **C9 is a second analysis, not a second experiment.**
5. **A population off-target rate.** 83 candidates over 34 constructs and 17 intended targets cannot
   support one.

### 5.2 What C9 changed relative to C8, honestly stated

- It **removed the two calibration defects** C8 identified, and the estimate survived a repair that
  works against it: −9.5 % on the point estimate, a null of essentially unchanged width, p from
  1.0e-5 to 2.0e-5.
- It **did not** add data, add constructs, narrow the biological question, or make the association
  causal. The interval is set by n = 34 and is barely different from C8's.
- It **weakened one of C8's own checks**: the wrong-recipient permutation, which in C8 was an
  independent reading immune to the matched-pool construction, is under repair (a) the same object
  as KA-4R. That loss is a real cost of the repair and is why Sec-1's real-spacer null — which
  resamples *sequences*, not genes, and agrees at its floor both centred and uncentred — carries
  more of the corroboration burden in C9 than it did in C8.

### 5.3 The non-targeting constructs are now nominally positive, and that is a consistency check, not a control failure

Under the repaired machinery the 31 non-targeting constructs read `S_NTC` = **+0.435
[+0.036, +0.828]**, p = 0.0067 — 61 % of the targeting aggregate, and nominally separated from its
own null where C8 read +0.286, p = 0.054. Read this correctly: **these are not negative controls**
(305's coordinator review, and both preregistrations, say so). A non-targeting construct has no
intended gene but still carries two protospacers that can seed-match promoters, and its candidate
genes are selected by the identical rule. If the association is carried by the spacer sequence, a
positive here is what the hypothesis predicts, on constructs for which **intended-target downstream
biology cannot be the explanation at all**. Each NTC is centred against the other 30, so the reading
is a genuine within-group contrast. It is corroborating, not contradicting — and it is not a verdict
cell.

### 5.4 What would resolve it, and what it would cost (unchanged from C8)

- **A sequence-mutant control — the decisive cheap experiment.** Re-clone the ~8 candidate-bearing
  constructs with 2–3 silent seed mutations that abolish the predicted off-target match while
  preserving the on-target site, and read the candidate genes. If repression tracks the seed and not
  the construct, the mechanism is established; if it survives the mutation, the association is a
  construct artefact. ~24 constructs × ~500 cells, one 10x lane, order **$3–6 k** plus cloning.
  This is the only design that separates the two surviving explanations.
- **An occupancy assay** (dCas9 ChIP-seq / CUT&RUN on 3–4 constructs), order **$5–10 k**, answers
  "is it bound", not "is it repressed".
- **An independent screen carrying the same library** would let the identical statistic be
  recomputed with cell line as a second unit — free if one exists (306's remit).
- **More constructs, not more cells.** With the construct as the unit the interval is set by n = 34.

## 6. Deviations from the preregistration

| # | deviation | before or after any repaired outcome |
|---|---|---|
| 1 | None of substance. Every statistic, threshold, seed, draw count and rule in `RESULTS.md` is the one in `PREREGISTRATION.md`. | — |
| 2 | `sens309.py` imports `c9.py`, so the primary is re-run in-process when the secondaries run. It is deterministic (fixed seeds) and reproduces the same values; no number was recomputed with a different seed after being seen. | stated for the record |
| 3 | Sec-1's pseudo-candidate filter is C8's S5 filter verbatim (library targets and the construct's own hits excluded), so a pseudo-candidate is not additionally screened against the new global `L*`-candidate exclusion, while the eligible universe it is scored against is. The asymmetry shifts a percentile by ≤ 1/14,000 and is the same asymmetry C8's S5 had. | fixed in the code before the row was run; noted here |
| 4 | The KA-4 (unchanged) row was run with 5,000 derangements at C9's seed rather than replaying C8's exact RNG stream; the preregistration set the rule as reproduction within ±0.010 for this reason. Measured +0.07168 vs +0.0717. | before |

## 7. Files

| file | what |
|---|---|
| `PREREGISTRATION.md`, `work/PREREGISTRATION.sha256` | the frozen design, hashed 16:27:26Z, in the ledger as C9 before any repaired value existed |
| `work/lib309.py` | 305's sums reader and the frozen response formula |
| `work/pools309.py` | the pool builder in both configurations (C8's 50-cell; C9's repaired 150-cell + global exclusion) and the balance/null helpers |
| `work/ka309.py`, `ka309.json`, `ka309.log` | KA-0, KA-1, KA-2(a)(b), KA-3, KA-7 |
| `work/c9.py`, `results309.json`, `c9.log`, `state309.npz` | the repaired primary, its joint null, KA-4 / KA-4R / KA-5 / KA-6 / KA-6R, Sec-3 … Sec-7 |
| **`work/enrichment309.csv`** | the per-construct table for all 65 constructs |
| `work/sens309.py`, `sens309.json`, `sens309.log` | Sec-1 (real-spacer null, uncentred and centred) and Sens-1 … Sens-4 |

No `vcc`, no `cell-eval2`, no build, no submission package, nothing under `submissions/` or `data/`,
no download, no model, no production change. Nothing was written outside this directory except the
C9 entry in `docs/CONFIRMATION-LEDGER.md`.
