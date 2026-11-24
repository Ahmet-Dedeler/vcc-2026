# 309 — PREREGISTRATION: Claim C9, the calibration repair of C8, re-tested on the SAME hypothesis at the SAME α

Written and hashed **before any repaired primary value was computed**. The hash of this file is
`work/PREREGISTRATION.sha256`; it is appended to `docs/CONFIRMATION-LEDGER.md` as
**Claim C9 — 309 (calibration repair of C8; no reserve)** before any number produced by the repaired
machinery exists. Every repaired number in `RESULTS.md` is produced by code written after this file
was hashed.

---

## 0. What this is, in one paragraph

307 measured one aggregate statistic — are the genes that 301's seed-match rule selects for a
construct more repressed in that construct's cells than sequence-, promoter- and expression-matched
genes the rule did not select? — with the **targeting construct** as the unit, and got
**S = +0.786 [+0.387, +1.166]** (n = 34), joint-null p at its floor 1.0e-5. It then **aborted to
`UNRESOLVABLE-AT-THIS-RESOLUTION`** because 2 of its 7 preregistered known-answer rows failed:
**KA-4** (the no-signal row) measured a residual construct-independent bias of **+0.0683** against a
`|mean| ≤ 0.050` rule (extended post-hoc to 5,000 derangements: **+0.0717 ± 0.0022**, i.e. **9 % of
the observed +0.786**), and **KA-6** (matching balance) missed on one axis of nine, promoter-window
count SMD **+0.2905** against a `≤ 0.250` rule. **Both are defects of the null/matching machinery,
not of the hypothesis.** C9 repairs exactly those two defects and re-tests the identical hypothesis
at the identical α with the identical single verdict cell. **C8's bounded estimate stands and is
reported alongside C9 whatever C9 says. C9 does not replace or retract it.**

## 1. The three rules this document is bound by, stated before anything is computed

1. **Same hypothesis, same α, same ONE verdict cell, same statistic family.** The hypothesis is
   C8's, verbatim. α = **0.05, one-sided** (candidates *more* repressed). Exactly **one** verdict
   cell. The statistic remains the within-construct matched-pool rank statistic
   `z_c = sqrt(12 k_c)(0.5 − r_c)` aggregated over the 34 targeting constructs. Nothing about the
   effect is redefined to make it easier to detect: the candidate rule (301's `L*`), the response
   formula (305's), the population (all 34 targeting constructs), the unit (the construct) and the
   direction convention are all frozen exactly as C8 had them.
2. **C8 stands.** C9 is reported *beside* C8, never in place of it. C8's line and interval are
   quoted in full in `RESULTS.md` and in the ledger entry.
3. **If a known-answer row fails again, abort again.** §8's abort clause is the same clause, not a
   weakened one, and it is not weakened by the fact that it has already fired once.

### 1.1 This is a SECOND TEST OF THE SAME HYPOTHESIS ON THE SAME DATA — and what that costs

C9 uses **the same 34 constructs, the same cells, the same integer sums, the same candidate genes**
as C8. It is **not an independent replication**, and no sentence in `RESULTS.md` may imply that it
is. Three consequences, stated now:

- **Multiplicity.** Across C8 and C9 there are **two verdict cells testing one hypothesis on one
  data set**. Family-wise control over the two attempts at α = 0.05 requires the smaller p ≤ 0.025
  (Bonferroni/Holm on two). Both designs have an attainable floor of 1/(B+1) = 1.0e-5, so the
  correction cannot change either verdict arithmetically — **but the arithmetic is not the point**.
  The real multiplicity is *analytic*: C9's design was chosen **after** seeing C8's estimate, its
  nulls and its diagnostics. A p-value counts the tail of a null, not the number of analyses that
  could have been run on the way to it.
- **What C9 can therefore buy, and what it cannot.** It can establish that the estimate is **not an
  artefact of the two specific calibration defects C8 found**. It cannot make one experiment into
  two, cannot narrow the interval (the unit is still 34 constructs), and cannot upgrade an
  observational association into a mechanism.
- **The point estimate is largely predictable in advance, and this is disclosed rather than hidden.**
  C8 already published, post-hoc, both repairs' approximate effects: subtracting the wrong-recipient
  level (its P5) gave **+0.71399 [+0.31501, +1.09436]**, and repairing the KA-6 axis (its P2) gave
  **S = +0.78813, p = 1.0e-5**. So C9's repaired point estimate is expected near **+0.71**, and this
  preregistration is **not blind to the outcome's magnitude**. Its value lies in the calibration —
  a properly repaired null, a repaired construct-level interval, a repaired leave-one-out condition
  and nine known-answer rows — not in surprise.

## 2. What a positive would and would not establish (C8's §1, unchanged)

A positive result says: *genes selected by a 12-nt-seed + NGG promoter match are, in aggregate and
with the construct as the unit, more repressed in the cells carrying that construct than sequence-,
promoter- and expression-matched genes that the same rule did not select, over and above the level
those same genes read in the cells of the other constructs of the same split.* It is
**observational corroboration of a sequence-linked association**.

It would **not** establish **dCas9-KRAB off-target binding**. Expression here is observational;
nothing measured is an occupancy assay; the 34 constructs differ in on-target knockdown efficacy
(`kd_lfc` −3.04 … +0.03), in cell number (95 … 3,415), in population composition and in handling;
and there is no biological replicate of any construct. C8's partner-construct test — the same
candidate genes read in the cells of the construct that knocks down the **same intended gene with
different protospacers** (+0.041 [−0.295, +0.369] versus +0.786) — excludes intended-target
downstream biology and generic candidate-gene properties, and narrows what is left to **something
carried by the protospacer sequence**; it does **not** separate seed-mediated binding from a
construct artefact correlated with spacer identity. C9 inherits that limit exactly, and in fact
builds the partner logic into its primary statistic (§6.1), which sharpens the narrowing and does
nothing at all for the separation.

A negative means "no aggregate expression corroboration under this test", never "off-targets do not
occur" and never "off-targeting does not explain failing constructs".

## 3. Population, unit, candidate sets, response formula — all inherited unchanged from C8

| | |
|---|---|
| **Unit** | the gene-targeting construct (programme rule 1); every interval and p is labelled with it |
| **Primary population** | the **34** gene-targeting constructs of the 2025 library (13 failing + 21 working by 293's `\|kd_lfc\| < 1.0`, frozen in 298's `seed298.csv`), covering 17 paired intended targets, each measured in its own split |
| **Candidate set `C_c`** | 301's frozen rule, unchanged: PAM-anchored exact match of the PAM-proximal L nt of either spacer, canonical NGG, inside a CAGE-TSS promoter window at offset −50 … +300, own-target symbol excluded, 10 kb distal filter; `C_c` = the hit set at `L*(c)`, the longest L ∈ {12, 11, 10, 9} yielding ≥ 1 hit. All 34 have a non-empty set (L\* = 12 for 8, 11 for 16, 10 for 7, 9 for 3); 83 construct–gene pairs |
| **Response** | 305's formula, unchanged: `e_c(g) = prof(S_c)[g] − prof(S_ctrl(split(c)))[g]`, `prof(S)[g] = log1p(50000·S[g]/Σ_j S[j])` over all 18,080 genes; `S_ctrl` = the pooled sums of all 31 non-targeting constructs of the split (leave-that-construct-out for a tested NTC construct). Read from 305's saved `sums305_<split>.npz`; **no h5ad is opened, nothing is re-streamed, no download, no reserve access** |
| **Secondary population** | the 31 non-targeting constructs, run identically, **labelled candidate-bearing observations selected by the same rule, NOT negative controls** (305's coordinator review) |

Everything in this table was fixed by C8 and is not touched by C9.

## 4. What has been read at the moment of hashing, and what has not

**Read:** C8's `PREREGISTRATION.md`, `RESULTS.md` and its published numbers (including its post-hoc
P1–P6), its source files (`lib307.py`, `seq307.py`, `pool307.py`, `probe307.py`, `enrich307.py`,
`enrich_pools307.py`, `sens307.py`, `extra307.py`, `ka307.py`) and its saved metadata
(`pool307_meta.json`, `probe307.json`, `ka307.json`, `constructs307.csv`), plus the *structure* of
`poolhits307.pkl` and `seq307_lookup.pkl` (keys, lengths, the first two pool rows). 305's saved
`sums305_<split>.npz`/`.json` have been opened only to confirm they load and to count constructs;
no expression value has been used in any C9 calculation.

**Not computed at hashing time:** no repaired pool, no repaired `z_c`, no wrong-recipient level
`b_c`, no repaired `S`, no repaired null draw, no repaired SMD, no repaired KA row. One
sequence-only fact was measured before hashing and is used in §6.2: the union of the `L*` candidate
genes of all 65 constructs is **134 genes** (83 over the 34 targeting constructs).

**No new data access of any kind is declared or consumed.** This is an extension on data already
used by 293, 305 and 307. **No reserve.**

## 5. Repair 0 — the inherited pool must be the CORRECTED one, and this is a known-answer row

C8's `pi_off` propensity comes from a rescan of all **209,070** hCRISPRi-v2 protospacers through
301's exact rule. That rescan was corrected during C8 for the Horlbeck supplementary sheet's
Excel-escaping of text cells (`="SYMBOL"`), which had stopped the own-target exclusion from firing
and inflated the per-spacer L = 12 hit rate (documented in `pool307.py:123-127` as 0.597 buggy vs
0.107 corrected in development). C9 **inherits `poolhits307.pkl` verbatim** and carries **KA-2(b)**
(§9) to prove it inherited the corrected object: the escape must be absent from every pool label,
no spacer's own target symbol may appear in its own hit set at any L, and the pool and real-spacer
hit rates must reproduce C8's saved values exactly.

## 6. THE TWO REPAIRS — specified exactly, before any repaired primary value is computed

### 6.1 KA-4 / residual bias — variant **(a)**: per-construct centring on the construct's own wrong-recipient distribution

**The repair.** For each construct `c`, let `retained(c)` = every **other targeting construct in
`c`'s own split** (train 16 / test 10 / validation 8 targeting constructs ⇒ 15 / 9 / 7 wrong
recipients; the split composition is itself reported, see §9 KA-0). Score `c`'s candidate set in
each wrong recipient's cells, **using `c`'s own matched pools, unchanged**:

```
q_d(G) = ( #{g in P_G : e_d(g) < e_d(G)} + 0.5 · #{g in P_G : e_d(g) = e_d(G)} ) / |P_G|
z_{d<-c} = sqrt(12 k_c) · (0.5 − mean over G in C_c of q_d(G))          for d in retained(c)
b_c      = mean over d in retained(c) of z_{d<-c}                        (the wrong-recipient level)
```

**The repaired per-construct score, and the primary statistic:**

```
z_c^rep = z_c − b_c                     z_c = sqrt(12 k_c)(0.5 − r_c) exactly as in C8
S       = (1 / 34) · Σ over the 34 targeting constructs of z_c^rep
```

**Why (a) and not (b), argued from the diagnostic and not from its effect on the answer.**

1. **The pedestal is real but it is not a clean constant.** C8's KA-4 gives its *mean*
   (+0.0717 ± 0.0022 at 5,000 draws), but C8's P3 gives its *spread*: the wrong-recipient level of
   the partner construct reads +0.041 with a construct-level 95 % interval **[−0.295, +0.369]** —
   a per-construct spread four to five times the mean it is scattered around. A single scalar
   (variant b) removes the pedestal from the point estimate and leaves every construct's own null
   level inside its contribution.
2. **The interval is the deliverable here, and (b) does not repair it.** With the construct as the
   unit, the reported uncertainty is a construct-level bootstrap over 34 numbers. Under (b) those
   34 numbers are still uncentred; under (a) each is referenced to its own null. (a) repairs the
   estimate and the interval together, and it **reduces to (b) if the per-construct levels turn out
   homogeneous** — so (a) is the weaker assumption, not the stronger one.
3. **Better covariate balance does not remove the pedestal.** C8's P2 measured the residual bias at
   **+0.081** after repairing the KA-6 axis. So the bias is not a covariate-imbalance artefact that
   the §6.2 repair would fix as a side effect; it needs its own repair.
4. **(a) is C8's own strongest structural check, generalised.** The partner-construct test is
   exactly `z_c` minus the same statistic in **one** wrong recipient (the partner). (a) is that
   contrast estimated against **every** wrong recipient in the split.

**Direction, stated before running: this repair works AGAINST the effect.** It subtracts a positive
pedestal of ≈ +0.072 from an estimate of +0.786, i.e. it is expected to move the point estimate
**down by roughly 9 %**, to ≈ +0.71. It makes passing harder on the point estimate. Its effect on
the *null's spread* is not signed in advance: the centred statistic is a paired difference, so its
null may be tighter than C8's (shared variance across recipients cancels) or wider (two estimated
terms); **both the repaired point estimate and the repaired null's sd are reported, and the p is
whatever it is.**

**What (a) costs, stated before running.** Under (a) the repaired KA-4 row (§9, KA-4R) is zero in
expectation **by construction**: the statistic subtracts precisely the level the row measures. Its
evidential content is therefore reduced to (i) verifying the correction is correctly implemented
and (ii) confirming by Monte Carlo that the interval brackets zero. That is an unavoidable property
of any bias-subtraction repair — variant (b) has it too — and it is the reason the structurally
different checks (KA-5's null calibration, the real-spacer resampling null, the partner test) are
kept and reported.

**Two conventions, fixed now, with their directions.**
- If a candidate gene `G ∈ C_c` is also in wrong recipient `d`'s own hit set at any L, `d` is
  **kept** and `G` is scored in `d`'s cells anyway. If that gene is genuinely repressed by `d`, the
  wrong-recipient level rises, `b_c` rises, and the correction gets **larger** — conservative.
- Pool members are never any construct's `L*` candidate (§6.2's global exclusion), but may be in
  some recipient's non-`L*` hit set. The expected overlap is ≈ 0.4 genes of 200 and the induced
  shift in a percentile is ≈ 0.002; it is reported, not surgically removed, so that **pools are
  recipient-independent** and the observed statistic and every null draw traverse one code path.

### 6.2 KA-6 / matching — promoter-window count added to the matched axes

**The repair.** The hard strata become **10 basal-expression deciles × 5 `pi_off_L` quintiles × 3
promoter-window classes (n_windows = 1, 2, ≥ 3) = 150 cells** (C8: 50 cells). Soft matching inside
the cell is unchanged: the **M = 200** genes nearest to the candidate in the z-scored space
`[gc, log1p(tss_density_10kb), map_mult, lowcomp, mean_offset/100, log1p(n_windows),
log1p(n_sites)]`. The widening rule is unchanged in form: if a cell holds fewer than **20** eligible
genes the pool widens to the same (`pi_off` quintile, window class) ∩ adjacent basal deciles, and
each widening is recorded per candidate; if it holds fewer than **50**, the whole cell is used
without the nearest-200 trim.

**Eligibility (C8's §5.2, unchanged) plus one new global exclusion.** A gene is eligible for
construct `c` iff `pi_off_{L*(c)}(g) > 0`; it is not one of the 300 genes any 2025 construct intends
to target; it is not `c`'s own target; it is not in `c`'s hit set at any L; and it is measured.
**New:** it is not the `L*` candidate of **any** of the 65 constructs (**134 genes**, measured
pre-hash, ≈ 0.9 % of the ≈ 14,132-gene eligible universe at L = 12).

**Why the new exclusion, and its direction — the one repair element that could favour the effect.**
It is required by §6.1: with it, a pool contains no gene that is scored as a candidate anywhere, so
**one pool serves every recipient** and the observed statistic, the wrong-recipient scores and every
null draw use identical reference sets. Its direction is **not** conservative: if those 134 genes
are themselves repressed, removing them lifts the pools slightly and can move the estimate **up**.
Expected magnitude: ≈ 1.9 of 200 pool members. **It is therefore reported as a declared sensitivity
(§10, Sens-3): `S` recomputed with the global exclusion switched off.** The primary keeps the
exclusion, decided here, before any value exists.

**Reported for the repaired pools (§9, KA-6R):** the **full nine-axis SMD table** — `pi_off_L*`,
`basal`, `gc`, `log1p_tss_density`, `map_mult`, `lowcomp`, `mean_offset/100`, `log1p_n_windows`,
`log1p_n_sites` — every axis must meet `|SMD| ≤ 0.250`. Also reported: mean and **minimum** pool
size, the number of pools that needed widening, and the number of candidates whose stratum cell
held < 20 eligible genes. **If adding the axis makes the matched pools too small for some
construct, that construct is NOT dropped**: it is reported by name with its attainable pool size,
and the resulting attainable resolution is stated.

### 6.3 What is NOT changed

The candidate rule, the response formula, the population, the unit, the direction convention, the
null's *family* (matched candidate-set resampling), `B = 100,000`, `M = 200`, the widening
thresholds, α, the number of verdict cells, the leave-one-out condition, and the abort clause.

### 6.4 No post-hoc choice between variants

The repair variants are fixed by this document. **No variant may be selected, swapped or blended
after any repaired primary value has been seen.** Variant (b) — scalar subtraction of the measured
bias with a propagated interval — is computed and reported as a declared secondary (§10, Sens-1)
for comparison only, and carries no verdict.

## 7. THE PRIMARY STATISTIC AND ITS NULL

**Primary statistic.** `S = (1/34) Σ_c (z_c − b_c)` over the 34 targeting constructs, on the
repaired pools of §6.2. Positive = candidates more repressed than their matched non-candidates,
over and above the level the same genes read in the other constructs of the split.

**The joint null (repaired).** `B = 100,000` joint draws. In draw `b`, **every** construct's whole
candidate set is redrawn — each candidate `G` replaced by a gene drawn uniformly from `G`'s own
matched pool `P_G`, without replacement across the `k_c` slots — and **both** terms are recomputed
from the drawn genes through the identical code path:

```
z'_c = sqrt(12 k_c)(0.5 − mean_G q_c(m_G))        b'_c = mean over d in retained(c) of sqrt(12 k_c)(0.5 − mean_G q_d(m_G))
S^(b) = (1/34) Σ_c ( z'_c − b'_c )
p_joint = (1 + #{ b : S^(b) >= S_obs }) / (B + 1)              one-sided, repression direction
```

`retained(c)` is fixed at its observed value and is identical for the observed statistic and every
null draw.

**Interval.** Construct-level bootstrap: **20,000** resamples of the 34 values `z_c^rep` with
replacement, percentile 95 % CI. It is labelled *n = 34 targeting constructs*, it measures
construct-to-construct variability, **not** the null, and the verdict rests on `p_joint`.

**Interpretable effect size (reported with the primary, no separate verdict).** The same centring on
the log1p-CP50k scale:
`Δ_c = mean_{G∈C_c}[ e_c(G) − median(e_c over P_G) ]`,
`Δ^wr_c = mean over d in retained(c) of mean_{G∈C_c}[ e_d(G) − median(e_d over P_G) ]`,
`Δ = (1/34) Σ_c (Δ_c − Δ^wr_c)`, with its own construct-level bootstrap CI.

**Reported diagnostics, no verdict.** The distribution of `b_c` over the 34 constructs (mean, sd,
range) — the heterogeneity that §6.1's argument turns on; the number of wrong recipients per
construct; the repaired null's mean and sd against C8's (+0.00005, 0.17173).

## 8. PASS RULE, VERDICT CELLS, ABORT CLAUSE — fixed before computing anything

- **Number of verdict cells: ONE.** No multiplicity correction is required *within* C9; §1.1's
  across-claim multiplicity is a reporting obligation, not a threshold change. Every quantity in
  §10 is reported without a verdict and cannot change this one.
- **α = 0.05, one-sided** (candidates more repressed).
- **C9 PASSES iff both hold:**
  1. `p_joint ≤ 0.05`; **and**
  2. **leave-one-construct-out robustness:** recomputing `S` and `p_joint` with each of the 34
     constructs deleted in turn (denominator 33), the **maximum** of the 34 leave-one-out p-values
     is `≤ 0.05`.
- **C9 FAILS** if condition 1 or condition 2 does not hold.
- **C9 is `UNRESOLVABLE-AT-THIS-RESOLUTION` iff any of the nine known-answer rows of §9 fails its
  stated rule.** This clause is C8's clause verbatim. It is **not** weakened because it has already
  fired once. A FAIL will be reported as a bounded estimate with its interval, never as evidence of
  absence.
- **Attainable resolution.** The joint null's smallest attainable p is `1/(B+1) = 1.0e-5`, three
  orders of magnitude below α, so a FAIL here would be informative.
- **Sign convention frozen:** `S > 0` = candidates repressed relative to matched non-candidates,
  net of the wrong-recipient level.

## 9. THE NINE KNOWN-ANSWER ROWS — run before the verdict; a failure aborts the conclusions

C8's seven rows are re-run unchanged; the two that failed are additionally recomputed **on the
repaired machinery**, where they must meet their **original** rules. Rows KA-4 and KA-6 on the
*unchanged* machinery govern by **reproduction of C8's published values** — they are regression
checks that C9 inherited C8's machinery correctly, not a re-litigation of the thresholds C8 aborted
on.

| id | row | pass rule |
|---|---|---|
| KA-1 | 305's 23 published rows (pooled effect, both half effects, control CP50k) recomputed from `sums305_<split>.npz` | max \|diff\| ≤ 1e-9 |
| KA-2(a) | 301's selection reproduced from its saved hits with 301's exact filters | failing **4/13**, working **4/21**, NTC **4/31**; **8** L = 12 targeting constructs; 83 `L*` pairs over 34 constructs |
| KA-2(b) | **the inherited pool is the CORRECTED one**: the Excel escape `="…"` absent from all 209,070 pool gene labels; no spacer's own target symbol in its own hit set at any L; the pool's and the 68 real spacers' L = 12 hit rates reproduce C8's saved values; the *unfiltered* (bug-equivalent) rate recomputed from `lmer_genes` for contrast | 0 escaped labels; **0** own-target violations at all four L; pool rate = 0.160879 and real rate = 0.132353 to ≤ 1e-6; bug-equivalent rate materially higher |
| KA-3 | the 31 NTC constructs' per-half integer sums across the three splits | identical names, max \|sum diff\| 0, max \|cell-count diff\| 0, 38,176 control cells |
| KA-4 | **unchanged machinery** — the no-signal row on C8's pools and C8's uncentred statistic, 5,000 within-split derangements | **reproduction**: mean within **±0.010** of C8's 5,000-draw **+0.0717** |
| KA-5 | null calibration — 200 draws of the **repaired** primary null, each p-valued against the rest; additionally reproduced on C8's null | KS uniformity p ≥ 0.05 **and** fraction ≤ 0.05 inside [0.02, 0.10], on **both** nulls |
| KA-6 | **unchanged machinery** — the nine-axis SMD table on C8's pools | **reproduction**: every axis within **0.002** of C8's published SMDs, including `log1p_n_windows` **+0.2905** |
| KA-7 | 293's on-target `kd_lfc` for all 34 constructs recomputed from 305's sums | max \|diff\| ≤ 1e-9 vs `seed298.csv` |
| **KA-4R** | **repaired** — the no-signal row recomputed on the repaired statistic `z^rep`, 5,000 within-split derangements | **C8's original rule**: \|mean\| ≤ 0.050 **and** the 2.5–97.5 % interval brackets 0 |
| **KA-6R** | **repaired** — the nine-axis SMD table on the repaired (150-cell) pools | **C8's original rule**: \|SMD\| ≤ 0.250 on **every** one of the nine axes |

(The table lists KA-2 as two clauses of one row, so the count of rows that govern is **nine**:
KA-1, KA-2, KA-3, KA-4, KA-5, KA-6, KA-7, KA-4R, KA-6R.)

**KA-0, reported not a verdict row:** the split composition of the 34 targeting constructs is
re-derived and printed. C8's `RESULTS.md` §1 states "21 train / 7 validation / 6 test"; C8's own
`constructs307.csv` reads train 16 / test 10 / validation 8. Whichever is right, C9's within-split
centring depends on it, so it is measured and reported, and any discrepancy with C8's write-up is
recorded as an audit finding.

## 10. DECLARED SECONDARIES — reported, no verdict, cannot change §8

All are computed on the repaired machinery unless stated; C8's corresponding values are quoted
beside each.

1. **The real-spacer resampling null** (C8's S5, kept): a pseudo-construct is a random **pair of
   stringency-0 hCRISPRi-v2 protospacers** rescanned through 301's exact rule, kept by rejection
   only if it yields `L* = L*(c)` and exactly `k_c` candidates after `c`'s exclusions; scoring is
   matching-free (percentile in the whole eligible universe). Reported **uncentred** (reproducing
   C8's `W_obs` = +0.74470, p = 0.00149) **and centred** by the same wrong-recipient construction.
2. **The wrong-recipient permutation** (C8's P1, kept): under the repaired statistic this is the
   same object as KA-4R; its draw distribution also yields a permutation p for the observed
   repaired `S`, which is reported. **Stated plainly:** under repair (a) this is no longer an
   independent check of the estimate, because the repair subtracts exactly the level it measures;
   C8's p = 2.0e-4 stands as the reading on C8's machinery.
3. **Partner-construct specificity** (C8's P3, kept): each construct's candidate genes re-scored in
   the cells of its **partner construct — same intended target, different protospacers** — on the
   **repaired pools**, **uncentred**, so it remains the same contrast C8 reported; own − partner
   with a paired Wilcoxon and a construct bootstrap.
4. **The failing-vs-working contrast** (301's original T7 question, kept): `z^rep` for the 13
   failing vs the 21 working constructs, difference with a construct-level bootstrap CI and a
   Mann–Whitney p.
5. **The 31 non-targeting constructs** run identically through the repaired pipeline, **labelled
   candidate-bearing observations selected by the same rule, not negative controls**.
6. **By seed level:** `S` on the L\* = 12 / 11 / 10 / 9 subsets.
7. **Sens-1 — repair variant (b):** `S_C8 − b̄` with an interval propagated from the
   wrong-recipient mean's standard error, for comparison with the chosen variant (a).
8. **Sens-2 — centring over all 33 other targeting constructs** instead of within-split only.
9. **Sens-3 — the global `L*`-candidate exclusion switched off** (§6.2's direction check).
10. **Sens-4 — the two repairs separately:** the KA-6 repair alone (uncentred, 150-cell pools) and
    the KA-4 repair alone (centred, C8's 50-cell pools), so the contribution of each is visible.
11. **Sign count, measurement floor and per-construct table**, as in C8.

## 11. Seeds, determinism, resources, and what is not done

`SEED309 = 20260907`; `numpy.random.default_rng`. `B = 100,000` joint null draws, 20,000 construct
bootstrap resamples, 5,000 derangements for KA-4 and KA-4R, 200 calibration draws — all fixed here.
Everything is read from 305's saved `.npz` sums and 301's/307's saved sequence products; **no h5ad
is opened, nothing is re-streamed, no download, no reserve access, no new seed search.** Peak RSS
budget **< 4 GiB**. `uv run python` throughout.

**Nothing is written outside `experiments/309-c9-calibration-repair/`** except the C9 entry in
`docs/CONFIRMATION-LEDGER.md`. **No `vcc`. No `cell-eval2`. No build. No submission package,
nothing under `submissions/` or `data/`. No model. No production change.**

## 12. Deviations

Every deviation from this document will be listed in `RESULTS.md` §Deviations, each marked as fixed
**before** or **after** any repaired outcome was seen.
