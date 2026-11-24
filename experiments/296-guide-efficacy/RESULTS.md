# 296 — RESULTS: is F14's "construct term" guide EFFICACY?

**GUIDE EFFICACY EXPLAINS THE CONSTRUCT TERM: PARTLY** — K562-GWPS both-≥cut ρ_guide **+0.340
[+0.113, +0.548]** vs weaker-<cut **+0.001 [−0.082, +0.105]** (n = 50, library-relative cut; at 293's
absolute cut +0.610 vs +0.073 on only 5 pairs); Nadig **Jurkat +0.458 [+0.245, +1.180]** vs
**+0.096 [−0.175, +0.210]** (PASSES at all three cuts and both rungs); Nadig **HepG2 +0.645
[+0.185, +1.016]** vs **+0.307 [−0.023, +0.498]** (difference +0.338 [−0.221, +0.808], not separated);
**Broad-16 (Cas9) +0.951 [+0.710, +1.184]** vs **+0.763 [+0.621, +0.886]** — no split, P3 confirmed;
**CD4 +0.018 [−0.006, +0.057]** vs **+0.014 [−0.017, +0.056]** — flat at every condition and budget, a
well-powered null. **FAILING CONSTRUCTS REPRODUCIBLE: YES** (own split-half above zero in all five
libraries: K562 +0.0069 [+0.0062, +0.0077], Jurkat +0.0348 [+0.0270, +0.0401], HepG2 +0.1375
[+0.1275, +0.1566], Broad-16 +0.0907 [+0.0866, +0.0945], CD4 +0.1042 [+0.0963, +0.1067] — but *less*
reproducible than working constructs, by +0.007 / +0.099 / +0.124 / +0.019 / +0.007; **shared vector
signature: NO** — a failing construct is no closer to the library mean than a working one (K562 +0.1029 vs
+0.1044) and failing constructs of different genes agree no better than working ones in un-centred space
(K562 +0.0099 vs +0.0099; Jurkat +0.0112 vs +0.0130, i.e. *worse*)). **FAILURE RATES (|lfc| < 1.0 / < 0.5,
own-gene column, full depth): Arc 2025 panel 18.7 % / 4.3 %; K562-GWPS 2026-panel targets 70.1 % / 24.0 %,
non-panel 75.0 % / 30.9 %; Nadig Jurkat 53.4 % / 18.8 %, HepG2 55.5 % / 19.0 %; Broad-16 (Cas9) 97.7 % /
82.0 %; CD4 GWCD4i 63.9 % / 29.0 % (the authors' own call: 26.7 % of 73,765 guide × condition rows are not
a significant knockdown)** — the absolute cut is **not** comparable across libraries, because Arc's panel
was selected for high basal expression of the target gene and the `log1p` scale compresses low-expressed
genes; the library-relative cut and the authors' own CD4 call are the transferable numbers.

Preregistration `PREREGISTRATION.md`, sha256 `ac84c52bd64b1ffd104363e7b3dd68d0c239c49f09bf7a0eed58ba15a1022ffa`
(`work/PREREGISTRATION.sha256`), hashed **2026-09-06 16:41:42Z**, before any stratified number. Nothing was
written outside `experiments/296-guide-efficacy/`; the write audit at the end lists only the two sibling
agents' own trees.

---

## 1. Known-answer rows (`work/ka296.{py,json,log}`) — 5/5 PASS

| id | row | required | got |
|---|---|---|---|
| **KA-1** | 285's K562 guide row at n = 25 (748 pairs) recomputed from `k562_guide.npz` | max \|Δ\| = 0.0e+00 | **0.0e+00** — (b) +0.000567, (a) +0.006320, ρ +0.089643, n 748; `pds` 0.5093 / 0.6330 |
| **KA-2** | 285 §6b's three knockdown terciles re-derived from `k562_kd.npz` by this code | max \|Δ\| = 0.0e+00 across n / weaker-KD / guide / cells / ρ | **0.0e+00** — T1 180 / +0.015 / +0.0009 / +0.0072 / 0.124; T2 179 / −0.119 / +0.0003 / +0.0064 / 0.041; T3 180 / −0.300 / +0.0021 / +0.0089 / 0.236 |
| **KA-3** | 293's ALL / both-≥1.0 / weaker-<1.0 rows re-derived from `construct293.csv` by this code | every cell of `strat293.json` to ≤ 1e-12 | **PASS** — ρ +0.074662 / +0.687417 / −0.002632, (a) +0.312262 / +0.149456 / +0.316801, (b) +0.023314 / +0.102738 / −0.000834, k 16 / 5 / 11 |
| **KA-4** | 268's halves satisfy `A + B == full` and `v1 + v2 ≤ full` on 12 lines | 0 / ≤ 0 | **0.0 / 0.0** |
| **KA-5a** | the CD4 bank's 590 targeting guides all appear in the authors' `guide_kd_efficiency` table | 0 missing | **0 missing** |
| **KA-5b** | my bank-computed on-target \|lfc\| vs the authors' own `−t_statistic`, 1,753 (condition, guide) rows | Spearman > +0.3 | **+0.640** (p 2.1e-202) |

*Two facts found on the way, both confirmatory.* (i) 293's printed §5a table shows (a) **+0.150** for the
both-≥1.0 stratum where the stored value is **0.1494559** — a display rounding, not a discrepancy.
(ii) `rates296.py` recomputes 293's 2025-panel failure rate from 293's own saved `kd_all` and reproduces
**18.7 % / 4.3 %** and median |lfc| **1.485** exactly.

**Also reproduced, not preregistered as a KA** (different construct caps and control caps, so exact equality
was not required): 285's published F14 rows. Nadig Jurkat n = 25 — (a) +0.0378 vs 285's +0.0416, (b) +0.0080
vs +0.0070, ρ_guide **+0.210** vs 285's **0.17**, on the same 124 pairs. Nadig HepG2 n = 25 — (a) +0.1608 vs
+0.1622, (b) +0.0640 vs +0.0659, ρ_guide **+0.398** vs **0.41**, on the same 70 pairs. Broad-16 at
T = 2×10⁵ — ρ_guide **+0.793** vs 285's **0.81**. CD4 at T = 1–2×10⁵ — ρ_guide **0.016–0.056** vs 285's
**0.02–0.06**. **The one row that does not reproduce stably is K562-GWPS**: ρ_guide **+0.027 [−0.076, +0.144]**
at n = 25 here against 285's **+0.090 [−0.038, +0.232]**, and +0.079 vs +0.062 at n = 50. The reason is
structural, not a bug: (b) is +0.0002…+0.0008 on a denominator of +0.0063…+0.0108, so **ρ_guide's ALL row in
K562-GWPS is a ratio of a number indistinguishable from zero to a small number** and moves by 3× between two
legitimate draws of the same design. Every interval overlaps. **The ALL row for K562 carries no information;
the stratified contrast does, and that is what §3 reports.**

---

## 2. What was measured

Five libraries, 13 (library × budget) cells, one frozen contrast set per cell (`PREREGISTRATION` §2.2):
**(a)** the larger construct's own two halves (the matched cell-sampling reference); **(a2)** the smaller
construct's own halves; **(b)** the two constructs against each other; **(x)** the cross-**gene** null
(side A against a *different* gene's second construct); ρ_guide = med(b)/med(a); ρ_dis = med(b)/√(med(a)·med(a2));
the scorer's `pds` for each. 8 draws, seed 20260906, paired target bootstraps 20,000 draws.

| library | modality | construct unit | pairs | budget | cells extracted / peak RSS |
|---|---|---|---:|---|---|
| K562-GWPS (Replogle 2022) | CRISPRi | `guide_id` inside one `perturbation` | 748 / 653 | n = 25 / 50 cells/side | 280,071 cells, **9.82 GiB**, 219 s |
| Nadig Jurkat | CRISPRi | `sgID_AB` inside one `gene` | 124 / 60 | n = 25 / 50 | 90,711 cells, 6.02 GiB, 13 s |
| Nadig HepG2 | CRISPRi | as above | 70 | n = 25 | 58,958 cells, 6.12 GiB, 8 s |
| Broad-16, 12 lines | **Cas9 KO** | 268's `v1` / `v2` sums | 803 / 404 (line, target) | T = 2×10⁵ / 5×10⁵ UMIs | sums only, 2.56 GiB, 122 s |
| CD4 GWCD4i | CRISPRi | guide slot 0 / 1 in one (donor, condition) | 158–231 per condition | T = 1×10⁵ / 2×10⁵ UMIs | bank only, 2.93 GiB, 128 s |

Every pass stayed under the 10 GiB cap (the K562 pass preallocates its output instead of concatenating
blocks, which is what took 285's equivalent to 11.84 GiB). The stratifier — each construct's own gene
column at full depth — is **excluded from every scored statistic**, so nothing here is circular.

### 2a. A structural fact that qualifies F14 before any stratification

**In K562-GWPS and Nadig, the two "constructs" of a gene are almost always two different PROMOTERS**, not
two guides against one transcript: the construct ids carry `-P1` / `-P2` suffixes
(`A1BG_+_58858964.23-P1|A1BG_-_58858788.23-P1`), and **713 of 748** K562 pairs, **103 of 120** Jurkat pairs
and **52 of 69** HepG2 pairs are P1-vs-P2. Arc's own 2025 library has the same design (293's `CAST_P1` /
`CAST_P2`). So F14's "guide construct" unit is, in every CRISPRi library on disk, a **transcription-start-site
unit**. §7 reports the same-promoter comparison; it is underpowered in K562 (3 pairs) and null in Jurkat.

---

## 3. The stratified table (`work/strat296.{py,json,log}`, `work/efficacy296.csv`)

Cuts, all frozen before the measurement: **C1** = 293's absolute |lfc| ≥ 1.0; **C2** = the library's own
median single-construct |lfc|; **C3** = knockdown fraction ≥ 0.5 (expression-free). Strata: both constructs
clear the cut / the weaker one does not. `k` = pairs in the stratum.

### 3a. Headline cells

| library (rung) | cut | k both / weak | ρ_guide **both** | ρ_guide **weaker fails** | ρ(both) − ρ(weak), stratum bootstrap | P1 |
|---|---|---:|---|---|---|---|
| **K562-GWPS** n = 25 | C1 (1.0) | 6 / 533 | +0.444 [−0.970, +3.068] | +0.041 [−0.078, +0.191] | +0.404 [−1.004, +2.909] | **UNDERPOWERED** |
| K562-GWPS n = 25 | C2 (0.320) | 81 / 458 | +0.123 [−0.059, +0.427] | +0.011 [−0.119, +0.198] | +0.112 [−0.138, +0.450] | fail |
| K562-GWPS n = 25 | C3 (50 % KD) | 81 / 458 | +0.288 [−0.035, +0.525] | −0.003 [−0.127, +0.165] | +0.291 [−0.074, +0.549] | fail |
| **K562-GWPS n = 50** | C1 (1.0) | 5 / 464 | +0.610 [−0.576, +4.541] | +0.073 [−0.046, +0.150] | +0.537 [−0.651, +4.472] | **UNDERPOWERED** |
| **K562-GWPS n = 50** | **C2 (0.324)** | **70 / 399** | **+0.340 [+0.113, +0.548]** | **+0.001 [−0.082, +0.105]** | **+0.339 [+0.089, +0.568]** | **PASS** |
| K562-GWPS n = 50 | C3 | 73 / 396 | +0.203 [+0.049, +0.388] | +0.018 [−0.080, +0.113] | +0.185 [−0.005, +0.391] | fail (marginal) |
| **Nadig-Jurkat n = 25** | **C1** | **8 / 94** | **+0.776 [+0.174, +1.544]** | **+0.147 [−0.036, +0.280]** | **+0.629 [+0.016, +1.340]** | **PASS** |
| **Nadig-Jurkat n = 25** | **C2 (0.548)** | **21 / 81** | **+0.458 [+0.245, +1.180]** | **+0.096 [−0.175, +0.210]** | **+0.362 [+0.136, +1.202]** | **PASS** |
| **Nadig-Jurkat n = 25** | **C3** | **31 / 71** | **+0.302 [+0.137, +0.765]** | **+0.095 [−0.227, +0.217]** | **+0.208 [+0.014, +0.826]** | **PASS** |
| Nadig-Jurkat n = 50 | C1 | 5 / 43 | +0.627 [−0.248, +1.119] | +0.026 [−0.037, +0.249] | +0.601 [−0.290, +1.093] | UNDERPOWERED |
| **Nadig-Jurkat n = 50** | **C2** | **12 / 36** | **+0.686 [+0.258, +0.994]** | **−0.001 [−0.101, +0.115]** | **+0.687 [+0.251, +1.029]** | **PASS** |
| **Nadig-Jurkat n = 50** | **C3** | **22 / 26** | **+0.283 [+0.029, +0.732]** | **−0.024 [−0.322, +0.159]** | **+0.307 [+0.001, +0.821]** | **PASS** |
| Nadig-HepG2 n = 25 | C1 | 4 / 57 | +0.910 [−0.992, +1.597] | +0.370 [+0.070, +0.496] | +0.541 [−1.356, +1.375] | UNDERPOWERED |
| Nadig-HepG2 n = 25 | C2 (0.568) | 11 / 50 | +0.645 [+0.185, +1.016] | +0.307 [−0.023, +0.498] | +0.338 [−0.221, +0.808] | fail |
| Nadig-HepG2 n = 25 | C3 | 18 / 43 | +0.688 [+0.320, +0.880] | +0.144 [−0.222, +0.492] | +0.544 [−0.018, +0.920] | fail (marginal) |
| **Broad-16** T = 2×10⁵ | C1 | 6 / 785 | +0.442 [−1.046, +12.6] | +0.802 [+0.706, +0.921] | −0.360 [−1.871, +11.8] | UNDERPOWERED |
| **Broad-16** T = 2×10⁵ | **C2 (0.242)** | **280 / 511** | **+0.951 [+0.710, +1.184]** | **+0.763 [+0.621, +0.886]** | **+0.187 [−0.093, +0.452]** | **fail → P3 PASS** |
| Broad-16 T = 5×10⁵ | C2 (0.250) | 152 / 246 | +0.931 [+0.748, +1.306] | +0.829 [+0.686, +1.045] | +0.102 [−0.181, +0.495] | fail → P3 PASS |
| **CD4** T = 1×10⁵ Rest | C1 | 53 / 156 | +0.019 [−0.016, +0.054] | +0.014 [−0.014, +0.054] | +0.004 [−0.048, +0.050] | **fail (powered)** |
| CD4 T = 1×10⁵ Stim8hr | C1 | 57 / 174 | +0.062 [+0.016, +0.100] | +0.040 [+0.005, +0.056] | +0.021 [−0.029, +0.068] | fail |
| CD4 T = 1×10⁵ Stim48hr | C1 | 49 / 181 | +0.030 [−0.017, +0.098] | +0.046 [+0.025, +0.063] | −0.016 [−0.069, +0.054] | fail |
| CD4 T = 2×10⁵ (3 conds) | C1/C2/C3 | 41–163 / 23–154 | +0.017…+0.053 | +0.028…+0.056 | −0.016…+0.026, every interval through 0 | fail |

The complete 13 × 3 grid, with (a), (a2), (b), (x), ρ_dis and `pds` for every stratum, is in
`work/strat296.json` and printed in `work/strat296.log`. Per-(library, target, construct) rows —
8,144 of them — are in **`work/efficacy296.csv`**.

### 3b. The ALL rows, and what the cross-gene null does to F14

| library (rung) | ρ_guide ALL | (b) − (x), paired | `pds` (b) / (x) | ρ null-corrected = (med b − med x)/med a |
|---|---|---|---|---:|
| K562-GWPS n = 25 | +0.027 [−0.076, +0.144] | +0.0004 [−0.0005, +0.0013] | 0.506 / 0.496 | +0.056 |
| K562-GWPS n = 50 | +0.079 [−0.026, +0.150] | +0.0004 [−0.0011, +0.0016] | 0.513 / 0.509 | +0.038 |
| Nadig-Jurkat n = 25 | +0.210 [+0.108, +0.331] | +0.0030 [−0.0025, +0.0089] | 0.552 / 0.475 | +0.079 |
| Nadig-HepG2 n = 25 | +0.398 [+0.211, +0.511] | **+0.0049 [−0.0282, +0.0330]** | 0.434 / 0.464 | **+0.030** |
| **Broad-16** T = 2×10⁵ | +0.793 [+0.706, +0.919] | **+0.0092 [+0.0075, +0.0109]** | **0.696 / 0.500** | **+0.792** |
| Broad-16 T = 5×10⁵ | +0.878 [+0.748, +1.029] | +0.0186 [+0.0153, +0.0240] | 0.764 / 0.489 | +0.910 |
| CD4 T = 2×10⁵ Stim8hr | +0.056 [+0.030, +0.076] | +0.0106 [+0.0040, +0.0152] | 0.607 / 0.509 | +0.053 |

**This is new, and it corrects F14's HepG2 number.** In the four CRISPRi libraries the construct term at
the median residual cosine is **not separated from the cross-gene null**: every (b) − (x) interval crosses
zero except CD4's. HepG2's ρ_guide **0.41** — F14's highest CRISPRi value — comes almost entirely from a
response shared across *different* genes (its targets are DepMap-essential genes and share a growth-arrest
programme that the 507-target panel centring does not remove): (x) reads **+0.0591** against (b) **+0.0640**,
so the null-corrected ratio is **+0.030**, not +0.40. The same correction leaves Jurkat at +0.079 and K562
at +0.038–0.056. **Only Broad-16 (Cas9) has a construct term that survives its own cross-gene null**
(+0.792 corrected vs +0.793 raw), and only Broad-16's `pds` for (b), 0.696, sits far above the null's 0.500.
This is exactly the pattern 293 found in Arc's 2025 library, where every (b) − (x) interval also crossed
zero while `pds` separated them. It does **not** change the stratified result in §3a — the strata share one
(x) — but it means **F14's cross-library ordering of ρ_guide (0.09 / 0.17 / 0.41 / 0.02–0.06 / 0.76–0.81)
partly ranks how much shared biology each library's target set carries, not how much construct-specific
biology there is.** Null-corrected, the ordering is K562 ≈ HepG2 ≈ CD4 ≈ 0.03–0.06, Jurkat 0.08, Broad-16 0.79.

---

## 4. The five frozen predictions

| | prediction | rule | **outcome** |
|---|---|---|---|
| **P1** | in every CRISPRi library ρ_guide(both ≥ cut) > ρ_guide(weaker < cut), interval above 0 | 2.5th stratum-bootstrap percentile > 0 at C1; < 8 pairs in either stratum ⇒ UNDERPOWERED | **PARTLY.** At **C1** only Jurkat n = 25 (8/94) and CD4 (41–57 per stratum) are testable: **Jurkat PASSES**, **CD4 FAILS**; K562 (5–6 pairs) and HepG2 (4 pairs) are UNDERPOWERED, and Broad-16 (6 pairs) too. At **C2**, where every library is powered: **K562 n = 50 PASSES**, **Jurkat PASSES at both rungs**, **HepG2 FAILS** (+0.338 [−0.221, +0.808]; the point estimate is large and positive), **CD4 FAILS flat**. So **2 of 4** CRISPRi libraries at the powered cut |
| **P2** | K562-GWPS both-≥1.0 ρ_guide ≥ 0.40 | point estimate at n = 25, cut C1 | **PASS on the letter, not on the power**: +0.444 (n = 25) and +0.610 (n = 50) — but on **6 and 5 pairs**, intervals [−0.970, +3.068] and [−0.576, +4.541]. The powered version of the same statement is C2's +0.340 [+0.113, +0.548] at n = 50 |
| **P3** | Broad-16 (Cas9, validated guides) shows no such split | P1's rule fails on Broad-16 at C1 | **CONFIRMED** — and more strongly than the rule asked: at the powered C2 cut, ρ_guide is +0.951 in the both-clear stratum and +0.763 where the weaker construct does not clear (difference +0.187 [−0.093, +0.452]); at C3, +0.834 vs +0.804. Broad-16's construct term is high *everywhere*, including where the target's mRNA does not drop at all |
| **P4** | the weaker construct's knockdown predicts ρ_dis at Spearman > +0.3 in each CRISPRi library | Spearman(kd_weak as strength, per-target ρ_dis) > +0.30 | **1 of 4.** Jurkat **+0.476** (p 6.2e-05, n = 80) and **+0.424** (n = 30) → **PASS**; K562 **+0.055 / +0.056** → FAIL; HepG2 **+0.025** (n = 26) → FAIL; CD4 **+0.003…+0.135** → FAIL. Broad-16 (not part of P4) +0.039 / +0.086 |
| **P5** | failing constructs are (i) reproducible and (ii) share a vector signature | (i) split-half interval above 0; (ii) failing × failing cross-gene agreement exceeds working × working | **(i) CONFIRMED in all five libraries. (ii) REFUTED** — see §5 |

**The single most informative number in P4** is not the coefficient but the contrast inside it. In Jurkat,
knockdown predicts the **construct** agreement (b) at **+0.446** while predicting the **cell-split**
agreement (a) at **−0.022**: the effect is specific to the construct unit, which is what 293 claimed. In
K562 it predicts both about equally (+0.062 vs +0.059 at n = 25; +0.141 vs +0.085 at n = 50) — **285 §6b's
finding reproduced exactly**, i.e. in K562 the tercile is selecting response magnitude, not construct
cleanliness. In CD4 it predicts neither. So the mechanism 293 found in Arc's library is present in one
public CRISPRi library (Jurkat), present but confounded with magnitude in another (K562), and absent in
the other two.

---

## 5. The failing-construct phenomenon (`PREREGISTRATION` §3)

"Failing" = own-gene |lfc| < 0.5 at full depth; "working" = ≥ 1.0. Counts are constructs, not targets.

| library (rung) | failing / working | own split-half **failing** | **working** | working − failing | cos to its partner, failing | cos to the LIBRARY MEAN (un-centred) failing / working |
|---|---:|---|---|---|---|---|
| K562-GWPS n = 25 | 716 / 131 | **+0.0069 [+0.0062, +0.0077]** | +0.0136 [+0.0095, +0.0172] | +0.0067 [+0.0023, +0.0102] | +0.0004 [−0.0003, +0.0013] | +0.1029 / +0.1044 (diff −0.0015 [−0.0123, +0.0051]) |
| K562-GWPS n = 50 | 620 / 113 | **+0.0114 [+0.0103, +0.0132]** | +0.0248 [+0.0197, +0.0313] | +0.0134 [+0.0073, +0.0196] | +0.0008 [−0.0003, +0.0016] | +0.0843 / +0.0863 |
| Nadig-Jurkat n = 25 | 101 / 56 | **+0.0348 [+0.0270, +0.0401]** | +0.1339 [+0.0660, +0.1927] | +0.0991 [+0.0351, +0.1576] | +0.0047 [−0.0052, +0.0079] | +0.1087 / +0.1511 (**−0.0424 [−0.0589, −0.0176]**) |
| Nadig-HepG2 n = 25 | 54 / 34 | **+0.1375 [+0.1275, +0.1566]** | +0.2616 [+0.1556, +0.3805] | +0.1241 [+0.0117, +0.2437] | +0.0581 [+0.0123, +0.0778] | +0.0929 / +0.3234 (**−0.2305 [−0.3263, −0.0896]**) |
| Broad-16 T = 2×10⁵ | 1,295 / 32 | **+0.0907 [+0.0866, +0.0945]** | +0.1094 [+0.0777, +0.1523] | +0.0187 [−0.0131, +0.0616] | +0.0093 [+0.0081, +0.0105] | — |
| CD4 T = 1×10⁵ Stim48hr | 318 / 388 | **+0.1042 [+0.0963, +0.1067]** | +0.1116 [+0.1082, +0.1165] | +0.0074 [+0.0034, +0.0171] | +0.0032 [+0.0015, +0.0048] | +0.1765 / +0.1737 |

**P5(i) CONFIRMED, everywhere.** A construct that misses its target still gives a reproducible
transcriptome: every failing-construct split-half interval is above zero, in all five libraries and at both
budgets. This generalises 293's two anecdotes (HIRA_P1 +0.268, MED15_P2 +0.047) to 2,400+ failing
constructs. It is the fact behind the "one truth in five" concern, and it holds.

**P5(ii) REFUTED. There is no shared vector signature.** Three independent readings agree:

* **Cosine to the library's mean response** (un-centred, so a common vector/stress component *would* show):
  failing constructs are **not** more aligned to it than working ones — K562 +0.1029 vs +0.1044 (difference
  −0.0015 [−0.0123, +0.0051]), CD4 +0.1765 vs +0.1737. In Jurkat and HepG2 they are **markedly less**
  aligned (−0.0424 and −0.2305, both intervals below zero).
* **Cross-gene agreement between failing constructs, un-centred**: K562 +0.0099 vs working × working
  +0.0099 (difference −0.0001 [−0.0003, +0.0003]); Jurkat +0.0112 vs +0.0130 (**failing lower**);
  HepG2 +0.0115 vs +0.0439 (failing far lower); Broad-16 +0.0236 vs +0.0510. Only CD4 shows a positive
  difference and it is +0.0011 [+0.0005, +0.0018] — a thousandth of a cosine.
* **Cross-gene agreement in the residual space** does favour failing constructs in K562 (+0.0002
  [+0.0000, +0.0005]), Jurkat (+0.0050) and HepG2 (+0.0253) — but this is an artefact of centring, not
  evidence: working constructs are *actively anti-correlated* there (−0.0052 Jurkat, −0.0072 HepG2, −0.0007
  K562) because real, strong, different responses must point apart once the panel mean is removed, while
  failing constructs sit near zero. Reading it as "failing constructs agree with each other" inverts the
  sign of the fact: failing constructs carry **less** gene-specific signal, they do not carry **more** shared
  signal.

So the honest statement is: **a failing construct's response is reproducible and it is not the vector.** It
is reproducible, systematically *less* reproducible than a working construct's, and it points somewhere that
is neither the gene's response nor a generic library-wide direction. What it is — an off-target effect of
that specific spacer, a partial knockdown below the detection threshold on the log1p scale, or a
clone-of-the-transduced-population effect — is not distinguished here.

---

## 6. Library failure rates (`work/rates296.{py,json,log}`)

On-target |lfc| at the gene's own column, full depth, `log1p(5e4·p)` scale, against each library's own
control pool.

| library / cohort | n | median \|lfc\| | \|lfc\| < 1.0 | < 0.5 |
|---|---:|---:|---:|---:|
| **Arc VCC-2025 panel** (293's `kd_all`, 3 splits) | 300 | 1.485 | **18.7 %** | **4.3 %** |
| K562-GWPS single-construct, **2026-panel targets** | 221 | 0.763 | 70.1 % | 24.0 % |
| K562-GWPS single-construct, non-panel | 6,897 | 0.680 | 75.0 % | 30.9 % |
| K562-GWPS, the individual two-construct constructs | 1,116 | 0.320 | 88.0 % | 66.4 % |
| Nadig Jurkat, single-construct targets | 2,051 | 0.949 | 53.4 % | 18.8 % |
| Nadig Jurkat, two-construct constructs | 240 | 0.561 | 72.9 % | 47.9 % |
| Nadig HepG2, single-construct targets | 2,051 | 0.923 | 55.5 % | 19.0 % |
| Nadig HepG2, two-construct constructs | 200 | 0.548 | 75.5 % | 44.5 % |
| **Broad-16 (Cas9)** per (line, target) | 1,032 | 0.237 | 97.7 % | 82.0 % |
| Broad-16 (Cas9) per (line, target, construct) | 2,035 | 0.254 | 96.5 % | 79.0 % |
| CD4 GWCD4i per (condition, guide), my |lfc| | 1,753 | 0.806 | 63.9 % | 29.0 % |
| CD4 GWCD4i, **the authors' own `signif_knockdown`** | 73,765 | — | **26.7 % not significant** (Rest 27.6 %, Stim8hr 25.8 %, Stim48hr 26.6 %); 7,002 rows flagged "high confidence no effect" | — |

**Read this table with three corrections, in order of importance.**

1. **The absolute cut is not comparable across libraries.** Arc's 2026/2025 panels were built with
   *"sufficiently high basal expression of the target gene"* as an inclusion criterion (AGENTS.md §2.1), and
   the `log1p(5e4·p)` scale compresses low-expressed genes: at low `p`, even a complete knockdown cannot
   reach |lfc| 1.0. That is most of the gap between Arc's 18.7 % and K562-GWPS's 70–75 %. **293's absolute
   cut is a property of Arc's expression-selected panel; it should not be exported.** This is exactly why
   C2 and C3 were preregistered, and why §3's headline uses C2.
2. **Broad-16 is Cas9.** Editing produces indels; mRNA drops only when nonsense-mediated decay follows.
   97.7 % of its (line, target) rows read |lfc| < 1.0, which cannot mean 97.7 % of its constructs failed —
   its ρ_guide is 0.79–0.88, the highest of any library. **For Cas9 the on-target mRNA column is a poor
   efficacy readout**, and this is the honest reason P3's C1 stratum has 6 usable pairs.
3. **The transferable number for a CRISPRi library is the authors' own CD4 call: ~26.7 % of guides are not
   a significant knockdown**, on 73,765 guide × condition rows with the study's own statistics, in a
   genome-wide CRISPRi library in a primary human cell type. That sits above 293's 18.7 % for Arc's
   validated 300-gene panel and well below the naive |lfc| reading of any public library. **The "one truth
   in five" estimate survives as a range: roughly one target in five to one in four carries a truth whose
   construct did not knock the gene down**, with Arc's own validated panel at the optimistic end.

*Indicative only* (`work/rates296.json`): of the **264** 2026-panel targets present in K562-GWPS as
single-construct perturbations, **221** have a measurable own-gene column and **70.1 %** read |lfc| < 1.0,
**24.0 %** < 0.5. The 2026 challenge uses different vectors, a different chemistry and different lines;
this is a base rate for CRISPRi against these particular genes, not a prediction about the challenge's
truths.

---

## 7. Descriptive tail (POST-HOC, labelled)

* **Promoter / transcript-level constructs** (`work/promoter296.json`). As §2a records, the construct pairs
  are overwhelmingly P1-vs-P2. Same-promoter pairs vs different-promoter pairs: **K562** ρ_guide +0.193
  (3 pairs) vs +0.002 (713) at n = 25 and +0.605 (2) vs +0.076 (625) at n = 50 — **the same-promoter cells
  have 2–3 pairs and are uninformative**, the bootstrap interval notwithstanding; **Jurkat** +0.172 (17) vs
  +0.197 (103), difference −0.025 [−0.212, +0.606] — **no effect**; **HepG2** +0.641 (17) vs +0.303 (52),
  +0.338 [−0.116, +0.695] — positive, not separated. In Nadig the same-promoter pairs also have a much
  stronger weaker-construct knockdown (median |lfc| 0.634 vs 0.133), so promoter identity and efficacy are
  confounded in exactly the cells that could test it. **Verdict: transcript-level construct choice is not
  shown to explain the split, and the design on disk cannot cleanly test it.**
* **The 20 most and least efficient constructs per library** are in `work/rates296.json`
  (`descriptive_top20`), for all five libraries.
* 293's `KDM2B` pattern (a pair whose second construct is *more* reproducible than the first and points
  elsewhere) recurs at scale: in K562 at n = 50, (a2) +0.0155 exceeds (a) +0.0108 for the ALL row — but the
  (a2) cohort is restricted to pairs whose *smaller* construct has ≥ 2n cells, so it is a larger-construct
  subset, not a paradox. The same restriction applies in 293.

---

## 8. What this does to F14

**F14 stands as a measurement and changes as a claim.** Three amendments follow from §3–§6:

1. **The construct term is partly guide efficacy, in some libraries.** In Nadig Jurkat it splits cleanly at
   all three cuts and both rungs (ρ_guide +0.46…+0.69 where both constructs knock down, ≈ 0 where the weaker
   does not) and knockdown predicts the construct agreement specifically (+0.446) and not the cell split
   (−0.022). In K562-GWPS it splits at the library-relative cut and n = 50 (+0.340 [+0.113, +0.548] vs
   +0.001), but knockdown there predicts the cell split as strongly as the guide split — **285 §6b's caution
   reproduced**. In HepG2 the point estimates split (+0.645 vs +0.307) with an interval that does not
   exclude zero. **In primary CD4 T cells there is no split at all**, at any condition, budget or cut, with
   41–199 pairs per stratum: ρ_guide is 0.02–0.06 whether both constructs knock down or neither does.
   So the headline **"guide efficacy is the dominant irreducible term" is not supported as a general claim**;
   "in immortalised-line CRISPRi libraries, a validated construct pair agrees much better than an
   unvalidated one, and in primary T cells it does not" is.
2. **Arc's own library is the favourable case, not the typical one.** 293 measured +0.687 vs −0.003 there.
   The closest public analogue (Jurkat, C1) reads +0.776 vs +0.147; the two large libraries (K562 748 pairs,
   Broad-16 803 rows) do not reproduce the size of the effect. 293's estimate remains the best number *for
   the challenge's own library*, and this experiment does not extend it to the public sources the pipeline
   actually trains on.
3. **F14's cross-library ρ_guide ordering is partly a cross-gene-null artefact** (§3b). Once (x) is
   subtracted, HepG2 falls from 0.41 to 0.03 and only Broad-16's Cas9 term survives. Any future statement of
   F14 should quote the null-corrected value or state that (b) − (x) is not separated from zero in the
   CRISPRi libraries.
4. **A construct that fails still gives a reproducible transcriptome — and it is not the vector.** P5(i)
   confirmed at scale, P5(ii) refuted (§5). The prediction problem this creates is therefore *harder* than
   293's framing suggested: a failing construct's truth is not a shared, learnable stress direction that a
   model could emit as a fallback.

**Nothing here is build-worthy.** No model, no arm, no slot. The competition-facing consequence is only the
one already recorded: the construct/efficacy ceiling caps every zero-shot entry equally and, at the
incumbent's ≈ 0.04 median transfer cosine, is not the binding constraint (R293 §6).

---

## 9. Honest limits

* **Multiplicity.** 13 (library × budget) cells × 3 cuts × 3 strata = 117 stratum rows, plus 13 P4
  coefficients, 13 failing-construct blocks and 5 promoter cells. The cuts, the strata, the contrast set and
  the five predictions were fixed in `PREREGISTRATION.md` before any of them was computed; no cell was added
  after a look. Individual C2/C3 passes at P ≈ 0.02 (Jurkat n = 50 C3, +0.001 lower bound) should be read
  against that count. The results that survive multiplicity comfortably are Jurkat's C2 at both rungs
  (P ≤ 0.001), K562 n = 50 C2 (P = 0.006), and CD4's null (13/13 cells flat).
* **The two-construct targets are not a random sample.** 293 measured this directly in Arc's library
  (the two-construct constructs are its weakest, p = 0.0015); the same holds in every library here — K562's
  two-construct constructs read median |lfc| 0.320 against 0.680–0.763 for single-construct targets,
  Jurkat's 0.561 against 0.949, HepG2's 0.548 against 0.923. **Every ρ_guide reported here is therefore a
  lower bound on the library-wide value**, and the both-clear stratum is the estimate that transfers.
* **The K562 ALL row is not stably estimated** (§1): (b) is not distinguishable from zero, so ρ_guide's
  point estimate moves 3× between legitimate draws. Only the stratified contrast is quoted as a result.
* **(a2) is measured on a subset.** The smaller construct's own split requires ≥ 2n of its cells, which
  selects the larger "small" constructs; ρ_dis inherits that selection. 293 has the same property.
* **Broad-16 is Cas9 and its cell split mixes the two constructs** (268 built A/B from all assigned cells),
  both of which push its ρ_guide up; its on-target mRNA column is a weak efficacy readout (§6, correction 2).
  Its (a1)/(a2) are binomial **count** splits, which are Poisson-only and read far above a cell split
  (+0.093 vs (a) +0.012 at T = 2×10⁵), so Broad-16's ρ_dis (+0.28–0.32) is biased **down** and is not
  comparable to the CRISPRi libraries' ρ_dis.
* **CD4 has no cells.** Its (a)/(a2) are count splits of pseudobulk rows (285 §6 measured the
  count-minus-cell offset at +0.0438 on the two sources that have both; it is reported, not applied). Its
  guide slots are the bank's two guides per (target, donor, condition), ~8–21 cells per side at these
  budgets. The CD4 null is powered in target count but thin in cells per target.
* **The knockdown measure is one number from one column.** It cannot separate "the guide never repressed"
  from "the guide repressed a transcript this assay does not resolve" from "the gene is too lowly expressed
  for the log1p scale to show it". C3 (knockdown fraction) is the expression-free companion and it agrees
  with C1/C2 on every qualitative call.
* **Nothing here is measured on A/B/C or on the challenge's own reference**; every source is a public
  screen, and the promoter fact (§2a) means "construct" is a TSS unit throughout.
* **One design difference between libraries**: K562 centres and scores on its own 748-pair cohort (285's
  `k562_guide.py` construction, so KA-1 is exact); Nadig centres on the 507–534-target pool (285's
  `nadig_285.py` construction); Broad-16 and CD4 centre per side on their own cohort. Strata within a
  library always share one centring, so §3a is unaffected; §3b's cross-library comparison is not.

---

## 10. Files, memory, provenance

| file (in `work/`) | what | peak RSS / time |
|---|---|---|
| `lib296.py` | the stratification, the ratio-of-medians bootstrap, the **stratum** bootstrap between disjoint strata; everything else imported from 285's `lib285.py` (which imports 218's `_per_target` verbatim) | — |
| `ka296.py` → `ka296.{json,log}` | KA-1…KA-5a | 1.6 s |
| `k562_296.py` → `k562_296.{npz,json,log}` | K562-GWPS: one pass = capped cells + exact full-depth on-target statistics for all 9,866 perturbations | **9.82 GiB** / 219 s |
| `nadig_296.py jurkat\|hepg2` → `nadig_*_296.{npz,json,log}` | Nadig, both lines | 6.02 / 6.12 GiB; 13 / 8 s |
| `broad_296.py` → `broad_296.{npz,json,log}`, `broad_296_cross_*.npz` | Broad-16, 12 lines, two UMI budgets | 2.56 GiB / 122 s |
| `cd4_296.py` → `cd4_296.{npz,json,log}`, `cd4_296_cross_*.npz` | CD4 GWCD4i, 3 conditions × 2 budgets, KA-5b | 2.93 GiB / 128 s |
| `strat296.py` → `strat296.{json,log}`, **`efficacy296.csv`** | the frozen strata, P1/P4, the failing-construct blocks, the 8,144-row csv | 101 s |
| `rates296.py` → `rates296.{json,log}` | §6's failure rates, the descriptive top-20s, the construct-id shapes | 0.1 s |
| `promoter296.json` | §7's same-promoter vs different-promoter cells | — |

Read-only inputs: `data/external/scperturb/ReplogleWeissman2022_K562_gwps.h5ad`,
`data/external/nadig2025/{jurkat,hepg2}_raw_singlecell.h5ad`,
`experiments/268-broad16-context-residuals/work/{halves/*.npz, inventory.json}`,
`experiments/280-gwcd4i-donor-replicates/work/{bank280_X.npy, bank280_rows.npy, obs280.npz, axis280.npz, lib280.py}`,
`data/external/cd4_cell2026/suppl_tables_2026-08-24/guide_kd_efficiency.suppl_table.csv`,
`data/challenge/pert_counts.csv`, and 285's / 293's saved arrays and csv.

**Write audit.** `find experiments data scripts src docs submissions -newer work/.marker-start -type f`
returns, outside this directory, only the two sibling agents' own trees (`experiments/294-writeup/`,
`experiments/297-audit-extension/`) — 29 files, none of them touched by any script here. No `vcc`, no
`cell-eval2`, no build, no download, nothing under `submissions/` or `data/`.
