# 284 — 281 rerun with EVERY membership of a held HGNC family excluded from training, target-cluster bootstraps, and the test that separates "the features cannot predict the remainder" from "the remainder is noise"

**FAMILY-HELD (ALL MEMBERSHIPS): ρ +0.3338 [+0.2164, +0.3545], R² +0.1121 [+0.0472, +0.1247]
(281: +0.320 / 0.111); vs reliability alone −0.0348 [−0.0815, +0.0330]; REMAINDER REPRODUCIBLE
ACROSS DISJOINT SOURCE SETS: YES (ρ +0.3277 [+0.2413, +0.3730]; ceiling +0.3816 [+0.2667, +0.4098])**

Three statements, all measured:

1. **281's split defect is real and large — and it did not cause 281's result.** 43.7 % of 281's
   test genes share an HGNC family with a training gene; 784 training rows carry one of the 359 held
   families as a non-first group, touching 130 of them. Closing the hole completely (S1: a gene is
   TEST if ANY of its groups is held; 0 training rows carry a held family) moves the primary cell
   **up**, from +0.3204 / R² +0.1113 to **+0.3338 / +0.1121**. The leak was not inflating the number.
2. **281's family-level result survives; its zinc-finger figure does not.** Under folds in which a
   test gene shares no family at all with its training rows, out-of-fold predicted family medians
   still reproduce 279's axis at **+0.8047 [+0.6786, +0.8828]** over 78 families (+0.8039 over 281's
   51 first-group families). But 279's zinc-finger deficit is **57 % predicted, not 84 %**
   (−0.01009 of −0.01785); 281's 84 % is recovered exactly from its own out-of-fold column and falls
   to 80 % / 57 % as the folds are tightened.
3. **281's conclusion "the gene-level residual is mostly measurement" is REFUTED at this depth.**
   Split the 28 sources into two halves sharing no source and no cell line, and the gene-level
   remainder — what is left of each half's conservation map after the 21 build-time features have
   predicted it out-of-fold — reproduces across the two halves at **ρ +0.3277 [+0.2413, +0.3730]**
   against a ceiling of +0.3816 and a within-component permutation null of +0.0435 ± 0.0181 (0/200
   draws reach it, z +15.7). It is positive on **20/20** admissible partitions (median +0.3191), it
   survives subtracting measured reliability and response magnitude as well (+0.2774 [+0.1747,
   +0.3036]), and it is not the perturbed gene's own knockdown (genes OFF the 2,739-gene axis:
   +0.3224 [+0.2391, +0.4110]). **The remainder is a gene property the features do not predict.**

**Nothing here is buildable and none is proposed.** Measurement only: no build, no `vcc`, no
`cell-eval2`, no download, nothing under `submissions/` or `data/`, no production file touched, no
file written outside this directory (`find` over the repo for files newer than `PREREGISTRATION.md`
returns only another agent's `experiments/287-*`). The spent reserves were not read: neither
`experiments/267-*/work/reserve/` nor JVE127 / MG63 / SLR23 / UMRC3. `uv run python` throughout.
Peak RSS **0.42 GiB** (`ka284.py`, `/usr/bin/time -l` 452,886,528 B; `extra284` 0.30 GiB, `fam284`
0.25 GiB) — the ~6 GiB cap was never approached.

**Date 2026-09-06, 13:46–14:15 UTC.** Deterministic, seed **20260906**. `PREREGISTRATION.md` sha256
**`51f30071220d21529161e67a389d35f0b412bb9dd792efa1e754449390e3e017`**, frozen **2026-09-06
13:46:19 UTC** (`work/PREREGISTRATION.sha256`) — before any fit on a new split and before any
number of S1, S2 or the remainder test existed. Everything labelled POST-HOC was written after the
preregistered numbers were read and cannot change a verdict.

---

## 0. Known answers

| row | required | measured | verdict |
|---|---|---|---|
| **KA-1** 281's primary cell, on 281's OWN split, from the copied table and copied functions | test set = `pred281.csv` `split_family` (3,422 genes, 481 blocks); ridge ρ **+0.320420**, R² **+0.111332**; gbm **+0.315974**; reliability alone **+0.351664**; all at 0.0e+00 | split matches on **17,106/17,106** genes, 3,422 / 481; **+0.320420 / +0.111332 / +0.315974 / +0.351664**; max \|Δ\| **0.0e+00** | **PASS** |
| **KA-2** `c_held` re-derived over all 364 cross-line pairs by the half-map code (`lib284.cmap`) | max \|Δ\| 0.0e+00, n 17,106, median +0.0186 | **0.0e+00**, 17,106, +0.018605 (a float64 recompute differs by 1.1e-07 — 279 stores `c_held` as float32) | **PASS** |
| **KA-3** 279's balanced core off `znf279.csv` | n 8,541; ZNF +0.011813, non-ZNF +0.029663, Δ −0.017850 | **8,541 / +0.011813 / +0.029663 / −0.017850** | **PASS** |
| **KA-4** shuffled training target on S1, 200 independent shuffles (281's corrected KA-4′ form) | the observed ρ outside the null's spread | null **−0.0117 ± 0.1075**; **0/200** draws reach +0.3338; z **+3.2** | **PASS** |
| **regression** 281's own out-of-fold core column re-read here | R281 §5 ρ +0.3398; §6 family axis +0.8115; ZNF predicted −0.01506 = 84 % | **+0.3398**, **+0.8115**, **−0.01506 = 84.4 %** | **PASS** |

---

## 1. The leak, quantified before it was fixed (mandate item 2)

Held genes by number of HGNC groups (`gene_group_all`, pipe-separated): **0 → 3,135 · 1 → 9,237 ·
2 → 3,557 · 3 → 893 · 4 → 211 · 5 → 56 · 6 → 10 · 7 → 4 · 8 → 3**. **4,734 of 17,106 (27.7 %) carry
two or more.**

Under 281's split (blocks = first group, else `NOFAM:<cytoband>`; 481 of 2,425 blocks held, of which
**359 are real families**):

| quantity | count |
|---|---:|
| TEST genes sharing at least one HGNC family with a TRAIN gene | **1,497 of 3,422 (43.7 %)** |
| TRAIN genes carrying one of the 359 held families as a **non-first** group | **784** |
| held families with at least one training row (i.e. not actually held) | **130 of 359 (36.2 %)** |

So "whole HGNC families held out" was true of 229 of the 359 held families and of 56.3 % of the test
genes. For scale, the random-gene split leaks 2,762 of 3,421 (80.7 %).

---

## 2. The splits (PREREG §4)

| split | test genes | units held | TRAIN genes carrying a held family | TEST genes sharing ANY family with TRAIN | test spans |
|---|---:|---|---:|---:|---|
| 281, first group | 3,422 (20.005 %) | 481 of 2,425 blocks | **784** | **1,497** | 273 components / 481 blocks |
| **S1, all memberships (PRIMARY)** | **3,424 (20.016 %)** | 374 blocks taken, **279 held families** | **0** | 1,035 | 207 components / 498 blocks |
| S2, component | 3,423 (20.011 %) | 717 of 1,393 components | — | **0** | 717 components |
| random, seed 20260907 | 3,421 (19.999 %) | — | — | 2,762 | 694 components |

S1's residual 1,035 is genes that carry a **non-held** family also present in training — a different
object from the leak above (no held family is contaminated), and it is what S2 removes entirely.
S1 and 281 agree on 2,878 test genes; 546 are only in S1 and 544 only in 281 (S1 takes fewer blocks
because each held family now pulls in more genes).

**The component structure (the reason S2 is not simply "a stricter S1").** The gene–family graph on
the 17,106 held genes has **1,393 connected components**, and **one of them holds 10,817 genes
(63.2 %)** and spans 1,309 distinct HGNC families — multi-family genes chain the families together.
Sizes: 281 singletons, 698 of 2–4, 390 of 5–19, 22 of 20–99, 2 of ≥ 100 (median 3). Removing the
cytoband edges for family-less genes changes nothing at the top (729 components, same giant).

Two consequences, both measured and both stated with every S2 number below:

* **S2's test set cannot contain the giant component** (it would be 63 % of the genes), so S2 holds
  out only peripheral clusters: its largest families are keratin-associated proteins (68), olfactory
  receptor families 4/10/6/52 (127 together), solute carrier family 22, gap junctions. Its genes are
  systematically less measured and less essential than S1's — median `n_pairs` **3 vs 5**, mean
  `r_bar` **0.0449 vs 0.0579**, mean `replogle_essential` **0.086 vs 0.133**. **S2 is not a stricter
  estimate of the same quantity; it estimates the model on a different, peripheral gene population.**
* **The frozen component bootstrap is near-degenerate on S1**: 63.2 % of S1's test genes sit in the
  single giant component, so the "cluster" resample effectively has one large unit plus 206 small
  ones. It is reported as frozen and next to 281's first-family block bootstrap throughout; see §3.1.

---

## 3. The preregistered tables (mandate items 3–4)

### 3.1 `C_held`, the primary target — S1, S2, 281's split, random

Point estimates are identical under both bootstraps (same fit); the intervals differ.

**Component bootstrap (frozen, PREREG §5), 2,000 draws:**

| arm | **S1 (primary)** | S2 (component) | 281 (first group) | random |
|---|---|---|---|---|
| **ridge, 21 features** | **+0.3338 [+0.2164, +0.3545]** | +0.2155 [+0.1730, +0.2550] | +0.3204 [+0.2014, +0.3495] | +0.2522 [+0.1218, +0.2839] |
| R² | **+0.1121 [+0.0472, +0.1247]** | +0.0470 [+0.0286, +0.0638] | +0.1113 [+0.0473, +0.1269] | +0.0763 [+0.0285, +0.0888] |
| AUC bottom decile | **0.7027 [0.6541, 0.7231]** | 0.6945 [0.6709, 0.7193] | 0.7030 [0.6667, 0.7253] | 0.7107 [0.6510, 0.7250] |
| gbm (depth 3) | +0.3296 [+0.2032, +0.3616] | +0.2242 [+0.1670, +0.2810] | +0.3160 [+0.1758, +0.3515] | +0.2569 [+0.1586, +0.2812] |
| best single (`replogle_essential`, calibrated) | +0.3721 [+0.2347, +0.4082] | +0.2850 [+0.2280, +0.3417] | +0.3669 [+0.2540, +0.3990] | +0.3095 [+0.2444, +0.3243] |
| **reliability alone** (`r_bar`) | **+0.3722 [+0.2418, +0.3969]** | +0.2818 [+0.2415, +0.3217] | +0.3517 [+0.2255, +0.3824] | +0.3314 [+0.2433, +0.3538] |
| features + reliability, ridge | +0.4250 [+0.2922, +0.4493] | +0.3284 [+0.2851, +0.3719] | +0.4075 [+0.2762, +0.4395] | +0.3744 [+0.2764, +0.3957] |
| features + reliability, gbm | +0.4346 [+0.2953, +0.4591] | +0.3294 [+0.2818, +0.3778] | +0.4148 [+0.2823, +0.4484] | +0.3866 [+0.2861, +0.4083] |
| mean predictor (R²) | −0.0004 | −0.0028 | −0.0007 | −0.0000 |

**Paired margins on the same component resample:**

| Δρ | S1 | S2 | 281 | random |
|---|---|---|---|---|
| features − reliability alone (ridge) | **−0.0348 [−0.0815, +0.0330]**, P(Δ≤0) 0.913 | −0.0668 [−0.1135, −0.0206], 0.998 | −0.0309 [−0.0798, +0.0170], 0.926 | −0.0947 [−0.1916, −0.0506], 1.000 |
| features − reliability alone (gbm) | −0.0388 [−0.0924, +0.0276], 0.898 | −0.0483 [−0.1033, +0.0018], 0.970 | −0.0378 [−0.1014, +0.0187], 0.932 | −0.0730 [−0.1478, −0.0264], 0.995 |
| (features + reliability) − reliability, ridge | **+0.0528 [+0.0400, +0.0657]**, 0.000 | +0.0465 [+0.0320, +0.0598], 0.000 | +0.0547 [+0.0393, +0.0657], 0.000 | +0.0396 [+0.0158, +0.0505], 0.000 |
| (features + reliability) − reliability, gbm | +0.0594 [+0.0253, +0.0760], 0.003 | +0.0575 [+0.0333, +0.0803], 0.000 | +0.0586 [+0.0265, +0.0742], 0.002 | +0.0582 [+0.0181, +0.0763], 0.002 |
| features − best single | −0.0400 [−0.1095, +0.0186], 0.946 | −0.0673 [−0.1066, −0.0264], 1.000 | −0.0515 [−0.1159, −0.0138], 0.991 | −0.0747 [−0.1533, −0.0338], 1.000 |

**The two bootstraps compared (S1, ridge on features).** Component (frozen): sd **0.0388**,
interval **[+0.2164, +0.3545]** — strongly left-skewed, the point estimate sitting near the upper
end, because 63 % of the test genes are one resampling unit and draws that under-represent it land
on the peripheral population. 281's first-family block bootstrap on the same fit: sd **0.0318**,
**[+0.2713, +0.3921]**, near-symmetric. On the S2 split, where the clusters are all small, the two
agree almost exactly (component sd 0.0216, [+0.1730, +0.2550]; first-family sd 0.0216, [+0.1711,
+0.2578]). **Reading: the component bootstrap is the honest cluster unit but it is degenerate on any
split that keeps the giant component in the test set — its width there is a statement about the
family graph, not about the estimate.** Neither interval touches zero on any fitted arm (the mean
predictor has no rank), and no paired margin changes sign between the two.

**What the frozen rule says.** 281's PREREGISTRATION §6 rule, applied to S1: the better of
ridge/gbm has ρ ≥ 0.20 with its interval entirely above 0 and R² > 0 with its interval entirely
above 0 → **YES**, under both bootstraps. **281's verdict is not changed by the fix.**

### 3.2 `C_rel` (reliability-corrected) and the balanced core

Component bootstrap; ridge/gbm on the 21 features.

| target | split | ridge ρ | R² | reliability alone | Δρ features − reliability |
|---|---|---|---|---|---|
| `C_rel` | **S1** (n 2,119) | +0.1959 [+0.1376, +0.2450] | +0.0209 [+0.0047, +0.0273] | +0.2197 [+0.0986, +0.2433] | −0.0060 [−0.0543, +0.1097] |
| `C_rel` | S2 (n 1,986) | +0.1576 [+0.1122, +0.2003] | +0.0112 | +0.1302 [+0.0829, +0.1722] | +0.0276 [−0.0334, +0.0908] |
| `C_rel` | 281 (n 2,098) | +0.2117 [+0.1567, +0.2529] | +0.0227 | +0.1930 [+0.0885, +0.2199] | +0.0319 [−0.0154, +0.1271] |
| `C_rel` | random (n 2,163) | +0.1591 [+0.0142, +0.1920] | +0.0095 | +0.1556 [+0.0507, +0.1813] | −0.0099 [−0.1202, +0.0447] |
| **core** | **S1** (n 1,633) | **+0.4133 [+0.2103, +0.4423]** | **+0.2427 [+0.1240, +0.2583]** | **+0.5294 [+0.4024, +0.5554]** | **−0.1429 [−0.2859, −0.0913]** |
| core | S2 (n 1,385) | +0.2125 [+0.1557, +0.2644] | +0.1327 | +0.3978 [+0.3460, +0.4510] | −0.1860 [−0.2466, −0.1251] |
| core | 281 (n 1,600) | +0.3816 [+0.1943, +0.4178] | +0.2406 | +0.5200 [+0.3938, +0.5488] | −0.1569 [−0.2740, −0.1114] |
| core | random (n 1,731) | +0.3719 [+0.2176, +0.3982] | +0.2173 | +0.4435 [+0.3089, +0.4672] | −0.0814 [−0.1725, −0.0290] |

On the measurement-matched core, S1 reads **higher** than 281 (+0.4133 vs +0.3816, R² +0.2427 vs
+0.2406) and reliability alone still wins by ρ +0.14; the features add
**+0.0148 [−0.0089, +0.0290]** on top of it (281: +0.0119 [−0.0029, +0.0257]). **281 §5's
conclusion — once the measurement is matched the gene features add essentially nothing to measured
reliability — survives the fix unchanged.**

---

## 4. The family axis under S1 (mandate item 5)

Grouped 5-fold out-of-fold predictions in which a test gene shares **no** family with its training
rows: the folds are whole components, so fold sizes are **[10,817 · 1,573 · 1,572 · 1,572 · 1,572]**
— the giant component is necessarily one fold, and its 10,817 genes are predicted by a model trained
on the other 6,289 only. Overall out-of-fold quality is nevertheless unchanged: core ρ **+0.3337**
(281-style first-family folds +0.3337; 281's own column +0.3398), raw ρ **+0.2630** (vs +0.2612 /
+0.2641).

| fold construction | families ≥ 20 genes | Spearman(observed family median, predicted family median) |
|---|---:|---|
| **components, ALL memberships (primary)** | **78** | **+0.8047 [+0.6786, +0.8828]** |
| components, first group (281's 51 families) | 51 | +0.8039 [+0.6386, +0.8961] |
| first-family folds, my code | 51 | +0.7929 |
| 281's own out-of-fold column (regression) | 51 | +0.8115 (R281 §6: +0.8115) |

Top and bottom of the axis (balanced core, all memberships): Mediator complex +0.2748 (predicted
+0.1028), SSU processome +0.1925 (+0.1192), spliceosomal C complex +0.1397 (+0.1151), spliceosomal
P complex +0.1200 (+0.1063), nucleoporins +0.1189 (+0.1190), DEAD-box helicases +0.1143 (+0.1139) …
NKL homeoboxes +0.0080 (+0.0188), **KRAB domain containing +0.0079 (+0.0224)**, TBC1 domain +0.0071
(+0.0295), V-set domain +0.0007 (+0.0282). **279's core-machinery → cell-identity ordering is
reproduced from build-time biology under folds that hold every family membership out.**

**The zinc-finger deficit is where 281's number moves.**

| prediction | observed Δmedian | predicted Δmedian | share |
|---|---:|---:|---:|
| 281's own out-of-fold column (regression) | −0.01785 | −0.01506 | **84 %** |
| first-family folds, my code | −0.01785 | −0.01436 | 80 % |
| **component folds (no family shared with training)** | **−0.01785** | **−0.01009** | **57 %** |
| component folds, raw target | −0.01047 | −0.00352 | 34 % |

**Stated with its caveat:** under the strict folds a zinc finger is predicted by a model that has
seen no member of any family it belongs to — which is the point — but that model is also trained on
the 6,289 non-giant genes only, so the drop from 84 % to 57 % mixes "no family leak" with "a smaller,
differently composed training set". The direction is unambiguous and the honest reading of F12/F-C
is **"about half to four fifths of 279's ZNF deficit is predictable from essentiality and
expression"**, not 84 %.

---

## 5. The remainder test (mandate item 6) — the distinction the review asked for

**Construction.** The 28 sources carry 21 distinct lines; three lines repeat (**K562 ×5, A549 ×3,
H1 ×2**), so a partition assigns whole line groups. 20 admissible partitions (each half ≥ 8 lines
and ≥ 20 cross-line pairs) were drawn in 22 tries with seed 20260906. **Partition #1 (primary):**
half 1 = DLD-1, HCT116, HEK293T, HT-29, K562, MCF7, RPE1, THP-1, iPSC-Tian19, melanoma (10 lines,
**81 pairs**); half 2 = A549, BxPC-3, Calu-3, H1, HAP1, HEK293, HepG2, Jurkat, TeloHAEC,
cardiomyocyte, iPSC-neuron (11 lines, **87 pairs**). **The two halves share no source and no cell
line.** 2,848 genes have ≥ 1 pair in both halves.

| quantity (primary partition, component block bootstrap, 2,000 draws) | value |
|---|---|
| **ceiling** — raw Spearman(C⁽¹⁾, C⁽²⁾) | **+0.3816 [+0.2667, +0.4098]** |
| **remainder** — features subtracted out-of-fold on each half's own map | **+0.3277 [+0.2413, +0.3730]** (sd 0.0296, n 2,848) |
| remainder, family-median-only prediction subtracted | +0.3807 [+0.2705, +0.4133] |
| null — remainder⁽²⁾ permuted **within components**, 200 draws | +0.0435 ± 0.0181; **0/200** reach the observed; z **+15.7** |
| attenuation-corrected remainder (sub-halves of each half's pairs, Spearman–Brown; half 1 rel +0.833, half 2 +0.667) | **+0.4396** (conservative: sub-halves share sources, so the reliabilities are biased up) |
| attenuation-corrected ceiling, same construction | +0.5050 |
| the feature model's own out-of-fold ρ on each half's map | +0.238 / +0.239 |

**Over 20 admissible partitions:** remainder median **+0.3191**, range [+0.0981, +0.4641],
**20/20 > 0** and **20/20 ≥ 0.05**; ceiling median +0.3586, range [+0.0742, +0.5658]; the
family-median remainder median +0.3586. The one weak partition (#4, remainder +0.0981, ceiling
+0.0742) is the one whose halves overlap on only 1,345 genes.

**Strata (descriptive, same partition).** The remainder correlation rises with how well each gene is
measured, which is the signature of attenuation by noise, not of a noise artefact:

| stratum | n | remainder ρ | ceiling ρ |
|---|---:|---:|---:|
| `r_bar` low tercile | 948 | +0.2067 | +0.2752 |
| `r_bar` mid | 948 | +0.3176 | +0.3477 |
| `r_bar` high | 948 | +0.3402 | +0.3488 |
| min pair count 1 | 1,013 | +0.1381 | +0.2458 |
| min pair count 2–4 | 1,013 | +0.4106 | +0.4314 |
| min pair count ≥ 5 | 822 | +0.4476 | +0.4465 |

**The frozen pass rule (PREREG §7).** The interval **excludes 0** and the point estimate
**+0.3277 ≥ 0.05**; the 20-partition median +0.3191 agrees in sign and stays above 0.05, so the
result is not EQUIVOCAL. → **"THE REMAINDER IS A GENE PROPERTY THE FEATURES DO NOT PREDICT."**

### 5.1 Which of 281's two readings survives

281 wrote: *"Conservation is a family-level biological property; the gene-level residual is mostly
measurement."* The first clause survives (§4). **The second is refuted at this depth**: what the 21
build-time features fail to predict about a gene's conservation is reproduced by an independent set
of cell lines and datasets at ρ +0.33 against a +0.38 ceiling — the remainder retains **86 % of the
ceiling's cross-half rank correlation** (+0.3277 of +0.3816; a ratio of Spearmans, not a variance
share). 281 was entitled to say
"our features cannot predict the gene-level remainder"; it was not entitled to say the remainder is
noise, and it is not.

### 5.2 Three controls the conclusion needs (POST-HOC, labelled)

| control | remainder ρ (component bootstrap) |
|---|---|
| frozen: 21 features subtracted | +0.3277 [+0.2496, +0.3704] |
| **+ measured reliability `r_bar` subtracted** | **+0.2774 [+0.1747, +0.3036]** |
| + `r_bar` + response magnitude subtracted | +0.2768 [+0.1715, +0.3005] |
| `r_bar` alone subtracted | +0.3171 [+0.2166, +0.3513] |
| `r_bar` + magnitude alone subtracted | +0.3211 [+0.2184, +0.3529] |
| outside the giant component (n 734) | +0.3112 [+0.2265, +0.3845] (ceiling +0.3468) |
| inside the giant component (n 2,114) | +0.3373 (one cluster: no component interval exists) |
| **genes OFF 279's 2,739-gene axis** (self-knockdown impossible), n 1,816 | **+0.3224 [+0.2391, +0.4110]** (ceiling +0.4055) |
| genes ON the axis (self-knockdown possible), n 1,032 | +0.2920 [+0.1418, +0.3234] (ceiling +0.2685) |

(The first row is the frozen statistic recomputed in `extra284.py` with a different bootstrap seed:
the same point estimate +0.3277, interval [+0.2496, +0.3704] against the frozen script's [+0.2413,
+0.3730] — a Monte-Carlo difference of the 2,000-draw bootstrap, not of the estimate. The frozen
interval is the one reported everywhere else.)

`r_bar` and the magnitude are measured on **all 28 sources**, i.e. on both halves, so subtracting
them removes information the two halves share — the test is conservative in that direction, and the
remainder still reads +0.28. The axis control matters because 279's statistic is a cosine of
residual response rows on a common 2,739-gene axis and does **not** drop the perturbed gene's own
column: a gene on that axis contributes its own on-target knockdown to every pair. Only 15.3 % of
held genes are on it, and the off-axis genes — for whom the cosine is downstream response only —
carry the effect at full strength. **The reproducible remainder is neither reliability, nor response
size, nor the knockdown of the gene itself.**

---

## 6. The reproducible remainder, described (POST-HOC, item 7)

The 20 genes furthest above and below prediction **in both halves independently** (sorted by the
weaker of the two halves' residuals; `r1 / r2`, pair counts, Replogle-essential flag):

**Above prediction — conserved more than essentiality and expression predict.** `EIF3H` (+0.406 /
+0.436), `ZC3H13` (+0.384 / +0.297), `OGFOD1`, `NARS`, `TANGO6`, `RPAP1`, `C14orf178`, `TRRAP`,
`SMG7`, `SUDS3`, `SMG5`, `RBMX2`, `BOD1L1`, `ZFC3H1`, `PPRC1`, `TCERG1`, `CASC3`, `CMTR1`, `UPF2`,
`MED7`. **Four are core nonsense-mediated-decay / exon-junction factors** (`SMG5`, `SMG7`, `UPF2`,
`CASC3`) and six more decide an mRNA's fate at another step (`ZC3H13` m6A, `ZFC3H1` PAXT nuclear
decay, `RBMX2` spliceosome, `CMTR1` cap 2′-O-methylation, `TCERG1` elongation–splicing coupling,
`EIF3H` translation initiation): perturbing the machinery that decides what happens to an mRNA
produces the same downstream response in every line measured. The rest of the tail is transcription
and translation machinery (`TRRAP`/SAGA, `MED7`/Mediator, `RPAP1`, `SUDS3`/SIN3, `NARS`, `OGFOD1`).

**Below prediction — conserved less than their essentiality predicts.** `CDC37` (−0.131 / −0.127),
`SAP30BP`, `WDR74`, `NCL`, `TUBGCP6`, `ATP5MF`, `TBCC`, `TAF1C`, `STX5`, `VCP`, `DTYMK`, `VPS37A`,
`NCAPG`, `PSMB4`, `NCAPD2`, `BOP1`, `POP5`, `HAUS6`, `GPSM2`, `IRS2`. **Chaperone / proteostasis and
chromosome-segregation / ribosome-biogenesis machinery** — `CDC37` (HSP90 kinase co-chaperone),
`VCP`, `PSMB4`, `NCAPG` + `NCAPD2` (condensin I), `HAUS6` (augmin), `TUBGCP6` (γ-TuRC), `BOP1`,
`WDR74`, `NCL`, `POP5`. Their response is large enough to be measured and is *not* the same across
lines: the recipient's proliferative and proteostatic state conditions it. This is the same class
281 §6.1 found in its raw-target tail (`CCNA2`, `CHEK1`, `CDC16`, `CDCA5`, `NCAPD3`), now recovered
from two disjoint source sets that agree.

Both tails are essential (mean `replogle_essential` 0.75 above, 0.85 below) and complex-bound (mean
hu.MAP complexes 2.65 / 3.40): **the reproducible remainder separates within the essential genes,
which is exactly where the essentiality × expression axis of §4 has no resolution.** Caveat: several
above-tail genes have only one pair in half 2 (the "min pairs 1" stratum reads ρ +0.14, the weakest),
so individual rows are noisier than the aggregate; the asymmetry in magnitude (+0.18…+0.44 above
vs −0.08…−0.22 below) is the floor of a bounded statistic whose prediction sits at +0.02…+0.10.

---

## 7. Limitations, stated plainly

1. **The component bootstrap is degenerate wherever the giant component is in the test set** (63 %
   of S1's test genes, 0 % of S2's). It is reported as frozen; the first-family block bootstrap is
   reported beside it; they agree where the clusters are small. No verdict in this document depends
   on which of the two is used.
2. **S2 is not a stricter S1.** Holding whole components necessarily excludes the giant one from the
   test set, so S2 measures the model on keratin-associated, olfactory-receptor and solute-carrier
   clusters that are less measured and less essential than average. Its lower ρ (+0.2155) is a
   population difference, not evidence that S1 is inflated.
3. **The strict 5-fold is unbalanced by construction** (one fold of 10,817). §4's ZNF number carries
   that caveat explicitly.
4. **The remainder test is about reproducibility across source sets, not about a new experiment.**
   Both halves come from the same 28 banks on disk and the same 2,739-gene axis; a processing
   artefact shared by, say, all Replogle-derived banks would appear in both halves. The line-group
   partition removes shared sources and shared lines, not shared pipelines.
5. **"Gene property" is a statement about reproducibility, not about mechanism.** The controls in
   §5.2 exclude measured reliability, response magnitude and the self-knockdown column; they do not
   exclude every gene-level technical covariate (e.g. guide-library composition correlated with gene
   identity across studies).
6. **The remainder's absolute size is small.** The ceiling itself is +0.38 (+0.51 attenuation-
   corrected): most of the variance of a single-gene conservation value is still sampling noise, and
   the remainder is a rank statistic over a bounded quantity whose median is +0.019.
7. **Multiplicity.** Verdict cells: 3 (S1 ridge, S1 gbm, the remainder test's primary partition).
   Reported in total: ~90 cells across 4 splits × 3 targets × 8 arms × 2 bootstraps, the family axis
   in four fold constructions, 20 partitions, 5 stronger-predictor controls, 2 axis strata, 6 strata
   and the two tails — all but the three labelled secondary, comparative or post-hoc.
8. **One map, one statistic.** Everything is 279's `C_held` on its 28 banks. The conclusions are at
   that scope.

---

## 8. What this is worth

**Scientific track — three amendments to F12.**

* **F12-a (correction of the split, no change of verdict).** 281's family-held split leaked
  (43.7 % of test genes shared a family with training; 784 training rows carried a held family).
  With every membership of a held family excluded, the primary cell is **ρ +0.3338 [+0.2164,
  +0.3545] / R² +0.1121**, marginally **above** 281's +0.3204 / +0.1113. The reliability comparison
  is unchanged (features − reliability −0.0348 [−0.0815, +0.0330]; +0.0528 [+0.0400, +0.0657] added
  to it; +0.0148 [−0.0089, +0.0290] on the measurement-matched core). **F12's numbers stand; its
  "whole families held out" claim is now true.**
* **F12-b (the family axis survives, the ZNF figure moves).** Out-of-fold predicted family medians
  reproduce 279's axis at **+0.8047 [+0.6786, +0.8828]** over 78 families under folds that share no
  family with training. But 279's zinc-finger deficit is **57 % predicted (−0.01009 of −0.01785)**,
  not 84 %; the 84 % is exactly recovered from 281's own column and is a property of its leakier
  folds. **F12 and F9's headline "84 % predicted from biology alone" should read "57–80 %,
  depending on how completely the gene's families are held out".**
* **F12-c (new, and the reason the experiment was run).** The gene-level remainder is **not
  measurement noise**: across two halves of the source set that share no source and no cell line it
  reproduces at **ρ +0.3277 [+0.2413, +0.3730]** (ceiling +0.3816; null +0.044 ± 0.018, 0/200;
  20/20 partitions positive; +0.2774 after reliability and magnitude are also subtracted; +0.3224 on
  genes whose own knockdown cannot contribute). **There is real, reproducible, per-gene conservation
  biology beyond the family axis — 21 build-time annotation features simply do not carry it.**
  The described tail says what it looks like: mRNA-surveillance machinery above the essentiality
  prediction, proteostasis and chromosome-segregation machinery below it.

**What it changes for the programme.** F12's "for a model, this is a family-level prior, not a
per-gene one" was a statement about *these features*, not about the quantity. The per-gene signal
exists and is worth about as much rank information as the whole feature model (+0.33 reproducible
vs +0.33 predicted). Getting at it needs different features — response-derived or pathway-level ones
— or the depth axis 276/278/281 all pointed at; it does not need more annotation blocks (regulon
size and TF status read at chance in 281 and nothing here contradicts that).

**Competition track: nothing.** No arm, no dose, no rendering, no model use is proposed. The
per-gene R² is 0.11 (0.24 measurement-matched) and 277 already closed the one action a
conservation-keyed rendering could take.

---

## Files

| path | what |
|---|---|
| `PREREGISTRATION.md` + `work/PREREGISTRATION.sha256` | the frozen splits, models, metrics, uncertainty, remainder test and pass rule (sha256 `51f30071…`, 2026-09-06 13:46:19 UTC) |
| `lib284.py` | the functions COPIED from 281 (with their source lines) + the S1 draw, the gene–family graph, components, component folds |
| `ka284.py` → `work/ka284.json`, `work/ka284.log`, `work/splits284.npz` | KA-1…KA-3, the leak quantification, the four splits and the component structure |
| `fit284.py` → `work/fit284.json`, `work/fit284.log`, `work/preds284.npz` | the 4 splits × 3 targets × 8 arms tables, both bootstraps, KA-4 |
| `fam284.py` → `work/fam284.json`, `work/oof284.npz` | the family axis under the strict folds, the four fold constructions, the ZNF deficit |
| `rem284.py` → `work/rem284.json`, `work/rem284.npz` | the remainder test: partitions, ceiling, remainder, attenuation, null, 20 partitions, the pass rule |
| `extra284.py` → `work/extra284.json`, **`work/pred284.csv`** | the stronger-predictor controls, the giant-component split, the tails, the deliverable |
| `axis284.py` → `work/axis284.json` | the self-knockdown control (on/off 279's axis) |
| **`work/pred284.csv`** | 17,106 rows × 29 columns — `gene`, observed (`C_held` / `C_rel` / core), `r_bar`, `n_pairs`, `gene_group_all`, `component`, `component_fold`, the four split labels, the S1 and S2 predictions (ridge / gbm / reliability / features+reliability), the strict out-of-fold predictions, `C_half1` / `C_half2`, their pair counts and `remainder_half1` / `remainder_half2` |

Reproduce, from this directory: `PYTHONDONTWRITEBYTECODE=1 uv run python {ka,fit,fam,rem,extra,axis}284.py`
in that order.

---

*Coordinator addendum, 2026-09-06 15:15Z (not written by 284).* 290 (`experiments/290-remainder-response-features/RESULTS.md` §6) corrects one sentence of §5.2: `r_bar` is measured on FOUR banks (`k562_gwps`, `orion_hct116`, `orion_hek293t`, `dld1_lowmoi`; `rel279.py:44-49,222`), all in the primary partition's half 1 — not "on all 28 sources"; 279's `mag` is. No verdict here changes (`r_bar` is subtracted from both halves symmetrically). 290 then showed the remainder IS carried by response-derived features measured on the other half (ρ +0.317, at this experiment's +0.328 ceiling) — the test §8 asked for.
