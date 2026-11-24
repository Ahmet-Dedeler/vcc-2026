# 290 — what carries 284's reproducible per-gene conservation remainder: the responses themselves, measured on one half of the source set

**RESPONSE-DERIVED FEATURES PREDICT THE HELD-HALF REMAINDER: YES (ρ +0.3174 [+0.2897, +0.3914];
annotation −0.0566 [−0.1380, −0.0205]; C⁽¹⁾ alone +0.2948 [+0.1989, +0.3527]); BEAT ANNOTATION:
YES (paired Δ +0.3948 [+0.3442, +0.4896]); PROLIFERATION / STRESS ENGAGEMENT → LESS CONSERVED:
NO (ρ +0.1184 [+0.0999, +0.2327] — the opposite sign); GENERICNESS → LESS CONSERVED: NO
(ρ +0.1147 [+0.0694, +0.1878] — the opposite sign).**

Four statements, all measured:

1. **284's remainder is predictable from the responses of a disjoint set of cell lines.** Features
   computed on half 1's 14 sources predict half 2's remainder at **ρ +0.3174**, against an
   annotation model at **−0.0566** and the empirical ceiling — half 1's OWN remainder, 284's
   +0.3277 — which the response features **match** (paired Δ +0.0048 [−0.0329, +0.1010]). The
   symmetric direction reads **+0.3411**, and the effect is positive on **40/40** partition ×
   direction cells (median +0.2479, 39/40 ≥ +0.10).
2. **It is not magnitude, not reliability, not coverage, and not a leak.** `logmag_h1` alone reads
   **+0.0453 [−0.0558, +0.0623]** — not distinguishable from 0; partialling magnitude out of both
   sides leaves +0.3162, and magnitude + `r_bar` leaves +0.2832. `r_bar` alone +0.1633, coverage
   alone +0.1381; the model beats both (Δ +0.1722 / +0.1924). 200 target shuffles read
   −0.0070 ± 0.0216 with 0/200 reaching the observed (z +15.0).
3. **The preregistered mechanism is REFUTED with its sign inverted.** The hypothesis was that
   responses engaging the recipient's proliferation / stress axis are the state-dependent ones and
   generic responses the less conserved ones. Both read **positive**: proliferation engagement
   +0.1184, genericness (top-5-PC energy) +0.1147, top-20-PC energy +0.1578. **A response that is
   broad, that lies in the directions the whole bank varies along, and that loads on the
   proliferation programmes is MORE conserved across lines, not less.**
4. **What actually carries it is response GEOMETRY, not programme identity.** Block ablation:
   genericness alone **+0.2352** (the whole model is +0.3174; without it +0.2748), programmes alone
   +0.1576, on-target alone +0.1400, breadth alone +0.0853, magnitude alone +0.0416. Permutation
   importance is led by `gen_pc20` (+0.2430) and `breadth` (+0.1833); `pe_prolif_abs` contributes
   **+0.0059**. Replacing the 9 hallmark sets by random gene sets of matched size costs only
   **+0.023** (model +0.2943 ± 0.0074 over 20 draws), and the |projection| onto a *random* matched
   228-gene direction already reads **+0.0792 ± 0.0119** of the observed +0.1184. **The hallmark
   labels add little; the energy geometry is the signal.** It is not one dataset either: half 1
   **without** its two dominant Orion banks still reads +0.3084, and eight individual banks each
   predict the held half's remainder on their own (+0.15 … +0.41).

**Nothing here is buildable and none is proposed.** Analysis only: no build, no `vcc`, no
`cell-eval2`, no download, nothing under `submissions/` or `data/`, no production file touched, no
file written outside this directory. The spent reserves were not read: neither
`experiments/267-*/work/reserve/` nor JVE127 / MG63 / SLR23 / UMRC3. `uv run python` throughout.
Peak RSS **5.65 GiB** (`build290.py`, `/usr/bin/time -l` 6,065,569,792 B; `ctrl290` 4.84 GiB, all
others ≤ 0.35 GiB) — under the ~8 GiB cap.

**Date 2026-09-06, 14:23–15:10 UTC.** Deterministic, seed **20260906**. `PREREGISTRATION.md` sha256
**`adce9656097f534d4ef43770d77e5a91b8d7ae0f2b2eb9f95f59fe6ede5bfd0f`**, frozen **14:23:29 UTC**
(`work/PREREGISTRATION.sha256`) — before any feature was computed and before any feature–remainder
number existed. Everything labelled POST-HOC was written after the preregistered numbers were read
and cannot change a verdict.

---

## 0. Known answers

| row | required | measured | verdict |
|---|---|---|---|
| **KA-1** 279's `c_held` over all 364 cross-line pairs, rebuilt from the banks with the copied `residual_rows` / `load_unit` / `pair_cos` | max abs diff **0.0e+00** vs `held279.npz`, n 17,106, median +0.018605 | `Ccross` **0.0e+00**, `c_held` **0.0e+00**, n 17,106, median +0.018605; 28 sources / 21 lines / 2,739-gene axis / 364 pairs all matched | **PASS** |
| **KA-2** 284's primary-partition maps and remainders, rebuilt from `Ccross` + `feat281` + the component folds | `C1, C2, r1, r2, ok` at **0.0e+00** vs `rem284.npz`; n 2,848; ceiling +0.3816; remainder +0.3277 | all four at **0.0e+00**, `ok` identical, n **2,848**, ceiling **+0.3816**, remainder **+0.3277**; the partition draw reproduced 20 partitions in **22 tries**, as `rem284.py` did | **PASS** |
| **KA-3** the two declared replication partitions rebuilt the same way | #2 ceiling +0.4622 / remainder +0.3389; #3 +0.5518 / +0.4043 | **+0.4622 / +0.3389** and **+0.5518 / +0.4043** | **PASS** |
| **KA-4** shuffled target, 200 independent shuffles (281's KA-4′ form), ridge | the observed ρ outside the null's spread | null **−0.0070 ± 0.0216**; **0/200** reach +0.3174; z **+15.0** | **PASS** |
| **KA-5** (new, internal) the per-source feature store re-aggregated over the primary halves reproduces `feat290.csv` | max abs diff < 1e-6 and an identical NaN pattern | **2.4e-07** on both halves, NaN patterns identical | **PASS** |

Two further reproductions carried inside the later screens (AGENTS.md:1375-1382, a known-answer row
in every screen): `ctrl290.py` re-derives 284's remainder ρ **+0.3277** from `rem284.npz` before its
first control draw, and `desc290.py`'s per-fold refit of the frozen ridge reproduces the frozen
out-of-fold ρ at **+0.3174**.

KA-5's bar is `< 1e-6`, not `0.0e+00`: the frozen artifact is a CSV of float32 aggregates and the
decimal round-trip carries ~1e-7. 279's own KA-0 / KA-1 use the same bar for the same reason
(`build279.py:131,58`, measured 1.17e-07 / 4.99e-07). It was first written as `0.0e+00`, failed at
2.4e-07, and the bar — not the number — was changed, with the reason recorded in the code.

---

## 1. The design

**The unit.** 284's primary partition. Half 1 = DLD-1, HCT116, HEK293T, HT-29, K562, MCF7, RPE1,
THP-1, iPSC-Tian19, melanoma — **14 sources** (`dld1_lowmoi`, `frangieh_melanoma`,
`gasperini_k562`, `jiang_HT29`, `jiang_K562`, `jiang_MCF7`, `k562_essential`, `k562_gwps`,
`orion_hct116`, `orion_hek293t`, `perturbmatch`, `rpe1`, `tian2019_ipsc`, `wessels_thp1`),
81 cross-line pairs. Half 2 = A549, BxPC-3, Calu-3, H1, HAP1, HEK293, HepG2, Jurkat, TeloHAEC,
cardiomyocyte, iPSC-neuron — **14 sources** (`altair_a549`, `gse295214_h1`, `h1_2025`,
`igvf_cardiomyocyte`, `igvf_telohaec`, `jiang_A549`, `jiang_BXPC3`, `jiang_HAP1`, `myllia_a549`,
`nadig_hepg2`, `nadig_jurkat`, `sunshine_calu3`, `tian_ipsc_neuron`, `xucao_hek293`), 87 pairs.
**No shared source, no shared line.** **2,848 genes** have ≥ 1 cross-line pair in both halves.

**The targets.** `remainder_half2` = C⁽²⁾ minus the 21-feature ridge's out-of-fold prediction of
C⁽²⁾ (primary); C⁽²⁾ raw (secondary); and the symmetric pair.

**The splits.** 284's grouped 5-fold over connected components of the HGNC-family graph, so a test
gene shares no family with its training rows. On the unit the folds are **[2114, 201, 182, 162,
189]** across 437 components and 941 first-family blocks — the giant family component is fold 0, so
74 % of the unit is predicted by a model trained on 734 genes. That unbalance is 284's (§7.3) and is
kept as frozen; §8.1 reports the per-fold reading and a family-blind upper bound.

**A correction to the preregistration's own arithmetic.** PREREG §5 labels the frozen model feature
set "17 columns" and then enumerates 18 (6 genericness/magnitude/breadth + 9 `pe_*` +
`pe_prolif_abs` + `ontarget` + `on_axis`). The **enumerated list is what `lib290.FEATS290` contains
and what every fit used**; the "(17)" is a miscount in the prose and nothing depends on it.

---

## 2. The feature table — formulas as frozen, and coverage

Per source `s`: `blk_s` = the residual rows on the common 2,739-gene axis (each bank centred on its
own mean over the held rows it carries — DEFINITION-276 §2); `nrm_s(g) = ||blk_s(g)||`;
`u_s(g) = blk_s(g)/nrm_s(g)` is the unit row that enters 279's cosine; `mu_s` = the mean raw row the
centring subtracts. Aggregation is the unweighted mean over the half's sources that measure `g`.

| column | formula | coverage on the unit (h1 / h2) |
|---|---|---:|
| `logmag` (a) | `log mean_s [ nrm_s(g) / median_g' nrm_s(g') ]` | 1.000 / 1.000 |
| `breadth` (b) | `mean_s log1p #{ j : |blk_s(g)[j]| > τ_s }`, τ_s = the 99th percentile of \|blk_s\| over the bank's whole held block | 1.000 / 1.000 |
| `gen_cos` (c1) | `mean_s [ u_s(g) · mu_s/||mu_s|| ]` | 1.000 / 1.000 |
| `gen_abscos` (c2) | `mean_s | u_s(g) · mu_s/||mu_s|| |` | 1.000 / 1.000 |
| `gen_pc5` (c3) | `mean_s Σ_{k≤5} ( u_s(g) · v_k^s )²`, `v^s` = the bank's Gram eigenvectors (exact, not randomized); banks with < 200 held rows contribute NaN | 1.000 / 0.999 |
| `gen_pc20` (c4) | the same with `k ≤ 20` | 1.000 / 0.999 |
| `pe_<SET>` (d), 9 columns | `mean_s [ u_s(g) · w_H ]`, `w_H` = the unit indicator of the hallmark set on the axis | 1.000 / 1.000 |
| `pe_prolif_abs` (d-summary) | `mean_s | u_s(g) · w_P |`, `w_P` = the unit indicator of E2F ∪ G2M ∪ MYC_V1 (228 axis genes) | 1.000 / 1.000 |
| `ontarget` (f) | `mean_s u_s(g)[g]` where `g` is on the axis | 0.362 / 0.362 |
| `on_axis` (f-ind) | 1 if `g` is on the 2,739-gene axis | 1.000 / 1.000 |
| *(reported, not in the model)* `pe_stress_abs`, `n_src`, `log_n_src` | stress axis = P53 ∪ UPR ∪ APOPTOSIS (126 axis genes); coverage counts | 1.000 / 1.000 |

Hallmark sets on the axis: E2F 104/200, G2M 102/200, MYC_V1 96/200, P53 44/200, UPR 43/113,
APOPTOSIS 49/161, MTORC1 101/200, OXPHOS 79/200, IFN-γ 38/200. **2,624 of 17,106 held genes (15.3 %)
sit on the 2,739-gene axis**, so `ontarget` is missing for 64 % of the unit by design and the
standardiser imputes the training median.

Note on (c1): the bank's mean *residual* row is identically zero by construction, so the generic
direction is the bank's mean RAW row — the same object before centring. This was stated in the
preregistration before any number existed.

---

## 3. The primary result

Out-of-fold Spearman with `remainder_half2` on the 2,848-gene unit, ridge primary / GBM secondary,
component block bootstrap (2,000 draws) and 281's first-family block bootstrap.

| predictor | ρ, component bootstrap | ρ, first-family bootstrap | R² |
|---|---|---|---:|
| **response features, ridge (frozen cell)** | **+0.3174 [+0.2897, +0.3914]** | **+0.3174 [+0.2737, +0.3563]** | +0.0688 |
| response features, GBM | +0.3112 [+0.2835, +0.4024] | +0.3112 [+0.2692, +0.3499] | +0.0662 |
| B6 `remainder_half1` alone — **the empirical ceiling** (284's statistic) | +0.3277 [+0.2430, +0.3670] | +0.3277 [+0.2892, +0.3674] | −0.3316 |
| B5 `C_half1` alone — within-half conservation | +0.2948 [+0.1989, +0.3527] | +0.2948 [+0.2488, +0.3388] | −1.5001 |
| B3 reliability `r_bar` | +0.1633 [+0.0754, +0.1905] | +0.1633 [+0.1142, +0.2060] | −3.7703 |
| B8 coverage (`log n_src_h1`) | +0.1381 [+0.0486, +0.2083] | +0.1381 [+0.0855, +0.1861] | −460.07 |
| B4 magnitude (`logmag_h1`) | +0.0453 [−0.0596, +0.0657] | +0.0453 [−0.0054, +0.0914] | −17.97 |
| B2 the annotation model's prediction of C⁽²⁾ as a covariate | −0.0177 [−0.0847, +0.0252] | −0.0177 [−0.0616, +0.0226] | −0.6175 |
| **B1 the annotation model (284's 21 features, refit on this target)** | **−0.0566 [−0.1380, −0.0205]** | −0.0566 [−0.1007, −0.0128] | −0.0059 |
| B7 mean predictor | ρ undefined (constant); R² 0 by construction | — | 0 |

**Paired deltas from the same resample** (component / first-family):

| Δρ vs | component | first-family | P(Δ ≤ 0) |
|---|---|---|---:|
| annotation ridge (B1) | **+0.3948 [+0.3442, +0.4896]** | +0.3742 [+0.3172, +0.4289] | 0.000 |
| annotation prediction (B2) | +0.3475 [+0.3080, +0.4305] | +0.3356 [+0.2876, +0.3840] | 0.000 |
| magnitude (B4) | +0.2955 [+0.2511, +0.4047] | +0.2733 [+0.2283, +0.3197] | 0.000 |
| coverage (B8) | +0.1924 [+0.1150, +0.3184] | +0.1807 [+0.1224, +0.2410] | 0.000 |
| reliability (B3) | +0.1722 [+0.1344, +0.2705] | +0.1547 [+0.1100, +0.2012] | 0.000 |
| `C_half1` (B5) | +0.0372 [−0.0062, +0.1436] | +0.0227 [−0.0196, +0.0624] | 0.043 / 0.136 |
| **`remainder_half1` (B6, the ceiling)** | **+0.0048 [−0.0329, +0.1010]** | −0.0111 [−0.0584, +0.0350] | 0.639 / 0.679 |
| GBM | +0.0038 [−0.0423, +0.0403] | +0.0064 [−0.0265, +0.0397] | 0.277 |

**Reading.** The response features reach the ceiling: predicting half 2's remainder from half 1's
*response geometry* is as good as predicting it from half 1's *own measured remainder*, and no
better. Ridge alpha per fold: 3.162 / 0.316 / 0.001 / 0.001 / 0.100; the GBM chose
`lr 0.03, 200 iters, depth 3` in all five folds.

### 3.1 Per-feature Spearman with the primary target (component bootstrap)

| feature | ρ [95 %] | | feature | ρ [95 %] |
|---|---|---|---|---|
| `ontarget` (n 1,032) | **−0.1957 [−0.2693, −0.0748]** | | `pe_G2M_CHECKPOINT` | +0.0684 [+0.0353, +0.1472] |
| `gen_pc20` | **+0.1578 [+0.1152, +0.2346]** | | `pe_UNFOLDED_PROTEIN_RESPONSE` | +0.0649 [+0.0361, +0.1826] |
| `on_axis` | +0.1387 [+0.0539, +0.1745] | | `gen_abscos` | +0.0627 [+0.0445, +0.1714] |
| `n_src` / `log_n_src` | +0.1381 [+0.0461, +0.2042] | | `pe_E2F_TARGETS` | +0.0615 [+0.0159, +0.1233] |
| **`pe_prolif_abs` (P3)** | **+0.1184 [+0.0999, +0.2327]** | | `pe_MTORC1_SIGNALING` | +0.0447 [+0.0226, +0.1597] |
| **`gen_pc5` (P4)** | **+0.1147 [+0.0694, +0.1878]** | | `logmag` | +0.0453 [−0.0558, +0.0623] |
| `breadth` | +0.0869 [+0.0097, +0.1122] | | `pe_INTERFERON_GAMMA_RESPONSE` | +0.0401 [−0.0423, +0.0752] |
| `pe_stress_abs` | +0.0367 [−0.0115, +0.0824] | | `pe_APOPTOSIS` | +0.0351 [−0.0158, +0.0977] |
| `pe_P53_PATHWAY` | −0.0375 [−0.0877, +0.0284] | | `pe_OXIDATIVE_PHOSPHORYLATION` | −0.0399 [−0.1472, −0.0188] |
| `gen_cos` | −0.0368 [−0.0913, +0.0123] | | `pe_MYC_TARGETS_V1` | +0.0207 [−0.0285, +0.0837] |

`ontarget` is the *signed* unit-row entry at the gene's own column and is normally negative
(knockdown): ρ −0.196 means **a deeper relative on-target knockdown goes with a MORE conserved
response**, on the 1,032 unit genes that sit on the axis.

---

## 4. The frozen pass rules

| rule | statistic | frozen bar | measured | verdict |
|---|---|---|---|---|
| **P1** response-derived features predict the held half's remainder | OOF Spearman, ridge, component bootstrap | ≥ +0.10 **and** interval excludes 0 | **+0.3174 [+0.2897, +0.3914]** | **PASS** |
| **P2** they beat the annotation model on the same target | paired Δρ, same resample | interval entirely above 0 | **+0.3948 [+0.3442, +0.4896]** | **PASS** |
| **P3** proliferation / stress engagement → LESS conserved | ρ(`pe_prolif_abs_h1`, `remainder_half2`) | **< 0**, interval excludes 0 | **+0.1184 [+0.0999, +0.2327]** | **FAIL — significant in the OPPOSITE direction** |
| **P4** genericness → LESS conserved | ρ(`gen_pc5_h1`, `remainder_half2`) | **< 0**, interval excludes 0 | **+0.1147 [+0.0694, +0.1878]** | **FAIL — significant in the OPPOSITE direction** |

Secondary readings declared with P3/P4 and reported without re-tuning: `pe_stress_abs`
+0.0367 [−0.0115, +0.0824] (null); `gen_pc20` +0.1578, `gen_abscos` +0.0627, `gen_cos`
−0.0368 [−0.0913, +0.0123]. Of the nine individual `pe_*` columns, **one** is significantly
negative (`OXIDATIVE_PHOSPHORYLATION`, −0.0399 [−0.1472, −0.0188]) and **four** are significantly
positive (`E2F_TARGETS`, `G2M_CHECKPOINT`, `UNFOLDED_PROTEIN_RESPONSE`, `MTORC1_SIGNALING`);
`P53_PATHWAY`, `MYC_TARGETS_V1`, `APOPTOSIS` and `INTERFERON_GAMMA_RESPONSE` are null.

**The hypothesis the mandate froze — "a response is conserved when it engages gene-specific
downstream programmes, and state-dependent when it engages the recipient's proliferation / stress
axis" — is refuted on this map.** The association runs the other way.

---

## 5. The three controls

| control | result |
|---|---|
| **C1 — genes shuffled** (200 independent shuffles of the target, ridge, 281 KA-4′ form) | null **−0.0070 ± 0.0216**; **0/200** draws reach +0.3174; **z +15.0** |
| **C2 — features on the SAME half as the target** (the leakage reading, labelled) | **+0.4878 [+0.4741, +0.5666]** vs the honest +0.3174; paired Δ **+0.1747 [+0.1259, +0.2449]**. A careless design would have bought +0.17 of Spearman for free |
| **C3 — programmes replaced by random matched gene sets** (20 draws) | model **+0.2943 ± 0.0074 [+0.2814, +0.3063]**; observed +0.3174; **0/20** reach it. The 9 real hallmark sets are worth **+0.023**. And the \|projection\| onto a *random* 228-gene direction already reads **+0.0792 ± 0.0119** against `pe_prolif_abs`'s +0.1184 — **two thirds of the "programme engagement" signal is generic energy geometry, not programme identity** |

---

## 6. Symmetric direction and the secondary target

| cell | response ridge | response GBM | annotation | the other half's own map | `r_bar` |
|---|---|---|---|---|---|
| **primary** h1 features → `remainder_half2` | **+0.3174 [+0.2897, +0.3914]** | +0.3112 | −0.0566 | C⁽¹⁾ +0.2948 · r⁽¹⁾ +0.3277 | +0.1633 |
| **symmetric** h2 features → `remainder_half1` | **+0.3411 [+0.2885, +0.3916]** | +0.3612 | +0.1046 [+0.0768, +0.2442] | C⁽²⁾ +0.3321 · r⁽²⁾ +0.3277 | **+0.5143** |
| secondary h1 features → C⁽²⁾ raw | **+0.3878 [+0.3658, +0.4628]** | +0.3865 | +0.2364 [+0.1691, +0.2740] | C⁽¹⁾ +0.3816 | +0.2381 |
| symmetric secondary h2 features → C⁽¹⁾ raw | **+0.5133 [+0.4371, +0.5309]** | +0.5403 | — | C⁽²⁾ +0.3816 | — |

Paired: response − annotation **+0.2155 [+0.0890, +0.2717]** (symmetric) and **+0.1659 [+0.1277,
+0.2537]** (secondary, raw C⁽²⁾). On the raw conservation map the response features **beat half 1's
own conservation map** at predicting half 2's (+0.3878 vs +0.3816, paired +0.0300 [−0.0145,
+0.1445], i.e. matched), and beat it clearly in the symmetric direction (+0.5133 vs +0.3816).

**One baseline must be read with care, and it is the one that wins.** In the symmetric direction
`r_bar` alone reads **+0.5143** and beats the response model (paired −0.1446 [−0.1935, −0.0274]).
`r_bar` is **not** a half-independent quantity: it is the median split-half residual cosine over
exactly four banks — `k562_gwps`, `orion_hct116`, `orion_hek293t`, `dld1_lowmoi`
(`rel279.py:44-49,222`) — and **all four are half-1 sources**, three of them the largest banks in
the study. Predicting `remainder_half1` from `r_bar` therefore shares measurement noise with the
target. In the primary direction, where `r_bar` is an independent half-1 measurement predicting a
half-2 target, it reads only +0.1633 and the response model beats it by +0.1722. **This also
corrects R284 §5.2's parenthetical that "`r_bar` and the magnitude are measured on all 28
sources": the magnitude is (279's `mag`, a median over every source), `r_bar` is measured on four
banks, all in half 1.** It does not change any 284 verdict — 284 subtracts `r_bar` from both halves
symmetrically — but the sentence as written is not accurate.

---

## 7. Replication

**The two declared partitions** (#2 and #3, named in the preregistration before any number was
read; their halves mix the primary partition's source halves, so their feature blocks were
re-aggregated from the per-source store — KA-5).

| partition | n | direction | response | annotation | paired Δ | the other half's remainder | `pe_prolif_abs` | `gen_pc5` |
|---|---:|---|---|---:|---|---:|---|---|
| **#2** (ceiling +0.4622, remainder +0.3389) | 5,623 | h1→r2 | +0.2443 [+0.1856, +0.2725] | −0.0838 | +0.3434 [+0.2970, +0.4555] | +0.3389 | +0.1088 [+0.0626, +0.1624] | +0.1321 [+0.0542, +0.1679] |
| | | h2→r1 | +0.2552 [+0.1810, +0.2722] | −0.0550 | +0.3078 [+0.2532, +0.3553] | +0.3389 | −0.0104 [−0.0792, +0.0070] | +0.0043 [−0.0612, +0.0245] |
| **#3** (ceiling +0.5518, remainder +0.4043) | 4,667 | h1→r2 | +0.2469 [+0.2323, +0.3298] | −0.0861 | +0.3346 [+0.2425, +0.4488] | +0.4043 | −0.0081 [−0.0293, +0.1001] | +0.0351 [+0.0170, +0.1314] |
| | | h2→r1 | +0.2371 [+0.1462, +0.2757] | −0.0823 | +0.3475 [+0.2763, +0.5214] | +0.4043 | +0.0959 [+0.0603, +0.1772] | +0.1285 [+0.0494, +0.1805] |

**P1 and P2 replicate on all four cells. P3 and P4 do not reverse:** `pe_prolif_abs` and `gen_pc5`
are positive on 2/4 and 3/4 and never significantly negative.

**POST-HOC sweep over all 20 admissible partitions × both directions (40 cells, point estimates):**

| quantity | median | range | cells |
|---|---|---|---|
| response ridge | **+0.2479** | [+0.0891, +0.3700] | **40/40 > 0**, 39/40 ≥ +0.10 |
| annotation ridge | −0.0735 | [−0.2474, +0.1417] | 31/40 < 0 |
| paired response − annotation | — | — | **40/40 > 0** |
| the other half's own remainder (284's statistic) | +0.3191 | [+0.0981, +0.4641] | response higher on 11/40 |
| `logmag` alone | +0.1174 | [−0.0118, +0.2711] | 38/40 > 0 |
| `pe_prolif_abs` | +0.0420 | [−0.0380, +0.1512] | 30/40 > 0, **10/40 < 0** |
| `gen_pc5` | +0.0635 | [−0.0443, +0.2280] | 36/40 > 0, **4/40 < 0** |

The weakest cell is partition #4 (the one whose halves overlap on only 1,345 genes; 284 flagged it
for the same reason).

---

## 8. What actually carries the remainder (POST-HOC, labelled)

### 8.1 Fold-by-fold, and the family-held cost

| fold | n (train) | response | annotation | `remainder_half1` |
|---|---|---|---|---|
| 0 (the giant family component) | 2,114 (734) | +0.3096 | −0.0773 | +0.3373 |
| 1 | 201 (2,647) | +0.3875 | −0.1331 | +0.3945 |
| 2 | 182 (2,666) | +0.3474 | −0.1456 | +0.3834 |
| 3 | 162 (2,686) | +0.3766 | +0.0131 | +0.3499 |
| 4 | 189 (2,659) | +0.2733 | −0.0116 | +0.1268 |

A family-blind random 5-fold (families NOT held; an upper bound, not a verdict) reads **+0.3455** —
so holding whole family components costs about **0.03**, and the frozen fold unbalance is not what
produces the result.

### 8.2 Attribution

**Permutation importance** (drop in Spearman when one column is permuted within its held fold, 20
repeats, refit base +0.3174):

`gen_pc20` **+0.2430 ± 0.0126** · `breadth` +0.1833 ± 0.0112 · `logmag` +0.1037 ± 0.0081 ·
`pe_G2M_CHECKPOINT` +0.0857 · `gen_pc5` +0.0848 · `pe_P53_PATHWAY` +0.0846 · `pe_MYC_TARGETS_V1`
+0.0808 · `on_axis` +0.0400 · `pe_APOPTOSIS` +0.0356 · `ontarget` +0.0343 · `pe_UPR` +0.0312 ·
`pe_OXPHOS` +0.0235 · `pe_E2F` +0.0126 · `pe_MTORC1` +0.0063 · **`pe_prolif_abs` +0.0059** ·
`gen_abscos` +0.0028 · `pe_IFN-γ` −0.0011 · `gen_cos` −0.0027.

**Block ablation** (each block alone, and the model without it):

| block | alone | model without it |
|---|---:|---:|
| (c) genericness — 4 columns | **+0.2352** | **+0.2748** |
| (d) programmes — 10 columns | +0.1576 | +0.2941 |
| (f) on-target — 2 columns | +0.1400 | +0.2664 |
| (b) breadth — 1 column | +0.0853 | +0.3130 |
| (a) magnitude — 1 column | +0.0416 | +0.3162 |

### 8.3 It is not magnitude, reliability, coverage, or a restatement of C⁽¹⁾

| covariate partialled out of BOTH sides (rank residuals) | response model | `remainder_half1` |
|---|---:|---:|
| nothing | +0.3174 | +0.3277 |
| `logmag_h1` | +0.3162 | +0.3357 |
| `n_pairs_half2` | +0.3244 | +0.3264 |
| `r_bar` | +0.2779 | +0.2901 |
| `logmag_h1` + `r_bar` | +0.2832 | +0.2970 |

Magnitude-matched terciles of `logmag_h1`: +0.2916 (n 950) / +0.3665 (949) / +0.3027 (949) —
present at every magnitude. **The precise claim, since two readings of magnitude differ:** its
MARGINAL correlation with the primary target is +0.0453 and not distinguishable from 0 (and its
block alone reads +0.0416), while its CONDITIONAL contribution inside the fitted model is real
(permutation drop +0.1037) — magnitude is a useful adjustment to the other columns, not a predictor
on its own. On the 40-cell sweep its marginal reading is larger (median +0.1174, 38/40 > 0), so
"magnitude carries nothing" would be too strong; "magnitude does not explain the +0.3174" is what
the partials and terciles support. Adding coverage as a 19th column changes nothing (+0.3137, paired
−0.0033 [−0.0066, +0.0024]).

**And it is not the same information half 1's own remainder carries.** ρ(prediction, r⁽¹⁾) = +0.393,
ρ(prediction, C⁽¹⁾) = +0.472; partialling r⁽¹⁾ out of both sides still leaves
ρ(r⁽²⁾, prediction | r⁽¹⁾) = **+0.2173**, and symmetrically ρ(r⁽²⁾, r⁽¹⁾ | prediction) = **+0.2319**.
Adding r⁽¹⁾ as a 19th column lifts the model to **+0.3834** (paired +0.0580 [−0.0081, +0.0870]).
**The response geometry and the measured remainder are two partly independent views of the same
gene property.**

### 8.4 The 20 genes at each end

Most **under**-predicted (the actual remainder far above the model): `HMGB2`, `NFX1`, `ZNF286A`,
`ZNF92`, `ZNF426`, `TCF3`, `CENPB`, `MTF2`, `CAMTA1`, `CREM`, `EIF3CL`, `CNBP`, `HDGF`, `SMYD3`,
`GPBP1`, `TET1`, `TSC22D4`, `ZNF146`, `TGIF2`, `ZNF644` — **overwhelmingly sequence-specific DNA
binders (four C2H2 zinc fingers, a homeobox-adjacent TGIF2, TCF3, CAMTA1, CREM, MTF2, TET1)**, i.e.
exactly the cell-identity families 279 placed at the state-dependent end of the axis. The model
under-predicts how conserved *their* responses are.

Most **over**-predicted: `ARCN1`, `DDX27`, `PSMA4`, `SKP1`, `TMEM214`, `KDM4A`, `BUB1`, `TMED2`,
`PRORP`, `ATP5MF`, `PSMA1`, `ZMYM3`, `RHEB`, `WDR24`, `FTL`, `PSMA3`, `COPB2`, `POLR2H`, `ISG20L2`,
`PSMD2` — **proteasome (PSMA1/3/4, PSMD2), COPI (ARCN1, COPB2, TMED2), mitotic and growth machinery
(BUB1, RHEB, WDR24)**: large, broad, generic responses that the model reads as conserved and that
in fact are not.

### 8.5 It is not the two Orion banks (POST-HOC)

Half 1's 14 sources are wildly unequal — `orion_hct116` (16,594 held rows) and `orion_hek293t`
(16,721) against 23 for `tian2019_ipsc` — so the aggregate could be a property of Orion rather than
of perturbation responses in general. Re-aggregating the same 18 columns over source subsets (the
full half-1 subset reproduces the frozen cell at **+0.317385**, an internal known-answer row):

| half-1 subset | sources | unit genes covered | model ρ | `gen_pc20` | `breadth` |
|---|---:|---:|---|---:|---:|
| **all 14 (frozen)** | 14 | 2,848 | **+0.3174 [+0.2931, +0.3943]**† | +0.1578 | +0.0869 |
| **without the 2 Orion banks** | 12 | 2,818 | **+0.3084 [+0.2767, +0.3781]** | +0.1940 | +0.1013 |
| the 2 Orion banks only | 2 | 2,690 | +0.2537 [+0.1832, +0.2977] | +0.0196 | +0.0441 |
| the 7 banks with ≥ 1,000 rows | 7 | 2,848 | +0.3018 [+0.2581, +0.3615] | +0.1528 | +0.0857 |
| the 7 banks with < 1,000 rows | 7 | 314 | +0.2289 [+0.1009, +0.3852] | +0.0673 | −0.0286 |

**Single sources** (each one bank, predicting a remainder measured on 14 disjoint sources in 11
other lines): `jiang_K562` +0.4074 (n 206) · `jiang_HT29` +0.3130 (206) · `rpe1` +0.2768 (1,950) ·
`orion_hct116` +0.2581 (2,655) · `jiang_MCF7` +0.2438 (206) · `perturbmatch` +0.2440 (2,658) ·
`k562_gwps` +0.2383 (2,751) · `orion_hek293t` +0.1885 (2,658) · `k562_essential` +0.1473 (1,754) ·
`dld1_lowmoi` +0.0155 (725, the one null). Four half-1 banks cover < 200 unit genes and were
skipped.

† the point estimate is the frozen +0.317385 to six decimals; the interval differs from §3's
[+0.2897, +0.3914] only because this table's bootstrap uses a different seed (a 2,000-draw
Monte-Carlo difference, not a different estimate). The intervals in this table are therefore
comparable with each other, and §3's is the one quoted everywhere else.

**A single bank's response geometry already predicts the held half's remainder**, and dropping the
two banks that dominate coverage does not weaken the model — it strengthens the genericness columns
(`gen_pc20` +0.1940 without Orion vs +0.1578 with it). The result is a property of perturbation
responses, not of one dataset.

### 8.6 284's two described tails

| feature (half-1) | 284's ABOVE tail (n 20) | 284's BELOW tail (n 20) | rest | MWU p |
|---|---:|---:|---:|---:|
| `logmag` | **+0.473** | +0.115 | +0.186 | 0.0001 |
| `breadth` | **+3.863** | +1.904 | +2.333 | < 0.0001 |
| `gen_pc20` | **+0.407** | +0.282 | +0.279 | 0.0071 |
| `gen_pc5` | +0.274 | +0.224 | +0.206 | 0.209 |
| `pe_MYC_TARGETS_V1` | −0.056 | +0.004 | −0.010 | 0.0015 |
| `pe_prolif_abs` | +0.115 | +0.096 | +0.090 | 0.060 |
| `pe_APOPTOSIS` | +0.022 | +0.004 | +0.004 | 0.017 |
| `pe_OXIDATIVE_PHOSPHORYLATION` | −0.013 | +0.009 | −0.001 | 0.016 |
| `pe_INTERFERON_GAMMA_RESPONSE` | +0.017 | −0.002 | +0.001 | 0.014 |

**284's mRNA-surveillance tail is separated from its proteostasis / condensin tail mainly by size,
breadth and PC-energy — not by which hallmark programme the response engages.** The one programme
that does separate them is MYC targets, and it points the *other* way (the conserved tail loads
*less* on MYC). Over the 40 tail genes the model orders them at ρ +0.5677 (predicted +0.0250 above
vs −0.0139 below, against actual +0.2605 vs −0.0941): the ordering is right, the magnitude is
compressed ~10×.

---

## 9. What is new against 284 and 281 (programme rule 10)

* **NEW.** 284 showed the remainder reproduces across disjoint source sets and that 21 annotation
  features do not carry it, and closed with *"getting at it needs different features — response-
  derived or pathway-level ones"*. **This is that test, and response-derived features do carry it,
  at the ceiling** (ρ +0.3174 vs the +0.3277 that half 1's own remainder achieves; annotation
  −0.0566). 284 could not know whether the remainder was reachable from anything except itself.
* **NEW.** The mechanism is **response geometry, not programme identity**: genericness (top-5 /
  top-20 PC energy) and breadth carry it; random gene sets of matched size recover 93 % of the
  model (+0.2943 of +0.3174) and two thirds of the "proliferation engagement" statistic. **A
  pathway-level feature set is not what is needed; a spectral one is.** And it is a property of
  perturbation responses rather than of any one bank: half 1 minus its two Orion banks reads
  +0.3084, and a single bank of 206–2,751 measured genes already reaches +0.15 … +0.41 (§8.5).
* **NEW, and the sign is the finding.** More generic and more proliferation-loaded responses are
  **more** conserved across lines — the opposite of what R284 §6's described tail suggested and of
  what this experiment preregistered. It fits 279's axis read from the other side: core-machinery
  perturbations produce large, broad, spectrally central responses and transfer; cell-identity
  perturbations produce small, narrow, idiosyncratic ones and transfer less. **With one systematic
  exception the geometry rule gets wrong** (§8.4): the 20 most under-predicted genes are
  overwhelmingly sequence-specific DNA binders (four C2H2 zinc fingers, `TCF3`, `CAMTA1`, `CREM`,
  `MTF2`, `TET1`, `TGIF2`) — their responses are narrow and non-generic, the model calls them
  state-dependent, and they are more conserved than that. Whatever carries the last of the
  remainder is not geometry.
* **EXTENDED.** 281's `r_bar` is a four-bank quantity (`k562_gwps`, `orion_hct116`,
  `orion_hek293t`, `dld1_lowmoi`), all in half 1 — which is why it reads +0.16 predicting half 2 and
  +0.51 predicting half 1. R284 §5.2's parenthetical ("measured on all 28 sources") is inaccurate
  for `r_bar`; no 284 verdict changes.
* **EXTENDED.** 284's per-gene remainder is not a single quantity: the response model and the
  measured remainder each add ≈ +0.22 of rank information beyond the other, and together reach
  +0.3834.
* **CONTRADICTED (a frozen prediction of this experiment, not of 284).** P3 and P4 both fail with
  the opposite sign. Recorded as a negative; nothing was re-tuned after they failed.

**Competition track: nothing.** No arm, no dose, no rendering, no model use is proposed. The R² is
0.069 on a rank statistic whose ceiling is itself +0.38, 277 already closed the one action a
conservation-keyed rendering could take, and none of these features is available at build time for
a challenge target in an unmeasured line.

---

## 10. Limitations, stated plainly

1. **Both halves come from the same 28 banks and the same 2,739-gene axis.** The line-group
   partition removes shared sources and shared lines, not shared pipelines. A processing artefact
   common to, say, all Replogle-derived banks would appear in both halves (284 §7.4, unchanged).
2. **"Genericness" is defined against each bank's own row population.** A gene's own row
   contributes to its bank's PCs and to `mu_s`; with 200–16,721 rows per bank the self-contribution
   is small but not zero, and it is a within-half quantity either way, so it cannot see the target.
3. **The unit is 2,848 genes and the fold structure is unbalanced** (74 % in one family component,
   trained on 734 genes). §8.1 shows the result is present in every fold and that a family-blind
   split gains only 0.03.
4. **`ontarget` is missing for 64 % of the unit** and is median-imputed. Its −0.196 correlation is
   measured on the 1,032 on-axis genes only.
5. **The breadth threshold is a bank-level quantile, not a DE test.** The compacted banks carry a
   control *mean* (`b_ctrl`) and no control variance, so a real DE count is not computable from
   them; τ_s fixes the expected count at 1 % of the axis per bank, which makes `breadth` a relative,
   not an absolute, DE count.
6. **The remainder is defined by 284's annotation model.** "The features do not carry it, the
   responses do" is a statement about those 21 features and this map.
7. **Multiplicity.** Verdict cells: 4 (P1–P4). Reported in total: ~140 cells across 2 targets ×
   2 directions × 2 models × 2 bootstraps, 8 baselines, 3 controls, 20 partitions × 2 directions,
   21 per-feature correlations, 18 permutation-importance rows, 5 block ablations, 5 partials,
   3 strata, 2 tails, 15 source-subset cells — all but the four labelled secondary, comparative or post-hoc.
8. **One map, one statistic.** Everything is 279's cross-line residual cosine on its 28 banks.

---

## Files

| path | what |
|---|---|
| `PREREGISTRATION.md` / `work/PREREGISTRATION.sha256` | frozen 14:23:29 UTC, `adce9656…`; every script aborts on a digest mismatch |
| `lib290.py` | copied loaders / metrics / models / bootstraps (provenance line-cited in its docstring) |
| `build290.py` → `work/feat290.csv`, `work/persrc290.npz`, `work/build290.json` | KA-1 and the response-derived feature block; per-gene half-1 and half-2 columns, and the per-source store |
| `fit290.py` → `work/pred290.csv`, `work/fit290.json` | KA-2 / KA-3, the preregistered fits, baselines, controls C1 and C2, P1–P4 |
| `ctrl290.py` → `work/ctrl290.json` | control C3, random matched programmes, 20 draws |
| `rep290.py` → `work/rep290.json` | KA-5, the two declared replication partitions, the 20-partition sweep |
| `desc290.py` → `work/desc290.json` | POST-HOC: per-fold, permutation importance, block ablations, partials, tails, best/worst genes |
| `srcabl290.py` → `work/srcabl290.json` | POST-HOC: the half-1 source-subset ablation (§8.5), with the frozen cell as its known-answer row |
| `work/*.log`, `work/*.nohup`, `work/rss_build.txt` | run logs and `/usr/bin/time -l` peak RSS |
