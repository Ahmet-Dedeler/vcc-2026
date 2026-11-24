# 291 — F13 (response geometry predicts cross-line conservation) on two surfaces OUTSIDE the 28-bank map

**C3 (Broad-16, disjoint lines, Cas9): FAIL (ρ +0.2099 [−0.0130, +0.4263]; pc20 +0.2952
[+0.0695, +0.5065]; annotation +0.2996 [+0.0826, +0.4894]; C_held +0.2104 [−0.0043, +0.4144]);
C4 (GWCD4i CD4-vs-K562, Stim48hr): FAIL (ρ +0.2032 [+0.1296, +0.2517]; partial given magnitude
+0.0610 [+0.0296, +0.1563]; annotation +0.2702 [+0.1922, +0.3111]; pc20 +0.1858 [+0.1189, +0.2324])**

Both claims FAIL, and they fail on **different halves of the same rule for the same reason**. What
survives, measured on two surfaces the 28-bank map never touched:

1. **F13's transfer is real but MUCH weaker off-map, and on Broad-16 it is not separated from
   zero.** The 290 ridge fitted on the map's own conservation reads +0.2032 [+0.1296, +0.2517] on
   the GWCD4i authors' published CD4-vs-K562 correlation (1,773 genes; 0/200 gene shuffles reach
   it, z +8.4) and +0.2099 [−0.0130, +0.4263] on Broad-16's 12-line conservation (82 targets;
   10/200 target shuffles reach it, z +1.8). On the map itself the same model reads +0.383
   (grouped-CV Spearman against `C_held`). **Off-map, roughly half the rank information is gone;
   on 82 targets what remains does not clear its own interval.**
2. **F13's headline — response geometry BEATS annotation — does not transfer at all.** Inside the
   map the response ridge beat the annotation ridge by **+0.3948** (R290 §3). On both surfaces here
   the ordering **reverses**: Broad-16 annotation +0.2996 vs response +0.2099 (paired −0.0837
   [−0.3570, +0.2062]); GWCD4i Stim48hr annotation +0.2702 vs response +0.2032 (paired −0.0641
   [−0.1236, +0.0112], P(Δ ≤ 0) = 0.964); the same sign in all three CD4 conditions and on 268's
   shared-energy fraction (−0.2391 [−0.5117, +0.0363]). **290's +0.39 was a property of its target
   — the remainder AFTER the annotation model — not a property of conservation.** C4's condition 2
   is the only one that fails, in every condition.
3. **The inverted-sign finding — broad, spectrally central responses are MORE conserved — DOES
   transfer, on both surfaces, as a single feature.** `gen_pc20` reads **+0.2952 [+0.0695, +0.5065]**
   on Broad-16 and **+0.1858 [+0.1189, +0.2324]** on GWCD4i-Stim48hr; `breadth` +0.1857 / +0.1275
   [+0.0160, +0.1555]. On Broad-16 `gen_pc20` ALONE beats the 18-column ridge that contains it.
   The one preregistered condition that passes on both surfaces is the sign of the geometry.
4. **Against the authors' own magnitude the effect shrinks but stays.** Partialling out
   `n_degs_MASH_K562` and `n_degs_MASH_Stim48hr` leaves +0.0610 [+0.0296, +0.1563] for the model and
   +0.0720 [+0.0308, +0.1406] for `gen_pc20` — the K562 DE count alone reads +0.2959 and the
   condition's own DE count +0.5376, so **most of the published conservation statistic is
   magnitude**, and about a quarter of the geometry signal survives it.
5. **The annotation and response models are two partly independent views, on both surfaces.**
   ρ(M_resp, outcome | M_ann) = +0.1621 (Broad-16) and +0.1113 (GWCD4i); ρ(M_ann, outcome | M_resp)
   = +0.2699 / +0.2122; the two predictions correlate only +0.194 / +0.386. But the combined 39-column
   model does **not** beat annotation alone (GWCD4i-Stim48hr +0.2166, paired −0.0514 [−0.1060,
   +0.0170]).
6. **290's DNA-binder exception does NOT reproduce on the CD4 surface, and inverts.** 290 §8.4 found
   the most under-predicted genes were sequence-specific DNA binders. Here ZNFs sit at the model
   (mean rank-residual −0.0039, p = 0.44) and annotated TFs sit **below** it (−0.0550, i.e. LESS
   conserved between CD4 and K562 than their geometry predicts).

**Nothing is buildable and none is proposed.** Analysis only: no build, no `vcc`, no `cell-eval2`,
no download, nothing under `submissions/` or `data/`, no production file touched, no file written
outside this directory except the two claim entries appended to `docs/CONFIRMATION-LEDGER.md`. The
spent Broad-16 reserve (`JVE127`, `MG63`, `SLR23`, `UMRC3`) was not read — `load268.py`'s assertion
is kept and fired on every load. `uv run python` throughout. Peak RSS **1.69 GiB** (`surfA291.py`,
`/usr/bin/time -l` 1,811,349,504 B; `surfB291.py` 0.45 GiB, others ≤ 0.3 GiB) — far under the ~8 GiB cap.

**Date 2026-09-06, 15:07–15:30 UTC.** Deterministic, seed **20260906**. `PREREGISTRATION.md` sha256
**`3a4f5f8534469fba203268ab1a0d96afb2d9ce2e310ba898156ccb9edffa78da`**, frozen **15:13:32 UTC**
(`work/PREREGISTRATION.sha256`), and the two claims were appended to `docs/CONFIRMATION-LEDGER.md`
before any feature–outcome number on either surface existed. Every script aborts on a digest
mismatch. Everything labelled POST-HOC was written after the two verdicts were read.

---

## 0. Known answers

| row | required | measured | verdict |
|---|---|---|---|
| **KA-1** 290's primary cell rebuilt from the per-source store (`persrc290.npz` → half-1 aggregation → the frozen out-of-fold ridge against `remainder_half2`) | ρ = **+0.317385** on 2,848 genes, folds [2114, 201, 182, 162, 189] | ρ **+0.31738512** (six decimals **+0.317385**, diff **1.2e-07**), n **2,848**, folds **[2114, 201, 182, 162, 189]**, per-fold α 3.162 / 0.316 / 0.001 / 0.001 / 0.100 — R290 §3's list | **PASS** |
| **KA-2** 279's `c_held`, READ from `held279.npz` (not recomputed) | n 17,106, median +0.018605, equal to `znf279.csv`'s `C_held` | n **17,106**, median **+0.0186052**, max abs diff vs `znf279.C_held` **5.0e-07** (the CSV's decimal round-trip; 279's own KA bar) | **PASS** |
| **KA-3** 287's GWCD4i row — the random controls at chance before anything else on surface B | \|median `random_r*`\| < 0.01 and real median in (0.25, 0.36), all three conditions | Rest real **+0.3030** / random **+0.00396** (p95 +0.0546); Stim8hr **+0.3147** / **+0.00365** (+0.0550); Stim48hr **+0.3117** / **+0.00382** (+0.0555); n = 1,880 rows per condition | **PASS** |
| **KA-4** 268's per-target decomposition recomputed from `work/halves/<line>.npz` on the 12 working lines | `shared_full`/`line_full`/`noise_full` reproduce `per_target_decomp_cpm5.csv` (82 targets); `decomp268_cpm5.json`'s `cell_halves.pooled_full_depth` = **0.090681 / 0.115379 / 0.793940** | max abs diff **5.5e-08 / 2.5e-08 / 8.2e-08** (the CSV is float32); pooled **0.090681 / 0.115379 / 0.793940**, max abs diff **3.5e-08** — the row reproduced | **PASS** |
| **KA-5** (internal) the aggregation over an arbitrary source subset reproduces `feat290.csv` | < 1e-6, identical NaN pattern, both halves | **2.38e-07** (h1) / **2.37e-07** (h2), NaN patterns identical | **PASS** |

KA-1's and KA-5's bar is a float32/decimal round-trip (290's own KA-5 records the same reason);
KA-4's is `per_target_decomp_cpm5.csv` being a float32 CSV. **`per_target_decomp.csv` and
`decomp268.json` without the `_cpm5` suffix are an EIGHT-line partial run** (`genes_kept` 12,003,
`lines` 8) — the mandate's file pointer resolves to the 12-line objects `per_target_decomp_cpm5.csv`
/ `decomp268_cpm5.json` (12,199 genes, 12 lines), which is what R268 §1–§2 reports and what is used
here.

---

## 1. The two surfaces, and what "outside the map" means

| | surface A — Broad-16 | surface B — GWCD4i CD4-vs-K562 |
|---|---|---|
| source | DepMap/Broad 16-line Perturb-seq, 12 working lines, **Cas9 knockout** | the GWCD4i authors' own `K562_comparison.suppl_table.csv` |
| outcome | mean over 66 line pairs of cos(line-centred residual rows) — DEFINITION-276 §3 on 268's full responses, `min_cpm` 5.0, **12,199 genes** | `logfc_pearson_r`: their Pearson correlation of DESeq2 LFCs between primary CD4 T cells and K562 |
| unit | **82** non-OR targets available on all 12 lines and in the held map (+7 OR controls as a secondary; `KIF11`, `PSMA1`, `RRM1` unavailable; 5 non-OR available targets are outside the held map) | **1,773** genes per condition (of 1,880 in the table), 408 family components |
| overlap with the 28 banks | **none** — no Broad-16 line and no Broad-16 dataset is among the 28 | CD4 is in **none** of the 28; the K562 side is excluded from the features by construction |
| features from | **all 28 banks** | the **23 K562-free banks** (excluding `gasperini_k562`, `jiang_K562`, `k562_essential`, `k562_gwps`, `perturbmatch`) |
| perturbation modality | Cas9 KO, 2 guides/target, one guide per cell | CRISPRi in primary human T cells, 4 donors × 3 states |

**Neither surface is untouched and neither is a reserve** (268 / 283 / 289 analysed Broad-16; 287
read this exact GWCD4i table against ZNF labels, donor reproducibility and MASH DE counts). What is
new is that neither has ever been compared with 290's response-geometry features or with any model
fitted on the 28-bank map, and both claims were frozen before that comparison. This is an out-of-map
**generalisation** test, not an untouched confirmation.

**The model.** 290's 18 response-geometry columns aggregated over the stated bank subset; ridge
fitted to predict 279's `C_held` on all 17,106 held genes, α by grouped CV over 284's
family-component folds, then applied. Fitted quality on the map itself:
**M_resp α 31.6, grouped-CV ρ +0.3833** (23-bank version +0.3692); **M_ann α 3,162, grouped-CV
ρ +0.2263**. No surface outcome enters any fit. The leak-free out-of-fold variant is reported in
every table and changes nothing (Broad-16 +0.2152; GWCD4i-Stim48hr +0.1920).

---

## 2. C3 — Broad-16

`C_broad16` on the 82-target unit: median **+0.02957**, IQR [+0.01585, +0.05701], range
[+0.00200, +0.37925], **0 targets ≤ 0**, every target on all 66 pairs.

### 2.1 The primary outcome (paired target bootstrap, 2,000 draws)

| predictor | ρ [95 %] |
|---|---|
| **M_resp — the 290 ridge, 28 banks (the frozen cell)** | **+0.2099 [−0.0130, +0.4263]** |
| **M_ann — 284's annotation ridge (baseline)** | **+0.2996 [+0.0826, +0.4894]** |
| `C_held` — the map's own conservation (the reference) | +0.2104 [−0.0043, +0.4144] |
| **`gen_pc20` (the frozen single feature)** | **+0.2952 [+0.0695, +0.5065]** |
| `breadth` | +0.1857 [−0.0439, +0.3852] |
| `logmag` (magnitude alone) | +0.1518 [−0.0685, +0.3639] |
| M_resp, out-of-fold variant | +0.2152 [−0.0058, +0.4324] |
| M_ann, out-of-fold variant | +0.3000 [+0.0813, +0.4927] |

Paired, same resample: M_resp − M_ann **−0.0837 [−0.3570, +0.2062]** (P(Δ ≤ 0) 0.715);
M_resp − `C_held` +0.0021 [−0.1700, +0.1647]; M_resp − `logmag` +0.0611 [−0.1983, +0.3124].

**Null**: 200 target shuffles read **+0.0157 ± 0.1068**; **10/200** reach +0.2099; z **+1.8**.

**Secondary unit** (the 82 + the 7 held-map OR controls, n 89): M_resp +0.2083 [−0.0041, +0.4195],
M_ann +0.2850 [+0.0836, +0.4636], `C_held` +0.1537 [−0.0713, +0.3621], `gen_pc20` +0.2597
[+0.0476, +0.4615]. No verdict changes.

### 2.2 The second outcome — 268's per-target shared-energy fraction (n 77)

| predictor | ρ [95 %] |
|---|---|
| M_ann | **+0.3528 [+0.1347, +0.5474]** |
| `logmag` | +0.2766 [+0.0499, +0.4674] |
| `breadth` | +0.2617 [+0.0338, +0.4642] |
| `gen_pc20` | +0.1976 [−0.0152, +0.4151] |
| M_resp | +0.1118 [−0.1209, +0.3190] |
| `C_held` | +0.0878 [−0.1393, +0.3091] |

Paired M_resp − M_ann **−0.2391 [−0.5117, +0.0363]** (P(Δ ≤ 0) 0.955). The shared-energy fraction is
a *variance* statistic, not a cosine, and it is the one place where magnitude and breadth read
clearly positive while the fitted ridge does not.

### 2.3 Depth (programme rule 9 — reported, never in the model)

| statistic | value |
|---|---|
| ρ(`C_broad16`, target's median cells/KO over 12 lines) | **−0.081** |
| ρ(`C_broad16`, target's minimum cells/KO) | −0.035 |
| ρ(pair mean cosine, pair geometric-mean UMI/KO), over the **66 line pairs** | **+0.400** |
| C3's ρ after partialling log median cells/KO out of both sides | +0.2106 |
| `gen_pc20`'s ρ after the same partial | +0.2927 |

**Depth acts on the LINE, not on the target.** Deeper line pairs agree more (+0.400 over 66 pairs —
the Broad-16 paper's own UMI–power law, and the mechanism 283 found behind C2a's miss), while
cells-per-KO carries nothing at the target level and partialling it out moves C3 by 0.0007. The
frozen rule is unconditioned, as declared; the depth-adjusted reading is identical.

### 2.4 The frozen rule

| condition | measured | verdict |
|---|---|---|
| 1. M_resp interval entirely above 0 | +0.2099 **[−0.0130, +0.4263]** | **FAIL** |
| 2. `gen_pc20` > 0 with its interval entirely above 0 | +0.2952 **[+0.0695, +0.5065]** | **PASS** |

**C3 FAILS.** It fails on condition 1 only, and narrowly (lower bound −0.013 on 82 targets). The
honest reading: on 82 Cas9-knockout targets in 12 disjoint cancer lines the response-geometry model
transfers at the same level as the 28-bank map's own conservation score (+0.2099 vs +0.2104,
paired +0.002) — and **neither is separated from zero at this n**, while the annotation model and
the single spectral feature both are. A negative on 82 targets with a ±0.22 interval is a
low-power negative; it is not evidence that the association is absent, and it is written as a FAIL
because that is the rule that was frozen.

---

## 3. C4 — GWCD4i, the authors' own CD4-vs-K562 conservation

n = 1,773 genes, 408 component blocks, component-block bootstrap, 2,000 draws. **Primary condition
`Stim48hr`**; the other two are secondary and read the same.

| predictor | Rest | Stim8hr | **Stim48hr (primary)** |
|---|---|---|---|
| `n_degs_MASH_<cond>` — the condition's own DE count | +0.5065 [+0.4288, +0.5289] | +0.5251 [+0.4582, +0.5508] | **+0.5376 [+0.4723, +0.5633]** |
| `n_degs_MASH_K562` — magnitude in K562 | +0.2802 [+0.1822, +0.3061] | +0.2975 [+0.1937, +0.3239] | **+0.2959 [+0.1921, +0.3224]** |
| **M_ann — the annotation ridge (baseline)** | +0.2847 [+0.1899, +0.3169] | +0.2813 [+0.2019, +0.3245] | **+0.2702 [+0.1922, +0.3111]** |
| `C_held` — the map's own conservation | +0.2682 [+0.2003, +0.3209] | +0.2620 [+0.1860, +0.3071] | **+0.2739 [+0.2104, +0.3319]** |
| M_resp, 28-bank sensitivity arm (**labelled**) | +0.2332 [+0.1630, +0.2864] | +0.2295 [+0.1654, +0.2857] | **+0.2279 [+0.1637, +0.2776]** |
| **M_resp — 23 K562-free banks (the frozen cell)** | +0.2098 [+0.1267, +0.2547] | +0.2047 [+0.1238, +0.2474] | **+0.2032 [+0.1296, +0.2517]** |
| M_resp, out-of-fold variant | +0.1995 [+0.1286, +0.2551] | +0.1937 [+0.1247, +0.2476] | +0.1920 [+0.1293, +0.2514] |
| M_resp, K562-free training target (**labelled**) | +0.1907 [+0.1117, +0.2344] | +0.1843 [+0.1036, +0.2286] | +0.1846 [+0.1153, +0.2352] |
| **`gen_pc20`** | +0.1900 [+0.1055, +0.2222] | +0.1906 [+0.1097, +0.2264] | **+0.1858 [+0.1189, +0.2324]** |
| **`breadth`** | +0.1344 [+0.0316, +0.1662] | +0.1382 [+0.0363, +0.1701] | **+0.1275 [+0.0160, +0.1555]** |
| `logmag` | +0.1037 [+0.0001, +0.1333] | +0.1041 [−0.0018, +0.1322] | +0.1005 [−0.0139, +0.1270] |
| the authors' `donor_correlation_mean` (their own cross-donor reproducibility; **subset**) | +0.1689 [+0.0207, +0.2431] (n 524) | +0.1047 [−0.1648, +0.1541] (n 587) | +0.1009 [−0.0615, +0.1819] (n 590) |

Paired deltas from the same resample (Stim48hr): M_resp − M_ann **−0.0641 [−0.1236, +0.0112]**
(P(Δ ≤ 0) 0.964); M_resp − `C_held` −0.0739 [−0.1344, −0.0267]; M_resp − `n_degs_MASH_K562`
−0.0831 [−0.1395, +0.0121]; M_resp − `gen_pc20` +0.0163 [−0.0372, +0.0691].

**Null**: 200 gene shuffles read +0.0003 ± 0.0241; **0/200** reach +0.2032; z **+8.4** (Rest z +8.5,
Stim8hr z +9.1).

### 3.1 Partial, given the authors' own magnitude (`n_degs_MASH_K562` + the condition's count)

| predictor | Rest | Stim8hr | **Stim48hr** |
|---|---|---|---|
| M_ann | +0.1499 [+0.0684, +0.1897] | +0.1375 [+0.0714, +0.1934] | **+0.1168 [+0.0543, +0.1794]** |
| `C_held` | +0.1135 [+0.0820, +0.2139] | +0.0976 [+0.0595, +0.1894] | **+0.1100 [+0.0796, +0.2084]** |
| `gen_pc20` | +0.0879 [+0.0104, +0.1287] | +0.0826 [+0.0127, +0.1242] | **+0.0720 [+0.0308, +0.1406]** |
| **M_resp** | +0.0816 [+0.0393, +0.1650] | +0.0678 [+0.0282, +0.1560] | **+0.0610 [+0.0296, +0.1563]** |
| `breadth` | +0.0492 [−0.0512, +0.0846] | +0.0520 [−0.0471, +0.0855] | +0.0361 [−0.0797, +0.0646] |

### 3.2 The frozen rule (`Stim48hr`)

| condition | measured | verdict |
|---|---|---|
| 1. ρ ≥ +0.15 with the interval above 0 | +0.2032 [+0.1296, +0.2517] | **PASS** |
| 2. beats M_ann, paired Δ interval above 0 | **−0.0641 [−0.1236, +0.0112]** | **FAIL** |
| 3. partial given magnitude > 0, interval above 0 | +0.0610 [+0.0296, +0.1563] | **PASS** |
| 4. `gen_pc20` AND `breadth` positive, intervals above 0 | +0.1858 [+0.1189, +0.2324] and +0.1275 [+0.0160, +0.1555] | **PASS** |

**C4 FAILS** — on condition 2 only, in all three conditions (Rest and Stim8hr identical:
1 PASS / 2 FAIL / 3 PASS / 4 PASS). The response-geometry model does predict the authors' published
conservation, survives their own magnitude, and its two geometry columns carry the predicted
(inverted) sign; it does **not** beat annotation, which is what 290's headline claimed.

---

## 4. What this does to F13 (programme rule 10)

* **CONFIRMED, on two independent surfaces.** The *sign* of F13's mechanism — broad, spectrally
  central responses are MORE conserved across contexts — reproduces off-map on both a Cas9-knockout
  cancer-line panel and a primary-T-cell-vs-K562 comparison published by other authors, as a single
  feature with its interval above 0 in both (`gen_pc20` +0.2952 and +0.1858), and it survives the
  authors' own magnitude covariates (+0.0720). This was the part of F13 that was most surprising and
  most at risk (290's own P3/P4 were preregistered with the opposite sign and failed).
* **NARROWED, and this is the correction.** F13's headline number — response features beat
  annotation by +0.3948 — **is specific to 290's target**, the remainder after the annotation model,
  where annotation is orthogonal to the target by construction. Against raw conservation on new
  data the annotation ridge is as good or better (Broad-16 +0.2996 vs +0.2099; GWCD4i +0.2702 vs
  +0.2032; 268's shared-energy fraction +0.3528 vs +0.1118). F13's sentence "a pathway-level feature
  set is not what is needed; a spectral one is" does not hold outside the map.
* **NARROWED, magnitude.** Off-map, the published conservation statistic is mostly magnitude
  (`n_degs_MASH_<cond>` +0.5376, `n_degs_MASH_K562` +0.2959), and the geometry model keeps about a
  quarter of its rank correlation after both are partialled out. R290 §8.3's "it is not magnitude"
  is a statement about the map's remainder; on a published cross-cell-type statistic magnitude is
  the dominant term and geometry is an increment.
* **CONTRADICTED (a 290 post-hoc, not a 290 verdict).** 290 §8.4's DNA-binder exception — narrow,
  non-generic responses of sequence-specific DNA binders being MORE conserved than the geometry says
  — does not reproduce on the CD4 surface (§5 B1): ZNFs sit at the model and annotated TFs sit
  *below* it.
* **EXTENDED.** The 28-bank map itself (`C_held`) predicts both surfaces about as well as the model
  fitted on it (Broad-16 +0.2104 vs +0.2099; GWCD4i +0.2739 vs +0.2032). **The map transfers; the
  fitted feature model adds nothing over the map on either surface.** That is the strongest single
  argument for the map (`E279/work/znf279.csv`) being the object worth keeping, which is what 284's
  synthesis already said for a different reason.

---

## 5. Descriptive tail (POST-HOC, labelled; written after both verdicts were read)

**A1 — is `C_broad16` above chance, and where do the OR controls sit?** Wrong-target permutation
(20 independent draws, target t's row in line i paired with another target's row in line j, same
centring): null median **+0.00313**, p95 +0.02798, against an observed median of **+0.02957**;
**52.4 %** of the 82 targets sit above the null's p95. The statistic is target-specific.
**The 7 olfactory-receptor controls read +0.02672 (MWU targets > OR, p = 0.253) — NOT separated
from the real targets, and this row is not a clean negative control**: OR rows are built against the
other nine controls and then projected on a centring computed from the *target* rows, so their
common component is not removed, while the targets' is. The wrong-target permutation is the
meaningful null here and the OR row should not be read as one.

**A2 — the Broad-16 extremes** (rank residual, observed minus predicted). Most **under**-predicted
(more conserved than the geometry says): `CAB39`, `RAB7A`, `SMG6`, `AURKB`, `SEC23IP`, `YTHDF2`,
`KCTD5`, `TOP2A`, `SLC39A10`, `UBE2S`, `BRD4`, `DNAJC13`. Most **over**-predicted: `SLC25A3`,
`CD46`, `LRP8`, `KIF2C`, `RAD51`, `CDK4`, `IER3IP1`, `TFRC`, `USP37`, `SCYL1`, `UFL1`, `MDM2`. No
family pattern of 290's kind is visible; the under-predicted tail is trafficking / mRNA-surveillance
/ mitotic-exit rather than DNA binders. **290's exception is barely testable on this panel**: of the
82 unit targets, 279's classes are 77 `NONZNF`, 4 `OTHER_ZF`, **1 `ZNF`**, and only 5 are annotated
TFs — a limitation, not a refutation.

**A3 — two views on Broad-16.** ρ(M_resp, `C_broad16` | M_ann) = **+0.1621**;
ρ(M_ann, `C_broad16` | M_resp) = **+0.2699**; ρ(M_resp, M_ann) = +0.1941. Both models carry rank
information the other does not.

**B1 — the DNA-binder exception on GWCD4i** (mean rank residual, outcome minus M_resp; one-sided
MWU against the rest):

| group | Rest | Stim8hr | Stim48hr |
|---|---|---|---|
| ZNF (279's label), n 79 | −0.0048 (p 0.46) | −0.0041 (p 0.48) | −0.0039 (p 0.44) |
| TF (CollecTRI ∪ DoRothEA), n 148 | −0.0667 (p 0.98) | −0.0541 (p 0.96) | −0.0550 (p 0.96) |
| other, n 1,571 | +0.0075 (p 0.020) | +0.0063 (p 0.035) | +0.0063 (p 0.044) |

**B2 — do the two feature sets combine?** A 39-column ridge (18 response + 21 annotation), same
protocol: Stim48hr **+0.2166 [+0.1417, +0.2583]**, against M_ann's +0.2702 — paired
**−0.0514 [−0.1060, +0.0170]**, i.e. the combination does not beat annotation alone even though
ρ(M_resp, y | M_ann) = +0.1113 and ρ(M_ann, y | M_resp) = +0.2122 (Rest and Stim8hr the same).
The two predictions correlate +0.3864. **The response columns add rank information that the ridge
cannot convert into a better fit** — the fit is chosen on the map, where the response block is much
the stronger (grouped-CV +0.383 vs +0.226), and that choice does not survive the transfer.

---

## 6. Limitations, stated plainly

1. **Neither surface is untouched data and neither is a reserve.** 268 / 283 / 289 read Broad-16;
   287 read this GWCD4i table. Only the comparison with these features is new. A future claim on
   either surface is a third look at it.
2. **C3 has 82 targets and 21 family components.** The paired-target bootstrap interval is ±0.22
   wide; the test could not have detected anything below ρ ≈ 0.22 with confidence. The frozen rule
   was set knowing n ≈ 92. Reported as a FAIL; read as low-powered.
3. **The Broad-16 panel is a growth/essentiality panel** (90 analysis targets, mostly core
   machinery): of the 82 unit targets, 279's class labels are 77 `NONZNF` / 4 `OTHER_ZF` / 1 `ZNF`
   and only 5 are annotated TFs, so the family axis 279 measured — core machinery to cell identity —
   is largely compressed out of it. That is the most likely reason both models read low there.
4. **The GWCD4i outcome is the authors' statistic, computed by their pipeline** (DESeq2 LFC
   correlation vs K562, MASH DE sets). It is not 276's cosine and not corrected for reliability;
   its dominant covariate is their own DE count (+0.54).
5. **`n_degs_MASH_<cond>` is not an independent covariate.** Partialling it out removes part of the
   outcome's own construction as well as the confound; §3.1's partial is a floor, not an unbiased
   effect.
6. **The training target `C_held` includes K562 pairs** even in the K562-free arm; the labelled
   sensitivity that also removes them from the target reads +0.1846 — no verdict changes.
7. **Both models are fitted to predict `C_held` and applied without recalibration.** A model refit
   on either surface would read higher; that would not be a test of F13.
8. **Multiplicity.** Verdict cells: **2** (C3, C4). Reported in total: ~165 cells — 20 Spearman
   rows across surface A's three units, 3 paired deltas, 1 null, 6 depth rows; 33 Spearman rows
   across surface B's three conditions, 12 paired deltas, 6 donor rows, 15 partials, 3 nulls; ~66
   post-hoc cells. All but the two verdict cells are secondary, comparative, or post-hoc.

---

## Files

| path | what |
|---|---|
| `PREREGISTRATION.md` / `work/PREREGISTRATION.sha256` | frozen 15:13:32Z, `3a4f5f85…`; every script aborts on a digest mismatch |
| `lib291.py` | copied metrics / models / bootstraps (provenance line-cited in its docstring), the per-source aggregator, the fit-then-apply protocol, the partial Spearman |
| `ka291.py` → `work/ka291.json` | KA-1, KA-2, KA-5 |
| `surfA291.py` → `work/pred291_broad16.csv`, `work/surfA291.json` | KA-4, `C_broad16`, C3's table, the second outcome, the null, depth |
| `surfB291.py` → `work/pred291_gwcd4i.csv`, `work/surfB291.json` | KA-3, C4's table per condition, the partials, the sensitivity arms |
| `desc291.py` → `work/desc291.json` | POST-HOC: A1–A3, B1–B2 |
| `work/recon291.py`, `work/recon291b.py` | the pre-freeze reconnaissance (counts and column availability only) |
| `work/*.log` | run logs |
| `docs/CONFIRMATION-LEDGER.md` | claims C3 / C4 appended before the computation, results appended after |
