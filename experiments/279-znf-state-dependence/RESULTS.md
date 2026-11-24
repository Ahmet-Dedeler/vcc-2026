# 279 — The preregistered zinc-finger claim, tested on 17,106 HELD (non-panel) perturbed genes

**CONFIRMED on held targets (median difference −0.01047 [−0.01269, −0.00844], AUC 0.564 [0.548, 0.581]; matched median difference −0.00635 [−0.01172, −0.00079]).** Both arms of the frozen pass rule fire in the predicted direction.

**But the claim's own framing does not survive its controls.** The deficit is *not* zinc-finger-specific: homeobox families read **larger** (Δ −0.01476, AUC 0.603) and bHLH the same (Δ −0.01169, AUC 0.574); among the 93 HGNC gene groups with ≥ 30 held genes the two zinc-finger groups rank **13th** (KRAB, Δ −0.01205) and **21st** (C2H2, Δ −0.01060) from the bottom, behind interleukin receptors, olfactory receptors, PAS-domain proteins, PRD/HOXL/NKL homeoboxes and Ig-like cell-adhesion molecules. What the map actually separates is **core cellular machinery** (spliceosome +0.107, SSU processome +0.107, proteasome +0.095, ribosome biogenesis +0.074, DEAD-box helicases +0.073 — all *conserved*) from **cell-identity and signalling families** (*state-dependent*), and the zinc fingers sit on the identity side without being special within it. Post hoc, roughly half the ZNF deficit is generic "TF-ness" (CollecTRI TFs vs non-TFs, Δ −0.00520) and half is zinc-finger-specific relative to other curated TFs (Δ −0.00562 [−0.00906, −0.00271], matched −0.00406 [−0.00815, −0.00028]).

**The one place the held test overturns 276.** 276 read the in-sample effect as "mostly a weak, poorly-measured response" (R276 §5). On held genes with four sources' *measured* split-half reliability that reading fails: the rank effect is unchanged by the reliability correction (AUC 0.561 uncorrected → **0.567** corrected on the same genes) and the deficit is **largest among the best-measured genes** (AUC 0.517 / 0.551 / **0.614** across reliability terciles; Δ −0.0037 / −0.0076 / **−0.0267**). Noise cannot produce that ordering. The zinc-finger deficit is **state dependence, not measurement** — with the caveat that the two supporting statistics are *unmatched*: inside the best-measured tercile the expression/magnitude-matched interval is −0.0143 [−0.0402, +0.0083] on n = 161 ZNF genes, i.e. the ordering is clear and the confound-adjusted version of its top cell is underpowered (§6iii).

**Date 2026-09-06, 06:21–07:05 UTC (09:21–10:05 Ankara).** Deterministic (seed 20260906). **Measurement only: no build, no `vcc`, no `cell-eval2`, no model, no submission, nothing under `submissions/` or `data/`, no download, no production file touched** — everything was already on disk. `uv run python` throughout. Peak RSS 4.44 GiB (`build279`), 6.65 GiB (`rel279`), 3.99 GiB (`analyse279`) — measured both as `ru_maxrss` and by sampling the live process (5.97 GiB concurrent for `rel279`). `find` over the repo for files newer than `PREREGISTRATION.md` returns, outside this directory, only other agents' work (`docs/LOG.md`, `docs/READ-FIRST.md`, `experiments/278-*`, `experiments/280-*`) — plus **one file of mine, now removed**: a mistyped shell line ran `rel279.py` from the repo root and left a 176-byte `work/rel279.log` there holding nothing but the interpreter's "can't open file" error. It was deleted (`rm work/rel279.log`); the pre-existing root `work/` directory and its other files were untouched, and no other path outside `experiments/279-znf-state-dependence/` was written.

**PREREGISTRATION.md** sha256 `a3574d7a92c1bfc3182ffa8d51b625f5ff08ca7ea8ee30cd14bac19e1ddfced9`, frozen **2026-09-06 06:25:54 UTC** (`work/PREREGISTRATION.sha256`), after the metadata-only inventory and before a single held-gene cosine existed. Everything labelled POST-HOC below was written after the preregistered numbers were read and cannot change the verdict.

---

## 0. Known answers (all read before any held number)

| row | required | measured | verdict |
|---|---|---|---|
| **KA-0** the block reader reproduces `lib276.load_source_path` | max abs diff < 1e-6, identical availability | identical index, **1.17e-07** | PASS |
| **KA-1** 276's panel map rebuilt from its 26 sources | axis 2,994 / 315 cross-line pairs / median +0.0198672 / 16.0 % ≤ 0, and `conservation276.csv`'s `C_raw` to < 1e-6 | 2,994 / 315 / **+0.0198672** / **0.160**; max abs diff **4.99e-07** | PASS |
| **KA-2** 276 §5's post-hoc panel ZNF observation (symbol rule) | 35 targets, median +0.0044 vs +0.0223, MWU p 1.1e-05, 23.3 % of the bottom decile / 6.7 % of the top | **35**, **+0.0044** vs **+0.0223**, p **1.08e-05**, **23.3 %** / **6.7 %** | PASS |
| **KA-3** DLD-1 replicate-pair reliability median on the panel (rel276.json +0.086321) | ±1e-5 | **+0.0863211** | PASS |
| **KA-4** the rebuilt K562/HCT116/HEK293T half rows vs 271's | 271's medians +0.04706 / +0.02756 / +0.02825 | **+0.05262 / +0.03113 / +0.03543** — NOT the same number | see below |
| **KA-4b** per-target agreement with 271's own `r_half_raw` | — | Spearman **+0.968 / +0.974 / +0.980** (Pearson +0.985/+0.984/+0.991) on n = 269/290/294 | PASS |

KA-4 as written in the preregistration was mis-specified by me, and I record that rather than paper over it: 271's `r_half_raw` lives on the **18,533-gene official axis**, centred on the bank's *every-perturbation* mean with the panel's own gene columns zeroed, while this experiment's reliability is computed on the **7,543-gene half-map axis** centred on the rows of *this* map (DEFINITION-276 §2). The two are not the same statistic and their medians differ by 12–25 %. KA-4b is the check that was actually available — the per-target ordering agrees at ρ ≥ 0.968 — and KA-3, which *is* an identical-convention number, reproduces to seven decimals. The rebuild is correct; the preregistered tolerance was not.

---

## 1. The held set (item 1)

**Held test set: 17,106 genes** — every non-panel perturbed gene that appears in ≥ 2 included sources from **two different cell lines** and is one of the challenge's 18,533 symbols. The 300 panel targets (the hypothesis-forming set) are excluded by construction, so nothing that produced the observation is in the test.

| class (HGNC `gene_group`, frozen §2) | n |
|---|---:|
| **ZNF** = "Zinc fingers C2H2-type" ∪ "KRAB domain containing" | **696** (331 KRAB, 365 other C2H2) |
| **non-ZNF** = no zinc-finger group of any type | **16,060** |
| other zinc-finger families (CCCH, FYVE, matrin, CCHC, RANBP2, DHHC, MYND, ZZ, GATA …) — excluded from both | 350 |

Membership source: `data/external/xcell/hgnc_protein_coding.txt` (the HGNC complete protein-coding set, `symbol` + `gene_group`; genenames.org gene-group taxonomy). 276's post-hoc symbol rule (`ZNF*`/`ZKSCAN*`/`ZSCAN*`) is carried as a declared secondary; it selects 475 held genes.

**28 sources / 21 cell lines**, admitted by DEFINITION-276 §1's rule with "panel targets" replaced by "held genes" (≥ 5): 276's 26 minus the five that carry **0** held genes (`kolf2`, `feng_ipsc` and the three CD4 conditions are panel-only banks: every perturbation they contain is one of the 300) plus the seven 276 excluded only for carrying 0–1 *panel* targets (`k562_essential`, `rpe1`, `nadig_jurkat`, `nadig_hepg2`, `tian2019_ipsc`, `gasperini_k562`, `wessels_thp1`). Full inventory in `work/inventory279_sources.csv`. The widest: `orion_hek293t` 16,721 held genes, `orion_hct116` 16,594, `k562_gwps` 9,036, `perturbmatch` 8,486, `dld1_lowmoi` 4,459, `gse295214_h1` 2,128, `nadig_jurkat` 2,059, `rpe1` 1,997, `nadig_hepg2` 1,967, `igvf_telohaec` 1,942, `k562_essential` 1,821, `igvf_cardiomyocyte` 1,707.

*(One correction against the preregistration: its §1 says "28 sources / **24** distinct cell lines". The source **set** is exactly the one the rule specifies — 28 banks — but they span **21** lines, not 24; the count in the frozen text was mine and was wrong. Nothing downstream uses it.)*

**Common gene axis: 2,739 genes** (the intersection over the 28 sources ∩ the challenge's 18,533; 276's panel axis was 2,994 over 26 sources). **364 cross-line pairs**, 14 same-line pairs.

### 1.1 The held map itself (a new scientific object)

| statistic | mean | median | p10 | p90 | share ≤ 0 | n |
|---|---:|---:|---:|---:|---:|---:|
| **`C_held`** cross-line, 364 pairs (median 5/gene, max 164) | +0.0308 | **+0.0186** | −0.0379 | +0.1259 | **30.6 %** | 17,106 |
| `C_within` same-line (14 pairs) | — | +0.5965 | — | — | — | 8,486 |
| target-shuffled null, 50 draws | +0.00023 | +0.00028 | — | — | — | 17,106 |
| (276's panel map, for reference) | +0.0263 | +0.0199 | −0.0028 | +0.0676 | 16.0 % | 300 |

The held map's median conservation (+0.0186) is within a hair of the panel's (+0.0199) even though the held genes were never selected for measurability — the panel's advantage shows up in the **tails**, not the centre: 30.6 % of held genes read ≤ 0 against 16.0 % of panel targets, and the shuffled null's spread over genes is 0.0107 against the panel's 0.0013, because the median held gene has 5 pairs and the median panel target had 33. **Per-gene scores on the held map are far noisier than on the panel map; only group statistics are read here.**

`nadig_jurkat` is Jurkat = challenge context A. It carries **0/300 panel targets**, so no panel response can leak through it; dropping it moves the primary effect from −0.01047 to −0.01012 (§5).

---

## 2. The frozen prediction and pass rule (item 2, quoted from PREREGISTRATION §4)

> Predicted direction (one-sided): `median(C_held | ZNF) − median(C_held | non-ZNF) < 0`.
> Unit of resampling: **the perturbed gene**. Confounds named in advance: expression level (CCLE median over 1,600 lines) and response magnitude (per-source scale-normalised row norm), matched 1:1 without replacement inside the same (expr-decile × mag-decile) cell, 50 draws.
> **CONFIRMED iff both** (1) the raw Δmed is negative with its 95 % bootstrap interval entirely below 0, **and** (2) the matched Δmed is negative with its 95 % interval entirely below 0.

---

## 3. The result (item 4)

| test | ZNF median | non-ZNF median | **Δ median** | 95 % CI | AUC (P[ZNF < non-ZNF]) | n |
|---|---:|---:|---:|---|---:|---|
| **PRIMARY, raw** | +0.00878 | +0.01925 | **−0.01047** | **[−0.01269, −0.00844]** | **0.564** [0.548, 0.581] | 696 / 16,060 |
| **PRIMARY, matched** (expr × mag deciles, 50 draws, 0 dropped) | — | — | **−0.00635** | **[−0.01172, −0.00079]** | 0.530 | 696 pairs |

Mann–Whitney one-sided p = **4.06e-09**. Mean difference −0.01412. **Both pass-rule conditions fire → CONFIRMED.**

Three things about the matched arm, stated because they are the weak link:

* the per-draw spread over the 50 matching draws is [−0.00959, **+0.00006**] — one draw touches zero. Repeated post hoc with **200 draws**: Δ −0.00664 [−0.01205, −0.00157], per-draw [−0.01141, +0.00006]. The interval clears, the draw distribution grazes zero;
* the matching-free version of the same adjustment (POST-HOC): the **stratified AUC** — P(a ZNF gene sits below a non-ZNF gene *from the same expr × mag cell*) — is **0.5299 [0.5124, 0.5478]** over 1.26e+05 within-cell gene pairs, and **0.5536 [0.5247, 0.5825]** on the balanced core of §4(c3). Confound-adjusted, the effect is real and small;
* the named confounds are small but real and in the expected direction: ZNF CCLE median expression 2.296 vs 2.405 (p 0.031), response magnitude 0.966 vs 0.992 (p 0.0040). The one large imbalance runs the *other* way — ZNF genes have **more** measurements (median 9 cross-line pairs vs 5, 4 lines vs 3, p ≈ 1e-106) — which §4(c3) removes by construction.

---

## 4. Controls (item 3)

**(a) Size-matched random families.** 100 draws of 696 genes from the held non-ZNF pool: Δmed null band **[−0.0030, +0.0045]**, and **0/100** draws reach the observed −0.01047. The observed effect is far outside the sampling band for a set of this size.

**(a2) Real HGNC families — the family-level null (a stronger control than (a)).** All 93 HGNC gene groups with ≥ 30 held genes, ranked by their own Δmed against the non-ZNF rest:

| rank | family | n | Δ median |
|---:|---|---:|---:|
| 1 | Interleukin receptors | 37 | −0.02131 |
| 2 | Olfactory receptors, family 2 | 58 | −0.02082 |
| 3 | PAS domain containing | 32 | −0.01650 |
| 4 | PRD class homeoboxes | 48 | −0.01518 |
| 5 | Ig-like cell adhesion molecules | 61 | −0.01510 |
| 6 | HOXL subclass homeoboxes | 51 | −0.01435 |
| 7 | RAB GTPases | 61 | −0.01344 |
| 8 | SCAN domain containing | 47 | −0.01343 |
| 9 | TBC1 domain family | 40 | −0.01331 |
| 10 | NKL subclass homeoboxes | 42 | −0.01313 |
| 11 | GPCRs, class A orphans | 66 | −0.01302 |
| 12 | C-type lectin domain containing | 80 | −0.01234 |
| **13** | **KRAB domain containing** | **331** | **−0.01205** |
| 15 | Basic helix-loop-helix proteins | 106 | −0.01169 |
| **21** | **Zinc fingers C2H2-type** | **691** | **−0.01060** |
| … | … | | |
| 87–93 | AAA ATPases +0.028, snoRNA host genes +0.036, RNA-binding-motif +0.047, RNA polymerase subunits +0.053, DEAD-box helicases +0.073, ribosome biogenesis +0.074, proteasome +0.095, spliceosomal P complex +0.096, SSU processome +0.107, spliceosomal C complex +0.107 | | |

At the family level the zinc fingers are in the low tail but **not exceptional** (empirical one-sided family rank 13/93 and 21/93). The two zinc-finger groups are by far the largest families in the table, so their Δ is the most precisely estimated one — a like-for-like family null is not available at n ≈ 700 — but the ranking is the honest multiplicity statement. (SCAN at rank 8 overlaps the ZNF class: SCAN-domain genes are mostly C2H2 zinc fingers.)

**(b) The other large TF families — is it ZNF-specific?** **No.**

| family | n | median | Δ median | 95 % CI | AUC | matched Δ |
|---|---:|---:|---:|---|---:|---|
| Homeoboxes (all HGNC subclasses) | 219 | +0.00483 | **−0.01476** | [−0.01713, −0.00961] | **0.603** | −0.00724 [−0.01779, +0.00412] |
| Basic helix-loop-helix | 106 | +0.00763 | **−0.01169** | [−0.01803, −0.00627] | 0.574 | −0.00776 [−0.02560, +0.00294] |
| (zinc fingers, for comparison) | 696 | +0.00878 | −0.01047 | [−0.01269, −0.00844] | 0.564 | −0.00635 [−0.01172, −0.00079] |

**(c) Single-dataset replications — is it a dataset mix?** **No, where it is powered.**

| variant | pairs used | ZNF median | non-ZNF median | Δ median | 95 % CI | AUC | matched Δ |
|---|---|---:|---:|---:|---|---:|---|
| **(c1)** Orion HCT116 × HEK293T only | 1 | +0.01906 | +0.02808 | −0.00901 | [−0.01573, −0.00091] | 0.521 | −0.00321 [−0.01513, +0.00971] |
| **(c2)** Replogle K562-GWPS × RPE1 only | 1 | +0.09342 | +0.09585 | −0.00243 | [−0.03252, +0.02678] | 0.513 | +0.02790 [−0.01880, +0.07539] |
| **(c3)** balanced core {K562-GWPS, Orion HCT116, Orion HEK293T}, genes in all three | 3, identical for every gene | +0.01181 | +0.02966 | **−0.01785** | **[−0.02197, −0.01231]** | **0.604** | **−0.00792 [−0.01536, −0.00015]** |

(c1) reproduces the effect inside one dataset. (c2) is a null with n = 47 ZNF genes — Replogle's RPE1 bank is an *essential-gene* screen and zinc fingers are almost never essential, so this cell is underpowered rather than contradictory. **(c3) is the cleanest reading in the whole experiment**: with every gene carrying exactly the same three cross-line pairs, the "which sources measured this gene" confound is gone and the effect is *larger* (Δ −0.01785, AUC 0.604), with the matched interval also excluding zero. On the balanced core the zinc fingers are **10.8 % of the bottom decile and 1.9 % of the top** (7.9 % of the set), and **38.2 %** of them read ≤ 0 against 25.4 % of non-ZNF — the same signature 276 reported in-sample (23.3 % / 6.7 %).

**(d) The map's own shuffled null**: target-shuffled pairings, 50 draws, mean +0.00023 (median +0.00028) against the map's +0.0186 — the held map is above chance, as the panel map was.

---

## 5. Declared secondaries

| variant | n ZNF | Δ median | 95 % CI | AUC | matched Δ |
|---|---:|---:|---|---:|---|
| KRAB-domain only | 331 | −0.01205 | [−0.01544, −0.00948] | 0.589 | −0.00806 [−0.01786, −0.00095] |
| 276's symbol rule `ZNF*`/`ZKSCAN*`/`ZSCAN*` | 475 | −0.01179 | [−0.01454, −0.00903] | 0.575 | −0.00792 [−0.01436, −0.00072] |
| without `nadig_jurkat` (challenge line A) | 696 | −0.01012 | [−0.01249, −0.00810] | 0.564 | — |

Both alternative membership rules give the same answer, and the KRAB subfamily — the one the mechanistic story names — reads slightly *stronger* than the C2H2 family as a whole.

## 5.1 In-sample reference (item 4, **labelled: this is the hypothesis-forming set**)

The same test on the 300 panel targets with 279's sources and axis:

| rule | n ZNF | ZNF median | rest | Δ median | 95 % CI | AUC |
|---|---:|---:|---:|---:|---|---:|
| HGNC (this experiment's) | 48 | +0.01209 | +0.03207 | −0.01998 | [−0.03033, −0.00492] | 0.648 |
| 276's symbol rule | 35 | +0.00498 | +0.03172 | −0.02675 | [−0.03547, −0.01561] | 0.696 |
| 276's own published panel numbers (26 sources, 2,994-gene axis) | 35 | +0.0044 | +0.0223 | −0.0179 | — | — |

**The held effect is 0.44–0.52× the in-sample effect** (−0.01047 vs −0.01998 on the HGNC rule; −0.01179 vs −0.02675 on the symbol rule) and the AUC falls from 0.65–0.70 to 0.56–0.58. That is the expected shrinkage of an effect discovered post hoc in one screen and re-measured on data that did not select it — the claim survives, at about half its discovered size.

---

## 6. Noise or state dependence? (item 5, POST-HOC as declared)

Reliability comes from four sources with independent halves on non-panel genes: 267's batch-parity accumulators for `k562_gwps` / `orion_hct116` / `orion_hek293t` (half rows rebuilt here with `truth267.write_half_bank`'s arithmetic, KA-3/KA-4b) and DLD-1's two low-MOI replicates. Half-map axis 7,543 genes, 6 cross-line pairs.

**(i) Are ZNF responses measured worse?** A little, and only in one source:

| source | r½ ZNF | r½ non-ZNF | AUC | one-sided p | share r ≤ 0.02 |
|---|---:|---:|---:|---:|---|
| `k562_gwps` | +0.0275 | +0.0438 | 0.624 | 1.5e-27 | 37 % vs 25 % |
| `orion_hct116` | +0.0145 | +0.0210 | 0.530 | 3.9e-03 | 55 % vs 49 % |
| `orion_hek293t` | +0.0272 | +0.0359 | 0.523 | 2.1e-02 | 46 % vs 43 % |
| `dld1_lowmoi` | +0.0609 | +0.0632 | 0.511 | 0.31 | 13 % vs 12 % |

**(ii) Does the correction remove the deficit?** No — the *rank* effect is untouched:

| score (same four sources) | ZNF | non-ZNF | Δ median | 95 % CI | **AUC** |
|---|---:|---:|---:|---|---:|
| `C_full4` uncorrected (n 687/15,831) | +0.0094 | +0.0183 | −0.0089 | [−0.0119, −0.0056] | **0.561** |
| `C_rel` reliability-corrected (n 502/9,941) | +0.1438 | +0.2300 | −0.0862 | [−0.1216, −0.0575] | **0.567** |
| `C_rel` without `orion_hek293t` (n 359/6,093) | +0.1273 | +0.1957 | −0.0684 | [−0.1098, −0.0293] | 0.567 |

Dividing by √(r·r′) inflates the whole map ~8.8× (non-ZNF median +0.0183 → +0.2300); if the ZNF deficit were measurement, it would shrink toward zero under that operation. It does not: the corrected Δ is −0.0862 against −0.1164 predicted by a purely-biological deficit (Δ_uncorrected × 8.79), i.e. **74 % of the deficit survives**, and the AUC — which is scale-free — moves from 0.561 to 0.567.

**(iii) The decisive ordering.** Split all held genes into terciles of measured reliability (median r ≈ −0.007 / +0.036 / +0.105) and repeat the contrast inside each:

| tercile | r median | n ZNF / non-ZNF | Δ median | 95 % CI | **AUC** | matched Δ |
|---|---:|---|---:|---|---:|---|
| 0 (unmeasurable) | −0.007 | 258 / 5,344 | −0.00367 | [−0.00599, −0.00089] | **0.517** | −0.00294 [−0.00869, +0.00243] |
| 1 | +0.036 | 277 / 5,295 | −0.00756 | [−0.01079, −0.00486] | **0.551** | −0.00568 [−0.01290, +0.00089] |
| 2 (best measured) | +0.105 | 161 / 5,413 | −0.02666 | [−0.03302, −0.01904] | **0.614** | −0.01426 [−0.04017, +0.00827] |

A noise explanation predicts the deficit shrinks as measurement improves. It **grows** — and it grows on the rank scale (AUC 0.517 → 0.551 → 0.614), where the mechanical inflation of cosines with reliability cannot act. Two honest qualifications: the tercile contrasts are **not** expression/magnitude-matched (the matched column is shown and only tercile 2's point estimate is large, with an interval that includes zero at n = 161), and the ZNF genes that reach the top tercile are a subset that this experiment cannot show is compositionally the same as the ZNF genes in the bottom one. **Where the response is measured well enough to resolve a direction, zinc-finger responses point in genuinely different directions in different cell lines.** This overturns 276 §5's in-sample reading ("mostly a weak, poorly-measured response"), which was made at n = 28 corrected panel targets.

---

## 7. What the effect actually is (POST-HOC, `work/tf279.json`)

**(1) Paralogue clustering does not inflate the interval.** KRAB-ZNFs sit in tandem arrays (the six largest held ZNF cytobands are 19q13.43 with 58 genes, 19p13.2 with 35, 19q13.41 with 26, 19q13.12 with 26, 19p12 with 26, 19q13.31 with 22), so the preregistered gene-level bootstrap treats correlated genes as independent. Resampling **cytobands** instead (261 ZNF bands / 1,425 non-ZNF bands, 5,000 draws) gives Δ **−0.01047 [−0.01285, −0.00785]** — indistinguishable from the gene-level interval.

**(2) Half of it is generic TF-ness, half is zinc-finger-specific.** Using CollecTRI's curated TF catalogue on disk (`data/external/regulon/collectri_omnipath.tsv`, 1,201 TFs, 1,116 in the held set of which 225 are HGNC zinc fingers):

| contrast | n | Δ median | 95 % CI | AUC | matched Δ |
|---|---|---:|---|---:|---|
| TFs (ZNFs removed) vs non-TFs | 891 / 15,519 | −0.00520 | [−0.00787, −0.00235] | 0.513 | −0.00282 [−0.00823, +0.00205] |
| **ZNF vs other TFs** | 696 / 891 | **−0.00562** | **[−0.00906, −0.00271]** | 0.561 | **−0.00406 [−0.00815, −0.00028]** |
| ZNF vs non-TF non-ZNF | 696 / 15,519 | −0.01082 | [−0.01315, −0.00879] | 0.566 | −0.00636 [−0.01125, −0.00108] |

Being a transcription factor costs about −0.005 of conservation; being a zinc finger costs about another −0.006 **on top of** being a transcription factor. So the claim is not merely "TFs transfer badly" — but it is also not unique to zinc fingers, since homeobox families (§4b) read lower still.

**(3) The stated mechanism (KRAB) is directionally supported and weak.** Inside the HGNC C2H2 family, KRAB-domain members (n 331, median +0.00720) vs non-KRAB C2H2 members (n 365, +0.01027): Δ −0.00307 [−0.00800, **+0.00045**], AUC 0.555, one-sided p 0.006. The KRAB subfamily is lower, but the median difference's interval includes zero — this experiment cannot separate "KRAB-ZNF" from "C2H2-ZNF" as the operative class.

---

## 8. Limitations, stated plainly

1. **Confirmed, but not the claim as framed.** The preregistered statistic passes; the mechanistic framing ("KRAB-ZNFs are special because they are lineage-restricted TE repressors") is *not* what the data singles out. Homeobox families read lower, and the zinc fingers rank 13/93 and 21/93 among real HGNC families. The defensible statement is the one in §7(2): a TF-family cost plus a zinc-finger-specific cost of similar size.
2. **The confound-adjusted effect is small.** Stratified AUC 0.530 [0.512, 0.548] (0.554 on the balanced core). This is a group-level property, not a per-gene predictor: a single held gene's `C_held` has a shuffled-null spread of 0.0107, half the size of the whole effect.
3. **One matching draw in 50 (and in 200) reads +0.00006.** The matched arm clears its interval, but it is close to the boundary; the matching-free stratified AUC is the number to quote.
4. **Held ≠ independent.** The held genes are new, but the *sources* are the same banks 276 used (plus seven). A shared analysis artefact in Orion or K562-GWPS would appear in both maps. (c1) and (c3) show it is not one dataset, and the cytoband bootstrap shows it is not one chromosome — but no genuinely new experiment was run.
5. **The 2,739-gene common axis** costs absolute cosine (276 §2.3 measured Spearman +0.85…+0.96 against the donors' own wider axis and a loss of ρ in *prediction*). Every number here is computed on the same axis for both groups, so the comparison is fair; the absolute values are not the widest-axis ones.
6. **Reliability coverage.** `C_rel` exists for 10,695 held genes (502 ZNF) after the r ≤ 0.02 floor drops 22,420 pairs; the tercile analysis uses the median reliability over the four halved sources, which for most genes is Orion's.
7. **The reserve was not touched.** Programme rule 2 reserves the `orion_hek293t` split-half **B truths on the 300 panel targets**. This experiment reads that line's batch-B accumulator only for **non-panel** genes, and every reliability statistic is also reported without `orion_hek293t` (§6ii, same AUC 0.567). No entry in `docs/CONFIRMATION-LEDGER.md` is required.
8. **Multiplicity.** Cells screened here: 1 preregistered primary + 1 matched arm; 4 controls; 3 secondaries; 3 reliability variants; 3 terciles; 4 post-hoc contrasts ≈ 19 reported cells, of which the verdict depends on 2 named in advance.

---

## 9. What this is worth

**Scientific track.** This is the programme's first *preregistered biological claim* — a statement about perturbation biology rather than about a predictor — tested on targets that did not exist in the map when it was formed. (Claim C1 in `docs/CONFIRMATION-LEDGER.md` was the first preregistered claim confirmed on *reserved* data; that reserve is untouched here, §8.7.) Three findings that outlive the verdict:

* **F-A. A 17,106-gene cross-line conservation map exists** (`work/znf279.csv`), built with 276's frozen statistic over 28 sources and 21 lines — 57× the panel map, and the first version of the map that is not restricted to the challenge's own designed panel. Its family structure is a clean axis from core machinery (spliceosome/proteasome/ribosome biogenesis, +0.07…+0.11) to cell-identity families (receptors, homeoboxes, zinc fingers, −0.01…−0.02).
* **F-B. Zinc-finger responses are state-dependent, not badly measured** — the deficit is *largest* among the best-measured genes (AUC 0.517 → 0.614 across reliability terciles) and survives the reliability correction (74 % of it). This corrects 276 §5's in-sample reading.
* **F-C. The cost decomposes**: −0.005 for being a transcription factor, another −0.006 for being a zinc finger, on the same scale as a map whose median is +0.019.

**Competition track: nothing to build.** 277 already closed the one competition-facing use of this class (rendering the identifiable low-conservation rows differently does not clear: −0.0043 [−0.0092, +0.0006] `pds`), and the effect here is a group-level property with a stratified AUC of 0.53 — far too weak to key an action on. This experiment changes no production file and proposes no arm.

**If it earns a stage 2**, the measurement that would move it is the depth axis 278 is inventorying: the tercile result says the effect is only visible where the response is resolved, so *more cells per condition on a few lines* — not more lines — is what would sharpen it, exactly as R276 §6 concluded on independent grounds.

---

## Files

| path | what |
|---|---|
| `PREREGISTRATION.md` + `work/PREREGISTRATION.sha256` | the frozen claim, sets, statistic, pass rule, controls (sha256 `a3574d7a…`, 06:25:54 UTC) |
| `inv279.py` → `work/inventory279.json`, `work/inventory279_sources.csv` | the metadata-only held-set inventory (no `delta.npy` opened) |
| `build279.py` → `work/held279.npz`, `work/build279.json`, `work/build279.log` | KA-0/KA-1 and the held map (per-pair cosines, magnitude, CCLE expression, shuffled null) |
| `analyse279.py` → `work/analyse279.json`, `work/znf279.csv` | the preregistered test, the controls, the secondaries, the in-sample reference |
| `rel279.py` → `work/rel279.npz`, `work/rel279.json` | split-half reliability on held genes, `C_rel`, KA-3/KA-4b |
| `desc279.py` → `work/desc279.json` | the noise-vs-state-dependence descriptor (§6) |
| `extra279.py` → `work/extra279.json` | POST-HOC: decile shares, stratified AUC, 200 matching draws |
| `tf279.py` → `work/tf279.json` | POST-HOC: cytoband block bootstrap, TF-class decomposition, KRAB vs other C2H2 |
| **`work/znf279.csv`** | the per-gene deliverable: 17,106 rows × 19 columns (`class`, `krab`, `symbol_rule`, `C_held`, `n_pairs`, `n_lines`, `C_within`, `mag`, `ccle_median`, deciles, `null_mean`, the two single-dataset scores, `r_bar`, `C_full4`, `C_rel`, `rel_tercile`) |

Reproduce: `PYTHONDONTWRITEBYTECODE=1 uv run python {inv,build,analyse,rel,desc,extra,tf}279.py` in that order, from this directory.
