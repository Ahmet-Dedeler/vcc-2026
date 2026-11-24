# 289 — Do the source paper's own state features predict Broad-16's line-specific response component?

Date 2026-09-06. **Analysis only**: no build, no `vcc`, no `cell-eval2`, nothing written outside
`experiments/289-broad16-state-features/`, nothing under `data/` or `submissions/`. The four reserved Broad-16 lines
(JVE127, MG63, SLR23, UMRC3) were **not read**: `load268.py`'s assertion is unmodified (`grep -c` = 1) and no reserved
line name appears in any output file of this experiment. Every number is on the **12 working lines** under
leave-one-line-out. Wall clock ≈ 1 h 40 min; peak RSS < 2 GB; 200.9 MB of public DepMap tables fetched into
`work/depmap/` (sha256 in `PREREGISTRATION.md` §2 F-CN).

**Preregistration**: `PREREGISTRATION.md`, sha256
`4d84c53178701fa995d9587e0082f072f5909eed7631ec4bdee21c4f46d4ed39`, hashed **2026-09-06 14:20:21Z**
(`work/PREREGISTRATION.sha256`) — before any state-feature value, any LOLO number, any P1/P2 statistic and any gating
number existed. The only number that predated it is KA-1, a bit-identical re-run of 268's published arm.

---

## 0. Verdict

> **STATE FEATURES PREDICT THE LINE-SPECIFIC RESIDUAL: NO** (S1 predicted-residual cosine **+0.0110 [+0.0060, +0.0162]**
> per gene and **+0.0091 [+0.0015, +0.0169]** programme-level, vs **wrong-recipient +0.0069 / +0.0043** and
> **shuffled +0.0002 / −0.0005**; the paired gap to the wrong-recipient control crosses zero —
> +0.0041 [−0.0004, +0.0088] and +0.0048 [−0.0026, +0.0123] — and the full-response Δcos is
> **−0.0588 [−0.0714, −0.0467]** per gene, **−0.0236 [−0.0302, −0.0172]** programme-level, both intervals below 0).
> **PAPER'S P1 REPRODUCES ON OUR DATA: YES** (P53_PATHWAY score of the line-specific residual after UBE3A/MDM2:
> SKGII +0.1025 and C4I +0.0684 are the top 2 of 11 lines, every other line negative; exact permutation over all 55
> line pairs, the HPV pair ranks **1st, p = 0.018**; the E2F/G2M direction shows the opposite ordering, p = 1.000).
> **P2: NO** (testable — 3 co-deletion lines — and the effect runs the **wrong way**: IER3IP1 dependency −0.159 in the
> co-deleted lines vs −0.476 in the rest, one-sided Mann–Whitney p = 0.93; badly underpowered at n = 12).
> **PER-TARGET GATING beats chance on 17.1 % of targets** on the p53 direction (vs 7.3 % binomial expectation and
> 7.8 % shuffled; pooled sign accuracy 0.559 vs 0.506 shuffled, paired +0.053 [+0.015, +0.089], P(≤0) = 0.003), and on
> 20.7 % on the E2F/G2M direction but there **not separated from its own control** (10.0 % shuffled; paired
> +0.015 [−0.034, +0.064], P(≤0) = 0.28).

**In one sentence.** The source paper's mechanism is *in* our residual object and is recovered blind by the
leave-one-line-out machinery — UBE3A is the single best-predicted target of all 82 and the feature that carries it is
HPV status — but a state-feature model is still, transcriptome-wide, a **worse** transfer prediction than the shared
response, exactly as two basal expression PCs were. **The two statements in 287's contradiction C-1 are both true**;
what was wrong was reading 268's aggregate null as a statement about mechanism.

No model use for the competition is proposed.

---

## 1. Known-answer rows

| KA | what it checks | result |
|---|---|---|
| **KA-1** | `work/ka289.py` = `lolo268.py` with output paths redirected; must reproduce R268 §4 | **PASS at 0.0e+00** — all 788 numeric entries of the summary JSON and all 984 × 45 CSV cells identical to `E268/work/lolo268_K2.{json,csv}`; (c1) cos_res **+0.0181 [+0.0134, +0.0229]**, Δcos_full **−0.0071**, Δretrieval **+0.0003**; ceiling median 0.0455 |
| **KA-2** | the p53 dependency proxy, on all 1,150 DepMap 24Q2 models: TP53 and MDM2 Chronos effects anti-correlated (≤ −0.3) and MDM2-dependent models with higher TP53 effect | **PASS** — Pearson **−0.7408**; mean TP53 effect **+0.8115** in the 400 models with MDM2 ≤ −0.5 vs **+0.1126** in the rest |
| **KA-3** | the generalised fold loop with the 2 basal PCs must reproduce KA-1's (c1)/(c2) to ≤ 1e-12 | **PASS at 0.0e+00** on every readout of `pc2.g`, `pc2.p`, `pc2.b`, `nn`, `oracle` (all 984 rows). One control mean, `nn__WRONG_recipient|dcos_full`, differs at **5.96e-08** — float32 summation order only (268 subtracts the baseline after averaging the 11 wrong recipients, 289 before); outside KA-3's stated scope and reported for the record |
| **KA-4** | hallmark parse | **PASS** — HALLMARK_P53_PATHWAY / E2F_TARGETS / G2M_CHECKPOINT are 200 symbols each; on the 12,199-gene axis 176 p53-pathway and 303 E2F∪G2M genes are present |

---

## 2. The 12 lines' frozen state features (`work/state289.csv`)

Rules from `PREREGISTRATION.md` §2. Depth columns are **nuisance** (rule 9) and entered no model.

| line | lineage | batch | HPV-E6 | TP53 eff. | MDM2 eff. | **p53_wt** (dep.) | p53_score | TP53 dam./hot. | **p53_wt** (mut.) | e2f_g2m | p53_path | cn 18q21 (log2) | cn 11q23-24 (log2) | **codel** | UMIs/cell | UMIs/KO | cells/KO |
|---|---|---|---:|---:|---:|:--:|---:|---:|:--:|---:|---:|---:|---:|:--:|---:|---:|---:|
| C4I | Cervix | apollo33 | **1** | +0.268 | −0.315 | 0 | 0.584 | 0 / 0 | 1 | 0.542 | 0.228 | −0.987 | −0.115 | 0 | 15,055 | 1.04 M | 68.9 |
| HCT15 | Bowel | apollo31 | 0 | +0.008 | −0.391 | 0 | 0.399 | 2 / 1 | 0 | 0.719 | 0.015 | −0.004 | −0.010 | 0 | 47,408 | 2.56 M | 53.5 |
| HKGZCC | Biliary Tract | apollo38 | 0 | +0.232 | −0.463 | 0 | 0.696 | 2 / 2 | 0 | 0.657 | 0.253 | −0.727 | +0.500 | 0 | 32,599 | 3.88 M | 122.7 |
| KMRC20 | Kidney | apollo35 | 0 | +0.295 | −0.560 | **1** | 0.855 | 2 / 0 | 0 | 0.638 | 0.078 | −0.717 | −0.172 | 0 | 12,503 | 1.34 M | 108.2 |
| KYSE140 | Esophagus/Stomach | apollo31 | 0 | +0.159 | −0.133 | 0 | 0.292 | 2 / 2 | 0 | 0.676 | 0.027 | −0.649 | −0.545 | **1** | 45,577 | 3.58 M | 79.0 |
| KYSE450 | Esophagus/Stomach | apollo31 | 0 | +0.298 | −0.326 | 0 | 0.625 | 1 / 1 | 0 | 0.587 | −0.002 | −0.978 | −0.576 | **1** | 35,964 | 6.05 M | 164.3 |
| ONE58 | Pleura | apollo34 | 0 | +0.245 | −0.417 | 0 | 0.662 | 2 / 2 | 0 | 0.758 | 0.029 | −1.403 | +0.079 | 0 | 20,071 | 2.79 M | 141.4 |
| RPMI7951 | Skin | apollo31 | 0 | +0.121 | −0.343 | 0 | 0.464 | 2 / 0 | 0 | 0.540 | 0.051 | −0.119 | +0.218 | 0 | 51,862 | 3.44 M | 67.5 |
| SKGII | Cervix | apollo32 | **1** | +0.266 | −0.439 | 0 | 0.705 | 0 / 0 | 1 | 0.612 | 0.126 | −0.769 | −0.722 | **1** | 75,928 | 4.81 M | 63.3 |
| SNU8 | Ovary/Fallopian Tube | apollo35 | 0 | +0.127 | −0.376 | 0 | 0.504 | 2 / 0 | 0 | 0.482 | 0.092 | −0.373 | −0.168 | 0 | 11,357 | 1.17 M | 102.3 |
| T3M5 | Head and Neck | apollo34 | 0 | +0.537 | −0.444 | 0 | 0.981 | 2 / 0 | 0 | 0.429 | 0.237 | −0.634 | −0.033 | 0 | 29,515 | 1.80 M | 60.9 |
| U343 | CNS/Brain | apollo32 | 0 | +0.754 | −1.407 | **1** | 2.160 | 0 / 0 | 1 | 0.416 | 0.050 | −0.610 | −0.090 | 0 | 29,100 | 1.48 M | 51.7 |

- **p53 (dependency rule)**: `MDM2 ≤ −0.5` **and** `TP53 ≥ 0` ⇒ wild-type. **KMRC20 and U343** only (2 wt / 10 def) — not
  degenerate by the frozen contingency (≥ 2 per class), so the binary stands.
- **p53 (mutation annotation, agreement check)**: no damaging and no hotspot TP53 mutation ⇒ **C4I, SKGII, U343** (3/12).
  Agreement of the two calls **0.75 (9/12)**. The three disagreements are informative rather than noise: C4I and SKGII
  are TP53 **wild-type by sequence** and p53-**deficient by function** — precisely the HPV-E6 mechanism, where the viral
  oncoprotein degrades p53 protein without a mutation — and KMRC20 carries a damaging TP53 call yet retains an
  MDM2 dependency. Neither annotation is a gold standard here; both were carried as separate arms (S1 vs S4).
- **HPV-E6**: SKGII, C4I (the paper, quoted through E287 row 10). **Declared confound** (prereg §2): these are the only
  two Cervix lines among the 12, so `hpv` and "Cervix lineage" are the same column on this data.
- **Copy number**: `cn` = log2 of the length-weighted mean relative segment CN over the region (DepMap 24Q2 relative
  segments; 18q21 = 197 genes, 45.80–63.90 Mb; 11q23-24 = 390 genes, 110.58–130.87 Mb). `codel` (both regions
  < 0.85 relative) = **KYSE140, KYSE450, SKGII**.
- **Depth spans 6.7× in UMIs per cell and 5.8× in UMIs per perturbation**; four lines sit below the source paper's
  ≈ 5 M UMI-per-perturbation elbow by more than 3×.

---

## 3. The LOLO table (82 targets × 12 folds, paired target bootstrap, 4,000 draws)

Everything is refit inside the 11 training lines. `.g` = per-gene ridge (268's c1 form), `.p` = programme-level, J = 5
(268's c2 form). Split-half **ceiling**: median 0.0455, mean 0.0607.

| arm | features | cos(pred. resid., held resid.) | − WRONG recipient | − SHUFFLED labels | Δcos(full response) | folds cos>0 |
|---|---|---|---|---|---|---|
| (a) target-only | — | 0 | — | — | 0 | — |
| **(c1) `pc2.g`** | 2 basal PCs | +0.0181 [+0.0134, +0.0229] | +0.0169 [+0.0122, +0.0218] **P 0.000** | +0.0241 [+0.0186, +0.0297] **P 0.000** | −0.0071 [−0.0093, −0.0049] | 9/12 |
| **(c2) `pc2.p`** | 2 basal PCs | +0.0174 [+0.0108, +0.0239] | +0.0167 [+0.0099, +0.0233] **P 0.000** | +0.0206 [+0.0130, +0.0282] **P 0.000** | +0.0000 [−0.0016, +0.0017] | 8/12 |
| **S1.g** (headline) | p53, HPV, E2F/G2M, p53-path, 18q21, 11q23-24 | +0.0110 [+0.0060, +0.0162] | **+0.0041 [−0.0004, +0.0088] P 0.038** | +0.0107 [+0.0054, +0.0161] P 0.000 | **−0.0588 [−0.0714, −0.0467]** | 6/12 |
| **S1.p** (headline) | same | +0.0091 [+0.0015, +0.0169] | **+0.0048 [−0.0026, +0.0123] P 0.103** | +0.0096 [+0.0017, +0.0176] P 0.010 | **−0.0236 [−0.0302, −0.0172]** | 8/12 |
| S1n.g | S1 without CN | +0.0227 [+0.0172, +0.0289] | +0.0092 [+0.0045, +0.0141] P 0.000 | +0.0230 [+0.0170, +0.0301] P 0.000 | −0.0445 [−0.0540, −0.0352] | 7/12 |
| S1n.p | S1 without CN | +0.0113 [+0.0036, +0.0195] | +0.0046 [−0.0033, +0.0123] P 0.127 | +0.0115 [+0.0031, +0.0203] P 0.003 | −0.0177 [−0.0233, −0.0121] | 9/12 |
| S3.g | S1 + 2 basal PCs | +0.0150 [+0.0099, +0.0202] | +0.0065 [+0.0019, +0.0111] P 0.004 | +0.0228 [+0.0172, +0.0285] P 0.000 | −0.0580 [−0.0710, −0.0459] | 6/12 |
| S3.p | S1 + 2 basal PCs | +0.0119 [+0.0041, +0.0198] | +0.0066 [−0.0009, +0.0138] P 0.042 | +0.0157 [+0.0075, +0.0241] P 0.000 | −0.0239 [−0.0309, −0.0172] | 9/12 |
| S4.g | S1, p53 by mutation | +0.0105 [+0.0052, +0.0165] | +0.0008 [−0.0039, +0.0058] P 0.388 | +0.0133 [+0.0073, +0.0196] P 0.000 | −0.0516 [−0.0629, −0.0406] | 7/12 |
| S4.p | S1, p53 by mutation | +0.0103 [+0.0019, +0.0189] | +0.0049 [−0.0031, +0.0131] P 0.108 | +0.0073 [−0.0014, +0.0161] P 0.046 | −0.0190 [−0.0248, −0.0130] | 9/12 |
| (nn) nearest basal line | — | +0.0353 [+0.0305, +0.0400] | +0.0241 (vs +0.0112) | — | −0.0494 [−0.0603, −0.0391] | 8/12 |
| **(nnL) nearest same-lineage line** | — | +0.0439 [+0.0345, +0.0531] | — | — | −0.0526 [−0.0646, −0.0404] | 3/**4** folds |
| oracle line offset (held data) | — | +0.2968 [+0.2910, +0.3027] | — | — | +0.1690 [+0.1574, +0.1806] | 12/12 |

**Single features (S2), the screen — 8 features × 2 forms, read with the multiplicity of §3.1:**

| feature | `.g` cos_res | − WRONG | − SHUF | `.g` Δcos(full) | `.p` cos_res | − WRONG |
|---|---|---|---|---|---|---|
| **p53_wt** | **+0.0390 [+0.0318, +0.0468]** | +0.0088 [+0.0061, +0.0115] P 0.000 | +0.0388 P 0.000 | −0.0117 [−0.0150, −0.0087] | +0.0226 [+0.0159, +0.0300] | +0.0076 [+0.0029, +0.0123] P 0.001 |
| p53_score | +0.0283 [+0.0236, +0.0330] | +0.0068 [+0.0041, +0.0094] P 0.000 | +0.0388 P 0.000 | −0.0121 [−0.0158, −0.0086] | +0.0145 [+0.0086, +0.0207] | +0.0020 [−0.0027, +0.0067] P 0.199 |
| **hpv** | +0.0216 [+0.0143, +0.0310] | +0.0059 [+0.0022, +0.0113] P 0.000 | +0.0241 P 0.000 | −0.0120 [−0.0149, −0.0090] | +0.0019 [−0.0051, +0.0098] | +0.0031 [−0.0021, +0.0092] P 0.137 |
| e2f_g2m | +0.0212 [+0.0168, +0.0254] | +0.0147 [+0.0106, +0.0188] P 0.000 | +0.0187 P 0.000 | −0.0166 [−0.0201, −0.0132] | +0.0041 [−0.0028, +0.0113] | +0.0016 [−0.0053, +0.0088] P 0.334 |
| p53_path | +0.0192 [+0.0134, +0.0254] | +0.0079 [+0.0035, +0.0128] P 0.000 | +0.0257 P 0.000 | −0.0135 [−0.0170, −0.0100] | +0.0076 [−0.0006, +0.0165] | +0.0057 [−0.0017, +0.0137] P 0.077 |
| cn_18q21 | −0.0012 [−0.0050, +0.0025] | −0.0008 P 0.681 | +0.0040 P 0.022 | −0.0145 [−0.0179, −0.0112] | +0.0007 [−0.0068, +0.0080] | +0.0001 P 0.491 |
| cn_11q23_24 | −0.0098 [−0.0136, −0.0057] | −0.0042 P 0.975 | −0.0126 P 1.000 | −0.0164 [−0.0196, −0.0130] | +0.0023 [−0.0047, +0.0098] | +0.0047 P 0.092 |
| codel | −0.0360 [−0.0400, −0.0321] | +0.0057 P 0.003 | −0.0392 P 1.000 | −0.0169 [−0.0201, −0.0137] | −0.0178 [−0.0238, −0.0119] | −0.0008 P 0.585 |

### 3.1 Reading

1. **The headline arm fails the frozen pass rule on both conditions.** S1's predicted-residual cosine is above the
   shuffled-label control on both forms, but its paired gap to the **wrong-recipient** control includes zero
   (+0.0041 [−0.0004, +0.0088] per gene; +0.0048 [−0.0026, +0.0123] programme-level), and the full-response Δcos is
   **negative with the interval below zero** on the programme-level readout (−0.0236). Condition (i) fails, condition
   (ii) fails. **NO.**
2. **The wrong-recipient control is the discriminating one, and it is much larger for state features than for basal
   PCs** (+0.0069 vs +0.0012 on the `.g` form). Handing the model *any* line's HPV/p53/CN vector recovers most of the
   apparent cosine, because the frozen features are low-cardinality (two are binary; ten of twelve lines share the
   same p53 call, ten share the same HPV call), so a wrong recipient usually carries the *same* value as the right
   one. 268's continuous PCs do not have this property. **A state-feature model must be scored against the
   wrong-recipient control specifically; the shuffled-label control alone would have passed every state arm here.**
3. **Adding the paper's features to 268's PCs makes the prediction worse, not better** (S3 vs pc2: cos_res
   +0.0119 vs +0.0174 programme-level, Δcos −0.0239 vs +0.0000). With 11 training rows, six extra columns cost more in
   variance than they buy in signal.
4. **The CN features are inert or harmful** on every readout — the only two features whose single-feature arm is at or
   below its own controls, and `codel` is strongly negative (−0.0360, below its shuffled control at P 1.000).
5. **`p53_wt` alone is the strongest single state feature** and clears both controls in both forms
   (+0.0390, gap to wrong-recipient +0.0088 P 0.000; +0.0226, gap +0.0076 P 0.001) — larger than 268's two basal PCs
   (+0.0181). It still buys **nothing** on the full response (−0.0117 / −0.0043, both intervals below 0), so it fails
   condition (ii) like everything else. This is a **cell of a 25-cell preregistered screen** (12 feature sets × 2 forms
   + nnL; the 13 recipient-offset `.b` forms were additionally computed and are descriptive, not in that count), and
   its per-cell interval ignores the screen.
6. **Nearest same-lineage line** has the highest non-oracle residual cosine of any arm (+0.0439) — but on only the
   **4 folds** that have a same-lineage partner among the 12 working lines (C4I↔SKGII, KYSE140↔KYSE450; Kidney and
   Bowel lost their partners to the reserve), and it damages the full response by −0.0526. Same-lineage similarity is
   real and is not a usable transfer prediction; this is 268's (nn) result with a sharper label.
7. **The oracle is unchanged** (+0.2968 residual cosine, +0.1690 on the full response): the headroom exists, but it is
   the held line's *own* mean residual, which no state feature reaches.

---

## 4. P1 — the paper's p53/HPV prediction: **REPRODUCES**

Frozen rule, prereg §5. Full-data responses, residual against the included-line mean, HALLMARK_P53_PATHWAY score
minus the all-gene mean, 176 p53-pathway genes on the 12,199-gene axis. MDM2 lost **U343** (5 cells < 6, which is why
MDM2 is not in 268's 82); **n = 11 common lines**. Cells per line: UBE3A 21–192, MDM2 26–114.

| line | HPV | p53 (dep.) | **P53_PATHWAY score of the residual, mean(UBE3A, MDM2)** | UBE3A alone | MDM2 alone | E2F/G2M score (specificity) |
|---|:--:|:--:|---:|---:|---:|---:|
| **SKGII** | **1** | 0 | **+0.1025** | **+0.1906** | +0.0144 | −0.0495 |
| **C4I** | **1** | 0 | **+0.0684** | **+0.1385** | −0.0017 | −0.0523 |
| HCT15 | 0 | 0 | −0.0036 | −0.0187 | +0.0115 | +0.0185 |
| KYSE140 | 0 | 0 | −0.0118 | −0.0288 | +0.0052 | −0.0073 |
| HKGZCC | 0 | 0 | −0.0146 | −0.0305 | +0.0012 | −0.0160 |
| KMRC20 | 0 | **1** | −0.0151 | −0.0263 | −0.0040 | +0.0207 |
| SNU8 | 0 | 0 | −0.0177 | −0.0356 | +0.0002 | +0.0076 |
| ONE58 | 0 | 0 | −0.0235 | −0.0447 | −0.0023 | +0.0236 |
| KYSE450 | 0 | 0 | −0.0257 | −0.0501 | −0.0013 | +0.0109 |
| RPMI7951 | 0 | 0 | −0.0299 | −0.0516 | −0.0082 | −0.0087 |
| T3M5 | 0 | 0 | −0.0367 | −0.0584 | −0.0150 | +0.0482 |

- **Both orderings hold**: HPV mean **+0.0854** > p53-wt-non-HPV (KMRC20) **−0.0151** > p53-deficient mean **−0.0204**.
- **Exact permutation over all 55 line pairs: the observed {SKGII, C4I} pair ranks 1st, one-sided p = 0.0182.**
  **P1 PASSES its frozen rule.**
- **Driven by UBE3A** (HPV mean +0.1645, p = 0.0182, pair ranks 1st). **MDM2 alone does not reach it**
  (HPV mean +0.0064, p = 0.127; SKGII is 1st but C4I is 7th of 11) — on our pseudobulks, at 26–114 cells per line,
  the MDM2 leg of the paper's claim is not resolved.
- **Specificity check passes**: on the E2F/G2M direction the two HPV lines are the **bottom two** of 11 (p = 1.000).
  The effect is p53-directional, not a general "the cervical lines respond more" offset.
- **Depth conditioning (rule 9).** The P1 score is uncorrelated with depth (Spearman with UMIs per perturbation
  −0.009, p = 0.98; with UMIs per cell +0.191; with cells per KO −0.245), and the two HPV lines sit at **opposite
  ends of the depth range** (C4I is the shallowest of the 12 at 1.04 M UMIs/KO, SKGII the 11th at 4.81 M). *Post-hoc,
  labelled*: restricting the permutation null to the 9 line pairs whose depth-rank gap is ≥ 8 leaves the HPV pair
  ranked 1st of 9, p = 0.111 — the smallest value that null can produce.
- **What P1 does and does not license.** It confirms that HPV-E6 status predicts the *direction* of a specific
  target's line-specific residual, in our own residual object, on the paper's own lines. It does **not** separate HPV
  status from cervical lineage (§2), and the two positive lines are two lines.

### 4.1 The LOLO finds the same mechanism blind

Not preregistered as a test, but it is the same machinery and it is decisive for the interpretation. Ranking all 82
targets by S1's predicted-residual cosine:

| rank | target | S1.g cos_res | − own wrong-recipient | S1.p cos_res | pc2.g cos_res | S1.g Δcos(full) | ceiling | best single feature (cos_res) |
|---:|---|---:|---:|---:|---:|---:|---:|---|
| **1** | **UBE3A** | **+0.1153** | **+0.0717** | +0.1027 | +0.0082 | −0.0194 | 0.0660 | **hpv (+0.2908)** |
| 2 | HDAC3 | +0.0616 | +0.0556 | +0.0981 | +0.0418 | −0.1160 | 0.1386 | p53_score (+0.0787) |
| 3 | MYC | +0.0546 | +0.0478 | +0.0723 | +0.0852 | −0.1762 | 0.1745 | p53_path (+0.1031) |
| 4 | IER3IP1 | +0.0520 | +0.0419 | +0.0654 | −0.0004 | −0.0129 | 0.0348 | hpv (+0.0661) |
| 5 | USP37 | +0.0488 | +0.0087 | +0.0224 | +0.0168 | −0.0286 | 0.1010 | p53_wt (+0.1955) |
| 6 | POLR2D | +0.0483 | +0.0316 | +0.0558 | +0.0822 | −0.1562 | 0.1152 | p53_wt (+0.0844) |
| 7 | SMG6 | +0.0472 | +0.0376 | +0.0787 | +0.0345 | −0.1807 | 0.1425 | p53_path (+0.0936) |
| 8 | SLC39A9 | +0.0462 | +0.0439 | +0.0772 | +0.0148 | +0.0108 | 0.0380 | p53_path (+0.0447) |
| 9 | PRKRA | +0.0406 | +0.0361 | +0.0666 | +0.0234 | −0.0082 | 0.0423 | p53_path (+0.0423) |
| 10 | SCYL1 | +0.0390 | +0.0276 | +0.0655 | −0.0064 | −0.0143 | 0.0145 | p53_path (+0.1136) |

(full table `work/descriptive289.csv`; ranks 11–20 UBA5, WDR48, IBA57, ADAR, TXN, C3orf38, NHLRC2, UBE2H, N6AMT1,
PITRM1, all +0.024…+0.033.)

**UBE3A is the best-predicted target of all 82, its carrier is HPV status, and the HPV arm predicts it on 12/12 folds**
(per-fold cosine +0.178…+0.402, including **+0.287 on the held-out C4I fold and +0.241 on the held-out SKGII fold**,
where the coefficient was estimated from the single remaining HPV line). Its wrong-recipient control is high
(+0.1211 — see §3.1 reading 2) and its shuffled control is −0.0247, so the recipient-specific part is +0.170. The
mechanism the paper reports is therefore not merely present in our data: it is **the strongest thing a
state-conditioned LOLO finds on Broad-16**, and it is a single target out of 82.

Best-carrier counts over the 82 targets: p53_wt 36, p53_path 14, e2f_g2m 13, hpv 10, p53_score 9, **CN features 0**.

---

## 5. P2 — the IER3IP1 / 18q–11q co-deletion prediction: **DOES NOT REPRODUCE**

Testable here (3 co-deletion lines ≥ 2), and it fails in the predicted direction on every readout.

| readout | co-deleted (KYSE140, KYSE450, SKGII) | the other 9 | one-sided MW p | Spearman vs continuous co-deletion score |
|---|---:|---:|---:|---|
| DepMap Chronos IER3IP1 gene effect (more negative = stronger) | **−0.159** | **−0.476** | **0.927** | ρ −0.035, p 0.914 |
| IER3IP1 response norm / line median response norm | **0.784** | **0.908** | 0.895 | ρ +0.147, p 0.649 |
| IER3IP1 residual norm / line median residual norm | 0.809 | 0.904 | 0.927 | ρ +0.294, p 0.354 |

Per line, the strongest IER3IP1 dependencies are **C4I (−1.198), HKGZCC (−0.893), ONE58 (−0.651), U343 (−0.533)** —
none of them co-deleted by our call; the three co-deleted lines are −0.027, −0.249, −0.202. **P2 = NO.**

Honest scope, in order of importance:
1. **Power.** n = 12 lines, 3 vs 9; the smallest attainable one-sided Mann–Whitney p is 0.0045, and the paper's own
   p = 3.28 × 10⁻⁶ comes from a genome-scale comparison across hundreds of DepMap models, not from these 12. This is
   a failure to reproduce **on 12 lines**, not a refutation of the paper's association.
2. **Our co-deletion call is a proxy.** The paper writes "hemizygous 18q21.1-32 / 11q23-4"; we froze
   band-prefix regions (18q21·*, 11q23·*∪11q24·*) and a length-weighted relative-CN threshold of 0.85 before looking.
   A different region or threshold could give a different 3 lines.
3. **The direction of the miss is not subtle** (co-deleted lines are *less* IER3IP1-dependent by 0.32 Chronos units),
   which argues against a marginal-threshold explanation and for a genuine absence of the association in this
   12-line slice.
4. Notably, IER3IP1 *is* the 4th best-predicted target under S1 (§4.1) — but its carrier is **hpv**, and its hpv arm
   sits **at its own wrong-recipient control** (+0.0661 vs +0.0686), i.e. that ranking is not recipient-specific
   information.

---

## 6. Per-target gating (prereg §7): state predicts the **sign** of the p53 projection above chance

Fixed hallmark directions on each fold's gene axis; ridge on (hpv, p53_wt) for the p53 direction and on the basal
E2F/G2M score for the E2F direction; sign of the held line's projection; 5-permutation shuffled control.

| direction | pooled sign accuracy | shuffled | paired (target bootstrap) | pooled binomial p | targets ≥ 9/12 | shuffled ≥ 9/12 | binomial expectation |
|---|---:|---:|---|---:|---:|---:|---:|
| **P53_PATHWAY** | **0.5589** | 0.5063 | **+0.0526 [+0.0146, +0.0888], P(≤0) 0.0028** | **1.2 × 10⁻⁴** | **17.1 %** (14/82) | 7.8 % | 7.3 % |
| E2F ∪ G2M | 0.5264 | 0.5116 | +0.0148 [−0.0342, +0.0640], P(≤0) 0.28 | 0.052 | 20.7 % (17/82) | 10.0 % | 7.3 % |

**Reading.** Target-by-target, two binary state features predict *which way* a held line's residual moves along the
p53 axis better than chance and better than their own shuffled control — 55.9 % against 50.6 %, and 2.3× the expected
number of targets clearing 9/12. This is the smallest true statement that survives: **state predicts a sign, on one
programme, at an effect size of five percentage points.** It is not the same object as a transfer prediction — the
same features, given a magnitude and a direction to write into a response vector, make the full-response cosine
*worse* (§3). The E2F/G2M direction shows a similar target count at the 9/12 bar but is **not separated from its own
control**, so its 20.7 % is a multiplicity artefact of 82 correlated targets, not signal.

---

## 7. Depth as a covariate (programme rule 9) — the negative is not a shallow-data artefact

`work/folds289_depth.csv`, `work/depth289.json`; n = 12 folds, so every Spearman below is a 12-point statistic.

- **The split-half ceiling tracks depth**, as the source paper's UMI–power law says it should: Spearman(ceiling median,
  UMIs per perturbation) = **+0.594, p = 0.042**. Rule 9's premise is visible in our own folds.
- **S1's residual cosine is not concentrated in the deep lines**: Spearman with UMIs per perturbation −0.063 (`.g`) and
  +0.294 (`.p`), both p > 0.35. Median by depth half: `.g` +0.0078 deep vs −0.0005 shallow; `.p` +0.0107 vs −0.0003.
- **S1's damage to the full response is *worse* in the deep lines**: Spearman(Δcos_full, UMIs per perturbation)
  **−0.720, p = 0.008** (`.g`) and **−0.587, p = 0.045** (`.p`). Median Δcos: −0.0743 deep vs −0.0538 shallow (`.g`).
- Therefore **the failure of condition (ii) cannot be explained away as noise in shallow lines**: where the truth is
  best resolved, the state-conditioned prediction hurts most. (Mechanically this is expected — a deeper held line has a
  larger, better-measured true response, so a mis-directed additive correction of fixed relative size costs more
  cosine — but it removes the "we'd see it with more UMIs" reading of this particular negative.)
- The two HPV lines straddle the depth range (§4), so P1 is not a depth effect either.

---

## 8. What this closes, and what it does not

**Closes (at this construction and this n).** A linear model on the *source paper's own state features* —
p53 functional status (by DepMap dependency or by mutation), HPV-E6 status, basal E2F/G2M and P53_PATHWAY programme
scores, and arm-level 18q21 / 11q23-24 copy number — fitted on 11 Broad-16 lines and predicting the held line's
line-specific residual over 82 targets and 12,199 genes, **does not beat the shared response**, and does not beat the
wrong-recipient control on the residual either. Adding those features to 268's basal PCs makes it worse. This extends
F5's amended headline from "not predictable from basal **expression** by a linear model at 11 training lines" to
**"not predictable, by a linear transcriptome-wide model at 11 training lines, from the paper's own genotype, viral
and programme features either"**. 25 preregistered arm-cells; 0 pass.

**Does not close, and 287's C-1 is now settled the other way round.** The paper's mechanistic claims are **not** in
conflict with 268's null, and one of them is confirmed here: **P1 reproduces in our residual object** (p = 0.018,
with a passing specificity check), **the LOLO recovers it blind** (UBE3A is the best-predicted of 82 targets and HPV
is its carrier, 12/12 folds), and **per-target sign gating on the p53 axis beats chance and its shuffled control**.
The correct statement is therefore:

> a state feature can predict *which way a particular target's response bends in a particular line* — and
> transcriptome-wide, target-averaged, additive prediction of the whole residual from 11 training lines still buys
> nothing. The two live at different granularities, and 268's aggregate readout could never have detected the
> mechanism it was said to contradict: **UBE3A is 1 of 82 targets, and the P53_PATHWAY genes are 176 of 12,199.**

**Recommended amendment to `docs/FINDINGS-2026-09.md` F5** (for the coordinator; not applied here — this experiment
writes only inside its own directory): keep F5's negative with its scope, and add that 289 tested the paper's own
state features under the same protocol with the same result, **while confirming the paper's P1 on our data** and
finding that the mechanism is visible per target (UBE3A/HPV) and per direction (p53 sign gating, 17.1 % of targets)
but not in the transcriptome-wide aggregate. C-1 should be re-graded from CONTRADICTED to **RECONCILED**: both claims
hold, at different granularities. P2 did not reproduce on these 12 lines and is underpowered by construction.

**What a stage 2 would have to be** (not proposed as work, stated as scope): the object that has signal here is a
**per-target, per-programme sign/gate**, not a transcriptome-wide additive residual model — and its effect size
(+0.05 accuracy, one programme) is far below what would change a transfer prediction. Buying more *lines* is what
rule 7 would say next (11 training rows cannot support 6 features), not more capacity.

---

## 9. Honest limits

1. **n = 12 lines, 11 training rows.** Every LOLO interval is a target bootstrap; the *line* is the unit that limits
   what can be learned, and it is not resampled. Six features on 11 rows is the regime where S3 loses to pc2.
2. **HPV-E6 is confounded with cervical lineage** on these 12 lines (§2), and with nothing else measured. P1's
   two positive lines are two lines, of one lineage, in two different batches (apollo32, apollo33) — which is the one
   thing that keeps it from being a batch effect, not a proof that it is not one.
3. **268's limits carry over**: split halves inside one experiment reproduce batch and line-level technical structure
   as well as biology (268 §7.1, batch permutation p = 0.001); cells per (line, KO) are small so the per-(t,l)
   residual is barely reproducible (ceiling median 0.0455); the five most lethal targets are excluded by the 6-cell
   rule, and MDM2 is one of them — which is why P1's MDM2 leg rests on 11 lines at 26–114 cells.
4. **The p53 feature is a proxy in both forms.** Neither the DepMap dependency rule nor the mutation annotation is
   ground truth for p53 *function*; they disagree on 3 of 12 lines, and the disagreements are mechanistically
   interpretable (§2) rather than random, which means neither call is neutral with respect to P1.
5. **The CN region definition was frozen from a second-hand quotation** of the paper's bands ("18q21.1-32 / 11q23-4",
   via E287 row 10; the preprint was not re-read). A different region or a different threshold gives a different
   `codel` set, and P2's failure is conditional on ours.
6. **Multiplicity**: 25 preregistered arm-cells × 3 readouts, plus 13 descriptive `.b` cells; the S2 single-feature
   rows are a screen and their per-cell intervals ignore it. §4.1's ranking and the depth-matched P1 permutation are
   labelled post-hoc.
7. **Zero panel overlap** between Broad-16's 90 targets and the 300-target VCC panel (268 §5): nothing here is a
   measured statement about the competition, and none is made.

---

## 10. Files

| file | what |
|---|---|
| `PREREGISTRATION.md` | the frozen document (sha256 `4d84c531…`, 14:20:21Z) |
| `work/PREREGISTRATION.sha256` | the hash + timestamp |
| `work/ka289.py`, `work/ka289_K2.log`, `work/ka289_lolo268_K2.{csv,json}` | KA-1 (268's arm re-run, 0.0e+00) |
| `work/features289.py` → `work/state289.csv`, `work/features289.json` | the frozen state features, KA-2, KA-4, CN regions, p53 call agreement |
| `work/lolo289.py` → `work/lolo289.{csv,json}`, `work/folds289.csv`, `work/lolo289.log` | the state-feature LOLO, KA-3, all arms and controls |
| `work/p1p2_289.py` → `work/p1p2_289.json`, `work/p1_scores289.csv`, `work/gating289b.csv` | P1, P2, gating with per-permutation shuffles |
| `work/descriptive289.csv`, `work/folds289_depth.csv`, `work/depth289.json` | §4.1 tail, depth covariate |
| `work/depmap/` (200.9 MB) | DepMap 24Q2 Public segment + mutation tables, figshare 25880521; sha256 in `PREREGISTRATION.md` §2 |
| `work/gating289.csv` | the in-loop gating pass (superseded by `gating289b.csv`, kept for the record) |

Read-only inputs outside this directory: `E268/work/{halves/*.npz, inventory.json, reserve.json, load268.py,
lolo268.py, lolo268_K2.*}`, `E269/work/depmap/CRISPRGeneEffect.csv`,
`data/external/broad16line/{metadata/cell_line_metadata.csv, metadata/hallmark_symbols.gmt, data/gene_locations.csv}`,
`data/external/depmap26q1/extracted/processed/cell_quality_covariates.csv`. Nothing outside
`experiments/289-broad16-state-features/` was written.
