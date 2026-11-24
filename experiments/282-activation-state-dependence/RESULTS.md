# 282 — State dependence within one individual: **the preregistered across-state statistic is a construct artefact; with the guide broken, a CD4 perturbation response keeps only ~half of its across-donor agreement when the activation state changes**

**The line the mandate asks for.**

> **STATE DEPENDENCE WITHIN INDIVIDUALS: median SD = +0.49 [+0.33, +0.61]**
> (deconfounded index, 74 targets, ~107 cells/side, full depth; depth-matched +0.43 [+0.28, +0.53];
> pooled ratio +0.385). **The preregistered index reads −2.68 [−3.56, −2.07] and is an ARTEFACT** — it
> says a response agrees across states 5.6× better than across donors, and the crossed control shows that
> agreement belongs to the **(donor × guide) construct**, not to the gene.
> **PREDICTED BY CROSS-LINE CONSERVATION: NO** — preregistered ρ(SD, `C_noCD4`) = **+0.262
> [+0.120, +0.399]**, interval entirely on the *wrong side of zero* (the frozen prediction was ρ < 0), and
> on the deconfounded index ρ = **−0.144 [−0.403, +0.126]**, which covers 0. What conservation predicts is
> the across-**donor** ceiling (ρ +0.301 [+0.177, +0.419]), not across-**state** agreement (ρ +0.084
> [−0.045, +0.209]) — S1 is flat across conservation terciles (+0.258 / +0.241 / +0.250) while S2 rises
> 4× (+0.017 / +0.046 / +0.067).

**And the finding that outlives the mandate.**

> **The reproducible part of a small CRISPRi pseudobulk in primary CD4 T cells is a *donor × guide*
> interaction, not a gene effect.** At ~66 cells/side and matched depth: same donor + same guide, state
> changed **+0.0424**; same donor, *other guide for the same gene*, same state **+0.0116**; same guide,
> different donor, same state **+0.0141**; change any two **+0.008 … +0.011**. Only the cell that holds
> **both** the person and the construct is large, and the two guides of one gene in one person agree no
> better than two different people. **Any replicate design that holds the construct and the donor fixed —
> two culture arms, two timepoints, two libraries of one transduced pool — overestimates response
> reproducibility by 3.7× (depth-matched) to 45× (full depth).** The three GWCD4i states are twelve independent libraries (one per donor ×
> condition), so this is not a shared-library artefact.

**Date** 2026-09-06, 07:57–08:25 UTC (10:57–11:25 Ankara). **Preregistration** `PREREGISTRATION.md`, sha256
`7d75116b8c5d17943558b484f0e1f07a52a851671109661580bbad5755ead13f`, hashed **08:03:50Z**
(`work/prereg282.sha256`) — before any across-state cosine, index, contrast or decomposition existed; only
`.obs` design metadata, HGNC groups and column names had been read. **No download, no new extraction, no
model, no build, no `vcc`, no `cell-eval2`; nothing written outside this experiment directory.**
`uv run python`, `PYTHONDONTWRITEBYTECODE=1`. Peak RSS **4.53 GB** (`/usr/bin/time -l` on `run282.py`),
cap 8 GB. Total compute ≈ 35 s across five scripts; 280's 1.2 GB bank read memmapped, nothing new on disk
above 1 MB.

---

## 0. Known-answer rows

| row | required | measured | verdict |
|---|---|---|---|
| **KA-1 (the ceiling reproduces 280)** | 280's §2.1 across-donor medians and intervals, its filter, its per-condition target sets, its code path. Bar max\|Δ\| ≤ 5e-5 | compared against **280's own full-precision table** (`q1_halves_keepDE.csv`, not R280's 4-d.p. print): 12/12 rows, **max\|Δ\| = 9.0e-17** (≈ 1 ulp), every `n` identical, Stim48hr POOLED `0.0538943254` and its interval `[0.0426592578, 0.0619975880]` reproduced to the last digit | **PASS** (floating-point identity) |
| **KA-2 (extraction integrity)** | 280's KA-0 on a 4,000-row sample: row gene-sum == `obs.total_counts` | **max\|Δ\| = 0** | **PASS** |
| **KA-3 (conservation join)** | every analysis target carries a `conservation276.csv` row; `C_raw` re-derived from `cons276.npz`'s pair matrix matches the published column (bar 1e-6) | 241/241 targets present; **max\|Δ\| = 4.98e-07**; 69 of 315 source pairs dropped as CD4-derived (`cd4_rest`, `cd4_stim8hr`, `cd4_stim48hr`) | **PASS** |
| **KA-4 (null)** | target-shuffled across-state pairing within ±0.02 of 0 | mean **+0.00022**, median −0.00048, sd 0.0402, against a real S1 median of **+0.2480** | **PASS** |
| **KA-5 (design)** | pooled 4-donor cells per (target, state) reconciles with 280 §1.4 / KA-5 | state-common guides **536 / 537 / 587** (n = 241); 280's filter 537 / 536 / 587 (n = 242); all guides 518 / 502 / 523 (n = 273) on the **common** set vs R280's 500 / 476 / 481 on its own per-condition sets | **PASS** |

Everything below was computed after these five rows printed.

---

## 1. The preregistered result, exactly as frozen

Objects, geometry and arithmetic are 280's, reused by import (`lib280.Bank`): `p = log1p(5e4·S/ΣS)`, control =
the same transform of the side's own donors' NTC pool in the same state, `delta = p − b`, residual =
`delta` minus its own across-target mean over the analysis set, panel positions zeroed, cosine on the
17,471 scored columns. Sides are **donor pairs** (~215 cells). Analysis set: the **241** panel targets whose
guides pass `keep_for_DE` in all three states, in all four donors (the state-common guide filter; §2.4 of
the preregistration explains why, and the two sensitivity arms are in the table).

### 1.1 The 2×2 (241 targets, ~215 cells/side)

| cell | what changes | median | 95 % (paired target bootstrap) | ≤ 0 |
|---|---|---:|---|---:|
| **W** counting floor (same donors, same guides, same state, disjoint count halves) — **half depth**, context only | nothing | **+0.4690** | [+0.4522, +0.4866] | 0 % |
| **S1** across **STATE**, donors + guides held | state | **+0.2480** | [+0.2375, +0.2617] | **0.0 %** |
| **S2** across **DONORS**, state held — 280's statistic, the preregistered **ceiling** | donors | **+0.0440** | [+0.0318, +0.0513] | 25.3 % |
| **S12** across state **and** donors — the floor | both | **+0.0278** | [+0.0204, +0.0328] | 25.7 % |

| contrast | value |
|---|---|
| `S2 − S1` (the preregistered contrast) | **−0.2040 [−0.2214, −0.1918]** |
| `S2 − S12` (value of sharing the STATE) | +0.0163 [+0.0090, +0.0239] |
| `S1 − S12` (value of sharing the DONORS) | **+0.2202 [+0.2089, +0.2375]** |
| **`SD = 1 − S1/S2`** (n = 180 with S2 > 0) | **−2.677 [−3.562, −2.067]**; pooled −4.63; SD > 0 for only **2.9 %** of targets |

**Read literally, the preregistered answer is: within one person a perturbation response survives an
activation-state change 5.6× better than it survives a change of person.** Every sensitivity arm agrees
(280's per-state filter S1 +0.2471 / S2 +0.0463 / SD −2.55; all guides +0.2449 / +0.0456 / −2.22; the
depth-matched N3 arm +0.1748 / +0.0291 / −2.23), the null is flat (KA-4), all six donor pairs give S1
between +0.179 and +0.321 (N4), and S1 is positive for **every one** of the 241 targets.

Per state pair: Rest·Stim8hr **+0.2558**, Rest·Stim48hr **+0.2105**, Stim8hr·Stim48hr **+0.2828** — the two
stimulated states are the closest pair and Rest·Stim48hr the furthest, the expected ordering.
Per state, S2: Rest +0.0446, Stim8hr +0.0304, Stim48hr +0.0463.

### 1.2 …and it is an artefact. Two facts kill it

1. **S1 tracks sequencing depth; S2 does not.** Spearman(S1, cells per target) = **+0.753**
   (with total UMI +0.749), Spearman(S2, cells) = **+0.006**, Spearman(S12, cells) = −0.079. A statistic
   whose two sides come from the same donors and the same constructs inherits that (donor, guide)
   population's *abundance*, and abundance is the strongest predictor of S1 on the panel.
2. **Swapping to the target's other guide destroys it.** §2.

Depth alone does not explain it: with **every** cell thinned to one global budget (4.138e5 UMI, all 18
cells of all 217 surviving targets at exactly the same total), S1 = +0.0878 [+0.0827, +0.0925],
S2 = +0.0132 [+0.0096, +0.0170], S12 = +0.0082 — the ratio is 6.6×, larger than before. So the excess is
not a depth artefact; §2 identifies what it is.

---

## 2. The crossed control (POST-HOC, and decisive)

Added because the preregistered primary produced the opposite of what the preregistration anticipated, and
the interpretation turns entirely on **whether the shared component is the gene's response or the
construct's cell population**. GWCD4i carries the **identical two guides in all four donors and all three
states** for **74** panel targets (verified: 74/74 have one guide-id pair), so donor × guide × state is a
fully crossed 2×2×2 at ~66 cells/side.

| cell | donor | guide | state | full depth | **depth-matched** |
|---|---|---|---|---:|---:|
| `A_s` | same | **same** | diff | **+0.1693** [+0.1550, +0.1793] | **+0.0424** [+0.0387, +0.0607] |
| `A_g` | same | diff | same | +0.0133 [−0.0023, +0.0210] | +0.0116 [+0.0084, +0.0167] |
| `A_gs` | same | diff | diff | +0.0017 [−0.0023, +0.0092] | +0.0089 [+0.0063, +0.0119] |
| `B_g` | diff | **same** | same | +0.0245 [+0.0150, +0.0296] | +0.0141 [+0.0097, +0.0183] |
| `B_gg` | diff | diff | same | +0.0166 [+0.0047, +0.0214] | +0.0112 [+0.0067, +0.0133] |
| `B_s` | diff | **same** | diff | +0.0119 [+0.0088, +0.0162] | +0.0103 [+0.0076, +0.0138] |
| `B_all` | diff | diff | diff | +0.0068 [+0.0042, +0.0095] | +0.0083 [+0.0059, +0.0126] |

| contrast | full depth | depth-matched |
|---|---:|---:|
| **guide cost, state already differs** (`A_s − A_gs`) | **+0.1676 [+0.1508, +0.1777]** | **+0.0334 [+0.0297, +0.0506]** |
| **state cost, guide already differs** (`A_g − A_gs`) | +0.0116 [−0.0017, +0.0176] | +0.0026 [−0.0001, +0.0064] |
| donor cost, guide held, state differs (`A_s − B_s`) | +0.1574 [+0.1427, +0.1679] | +0.0321 [+0.0284, +0.0497] |
| donor cost, guide held, same state (`B_g − B_gg`) | +0.0079 [+0.0005, +0.0203] | +0.0029 [+0.0004, +0.0065] |

**Only one cell is large, and it is the cell that holds both the person and the construct.** Reading the
matched column: two *different guides for the same gene in the same person and the same state* agree at
**+0.0116** — no better than the *same guide in two different people* (+0.0141) and barely above the
all-different floor (+0.0083). But the *same guide in the same person across two activation states* reads
**+0.0424**, 3.7× higher. A gene-level response would be shared by the two guides; a guide's off-target
biology would be shared across donors (`B_g − B_gg` = +0.003). Neither is. What is shared is the
**(donor, guide) cell population** — and it repeats across the three culture arms.

At the donor-pair scale (~107 cells/side, same 74 targets) the same swap gives the size of the artefact
directly: **S1 (same donors, same guide, diff state) = +0.2260 → S1x (same donors, guide swapped, diff
state) = +0.0050 [−0.0069, +0.0133]**; `S1 − S1x` = **+0.2210 [+0.1996, +0.2496]**. The depth signature
goes with it: Spearman(S1, cells) +0.788 → Spearman(S1x, cells) **−0.064**.

**Most parsimonious mechanism** (interpretation, not a measurement): the readme calls the three states a
"culture condition **applied to the cells**", i.e. one transduced, expanded pool per donor split into three
arms, so a (donor, guide) population's clonal and subset composition is shared across the arms while the
cells, the library (12 libraries, one per donor × condition — checked in `sample_metadata.suppl_table.csv`)
and the sequencing are not. Whatever the mechanism, the *measured* fact is the donor × guide interaction,
and it is not the gene's response.

---

## 3. The deconfounded answer to the mandate's question

Both sides now differ in donor **and** guide, so the only factor separating the two rows of the table is the
**state**. 74 targets, ~107 cells/side.

| | same state | diff state | **state cost** | **SD = 1 − diff/same** |
|---|---:|---:|---:|---:|
| **full depth** | **+0.0206** | **+0.0127** | +0.0080 [−0.0026, +0.0186] | **+0.4948 [+0.3292, +0.6080]** (pooled +0.3854) |
| **depth-matched** | +0.0106 | +0.0065 | +0.0040 [−0.0010, +0.0070]\* | **+0.4284 [+0.2757, +0.5348]** (pooled +0.3831) |

\* the same contrast built cell-by-cell from the crossed 2×2×2 (`S2x − S1x` in `pair_guide282.csv`)
reads **+0.0065 [+0.0004, +0.0119]** and excludes 0; the pooled form above touches it. Both are the same
quantity at different aggregations, and at this n the state cost sits at the edge of resolution.

The variant that keeps the donors fixed on the across-state side (closest to the mandate's literal
`S1/S2`, with the construct still broken) reads **SD = 1 − S1x/S2x = +0.677 [+0.484, +1.038]** (matched
+0.548 [+0.351, +0.916]).

**The answer.** *Changing the recipient's activation state costs roughly half of a perturbation
response's already-small reproducibility*: SD **+0.43 … +0.68** depending on the aggregation, i.e. across
two states a response retains **~32–57 %** of the agreement two independent groups of donors reach within
one state. The absolute numbers are tiny (**+0.021 same state, +0.013 across states** at ~107 cells/side)
because the GWCD4i panel response is tiny — 280's F10 regime, not H1's.

**Coherence with 280 (F10).** The donor term is again marginal: `S1x − S12x` = −0.0077 [−0.0186, +0.0016]
(sharing the two donors buys nothing once the construct is broken), matching 280's `Δ_donor` = +0.0033
[+0.0004, +0.0072]. And the value of sharing the state, `S2x − S12x` = +0.0080 [−0.0026, +0.0186].

---

## 4. The shared vs state-specific decomposition (mandate item 3)

Cross-replicate inner products only, so every quantity is unbiased for the reproducible part.
**Replicates that share neither donor nor guide** (pair `P` guide 0 vs pair `P̄` guide 1, both orders, all
three pairings), 74 targets, energy-weighted:

| component | deconfounded (74 targets) | construct-confounded (241 targets, for contrast) |
|---|---:|---:|
| **shared across the three states** | **3.39 %** | 5.23 % |
| **state-specific & reproducible** | **1.45 %** | 1.74 % |
| noise | 95.16 % | 93.03 % |
| reproducible total | 4.84 % | 6.97 % |
| **state share of the reproducible part** | **29.9 %** | 25.0 % |

Per-target medians on the 241-target arm: shared +0.0309 [+0.0264, +0.0378], state-specific +0.0090
[+0.0051, +0.0136], noise +0.9595 [+0.9503, +0.9702].

**Scope note, preregistered:** the residual centring removes each state's across-target mean, so a *global*
activation programme common to every target is already gone. The 1.45 % is the **target-specific** part of
state dependence and is a lower bound on total state dependence.

**Effective rank** (participation ratio of the cross-replicate Gram in the (target, state) space; its
non-zero spectrum equals the gene-space one):

| part | rows | PR | top-1 | top-5 | top-20 |
|---|---:|---:|---:|---:|---:|
| deconfounded state-specific | 222 | **14.0** | 0.214 | 0.466 | 0.776 |
| deconfounded shared | 74 | **14.4** | 0.144 | 0.507 | 0.908 |
| confounded state-specific (3 pairings) | 723 | 20.4 / 20.6 / 24.3 | 0.14–0.18 | 0.37–0.39 | 0.62–0.63 |
| confounded shared | 241 | 15.8 / 17.2 / 18.3 | 0.17–0.20 | 0.43–0.45 | 0.71–0.72 |

**Both the shared and the state-specific reproducible parts are low-dimensional and of the same
dimensionality** — ≈ 14 effective directions over 74 targets, with the top 20 carrying 78–91 % of the
reproducible energy. The state-specific part is *not* one activation axis (top-1 is 0.21, and the global
axis is removed by the centring anyway); it is a handful of programmes, not a per-target private response.

By conservation tercile (241-target arm) the state share of the reproducible part is flat —
**0.2525 / 0.2465 / 0.2474** — while the absolute reproducible energy triples (`shared/T` 0.035 → 0.075).
More conserved targets have *bigger* reproducible responses, split the same way between shared and
state-specific.

---

## 5. S4 — does cross-line conservation predict within-individual state dependence? **No**

Predictor `C_noCD4`: 276's frozen statistic recomputed with all 69 CD4-derived source pairs removed
(KA-3), so nothing in it comes from the cells measured here; `C_raw` (the F8 object) as the declared
secondary.

| test | ρ | 95 % | n | verdict vs the frozen rule (ρ < 0, interval below 0) |
|---|---:|---|---:|---|
| **ρ(SD, `C_noCD4`)**, primary arm | **+0.262** | [+0.120, +0.399] | 178 | **FAIL — interval entirely on the wrong side of 0** |
| ρ(SD, `C_noCD4`), 280-filter arm | +0.275 | [+0.131, +0.412] | 179 | fails identically |
| ρ(SD, `C_noCD4`), all guides / depth-matched | +0.270 / +0.284 | [+0.135,+0.400] / [+0.140,+0.421] | 204 / 179 | fails identically |
| ρ(SD, `C_raw`) (secondary) | +0.267 … +0.289 | all above 0 | 180–206 | fails identically |
| **ρ(SD, `C_noCD4`) on the DECONFOUNDED index** (post-hoc) | **−0.144** | [−0.403, +0.126] | 49 | sign as predicted, **interval covers 0 — NOT SUPPORTED** |

**Why the preregistered test fires the wrong way, and what is actually true.** `SD`'s denominator is the
across-donor ceiling, and *that* is what conservation predicts:

| | ρ with `C_noCD4` | 95 % |
|---|---:|---|
| **S2** (across donors, state held) | **+0.301** | [+0.177, +0.419] |
| **S1** (across states, donors held) | +0.084 | [−0.045, +0.209] |
| S2 − S1 | +0.205 | [+0.079, +0.330] |

and by tercile of `C_noCD4` (low → high): **S1 +0.2575 / +0.2410 / +0.2501 — flat**, while
**S2 +0.0165 / +0.0460 / +0.0673 — a 4× rise**. So the ratio moves entirely with its denominator.

**The scientific reading.** A target's cross-line conservation predicts **how reproducibly its response can
be measured across people** (280's P3 / F8, re-derived here at ρ +0.30 on an independent statistic) and
does **not** predict how much of that response survives a change of activation state within one person.
On the deconfounded rows both sides scale together — ρ(same-state, `C_noCD4`) **+0.338 [+0.101, +0.550]**,
ρ(diff-state, `C_noCD4`) **+0.333 [+0.091, +0.544]** — i.e. conservation buys *magnitude*, equally in both
cells, and cancels in the ratio. **F8/F9's cross-line conservation axis and within-individual state
dependence are different properties.** That is a clean negative for the frozen prediction, on 178 targets
with matched intervals in four arms.

### 5.1 Families (preregistered: HGNC groups with ≥ 8 analysis targets)

| family | n | SD | SD (rest) | Δ [95 %] | S1 | S2 |
|---|---:|---:|---:|---|---:|---:|
| Zinc fingers C2H2-type | 41 | −3.98 | −2.57 | −1.42 [−3.36, +0.68] | +0.264 | +0.0374 |
| MicroRNA protein coding host genes | 16 | −3.61 | −2.63 | −0.98 [−9.46, +1.38] | +0.242 | +0.0313 |
| KRAB domain containing | 16 | −4.85 | −2.61 | −2.23 [−3.97, +0.63] | +0.266 | +0.0296 |
| RNA binding motif containing | 9 | −1.48 | −2.72 | +1.23 [−2.42, +2.22] | +0.248 | +0.0698 |
| Ring finger proteins | 9 | −2.29 | −2.75 | +0.47 [−1.12, +1.95] | +0.196 | +0.0492 |
| **PHD finger proteins** | 8 | −0.10 | −2.76 | **+2.66 [+1.27, +3.61]** | +0.255 | +0.1184 |

Only PHD fingers clear, and inspection of their `S2` (+0.118, 2.7× the panel) shows the same denominator
mechanism as §5 — nothing here is a state-dependence statement. On the deconfounded rows (74 targets,
15 C2H2 zinc fingers) the ZNF contrast is null in every column: state cost Δ −0.0157 [−0.0410, +0.0081],
SD Δ −0.108 [−0.951, +0.170]. **279's family axis does not reappear as within-individual state dependence
on this panel** — which is consistent with §5's main negative, and is a genuinely low-powered test
(6 families, 8–41 members, on a 300-gene panel).

---

## 6. Descriptive tail (POST-HOC, no pass rule, and noisy)

Per-target state cost = (same-state − diff-state) deconfounded cosine, 74 targets. Its split-half
reproducibility over the six replicate constructions is **Spearman +0.458 [+0.230, +0.643]**, so the
ranking carries real but modest signal.

**Most state-dependent:** `SUPT7L +0.218, DNTTIP1 +0.186, MED15 +0.151, CAND1 +0.099, BRPF1 +0.098,
ACLY +0.093, LPCAT3 +0.068, PRDX3 +0.055, FHL3 +0.055, PHF23 +0.052, NCLN +0.052, SRSF5 +0.049,
DLAT +0.044, RAB11FIP5 +0.039, TRPC4AP +0.038, RNF10 +0.036, HDAC8 +0.036, ANP32B +0.034, SIN3B +0.032,
ZNF444 +0.031`.

**Least:** `RBM4B −0.055, DTNBP1 −0.041, ZNF821 −0.040, MAPK12 −0.038, TMEM104 −0.036, MTM1 −0.036,
ZNF280C −0.030, ZNF627 −0.029, VIM −0.026, ZNF526 −0.024, ZNF496 −0.022, ZNF714 −0.018, TWF2 −0.015,
POR −0.015, MAN2A1 −0.012, ZNF251 −0.010, ZC3HAV1 −0.010, ZBTB48 −0.007, ZDHHC4 −0.006, HMGN3 −0.006`.

**This is a magnitude ranking, not a state-dependence ranking.** The top list is SAGA/ATAC (`SUPT7L`),
Sin3–HDAC (`DNTTIP1`, `SIN3B`, `HDAC8`), Mediator (`MED15`), HBO1 (`BRPF1`), `CAND1`, `ACLY` — i.e.
essentially 280 §5's *most donor-reproducible* targets; the bottom list is C2H2 zinc fingers, i.e. 280's
*most donor-variable* / near-zero-response targets. Spearman(state cost, DE-gene count) = **+0.354**, and
the same-state and diff-state cosines correlate with DE count at +0.466 / +0.458. A target with no
measurable response has no measurable state dependence. The **ratio** SD is the magnitude-free version, and
it correlates with DE count at only −0.123.

**Are the state-dependent responses the immune-programme genes?** *Unanswerable on this panel, and that is
the answer.* Labelling by HGNC gene group plus the preregistered keyword list on the HGNC `name` field
(`interleukin|chemokine|cytokine|interferon|T cell receptor|CD\d|tumor necrosis|JAK|STAT|NFKB|MAP kinase|
cyclin|proliferat`) marks **9 of the 241** analysis targets (3.7 %: `CDK8, DMTF1, IFI27L2, IFNAR2, IFNGR2,
IL10RB, MAPK12, MAPK7, STAT6`) and **1 of the 74** deconfounded targets. The VCC panel is a
core-machinery / chromatin / transcription panel, not an immune panel; 0/20 of the top list and 1/20 of the
bottom list are immune-labelled, which at a 1–4 % background is uninformative in either direction.

---

## 7. What this establishes, what it does not

**Established.**
1. **A construct-and-donor-held replicate is not a replicate.** In GWCD4i, the (donor × guide) interaction
   is the dominant reproducible component of a small pseudobulk: +0.042 matched / +0.169 full depth, against
   +0.008 … +0.014 (matched) and +0.002 … +0.025 (full depth) for every cell that breaks either factor, and the two guides of one gene in one person
   agree no better than two different people. This is measured across **twelve independent libraries**, so
   it is not a shared library, lane or run. **It sets a design rule for the programme: an across-timepoint,
   across-arm or across-run comparison that keeps the same donor and the same construct measures the cell
   population, not the perturbation.** 280's donor-half design and F1's guide-construct finding are the
   right constructions; this quantifies what the wrong one costs (4–45×).
2. **The deconfounded within-individual state dependence of a CD4 perturbation response is SD ≈ +0.43 …
   +0.68** — an activation-state change costs roughly half of the (small) reproducibility that two
   independent donor groups achieve within a state. Absolute cosines +0.021 (same state) → +0.013 (across
   states) at ~107 cells/side.
3. **The reproducible response energy splits 3.4 % shared across states / 1.5 % state-specific /
   95.2 % noise**, i.e. **30 % of the reproducible, target-specific response is state-specific**, and both
   parts are low-dimensional (PR ≈ 14 over 74 targets, top-20 ≥ 78 %). The state-specific part is a small
   set of programmes, not a private per-target response — and this is measured *after* the global
   activation programme is removed by the centring, so it is a lower bound.
4. **Cross-line conservation (F8/F9) does not predict within-individual state dependence.** It predicts the
   across-donor ceiling (ρ +0.30) and response magnitude, not the state ratio (deconfounded ρ −0.14
   [−0.40, +0.13]); across-state agreement is flat across conservation terciles. The organising question's
   two axes — "conserved across cell **lines**" and "depends on the recipient's **state**" — are separable
   properties, and the first is not a proxy for the second.
5. **The preregistered statistic S1/S2 as specified in the mandate is not usable on this design**, and the
   experiment says so with its own control rather than reporting the +0.248 headline the frozen rule
   produced. Reported both ways, hashed prereg first.

**Not established.** Nothing here is a competition claim and nothing was priced. The mechanism behind the
donor × guide interaction is inferred, not measured — separating "clonal/subset composition of the cells
carrying that construct" from "guide-by-donor genotype interaction" needs cell-level data (540 GB for one
condition, out of scope). The deconfounded arm rests on **74 targets** (the panel targets carrying the same
two guides in all 4 donors × 3 states); they are more reproducible across donors than the panel
(median S2 +0.052 vs +0.044) though not larger-responding (median DE-gene count 4 vs 5), so they are not a
random subsample. Every deconfounded cosine is +0.005 … +0.021,
close to the resolution of the instrument at ~107 cells/side, and the state-cost interval touches zero in
one of its two aggregations. Four donors, one cell type, one lab, one chemistry; the states are three
culture arms of one design, not three independently sampled states. The estimator is 280's pseudobulk
`log1p(CPM)` residual, which agrees with the study's DESeq2 rows at +0.35–0.53 (R280 §0.1), so absolute
values are not interchangeable with the production bank's.

**What it costs the programme.** 280 closed the depth axis; 282 closes the *within-individual replicate*
axis at this construction: GWCD4i's three states cannot serve as replicates of one another (the naive use
inflates reproducibility 5.6×), and once they are used correctly they measure a real but small state
effect on a reference that is itself near its noise floor. A stage 2 that wants to model state dependence
needs a **cell type with larger responses**, not more states of this one — the same conclusion 280 reached
on the depth axis, reached here from the state axis. The one directly reusable object is the design rule in
(1), which applies to every future replicate-validated measurement in the programme.

---

## 8. Files

`PREREGISTRATION.md` (sha256 `7d75116b…`, hashed 08:03:50Z) · `work/prereg282.sha256` —
**scripts**: `lib282.py` (rows and statistics; imports 280's `lib280.py` unchanged), `meta282.py` (design
metadata, pre-hash), `ka282.py` (KA-1…KA-5), `run282.py` (the preregistered analysis), `probe282.py`
(P-A crossed 2×2×2, P-B global equal depth, P-C depth diagnostics), `probe2_282.py` (the deconfounded
donor-pair statistics), `probe3_282.py` (pooled state contrast, deconfounded decomposition and S4, the
tail) —
**tables**: **`work/state282.csv`** (241 targets × 48 columns: S1/S2/S12/W/SD, per state pair and per
state, the energy decomposition, the three sensitivity arms, `C_noCD4`/`C_raw`/tercile, HGNC family and
name, DE count and cells, and the deconfounded columns for the 74-target subset),
`main282.csv`, `breakdown282.csv`, `s4_282.csv`, `families282.csv`, `decomp282.csv`, `rank282.csv`,
`pa282.csv`, `pb282.csv`, `pair_guide282.csv`, `pooled_state282.csv`, `s4_deconf282.csv`, `tail282.csv` —
**arrays**: `cnocd4_282.npz` (the CD4-free conservation), `pa282_matched.npz`, `pair_guide282.npz` —
**json/logs**: `ka282.json`, `run282.json`/`.log`, `probe282.json`/`.log`, `probe2_282.log`,
`probe3_282.json`/`.log`, `meta282.json`.
Inputs read and not modified: `experiments/280-gwcd4i-donor-replicates/work/{bank280_X.npy,
bank280_rows.npy, obs280.npz, axis280.npz, de_stats_panel280.csv, q1_halves_keepDE.csv, lib280.py,
meta/sample_metadata.suppl_table.csv}`, `experiments/276-conservation-map/work/{conservation276.csv,
cons276.npz}`, `data/external/xcell/hgnc_protein_coding.txt`, `data/challenge/pert_counts.csv`.
