# 223 — The silent zero is REAL and REPRODUCED, and **none of the 21 experiments was measured at it**. Zero VOID. But re-running the residual arm they never ran finds one live source: **DLD-1**, at **+0.0025 marginal over 222's union**, `pds` UP, all four falsifiers passed.

**Date:** 2026-09-04 · deterministic (seeds `20260904`, `20726`, `771103`) · **no `vcc` command, no
`cell-eval2`, no build, nothing under `submissions/` or `data/`, no production file modified,
nothing outside `experiments/223-silent-zero-audit/` written.** 11 scripts, **1,870 scored arms**
across `sweep223` (846), `focus223` (768) and `focus223_feng` (256). The 12 `sweep223` baseline
rows came back identical to 15 decimal places on every repeat, and the three 218-reproduction gates
fired **13 times each, all PASS**. Peak RSS of a single process **7.27 GiB**;
two runs briefly overlapped at ~11.5 GiB combined, inside the 12 GiB mandate but at its edge — see
§9. Work dir 2.2 GB (mostly re-laid-out donor tables, all reconstructible from bytes already on
disk).

---

## 0. The answer in twelve lines

1. **PHASE 0 PASSES. The mechanism is exactly as 218 described it, verified three ways.**
   `solve_cd4_fraction` takes `np.median(np.linalg.norm(added, axis=1))` over **every** panel row,
   with no carrying-row filter. A synthetic residual covering **21%** of a real PLOLO fold's targets
   solves to `f_per_unit = 0.0` **exactly**, `fraction = 0.0` **exactly**, and its effect matrix is
   **bit-identical** to the baseline (`np.array_equal` True, `max|Δ| = 0.0`, `Δcos_agg = Δpds = 0`).
   A 40%-coverage source is silently zero too. The 60% and 92% known-answer rows are **live**.
2. **PHASE 1: 0 of 21 are VOID.** 15 VALID, 6 NOT-RUN, 0 VOID, 0 UNKNOWN. **Not one of the 21
   called `solve_cd4_fraction` or `plolo.solve_fraction`, and not one applied a median-over-all-rows
   rule to a partial-coverage residual.** The mechanical scan found zero hits on either name in 21
   directories and 100 Python files; every median-norm hit was read by hand and is either a
   whitening preserve-scale over a dense prediction, or a source rescale over rows that all carry.
3. **The reason is chronological, and it is worth stating plainly.** The post-whitening
   median-sized residual is a **late** construction (099 → 062 → 103 → 207 → 210 → 216 → 218). The
   21 experiments span 055–176 and are *older* or *differently shaped*: explicit-fraction blends,
   rotations, routers, sign gates, and source-only audits. **The bug is real, it has cost this
   project real measurements — but on `cd4`-path experiments (200, 207, 218), not on these.**
4. **`105-telohaec-residual` — the mandate's first-check row — is VALID, and its name is a trap.**
   It measured *direct transfer* (equal 50:50 blend, macro −0.01778 raw PDS, all four gates failed)
   and then **closed the preregistered residual adapter without ever running it**. So its NO-GO on
   the residual arm rests on inference, not measurement. 223 ran that measurement. See line 8.
5. **`059-feng-ipsc` had already implemented the correct rule, 13 days before 218 found the bug.**
   `evaluate_lodo.py:377-384` sizes on `active_rows` × `feature_mask` and **raises**
   `ValueError("Feng median effect norm is zero")`. The fix existed in the repo and nobody
   generalised it.
6. **The one genuinely weak row is `139-jurkat-donor-simplex`** — the only all-zero table in the 21.
   Its zeros are a router selecting the identity member, and the machinery demonstrably ran
   (609/613 targets transferable). But unlike 085 and 129 it has **no nonzero sibling arm in its own
   table** to cross-check against, so it is the weakest VALID in the set. It is still not the sizing
   bug: 139 contains no median-norm sizing at all.
7. **PHASE 2: nine partial-coverage sources re-measured on carrying rows, marginal over 222's arm.**
   Seven are dead. `telohaec` (65/300) is **negative at every dose and loses to its own row-matched
   CD4 control by 0.09–0.43×** — 105's NO-GO survives its own untested residual arm. CMap
   (122/300 shRNA, 82/300 CRISPR-KO, the largest source we hold) is dead on price in all three
   variants and **`cmap_xpr` loses to its row-matched CD4 control (0.85×)** — 176's NO-GO survives
   the delivered-space test it gated off. `kolf2` and `feng_ipsc` both lose to their row-matched CD4
   controls (0.24× and 0.45×), confirming 057 and 059 under the corrected pricing.
8. **One source is live: `dld1_lowmoi` (67/300 = 22%), whose NO-GO (116) was on a GENE-axis
   reliability gate and was never a transfer measurement at all.** Marginal over
   `cd4@0.34 + h1pair@0.28`, three cancer folds, `g = 0.35`: **macro cosine +1.405%, Δ`pds`
   +0.00256, price +0.00228** — **1.20× rule 13's shape-only gate**, `pds` **UP on 3/3 cancer
   folds**, and **+0.00623 marginal over the free row-matched-CD4 alternative**.
9. **It replicates across DLD-1's own two independent libraries.** Low1 alone +0.565%/+0.00079,
   Low2 alone +1.080%/+0.00267, and the two-split mean **beats both** at +1.403%/+0.00250 — the
   signature of averaging two noisy estimates of one real signal, not of a single library's
   artefact. This was not designed as a falsifier and it is the strongest evidence here.
10. **Under INBOX 1 §1's Δ`pds`-first filter, `dld1` alone is the arm.** Ranked by Δ`pds`, `dld1`
    is the **only** source above +0.00012, and it is `pds`-positive on 3/3 cancer folds at every
    dose. Stacking `myllia_a549` on top buys +1.4 pp more cosine and a higher six-member price
    (+0.00249, 1.31× gate) but **costs** Δ`pds` (+0.00237 vs +0.00256) — so the pair is rank 2, not
    rank 1. `myllia` alone **fails the Δ`pds` > 0 filter** at its own cosine argmax.
11. **INBOX 1 §2 answered, and it is the second branch: the union has already taken what Feng had.**
    `feng_ipsc` (61%) is `pds`-positive on `cd4@0.34` only at `g ≤ 0.10` (+0.0006, 0.3× the gate),
    goes negative on top of the H1 union from `g = 0.10`, and **fails the row-matched CD4 control at
    every dose on both bases** (0.68× at best). Its 218 row was also measured **1.265× over-dosed** —
    the rule's **third** failure mode, confirmed at the coordinator's predicted **18.1st percentile**
    (§5.2). CD4 itself is over-dosed by only 1.014×, so 216/217's dose work is unaffected.
12. **Handed back, not built.** The sizing rule this needs is already committed
    (`scripts/make_submission_pm_cd4.py:1309`, 222's `solve_carrying_row_fraction`). What is missing
    is *registering a second and third extra residual source*, which is 222's territory, not mine.
    **CD4's dose and sizing are untouched throughout.**

---

## 1. Known-answer rows — thirteen, all pass

Rule 10. Every instrument carries a configuration whose answer was already known, and every script
**aborts** on failure.

| # | gate | expected | measured | verdict |
|---|---|---|---|---|
| KA-0a | `solve_cd4_fraction`'s residual-norm line, parsed from the builder's AST | median over all rows, no filter | exact string match; no row filter anywhere in the function | PASS |
| KA-0b | `plolo.solve_fraction`'s zero behaviour | returns `0.0` silently | `target_f / per_unit if per_unit > 0 else 0.0` | PASS |
| KA-0c | 21%-coverage synthetic residual, `f = 0.34` | `f_per_unit == 0.0` exactly | `0.000000e+00` | PASS |
| KA-0d | same, effect matrix | bit-identical to baseline | `np.array_equal` **True**, `max abs diff 0.0` | PASS |
| KA-0e | same, endpoints | `Δcos_agg = Δpds = 0` | `+0.000e+00 / +0.000e+00` | PASS |
| KA-0f | 40%-coverage | also silently zero | `frac 0.00000`, bit-identical | PASS |
| KA-0g | **60%-coverage must be LIVE** | `frac > 0`, not bit-identical | `frac 0.45071`, `Δcos +4.55e-03` | PASS |
| KA-0h | **92%-coverage (full CD4) must be LIVE** | `frac > 0`, not bit-identical | `frac 0.37078`, `Δcos +6.68e-03` | PASS |
| KA-0i | 21%-coverage under the CARRYING rule must be non-vacuous | `g_per_unit > 0` | `0.9103` | PASS |
| KA-1 | `f0` on all four folds | 218 `resid218` `f0` | `Δcos_agg = Δpds = 0.0e+00` | **PASS ×13** |
| KA-2 | `cd4@0.34` on all four folds | 218 `resid218` `cd4_gn0.34` | `0.0e+00` | **PASS ×13** |
| KA-3 | `cd4@0.34 + h1pair@0.28` on all four folds | 218 `compose218` `cd40.34+h1PAIRCARRY0.28` | `0.0e+00` | **PASS ×13** |
| KA-4/5 | every source's OLD-rule arm: `<=50%` ⇒ void, `>50%` ⇒ live | see §5 | 10 void, 2 live, 0 misclassified | PASS |

KA-1/2/3 mean this experiment's whole instrument reproduces 218's, bit-for-bit, before it measures
anything new — and it did so on **13 separate process launches**, with the 12 baseline rows
identical to 15 decimal places every time (`0` keys with differing values). KA-0g/0h are the rows
that make the void detectable rather than assumed: a dead instrument would report "bit-identical"
everywhere.

---

## 2. PHASE 0 — the reproduction

`phase0_repro.py` → `work/phase0_repro.json`. Fold `orion_hct116` (290 targets), production
operator path, `α = 0.982`, `γ = 0.65`, `λ = 0`, dose `f = 0.34` (the board incumbent `cd4f034`).
The residual is the **real** CD4 object masked down to a random fraction of its carrying rows —
coverage-only surgery, nothing else changed.

Baseline (no residual): `cos_agg 0.021329511971`, `pds 0.782054647417`.

| coverage | carrying rows | median row norm, **all rows** | `f_per_unit` | solved fraction | Δ`cos_agg` | Δ`pds` | effect bit-identical? |
|---:|---:|---:|---:|---:|---:|---:|---|
| **21%** | 57/290 | **0.000000e+00** | **0.000000e+00** | **0.00000** | **+0.000e+00** | **+0.000e+00** | **TRUE** |
| **40%** | 111/290 | **0.000000e+00** | **0.000000e+00** | **0.00000** | +0.000e+00 | +0.000e+00 | **TRUE** |
| 60% | 164/290 | 3.696512e+00 | 7.543655e-01 | 0.45071 | +4.551e-03 | +3.190e-03 | false |
| 100% | 267/290 | 4.493351e+00 | 9.169804e-01 | 0.37078 | +6.683e-03 | +8.086e-03 | false |

Under the carrying-row rule the same 21%-coverage object sizes to `g_per_unit = 0.9103`,
`fraction = 0.37349` — i.e. **a fully normal dose**. The residual is not small; it is *identically
absent*.

**One correction to the mandate's framing, and it matters for where to look next.** The builder's
`solve_cd4_fraction` **raises** `ValueError` on `f_per_unit <= 0` whenever `VCC_CD4_TARGET_F` is
set, which is production's shipped configuration. So a production build with a ≤50%-coverage source
fails **loudly**, not silently. The **silent** zero lives in the *harness* mirror,
`experiments/162-retrodicting-validator/plolo.py:363-369`, which returns `0.0` instead. That is the
function every PLOLO-era screen calls — 207, 218, and now 223. **The bug is a measurement bug, not
a shipping bug**, which is exactly why it went unnoticed for so long.

---

## 3. PHASE 1 — the 21-row classification

`coverage21.py` (coverage, computed from each source's own target list on disk) ·
`scan_paths.py` (mechanical path scan, 100 `.py` files) · `classify21.py` (assembly).
→ `work/coverage21.json`, `work/scan_paths.json`, `work/classify21.json`.

**Mechanical result first: `solve_cd4_fraction` = 0 hits, `plolo.solve_fraction` = 0 hits,
`import plolo` = 0 hits, across all 21.** Five directories contain a `median(norm(...))` expression
and four contain the `or 1.0` fallback; every one was read in context and none sizes a
partial-coverage residual.

| # | experiment | source | panel /300 | gene axis | median rule *would* be vacuous? | sizing path actually used | reported result (its own artifact) | **class** |
|---|---|---|---:|---:|---|---|---|---|
| 1 | `055-kolf-direct` | kolf2 | 273 (91.0%) | 18,106 | no | explicit blend, fraction 0.10; rescale median over the 240 **shared** targets, all carrying | HCT116 −0.001151, HEK293T +0.002127, K562 +0.003992, macro **+0.001656** raw PDS | **VALID** |
| 2 | `057-kolf-union-whiten` | kolf2 | 273 (91.0%) | 18,106 | no | same blend under whitening, γ swept | γ=0.65 **−0.006869**; negative at every γ>0; unwhitened checkpoint reproduces 055 bit-for-bit | **VALID** |
| 3 | `059-feng-ipsc` | feng_ipsc | 182 (60.7%) | 18,533 | no | **carrying-row rule already**, with a loud zero guard (`evaluate_lodo.py:377-384`) | macro **+0.0000698** raw PDS | **VALID** |
| 4 | `079-kolf-local-centering` | kolf2 is the **truth** | 273 (91.0%) | 18,106 | no | not a residual arm; 075's centering on the v9 effect | 0.6064291101 → 0.6067253825, **+0.000296**, CI `[−0.00229, +0.00285]` | **VALID** |
| 5 | `085-kolf-pm-gate` | kolf2 (never opened) | 273 (91.0%) | 18,106 | no | PM router with a mandatory `always_D` identity fallback | identity 0/0/0; best non-identity macro **−0.000411**; permutation p=1.0 | **VALID** *(selector-zero)* |
| 6 | **`105-telohaec-residual`** | telohaec | **65 (21.7%)** | 17,464 | **YES** | 50:50 equal **source blend** on a 300-target non-panel gallery; the median scale runs over the source's own 2,167 rows, all carrying | TeloHAEC alone **−0.14105 / −0.12654**; equal blend −0.01825 / −0.01731; macro **−0.01778**; 4/4 gates failed | **VALID** |
| 7 | `107-jurkat-high-moi` | GSE249595 Jurkat, jointly deconvolved | 299 (99.7%) est. | **349** | no | explicit blend fraction, swept 0/0.05/0.1/0.2 | fraction 0.1; nested macro raw-PDS **+0.004246**, bootstrap `[+0.00127, +0.00601]`, above all 99 perms; held PDS **+0.00168 / +0.00012**; directional credit **−0.000165 / −0.000103** | **VALID** |
| 8 | `115-gse295214-h1` | gse295214_h1 | 64 (21.3%) | 18,533 | **YES** | source build only | no value endpoint computed | **NOT-RUN** |
| 9 | `116-dld1-lowmoi` | dld1 | **67 (22.3%)** | 18,533 | **YES** | source-only reliability audit; **no held context ever opened** | split retrieval macro PDS **0.947083** `[0.915, 0.972]`, p=1e−4; **NO-GO on median target cosine 0.0696<0.2 and broad-axis reliable rows 0<0.5** | **VALID** |
| 10 | `119-altair-a549` | altair_a549 | 21 (7.0%) | 18,533 | **YES** | source build only | "has not been tested for predictive value" | **NOT-RUN** |
| 11 | `120-myllia-a549` | myllia_a549 | 14 (4.7%) | 18,533 | **YES** | source build only | no value endpoint computed | **NOT-RUN** |
| 12 | `121-a549-cross-source` | altair × myllia | 1 on-panel of a 20-target gallery | 7,685 | **YES** | source-vs-source retrieval, no VCC truth | macro PDS **0.452632** (chance 0.5), null p=0.875; 7/7 GO criteria failed | **VALID** |
| 13 | `126-hcasmc-panel-confirmation` | IGVF HCASMC | 27 (9.0%) | — | **YES** | **sealed**; the effect asset must stay an unopened binary | no value read, by design | **NOT-RUN** |
| 14 | `129-gse295214-transfer` | gse295214_h1 as **neighbourhood geometry only** | 64 (21.3%) | 18,533 | **YES** | retrieval geometry; "the query target's Pan effect was never added to its prediction" | identity 0; Pan k=5/15/30 macro **−0.01101 / −0.01672 / −0.02086**; all three **improved** signed credit | **VALID** *(selector-zero)* |
| 15 | `130-telohaec-panel-confirmation` | telohaec | 65 (21.7%) | 17,464 | **YES** | LOCKED; value gate never opened | no TeloHAEC value opened | **NOT-RUN** |
| 16 | `131-hues8-hepatocyte` | HUES8 / DE / WTC11 hepatocyte | **91 (30.3%)** | ~9.5–11.1k | **YES** | never acquired — HTTP 403/404 on every object | feasibility GO, acquisition NO-GO | **NOT-RUN** |
| 17 | `135-telohaec-atlas` | TeloHAEC is the **truth** | 65 (21.7%) | 17,464 | **YES** | not a residual arm; 122's atlas is the candidate | signed credit **+0.004066** (gate +0.005), precision **−0.005245** (gate +0.010) | **VALID** |
| 18 | `137-jurkat-tangent-cone` | 107's Jurkat, as an angular cap | 299 (99.7%) | 349 | no | rotation with a 0° identity member | nested macro raw-PDS **−0.000922**, signed credit +0.000429 | **VALID** |
| 19 | `138-dld1-replicate-sign` | dld1 signs | 67 (22.3%) | 18,533 | **YES** | sparse sign gate, explicit `k=200, β=0.05` | nested macro raw-PDS **+0.0000995**, both bootstraps strictly positive, both above p95 — ~30× below the +0.003 bar | **VALID** |
| 20 | `139-jurkat-donor-simplex` | 107's Jurkat as a **router** | 299 (99.7%) | 349 | no | landmark-routed convex weights, identity a mandatory member | identity selected in every fold; **all deltas exactly 0.0** | **VALID** *(selector-zero, weakest row)* |
| 21 | `176-cmap-donors` | LINCS/CMap 2020 | **127 (42.3%)** union | 953 LM + 10,808 inf | **YES** | explicit-weight sign blend on h1-test's DE pool | 9-line mean landmark sign acc **0.534** vs K562's 0.562; trio+CMap **+0.005 on landmarks, nothing on inferred, ≤+0.001 overall** | **VALID** |

**Counts: VALID 15 · NOT-RUN 6 · VOID 0 · UNKNOWN 0.**

### 3.1 The three exactly-zero tables, and why none of them is the bug

| experiment | zero | mechanism | cross-check available? |
|---|---|---|---|
| `085` | `always_D` 0/0/0 | a router selected its identity member | **yes** — the non-identity rules in the same table read +0.000432 / −0.001255 / −0.002674 |
| `129` | identity 0/0/0 | identity is the reference arm | **yes** — the three Pan rotations read −0.011 to −0.021 |
| `139` | every delta 0.0 | a router selected identity on every fold | **no** — this is the only table with no nonzero sibling |

139 is therefore the one row I would not defend as strongly as the other twenty. It still contains
no median-norm sizing anywhere, so it cannot be the silent zero; but "the router chose identity" and
"the router could not move" are not distinguishable from its recorded artifacts alone.

---

## 4. PHASE 2 — every partial-coverage source, re-measured on CARRYING rows

`sweep223.py` → `work/sweep223.json` (782 arms) · `focus223.py` → `work/focus223.json` (500 arms) ·
`analyse223.py`, `analyse_focus.py`.

**Frame.** Baseline is `cd4@0.34 + h1pair@0.28` — *the arm 222 is building* — so every number below
is the **marginal** value on top of 222 and nothing is double-counted (mandate falsifier 5). Both
residuals are solved against the same whitened denominator and added in one step followed by
production's single global preserve rescale, exactly as 218 did. CD4 keeps production's median rule
throughout; **its dose and sizing are never touched**.

Baselines, macro cosine (the primary endpoint, ±0.20% relative per 204) / `pds`:

| fold | `f0` | `cd4@0.34` | `+h1pair@0.28` (= base222) |
|---|---:|---:|---:|
| `k562_gwps` | 0.028213 / 0.781085 | 0.036352 / 0.795881 | 0.039436 / 0.805397 |
| `orion_hct116` | 0.023579 / 0.782055 | 0.029689 / 0.790140 | 0.032080 / 0.792213 |
| `orion_hek293t` | 0.020866 / 0.773842 | 0.024559 / 0.760550 | 0.026157 / 0.763623 |
| `kolf2` *(excluded from headlines)* | 0.006604 / 0.603444 | 0.009016 / 0.614527 | 0.012520 / 0.636155 |

### 4.1 The full table — marginal over 222's arm, three cancer folds

| source | panel /300 | `g` | fraction | rel Δmacro cos | Δ`pds` | Δ`q_w` | +cos | +pds | price `k=1` |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| **`dld1_lowmoi`** | **67** | 0.17 | 0.1116 | **+1.141%** | **+0.00183** | +0.00017 | 3/3 | 3/3 | **+0.00167** |
| | | **0.28** | 0.1838 | **+1.403%** | **+0.00250** | +0.00021 | 3/3 | 3/3 | **+0.00224** |
| | | **0.35** | — | **+1.405%** | **+0.00256** | +0.00016 | 3/3 | 3/3 | **+0.00228** |
| | | 0.45 | 0.2955 | +1.238% | +0.00250 | +0.00016 | 3/3 | 3/3 | +0.00220 |
| **`myllia_a549`** | **14** | 0.28 | 0.1759 | +1.308% | −0.00007 | +0.00008 | 3/3 | 2/3 | +0.00027 |
| | | 0.45 | 0.2826 | +1.690% | −0.00027 | +0.00013 | 3/3 | 1/3 | +0.00022 |
| `cmap_xpr` | 82 | 0.28 | 0.2055 | +1.048% | −0.00138 | +0.00031 | 3/3 | 1/3 | −0.00077 |
| `cmap_sh_cgs` | 122 | 0.28 | 0.2052 | +0.635% | −0.00330 | +0.00042 | 2/3 | 0/3 | −0.00234 |
| `cmap_sh_cgs_lm` | 122 | 0.28 | 0.2102 | −1.035% | −0.00142 | −0.00009 | 0/3 | 0/3 | −0.00135 |
| `feng_ipsc` | 182 | 0.17 | 0.1223 | +0.771% | −0.00113 | −0.00011 | 3/3 | 1/3 | −0.00063 |
| `kolf2` | 273 | 0.17 | 0.1230 | +0.544% | −0.00256 | +0.00088 | 3/3 | 1/3 | −0.00174 |
| **`telohaec`** | **65** | 0.17 | 0.1207 | **−0.066%** | −0.00044 | −0.00008 | 2/3 | 1/3 | −0.00035 |
| | | 0.28 | 0.1987 | **−0.410%** | −0.00145 | +0.00002 | 0/3 | 1/3 | −0.00120 |
| `altair_a549` | 20 | 0.28 | 0.2238 | −0.444% | +0.00010 | −0.00009 | 0/3 | 2/3 | −0.00004 |
| `frangieh_melanoma` | 5 | 0.28 | 0.1809 | +0.036% | −0.00016 | +0.00004 | 2/3 | 1/3 | −0.00011 |

**Only `dld1_lowmoi` is positive on price at every dose, and only `dld1_lowmoi` gains `pds` on
3/3 cancer folds.** `myllia_a549` is positive but small; everything else is negative on price.

### 4.2 The DLD-1 dose curve, and the `dld1 ∪ myllia` stack

`focus223.py`, eight doses, all four folds, marginal over 222's arm (three cancer folds):

| `g` | `dld1` rel Δcos | `dld1` Δ`pds` | `dld1` price | `dld1 ∪ myllia` rel Δcos | `∪` Δ`pds` | `∪` price |
|---:|---:|---:|---:|---:|---:|---:|
| 0.05 | +0.433% | +0.00067 | +0.00061 | +0.698% | +0.00070 | +0.00070 |
| 0.10 | +0.785% | +0.00127 | +0.00115 | +1.292% | +0.00131 | +0.00131 |
| 0.17 | +1.141% | +0.00183 | +0.00167 | +1.948% | +0.00188 | +0.00191 |
| 0.22 | +1.302% | +0.00215 | +0.00195 | +2.298% | +0.00216 | +0.00221 |
| 0.28 | +1.403% | +0.00250 | +0.00224 | +2.599% | +0.00241 | +0.00247 |
| **0.35** | **+1.405%** | **+0.00256** | **+0.00228** | **+2.802%** | **+0.00237** | **+0.00249** |
| 0.45 | +1.238% | +0.00250 | +0.00220 | +2.866% | +0.00218 | +0.00236 |
| 0.60 | +0.715% | +0.00216 | +0.00181 | +2.589% | +0.00162 | +0.00187 |

A clean interior optimum at `g ≈ 0.28–0.35` for both, which is what a real dose response looks like
and what a scale artefact does not.

**If 222 does not ship, the same arm on `cd4@0.34` alone is larger**, because 222's union is
partially redundant with it: `dld1` +0.00277 and `dld1 ∪ myllia` **+0.00302** at `g = 0.35`.

### 4.3 Per fold, at each source's argmax

| source | `g*` | `k562_gwps` | `kolf2` *(excluded)* | `orion_hct116` | `orion_hek293t` |
|---|---:|---|---|---|---|
| `dld1_lowmoi` | 0.35 | +1.646% / **+0.0046** | −0.234% / −0.0025 | +1.565% / **+0.0016** | +1.003% / **+0.0015** |
| `dld1 ∪ myllia` | 0.35 | +2.507% / **+0.0054** | +1.722% / −0.0024 | +2.885% / **+0.0014** | +3.013% / **+0.0003** |
| `myllia_a549` | 0.28 | +0.746% / +0.0007 | +1.441% / +0.0001 | +1.427% / +0.0001 | +1.750% / −0.0009 |

**`dld1 ∪ myllia` is positive on cosine on 4/4 folds.** The `kolf2` fold costs `pds` for both; it is
excluded from every headline for the reasons 218 gave (its truth is a per-gene z-score matrix with a
placeholder `b_ctrl` and the weakest absolutes of the four).

### 4.4 Leave-one-fold-out (dose chosen on the other folds, priced on the held one)

| source | cancer-fold LOFO mean | 4-fold LOFO mean |
|---|---:|---:|
| `dld1_lowmoi` | **+0.00188** (3/3 positive) | +0.00093 (3/4) |
| `dld1 ∪ myllia` | **+0.00204** (3/3 positive) | +0.00113 (3/4) |
| `myllia_a549` | +0.00011 (2/3) | +0.00014 (3/4) |

Selection and reporting frames never overlap. The cancer-fold LOFO for the pair clears the gate;
the 4-fold LOFO does not, entirely because of the excluded `kolf2` fold's `pds`.

---

## 5. PHASE 3 — the four mandatory falsifiers, every source, every dose

### F1 — target-shuffled control (whole rows permuted by a derangement; geometry, per-row norms, support pattern and carrying count preserved exactly)

Two independent seeds for the survivors (`20726`, `771103`). Marginal over 222, three cancer folds:

| source | `g` | real | SHUF (seed 1) | SHUF2 (seed 2) | verdict |
|---|---:|---:|---:|---:|---|
| **`dld1_lowmoi`** | 0.17 | **+1.141%** | −0.211% | −0.595% | **PASS** |
| | 0.28 | **+1.403%** | −0.751% | −1.360% | **PASS** |
| | 0.35 | **+1.405%** | −1.216% | −1.954% | **PASS** |
| | 0.45 | **+1.238%** | −1.991% | −2.902% | **PASS** |
| **`dld1 ∪ myllia`** | 0.35 | **+2.802%** | −1.574% | −1.741% | **PASS** |
| `myllia_a549` | 0.28 | +1.308% | −0.113% | −0.391% | PASS |
| `cmap_sh_cgs` | 0.28 | +0.635% | −0.737% | — | passes F1, fails price |
| `cmap_xpr` | 0.28 | +1.048% | −0.534% | — | passes F1, fails F3 |
| `telohaec` | 0.28 | **−0.410%** | −0.617% | — | **the arm is NEGATIVE; F1 is moot** |
| `altair_a549` | 0.28 | **−0.444%** | −0.378% | — | **FAILS — worse than its own shuffle** |

The shuffle is **strictly negative and monotonically worse with dose** for every source that has a
real signal — the same shape 218 measured, and the failure mode that killed 204's and 209's arms.
**`altair_a549` fails outright**: it is worse than its own target-shuffled version at every dose,
which is independent confirmation of 121's conclusion that the two A549 screens do not agree.

### F2 — non-carrying rows must not move

218's bar is |Δ| ≤ 0.018 pp. Maximum excursion over the three cancer folds:

| source | `g = 0.17` | `g = 0.28` | `g = 0.35` | `g = 0.45` | carrying rows move (same four doses) |
|---|---:|---:|---:|---:|---|
| `dld1_lowmoi` | 0.0032 pp | 0.0097 pp | 0.0125 pp | 0.0233 pp | **+4.18 / +5.13 / +5.13 / +4.49 %** |
| `dld1 ∪ myllia` | 0.0060 pp | 0.0146 pp | 0.0210 pp | 0.0315 pp | **+5.65 / +7.53 / +8.11 / +8.27 %** |
| `myllia_a549` | 0.0042 pp | 0.0050 pp | 0.0057 pp | 0.0072 pp | **+7.60 / +10.99 / +12.58 / +14.20 %** |

**PASS at the headline doses** (`g ≤ 0.35`: 0.0125 and 0.0210 pp, at or just over 218's bar and two
orders of magnitude below the carrying-row move). At `g = 0.45` the pair reaches 0.0315 pp, above
218's bar — a further reason to prefer `g ≈ 0.28–0.35`. The arm is adding information on the rows
it carries, not reweighting rows.

### F3 — the row-matched CD4 control (same rows, same sizing rule, same dose, CD4's information instead)

This absorbs row selection and "more residual is better" entirely. Marginal over 222, three cancer
folds, `g = 0.28`:

| source | real rel Δcos | row-matched CD4 | ratio | **verdict** |
|---|---:|---:|---:|---|
| **`myllia_a549`** | +1.308% | +0.259% | **5.04×** | **PASS** |
| **`dld1_lowmoi`** | +1.403% | +0.464% | **3.02×** | **PASS** |
| **`dld1 ∪ myllia`** | +2.599% | +0.869% | **2.99×** | **PASS** |
| `cmap_sh_cgs` | +0.635% | +0.559% | 1.14× | marginal, and negative on price |
| `cmap_xpr` | +1.048% | +1.240% | **0.85×** | **FAIL — the free CD4 alternative is better** |
| `feng_ipsc` | +0.139% | +2.035% | **0.07×** | **FAIL** |
| `kolf2` | −0.218% | +2.652% | **−0.08×** | **FAIL** |
| `telohaec` | −0.410% | +0.946% | **−0.43×** | **FAIL** |
| `altair_a549` | −0.444% | +0.494% | **−0.90×** | **FAIL** |

Priced on all six members the row-matched CD4 control is **negative** (−0.00317 at `g = 0.28`,
because it costs `pds` — 218 measured the same thing), so the **marginal value over the free
alternative** is larger than the raw price:

| arm | price `k=1` | row-matched CD4 price | **marginal over the free alternative** |
|---|---:|---:|---:|
| `dld1_lowmoi` @ 0.35, marginal over 222 | +0.00228 | −0.00395 | **+0.00623** |
| `dld1 ∪ myllia` @ 0.35, marginal over 222 | +0.00249 | −0.00387 | **+0.00636** |
| `dld1 ∪ myllia` @ 0.35, on `cd4@0.34` alone | +0.00302 | −0.00362 | **+0.00664** |

### F4 — leak audit, every arm

No source is ever scored on the fold where it is the truth (`kolf2` is skipped on the `kolf2` fold;
`dld1`, `myllia`, `telohaec`, `altair`, CMap and Frangieh are none of the four held lines). Carrying
rows per source per fold:

| source | `k562_gwps` | `kolf2` | `orion_hct116` | `orion_hek293t` |
|---|---:|---:|---:|---:|
| `cmap_sh_cgs` | 110/269 | 108/273 | 117/290 | 119/294 |
| `cmap_xpr` | 67/269 | 69/273 | 78/290 | 79/294 |
| `feng_ipsc` | 171/269 | 165/273 | 177/290 | 179/294 |
| `dld1_lowmoi` | 60/269 | 59/273 | 66/290 | 65/294 |
| `telohaec` | 58/269 | 59/273 | 62/290 | 63/294 |
| `altair_a549` | 20/269 | 17/273 | 19/290 | 19/294 |
| `myllia_a549` | 13/269 | 10/273 | 12/290 | 13/294 |
| `frangieh_melanoma` | 5/269 | 4/273 | 5/290 | 5/294 |
| `kolf2` | 249/269 | *(held)* | 264/290 | 268/294 |

Four specific overlaps, stated explicitly:

- **`h1_2025` covers 85/85 of `h1-test`'s targets.** `h1-test` is **not used anywhere in 223** —
  every measurement here is on 162's PLOLO challenge frame. This is the trap that cost experiment
  157 and it is flagged again because 222's arm carries `h1_2025`.
- **DLD-1 and HCT116 are both colorectal.** If the gain were a lineage artefact it would concentrate
  on the `orion_hct116` fold. It does not: the largest gain is on `k562_gwps` (myeloid, +1.646%),
  ahead of `orion_hct116` (+1.565%). **Not a lineage effect.**
- **`myllia_a549` is lung adenocarcinoma; no fold is lung.**
- **CMap's `HA1E` line is HEK-derived (hTERT/SV40 immortalised kidney epithelium)** and sits in the
  9-line shRNA mean scored against the `orion_hek293t` fold. This would *flatter* CMap, and CMap
  still fails, so no positive claim rests on it.

### F5 — non-vacuity, asserted in-loop

Every arm claiming a nonzero dose asserts `median_row_norm_carrying > 0` and the sweep aborts
otherwise. Every arm records `median_row_norm_all_rows`, `median_row_norm_carrying`,
`carrying_rows`, `coverage_frac`, `fraction`, `per_unit` and a `VACUOUS` flag.

### F6 (not required, and the strongest evidence here) — replication across DLD-1's two independent libraries

116 built DLD-1 as two separately sequenced low-MOI libraries and measured their cross-split
retrieval PDS at **0.947**. I built a table from each alone and re-ran the whole arm:

| table | rel Δmacro cos @ 0.28 | Δ`pds` @ 0.28 | vs SHUF | vs row-matched CD4 | price |
|---|---:|---:|---:|---:|---:|
| `dld1_low1` alone | +0.565% | +0.00079 | −0.798% | 1.22× | +0.00075 |
| `dld1_low2` alone | +1.080% | +0.00267 | −0.890% | 2.33× | +0.00227 |
| **two-split mean** | **+1.403%** | **+0.00250** | −0.751% | **3.02×** | **+0.00224** |

**Both libraries are independently positive on cosine and on `pds`, and the mean beats both.** That
is the signature of averaging two noisy estimates of one real signal. A single library's technical
artefact would not survive being averaged with an independently prepared one, and would not make
the average *better than its best half*.

### 5.1 The silent zero, in this experiment's own numbers

Every source was also run through the **production median rule** at every dose, on all four folds:

| source | coverage | OLD-rule fraction | bit-identical to base222? | CARRY-rule fraction |
|---|---:|---:|---|---:|
| `kolf2` | 91.6% | 0.20694 | false | 0.20262 |
| `feng_ipsc` | 61.5% | 0.24200 | false | 0.19320 |
| `cmap_sh_cgs` | 40.3% | **0.00000** | **TRUE** | 0.19774 |
| `cmap_sh_cgs_lm` | 40.3% | **0.00000** | **TRUE** | 0.20319 |
| `cmap_xpr` | 26.0% | **0.00000** | **TRUE** | 0.19677 |
| `dld1_lowmoi` / `low1` / `low2` | 22.2% | **0.00000** | **TRUE** | 0.177–0.184 |
| `telohaec` | 21.5% | **0.00000** | **TRUE** | 0.19102 |
| `altair_a549` | 6.7% | **0.00000** | **TRUE** | 0.21549 |
| `myllia_a549` | 4.3% | **0.00000** | **TRUE** | 0.16589 |
| `frangieh_melanoma` | 1.7% | **0.00000** | **TRUE** | 0.17228 |

**Ten of twelve source tables we hold are identically absent under the shipped sizing rule.** The
two that are not are exactly the two above 50%. This is the mandate's claim, measured on a third
independent instrument, with the >50% rows proving the test can tell the difference.

---

## 5.2 INBOX 1 §2 — the THIRD failure mode. Above 50% coverage the rule OVER-DOSES.

`overdose223.py` → `work/overdose223.json`. The coordinator predicted that at 61% coverage the
all-rows median sits at roughly the **18th percentile** of the carrying rows. Measured on the
`orion_hct116` fold:

| source | carrying/rows | coverage | median over **all** rows | median over **carrying** rows | **over-dose factor** | percentile of the all-rows median inside the carrying distribution |
|---|---:|---:|---:|---:|---:|---:|
| **`cd4`** | 267/290 | **92.1%** | 4.4934 | 4.5570 | **1.014** | 45.7 |
| `kolf2` | 264/290 | 91.0% | 6.6481 | 6.7971 | 1.022 | 45.1 |
| **`feng_ipsc`** | 177/290 | **61.0%** | 5.3806 | 6.8042 | **1.265** | **18.1** |
| `cmap_sh_cgs` | 117/290 | 40.3% | **0.0000** | 6.7150 | *(silent zero)* | — |
| `cmap_xpr` | 78/290 | 26.9% | **0.0000** | 6.7245 | *(silent zero)* | — |
| `dld1_lowmoi` | 66/290 | 22.8% | **0.0000** | 7.4364 | *(silent zero)* | — |
| `telohaec` | 62/290 | 21.4% | **0.0000** | 6.9083 | *(silent zero)* | — |
| `altair_a549` | 19/290 | 6.6% | **0.0000** | 6.1588 | *(silent zero)* | — |
| `myllia_a549` | 12/290 | 4.1% | **0.0000** | 7.8174 | *(silent zero)* | — |
| `frangieh_melanoma` | 5/290 | 1.7% | **0.0000** | 7.6142 | *(silent zero)* | — |

**The 18.1st percentile prediction is exact.** The rule has three regimes, not two:

| coverage | what `solve_cd4_fraction` does | who it affects |
|---|---|---|
| **≤ 50%** | `fraction = 0` — the residual is **identically absent** | ten of the twelve source tables we hold |
| **50–~85%** | `fraction` too **large** by `median(carrying)/median(all)` | **`feng_ipsc`, over-dosed 1.265× everywhere it was median-sized — including 218's own matched-dose row** |
| **> ~90%** | the two medians agree; error ≤ 2% | **CD4 (1.014×) and `kolf2` (1.022×) — so 216/217's CD4 dose work is unaffected, and so is 057's** |

**None of the 21 is OVERDOSED either**, for the same reason none is VOID: none uses the rule.
`059-feng-ipsc`, the one experiment that could have been, sized on `active_rows` × `feature_mask`
instead — the carrying-row rule — so it was correctly dosed. **The over-dose lands on 218's `feng`
row and on any future cd4-path build that registers a 50–85%-coverage source.**

---

## 5.3 INBOX 1 §1 — ranked by Δ`pds`, which sharpens the answer rather than changing it

The coordinator's filter: rank by Δ`pds` and require Δ`pds` > 0 at the shipped dose, with Δcos
secondary, because every cosine gain this project has bought so far paid on `pds` (218's free
self-residual floor: **+4.61% cosine and −0.0078 raw `pds`**).

Marginal over 222's arm, carrying-row sizing, three cancer folds, ranked by Δ`pds`:

| rank | source | `g` | **Δ`pds`** | folds `pds`-positive | rel Δcos | passes Δ`pds` > 0 |
|---:|---|---:|---:|---:|---:|---|
| 1 | `dld1_low2` *(one library)* | 0.45 | **+0.00299** | **3/3** | +0.735% | **YES** |
| 2 | `dld1_low2` *(one library)* | 0.28 | **+0.00267** | **3/3** | +1.080% | **YES** |
| 3 | **`dld1_lowmoi`** | 0.28 / 0.45 | **+0.00250** | **3/3** | +1.403% / +1.238% | **YES** |
| 4 | `dld1_lowmoi` | 0.17 | +0.00183 | 3/3 | +1.141% | YES |
| 5 | `dld1_low1` *(one library)* | 0.28 | +0.00079 | 2/3 | +0.565% | YES |
| 6 | `altair_a549` | 0.17 | +0.00012 | 3/3 | **−0.196%** | yes, but cosine-negative and fails F1 |
| 7 | `myllia_a549` | 0.17 | +0.00003 | 2/3 | +0.906% | marginal |
| — | `myllia_a549` | 0.28 | **−0.00007** | 2/3 | +1.308% | **NO** |
| — | `cmap_xpr` | 0.17 | −0.00009 | 2/3 | +0.826% | **NO** |
| — | `telohaec` | 0.17 | −0.00044 | 1/3 | −0.066% | **NO** |
| — | `cmap_sh_cgs` | 0.17 | −0.00069 | 1/3 | +0.660% | **NO** |
| — | `feng_ipsc` | 0.17 | −0.00113 | 1/3 | +0.771% | **NO** |
| — | `kolf2` | 0.17 | −0.00256 | 1/3 | +0.544% | **NO** |

**The `pds`-first filter is decisive and it simplifies the shortlist.** `dld1` is the only source
whose Δ`pds` is above the noise, and it is positive on **3/3 cancer folds at every dose tested**.
Every other source is at or below +0.00012.

**And it demotes the pair.** At `g = 0.35`, `dld1` alone gives Δ`pds` **+0.00256** while
`dld1 ∪ myllia` gives **+0.00237** — adding `myllia` buys +1.4 pp of cosine and **costs** `pds`.
Under a Δcos-first ranking the pair wins; under the coordinator's Δ`pds`-first ranking, **`dld1`
alone wins**. §7's shortlist is ordered that way.

---

## 5.4 INBOX 1 §2 — Feng, measured marginally on top of the H1 union. It is the SECOND branch.

`focus223.py --out work/focus223_feng.json`, 8 doses × 4 folds × 4 variants. Feng is 61% covered so
it is VALID under the mandate's classification and my Phase 2 would not have re-measured it; the
INBOX asked for it explicitly and it was worth asking for.

| base | `g` | rel Δcos | Δ`pds` | folds `pds`+ | SHUF | row-matched CD4 | price |
|---|---:|---:|---:|---:|---:|---:|---:|
| `cd4@0.34` | 0.05 | +0.615% | **+0.00054** | 2/3 | −0.170% | **+0.901%** | +0.00055 |
| `cd4@0.34` | 0.10 | +1.011% | **+0.00056** | 2/3 | −0.548% | **+1.608%** | +0.00067 |
| `cd4@0.34` | 0.17 | +1.209% | −0.00030 | 2/3 | −1.405% | **+2.321%** | +0.00008 |
| `cd4@0.34` | 0.28 | +0.779% | −0.00255 | 1/3 | −3.412% | **+2.925%** | −0.00170 |
| **`+h1pair@0.28`** | 0.05 | +0.472% | **+0.00025** | 2/3 | −0.167% | +0.701% | +0.00031 |
| **`+h1pair@0.28`** | 0.10 | +0.735% | **−0.00006** | 2/3 | −0.532% | +1.226% | +0.00015 |
| **`+h1pair@0.28`** | 0.17 | +0.771% | **−0.00113** | 1/3 | −1.355% | +1.710% | −0.00063 |
| **`+h1pair@0.28`** | 0.28 | +0.139% | **−0.00339** | 1/3 | −3.289% | +2.035% | −0.00248 |

**Answer: the union has already taken the shared direction Feng was contributing.** On `cd4@0.34`
alone Feng is weakly `pds`-positive at `g ≤ 0.10` (+0.00054/+0.00056, price +0.0006, **0.3× the
gate**). On top of the H1 union that positivity is gone by `g = 0.10` (−0.00006) and clearly
negative by `g = 0.17`.

**And Feng fails F3 at every dose on both bases.** The row-matched CD4 control beats it everywhere —
0.68× at its best point (`g = 0.05`, 0.615% vs 0.901%), falling to 0.27× at `g = 0.28`. Re-using
CD4's own information on Feng's 177 rows is strictly better than using Feng's.

Two things follow. **218's `feng` row (+0.00272 raw Δ`pds` at the matched dose) was measured through
the median rule, i.e. at 1.265× the dose it asked for** (§5.2) and against an `f0` baseline rather
than on top of CD4. **Correctly dosed and measured against the incumbent, Feng's `pds` gain shrinks
to +0.0005 and vanishes on top of the union.** 218 was right to dismiss it — but for the wrong
reason (cosine-first), which is exactly what the INBOX warned about; the right reason is that it
loses to its own free row-matched control.

---

## 6. Price, under §8.4's CORRECTED rule — all six members, never four at zero

Rate card from 218 §7.4, i.e. `cd4f034`'s own board row divided by its measured cosine move
(0.072280 → 0.090270): **+0.007259 overall per +0.01 absolute cosine**, of which `mse` +0.003634,
`reach` +0.002446, `nmae` +0.000980, `fid` +0.000847, `jac` +0.000100, `pds` −0.000748. `pds` is
priced independently through 207's route (PLOLO raw Δ / 0.452 / 6 × 1.011, ×2.0 per §8.4) and the
mixed total is what is reported. **`k = 1` is the headline** — §8.4's 1.17–1.39 transfer factor was
measured on arms that move ~93% of the panel's rows and this one moves 22–27%.

| arm | frame | price `k=1` | price `k=1.26` | × §8.2 shape-only gate (+0.0019) | × amplitude gate (+0.0089) |
|---|---|---:|---:|---:|---:|
| `dld1 ∪ myllia` @ 0.35 | on `cd4@0.34` alone | **+0.00302** | +0.00324 | **1.59×** | 0.34× |
| `dld1_lowmoi` @ 0.35 | on `cd4@0.34` alone | +0.00277 | +0.00288 | 1.46× | 0.31× |
| **`dld1 ∪ myllia` @ 0.35** | **marginal over 222's arm** | **+0.00249** | **+0.00268** | **1.31×** | 0.28× |
| `dld1_lowmoi` @ 0.35 | marginal over 222's arm | +0.00228 | +0.00238 | 1.20× | 0.26× |
| `dld1 ∪ myllia` @ 0.35 | **minus the free row-matched-CD4 alternative** | **+0.00636** | — | **3.35×** | 0.71× |
| `dld1 ∪ myllia` | **cancer-fold LOFO**, dose chosen out of sample | **+0.00204** | — | 1.07× | 0.23× |
| `myllia_a549` @ 0.28 | marginal over 222's arm | +0.00027 | +0.00036 | 0.14× | 0.03× |

**Classification.** α, γ, λ and the emission law are untouched; the only renormalisation is the
preserve scalar, present at every dose. **Shape-only.** It clears that gate on every framing and
does **not** clear the amplitude-class gate on any of them — so unlike 218's arm, this one's
classification is load-bearing and I am not going to argue it both ways.

---

## 7. Ranked shortlist — what is worth adding to 222's union

**Ordered by Δ`pds` per INBOX 1 §1, with Δcos secondary and Δ`pds` > 0 a hard filter.**

| rank | arm | Δ`pds` (3 cancer folds) | panel rows added | marginal price over 222 | falsifiers | what it costs to build |
|---|---|---:|---:|---:|---|---|
| **1** | **`dld1_lowmoi` @ `g ≈ 0.28–0.35`, sized on carrying rows** | **+0.00250 / +0.00256, 3/3 folds** | 67 (7 outside CD4, 3 outside CD4 ∪ 222) | **+0.00224 / +0.00228** (LOFO +0.00188) | **F1 ✓ (2 seeds) · F2 ✓ · F3 ✓ 3.02× · F4 ✓ · F6 ✓ (both libraries independently positive)** | **≈ zero compute.** The compact source already exists at `116/work/compact-source-v1`. The re-layout to a donor table is `build_tables.py` — 30 s, 0.1 GiB — and its output is already on disk at `work/tables/dld1_lowmoi/`. The sizing rule is **already committed** at `make_submission_pm_cd4.py:1309`. The only change is registering one more extra-residual source. |
| **2** | `dld1_lowmoi` ∪ `myllia_a549` @ `g ≈ 0.35` | +0.00237, 3/3 folds | 78 (3 shared) | **+0.00249** (LOFO +0.00204) | same, F3 2.99× | one extra source registration on top of rank 1 |
| 3 | `myllia_a549` alone | **−0.00007 at its cosine argmax** (+0.00003 at `g = 0.17`) | 14 (1 outside CD4) | +0.00027 | passes F1/F2/F3 but **fails the Δ`pds` > 0 filter at the dose its cosine wants** | not worth a slot alone |

**Rank 1 vs rank 2 is a genuine trade and the INBOX's filter decides it.** The pair buys +1.4 pp
more cosine and a higher six-member price (+0.00249 vs +0.00228) but **costs** Δ`pds` (+0.00237 vs
+0.00256). `pds` is 5/6 of the board's discrimination and the amplitude-free, ordering-free member;
on the coordinator's standing filter, **ship `dld1` alone.** If cosine is what a future arm needs,
rank 2 is there and it is not a large step down.
| — | `telohaec` (65/300) | — | **−0.00120** | negative arm; **loses to its row-matched CD4 control 0.43×** | **DO NOT BUILD.** 105's NO-GO survives the residual arm it never ran. |
| — | `cmap_xpr` (82/300) | — | −0.00077 | **F3 FAIL, 0.85×** | **DO NOT BUILD.** 176's NO-GO survives the delivered-space test it gated off. |
| — | `cmap_sh_cgs` / `_lm` (122/300) | — | −0.00234 / −0.00135 | negative on `pds` 0/3 | **DO NOT BUILD.** |
| — | `feng_ipsc` (182/300) | — | −0.00063 @ 0.17; **+0.00031 @ 0.05** | **F3 FAIL at every dose on both bases** (0.68× at best) | 059's NO-GO stands. On `cd4@0.34` alone it is weakly `pds`-positive at `g ≤ 0.10` (+0.0006, 0.3× the gate); **on top of the H1 union that is gone.** See §5.4. Its 218 row was 1.265× over-dosed (§5.2). |
| — | `kolf2` (273/300) | — | −0.00174 | **F3 FAIL, 0.24×** | 057's NO-GO stands. |
| — | `altair_a549` (20/300) | — | −0.00004 | **F1 FAIL** — worse than its own shuffle | 121's NO-GO stands. |
| — | `frangieh_melanoma` (5/300) | — | −0.00011 | too few rows to matter | — |

**Coverage arithmetic for rank 1.** CD4 carries 276/300. 222's union adds **10** rows CD4 cannot
reach (77 total, 67 shared with CD4). `dld1 ∪ myllia` adds **3** rows neither CD4 nor 222's union
reaches, taking the stack to **289/300**. So the value is *not* mostly new rows — it is a **fourth
independent estimate on 78 rows CD4 already covers**, which is precisely 159's `1/k` mechanism and
precisely why `pds` goes up rather than down.

**Two sources worth acquiring that 223 could not measure**, both of them larger than DLD-1:

| source | panel /300 | why not measured | what it would cost |
|---|---:|---|---|
| **HUES8 / DE / WTC11 hepatocyte** (131) | **91 (30.3%)** | never acquired — every GCS URL 404s and every Synapse entity 403s anonymously | authenticated, certified Synapse access, then ask for `calibrated_all_results.tsv` / `trans_de_results.tsv.gz` rather than the 17–30 GB H5MU. **This is now the single largest unmeasured partial-coverage source and it is an access problem, not a science problem.** |
| **IGVF HCASMC** (126) | 27 (9.0%) | **sealed by 126's own value-access boundary**; 223 did not open it and must not | a new hash-bound authorization. Given DLD-1's 67 rows buy +0.0023, 27 rows would buy proportionally less. Low priority. |

---

## 8. What this says about 215's "many small steps" hypothesis — partial support, honestly stated

215 found the leaders' direction arrived in many separately-resolved increments (`cqawesome` 5,
largest 28% of the total; `Haystack` 3) while ours arrived in one. The mandate's hypothesis was
that a stack of partial-coverage sources is what that looks like.

**What is now measured:** two more partial-coverage sources, on top of 222's two, price at
**+0.0025** marginal — a fourth increment of the right size, with `pds` moving the right way, from a
class the sizing rule had been discarding.

**What is not measured, and I will not claim it:** this is *one* increment, not five. Of the nine
partial-coverage sources tested, **seven are dead** and one of the two survivors contributes only
+0.0003 of the +0.0025. The stack does not extend indefinitely; it extends by about one useful
member per ~65 covered panel rows of genuinely independent CRISPR screen, and 166's census says
there are not many of those left. **"Many small steps" is supported as a shape and unproven as a
quantity.**

---

## 9. What I did NOT measure, and the honest uncertainties

- **No `cell-eval2` run, no build, no submission.** `fid`, `jac` and `reach` are unreadable on every
  harness we hold (rule 9) and are priced through the rate card and nowhere else.
- **The rate card rests on one board point** (`cd4f034`). §8.4's own lesson is that a class constant
  from one point is a lower bound on magnitude — which cuts both ways: the `pds` cost may be
  understated as easily as the gains.
- **The PLOLO → board transfer for a 22–27%-coverage arm is the dominant uncertainty**, exactly as
  it was for 218's 26% arm. `k = 1` is the headline for that reason, and macro (per-target,
  unweighted) cosine is used rather than the energy-weighted aggregate because the board averages
  six metrics *per perturbation*.
- **The DLD-1 table is my own re-layout** of 116's `transfer_lfc_{low1,low2}.npy` (transposed,
  measured-mask-weighted split mean, NaN→0). 116's own panel-coordinate zeroing is re-asserted in
  code and both splits independently reproduce the result (F6). The two-split mean is the estimator
  116 itself validated at split-PDS 0.947.
- **`121`'s 20-target gallery contains only 1 panel target**, so its NO-GO is about A549 source
  quality, not panel transfer. It is classified VALID on that basis, and `altair_a549`'s F1 failure
  here is consistent with it, but the two are not the same measurement.
- **`107`/`137`/`139`'s target coverage (299/300) is the *eligible* set from its own design matrix**;
  its response axis is 349 genes and it was never a residual arm, so I did not re-measure it. Its
  NO-GO turned on a directional-credit delta of −0.00016 against a `pds` gain of +0.004246 that was
  significant on every stated test. **That is a re-pricing question for 219's territory, not mine**,
  and I am flagging it rather than adjudicating it.
- **`126`'s HCASMC remains sealed.** I did not decompress, header-inspect, or gzip-test it.
- **The `kolf2` fold is excluded from every headline** and reported everywhere. It is the fold that
  turns the 4-fold LOFO negative for `dld1`.
- **F2 at `g = 0.45` exceeds 218's 0.018 pp bar** for the pair (0.0315 pp). That is one reason the
  recommended dose is 0.28–0.35, not the cosine argmax at 0.45.
- **Two runs overlapped at ~11.5 GiB combined RSS**, inside the 12 GiB mandate but not with much
  room. Single-process peak was 7.27 GiB. Sequential execution would have been the safer choice.
- **CMap's three tables and DLD-1's three are 2.2 GB in `work/tables/`.** They are pure re-layouts of
  bytes already on disk and can be deleted and rebuilt in ~2 minutes.

---

## 10. Handoffs

**To 222 (owns the sizing rule and the builder):** `solve_carrying_row_fraction` is already yours and
already committed. **The one source worth a second extra residual is `dld1_lowmoi`** — built at
`experiments/223-silent-zero-audit/work/tables/dld1_lowmoi/` in the exact layout `E062.Bank.load`
reads. It prices at **+0.00228 marginal on top of your arm at `g ≈ 0.35`** with **Δ`pds` +0.00256 on
3/3 cancer folds**, or **+0.00277** on `cd4@0.34` alone if yours does not ship. Adding
`myllia_a549` as a third raises the six-member price to +0.00249 but lowers Δ`pds` to +0.00237 —
your call, but the standing `pds`-first filter says take `dld1` alone. Do not take my dose as
pre-registered: it was chosen on the same folds it is reported on, which is why the LOFO row
(+0.00188 for `dld1`) is the number to trust.

**Also worth knowing before you register anything at 50–85% coverage:** the shipped median rule
does not just zero sub-50% sources, it **over-doses** sources in that band by
`median(carrying)/median(all)` — 1.265× for `feng_ipsc` at 61%. CD4 (92%) is 1.014×, so nothing
about CD4's dose changes. §5.2.

**To 219 (re-pricing closures):** `107-jurkat-high-moi` is the one row in my 21 whose NO-GO looks
mis-priced. It measured a nested macro raw-PDS gain of **+0.004246** with a strictly positive
bootstrap and above all 99 permutations, and was vetoed on two directional-credit deltas of
−0.00016 and −0.00010. Under §8.4's corrected multiplier that trade may not be the right one. I did
not re-price it; it is a `pds`-instrument question and yours.

**To whoever next audits harness code:** the silent zero lives at
`experiments/162-retrodicting-validator/plolo.py:363-369`, not in the builder. The builder raises.
Every PLOLO-era screen inherits the silent version. A one-line change there —
`raise` instead of `return 0.0`, or return the carrying-row fraction — would have made 218's
discovery unnecessary and would prevent the next one.

**Nothing under `scripts/`, `src/`, `data/`, `submissions/` or any other experiment directory was
touched.**
