# 322 — Does a BROADER source stack degrade better on a disjoint panel? Measured on 233's own instrument.

**322 (broad-source October insurance): NOT WORTH IT** — retention Δ on the channel that
actually loses is **≤ 0**: `kolf2`'s Δ`pds` contribution is **−0.01017 today, +0.00000 in October
as the bank is built, and −0.00006 / −0.00468 in October after the bank is re-extracted so it
carries 251 / 111–153 rows** — 0/3 folds positive at the amplitude-matched dose on **all three**
panels, in **both** constructions. Today's cost is Δ`pds` **−0.01017** (residual form) /
**−0.00758** (264's base-donor form, amplitude-matched), ≈ **−0.0038 / −0.0028 overall** by 233
§2's own `pds` route. The **cosine** half does retain (+0.84 % → +1.53 / +2.00 % relative
`cos_agg`, **210 %**) — but cosine is the channel 233 measured at **93 % retention**, i.e. the one
that was never losing. **Breadth is insurance for the half that is not in danger.**

**Date 2026-09-08, 11:19–12:35 UTC (~1 h 15 min of a 5 h box).** **Spend $0** — no cloud, no
`cell-eval2`, no `vcc`, no build, nothing written under `submissions/`, `data/`, `scripts/`,
`docs/` or any other experiment's directory. Network: **26.5 MB** of commit-pinned HTTP-Range
reads from `raw.githubusercontent.com` (the KOLF2.1J atlas, §4). One deviation to declare:
`from oct233 import deranged` caused CPython to write
`experiments/233-october-retention/__pycache__/oct233.cpython-312.pyc`; it was **deleted** and a
`find -newer` audit over `experiments/ data/ scripts/ docs/ submissions/ src/` now returns nothing
outside this directory. Memory: 233's gate (`memory_pressure` ≥ 40 % free **and** free+inactive
≥ 20 GiB) checked at every launch — 79–80 % free, 32.0–32.4 GiB free+inactive throughout; heavy
jobs strictly serial; `prepare_main_source`'s 8 GiB cap **fired once** and the run aborted rather
than swapping (fixed by preparing the bank before loading the 2.6 GiB operator and running one
`(panel, fold)` per process, 264's rule).

---

## 0. The mandate's premise is inverted, and that is the result

> *"`kolf2` covers **273 of 300** panel targets — far broader than `dld1_lowmoi` (67) or
> `gse295214_h1` (64)."*

Both halves of that sentence are true and they are about **different quantities**. Measured
(`census322.py`, from each bank's own `perts.txt`):

| | `kolf2` | `dld1_lowmoi` | `gse295214_h1` |
|---|---:|---:|---:|
| **panel** coverage — of TODAY's 300 | **273** | 67 | 64 |
| **source** breadth — perturbations the bank contains AT ALL | **273** | **4,715** | **2,354** |
| coverage on 233's `real` draw | **0** | 67 | 57 |
| coverage on 233's `pess` draw | **0** | 73 | 34 |

`kolf2` is the **broadest non-production source on today's panel and, of these three, much the
narrowest as a source**: all 273 of its rows *are* panel rows (it ranks 6th-narrowest of the 21
banks on disk by perturbation count, and 3rd-broadest by panel overlap). `dld1_lowmoi`, the
mandate's "narrow" contrast, is **17.3× broader** as a source and is one of only three banks on
disk whose coverage *rises* on a disjoint draw (the other two are the Orion production donors).

Over all 21 banks on disk with a non-empty perturbation list, **which axis predicts October
coverage is not close**:

| predictor of coverage on a disjoint draw | Spearman ρ, `real` | `pess` |
|---|---:|---:|
| coverage on **today's** panel | +0.476 | +0.482 |
| **log₁₀ (perturbations in the source)** | **+0.940** | **+0.933** |

`docs/DEF-RUNBOOK.md` §8.3 already records the `kolf2` 273 → 0 collapse and the rule *"never rank
a source by today's overlap again"*; 233's KA-F reproduced it. **What was never done is pricing
the arm** — the runbook rules `kolf2` out on today's price and on October's coverage, but no
experiment had measured what a broad source is worth *on a disjoint panel*. That is §3–§6.

---

## 1. Known-answer rows. **Six gates, six passes.**

Every number below is produced by 233's own modules (`panels233`, `fold233`, `score233`,
`oct233.run`'s body) reading 233's own artefacts — the nine re-fitted operators, the nine fold
caches and the three re-extracted CD4 bank sets, all still on disk (12 GB) and opened read-only.

| # | gate | expected | measured | verdict |
|---|---|---|---|---|
| **KA-1** | 233's coverage census, 15 sources × 3 panels (R233 KA-F) | 25 published rows | **15/15 source rows exact**, incl. `kolf2` 273→0/0, `feng_ipsc` 182→0/0, `k562_gwps` 269→300/150 | **PASS** |
| **KA-2** | the shipped stack's `cos_agg` / `pds` on all nine (panel, fold) cells | R233 §3.1's table | **9/9**, worst \|Δ\| **4.9e-06** (R233 quotes 5 d.p.) | **PASS** |
| **KA-3** | `controls322.py`'s own BASE == the screen's BASE | 0 | max \|Δ\| **0.0e+00**, 9 cells | **PASS** |
| **KA-4** | `basedonor322.py`'s two-donor consensus == the screen's BASE | 0 | max \|Δ\| **0.0e+00**, 9 cells | **PASS** |
| **KA-5** | **233 §5's own ARM rows** — `rpe1` at doses 0.20 and 0.34 plus its derangement, on both disjoint draws — re-run through 322's copy of `oct233.run()` | 0 | **20/20 rows at 0.0e+00** against `E233/work/oct233.json` | **PASS** |
| **KA-6** | the KOLF2.1J range-extractor, re-run against **today's** panel, must reproduce the frozen bank | byte-identical | **5/5 files identical** (`delta.npy` `a8ed5c7e…`, `b_ctrl.npy`, `genes.txt`, `perts.txt`, `gene_map.csv`) | **PASS** |

KA-5 is the load-bearing one: it checks the **arm** path, not only the base, and it is bit-exact
against 233's cached Phase-3 rows. KA-6 licenses §4.

**Not a known-answer row, and labelled as such.** 233 re-fits every operator with its own recipe
(no essential prune, no 12,000-row cap) and measures that recipe change at −3.5 / −5.8 / −1.8 %
`cos_agg` and +0.0018 / −0.0029 / −0.0046 `pds` on today's panel (R233 "Bridge"). So §5's `abc`
base-donor numbers are **264's construction on a different operator** and are expected to agree in
sign and order, not to the digit. They do: 264's (a) − (c) reads −0.056 `pds` by the letter and
**−0.0056 amplitude-matched, 0/3 folds**; mine read **−0.0802** and **−0.00758, 0/3**. The
whitening pathology 264 named reproduces independently: `kolf2`'s rows are amplified **49× / 76× /
166×** through an operator fitted without them, against **6.2–6.8×** for the production donors.

---

## 2. Instrument, unit, and what is *not* measured

**Form.** 218's un-whitened post-whitening residual on top of the shipped stack (`VCC_ALPHA`
0.982, γ 0.65, λ 0, CD4 `f` 0.34, H1 union `g` 0.28, PerturbMatch 50:50), dose solved on the
source's **carrying** rows, scored on 233's three cancer PLOLO folds after the renderer and the
scorer's read-back. This is the form 233 §7 item 3's argument is about: CD4 and the H1 union are
both post-whitening residual sources sized on rows, and they are the two components whose
contributions §4.1 of R233 tabulates. 264's alternative — `kolf2` as an extra **base donor** — is
§5.

**The resampled unit is the PANEL DRAW.** 233 built two (`real`, `pess`) plus the control (`abc`).
**n = 2** bounds every claim here. Three folds are three re-readings of one draw, not three draws:
fold spread is reported as a range and never as an interval on the draw. No interval is quoted on
a retention ratio.

**Not measured.** The *context* half of Oct 22 (D/E/F are new cell lines) — no local instrument
reaches it, as R233 §0 states. And 233's panels are **random** draws while Arc's is **selected**
(R233 §6: a random draw's model-free transferability ceiling is 1.3–2.0× below Arc's), so every
*cosine* and *`pds`* retention figure inherits 233's lower-bound status. **The coverage facts do
not** — `kolf2`'s frozen bank holds 273 rows that are all panel rows, so it carries exactly zero on
*any* panel disjoint from today's, selected or random.

---

## 3. The retention table — contribution on a disjoint draw / contribution today

Contribution = arm − the shipped stack on the same fold, same operator, same CD4 banks. `shuf` is
the same source **whole-row deranged** at the same dose (233's null: magnitudes held, recipients
scrambled). Dose g = 0.34, 233's headline dose. Full tables, per fold and per dose, in
`work/tables322.md`.

| source | perts in bank | panel | carrying rows (K562 / HCT116 / HEK293T) | rel `cos_agg` | its `shuf` | Δ`pds` | folds > 0 | its `shuf` |
|---|---:|---|---|---:|---:|---:|---:|---:|
| **`kolf2`** (as built) | 273 | abc | 249 / 264 / 268 | **+0.84 %** | −3.85 % | **−0.01017** | **0/3** | −0.01092 |
| | | real | **0 / 0 / 0** | **+0.00 %** | — | **+0.00000** | — | — |
| | | pess | **0 / 0 / 0** | **+0.00 %** | — | **+0.00000** | — | — |
| **`dld1_lowmoi`** | 4,715 | abc | 60 / 66 / 65 | **+2.38 %** | −1.22 % | **+0.00288** | **3/3** | −0.00302 |
| | | real | 67 / 67 / 67 | **+1.27 %** | −1.07 % | +0.00023 | 1/3 | +0.00031 |
| | | pess | 38 / 72 / 72 | **+3.03 %** | −0.84 % | +0.00009 | 1/3 | −0.00192 |
| `telohaec` | 2,167 | abc | 58 / 62 / 63 | −0.50 % | −0.42 % | −0.00209 | 1/3 | −0.00047 |
| | | real | 42 / 42 / 42 | −0.58 % | −0.45 % | −0.00185 | 1/3 | +0.00056 |
| | | pess | 26 / 39 / 39 | −0.53 % | −1.34 % | −0.00093 | 1/3 | −0.00372 |
| `igvf_cardiomyocyte` | 1,888 | abc | 86 / 86 / 86 | +1.12 % | −1.75 % | −0.00038 | 1/3 | −0.00184 |
| | | real | 60 / 60 / 60 | −1.05 % | −1.79 % | +0.00200 | 3/3 | −0.00053 |
| | | pess | 33 / 35 / 36 | +6.18 % | −3.02 % | +0.00117 | 2/3 | −0.00100 |

**Retention (disjoint mean / today):**

| source | retention of the **cosine** contribution | retention of the **`pds`** contribution |
|---|---:|---:|
| `kolf2` as built | **0 %** (contributes exactly nothing) | **0 %** |
| `dld1_lowmoi` | **90 %** (+2.38 % → +1.27 / +3.03 %) | **6 %** (+0.00288 → +0.00023 / +0.00009) |
| `telohaec` | 111 % (a loss retained as a loss) | 66 % (of a loss) |
| `igvf_cardiomyocyte` | not readable — today's contribution changes sign across folds and draws | not readable (today's Δ`pds` is −0.00038, a ratio on it is meaningless) |

**The two channels come apart, and the split is the finding.** `dld1_lowmoi` is the clean case: its
coverage is essentially **panel-invariant** (67 → 67 / 73), so anything that changes is not
coverage. Its **cosine** contribution keeps **90 %**. Its **`pds`** contribution keeps **6 %** and
lands on its own derangement (+0.00023 vs shuffle +0.00031 on `real`; 1/3 folds positive on both
draws against 3/3 today). This is 233 §4.2's attribution — *"coverage explains 78 % of the cosine
loss but only 18 % of the `pds` loss"* — reproduced from the opposite direction: **hold a source's
coverage fixed and its `pds` value still evaporates on a panel Arc did not select.**

---

## 4. The steel-man: the zero is an artefact of OUR EXTRACTION, and it is fixable

`data/external/kolf2` was built by `experiments/055-kolf-direct/extract_kolf.py` as a **range
extraction against today's panel** — its own `meta.json` records `range_request_count 273`,
`range_bytes_downloaded 10,539,438`, `full_matrix_bytes_avoided 392,507,202`. The upstream
KOLF2.1J atlas (commit `8f84531a…`, hash-pinned in `provenance.json`) holds **10,167
perturbations over 38,606 genes**, and it covers **251 / 300** of 233's `real` draw and **159 /
300** of its `pess` draw. So the runbook's *"`kolf2` goes to EXACTLY ZERO … #28 in October"* is a
property of the artefact on disk, **not of the source**, and undoing it costs **9.7 MB / 6.1 MB**
of range requests and **0.15 GiB** of RSS.

I re-extracted it for both draws (and for `abc`, which is **KA-6**: byte-identical to the frozen
bank), and re-ran both constructions.

| `kolf2_refit` | panel | carrying rows | rel `cos_agg` | its `shuf` | Δ`pds` | folds > 0 |
|---|---|---|---:|---:|---:|---:|
| residual form, g = 0.34 | abc | 249 / 264 / 268 | +0.84 % | −3.85 % | **−0.01017** | **0/3** |
| | real | **251 / 251 / 251** | **+1.53 %** | −2.19 % | **−0.00006** | 1/3 |
| | pess | **111 / 151 / 153** | **+2.00 %** | −3.05 % | **−0.00468** | **0/3** |

**Cosine retention 210 %** — a genuinely broad source is worth *more* in October than today, and
the effect is source-specific (its derangement reads −2.2 / −3.1 % where the source reads
+1.5 / +2.0 %). **`pds` retention 23 %, of a negative number**: the contribution is ≤ 0 on every
panel. And the `shuf` column says why the `pds` movement is not biology — on `abc` the source
reads −0.01017 against its derangement's −0.01092, and on `pess` −0.00468 against −0.00461.
**The `pds` cost is the generic amplitude cost of adding a residual at this dose; the cosine gain
is the source.**

---

## 5. `kolf2` as an extra BASE DONOR — 264's construction, on all three panels

264's construction transplanted operation for operation (`E062.prepare_main_source` → restrict to
the operator's columns → weight → `E062.build_baseline`), with 264's two weights: `equal` (1.0)
and `matched` (the fold-constant scalar equalising `kolf2`'s un-preserved whitened median norm
with the production donors', `lib264.stage3:209`). Means over the three folds:

| bank | panel | `mean3:equal` rel `cos_agg` / Δ`pds` | `mean3:matched` rel `cos_agg` / Δ`pds` | matched folds with Δ`pds` > 0 |
|---|---|---|---|---:|
| frozen (273 rows) | abc | −10.68 % / **−0.08020** | +1.73 % / **−0.00758** | **0/3** |
| | real | +0.00 % / +0.00000 | **identity** (0 rows) | 0/3 |
| | pess | −0.96 % / +0.00132 | **identity** (0 rows) | 0/3 |
| re-extracted (10,167-row atlas) | abc | −10.68 % / −0.08020 | +1.73 % / −0.00758 | **0/3** |
| | real | −16.46 % / **−0.05312** | +0.58 % / **−0.00706** | **0/3** |
| | pess | −22.98 % / **−0.05005** | −3.31 % / **−0.00844** | **0/3** |

**0/3 folds positive on `pds` in every panel × construction cell that has rows.** The
amplitude-matched loss is essentially *flat across panels* (−0.00758 today, −0.00706 / −0.00844
in October) — the arm does not degrade better tomorrow; it degrades identically.

Two by-products worth recording.

- **A zero-coverage donor is not always the identity.** `E062.build_baseline:311-314` fills rows no
  donor carries with `EFFECT_SCALE × mean(mean_on_axis) / ‖·‖`. Adding a donor changes that mean,
  so on `pess`/K562 — where **7 of 150** targets have no production carrier — the frozen,
  zero-coverage `kolf2` still moved `cos_agg` **−2.87 %** and `pds` **+0.00397**. On `real` (all
  300 targets carried by both donors) it is exactly the identity, max \|Δ\| **0.0e+00**. **Naming
  a vacuous source in `VCC_UNION_SOURCES` is a hazard on the uncovered rows, of uncontrolled
  sign** — a second face of 235's silent zero, which `make_submission_pm_cd4.py:864` already
  refuses for the union but which `build_baseline` does not see.
- The whitening pathology 264 measured is **worse, not better, in October**: amplification
  48–157× on `real` and 54–157× on `pess`, production donors 6.4–7.3×.

---

## 6. The controls — coverage or identity?

### 6.1 Control 1: the coverage ladder (`dld1_lowmoi`'s own rows cut at random to k targets, 2–3 seeds × 3 folds)

| k targets | abc rel `cos_agg` | real | pess | disjoint / abc | **abc Δ`pds`** | **real** | **pess** |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 20 | +1.04 % | +0.61 % | +1.32 % | 93 % | **+0.00079** | +0.00051 | −0.00108 |
| 40 | +1.40 % | +0.52 % | +0.69 % | 43 % | **+0.00148** | +0.00006 | +0.00051 |
| 60 | +2.01 % | +1.06 % | +2.73 % | 94 % | **+0.00197** | +0.00019 | +0.00001 |
| all (67 / 67 / 72) | +2.38 % | +1.27 % | +3.03 % | 90 % | **+0.00288** | +0.00023 | +0.00009 |

**This is the mechanism in one table.** On **today's** panel Δ`pds` is **monotone in coverage**
(+0.00079 → +0.00098 → +0.00148 → +0.00197 → +0.00288 across k = 20/26/40/60/all). On **both
disjoint draws** it is **flat at ≈ 0 at every rung** (largest value anywhere: +0.00051). The
cosine column, by contrast, rises with k on all three panels and retains 90–94 % at the top of the
ladder. **Coverage buys `pds` only on Arc's selected panel; it buys cosine on any panel.**

### 6.2 Control 2: coverage held fixed (both sources restricted to the rows BOTH carry)

Cells are `rel cos_agg / Δpds`, mean over three folds (and seeds for the random columns).

| panel | k | `dld1_lowmoi` **common** | `igvf_cardiomyocyte` **common** | `dld1` **random** k | `igvf` **random** k |
|---|---:|---:|---:|---:|---:|
| abc | 26 | **+1.02 % / +0.00104** | −0.24 % / −0.00102 | +0.65 % / +0.00098 | −0.02 % / +0.00007 |
| real | 12 | **+0.02 % / −0.00126** | −0.27 % / −0.00001 | −0.06 % / −0.00010 | −0.42 % / −0.00034 |
| pess | 9 | **+0.26 % / −0.00350** | +0.09 % / +0.00222 | +1.19 % / +0.00087 | +0.56 % / +0.00077 |

At **identical coverage on identical targets**, `dld1_lowmoi` beats `igvf_cardiomyocyte` on
`cos_agg` on **3/3 panels** (+1.26 / +0.29 / +0.17 points). So the cosine channel is **not** pure
reach — source identity carries a real part of it. And the **common** rows are not special: at
matched k the random restriction reads within ±0.9 points of the common one and is *larger* on
`pess`, so no "these are the targets everybody covers" effect exists to exploit.

---

## 7. What this says about October

1. **The specific question the mandate asked is answered NO, twice over.** `kolf2` as built
   delivers **exactly zero** in October while costing Δ`pds` −0.01017 (residual) / −0.00758
   (base-donor, matched) today. `kolf2` re-extracted so that it *does* carry October's rows still
   delivers Δ`pds` ≤ 0 on both draws and 0/3 folds in the base-donor form. A source that loses
   today does **not** degrade better tomorrow — this one degrades the same.

2. **"Buy breadth, not a new mechanism" (R233 §7.3) needs one word added: breadth buys COSINE.**
   233 measured cosine retention at 93 % and `pds` retention at 64 %, with the board-unit
   projection driven by `pds` (scaled 0.7454 → 0.4793 vs `mse` 0.1796 → 0.1560, R233 §3.3). Every
   measurement here says the same thing about *additions* to the stack: a broad source's cosine
   contribution retains 90–210 %, its `pds` contribution retains ≤ 6 % and sits on its own null.
   **The insurance is real and it is sold on the channel that is not at risk.**

3. **Rank sources by the size of the source, not by today's overlap** — the runbook already says
   this; here is the number: log₁₀(perturbations) predicts October coverage at Spearman **+0.94**,
   today's panel overlap at **+0.48**, over 21 banks.

4. **Two operational corrections, both cheap, neither competition-facing.**
   - `DEF-RUNBOOK` §8.3's *"`kolf2` … #28 in October"* is a property of our **extraction**, not the
     source. The atlas covers 251/159 of the two draws and `extract_kolf.py` reproduces the frozen
     bank byte-for-byte (KA-6), so re-extracting on Oct 22 is ~10 MB and 0.15 GiB. **It is still
     not worth doing** on these numbers — recorded so the option is not rediscovered under time
     pressure and then taken.
   - `E062.build_baseline`'s uncovered-row filler means a **zero-coverage donor is not inert**: on
     `pess`/K562 (7 uncovered targets) it moved `cos_agg` −2.87 % and `pds` +0.00397 with **zero**
     carrying rows. 235's guard covers the union sources; the donor list has no equivalent.

5. **What would actually be worth buying, and it is already priced.** 233 §5's four October-only
   sources are the mirror image of `kolf2` — 0/300 today, 30–48 in October, +5.5…+17.7 % relative
   cosine each, dose-monotone, shuffle-controlled. They are the breadth purchase; `kolf2` is not.
   Their `pds` legs (+0.00068…+0.00428) are of the same small size as everything else in this file.

**Nothing here is a candidate.** No arm is proposed, no slot is armed, no preregistration written,
no package built. Every magnitude above is a **local PLOLO reading on 233's synthetic panels**
(`docs/DEF-RUNBOOK.md` §5's two bundles) and must not be quoted as a delivered board value
(Correction 5, `docs/LOG.md`). The two board-unit figures in the verdict line are 233 §2's own
`pds` route (`Δpds / 0.452 / 6 × 1.011`) applied to a local Δ`pds`, and are indicative sizings, not
predictions.

---

## 8. Falsifiers

| # | falsifier | outcome |
|---|---|---|
| 1 | 233's KA must reproduce before any number of mine is trusted | KA-1…KA-5, six gates, **all pass**; KA-5 is bit-exact on 233's ARM rows, not only its base |
| 2 | panel disjointness asserted, not assumed | `census322.py` asserts \|abc ∩ real\| = \|abc ∩ pess\| = 0 at runtime; 300 unique symbols per panel |
| 3 | a measured null for every arm | every source is run whole-row **deranged** at the same dose; on `kolf2_refit` the derangement *matches* the `pds` movement (−0.01092 vs −0.01017 on abc; −0.00461 vs −0.00468 on pess) and is opposite on cosine — which is how the cosine gain is shown to be the source and the `pds` cost shown not to be |
| 4 | the "broad source" claim must survive the source's own extraction being fixed | §4: re-extracted from the 10,167-row atlas, `kolf2` carries 251/111–153 rows and STILL never buys `pds`; KA-6 licenses the re-extraction |
| 5 | coverage vs identity must be separated, not asserted | §6.1 ladder (coverage varied, source fixed) and §6.2 common-row restriction (coverage and targets fixed, source varied): cosine responds to both, `pds` responds to neither on a disjoint panel |
| 6 | two draws, not one; report both if they disagree | reported per draw throughout; they agree in direction on the `pds` verdict in every cell and disagree in magnitude on cosine (`real` +1.53 %, `pess` +2.00 % for `kolf2_refit`), which is why no interval is quoted on n = 2 |
| 7 | the narrow/broad contrast must not be a single-source anecdote | five candidate banks (four sources + the re-extracted `kolf2`) × 3 panels × 3 folds at two doses plus derangements, the two controls and the base-donor arm: **348 rows in total — 285 scored arms, 51 bases, 12 structural zero-row rows** — plus the 21-bank census correlation |
| 8 | leak audit | no candidate shares a cell line with any fold's truth (233's `FORBIDDEN` guard retained and empty for these four); `resid218.PANEL_SET` re-bound per panel through `score233.bind_panel` (R233 falsifier 8's silent leak); every panel uses its **own** re-fitted operator, verified by `axis_sha256` on load, and its **own** re-extracted CD4 banks |
| 9 | writes confined to this directory | `find -newer MANDATE.md` over `experiments/ data/ scripts/ docs/ submissions/ src/` returns nothing outside `experiments/322-october-breadth/` (one `.pyc` was created in 233's `__pycache__` and removed; declared at the top) |

---

## 9. Files

`lib322.py` (233's modules on the path, the extra source registrations, the two control
transforms, 233's memory gate) · `census322.py` (KA-1, the two-senses-of-breadth table, the
Spearman comparison) · `breadth322.py` (the retention screen — `oct233.run()`'s body with a
per-panel source registration and 252's `SilentZeroError` recorded as the 0.0 it is) ·
`basedonor322.py` (264's construction, panel-parameterised; `run_basedonor.sh` /
`run_basedonor_refit.sh` run one `(panel, fold)` per process) · `controls322.py` (the coverage
ladder and the common-row restriction) · `analyse322.py` (all tables → `work/tables322.md`) ·
`run_ka5.sh` (KA-5).

Artefacts in `work/`: `census322.json`, `breadth322.json` (105), `refit322.json` (36),
`basedonor322.json` (27), `basedonor_refit322.json` (27), `controls322.json` (129),
`ka5_322.json` (24), 
`contrib322.json`, `tables322.md`, and `kolf2-refit/{abc,real,pess}/` (the three re-extractions,
77 MB total including logs). 233's own `work/` was opened **read-only** throughout.
