# 207 — The CD4 dose was never swept against cosine. Swept, it wants to be **3.4× larger**, and the FORM should stay a residual, not become a library.

**Date:** 2026-09-03 · four PLOLO folds (162's instrument, imported verbatim) · **no `vcc`, no
`cell-eval2`, no build, no submission, no write outside `experiments/207-cd4-dose/`.** Peak RSS
**4.61 GB** (`/usr/bin/time -l`, single-fold reference run). 12 scripts; **744 scored
configurations** in `work/` (232 dose/control + 456 form + 28 form-falsifier + 28 whitened-residual).

---

## 0. The answer in ten lines

1. **206's H1 is CONFIRMED, and its three maximisers genuinely differ.** Over a 25-point dose
   sweep the **cosine** maximiser is `f = 0.80` (+31.5% over shipped), the **`pds`** maximiser is
   `f = 0.20`, and **`q_w`**'s is `f = 0.80`. The shipped `f = 0.10` is below all three.
2. **But 206's and the mandate's *premise* is a knob confusion, and correcting it makes the case
   stronger.** The residual does **not** ship at "effective fraction 0.0528". `0.0528` is a
   builder *diagnostic* (`effective_f_at_cd4_fraction`: the `f` you would get at the legacy
   `VCC_CD4_FRACTION = 0.05`). **a14 shipped `fraction = 0.05` — that is 062's dose. a15 shipped
   `VCC_CD4_TARGET_F = 0.10`, which the builder solves to `fraction = 0.09473`.** So the a14→a15
   board arm that 197 measured at **+6.8% cosine was itself a 1.89× DOSE INCREASE**, not the
   introduction of CD4. We already have a two-point board dose-response curve and its slope is
   positive.
3. **Every known-answer gate passes**, including two the mandate did not ask for: 062's three
   headline LODO numbers reproduce on this independent instrument at **0.99× / 1.14× / 1.18×**.
4. **Falsifier 1 (target-shuffled CD4) passes decisively at every dose.** Shuffled CD4 makes the
   cosine *monotonically worse* (−0.0047 at `f = 0.8`) and `pds` worse by −0.071. The gain is
   target-specific biology, not amplitude or geometry.
5. **The mechanism is that CD4 enters UN-WHITENED.** A CD4 residual pushed through the same
   operator first (`wres`) is **flat on cosine at every dose (+0.5 to +0.9%)**. 099's
   self-residual control recovers **44%** of the cosine gain but **loses 0.0078 `pds`**, where CD4
   gains 0.0070. So ~56% of the cosine and **100% of the `pds`** is CD4-specific information.
6. **The FORM should NOT change.** GISL's library construction (CD4 as a fourth equal-weight
   consensus member, whitened with the rest) **loses on every one of 15 weights and 42 composed
   settings**, out-of-sample price **−0.0081**, and its `nmae` blows up **+0.104 raw** because the
   delivered amplitude collapses to **0.52×** even at a re-solved α. The residual form is the only
   one of four that survives leave-one-fold-out selection.
7. **The gain is LINEAGE-GENERIC, which is the answer that matters for D/E/F.** It is positive and
   of the same order on four unrelated **non-T** recipient lines (+12.3 / +11.2 / +7.6 / +15.2% at
   the shipped dose). The mandate's per-context test **cannot be run at all** — `nadig_jurkat`
   covers **0/300** panel targets and the CD4 banks are panel-restricted, so no surface pairs a
   T-lineage recipient with the CD4 donor. There is a weak lineage *gradient* in CD4's raw
   alignment (highest to the one haematopoietic line), but the collectable component does not
   depend on the recipient's lineage.
8. **Chosen operating point: `VCC_CD4_TARGET_F = 0.34`** (3.4× shipped). Selected by three
   criteria that agree: the argmax of the measured six-member price, the modal leave-one-fold-out
   choice, and the largest dose at which `pds` is not resolvably worse (1.23σ < 1.645σ).
9. **Priced: `pds` −0.000685, `mse` +0.005206, `nmae` −0.000152 → +0.00437 overall**
   (board-calibrated cosine transfer: **+0.00572**). Out-of-sample leave-one-fold-out:
   **+0.00340 / +0.00502**. Ordering precision (the half of `reach`/`fid` that is readable)
   **improves +0.00757**.
10. **Gate: this is a SHAPE-ONLY arm (+0.0019), measured not assumed** — delivered median row norm
    **1.0009×** a15's, delivered `A` within **1.0%**, α unchanged and its re-solve a literal no-op.
    **It clears by 2.3× conservatively and by 1.8× out-of-sample. SHIP IT.**

---

## 1. Rule 11 — every operator, by path and sha256

| fold (held = truth) | sources | operator | sha256 | shape | wcols | targets |
|---|---|---|---|---|---:|---:|
| `k562_gwps` | hct116 + hek293t | `experiments/062-cd4-post-whiten/cloud-output/whiteners/held_k562_gwps.npz` | `13d988fc…3bc55` | 17,806² | 17,806 | 269 |
| `orion_hct116` | k562 + hek293t | `…/held_orion_hct116.npz` | `7cecbf7f…47cbc` | 17,942² | 17,942 | 290 |
| `orion_hek293t` | k562 + hct116 | `…/held_orion_hek293t.npz` | `97e486e0…c8475f2af1` | 17,942² | 17,942 | 294 |
| `kolf2` | all three | **`data/donors/whitener_union_pruned.npz`** | **`e2f65ecf…1db123`** | 17,942² | 17,942 | 273 |

The production operator is **`e2f65ecf`**, basis **29,394 × 17,942** (`whitener_union_pruned.json`:
`basis: [29394, 17942]`, `basis_prune: essential`, rows kept 7,784 / 10,807 / 10,803). This
confirms 193's correction; **`AGENTS.md` §4.2's "12,000 × 7,760" matches no file on disk.**
All **27** operators on disk and their shas are enumerated in §9 so no future screen has to
re-derive them.

The three held operators reproduce 062's own recorded hashes exactly (062 RESULTS §"Preservation
and hashes"), which is the provenance chain from 062's LODO to this sweep.

---

## 2. The dose frames — three of them, and the mandate inherited the wrong one

`scripts/make_submission_pm_cd4.py` carries **two** knobs with the word "fraction" in them plus a
diagnostic, and 206 §6.4 and this mandate both quote the diagnostic as the shipped dose.

| name | what it is | shipped value |
|---|---|---:|
| `VCC_CD4_FRACTION` (`fraction`) | **062's knob.** `add_post_whiten_residual` adds `fraction × cd4.mean_effect`. 062 selected **0.05**. | not used since a14 |
| `VCC_CD4_TARGET_F` (`f`) | **099's knob.** Ratio of the added residual's *median row norm* to the whitened base's. | **0.10** |
| `effective_f_at_cd4_fraction` | a **diagnostic** the builder prints: the `f` you would get at the legacy `fraction = 0.05`. | 0.05278 |

Read off the finished artifact (`experiments/180-def-readiness/work/abc-knee.h5ad.meta.json`,
`effects.cd4`), verbatim:

```json
"fraction": 0.09472706196346513,
"residual_frame": {"cd4_fraction_constant": 0.05,
                   "effective_f_at_cd4_fraction": 0.05278322684523277,
                   "f_per_unit_fraction": 1.0556645369046553,
                   "fraction_source": "VCC_CD4_TARGET_F",
                   "solved_fraction": 0.09472706196346513,
                   "target_f_requested": 0.1}
```

`experiments/191-rc-hybrid-build/build191.sh` line 22: `export VCC_CD4_TARGET_F=0.1`.

**Two consequences.**

- **062's selected dose is not what ships.** Production ships `fraction = 0.0947` — **1.89×**
  062's `0.05`. 062's `pds`-maximising selection was already overridden once, by 099, on `f`
  grounds. The mandate's sentence "062 selected that dose by maximising LODO `pds`" is therefore
  half true: 062 selected `0.05` on `pds`, and 099 moved it to `0.0947` for a reason unrelated to
  the objective.
- **197 §1.2's +6.8% is a dose DERIVATIVE, not an on/off measurement.** 162's own arm table
  (`run_arms.py` lines 41–42) has a14 at `cd4_frac = 0.05` and a15 at `cd4_f = 0.10`. Both carry
  CD4. So the board has already measured that **raising the dose 1.89× bought +6.8% cosine and
  +0.007257 scaled `pds`**, and the sweep below is a continuation of a curve the board has two
  points on, not an extrapolation from one.

Everything below reports the dose in the **`f`** frame (the shipped knob) with the solved
`fraction` alongside.

---

## 3. Phase 0 — known-answer gates. All pass.

`analyse207.py`, `work/analyse207.json`.

| gate | expected | measured | verdict |
|---|---:|---:|---|
| **KA-1** a14 → a15 PLOLO Δ`pds`, scaled, 4-fold mean | 162's **+0.006656** | **+0.006109** | **PASS** (−8.2%) |
| **KA-2** a15 → a23 mass-fix price, scaled | 162's **+0.004014** | **+0.004014** | **PASS (exact)** |
| **KA-3** a14 → a15 relative aggregate-cosine gain | board **+6.83%** (197 §1.2) | **+5.42%** | **PASS** (see §4.3) |
| **KA-4** `a15` arm ≡ `f0.1` arm, two code paths | identical | identical to 1e−12 | **PASS** |
| **KA-5** `lib0` (library at w = 0) ≡ `f0` | identical | identical, 4/4 folds | **PASS** |

**KA-2 lands exactly on 162's published +0.004014**, which is the strongest single statement that
this is 162's instrument and not a re-implementation of it. KA-1 is 92% of 162's value; the 8%
is the two experiments' bootstrap draw over the same folds, and 162's own realised value was
+0.007257, i.e. **both readings understate the board.**

### 3.1 The extra gate the mandate asked for — 062's own headline

062 reported HCT116 **+0.00482**, HEK293T **+0.00242**, frozen K562 **+0.00441** at `fraction`
0.05, on an *effect-level* LODO over 242 shared targets with no renderer. This experiment is a
*rendered*, challenge-frame, ~290-decoy retrieval. Measuring the identical configuration
(`a14` = `fraction` 0.05 exactly) against CD4-off:

| held line | 062 LODO Δ`pds` | this instrument | ratio |
|---|---:|---:|---:|
| Orion HCT116 | +0.00482 | **+0.005508** | **1.14×** |
| Orion HEK293T | +0.00242 | **+0.002856** | **1.18×** |
| K562-GWPS (062's post-selection holdout) | +0.00441 | **+0.004374** | **0.99×** |
| KOLF2 | — (no such fold in 062) | +0.004238 | — |

**All three reproduce in sign, in ordering, and to within 1–18% in magnitude across two
completely different scoring constructions.** The mandate asked which numbers the axis reproduces
and which it cannot: **it reproduces all three.** What it cannot reproduce is 062's *bootstrap
interval* (`[+0.00129, +0.00761]`), because that was a 20,000-replicate paired-target resample on
242 targets in the effect frame and the target sets differ per fold here.

---

## 4. Phase 1 — the dose sweep. Cosine is the primary endpoint.

`sweep207.py` → `work/sweep207.json` (59 arms × 4 folds). Baseline `f0` = CD4 **off**.
`cos_agg = sqrt(1 − sse)` is the **energy-weighted aggregate cosine**, the same functional 197
inverts off the board (`X / (√ΣD · mean_c √A_c)`). CIs are paired-**target** bootstraps
(4,000–8,000 replicates, seed 20260903) over all four folds jointly.

| `f` | `fraction`¹ | **cos_agg** | Δ vs `f0` | Δcos 95% CI | P>0 | `pds` | Δ`pds` | P>0 | `q_w` |
|---:|---:|---:|---:|---|---:|---:|---:|---:|---:|
| 0 | 0.0000 | 0.018091 | — | — | — | 0.735107 | — | — | 0.488256 |
| 0.02 | 0.0211 | 0.018510 | +0.000419 | [+0.00039,+0.00052] | 1.000 | 0.737156 | +0.002049 | 1.000 | 0.488537 |
| 0.05 | 0.0528 | 0.019121 | +0.001030 | [+0.00095,+0.00126] | 1.000 | 0.739579 | +0.004472 | 1.000 | 0.488636 |
| **0.10 ← SHIPPED** | 0.1057 | **0.020086** | +0.001995 | [+0.00178,+0.00238] | 1.000 | **0.742112** | +0.007005 | 1.000 | 0.488706 |
| 0.14 | 0.1479 | 0.020808 | +0.002717 | [+0.00237,+0.00317] | 1.000 | 0.743145 | +0.008039 | 1.000 | 0.488656 |
| 0.17 | 0.1796 | 0.021318 | +0.003227 | — | 1.000 | 0.743366 | +0.008260 | — | 0.488948 |
| **0.20 ← `pds` argmax** | 0.2113 | 0.021802 | +0.003711 | [+0.00311,+0.00419] | 1.000 | **0.743382** | **+0.008275** | 0.997 | 0.488788 |
| 0.24 | 0.2536 | 0.022404 | +0.004313 | — | 1.000 | 0.742702 | +0.007595 | — | 0.488605 |
| 0.28 | 0.2959 | 0.022958 | +0.004867 | [+0.00388,+0.00527] | 1.000 | 0.741874 | +0.006768 | 0.966 | 0.488877 |
| **0.34 ← CHOSEN** | 0.3593 | **0.023698** | +0.005607 | [+0.00251,+0.00354]² | 1.000 | 0.740275 | +0.005168 | — | 0.488959 |
| 0.40 | 0.4226 | 0.024330 | +0.006239 | [+0.00463,+0.00640] | 1.000 | 0.738481 | +0.003375 | 0.767 | 0.488986 |
| 0.55 | 0.5811 | 0.025490 | +0.007399 | [+0.00505,+0.00718] | 1.000 | 0.733126 | −0.001981 | 0.357 | 0.489247 |
| **0.80 ← cosine argmax** | 0.8453 | **0.026416** | +0.008325 | [+0.00503,+0.00756] | 1.000 | 0.723664 | −0.011443 | 0.039 | **0.489590** |
| 1.0 | 1.0566 | 0.026597 | +0.008506 | — | 1.000 | 0.716369 | −0.018737 | — | 0.489395 |
| 1.3 | 1.3736 | 0.026420 | +0.008329 | — | 1.000 | 0.707224 | −0.027883 | — | 0.489279 |
| 3.0 | 3.1699 | 0.024203 | +0.006112 | — | 1.000 | 0.681086 | −0.054020 | — | 0.488713 |
| 7.0 | 7.3964 | 0.022071 | +0.003980 | — | 1.000 | 0.664866 | −0.070240 | — | 0.488267 |

¹ the **solved** `fraction`, 4-fold mean. It is 11.6% larger here than in the challenge build for
the same `f` (`f_per_unit` 0.946 on PLOLO vs 1.0557 on the abc-knee artifact), because the two have
different targets and different whitened bases. **That difference is exactly what the `f` frame
normalises away, and it is why the dose must be swept and transferred in `f`, not in `fraction`.**
² the CI in that cell is against the SHIPPED point, not against `f0`; the full grid against both
baselines is in `work/analyse207.json`.

### 4.1 The three maximisers — the hypothesis under test, and it holds

| endpoint | argmax | value | vs the shipped `f = 0.10` |
|---|---|---:|---:|
| **aggregate cosine (197's route) — PRIMARY** | **`f = 0.80`** | 0.026416 | **+31.52%** |
| macro per-target cosine | `f = 0.80` | 0.026085 | +19.19% |
| **PLOLO `pds` (062's objective)** | **`f = 0.20`** | 0.743382 | +0.17% |
| **`|t|`-weighted `q_w` (187 §5)** | **`f = 0.80`** | 0.489590 | +0.18% |

**The maximisers differ by a factor of four in dose.** The member the dose was tuned on (`pds`)
turns over at `f = 0.20`; the member the direction actually pays through (cosine) keeps rising to
`f ≈ 0.8–1.0`. 206's H1 is confirmed as stated. The interior cosine maximum at `f ≈ 1.0` (the
residual's median row norm equal to the base's) rules out the trivial reading that "more is always
better" — at `f = 3.0` the cosine has fallen back below `f = 0.4`.

### 4.2 Why `pds` and cosine come apart — measured, not argued

It is **not** hubness. The mean off-diagonal row-row cosine *falls* with dose
(0.002256 → 0.000513 from `f = 0` to `f = 0.8`), i.e. the rows become *more* mutually dissimilar,
which is what `pds` normally wants. The shuffled control's `rr` falls identically, so the fall is
a property of adding any residual and carries no information. The turnover is the rank statistic
saturating: `pds` is a within-300 ranking and 197 §3 already measured that we extract 26% more
`pds` per unit of bulk cosine than the field, so `pds` under-rewards exactly this channel.

### 4.3 The transfer of the PLOLO cosine to the board — calibrated, not assumed

PLOLO's absolute cosine (0.0181–0.0264) is not the board's (0.079249); its truth is a single
donor line's pseudobulk and its prediction has two donors on three of four folds. The cosine is
carried across **relatively**, and the factor is measured on the one board point that exists for
this exact arm family:

```
a14 -> a15 :  board  +6.83%  (197 §1.2, 0.07418 -> 0.07925)
              PLOLO  +5.42%           (0.019053 -> 0.020086)
              => k_cos = 1.2598  (PLOLO UNDERSTATES the board by 26%)
```

Every price below is reported **twice**: at `k = 1.0` (uncalibrated, conservative — PLOLO's own
number carried straight across) and at `k = 1.2598`. **The conservative one is the headline.**
197's own row-to-row systematic on the board cosine is ±2% within an arm family and ±4% across
families, so `k_cos` itself carries roughly [0.9, 1.6]; `k = 1.0` sits at its low end.

---

## 5. Phase 2 — the FORM. The residual wins; the library form loses on every setting.

`library207.py` (15 library weights + 42 composed settings × 2 alphas × 4 folds = 464 rows),
`form207.py` (the whitened-residual decomposition), `libprice207.py`.

### 5.1 The three forms, priced against the shipped configuration

| form | best setting | Δcos_agg | Δ`pds` | **priced overall (k=1)** | **LOFO out-of-sample** |
|---|---|---:|---:|---:|---:|
| **post-whitening residual (incumbent form)** | **`f = 0.34`** | **+0.003612 (+18.0%)** | −0.001837 | **+0.00437** | **+0.00340** |
| CD4 as a 4th **library** member, whitened with the rest (GISL's form) | `w = 0.125` | −0.001780 | −0.003077 | −0.00372 | **−0.00806** |
| the two **composed** | `w = 0.25, f = 0.15` | +0.001033 | −0.003400 | +0.00022 | **−0.00195** |
| CD4 **whitened** then added as a residual | `f = 0.55` | −0.001975 | −0.003274 | −0.00407 | **−0.00560** |

**Not one library or composed configuration out of 57 beats the residual, and out-of-sample all
three alternative forms are negative.**

### 5.2 The mechanism, isolated — it is that CD4 must enter UN-WHITENED

The library form and the residual form differ in two ways at once (the member is whitened; the
member replaces consensus mass rather than being added and rescaled). `wres` isolates the first
by pushing CD4 through the *same* operator and *then* adding it at the same `f` solve:

| `f` | **un-whitened residual** cos_agg | rel | **WHITENED residual** cos_agg | rel |
|---:|---:|---:|---:|---:|
| 0.05 | 0.019121 | +5.7% | 0.018178 | +0.5% |
| 0.10 | 0.020086 | +11.0% | 0.018224 | +0.7% |
| 0.20 | 0.021802 | +20.5% | 0.018258 | +0.9% |
| 0.40 | 0.024330 | +34.5% | 0.018199 | +0.6% |
| 0.80 | 0.026416 | +46.0% | 0.017917 | −1.0% |

**Whitening CD4 destroys the entire cosine gain.** This is §4.2's own law read forwards:
whitening *lowers* the transfer cosine (measured there at −32%) while raising `pds`. The
post-whitening residual is the only place in the pipeline where information can be added without
paying that tax, and that is why the form matters more than the dose.

### 5.3 …and it is not merely whitening relaxation. 099's own control, run.

099 warned that adding back **the base mixture's own un-whitened effect** reproduces most of the
direction curve. Run here at matched `f`, with the same support mask and the same solve:

| arm at `f = 0.10` | Δcos_agg vs `f0` | Δ`pds` vs `f0` |
|---|---:|---:|
| **CD4 residual** | **+0.001995** | **+0.007005** |
| 099's self-residual (our own un-whitened effect) | +0.000873 (**44%** of CD4's) | **−0.007788** |
| target-shuffled CD4 | −0.000149 | −0.004641 |

**44% of the cosine gain is whitening relaxation, available from our own data. 56% of the cosine
and 100% of the `pds` is CD4-specific**, and the self-residual has the *opposite sign* on `pds`.
The two are not substitutes.

### 5.4 Why the library form fails, and why GISL's +0.01325 does not transfer

Three independent reasons, each measured:

1. **It gets whitened** (§5.2), which is where the cosine lives.
2. **It collapses the delivered amplitude.** At `w = 0.5` the *pre*-whitening median row norm
   falls to 0.855× and the re-solved α = 1.149 restores that — but the **delivered** median row
   norm is still **0.523×** a15's, because averaging in a fourth member shrinks rows
   non-uniformly and the renderer is nonlinear (170 §2 point 4's trap, in a new place).
3. **`nmae` therefore explodes:** raw 0.98476 → **1.08837** (+0.104), i.e. **−0.0288 overall**,
   with marginal refund `e = −1.195` (169 §5's threshold is `e > ½`). One member alone kills it.

**GISL's board-measured +0.01325 is not contradicted; it is not applicable.** They measured it on
a submission whose raw `mse` is ~5.0, i.e. with `mse` clamped to scaled 0 and — on the evidence of
that raw value — no spectral whitening at all. In a pipeline with no whitener there is no
distinction between "library member" and "un-whitened residual"; **their construction is closer to
our residual than to our library.** Read that way their number is corroboration for §5.2, not for
the library form.

### 5.5 The fold that dissents, and it is the one to expect

The library form is not uniformly bad: on `orion_hct116` `lib0.5` reads cos_agg 0.024304 against
`lib0`'s 0.021330 (+13.9%) and `pds` +0.0152. It collapses on `k562_gwps` (−36% cosine). The
asymmetry is explained by §7.2: on the K562 fold the two sources are the two Orion lines, which
share the X-Atlas long-gene technical axis (`AGENTS.md` §6.4: top-PC cosine 0.872), so that fold's consensus
is unusually well aligned already and the whitened CD4 member displaces it. The 4-fold mean is
162's pre-committed predictor and it is negative; the LOFO selection is negative on 3 of 4 held
folds. **Do not revive the library form on the HCT116 row alone.**

---

## 6. Phase 3 — the six-member price at the operating point

### 6.1 `nmae`, `reach`/`fid` ordering, and the delivered amplitude — and why NOT h1-test

**The mandate names h1-test. It cannot price this arm, and that is a measurement, not a
preference.** The CD4 banks on disk (`data/external/cd4t/cloud-compact/cd4_compact`) are
**panel-restricted** — 283 / 287 / 284 perturbations, all inside the 2026 300-gene panel.
h1-test's 85 targets intersect the panel in **5 genes**, and all three CD4 conditions cover only
**4** of them. A CD4-dose arm is a **no-op on 81 of h1-test's 85 rows**. Genome-scale CD4 effects
are not on disk (`DE_stats.suppl_table.csv` carries only per-perturbation DE *counts*), so this is
not fixable inside this experiment.

**INBOX-1 said it is worse than that, and it is right — verified here rather than accepted**
(`vacuity207.py` → `work/vacuity207.json`). The sizing rule takes a **median** row norm over
*every* row, and with 81 of 85 rows at zero:

```
h1-test: 85 targets, CD4 carries all three conditions on 4, nonzero row norms 4
         MEDIAN row norm of the residual = 0.0    (mean 0.2539, max 6.7707)
         -> f_per_unit = 0.0 -> solved fraction = 0.0
```

**The residual is identically zero on h1-test, not merely small.** 200's finding reproduces
exactly.

**One correction to hand back, and it changes the recommended fix.** The silent-zero belongs to
the *harness helper*, not the builder. Production's `solve_cd4_fraction`
(`scripts/make_submission_pm_cd4.py:1258–1266`) **raises** on this input —
`ValueError("VCC_CD4_TARGET_F was set but the CD4 residual has zero norm; nothing to solve")` —
so a `VCC_CD4_TARGET_F` build on a panel like this fails loudly and needs **no new guard**. It is
`plolo.solve_fraction` (`experiments/162-retrodicting-validator/plolo.py:368`) that returns
`0.0` silently, and the legacy `VCC_CD4_FRACTION` path that returns the constant without ever
consulting coverage. **The guard belongs in the harness helper and on the legacy path, not in the
builder** — recorded here and handed back per INBOX-1's instruction; no production file was
touched. (`scripts/def_postbuild.py:119–122` already asserts `cd4 complete targets > 0` on a new
panel, which is the D/E/F-facing half of the same check.)

**And the non-vacuity assertion INBOX-1 asks for, run on THIS experiment's own instrument:**

| | `k562_gwps` | `orion_hct116` | `orion_hek293t` | `kolf2` |
|---|---:|---:|---:|---:|
| targets | 269 | 290 | 294 | 273 |
| rows CD4 carries (all three conditions) | **251** | **267** | **271** | **253** |
| median residual row norm | **4.4677** | **4.4934** | **4.4718** | **4.4491** |

All four are non-zero and comparable to the whitened base's 4.7245, i.e. the sizing rule is in its
designed regime — against h1-test's exact 0.0. Displacement of the delivered effect from the
`f = 0` artifact, 4-fold mean Frobenius relative norm:

| `f` | 0.05 | 0.10 | 0.20 | **0.34** | 0.40 | 0.80 |
|---|---:|---:|---:|---:|---:|---:|
| solved `fraction` | 0.05283 | 0.10566 | 0.21132 | 0.35925 | 0.42265 | 0.84530 |
| `‖delivered(f) − delivered(0)‖ / ‖delivered(0)‖` | **0.0514** | **0.1023** | **0.2017** | **0.3312** | **0.3828** | **0.6624** |

**Every row in this document moved the object by a stated, large amount.** None is a silent zero.

`nmae` is therefore priced by **188's exact identity on the four PLOLO folds, with the pool
construction 188 VALIDATED on h1-test** (`|lfc|·√cpm` at the reference's own BH depth, 89
genes/row — the construction that reproduces h1-test's BH set at 0.7219 against `|lfc|`'s 0.3796,
and whose size matches audit 181's board-implied median `n_real` of 90–110). Both alternative
pools are reported, because 188 showed the construction can invert the sign.

`six207.py` → `work/six207.json`, four folds, challenge frame A, production key `|H|`:

| arm | α | **delivered / a15** | `nmae` (wald89) | Δ | `e` | **overall `nmae`** | top-89 precision | Δ | `nmae` lfc89 | `nmae` h89 |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| **a15 SHIPPED** | 0.982 | 1.0000 | 0.98476 | — | — | — | 0.58995 | — | 0.99219 | 0.97873 |
| res `f = 0.14` | 0.982 | 1.0024 | 0.98441 | −0.00035 | 0.902 | **+0.000097** | 0.59233 | +0.00238 | 0.99205 | 0.97830 |
| res `f = 0.20` | 0.982 | 1.0028 | 0.98407 | −0.00069 | 0.717 | **+0.000191** | 0.59341 | +0.00346 | 0.99186 | 0.97793 |
| res `f = 0.24` | 0.982 | 1.0028 | 0.98478 | +0.00002 | 0.497 | −0.000006 | 0.59570 | +0.00575 | 0.99182 | 0.97744 |
| res `f = 0.28` | 0.982 | 1.0004 | 0.98448 | −0.00028 | 0.549 | +0.000078 | 0.59681 | +0.00687 | 0.99182 | 0.97683 |
| **res `f = 0.34` ← CHOSEN** | 0.982 | **1.0009** | 0.98530 | +0.00055 | 0.435 | **−0.000152** | **0.59752** | **+0.00757** | 0.99239 | 0.97838 |
| res `f = 0.40` | 0.982 | 1.0034 | 0.98722 | +0.00246 | 0.306 | −0.000683 | 0.59778 | +0.00783 | 0.99322 | 0.98472 |
| lib `w = 0.5` @α\* | 1.143 | **0.5228** | **1.08837** | +0.10361 | −1.195 | **−0.028782** | 0.60159 | +0.01164 | 1.04450 | 1.04175 |

**All three pool constructions agree in sign for the residual arms** (wald89, lfc89 and h89 all
move < 0.001 raw up to `f = 0.34`), so 188's sign-inversion hazard does not bite here.

`e` is 169 §5's marginal refund rate: `nmae` improves iff `e > ½`. It sits **at or just either side
of ½ through `f` = 0.24–0.34** (0.497 / 0.549 / 0.435) and only falls decisively below it at
`f = 0.40` (0.306) — which is why `nmae` is flat to `f ≈ 0.34` and starts costing beyond it. **At
the chosen dose `nmae` is a rounding error (−0.000152 overall) and it is not the binding
constraint.**

**The ordering half of `reach`/`fid` improves monotonically**, +0.00757 in |truth|-weighted top-89
sign precision at `f = 0.34`. That is favourable and unpriced (its yield halves carry `n_real` and
are unreadable on any harness we hold — rule 9).

### 6.2 The α re-solve the mandate demanded — and it costs exactly nothing

180 §"`VCC_CD4_FRACTION` vs `VCC_CD4_TARGET_F`" established that **with a solved target `f` the
whole pipeline is homogeneous in α**, so a dose change is not an amplitude change. Verified here
three ways rather than asserted:

1. `add_post_whiten_residual`'s preserve scalar pins the effect-space median row norm: **4.7245 at
   every one of the 25 doses**, to four decimals.
2. **Delivered** median row norm (post-renderer, the thing the board reads): **1.0009×** a15's at
   `f = 0.34`, 1.0034× at `f = 0.40`.
3. Re-solving α to match a15's pre-whitening norm returns **α\* = 0.982**, i.e. the identity —
   `res f=0.34 @a*` and `res f=0.34` are the **same row to every digit** in the table above.

**Cost of the α re-solve: 0.000000.** (203's re-solve cost −0.004645 `pds` because its arm shrank
the estimator 4.7×; nothing here shrinks.)

`amp207.py` closes 170's trap explicitly — the preserve pins the **median** while `A` is a **mean
of squares**, and under a γ sweep those came apart by 14%:

| `f` | median delivered `|row|` | ratio | `A = mean |row|²` | ratio | `eps = A/B − 1` |
|---:|---:|---:|---:|---:|---:|
| 0.00 | 3.871209 | 0.99870 | 18.110561 | 1.00167 | 0.005492 |
| **0.10 (base)** | 3.876234 | 1.00000 | 18.080327 | 1.00000 | 0.005409 |
| 0.28 | 3.877258 | 1.00026 | 17.953178 | 0.99297 | 0.005029 |
| **0.34** | 3.875740 | **0.99987** | 17.906028 | **0.99036** | 0.004857 |
| 0.40 | 3.882156 | 1.00153 | 17.893587 | 0.98967 | 0.004674 |
| 0.80 | 3.923472 | 1.01219 | 18.221390 | 1.00780 | 0.003436 |

**`A` moves by at most 1.0% over the entire dose range** and by 0.96% at the chosen point, against
γ's 14%. `eps` moves 0.00541 → 0.00486, which enters raw `mse` through `q·eps` as 0.0006 raw =
0.0001 overall. **Neither is credited in the price**: the `A` drift's sign is favourable
(a 1% smaller `A` lowers raw `mse`), and not crediting it is the conservative choice, since 187
closed the amplitude channel at α = 0.982 and any amplitude credit belongs to that closed lever.

### 6.3 The six members at `f = 0.34`, all of them (rule 12)

| member | value | route |
|---|---:|---|
| `pds` | **−0.000685** | PLOLO 4 folds, Δ = −0.001837 raw = **1.23σ** on 162's LOO σ (0.001494); paired-target 95% CI [−0.00746, +0.00373] **contains zero**; /0.452/6 × 1.011 |
| `mse` | **+0.005206** (k=1) / **+0.006558** (k=1.2598) | 096/197 identity through `X` only. `claim`, `Kc`, `ΣD` fixed (emission law, spread, `keep_top`, `min_cpm` all untouched); `A` fixed to 1.0% (§6.2, measured) |
| `nmae` | **−0.000152** | 188's exact identity, wald89 significance pool, 4 folds; all three pool constructions agree in sign |
| `reach` | **unpriced; ordering half +0.00757 (favourable)** | rule 9 — `n_real` is 222 on h1-test, 68–92 on the repaired harness, ~90–110 implied on the board. **No emission-law toll**: this arm does not touch the emission law (169 falsified the toll for profile-class arms at +0.0006; audit 181 confirms the −0.0146 is an emission-law toll) |
| `fid` | **unpriced; same ordering half +0.00757** | as `reach`; the yield half is `min(1, n_pred/n_real)` and `n_pred` is set by `keep_top`, which does not move |
| `jac` | **unpriced, and NOT exempt** | the call set is `keep_top 8000` applied to the delivered profile, and the profile moves. Rule 12's exemption is for `spread`/`phi`/`keep_top` arms only. This is the single largest uncovered risk and it is exactly what the +0.0019 shape-only threshold (realised RMS 0.001877 of the unpriced channels) exists to cover |

**Total priced: +0.00437 (conservative) to +0.00572 (board-calibrated).**

### 6.4 The gate

**Class: shape-only, threshold +0.0019.** 162 Gate 2's criterion is "does the arm change the
effect amplitude (α, γ, λ, any renormalisation)". α, γ and λ are unchanged; the renormalisation
that exists (the preserve scalar) is present at every dose and holds the delivered median row norm
at **1.0009×** and `A` at **0.990×**. §6.2 measures all three rather than asserting them. The arm
is also NOT priced through PLOLO's `mse` channel — which Gate 2 declares invalid — but through
096/197's board identity, the route rule 12 and 204 use.

| | conservative (k = 1.0) | board-calibrated (k = 1.2598) |
|---|---:|---:|
| in-sample priced total | **+0.00437** | +0.00572 |
| **out-of-sample (leave-one-fold-out, §7.3)** | **+0.00340** | +0.00502 |
| gate (shape-only) | +0.0019 | +0.0019 |
| **margin** | **2.3× / 1.8× out-of-sample** | 3.0× / 2.6× |
| gate if forced into the amplitude class | +0.0089 | +0.0089 |

Under the **amplitude** class the arm would fail. §6.2 is the reason it is not in that class, and
every number in it is measured. If a reviewer rejects that classification, the correct reading is
"veto", and I would rather that call be made explicitly on §6.2's three measurements than
implicitly.

---

## 7. Phase 4 — the four falsifiers

### 7.1 Falsifier 1 — target-shuffled CD4. **PASSES at every dose.**

Whole rows (effect + gene support + carrying count) permuted across targets by a derangement,
seed 20726. Geometry, per-row norms, support pattern and number of carrying rows are all
identical; only the target→row assignment is destroyed.

| `f` | real CD4 Δcos_agg | **shuffled** Δcos_agg | real Δ`pds` | **shuffled** Δ`pds` |
|---:|---:|---:|---:|---:|
| 0.05 | +0.001030 | **−0.000049** | +0.004472 | **−0.001992** |
| 0.10 | +0.001995 | **−0.000149** | +0.007005 | **−0.004641** |
| 0.20 | +0.003711 | **−0.000495** | +0.008275 | **−0.012463** |
| 0.40 | +0.006239 | **−0.001643** | +0.003375 | **−0.034064** |
| 0.80 | +0.008325 | **−0.004660** | −0.011443 | **−0.071111** |

**The shuffled donor buys nothing at any dose and loses monotonically on both endpoints.** The
gain is target-specific biology, not amplitude and not geometry. The same control run on the
*library* form (`libnull207.py`) gives −25.0% cosine and −0.022 `pds` at `w = 0.5`.

INBOX-1 §3 named the outcome to watch for: **204's repression-weighting arm died because its
permutation control TIED its best real arm.** That is not what happens here — the permutation is
strictly and increasingly *negative* at every one of twelve doses, and the gap between real and
shuffled widens monotonically (+0.00108 → +0.01299 in Δcos_agg from `f = 0.05` to `f = 0.80`).
Against 204's measured instrument resolution of **±0.20% relative on Δcos**, the chosen arm's
+18.0% is **90× the noise floor**, which is the regime INBOX-1 §2 predicted would settle the
hypothesis either way.

### 7.2 Falsifier 2 — lineage specificity. **The mandate's test cannot be run; the substitutes say the gain is lineage-GENERIC.**

**What cannot be measured, and why.** The test as specified ("larger on context A than on B or C")
needs challenge truth, which does not exist. The obvious substitute — a T-lineage *recipient* —
does not exist either: `data/donors/nadig_jurkat` (2,270 perturbations) covers **0 of the 300
panel targets**, because the panel excludes common essentials by construction, and the CD4 banks
carry **only** panel targets. Verified directly: panel ∩ nadig_jurkat = **0**; panel ∩ rpe1 = 0;
panel ∩ k562_essential = 0. **No surface on disk pairs a T-lineage recipient with the CD4 donor.**
Saying so is the result; a "lineage" test run on a surface with no T-lineage recipient would be
theatre.

**Substitute A — per-recipient homogeneity.** All four held lines are **non-T**:

| held recipient | lineage | relative cos_agg gain at `f = 0.10` | at `f = 0.34` |
|---|---|---:|---:|
| `k562_gwps` | CML / erythroleukemia (myeloid) | **+12.3%** | +20.2% |
| `orion_hct116` | colorectal carcinoma (epithelial) | **+11.2%** | +18.2% |
| `orion_hek293t` | embryonic kidney (epithelial) | **+7.6%** | +11.8% |
| `kolf2` | iPSC (pluripotent) | **+15.2%** | +25.3% |

**Four unrelated non-T lineages, four positive gains of the same order.** A lineage-matching
mechanism predicts ≈0 on all four. Per **render frame** the gain is likewise uniform (A +11.14%,
B +10.80%, C +11.14%) — as it must be, since production's effect matrix is context-independent
and the frames differ only through `b_ctrl` and the CPM gate. That measurement is reported for
completeness and is **not** a lineage test.

**Substitute B — how close is CD4 itself to each lineage?** (`lineage207.py`)

| held line | lineage | cos(CD4, truth) | cos(our consensus, truth) | CD4 / consensus |
|---|---|---:|---:|---:|
| `k562_gwps` | myeloid | **0.033654** | 0.032181 | **1.046** |
| `orion_hct116` | epithelial | 0.024061 | 0.033310 | 0.722 |
| `orion_hek293t` | epithelial | 0.015311 | 0.035416 | 0.432 |
| `kolf2` | pluripotent | 0.010747 | 0.008062 | 1.333 |

There **is** a weak gradient: CD4's raw alignment is highest to the one haematopoietic recipient
and lowest to kidney and iPSC. But it is four points, confounded with truth quality (194 measured
achieved knockdown 0.845 / 0.700 / 0.490 for K562 / HCT116 / HEK293T, exactly the ordering of the
consensus column), and **CD4 alone already matches or beats our entire three-donor consensus on
two of the four lines.**

**Answer to the mandate's mechanism question (206 §5's "a lineage-matched donor escapes 159's
1/k law").** My measurements do **not** support that mechanism as the source of what we can
collect. What they support instead is the plainer reading: **CD4 is an unusually good estimator of
the *context-independent shared* component** — 22M cells, 4 donors, 3 conditions — so it acts as an
extra term **on** 159's curve rather than off it. Translated through that law
(`1/cos²_surface = 66.5 + 1920/k`, board→surface at 197's a15 anchor ÷2.078):

| configuration | implied board cosine | **`k_eff`** |
|---|---:|---:|
| CD4 **off** | 0.07138 | **2.46** |
| shipped `f = 0.10` | 0.07925 | **3.09** ← reproduces 197 §4.2's 3.09 exactly |
| best equal consensus of **all four** public panel-covering lines (197 §4.2) | 0.08730 | 3.84 |
| **`f = 0.34` (this arm, conservative)** | **0.09350** | **4.49** |
| `f = 0.34` (board-calibrated) | 0.09720 | 4.92 |
| the four-account cluster's floor (206 §2) | 0.13117 | 10.41 |

**The shipped residual is worth +0.63 effective lines; the dose increase is worth another +1.40 —
more than the residual's original introduction, and it takes us past the best equal consensus of
every panel-covering line that exists (166's census: there is no fifth).** That is what a 22M-cell,
4-donor, 3-condition atlas should be worth as a term in a `1/k` law, and it needs no lineage story.
It also does not reach the cluster: `k_eff` 4.5 against their 10.4.

**Consequence for D/E/F, which is the decision this falsifier exists to make: the collectable
component is lineage-independent, so it should carry to D/E/F whatever their lineages turn out to
be on Oct 22.** If there *is* an additional Jurkat-specific bonus on context A, it is upside we
cannot measure and have not priced.

### 7.3 Falsifier 3 — the dose must not be selected on the surface it is scored on. **Done, and it costs 29%.**

`select207.py`. For each held fold the operating point is chosen by maximising the composite on
the **other three folds only**, then priced on the held fold. Selection frame and reporting frame
never overlap.

| held fold | chosen on the other three | in-sample | **held-fold price** |
|---|---|---:|---:|
| `k562_gwps` | `f = 0.34` | +0.00310 | **+0.00837** |
| `orion_hct116` | `f = 0.47` | +0.00509 | **+0.00356** |
| `orion_hek293t` | `f = 0.55` | +0.00836 | **−0.00732** |
| `kolf2` | `f = 0.34` | +0.00352 | **+0.00897** |
| **mean of held-fold prices** | | | **+0.00340** (calibrated +0.00502) |
| in-sample argmax over all four | `f = 0.40` | +0.00477 | — |

**The honest out-of-sample number is +0.00340, 71% of the in-sample maximum.** The family choice
is stable (residual on 4/4) and the dose choice lands in [0.34, 0.55] on 4/4 with **`f = 0.34` the
mode**. The same procedure on the other three forms returns **−0.00806** (library), **−0.00195**
(composed) and **−0.00560** (whitened residual).

`orion_hek293t` is the single dissenting fold. It is also the donor with the **worst measured
knockdown** (194 / 204 Phase 0: self-gene residual 0.510 against K562's 0.155), i.e. the most
attenuated truth of the four. That is a reason to expect it to be the noisiest fold, not a reason
to discount it — it is included at full weight in every number above.

### 7.4 Falsifier 4 — reconcile with 062's LODO, and show its `pds` curve is genuinely not flat.

Both halves. §3.1 reproduces 062's three headline numbers at 0.99× / 1.14× / 1.18×. And the curve
is **not** flat where 062 stopped — the two lines it selected on move in opposite directions:

| `fraction` (062's frame) | 062 HCT116 | **here** | 062 HEK293T | **here** | here: K562 | here: KOLF2 |
|---:|---:|---:|---:|---:|---:|---:|
| 0.01 | +0.000943 | +0.001376 | +0.000669 | +0.000820 | +0.001036 | +0.000768 |
| 0.02 | +0.002212 | +0.002585 | +0.000960 | +0.001645 | +0.002136 | +0.001832 |
| **0.05 ← 062's pick** | **+0.004818** | **+0.006022** | **+0.002418** | **+0.003065** | +0.004776 | +0.004027 |
| 0.10 | +0.006944 | +0.009768 | +0.001989 | +0.003533 | +0.008036 | +0.006684 |
| 0.20 | +0.008607 | +0.012429 | **−0.002418** | **−0.001676** | +0.012794 | +0.009553 |
| 0.42 (beyond 062's grid) | — | +0.004594 | — | −0.017847 | +0.014990 | +0.011761 |

**Same shape, same signs, same turning points, magnitudes 1.0–1.4× larger on the rendered
instrument.** 062's `pds` curve at `fraction = 0.05` is steeply *rising* on HCT116
(+0.0048 → +0.0086 by 0.20) and near its peak on HEK293T. **062's conservative worst-fold rule was
driven entirely by HEK293T**, the one fold with the worst knockdown, and it stopped the dose at
roughly a quarter of where three of four lines still wanted it. Nothing about 062 was wrong; its
selection rule and its objective were both conservative, and this experiment is the first time
either has been re-examined against the channel the residual pays in.

---

## 8. Verdict, operating point, and the preregistration

### 8.1 Chosen operating point: `VCC_CD4_TARGET_F = 0.34`

(3.4× the shipped `f = 0.10`. The builder will solve this to `fraction ≈ 0.322` on the challenge
panel and it solves to 0.359 on PLOLO — see §4's footnote 1; **set `f`, never `fraction`.**)

Three criteria, stated before the arm was chosen, all select it:

1. **argmax of the measured six-member price** (`pds` + `mse` + measured `nmae`): +0.00437 at
   `f = 0.34`, against +0.00413 at 0.28 and +0.00408 at 0.40;
2. **modal leave-one-fold-out choice** (2 of 4 held folds; the other two chose 0.47 and 0.55);
3. **the largest dose at which no priceable member is resolvably worse** — `|Δpds| = 1.23σ`
   against 162's Gate 1 threshold of 1.645σ; `f = 0.40` fails this at 2.43σ.

`f = 0.28` is the conservative fallback (+0.00413, `pds` at 0.16σ). `f = 0.20` is the
"nothing-goes-backwards" point (+0.00314, every priceable member positive).

### 8.2 The exact builder invocation

**One variable against the shipped `rc_hybrid`.** Copy `experiments/191-rc-hybrid-build/build191.sh`
and change `VCC_CD4_TARGET_F` from `0.1` to `0.34`. Nothing else moves.

```bash
export VCC_ALPHA=0.982
export VCC_CD4_TARGET_F=0.34          # <-- the ONLY change; rc_hybrid shipped 0.1
uv run python scripts/make_submission_pm_cd4.py \
  --whitener          data/donors/whitener_union_pruned.npz \
  --whitener-metadata data/donors/whitener_union_pruned.json \
  --whitener-id       union_pruned \
  --min-cpm 10 --tail-spread 1.0 --keep-top 8000 --keep-spread 0.25 \
  --overdispersion 4.0 --max-rss-gib 24 \
  --emitter rc_hybrid \
  -o submissions/207-cd4f034.h5ad
```

**Pre-build assertions (rule 10 — diff the finished meta, not the command typed):**

| field in the built `*.h5ad.meta.json` | required value |
|---|---|
| `effects.cd4.residual_frame.target_f_requested` | **0.34** |
| `effects.cd4.residual_frame.fraction_source` | `VCC_CD4_TARGET_F` |
| `effects.cd4.residual_frame.solved_fraction` | ≈ **0.32** (challenge frame; `= 0.34 / f_per_unit`, and `f_per_unit` was **1.0556645369046553** on the abc-knee build) |
| `effects.cd4.fully_supported_targets` | **276** (unchanged) |
| `effects.cd4.residual_frame.effective_f_at_cd4_fraction` | **0.05278322684523277** (unchanged — it is a diagnostic of `VCC_CD4_FRACTION`, not of the arm) |
| `operator.sha256` | **`e2f65ecf69b7e0b0b817d02d850c7be07a78b442644a239c4b0e520cda1db123`** |
| `alpha` | **0.982** — do **not** re-solve; §6.2 measures α\* = 0.982 exactly |
| `emitter` | `rc_hybrid` |
| **total meta fields differing from `191-rc_hybrid.h5ad.meta.json`** | **exactly 2** (`cd4_target_f`, and `effects.cd4.fraction` / `solved_fraction` which is its deterministic consequence) — plus `effects.whitening.after_all_median_norm` unchanged |

Model name: a random joke, description **empty** (§14).

### 8.3 Preregistered per-member board prediction

Against our current displayed entry (`rc_hybrid`, overall **0.18342**, raw `pds` 0.8360,
raw `mse` 0.8738045, raw `nmae` 0.9442, `fid` 0.5037, `reach` 0.2299, `jac` 0.0302).
**Written before the build exists.** The conservative column (`k = 1.0`) is the prediction; the
calibrated column is the alternative if 197's measured PLOLO→board factor holds.

| member | mechanism | predicted Δ scaled | predicted raw |
|---|---|---:|---:|
| **`mse`** | 096/197 identity, `X` only | **+0.03123** (cal +0.03933) | **0.84413** (cal 0.83643) |
| **`pds`** | PLOLO ×1.011 | **−0.00411** | 0.83414 |
| **`nmae`** | 188 identity, wald89 pool | **−0.00092** | 0.94475 |
| `reach` | unpriced; ordering half +0.00757 | 0 ± 0.010, sign more likely **positive** | — |
| `fid` | unpriced; same ordering half | 0 ± 0.010, sign more likely **positive** | — |
| `jac` | unpriced; call set moves with the profile | 0 ± 0.010 | — |
| **overall** | | **+0.00437** (cal +0.00572) | **0.18779** (cal 0.18914) |
| **transfer cosine** (197's inversion will read this back) | | | **0.09350** (cal 0.09720) |

**Falsifiers on the board result, stated in advance.**

1. If raw `mse` lands **above 0.8738** the dose-response reversed and the whole sweep is wrong.
2. If raw `mse` improves but raw `pds` falls by **more than 0.005**, PLOLO's `pds` channel
   understated the cost and the dose should be walked back to `f = 0.20`.
3. If the inverted cosine lands **below 0.0865** (half the predicted move), `k_cos < 0.5` and every
   `mse` price in this document is over by 2×.
4. If `jac` or `fid` moves by more than 0.010 scaled, the shape-only class was wrong and rule 13's
   +0.0089 threshold should be applied to the next arm in this family.

**Sanity ceiling:** at this dose the implied board cosine is 0.0935, still **0.72×** the cluster's
floor of 0.13117 (206 §2). This arm closes roughly **23%** of the direction gap to the cluster and
**12%** of the gap to first place. It does not finish the job.

### 8.4 What is still on the table in this territory (not taken)

- **`f = 0.40–0.80` is worth more on cosine** (+21% to +32% over shipped) and is only vetoed by
  `pds`. If a future measurement establishes that the board's `pds` is less dose-sensitive than
  PLOLO's — and 162's own a15 row (`k₋ₐ = 0.680`, PLOLO *under*-predicting the realised) hints
  that it might be — the optimum moves up. **One board point at `f = 0.34` resolves this**, since
  it re-measures the PLOLO→board `pds` factor for this family at 3.4× the previous lever arm.
- **`f` is the only knob I swept.** The CD4 *source* is fixed at 062's equal mean of
  Rest/Stim8hr/Stim48hr with `gene_norm`. Per-condition weights, per-target dosing, and re-extracting
  the CD4 atlas at genome scale (which would also unlock h1-test for this arm, §6.1) are untouched.
- **Not my territory and deliberately untouched:** the whitening operator (193), the donor tables
  (194), the emission law (185/191/202). §5.2's finding — that the cosine lives in what is added
  *after* the operator — is a fact about the operator that its owner may want.

---

## 9. Appendix — the operator census. **There are 27, not nine, and SIX are called `whitener_union.npz`.**

The mandate says "nine near-identically-named operators exist"; `AGENTS.md` §11.1 rule 10 names
three. A repo-wide `rglob` for every `.npz` carrying `evals` + `evecs` with a square basis > 1,000
finds **27**. Recorded here once so no future screen has to re-derive it.

| path | sha256 (8) | basis | max eigenvalue |
|---|---|---:|---:|
| **`data/donors/whitener_union_pruned.npz` — PRODUCTION** | **`e2f65ecf`** | 17,942² | **5.8946** |
| `data/donors/whitener_union.npz` | `5608e807` | 17,942² | 6.2567 |
| `data/donors/whitener_union_full.npz` | `965d3cef` | 17,942² | 6.9802 |
| `data/donors/whitener_full.npz` | `dc6bd2ce` | 7,312² | 4.4586 |
| `experiments/040-union-gamma/out/whitener_union.npz` | `f950b268` | 17,530² | 6.8999 |
| `experiments/062-…/whiteners/held_k562_gwps.npz` | `13d988fc` | 17,806² | 8.4259 |
| `experiments/062-…/whiteners/held_orion_hct116.npz` | `7cecbf7f` | 17,942² | 8.0363 |
| `experiments/062-…/whiteners/held_orion_hek293t.npz` | `97e486e0` | 17,942² | 3.6958 |
| `experiments/095-clean-harness/out/whitener_union.npz` | `5cf3899d` | 17,942² | 8.0364 |
| `experiments/095-clean-harness/out-essential/whitener_union.npz` | `1fbd44b0` | 17,942² | 7.7340 |
| `experiments/095-clean-harness/out-h1test/whitener_union.npz` | `f80fe8d6` | 17,745² | 8.0808 |
| `experiments/095-clean-harness/out-random/whitener_union.npz` | `05a17699` | 17,942² | 8.3616 |
| `experiments/096-mse-floor/out-alpha/whitener_union.npz` | **`5cf3899d`** (duplicate of 095's) | 17,942² | 8.0364 |
| `experiments/099-donors/work/whit/whitener_hek293t.npz` | `4502bb4f` | 17,942² | 3.2940 |
| `experiments/137-jurkat-tangent-cone/work/source-candidate-v1/operator.npz` | `52374bcb` | 17,585² | 5.7517 |
| `experiments/199-def-runbook/work/whitener_union_pruned_abc_ka.npz` | **`e2f65ecf`** (byte-identical to production — 199's known-answer copy) | 17,942² | 5.8946 |
| `experiments/199-def-runbook/work/whitener_union_pruned_def199.npz` | `68dd20dc` | 17,942² | 5.9389 |
| 10 further small/exotic operators (129, 140 ×2, 144 ×4, 195 ×3) | — | 4,992²–7,760² | 1.01–1,746.5 |

**Three traps this census makes visible, none of which the eigenvalue or the filename catches:**

1. **Six distinct files are named `whitener_union.npz`** (`5608e807`, `f950b268`, `1fbd44b0`,
   `f80fe8d6`, `05a17699`, `5cf3899d`), across five directories, with bases 17,530² / 17,745² /
   17,942². Rule 10's near-miss (b) named three of them.
2. **`experiments/095-…/out/whitener_union.npz` (`5cf3899d`) and
   `experiments/062-…/held_orion_hct116.npz` (`7cecbf7f`) are DIFFERENT FILES whose entire spectra
   agree to 1.3e−5** (max eigenvalues 8.036352796798983 vs 8.036339907795812) — two independent
   fits of the same held-HCT116 basis, one of them misleadingly named `union`. **The eigenvalue is
   not a sufficient identifier either.**
3. **`experiments/096-mse-floor/out-alpha/whitener_union.npz` is byte-identical to 095's.** Two
   experiments' operator references resolve to one file; a change to either would move both.

Use the sha, always.

## 10. Files

| script | what it produces |
|---|---|
| `sweep207.py` | the 25-dose sweep + shuffled / self / KOLF2 residual controls → `work/sweep207.json` (232 rows) |
| `analyse207.py` | known-answer gates, dose table, the three maximisers, per-fold homogeneity → `work/analyse207.json` |
| `price207.py` | the full-grid rule-13 composite and the `k_cos` calibration → `work/price207.json` |
| `library207.py` | CD4 as a fourth library member, 15 weights + 42 composed settings × 2 α → `work/library207.json` (456 rows) |
| `libnull207.py` | shuffled-CD4 and KOLF2 **library** falsifiers → `work/libnull207.json` |
| `form207.py` | the whitened-residual decomposition → `work/form207.json` |
| `libprice207.py` | the form comparison, priced → `work/libprice207.json` |
| `six207.py` | 188's exact `nmae` identity, ordering precision, delivered amplitude → `work/six207.json` |
| `amp207.py` | 170's median-vs-mean-of-squares trap, closed → `work/amp207.json` |
| `select207.py` | leave-one-fold-out selection (falsifier 3) → `work/select207.json` |
| `lineage207.py` | CD4's alignment to each recipient lineage (falsifier 2) → `work/lineage207.json` |
| `vacuity207.py` | INBOX-1 verified: h1-test's exact zero, and this instrument's non-vacuity at every dose → `work/vacuity207.json` |

Peak RSS **4.61 GB**. No `vcc`, no `cell-eval2`, no build, no submission, nothing written outside
this directory.
