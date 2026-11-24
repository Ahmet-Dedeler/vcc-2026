# 177 — the resampling emitter: **the ordering defect is CLOSED** (98–110% of INBOX-4's counterfactual realised), the noise schedule is **refuted**, and the arm is **NO-GO on rule 13** because it pays in `n_pred`

Date: **2026-09-02**. Measurement only. No `cell-eval2`, no build, no submission, no slot spent.
Surfaces: `h1-test` (realisation, `reach` board-validated for emission-law changes) and the real
challenge panel (`mse` claim, zero-effect null). Donors `k562_gwps` + `orion_hek293t`; `h1_2025`
never enters. Peak RSS 5.8 GB.

## 0. The five things that changed

1. **The coordinator's first pass was a repair artefact, and it is now diagnosed exactly.**
   It rescaled every real control cell from ~20k to `EMIT_L = 100_000` while the pinned panel
   total stayed at the artifact's own ~20k/cell depth, so `emit_counts`'s exact-total repair had
   to peel **~80% of every gene's column, one count per cell per round** — a maximally
   variance-*destroying* operation. It called 10,776/10,779 genes significant.
   `signal177.py` reproduces that geometry and reads **10,658/10,779**. Fixed by passing raw
   counts and the cells' own libraries.
2. **`EMIT_L` was never the emitted depth**, so the mandated `quasi20` control is a no-op.
   `profile` is scaled by `S/(N_EMIT·EMIT_L)`, hence `totals = floor(N·L·profile) = floor(S·p̃)`
   for *any* `L`. The artifact's `S/400` is **20,250 counts/cell** against a real h1-test control
   median of **20,190** — **the shipped quasi arms were already at real depth.** A `quasi20` row
   is bit-identical to `base_mf` by construction (an identity, not a control), so the informative
   depth arm is the other direction: `quasi100`, a genuine 5×, which **loses** 0.0108 reach.
3. **The zero-effect null holds, and it is the largest single number in this experiment.**
   On all three challenge contexts, at the production depth, through the unmodified production
   `emit_counts`: **1 / 18 / 19 genes BH-significant at zero predicted effect** (0.0–0.2%),
   against **4,302 / 4,312 / 3,883 (43.3% / 44.8% / 38.4%)** for a15's shipped law.
4. **`rs_all` reverses INBOX-4's falsifier with no schedule at all.** Spearman(−log₁₀ p_adj, ·)
   inside our own call set moves from a15's **+0.485 expression / −0.391 |H|** to **+0.323 /
   +0.480**, and `hy_kt2000` reaches **+0.001 / +0.601**. Raw `reach` moves **+0.041 to +0.056**,
   which is **81–110% of INBOX-4's full `|H|` ceiling**.
5. **The mandated noise schedule is refuted, by a new mechanism.** Adding symmetric noise to a
   resampled sample does **not** suppress `|z|` — it *inflates* it, because the rank-sum's null
   variance is computed under exchangeability and over-dispersing our own sample breaks it
   upward. The permutation control **beats** the real schedule, exactly as in 171 §12b.

---

## 1. Instruments and known-answer gates — eleven, all passing

| gate | source | expected | measured |
| --- | --- | ---: | ---: |
| `KA_a15_raw` median `n_pred` | 157/164 | 6035 | **6034** ✓ |
| `KA_a15_raw` `reach` f=1 | 157/164 | 0.18283 | **0.18283** ✓ |
| `KA_a15_raw` `reach` key=\|H\| | INBOX-4 §2 | 0.23533 | **0.23533** ✓ |
| `KA_a15_raw` `reach` sig-then-\|H\| | INBOX-4 §2 | 0.21914 | **0.21914** ✓ |
| `base_mf` median `n_pred` | 173 | 6066 | **6066** ✓ |
| `base_mf` `reach` f=1 | 173 | 0.18432 | **0.18432** ✓ |
| recovery round trip | 173 | ≤ 3× seed | 0.00149 vs 0.01628 ✓ |
| `claim` flat 0.25/φ0 | 096/164 | 0.11860 | **0.11860** ✓ |
| `claim` flat 1.00/φ0 | 096/164 | 0.19216 | 0.19438 (+1.2%) ✓ |
| `claim` flat 0.25/φ4 (a8) | 096/164 | 0.16478 | 0.16652 (+1.1%) ✓ |
| `claim` a15 as shipped | 164 §5 | 0.24581 | **0.24581** ✓ |

Plus the zero-effect null's own gates against 157's closed form on challenge ctx A/B/C
(`work/null177.*.log`): quasi 0.25/φ4 **19.3 / 21.8 / 14.8 %** against 157's 19.0 / 22.8 / 16.2;
quasi 1.00/φ4 **77.3 / 72.1 / 74.2 %** against 78.9 / 73.9 / 76.2, all 0.0–0.2% UP; quasi
0.25/φ0 **44.8 / 49.7 / 36.3 %** against 44.9 / 48.5 / 36.3. The one row that does not reproduce
is quasi 1.00/φ0 (4.3% here vs 157's 7.7%), because `emit_counts` repairs the scattered part's
total exactly and 157's closed form does not; that row is not shipped anywhere.

**Rule-12 pseudobulk check, in 171 §6's non-vacuous form.** One sha256 of the *delivered*
per-gene panel total per (arm, target). **All twelve real-depth arms are identical on 85/85
targets**, digest `19842254b81a7f9b`; `KA_a15_raw` reads `b53d02aebfd0816a`, which is the digest
171 §6 recorded for a15. Without the exact-total repair the residual lands in the totals and
every arm differs, so the check reads differently when the thing it tests is absent. **So
Δ`pds` = 0.000000 is an identity, not an estimate**, and so is the profile half of Δ`mse`.

**One member the pin does NOT hold.** BRIEF §3 defines `lfc` as log₂ of the **arithmetic mean of
per-cell CPM**, not of the pseudobulk. The quasi law gives every cell the same library `L`; the
resample law gives each cell the real library it came with (7.3k–124k). So **`nmae` is genuinely
free to move on these arms and is measured, not asserted** — 157 §5a's "bit-identical" claim
covers `pds` and `mse` but not `nmae` once per-cell libraries change.

---

## 2. The zero-effect null — the property the whole experiment rests on

`null177.py`, production `emit_counts`, challenge controls, 400 emitted cells against the
context's own 18,400, exact mid-rank Mann-Whitney, BH at 0.05.

| emission law | ctx A | ctx B | ctx C | direction |
| --- | ---: | ---: | ---: | --- |
| quasi flat 0.25/φ4 | 28.4% | 33.1% | 18.8% | 84–91% UP |
| quasi flat 1.00/φ4 (a15's tail) | 82.4% | 77.7% | 80.1% | ~100% DOWN |
| **quasi SHIPPED (tail 1.0 / kt 8000 / ks 0.25)** | **43.3%** | **44.8%** | **38.4%** | mixed |
| **`rs_all` — resample, no jitter** | **0.0%** (1) | **0.2%** (18) | **0.2%** (19) | — |
| `rs_all`, `control_profile` = mean-CPM | 0.0% (3) | 0.1% (10) | 0.4% (36) | — |
| `rs` + multiplicative jitter σ 0.60 | 72.2% | 61.0% | 70.1% | ~100% DOWN |
| `rs` + log-normal jitter σ 0.60 | 63.5% | 54.8% | 62.9% | ~100% DOWN |
| `rs` + `fracpair` jitter σ 0.60 | 0.1% | 0.2% | 0.6% | — |
| `rs` + `fracpair` σ 1.00 | 39.2% | 13.9% | 36.6% | ~100% DOWN |
| `rs` + expression-compensating schedule | 25.9% | 19.9% | 27.4% | ~100% DOWN |

The first fact is the headline: **the resampling law is null under the null, to within one gene
in ten thousand, where the shipped law calls 38–45% of the DE universe at zero effect.** Real
control cells carrying a fold change are exchangeable with the reference controls, so `F̄ = ½`
identically rather than asymptotically.

The second fact kills the mandated schedule. `null177b.py` sweeps kernels × σ on the h1-test
controls and shows why: at a **pinned column total**, any multiplicative spread lowers the
median, and any *absolute* symmetric transfer drives the losing half onto the zero floor.
Only `fracpair` — a pair transfer of a uniform fraction of the **losing cell's own count**,
which preserves each gene's total exactly and can never clip — keeps the null, and only to
σ ≲ 0.6. The artefact concentrates in the **highest**-expression tercile (75.6% at absolute-pair
σ 0.8, against 0.0% in the lowest), i.e. exactly on the genes an expression-compensating
schedule must noise.

---

## 3. The signal response — why "add noise to suppress significance" is unavailable on ANY emitter

`signal177.py`, real predicted effects (169's `H`, production `minnorm` mass fix), 6 targets,
against the real h1-test controls.

| `fracpair` strength | median `n_pred` | \|z\| T0 | \|z\| T1 | \|z\| T2 | ρ(p, expr) in-set | ρ(p, \|H\|) in-set |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| **σ = 0 (plain resample)** | **695** | 0.61 | 0.88 | 1.31 | **+0.294** | **+0.462** |
| σ 0.50 ×1 | 654 | 0.60 | 0.92 | 1.35 | +0.307 | +0.385 |
| σ 1.00 ×1 | 2,684 | 0.94 | 1.84 | 1.99 | +0.288 | +0.551 |
| σ 1.00 ×2 | 5,147 | 1.04 | 2.57 | 3.24 | +0.386 | +0.382 |
| σ 1.00 ×4 | 6,780 | 1.07 | 3.22 | 4.96 | +0.604 | +0.126 |
| σ 1.00 ×8 | 7,226 | 1.01 | 3.66 | 6.48 | **+0.758** | **+0.058** |

`|z|` rises **monotonically** in the noise, in every expression tercile, and the expression
dependence the schedule was built to cancel is *doubled*. INBOX-4 §3/§3b's model — `z ≈ μh /
√(σ²_ctrl + σ²_ours)`, so raising `σ_ours` must lower `z` — is wrong for a **rank** statistic:
`z` is computed from `F̄` with the null variance `1/(12N)`, which assumes our sample and the
reference's are exchangeable. Over-dispersing our own sample makes `Var(F̄) > 1/(12N)`, so the
reported `|z|` inflates whatever the mean does.

**This generalises 171 §12c rather than escaping it.** 171 closed the spread-schedule family
because the quasi law's two noise channels both move the *median*. The resample law has no
median artefact — and the family is still closed, because a rank-sum has **no variance channel
at all** in the direction we need. The only way to lower `|z|` at a fixed mean is to make our
cells *tighter* than the controls, which is the even allocator, which is the median artefact.
**There is no setting of any per-gene noise vector, on any emitter, that suppresses unearned
significance.** The lever is the emission LAW, not a schedule on top of it.

---

## 4. The six-member table — h1-test, 85 perturbations, matched-seed twins

`realize177.py` / `six177.py`. `reach` is twin-averaged where a twin exists (the base pair's
same-policy spread is **0.01628**, 171 §7's known problem; the resample arms are far tighter —
`rs_all` 0.00329, `hy_kt2000` **0.00026**).

| arm | `n_pred` | `reach` (twin) | Δ raw | % of \|H\| ceiling | precision | `nmae` | ρ expr | ρ \|H\| |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| `KA_a15_raw` (known answer) | 6034 | 0.18283 | +0.0066 | 13% | 0.50701 | 0.96659 | **+0.485** | **−0.391** |
| **`base_mf` (incumbent, twin mean)** | **6066** | **0.17619** | — | 0% | 0.50789 | 0.96639 | +0.484 | −0.395 |
| `quasi100` (genuine 5× depth) | 8774 | 0.16542 | −0.0108 | −21% | 0.50561 | 0.95847 | −0.506 | −0.095 |
| **`rs_all`** | 1113 | **0.21751** | **+0.0413** | **81%** | 0.51842 | 0.96184 | +0.323 | **+0.480** |
| `rs_sched` (mandated delivery arm) | 2236 | 0.18228 | +0.0061 | 12% | 0.52111 | **1.02212** | +0.546 | +0.019 |
| **`rs_sched_perm` (CONTROL)** | 1911 | **0.19288** | **+0.0167** | 33% | 0.51354 | 0.97623 | +0.376 | +0.130 |
| `rs_hkey` (`rank_score` variant) | 1995 | 0.22584 | +0.0497 | 97% | 0.51211 | 0.98888 | +0.431 | +0.017 |
| `hy_kt2000` | 1579 | 0.22844 | +0.0523 | 102% | 0.51686 | 0.96749 | **+0.001** | **+0.601** |
| **`hy_kt4000`** | 2307 | **0.23211** | **+0.0559** | **110%** | 0.51320 | 0.97037 | +0.090 | +0.514 |
| `hy_kt6000` | 3020 | 0.21123 | +0.0350 | 69% | 0.51193 | 0.96873 | +0.127 | +0.506 |
| **`hy_kt8000`** | **3304** | 0.21142 | +0.0352 | 69% | 0.51305 | **0.96615** | +0.249 | +0.485 |

`hy_ktK` keeps the top-`K` genes by \|predicted effect\| on the parametric quasi law
(`param_mask`'s designed use) and resamples the rest. **`hy_kt8000` is therefore precisely
"a15 with its poisoned tail replaced"** — the `spread 1.0 / φ 4` tail block that 157 §5b showed
is 78.9% significant and 100% DOWN at zero effect, and that 171 §6 measured supplies **69.6% of
a15's reach head**.

### 4a. The ordering counterfactual — how much is realised

INBOX-4 priced reordering the *called* genes by `|H|` at **+0.036 raw** (0.21914 vs 0.18283) and
the full-pool key at **+0.0525** (0.23533). Both reproduce here to five decimals. On `base_mf`'s
own table the inside-block ceiling is **+0.03544** and the full ceiling **+0.05101**.

| arm | realised Δ`reach` | vs inside-block ceiling | vs full \|H\| ceiling | residual gap on its OWN table |
| --- | ---: | ---: | ---: | ---: |
| `rs_all` | +0.0348 (single-seed) | **98.3%** | 68% | 0.01055 (base: 0.03543) |
| `hy_kt2000` | +0.0523 | 147% | **102%** | **0.00026 twin spread**, gap 0.00722 |
| `hy_kt4000` | +0.0559 | 158% | **110%** | gap 0.00794 |
| `hy_kt8000` | +0.0352 | 99% | 69% | gap 0.01267 |

**The defect INBOX-4 found is closed.** The arms above the inside-block ceiling get there by
also reordering *outside* the call block, which is what the 0.23533 full-key row measures.

### 4b. The permutation control fires, and it fires against the schedule

`rs_sched_perm` uses the **same multiset of per-gene σ**, randomly permuted across the gated
axis — identical dispersion budget, identical σ histogram, only the gene↔σ pairing destroyed.
It reads **0.19288** against the schedule's twin mean **0.18228**, i.e. **+0.0106 HIGHER**, and
its `nmae` is 0.046 better. *(What would this read if the mechanism were present? The permuted
arm would score worse, having thrown away the compensation. It scores better — on both members.)*
This is 171 §12b's result reproduced on a different emitter.

---

## 5. `mse` — exactly zero, measured not asserted

`claim177.py`, real challenge panel, all four known-answer gates passing.

| policy | `Ĉ/C_p` | `claim` | branch | raw `mse` | Δ`mse` overall |
| --- | ---: | ---: | --- | ---: | ---: |
| a15 as shipped | 1.085 | 0.24581 | SATURATED | 0.88228 | — |
| **`rs_all`** | **1.094** | **0.24647** | SATURATED | 0.88228 | **0.000000** |
| `rs_sched` (+`fracpair` ×2) | 1.181 | **0.24653** = cap | SATURATED | 0.88228 | **0.000000** |
| `hy_kt2000` | 1.041 | 0.24502 | SATURATED | 0.88228 | **0.000000** |

The knee is `s²B = 0.23268` and every arm sits above it, so `mse` is exactly flat. The
`resample_mix` docstring's "`C_hat` rises to ~1.03 of the real-cell value" is confirmed in
direction and is slightly larger than claimed (1.085 → 1.094 for `rs_all`). **`mse` is free.**

---

## 6. `fid` and `jac` — where the arm actually dies, and it is rule 9

`nreal177.py`. Every previous emission arm held board `n_pred` near 4,424, where the yield cap
`min(1, n_pred/n_real)` is saturated at *both* ends of the inferred `n_real` range, so the cap
never mattered. **This experiment cuts `n_pred` by up to 5.4×, and the cap is now the verdict.**
Hard board-algebra bounds (153, harness-free): `jac ≤ n_real/n_pred` ⟹ **`n_real ≥ 191`**;
precision ≤ 1 on v14 ⟹ **`n_real ≤ 4,716`**. 153's enrichment meeting-point and a21 put it at
**~240**; 093's **~2,670** is flagged circular by rule 9 but not dead.

Sum of the four DE members in overall board units (`pds` and `mse` are exactly 0):

| arm | board `n_pred` | n_real 191 | **240** | 400 | 812 | 1200 | 2000 | **2670** | 4716 |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| `rs_all` | 812 | +0.0197 | **+0.0211** | +0.0239 | +0.0217 | −0.0807 | −0.1843 | **−0.2375** | −0.3280 |
| `hy_kt2000` | 1152 | +0.0121 | **+0.0125** | +0.0126 | +0.0074 | −0.1493 | −0.2139 | **−0.3001** | −0.3179 |
| `hy_kt4000` | 1683 | +0.0075 | **+0.0076** | +0.0073 | +0.0030 | −0.0696 | −0.1492 | **−0.2479** | −0.2687 |
| `hy_kt6000` | 2203 | +0.0046 | **+0.0045** | +0.0039 | −0.0001 | −0.0220 | −0.0879 | **−0.1991** | −0.2228 |
| **`hy_kt8000`** | **2410** | +0.0060 | **+0.0060** | +0.0056 | +0.0026 | −0.0158 | −0.0573 | **−0.1725** | −0.1968 |
| `rs_sched` | 1631 | +0.0003 | **+0.0020** | +0.0066 | +0.0133 | +0.0150 | −0.0444 | **−0.1139** | −0.2134 |
| `rs_hkey` | 1455 | +0.0061 | **+0.0073** | +0.0103 | +0.0125 | +0.0098 | −0.1496 | **−0.2355** | −0.2521 |

**Every arm is positive at the settled `n_real` and catastrophically negative at 093's.** The
swing on `rs_all` is +0.021 to −0.238 — a 0.26 range on a metric we score 0.173 overall. No
harness we hold can narrow it (rule 9); `hy_kt8000` is the most robust arm, staying positive to
`n_real ≈ 1,100`.

Two caveats on this table, both against the arm: `recall` is carried from h1-test where
`n_real` is 222, and `recall` is itself `n_real`-dependent; and the `jac` gain at small
`n_real` comes from the *union* shrinking with `n_pred`, which is the classic small-`n_real`
call-set-cut artefact that has already cost four board predictions.

---

## 7. GO / NO-GO

### NO-GO. Do not spend a slot on any arm here. a15 ships unchanged, zero flags different.

**Rule 13, applied literally.** This is an emission-law arm at a **verified** bit-identical
pseudobulk, so the **+0.0019** gate applies, and the composite is

```
predicted = 1.011 · Δpds/6   +   Δmse/6   −   0.0146/6   (the emission-law change penalty)
          = 1.011 · 0.000000 +   0.000000 −   0.00243     =  −0.00243
```

**−0.00243 < +0.0019 ⟹ VETO.** Rule 13 prices `reach`, `fid` and `jac` as *unpriceable* and
charges a flat 0.0146 for touching the emission law; an arm whose entire value lives in those
three members cannot clear it, by construction. The §6 table is why that is the right call
today and not merely a formality: the arm's value is `+0.021` or `−0.238` depending on a
quantity nobody has measured.

**The honest counter-argument, recorded so the next agent can weigh it.** Rule 13's "`reach`
2/6 signs, worse than chance" was fitted on arms whose `reach` deltas were *inside* the
same-policy seed spread (171 §7: 0.0075–0.0123 spread against deltas of the same size). These
deltas are **+0.035 to +0.056 raw**, 3–7× the base pair's spread and **200× `hy_kt2000`'s own
twin spread of 0.00026**. That is a different resolution regime, and 157 §1 / 165 §6 validate
h1-test for `reach` on exactly this arm class (emission law, fixed profile). If the reach
channel is real, `hy_kt8000` is worth `+0.0352 × 0.393 / 0.895 / 6 = +0.0026` on `reach` alone,
`hy_kt4000` `+0.0041`.

### What the mechanism buys, which is the real deliverable

The project has spent five experiments (093, 095, 157, 171, and the a15 tier) trying to deliver
the `|H|` ordering through a spread schedule, and 171 closed the family. **The emission law
delivers it outright**: `ρ(p, |H|)` inside our call set goes **−0.395 → +0.601**, `ρ(p, expr)`
**+0.484 → +0.001**, and 81–110% of INBOX-4's counterfactual is realised. That is the whole
+0.036–0.0525 that CURRENT_STATE called "our entire deficit below the field line," and it is
now a solved problem — **it just arrives attached to an `n_pred` cut nobody can price.**

### The production change, if the `n_real` question is ever settled

A new `emit_counts` path is a **builder** change, so it must be named and diffed (rule 10):

* `src/vcc2026/donor_transfer.py` — no code change; `mode="resample_mix"` already exists and is
  used unmodified here. **Two call-site contracts must be honoured or the arm inverts:**
  (a) `control_counts` = **raw** control counts with the cells' **own** `lib_sizes` — never
  rescaled to a nominal depth, or the exact-total repair peels the block flat (§0.1);
  (b) `control_profile` = the **pooled** control profile `expm1(b_ctrl)/Σ`, the space the
  rendered profile lives in — not the mean-per-cell-CPM profile, which differs by an sd of
  0.020 in log₂ and makes `fold ≠ 1` at zero predicted effect.
* `scripts/make_submission_pm_cd4.py` — new flags: `--emitter {quasi,resample_mix}` (default
  `quasi`), `--resample-keep-top K` (the `param_mask` cut; `K = keep_top = 8000` gives
  `hy_kt8000`, i.e. a15's keep block parametric and a15's tail resampled), and the control-cell
  source, which for the challenge is `data/challenge/context_{A,B,C}.h5ad` itself.
* Exactly **two** meta differences from a15 (`emitter`, `resample_keep_top`); everything else —
  `overdispersion 4.0`, `tail_spread 1.0`, `keep_top 8000`, `keep_spread 0.25`, `VCC_ALPHA
  0.982`, `whitener_id union_pruned` — unchanged. **`--min-cpm` must not move**; it changes the
  profile and leaves the regime h1-test can read (164 §7).
* Recommended point: **`hy_kt8000`** (most `n_real`-robust, `nmae` neutral at +0.0000,
  precision +0.0052, `reach` +0.0352). `hy_kt4000` if the 240 read is trusted.

### The one slot that would settle it

CURRENT_STATE already names the discriminator: **a fixed-profile, raise-`n_pred`-only
submission**, with `pds`/`mse`/`nmae` invariant by construction. This experiment adds that the
*converse* also works and is cheaper to read — `hy_kt8000` cuts board `n_pred` 4424 → 2410 at a
verified-identical pseudobulk, so its board `fid` alone reads `n_real` directly: `fid` rises to
≈ precision if `n_real ≤ 2410`, and falls to `precision × 2410/n_real` otherwise. That is a slot
that buys information rather than rank.

---

## 8. What I did NOT measure, and why

| not measured | why | can it move? |
| --- | --- | --- |
| `pds` on the delivered artifact | pinned by the sha256 identity on 85/85 targets (§1); BRIEF §1 defines it on exactly that pseudobulk | **No — identity** |
| `mse`'s profile half | same identity; the dispersion half IS measured (§5) and is exactly 0 | **No — identity** |
| challenge-side `n_pred` for any arm | needs a built `pred.h5ad`, i.e. a challenge count build — forbidden | It is *the* open variable; §6 brackets it |
| `n_real` on the challenge | rule 9; no harness we hold is in the regime, and 093's inversion is circular | Decides the sign of the whole arm |
| `hy_ktK` for K > 8000 | K = 10,779 is `base_mf` by construction; the frontier between 8000 and the gate is a two-arm interpolation | Would trade `reach` back for `n_pred` monotonically |
| a `φ` sweep on the parametric block of `hy_ktK` | rule 5, one variable per run; 171 §12c shows `φ` worsens the skewness pathology | Untested |
| the resample law at a genuine 5× depth | inflating real cells 5× is not a depth increase, it is a fabricated tightening — the confound §0.1 diagnoses. `quasi100` covers the depth question on the law where it is well defined | `quasi100` reads −0.0108 |
| Orion / HCT116 surfaces | 157 §1: repaired HCT116 misses a15's board `reach` by 2.6×; 165: h1-test is the validated surface for emission-law `reach` | Absolute level, not the deltas |
| `cell-eval2` on anything | forbidden (25 GB RSS, shared machine) | — |

### Integrity notes

* **Eleven known-answer gates, all passing** (§1), four of them bit-identical reproductions of
  157/164/173/INBOX-4 from fresh emissions on this machine.
* **The self-inclusion bias is real and quantified**: the 400 donor cells are drawn from the
  same 18,400 controls the Wilcoxon tests against, which shrinks `|F̄ − ½|` by a factor
  `(1 − N/m) = 0.978`. It is a 2.2% conservative bias on `|z|`, identical across all resample
  arms, and it makes the §2 null slightly *better* than the truth — so it is stated, not hidden.
* **Before trusting the permutation control**: if the schedule's mechanism were absent, the
  permuted arm would read the *same or better*. It reads better on `reach` **and** on `nmae`,
  on two members, which is not invariant by construction — a schedule that worked would have
  to beat it.
* **The `quasi20` control the mandate asked for was replaced, not skipped.** It is provably
  bit-identical to `base_mf` (§0.2) because `EMIT_L` never reached the emitted depth; an
  algebraic identity is not a control (this project has already been burned twice by exactly
  that), so `quasi100` was run instead and reports a real depth effect of −0.0108.

## 9. Files

| file | what |
| --- | --- |
| `null177.py` / `work/null177.A.log`, `work/null177.BC.log` | zero-effect null, all three challenge contexts, 157's known answers |
| `null177b.py` / `work/null177b.h1test.log` | noise-kernel × σ ceiling; the coordinator's confound as a KA row |
| `signal177.py` / `work/signal177.log` | signal response; the confound reproduced at 10,658/10,779 |
| `realize177.py` / `work/realize177.log`, `work/realize177b.log` | 16 realised arms, 85 perts, pseudobulk digests |
| `six177.py` / `work/six177.log`, `work/six177.json` | the six-member table, ordering diagnostics, `|H|` ceilings |
| `claim177.py` / `work/claim177.log` | `mse` through 096's `claim` route, four KA gates |
| `nreal177.py` / `work/nreal177.log` | the verdict across the admissible `n_real` range |
| `de177.*.npz` | 16 DE tables (85 × 10,779), real and predicted |

No slot spent. a15 unchanged.

---

## Coordinator addendum (22:50) — the per-perturbation K variant is closed before it is built

`hy_ktK`'s `jac` loss (−0.0142 on h1-test) lives entirely in the top `n_real` tercile (−0.042; the
other two terciles +0.0001 / +0.0005). The natural repair — a parametric block that grows with the
perturbation's predicted strength — requires predicted strength to track true `n_real`.
**Spearman(‖H_p‖, `n_real_p`) = −0.036; Spearman(a15 `n_pred_p`, `n_real_p`) = −0.085** (while
‖H_p‖ vs our own `n_pred_p` is +0.967). Our prediction carries no information about how many genes
truly respond, so no per-perturbation K can be keyed on anything we hold, and the `jac` loss is spread
across all predicted-strength terciles (−0.015 / −0.006 / −0.020). Only a global K remains, and at
K = 8000 the arm is ≈ neutral on expected score (see CURRENT_STATE 22:40). The manifest confirms the
regime match with h1-test (`ground_truth_cells` 138,400 = ~461 cells/perturbation; controls 18,400).

## Coordinator addendum 2 (23:00) — cross-check against 164's tail-φ0 arm

164's `armA` (a15 with the tail block at `spread 1.0 / φ 0`) has **the same DE profile as `hy_kt8000`**:

| | n_pred | reach | fid | jac | nmae | mse (overall) |
|---|---:|---:|---:|---:|---:|---:|
| a15 | 6034 | 0.1828 | 0.5070 | 0.0682 | 0.9666 | 0 |
| 164 armA tail φ0 | 3151 | 0.2076–0.2100 | 0.5046 | 0.0538 | 0.9666 | **−0.0185** (claim 0.246 → 0.171, below the knee) |
| 177 `hy_kt8000` | 3304 | 0.2114 | 0.5061 | 0.0541 | 0.9662 | **0.0000** (Ĉ/C_p 1.094, above the cap) |

Both remove the tail's zero-effect false calls and both pay the same `jac`/recall price. 164 rejected
armA because the quieter tail also *lowers the emitted dispersion* and forfeits the jackknife credit;
**resampled real cells carry the control's own dispersion, so 177 gets 164's DE gain with the `mse`
cost removed.** That is the whole advance, stated in one line — and it confirms the arm's expected
value is the DE-only sum: reach ≈ +0.0026, fid ≈ −0.001, jac ≈ −0.003 to −0.006 depending on the
challenge's `n_real` tail. A wash on score; the mechanism and the regime information are what it buys.
