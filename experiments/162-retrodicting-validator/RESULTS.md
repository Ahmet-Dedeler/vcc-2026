# 162 — a validator that retrodicts board `pds` to ±0.0015, and is blind to the other five

**Date:** 2026-09-02 · PLOLO (panel-restricted leave-one-library-out), 4 donor folds over the
real 300-target panel, deterministic pseudobulk truth, no emitter, scored after the renderer.

## Verdict

**GO — as a veto instrument on the `pds` channel of pure-geometry arms.** Across the six
reconstructible geometry arms, PLOLO gets **6 of 6 signs right**, a fitted scale of **0.899**
(i.e. it is very nearly 1:1, not 5.5× hot), a **leave-one-arm-out RMS residual of 0.001494
scaled `pds`**, against a predict-zero null of 0.003481. Exact permutation test over all 720
relabelings: **p = 0.0028**. That is the GISL-described capability, at the claimed accuracy
(~0.001–0.003), reproduced — **on one member of six.**

**NO-GO — as authorisation to spend a slot.** The five members it cannot price carry
**RMS 0.007350 on overall** across those same six arms, against a best-case priced signal of
**0.000997**. That is a **7.4× ratio**. Predicting board *overall* from PLOLO alone gets
**4 of 6 signs** — barely better than a coin. This validator kills bad arms; it does not
green-light good ones.

The previous runtime's headline that a23 was over-predicted **5.5×** was a **single-fold
artefact**. Over all four folds the a23 prediction is +0.004014 scaled against a realised
+0.001737 — 2.3× — and the residual is inside the instrument's own error bar.

---

## 1. What PLOLO is, and the fold record

Hold each donor library out as truth over the real 300-target panel; rebuild the production
effect matrix from the remaining libraries with a whitening operator fitted **without** the held
line; render through production's `profile_from_effect` into the **challenge** control frames
A/B/C; read back the way the scorer does; score against the held library's own pseudobulk.

| fold | sources | operator | sha | panel | scored genes | PM targets |
| --- | --- | --- | --- | ---: | ---: | ---: |
| `orion_hct116` | k562, hek293t | `062/held_orion_hct116` | `7cecbf7f` | 290/300 | 17 806 | 239 |
| `orion_hek293t` | k562, hct116 | `062/held_orion_hek293t` | `97e486e0` | 294/300 | 17 806 | 242 |
| `k562_gwps` | hct116, hek293t | `062/held_k562_gwps` | `13d988fc` | 269/300 | 7 448 | 0 |
| `kolf2` | all three | `data/donors/whitener_union_pruned` | `e2f65ecf` | 273/300 | 17 806 | 228 |

`kolf2` is not a production donor and contributed no row to any production operator, so that
fold is leak-free while using the exact shipped operator and the exact shipped three-source
mixture. `h1_2025` is excluded entirely (20/300 panel coverage, and it is the oracle-leakage
donor of the §6.4 harness warning). CD4 is a production input on every fold and never a truth.

**Predictor, fixed in advance:** unweighted mean over all four folds, challenge frame (mean of
A/B/C). No fold is selected. The frame was pre-committed by the previous runtime on the
mechanistic ground that it is the frame production renders into, and its render residual
(**2.88 / 2.90 / 2.87 %** on A/B/C, mass ratio 1.0018) reproduces `A23-PREDICTION.md`'s
independently measured **2.83–3.35 %**. See §7 for the disclosure on that choice.

## 2. The deliverable — predicted vs realised, per member

All values are **scaled** board units. `predicted` is the 4-fold mean before calibration;
`fold sd` is the between-fold scatter of that mean's terms; `LOO pred` applies a one-parameter
through-origin scale fitted on the *other* arms only.

### `pds` — RETRODICTS

| arm | predicted | fold sd | realised | LOO pred | LOO resid | k₋ₐ |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| a11 | +0.003722 | 0.001425 | +0.003328 | +0.003350 | −0.000022 | 0.900 |
| a12 | −0.002763 | 0.001033 | −0.002004 | −0.002532 | +0.000528 | 0.917 |
| a14 | +0.000354 | 0.000134 | +0.000729 | +0.000318 | +0.000411 | 0.897 |
| a15 | +0.006656 | 0.002766 | +0.007257 | +0.004527 | +0.002730 | 0.680 |
| a19 | −0.000994 | 0.005320 | −0.001183 | −0.000890 | −0.000293 | 0.895 |
| a23 | +0.004014 | 0.016538 | +0.001737 | +0.004060 | −0.002323 | 1.011 |

**k = 0.899 · sign 6/6 · LOO RMS 0.001494 vs null 0.003481 · permutation p = 0.0028.**

### `mse` — MEASURED AND REJECTED

`a11` is excluded: its base `a10` sits at scaled `mse` **exactly 0.0**, which is the clamp floor
(`mse` is clamped to [0,1]), so the realised +0.082041 is a censored lower bound, not a
measurement.

| arm | predicted | fold sd | realised | LOO pred | LOO resid |
| --- | ---: | ---: | ---: | ---: | ---: |
| a12 | −0.111697 | 0.043166 | +0.009122 | −0.221623 | +0.230745 |
| a14 | +0.020176 | 0.007874 | +0.019452 | +0.004461 | +0.014991 |
| a15 | +0.001328 | 0.000850 | +0.009442 | +0.000320 | +0.009122 |
| a19 | −0.043650 | 0.012409 | −0.095908 | +0.002078 | −0.097986 |
| a23 | +0.000257 | 0.000771 | +0.001524 | +0.000062 | +0.001462 |

**LOO RMS 0.112388 vs null 0.044162 — 2.5× WORSE than predicting zero. p = 0.958.**

### `nmae` — MEASURED AND REJECTED

sign **3/6**, LOO RMS 0.005112 vs null 0.004013, **p = 0.892**.
Fails at every pool size: K=120 → 0.00523, K=240 → 0.00511, K=480 → 0.00475, all above the
0.00401 null. The failure is not a bad choice of pool size.

### `reach` — MEASURED AND REJECTED

| arm | predicted | realised | | arm | predicted | realised |
| --- | ---: | ---: | --- | --- | ---: | ---: |
| a11 | +0.000215 | −0.011375 | | a15 | +0.000875 | +0.000826 |
| a12 | −0.000205 | +0.004570 | | a19 | +0.000370 | +0.002970 |
| a14 | +0.000033 | −0.002248 | | a23 | +0.000057 | **−0.007933** |

sign **2/6** (worse than chance), fitted scale **negative** (k = −2.108), LOO RMS 0.007839 vs
null 0.006161, **p = 0.614**. Fails at every K (0.00796 / 0.00784 / 0.00689 vs null 0.00616).
Predicted magnitudes are 10–100× too small.

### `fid`, `jac` — UNMEASURABLE, and they are not exempt

PLOLO emits no cells. Both members are functions of *call sets*: `fid` is precision on our own
BH-corrected significant set (`n_pred`), `jac` is the Jaccard of ours against the reference's.
With no cells there is no Wilcoxon, no p-value, no call set, and no `n_real`. This is **not** a
rule-12 "cannot move" exemption — a21 moved `fid` −0.029464 and `jac` −0.016582. They are
unmeasured members that must be priced by a separate instrument (experiment 164 owns that
surface).

## 3. Why `mse` fails — the mechanism, and it is structural

PLOLO's truth is a **different cell line**. The transfer cosine is 0.008–0.032, and the global
scale that minimises PLOLO's own squared error is `gstar` = **0.0215 / 0.0323 / 0.0458 /
0.0576** across the four folds — that is, PLOLO wants an amplitude **17–46× smaller** than
production ships. Its absolute `nmse_row` at production amplitude is **1.18–1.47**, i.e. *worse
than predicting no effect at all*, while the board's raw `mse` is 0.874. The two live in
disjoint regimes.

So PLOLO's gradient in α always points down, and every alpha arm is mis-priced: it reads a11
(α down) as an improvement and a12 (α up) as a regression, while the board says the optimum is
near 0.982. **A cross-line truth cannot locate an amplitude optimum.** This is the same fact
156 hit from the other side ("do not price a future profile arm by scaling a harness cosine
ratio").

Note the earlier endpoint set was *worse* than this: `run_arms.py` recorded only `sse = 1−cos²`
at the best global amplitude, which is scale-**invariant** by construction and read a19's
−0.0959 `mse` catastrophe as **+0.0000148**. Adding a scale-sensitive `nmse` (this experiment)
at least recovers a19's sign. It still fails the null.

### A flagged inconsistency in the board record itself

a11/a12/a14 nominally differ only in α (a12/a14 additionally carry the tier-8000 spread policy,
which rule 12 exempts as `mse`-neutral). Raw `mse` must then follow `A·α² − 2B·α + 1`, pinned at
1 for α=0. Fitting that law over the entire documented (b, span) range — 143 combinations —
leaves a best-case residual of **0.0116 raw `mse`, about 8× the size of a genuine board move**
(a23 moved 0.001457 raw). The three points are not on a common parabola. Either the tier-8000
spread policy is **not** `mse`-neutral, or the per-context b/span vary enough to matter. Worth
one follow-up; it also means there is no clean α→`mse` target for a validator to be fitted
against.

## 4. Why `reach` cannot be priced here — answering the mandate head-on

**The validator is profile-channel-only. It must be paired with a separate DE instrument.**
That is now a measurement, not an assertion:

1. **The spec forbids the substitution.** `BRIEF.md` §5: the reference's DE pool is ordered by
   *our significant calls first, then ascending predicted adjusted p-value*, and `|lfc|` is only
   the **fourth** tiebreak. PLOLO has no cells, hence no p-values. `reach_k` substitutes
   `|effect|` order — which prices the ordering's *sign content* only.
2. **The substitution measurably inverts.** 2/6 signs, negative fitted scale, fails the null at
   every K. On a23 specifically, PLOLO reads **+0.000057** where the board read **−0.007933**;
   156 independently measured the same inversion (its `|effect|`-keyed reach *rose* +0.00598).
3. **The dominant `reach` channel is invisible by construction.** PLOLO predicts exactly 0.0 for
   the four emission arms. The board moved:

| arm | `pds` | `mse` | `nmae` | `fid` | `reach` | `jac` | overall |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| a8 | +0.000000 | +0.000000 | −0.000038 | +0.006393 | +0.021398 | −0.003439 | +0.0040523 |
| a10 | +0.000000 | +0.000000 | +0.000037 | +0.000797 | −0.000754 | −0.001140 | −0.0001766 |
| a21 | −0.000000 | −0.000000 | −0.064214 | −0.029464 | −0.018482 | −0.016582 | −0.0214572 |
| a22 | −0.000000 | −0.000000 | +0.001520 | −0.004447 | −0.010656 | −0.001523 | −0.0025176 |

`pds` and `mse` are **exactly 0.0 to sixteen digits on all four** — the channel classification is
confirmed. But the RMS of what PLOLO calls zero is `nmae` **0.032116**, `fid` 0.015243, `reach`
0.015113, `jac` 0.008521. These four rows are deliberately **excluded from every RMS above**:
predicting zero and reading zero is a test of the classification, not of the instrument, and
folding them in would deflate the error bar with free points.

**A validator that cannot produce a `reach` number cannot have saved a23 on its own.** But
PLOLO *plus* a board-calibrated toll can, and that combination retrodicts leave-one-out (§5).

## 5. The gating rule

**Gate 0 — eligibility.** The arm must be pure geometry: no change to the emitter, the count
law, or the per-cell library sizes. Otherwise PLOLO predicts 0.0 and the board can move
−0.064 (`nmae`, a21). Read the *finished artifact's* meta, not the command typed (rule 10).

**Gate 1 — `pds` magnitude.** Require calibrated |0.899 · predicted Δ`pds`| ≥ **0.00246 scaled**
= 1.645 σ at σ = 0.001494. Under a Gaussian residual that is a one-sided **5 % false-sign
rate**. Empirically there were **0 sign errors in 6 arms at any magnitude**, but with n=6 that
only bounds the true sign-error rate at ≤39 % (Jeffreys, 95 %) — the permutation p = 0.0028 is
the stronger statement, and the 5 % figure is a model estimate, not an observed rate. Of the six
historical arms only a11 (+0.00333) and a15 (+0.00726) clear this bar.

**Gate 2 — amplitude exposure.** If the arm changes the effect amplitude (α, γ, λ, any
renormalisation), PLOLO's `mse` channel is invalid (§3). Unpriced RMS on overall is **0.008903**
for amplitude-changing arms vs **0.001877** for shape-only arms. Refuse the slot or price
amplitude elsewhere.

**Gate 3 — emitted-profile toll.** Any arm that changes the emitted profile pays a `reach` toll.
Calibrated on a21+a22 only: **−0.014569 scaled**, i.e. −0.002428 on overall.

**Gate 4 — the budget, and why nothing passes it yet.** Net expected Δoverall =
0.899·x/6 + toll. For a shape-only arm the unpriced risk is 0.001877 on overall, while the
largest priced move PLOLO has ever produced is **0.000997**. To clear its own noise a shape-only
arm needs predicted Δ`pds` > **+0.0126 scaled** — nearly **2× the largest value in our entire
history** (a15, +0.00666). **Under this gate none of the six historical geometry arms earns a
slot.** That is the honest answer, and it is why the verdict is a veto instrument.

### The a23 veto test — leave-one-out

| quantity | value |
| --- | ---: |
| `reach` toll from a21+a22 **only** (a23 held out) | −0.014569 |
| PLOLO Δ`pds` for a23, scale fitted **without** a23 (k=1.011) | +0.004060 → +0.000677 on overall |
| **composite prediction** | **−0.001752** |
| **realised overall** | **−0.000727** |

**VETOED, correct sign, 2.4× hot.** Nothing about a23 entered either term. The caveat: the
*decision* to apply an emitter-derived toll to a profile-only arm is a judgement that was not
made in advance, though the law was on paper before a23 (a22's write-up: "resampling a gene
changes its emitted distribution, which changes its p-value, which perturbs its rank"). Adding
this toll takes overall sign agreement from **4/6 to 5/6**; the one remaining error is a12,
which is an `mse` failure.

## 6. Controls — each one could have failed

| control | what it would read if the thing were absent | result |
| --- | --- | --- |
| board transcription: mean of 6 member deltas == overall | any mis-transcribed member breaks the identity | **PASS 12/12** |
| `a15` vs `a15-dup`, field-identical configs | state carried between arms → differing endpoints | **PASS, bit-identical, 20 rows** |
| pass-2 re-derivation of `pds`/`sse`/`gstar` vs pass 1 | fold-cache / operator / CD4 / config drift → 1e-3+ | **PASS, worst 2.1e-7; `pds` bit-identical** |
| predict-zero null on every member | a validator that explains nothing ties the null | **`pds` beats it 2.3×; `mse`/`nmae`/`reach` LOSE** |
| exact permutation, all 720 relabelings | destroys the arm↔prediction pairing, which is the claim | **`pds` p=0.0028; others p=0.61–0.96** |
| drop each arm and refit | a result carried by one arm collapses | **`pds` survives: −a23 → LOO 0.000859 p=0.0083; −a15 → 0.000997 p=0.0167** |
| drop each fold and refit | a result carried by one fold collapses | **`pds` 3 of 4 drops still beat the null; dropping `k562_gwps` does NOT (0.00435 > 0.00348)** |
| K ∈ {120, 240, 480} on the DE pools | a failure that is really a knob artefact | **`nmae`/`reach` fail at every K** |
| emission arms | if `pds`/`mse` were emitter-sensitive, they would move | **exactly 0.0, 16 digits, 4/4** |

Two honest failures inside the control set:

- **The `pds` result depends on `k562_gwps` being in the fold mean.** Three of the four
  three-fold subsets still beat the null; the one that drops `k562_gwps` does not. With only
  four folds the instrument is near the edge of its resolution.
- **The known-answer `intended` > `delivered` row (156 measured +0.00324/+0.00382) reproduces on
  only 3 of 4 folds** in the challenge frame: hct116 +0.0173, hek293t +0.0174, kolf2 +0.0055,
  **k562_gwps −0.0035**. Magnitudes on the Orion pair are 4–5× larger than 156's, consistent
  with the decoy pool going from 132 to ~290 (`pds` is a retrieval metric). Reported rather than
  laundered.

Also measured: `a23`'s `mass_consistent_effect` is a fixed-point iteration with an `rtol=1e-12`
exit, so a 2e-9 input shift changes which iteration exits and flips the sign of genes whose
delivered effect is ≈0. Drift on `sgn`/`sgn_w`: **6.1e-3 / 8.7e-3 on a23, exactly 0.0 on every
other arm.** Any sign-keyed statistic read off a mass-fixed prediction carries a ~1 %
reproducibility floor that has nothing to do with the model.

## 7. Disclosure on the frame choice

The challenge frame was pre-committed by the previous runtime, but that runtime had already seen
a23's board result and used it (challenge 5.5× vs library 14×) as part of the justification. So
the frame is not fully blind. Mitigations, all measured: the mechanistic justification (it is
the frame production renders into; its render residual reproduces an independent measurement) is
board-independent; the permutation test is within-frame; and the `pds` result survives dropping
a23 entirely (LOO 0.000859, p = 0.0083, 5/5 signs). Post-hoc, `library` fails the null (LOO
0.00507 > 0.00348) and `intended` has a lower RMS (0.00103) but only **4/6** signs, which is not
a usable gate.

## 8. What I did not measure, and why

- **`fid` and `jac`: not measured at all.** No cells → no call set. They are **not** exempt
  under rule 12 (a21 moved them −0.029464 / −0.016582). Owned by experiment 164.
- **`a13`** (donor-agreement per-gene + per-pert shrinkage): not reconstructible in this
  pipeline. Excluded from every fit.
- **`a2-prune` per-member:** v11's member decomposition was never recorded, so only its overall
  (+0.0001257) exists. Carried as an overall-only row; its predicted Δ`pds` is +0.001360 on the
  one fold (`kolf2`) that can express a `union → union_pruned` swap.
- **`a11`'s `mse`:** realised value is censored at the clamp floor. Excluded from the `mse` fit.
- **The emitter, end to end.** PLOLO has none. Every emission arm is predicted as exactly zero.
- **The emitted p-value order** (`|F̄_g − ½|/√tc_g`), which is what would actually fix `reach`.
  It is 164's mandate ("you own the CALL SET and the emitted p-value law") and I did not build
  it, to avoid two agents owning one surface.
- **Per-context truth.** PLOLO renders into A/B/C but scores against one donor library, so
  context-specific truth differences are not modelled; the A/B/C spread here is a render effect
  only.
- **`cell-eval2`: never run** (hard constraint). **Nothing was submitted; no `vcc submit`, no
  `vcc prep`, no challenge package, nothing written to `submissions/`.** No file outside
  `experiments/162-retrodicting-validator/` was modified.
- **Not attempted:** a same-line `mse` instrument. §3 argues a cross-line truth structurally
  cannot locate an amplitude optimum, which means fixing `mse` needs a *different* surface, not
  a better endpoint on this one. That is the highest-value follow-up.

## 9. Files

| file | what |
| --- | --- |
| `board.py` | board record transcribed from primary sources, each with a `src`; asserts mean(6 members) == overall |
| `plolo.py` | fold builder + production-verbatim renderer, library and challenge frames |
| `run_arms.py` | pass 1 — scale-invariant endpoints, 185 rows → `work/arms.json` |
| `rescore.py` | pass 2 — adds `nmse_row`, `nmae_k`, `reach_k`; gates against pass 1 → `work/arms2.json` |
| `retrodict.py` | the deliverable: predicted-vs-realised, LOO, permutation, gate arithmetic |
| `work/retrodict.log`, `work/retrodict.json` | the tables above, reproducible by `uv run python retrodict.py` |
