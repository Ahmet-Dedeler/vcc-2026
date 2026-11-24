# ADDENDUM to the `h1union` preregistration — written 2026-09-04 02:25 UTC, BEFORE the build finished and BEFORE submission

The preregistration is **not edited**. This sharpens three falsifiers and records one deterministic
consequence, all from artifacts already on disk. Prompted by a second Fable 5.1 session Ahmet asked
to review the day's work; I have verified each claim in the source rather than accepting it.

## 1. Falsifier 1 (`pds`) becomes THREE-WAY, not binary

As written — *"`pds` does not rise"* — the falsifier fires identically on a **transfer failure** and
on a **mechanism failure**, and we would not know which we had bought. `cd4f034` paid **−0.00365 raw
`pds`** for +0.018 cosine, and PLOLO `pds` transferred ~2× pessimistically on it. Against
`cd4f034`'s board raw `pds` of **0.83228**:

| board raw Δ`pds` | reading |
|---|---|
| **≥ +0.0010** | the **fourth-independent-estimate** mechanism is confirmed; the coverage-aware rule *and* the sources stand |
| **−0.0010 to +0.0010** | **transfer-limited**; the rule stands, the "4th estimate" claim is **not** adjudicated |
| **≤ −0.0015** | **dose-like** — the sources are adding *shared* direction on their rows, and the `g` ladder must go **down**, not up |

## 2. Falsifier 2 (`mse`) becomes two-sided, from 197's identity

`raw mse = 1.023611 − 2.074495 × cos`. `cd4f034` sits at cos **0.09027**, raw **0.83635**. 218's
macro **+7.68%** transferred 1:1 gives raw **0.822**; at 60% realisation **0.827**; at 219's 1.38×
transfer **0.817**.

**Point prediction: raw `mse` ≈ 0.822, band [0.817, 0.827].** "Must improve on 0.83635" is too weak
— a landing at 0.834 would pass while telling us the transfer factor for a 26%-row arm is a quarter
of the full-row one. **A pass inside the band and a pass at 0.834 are different results.**

## 3. The preserve scalar is a declared consequence, not a violation

From `work/stage-on028.json` and `work/stage-off.json`: `post_add_global_preserve_scale` is
**0.920359** for `cd4f034` and **0.909357** for the arm — ratio **0.9880**. So **every one of the
223 non-carrying rows, and CD4's own delivered contribution, is 1.2% smaller than in `cd4f034`.**

Through 187's identity at our `k ≈ 0.75` (slope −0.074 below `k = 1`) that is ≈ **+0.0005 raw
`nmae` on 74% of rows, ≈ −0.0001 overall** — inside noise.

**This is recorded now because the mandate says "identical elsewhere" while `profile222.py`
adjudicates it in the scale-corrected form** (non-carrying blocks differing by one shared scalar).
Those must not be allowed to disagree after the fact, and **the 1.2% shrink on CD4's delivered
amplitude must never be mistaken for CD4's dose moving** — `cd4_target_f` stays 0.34 and
`effects.cd4.*` stays untouched, which falsifier 5 checks.

## 4. Why `g = 0.28` and not LOFO's 0.45 — stated, since the preregistration did not

LOFO chose **0.45** on all four folds on a **cosine-weighted** price. `compose218.json`'s raw `pds`
per fold, base → 0.17 → 0.28 → 0.45:

| fold | base | 0.17 | 0.28 | 0.45 |
|---|---:|---:|---:|---:|
| `k562_gwps` | 0.79588 | 0.80339 | 0.80540 | 0.80642 *(still rising)* |
| `orion_hct116` | 0.79014 | 0.79202 | **0.79221** | 0.79131 *(turns over)* |
| `orion_hek293t` | 0.76055 | 0.76324 | **0.76362** | 0.76179 *(turns over)* |

**`pds` turns over between 0.28 and 0.45 on two of the three cancer folds** — the recipients that
resemble D/E/F. 0.28 is at or just before the discrimination optimum; 0.45 is past it on the
majority. This is 219's own lesson about local propagation factors, applied before the fact rather
than after.

**One ambiguity worth naming:** on the carrying rows at 0.28 the cosine gain is **+16.7% / +17.3% /
+12.7%** on the cancer folds, *above* 159's `1/k` prediction of **+13.7%** for k 3→4 (218 quoted the
law's agreement at `g = 0.17`, where the rows gain +8.8–11.7%). That is consistent with the H1
sources being **better-than-average estimators** (`h1_2025` runs 400 cells per perturbation) — and
it is *also* what an **over-dose** looks like. **§1's three-way `pds` reading is the discriminant.**

## 5. Pre-committed: the concentration check
220 found this morning that both published ordering effects were carried by **one to three
perturbations**. 218 never reported a leave-top-k-targets-out on the union's Δcos or Δ`pds`.
**Pre-commitment: if the top three of the 77 carried targets carry more than half of Δ`pds` on the
cancer folds, a positive board result is weaker evidence for 223's stack**, because the board would
be measuring a handful of perturbations rather than a source.
