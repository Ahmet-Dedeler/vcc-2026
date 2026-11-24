# a23 — mass-consistent rendering. Pre-registered 2026-09-02, before the build.

> ## AMENDMENT, 2026-09-02 01:20 UTC — written BEFORE the slot was spent
>
> Everything below was scaled off the two held Orion artifacts. I have since measured the
> mass ratio on the **challenge frame itself** (`challenge_frame.py`, contexts A/B/C, a15's
> CPM>10 gate, the effect matrix rebuilt and its `post_whiten_sha256` confirmed identical to
> a23's `b2426e15…`):
>
> | | mass ratio | render residual, % of row norm |
> | --- | ---: | ---: |
> | challenge A / B / C | 1.001339 / 1.001347 / 1.001273 | **2.83% / 3.35% / 2.94%** |
> | held Orion HCT116 / HEK293T | 1.00181 / 1.00168 | 3.96% / 4.53% |
>
> **The challenge's contamination is ~72% of the harness's**, so the prediction scales down.
> Revised: cosine gain **+9.8%** (pessimistic) to **+41%** (optimistic), giving
> **`mse` +0.0026 to +0.0107 overall**, `pds` +0.0001 to +0.0008, `nmae` ~+0.0004.
>
> **Revised total: `+0.004` overall (range `+0.003` to `+0.011`), landing ≈ `0.178`.**
>
> Falsifier 1 is correspondingly relaxed from raw `mse` < 0.8640 to **raw `mse` < 0.8680**
> (Δraw ≥ 0.0058, ~40% of the pessimistic prediction). This is an amendment to an *input*
> measurement made before any outcome was observed, not a post-hoc rescue; the original
> threshold is left in place below so the change is visible.
>
> **Verification of the shipped artifact** (the naive check was a tautology — see §8): the
> residual `r = (effect recovered from the emitted counts) − (intended effect)` points at the
> min-norm correction on **12/12** probed targets, mean cosine **+0.700** to `eps_minnorm`
> against **+0.338** to the renderer's implicit shift. The gap from 1.0 is integer-emission
> rounding, which adds an independent term.

**Arm.** a15 exactly, plus `VCC_MASS_FIX=minnorm`. Config diff vs a15 must be **exactly one
field** (`mass_fix: "minnorm"` added). Nothing upstream of the renderer changes, so the effect
matrix is bit-identical and `effects.whitening.after_all_median_norm` must read exactly
`4.2292956044543` and the operator sha exactly `e2f65ecf…`. If either drifts, the build is not
this arm and must not be submitted (AGENTS.md 11.1 rule 10).

## The mechanism

`profile_from_effect` renders the log-space effect as `expm1(b_ctrl + gated)` and divides by its
own sum. The scorer reads back `log1p(TS·p) − b_ctrl`. That division is the identity **iff**
`Σ_g expm1(b_ctrl + gated) = TS`. The control profile satisfies it exactly by construction, and
so does every reference perturbation — the truth lives on that manifold. Our effect does not:
measured mass ratio **1.00181 / 1.00168** on the two frozen Orion artifacts.

So the renderer applies its own correction, and it lands as a shift concentrated on the
highest-expressed genes that is nearly the same direction for every target: **3.96% / 4.53%** of
the row norm at **row-row cosine 0.303**, against an *intended* row-row cosine of **0.0040**.
Whitening exists to make those rows mutually dissimilar. The renderer was putting the shared
blob back, and every effect-level screen in this project since exp 017 scored the effect
*before* that happened.

`expm1` is convex with `S(0) = TS`, so `S` has no second root along `{γ·effect}` — scaling
cannot satisfy the constraint, the correction must move off the ray. But it is one scalar, so
the direction is free: `dS/dh_g = exp(b_g + h_g)`, and the minimum-L2 correction is
`ε = −(S−TS)·w/‖w‖²`, `w_g = exp(b_g+h_g)`, confined to the CPM gate. Cost: **1.76% / 2.29%**
of the row norm, and the read-back becomes exact to `1e-12`.

## Measured, two held Orion lines, scorer's own profile space (`experiments/156-*/`)

| | HCT116 | HEK293T |
| --- | ---: | ---: |
| Δ pds | `+0.00035` | `+0.00278` |
| transfer cosine (macro) | `0.02147 → 0.02460` **+14.6%** | `0.00784 → 0.01311` **+67.2%** |
| aggregate cosine (from `sse*`) | `0.02285 → 0.02596` **+13.6%** | `0.00900 → 0.01404` **+56.0%** |
| row-row cosine | `0.00564 → 0.00484` | `0.00480 → 0.00418` |

**The control that makes this a mechanism rather than a fit.** `massfix-const` satisfies the
*identical* scalar constraint with a uniform log shift instead of the minimum-norm one. It
**loses** `−0.00555 / −0.00312` pds and `−4.7% / −18.8%` cosine. Same constraint, different
direction, opposite sign of result — so the effect is geometric, not normalisation bookkeeping.
`massfix-minnorm` also beats the *infeasible* `intended` upper bound on HEK293T
(`+0.00278` vs `+0.00237`), i.e. the feasible object is closer to the truth than the one we
designed, which is what you expect if the truth's own manifold is the constraint set.

## All six members priced (AGENTS.md 11.1 rule 12)

| member | estimate | basis |
| --- | ---: | --- |
| `mse` | **+0.0036 to +0.0147**, point **+0.006** | `1.0237 − 2.077·cos`, our own board-fitted model (out-of-sample +0.32% on a8). cos `0.0722 → 0.0820` (worst) / `0.1126` (best) |
| `pds` | **+0.0001 to +0.0010** | worst-surface `+0.00035` raw / span 0.452 |
| `nmae` | **+0.0006** | debench raw `0.9501 → 0.9479` on HCT116, one surface |
| `reach` | **≥ 0**, not counted | `+0.00598` raw keyed on \|effect\|, **not** the shipped p-order; and rule 9 makes harness reach magnitudes unreadable |
| `fid` | **≥ 0**, not counted | sign accuracy on the reference DE pool `+0.0106`; our call set is a different cohort |
| `jac` | **~0** | call set moves only through near-threshold genes |

**Total predicted: `+0.005` overall (range `+0.003` to `+0.018`).**
Base is **a15 = `0.1740865087`**, not the displayed a22 (`0.1715689`), because this arm sets no
resample flag. So predicted landing ≈ **`0.179`**.

Supporting evidence on the DE side, HCT116, deterministic: sign accuracy `0.5311 → 0.5417`
(**+0.0106**), head purity `0.6751 → 0.6844`, and of the **318** (pert, gene) pairs in the
reference DE pool whose predicted sign moved, **68.9%** land on the *correct* side against 0.500
chance. The fix is removing a systematic bias, not reshuffling noise.

## Falsifiers, fixed before the slot is spent

1. **If raw `mse` is not below `0.8640`** (i.e. scaled `mse` gains less than `+0.010`), the
   cosine→mse transfer model is wrong for this lever and the central claim fails.
2. **If scaled `pds` moves outside `[−0.001, +0.005]`**, the harness pds reading does not
   transfer for profile-channel changes.
3. **If overall lands below a15's `0.1740865087`**, this is a NO-GO and an "all six positive on
   two surfaces" profile-channel screen is not trustworthy — which would be worth more than
   the arm.
4. **If scaled `nmae` falls**, the sign-flip accounting is wrong.

## What is deliberately NOT in this arm

`massfix + centroid d1.0` scores better on worst-surface pds (`+0.00237` vs `+0.00035`) but
*lower* aggregate cosine on HCT116, so it trades our largest gap (`mse`) for our smallest
(`pds`), and it is two variables. Held for a clean single-variable follow-up if a23 confirms.
Note `centroid` alone, pre-render, is worth exactly `+0.00000` on HCT116 — independently
matching a competitor's reported board failure of the same transform, and consistent with the
separate finding that our panel matrix already has only 0.35% shared energy.

## 8. The verification that failed, recorded because the failure is the lesson

My first artifact check recovered `δ̂ = log1p(TS·p) − b_ctrl` from the emitted counts and tested
`Σ expm1(b_ctrl + δ̂) = TS`. a23 read **exactly 1.00000000**. That looked like a clean
confirmation and it is **worthless**: `expm1(log1p(TS·p)) = TS·p`, which sums to `TS` for *any*
normalised profile. It was caught only because I ran the same check on **v9**, a pre-fix build,
which also read exactly 1.00000000.

This is AGENTS.md 11.1 rule 10's corollary operating exactly as written — a screen with no
known-answer row cannot detect this class of bug, because the wrong number looks precisely as
reasonable as the right one. The real check needs the intended effect matrix rebuilt
independently, which is what the amendment above reports.
