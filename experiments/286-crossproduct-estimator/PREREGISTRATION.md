# 286 — PREREGISTRATION (frozen before any held number)

Written 2026-09-06, before Stage A, B or C was run. `est286.py` (the estimator family and the controls) is frozen with
it; both sha256 are recorded in `work/prereg286.sha256` and asserted by every stage script. Nothing in this file or in
`est286.py` may change after the first held number; a change is recorded as an amendment with its own hash and the
reason, and every number measured before it is reported separately.

Binding: `docs/RESEARCH-PROGRAMME-2026-09.md` (freeze the family and the controls before any held number; name the
unit; strong baselines; the already-failed constructions are baselines here), `docs/READ-FIRST.md` §3 (pricing on
COHERENT PESSIMISTIC scenarios; uncertainty language), §4 (instruments), §5 (closed families). Reserves: BOTH SPENT —
`experiments/267-*/work/reserve/` and the Broad-16 lines JVE127 / MG63 / SLR23 / UMRC3 are never read.

## 0. The claim

A cross-product between two INDEPENDENT views of the same perturbation has expectation Σ_signal, so the DIRECTIONS of
the shared response can be estimated without the within-view sampling noise inflating them. A full-matrix (not
diagonal) combination built on those directions, with a floor so no direction is zeroed, is a DIFFERENT object from
271's scalar/diagonal shrinkage and from 272's re-weighting of the FITTED operator's eigenvalues. It could fail; it
has not failed because diagonal shrinkage failed.

## 1. The unit

* **Stage A (held VIEWS)** — the unit is a **panel target in one activation state of GWCD4i**; the views are
  the (donor, guide) pseudobulks of that target. Fit on the 500 NON-panel Stim48hr targets, score on panel targets.
* **Stage B (held CONTEXTS)** — the unit is a **panel target on a PLOLO fold**; the estimator is fitted per donor line
  on all its bank rows EXCEPT the 300 panel targets, and every held number is a per-target `pds` / cosine on the
  held line's cache truth. Resampled unit in every interval: the TARGET (265's paired target bootstrap).
* **Stage C (Flex reference)** — the unit is one of the 85 h1-test targets, union OFF (273's instrument).

## 2. The estimator family (frozen; `est286.py` is the definition)

Moments, per donor line j, from two independent views A, B of every fitted row (each centred on its own mean over the
bank's rows, converted to the production row's units):

    Ŝ = sym(AᵀB)/n = (UᵀU − DᵀD/4)/n     N̂ = DᵀD/(4n)     Ĉ = (AᵀA + BᵀB)/(2n) = Ŝ + 2N̂
    U = (A+B)/2, D = A−B.  Ŝ projected to PSD by `psd_eigh`; the clipped energy is reported.

| arm | construction | rôle |
|---|---|---|
| `E0` `mean` | equal-weight mean of the donor rows | the incumbent, KA at 0.0e+00 |
| `E1` `e1:phi` | G = Ŝ⁺(Ŝ⁺ + N̂/J + εI)⁻¹, gains raised to ≥ φ in G's eigenbasis (L Q), φ ∈ {0.10, 0.25}; φ = 0 = un-floored | **the arm** |
| `E2` `e2:r` | Ŝ truncated to rank r ∈ {20, 100}, no floor | the construction the review warned against |
| `E3` `e3:c` | 062's whitener recipe on Ŝ + floor (c·ν_i, c ∈ {0, 1}), γ = 0.65, in the operator's place | the corrected-directions object |
| `E4` `e4` | 271's diagonal per-gene Wiener (`est271.wiener_own`) | already-failed baseline |
| `E5` `e5` | 272's eigenvalue floor in the FITTED operator's basis, c = 1 | already-failed baseline |
| `E6` `e6` | 272's cross-fitted operator (its construction B) where its arrays are on disk | already-failed baseline; SKIPPED and said so if absent |

Every row arm is scored twice: **`:nm`** (norm-matched to the mean row's norm — shape-only, gate **+0.0019**) and
**`:lit`** (literal, amplitude-changing, gate **+0.0089**). Both are priced.

Constants frozen in `est286.py`: `EPS_REL = 1e-4`, `PHIS = (0.0, 0.10, 0.25)`, `RANKS = (20, 100)`,
`E3_CS = (0.0, 1.0)`, `E3_LMAX_FLOOR = 1e-6`, `GAMMA = 0.65`, seeds `20260906`, eigensolver float32 `evd`.

## 3. The controls (frozen; each built by the SAME code path as its arm)

| control | construction | what a gain there means |
|---|---|---|
| `Cs` `:shufpair` | Ŝ **and** N̂ rebuilt from view a of row t paired with view b of row π(t), π a fixed derangement | the estimator is not reading a shared response |
| `Cw` `:wrong` | another donor line's Ŝ on the shared columns (own N̂ kept) | the covariance is generic, not response-specific |
| `Cr` `:rand` | energy-matched random PSD covariance (Wishart, trace-matched) | the gain is the shrinkage, not the directions |
| `Cm` `:mirror` | G' = 2I − G_φ | the sign of the effect is not read |

## 4. Known-answer rows (every screen; a failure aborts that screen)

1. **KA-M** `‖Ĉ − Ŝ − 2N̂‖_F/‖Ĉ‖_F = 0` at machine precision (the exact algebraic identity). The mandate's literal
   `‖Ĉ − Ŝ − N̂‖_F/‖Ĉ‖_F` is reported beside it and equals `‖N̂‖_F/‖Ĉ‖_F` — it is NOT expected to be 0.
2. **KA-E1id** E1 with N̂ = 0 and φ = 1 equals E0 (the identity gain).
3. **KA-0** `lib263.stage`'s rebuilt per-donor rows, float32, equal the cached `pre` bit for bit; **KA-C** this
   experiment's per-line production rows equal `st["rows"][j]` on the fold's targets.
4. **KA-271** 271's `ebrow` and `relw` through this experiment's plumbing reproduce 271's per-fold `pds` at 0.0e+00.
5. **KA-BR249** the BR249 retrodiction: `cd4svd5` vs `STACK` through 261's pricing recipe reads +0.00472 at 0.0e+00
   against 257 (and the board's +0.0042634).
6. **KA-shuf** the shuffled-pairing Ŝ_shuf's top eigenvalue vs the real Ŝ's — the ratio is reported (an expectation-zero
   signal covariance must have a far smaller leading eigenvalue).
7. **Stage A KA** E0 reproduces 280's Q2 `k = 2` row at 0.0e+00 (`work/q2_<cond>_keepDE.npz`, `full_k2`/`matched_k2`).
8. **Stage C KA** 273's `base` and `cd4svd5` rows reproduce `R273`'s `pds` at 0.0e+00.

## 5. The gate (frozen, Stage B)

An arm **CLEARS** iff ALL of:

1. six-member price, **PESSIMISTIC** scenario (READ-FIRST §3; 261's recipe, per leg the band end that LOWERS the
   total for that cell's sign) ≥ **+0.0019** (`:nm`, shape-only) or ≥ **+0.0089** (`:lit`, amplitude-changing);
2. Δ`pds` > 0 on **3/3** cancer folds AND the pooled paired-target bootstrap interval excludes 0;
3. leave-top-3-out retention ≥ **50 %**;
4. **every control below half the arm's Δ`pds`** (Cs, Cw, Cr, Cm) and the arm beats E4 and E5 by more than the
   control's bootstrap half-width.

`kolf2` is a FOURTH, labelled fold if ≥ 1 hour remains; it is never one of the 3/3.

## 6. Predictions (recorded before the measurement)

* P1 Ŝ is far from PSD: > 30 % of its eigenvalues negative, clipped energy > 0.10 (the shared response is a small
  part of the row).
* P2 The Wiener gains are small in most directions (median g < 0.25) — this is why a floor is needed at all.
* P3 If the review's mechanism is right, E1 with a floor beats E2 (low rank) on `pds`, and E2 loses `pds` relative to
  E0 (it removes the very directions `pds` is bought in).
* P4 A well-supported NEGATIVE is a valid result. The pre-registered failure mode: E1's Δ`pds` is at its `:shufpair`
  or `:rand` control, i.e. the gain (if any) is the shrinkage, not the learnt directions.

## 7. Scope

No build, no `vcc`, no `cell-eval2`, nothing under `submissions/`, `data/harness/` or `data/`; no download. Everything
written lives under `experiments/286-crossproduct-estimator/`. A Stage B arm that clears the gate goes back to the
coordinator for a preregistered build mandate; **production rendering is not this mandate.**
