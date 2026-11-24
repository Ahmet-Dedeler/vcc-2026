# 232 — preregistration

Written **before** any Phase 1 or Phase 2 number was read. Phase 0's known-answer rows were run
first and are recorded in `work/phase0_232.log`; nothing below was conditioned on them beyond
"they passed".

## The claim under test

`AGENTS.md` §8.1d: a cosine gain earned by reshaping per-gene **magnitudes** is destroyed by
`nmae`; one earned by improving **signs** is not. The only untested construction in the consensus
is **how the donors' signs are combined**: production emits `sign(mean of donor effects)`, which is
a magnitude-weighted vote, not the donors' majority.

## Predictions, stated in advance

1. **Phase 1 terminal condition.** If the |t|-weighted **disagreement mass** (the share of the
   scored DE pool where the covering donors do not all agree in sign) is small, or if no
   implementable vote beats `sign(mean)` on that set, the experiment ends at Phase 1 and that is
   the complete result. No arm is built on a null.
2. **The arithmetic ceiling.** Where every covering donor agrees, every sign rule — mean, majority,
   any weighted vote, and the oracle — returns the same sign. So
   `prize_ceiling = mass_disagree × (1 − acc[sign(mean)] | disagreement)` bounds the |t|-weighted
   sign-accuracy gain of **any** sign-combination rule, including an oracle one. It is computed
   before any arm exists and it is the number that decides whether Phase 2 is worth running.
3. **If an arm is built and it works,** §8.1d predicts `nmae` moves the RIGHT way (down) — the
   falsifiable prediction this experiment exists to test. Efficiency `e = ΔM/ΔU` should be **above**
   the 0.5 break-even, not 228's 0.072.
4. **If an arm is built and `nmae` still rises,** §8.1d's mechanism is wrong as stated: the wall
   would not be magnitude-specific.

## Fixed in advance

- **Folds:** the three cancer PLOLO folds are the headline; `kolf2` is reported apart and is
  excluded from every headline number (pluripotent, `cos_agg` base 3–5× smaller, and it inflated a
  coordinator prediction on 2026-09-04 badly enough to force a band reissue).
- **Base:** the shipped `h1union` stack — three donors + PerturbMatch 50 + CD4 `f = 0.34`
  + H1 union `g = 0.28`, `alpha 0.982, gamma 0.65, lam 0`, a23 render for `pds`/`cos`.
- **Pool:** 188's `wald89` (`|lfc|·√cpm`, top 89/row) is primary; `lfc89` and `h89` carried.
- **Endpoint order:** `nmae` leads (§8.1d), with `M`, `U`, `e` alongside, then `pds`/cosine.
- **Arm form:** `|effect|` bit-identical to the incumbent's, asserted elementwise at runtime.
  Only signs move, and only on coordinates with **≥ 2 covering donors**; singletons are left alone.
- **Dose:** the vote's normalised confidence margin `tau ∈ {0, 0.1, 0.2, 0.3, 0.5, 0.7}`.
- **Seeds:** `20260904` for every random control.
- **Gate:** rule 13's shape-only `+0.0019`; noise floors per §8.1b.
- **Transfer:** PLOLO units only. `k_cos` is carried as the range `{1.000, 1.384}` and never as a
  point, per §8.1c's 7.84× spread finding. No PLOLO→board scalar conversion is quoted.

## Falsifiers, all mandatory, all fixed in advance

1. **Sign-shuffle null** — flip the same NUMBER of signs per target, drawn uniformly from the
   coordinates with ≥ 2 covering donors, seed `20260904`.
2. **Wrong-recipient null** — apply target `π(i)`'s flip set and vote signs to target `i` under a
   derangement. 184 died here (arm +0.0665, control +0.0667, null +0.0663).
3. **Leave-top-3-out retention, on cosine and `nmae` SEPARATELY** (222 §9: 14.6% vs 72.7% — the
   channels differ).
4. **Magnitude invariance** — `|effect|` elementwise identical, asserted, plus the delivered
   median row norm ratio reported.
5. **Coverage honesty** — a gene covered by one donor has no vote. The count of pool pairs with
   1 / 2 / 3 covering donors is reported and singletons are left alone.

Plus one control that is not a null: **`ctrl_mean`** imposes `sign(mean)` on the *identical*
coordinate set. It absorbs everything the coordinate selection does on its own, so the
`vote − ctrl_mean` difference is the vote and nothing else.
