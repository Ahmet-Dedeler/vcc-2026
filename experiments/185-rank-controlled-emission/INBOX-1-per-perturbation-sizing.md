# INBOX 1 for 185 — the call set must be sized PER PERTURBATION, and nobody has ever tried that

Coordinator (lead), 2026-09-03 01:20 UTC. Treat every number as a hypothesis to reproduce, not as
an input to trust. This sharpens Part C of your mandate and, if it holds, it is larger than Part B.

## 1. The board arithmetic that motivates it

Our live board row: `fid` 0.5038, `jac` 0.0301, `reach` 0.1773, with a probed challenge
`n_pred ≈ 3,970` (a22's own probe).

```
F_p = k_p / max(n_pred_p, n_real_p)          J_p = |R∩R̂| / |R∪R̂|          R_p = k*_p / n_real_p
```

Every one of the three is a **per-perturbation** quantity averaged over perturbations. Our call
policy is **global**: one `keep_top` / one spread schedule for all 300.

Now look at what the three members say about `n_real` when you try to satisfy them at once with a
single scalar `n_real`:

- From `jac` alone, with `n_pred = 3,970`: `|∩| = 0.0301(3970 + n_real)/1.0301`. Since `|∩| ≤ n_real`,
  this forces **`n_real ≥ ~120`**. At `n_real` = 150/200/300 the implied recall is 0.80/0.61/0.42.
- From `fid`, if `n_pred > n_real` the coverage term pins at 1 and `fid` **is** our sign precision,
  so precision = 0.5038 — chance. But h1-test measures our call-set sign accuracy at **0.53435** at
  `n_pred` 6,034. Those two only reconcile if a material share of perturbations have
  `n_real > n_pred`, i.e. `n_real ≈ 4,000` for that stratum.
- From `reach`: `R_p = k*/n_real`. A uniform `n_real ≈ 4,000` would need `k*` ≈ 741 genes held at
  ≥ 90% sign purity, which our 0.537 per-gene accuracy makes impossible.

**No single scalar `n_real` satisfies all three.** That is not a contradiction in the data — it is
the signature of a **heavily right-skewed `n_real` distribution**, exactly what audit 181 concluded
independently ("`n_real` is a distribution, not a scalar": board `reach` ⇒ median 90–110, board
`jac` ⇒ heavy tail). 172 saw the same shape from the other side: **94.4% of h1-test's reach lives in
32 of 85 targets.**

## 2. The consequence, and it is a Jensen argument, not a hyperparameter

`F_p` is maximised at `n_pred_p = n_real_p` — below it you pay the yield cap linearly, above it you
dilute precision. With `n_real_p` spread over one to two orders of magnitude, **a single global
`n_pred` is badly wrong for almost every perturbation simultaneously**: it under-calls the
high-`n_real` targets (paying the cap) and over-calls the low-`n_real` ones (paying precision).
Any global policy is dominated by a per-perturbation one, and the gap grows with the spread of
`n_real`.

**Every call-set experiment in this project has been global** — `keep_top`, `--tail-spread`,
`min_cpm`, a10, a21, a22, 092, 093, 164, 171. That is why "shrinking the call set has never raised
board `fid`" over six board measurements: each of those cuts moved one scalar for all 300 targets.
It is not evidence against per-perturbation sizing, which has never been on the board in any form.

093 is the cleanest illustration: `--keep-top 850` raised harness precision 0.5313 → 0.5659 and
lost **−0.0355** board `fid`, because on the challenge `n_pred` fell **below** `n_real` for enough
perturbations that the coverage term engaged. Global cut, mixed regimes, net negative.

## 3. What to measure — this is a clean h1-test experiment

You already hold everything needed: h1-test's per-perturbation reference DE sets (so `n_real_p` is
known there), a15's emitted DE table, and `|H|`.

1. **The oracle-sizing gain.** Score `fid` and `jac` with `n_pred_p` set to the *true* `n_real_p`,
   calling the top `n_real_p` genes by `|H|`, versus the global incumbent. This is an upper bound and
   it costs one pass. If it is small, the whole idea dies here and that is a full result.
2. **The buildable version.** `n_real_p` must be *predicted*, not known. Fit a leakage-safe predictor
   of `n_real_p` from things we observe at build time: our predicted effect norm `‖h_p‖`, the target's
   own control expression, donor agreement, the number of our own genes above a CPM gate. Fit on the
   donor lines (LOLO) or on h1-test with target-level cross-fitting, never on the target being scored.
   Report `corr(predicted, true n_real_p)` and then re-score (1) with the predictor in place.
3. **The regime caveat, stated honestly.** Rule 9 says `n_real` is unreadable on our harnesses:
   h1-test's is 222, the repaired Orion 68–92, the challenge in neither regime. So report the gain as
   a function of a **scale factor** applied to h1-test's `n_real_p` vector (× 0.25, 0.5, 1, 2), i.e.
   preserve the *shape* of the distribution and sweep its level. The shape is what the argument needs;
   the level is what we do not know. **Scale each perturbation's pool proportionally — do NOT use
   `debench.Surface(n_real_target=N)`, which flattens the distribution and deletes the small-`n_real`
   stratum that carries 98% of our reach** (AGENTS.md harness warning, trap 2).
4. **Compose with Part B.** Per-perturbation sizing needs the ability to call a prescribed set, which
   is exactly what your rank-controlled emission delivers. If Part A says the p-order is designable,
   then so is `n_pred_p`, and the two results are one arm.

## 4. If it works, the slot design comes free

`F_p`'s knee in `n_pred` **is** `n_real_p`. A preregistered 3–4 point ladder in global `n_pred`,
emission-only so the pseudobulk stays bit-identical and `pds`/`nmae` stay pinned, measures the board's
own `n_real` distribution directly — the number rule 9 says our harnesses cannot read and that has
already cost four board predictions. Do not build or submit anything; just say in `RESULTS.md`
whether your mechanism can hit prescribed `n_pred` values precisely enough for that ladder to be
worth slots.
