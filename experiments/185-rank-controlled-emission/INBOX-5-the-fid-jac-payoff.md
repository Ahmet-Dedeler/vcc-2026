# INBOX 5 for 185 — the payoff arithmetic for Part C, and why undershooting is the only real risk

Coordinator (lead), 2026-09-03 04:35 UTC. This is the number that justifies Part C. Reproduce it
before you believe it; every input is on the board or in the metric spec.

## 1. Solving our own call set from the board

Board: `fid` 0.5038, `jac` 0.0301, probed challenge `n_pred ≈ 3,970` (a22). Take `n_real ≈ 200`
(audit 181: board median 90–110 with a heavy tail; the answer below is not sensitive to this
between 150 and 300).

With `n_pred > n_real` the yield cap pins at 1, so **`fid` IS our sign precision over our own call
set: 0.5038 — chance.** And from `J = |∩| / (n_pred + n_real − |∩|)`:

```
0.0301 = |∩| / (3970 + 200 − |∩|)   ->   |∩| = 122.7
```

So of ~200 genes the reference calls, our 3,970 calls contain **123 (61% recall)** at a **precision
of 3.07%**. We are recovering most of the reference's set and burying it in 3,850 coin flips.

## 2. What the metric pays at `n_pred = n_real`

```
F_p = a(n) · min(1, n/n_real)        J_p = p(n)·n / (n + n_real − p(n)·n)
```
where `a(n)` is sign accuracy over the top `n` by our key and `p(n)` the fraction of them in the
reference's set. At `n = n_real` the cap is exactly 1, so **`fid` = `a(n_real)`** — nothing else.

| scenario at `n = n_real = 200` | `fid` raw | Δ overall | `jac` raw | Δ overall | total |
|---|---:|---:|---:|---:|---:|
| pessimistic `a` = 0.51, `p` = 0.08 | 0.510 | +0.0035 | 0.0417 | +0.0053 | **+0.0088** |
| central `a` = 0.56, `p` = 0.125 | 0.560 | +0.0315 | 0.0667 | +0.0166 | **+0.0481** |
| our measured top-decile-by-`|H|` accuracy 0.624 | 0.624 | +0.0674 | — | — | — |

**The asymmetry is the point.** `fid` at `n = n_real` equals the sign accuracy of our *top*
`n_real` genes, and our current 0.5038 is the accuracy of *all 3,970*. **Cutting to the head can
only lower `fid` if our ordering is anti-informative** — and 172 measured top-decile sign accuracy
on the DE pool at **0.624** against a pool mean of 0.537, so it is not.

## 3. The one real risk, and it has a board point

**Undershooting.** 093/v14 shipped `--keep-top 850` and lost **−0.0355** board `fid`: on the
challenge `n_pred` fell **below** `n_real` for enough perturbations that the coverage term engaged,
and the loss is linear in the shortfall. That is the entire historical case against a small call
set, and it is a case against a **global** cut, not against sizing.

`F_p` is per-perturbation and `max(n_pred, n_real)` is inside the average, so with `n_real` spread
over one to two orders of magnitude **a single global `n_pred` is simultaneously too small for the
big perturbations and too large for the small ones**. That is a Jensen loss, not a tuning error,
and it is why every global experiment in this family (a10, a21, a22, 092, 093, 164, 171) came back
flat or negative. **Per-perturbation sizing has never been on the board in any form.**

## 4. What to do with it

- Compute the table above on h1-test with the true per-perturbation `n_real_p`, once with **oracle**
  sizing and once with a **predicted** `n_real_p` (INBOX 1 §3 lists the build-time features). The
  oracle row bounds the idea; the predicted row is the arm.
- Size the safety margin explicitly. Because the penalty is one-sided — flat above `n_real_p`,
  linear below — **the expected-value optimum is deliberately ABOVE the predicted `n_real_p`**, by
  roughly the predictor's error sd. Report that margin; it is the whole design.
- Note this composes with, and does not depend on, Part B: `fid` and `jac` read *which* genes are
  called, `reach` reads the *order*. If Part A says only the set is controllable and not the order,
  Part C still stands on its own.

**Sizes for context:** the gap to first place is +0.0438; the whole ordering channel spans 0.0179
overall across the board. The central row above is **+0.0481** — larger than either. Treat that as
a reason for extra rigour, not for optimism: it is exactly the shape of result that turns out to be
a metric misreading, so re-derive `F_p` and `J_p` from
`docs/reference/vcc2026-metrics-BRIEF.md` §4 and §6 yourself before quoting any of it.

---

## 5. CORRECTION to §2, added 04:45 — `jac` can go DOWN, and the condition is sharp

I understated the downside on `jac`. Work the chance case: at `n_pred` 3,970 we hold `|∩|` = 123 of
~200. If the key carried **no** membership information, cutting to the top 200 would leave
`|∩| ≈ 123 × 200/3970 = 6.2`, so

```
J = 6.2 / (200 + 200 − 6.2) = 0.0157   against our current 0.0301 — WORSE.
```

**`jac` rises only if the key concentrates the reference's DE genes in the head by more than about
2×** (`p(200) > ~0.062`, i.e. `|∩|(200) > 12.4`). Below that the cut costs `jac`. So the `+0.0166`
in §2's central row is **conditional on a measurable property of the key that you must measure
first**, and it is a different property from the one `fid` needs:

| member | what it needs from the key | measurable on h1-test as |
|---|---|---|
| `fid` | **sign accuracy** in the head | `a(n)` — 172 measured 0.624 for the top `|H|` decile vs 0.537 pool mean, so this one is already known to hold |
| `jac` | **membership enrichment** in the head | `p(n) / (n_real/n_tested)` — the fold-enrichment of the reference's DE set. **Never measured here.** |

**Measure the enrichment curve before anything else in Part C.** It is one pass over h1-test's
reference DE sets and our own key, it costs nothing, and it decides whether Part C is a `fid`-only
arm (+0.003 to +0.032, one-sided upside because the yield cap pins at 1) or a `fid`+`jac` arm
(+0.009 to +0.048). Report `p(n)` and the fold-enrichment at `n` = 50, 100, 200, 400, 800, 1600 and
at the true `n_real_p` per perturbation.

Two notes that make the enrichment plausible but do not substitute for measuring it: our recall at
`n_pred` 3,970 is already 61%, so the reference's genes are not uniformly scattered through our call
set; and `jac`'s replicate anchor is only 0.3972, i.e. **real replicate experiments agree on DE
gene sets less than half the time**, so the achievable enrichment has a low ceiling even with a
perfect model.
