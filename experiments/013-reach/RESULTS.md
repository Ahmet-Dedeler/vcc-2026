# 013 — Gene ranking for `reach`, and the per-gene variance channel

**Date:** 2026-08-21 · LODO, 259 perturbations, top-200 genes by |true effect| standing in for
the reference's significant set, purity threshold 0.9 exactly as the spec defines.

## Why `reach` is worth attacking

Widest span of any metric — baseline 0.072, replicate 0.968 — and we sit at **0.042 scaled**,
against Vivai's 0.090. It is also the only metric that is purely a **ranking** question:
magnitude does not enter, only the order in which we stand behind our own sign calls.

Today that order is an accident. It falls out of whatever p-values the emitted counts happen to
produce, which we have never chosen deliberately.

## Result

| ranking signal | macro reach |
|---|---|
| **√(\|A\|·\|B\|) × sign-agreement** | **0.0419** |
| \|effect\| × agreement | 0.0367 |
| sign-agreement alone | 0.0367 |
| \|effect\| (what we effectively do now) | 0.0337 |
| min(\|A\|,\|B\|) | 0.0254 |
| random | 0.0067 |

**+24.2% over the magnitude ordering**, consistent across all three folds. The winner is the
geometric mean of the two donors' magnitudes, zeroed where they disagree on sign — it rewards
genes that are *both* strongly and *consistently* predicted, which is exactly the property
`reach` is testing.

Rough value: our raw reach 0.110 → ~0.136, scaled 0.042 → ~0.071, about **+0.005 overall**.
Comparable to the entire multi-donor gain, on a metric nobody has looked at.

## The part that matters more than the number

To use this I have to make the *scorer's* ordering follow the geometric signal. Per the spec,
`reach` orders genes by "significant calls first, then ascending adjusted p-value, then p-value,
then descending |lfc|". Those p-values come from the Wilcoxon on our emitted counts — and with
the even emitter, variance is near the integer floor, so the p-value ordering is essentially
monotone in |predicted effect|.

The obvious way to reorder is to reweight per-gene effect magnitudes. **That is exactly the
per-gene reweighting that lost in experiments 004 and 006**, because it distorts direction and
costs `pds`.

But there is a second channel. **A gene's p-value depends on its mean effect *and* on its
cell-to-cell variance, and we control both independently.** Emitting a gene with more variance
pushes it later in the ranking without touching its mean. So:

- genes we trust (high √(\|A\|·\|B\|), donors agree) → emit at the variance floor → significant,
  ranked first
- genes we do not trust → emit with inflated variance → weaker p-value, ranked later

The **mean profile is untouched**, so `pds`, `mse` and `nmae` should be unchanged, while `reach`
and possibly `fid` improve. Per-gene emission variance is a control channel that is orthogonal to
everything we have been tuning.

This is speculative until tested end to end — the count-level harness is the only thing that can
confirm it, and the interaction with `jac` (more variance → fewer significant calls → smaller
predicted set) is a real risk. But it is the first idea I have that improves one metric through a
degree of freedom the others cannot see.

## Caveat on the harness

LODO cannot compute `reach` exactly: I substituted "top-200 by |true effect|" for "genes the
reference called significant", and ranked donor effects directly rather than emitted-count
p-values. The *ordering* of candidate signals should be robust to both approximations, but the
absolute 0.0419 is not comparable to a leaderboard number.
