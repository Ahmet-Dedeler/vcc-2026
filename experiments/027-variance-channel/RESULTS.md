# 027 — The per-gene emission variance channel, and why the first version failed

**Date:** 2026-08-21 · real `cell-eval2` on `h1-fast`, whitened model (γ = 0.85, equal
donors, α = 1.08, CPM gate 5). Incumbent even emitter: `reach` 0.1519, `fid` 0.5108,
`nmae` 0.9534, `jac` 0.0748.

## The idea, and the piece experiment 023 supplied

Experiment 013 proposed it and could not test it. A gene's Wilcoxon p-value depends on its
mean effect **and** on its cell-to-cell spread, and we set those independently. So genes we
trust can be emitted at the variance floor — significant, and first in `reach`'s ordering —
while genes we do not get more spread and fall down the ordering, **with the mean profile
untouched**, so `pds`, `mse` and `nmae` cannot move.

Experiment 023 supplied the missing half: after whitening, `|predicted effect|` is a *worse*
per-gene confidence signal, while the un-whitened cross-donor geometric signal
(√ of the donors' magnitudes, zeroed on sign disagreement) holds up. The two channels want
**different objects**, so there is no conflict — whitened effects set the mean, un-whitened
agreement sets the variance.

## Version 1 — all the way to Poisson — is a wash

`mode="mixed"`: trusted genes get the even allocation, the rest get gamma-Poisson.

| tight k | fid | reach | nmae | jac |
|---|---|---|---|---|
| all (incumbent) | 0.5108 | 0.1519 | 0.9534 | 0.0748 |
| 1,000 | 0.5022 | **0.1994** | 0.9582 | 0.0680 |
| 3,000 | 0.5063 | 0.1827 | 0.9514 | 0.0680 |

`reach` +31%, and `fid` and `jac` take it straight back (≈ −0.002 overall net). The reason is
that full Poisson noise does not *reorder* a gene, it **removes it from the called set**.
`jac`'s recall collapses (exp 014: at our precision recall dominates Jaccard) and `fid` loses
its yield. What was wanted was a dial, and this was a switch.

## Version 2 — sub-Poisson, with the gene total preserved

`mode="quasi"`: split each gene's panel total `T` into an evenly-allocated part and a randomly
scattered part, with `spread` the fraction scattered. Per-cell variance is then ~`spread·T/n`
against Poisson's `T/n`, so **`spread` is literally the fraction of Poisson noise**, and the
gene's total is preserved — the pseudobulk profile is unchanged by construction, not by
approximation.

At k = 1,000, three seeds each:

| spread | reach | fid | nmae | jac |
|---|---|---|---|---|
| 0 (incumbent) | 0.1519 | 0.5108 | 0.9534 | 0.0748 |
| 0.10 | 0.1603 | 0.5123 | 0.9539 | 0.0728 |
| **0.15** | **0.1764 ± 0.013** | **0.5128** | 0.9528 | 0.0714 |
| 0.20 | 0.1751 | 0.5127 | 0.9534 | 0.0696 |

`reach` **+0.0245**, `fid` **+0.002** (slightly *better*, not worse), `nmae` unchanged as the
construction guarantees, `jac` −0.0034. On the official scale: reach +0.027, fid +0.007, jac
−0.009, nmae ~0 → **≈ +0.004 overall**.

## Seed replication: the gain is 3.4 sigma, and all the noise is on one side

The even emitter is nearly deterministic — `reach` 0.1519 / 0.1521 / 0.1511 over three seeds,
sd 0.0005 — because the only randomness in it is which control library sizes get drawn. The
quasi emitter adds real sampling, and its spread is 0.1621 / 0.1820 / 0.1850, sd 0.0125.

So the comparison is 0.1764 ± 0.0072 (SE) against a fixed 0.1517: **+0.0247, 3.4 sigma**.
`fid` is up in all three replicates (0.5133 / 0.5125 / 0.5127 vs 0.5108) and `nmae` is flat in
all three, which is what the construction guarantees.

Official scale: `reach` +0.028, `fid` +0.007, `jac` −0.009, `nmae` ~0 → **+0.004 overall**.

## Why it is not in tonight's submission

1. **The response is non-monotone at the low end.** `spread = 0.05` scored `reach` 0.1329,
   *below* doing nothing, on a single seed. A knob that dips below baseline before rising is
   either genuinely peaked or noisy, and one harness cannot tell which.
2. **`h1-fast` is the wrong biology for this.** It is 40 perturbations of an embryonic line on
   the 2025 common-essential-enriched panel — the exact gene class the 2026 panel excludes,
   and the confidence signal being tested is cross-donor agreement on *those* genes.

The somatic VCC-panel harness built alongside this (`data/harness/orion-hct116.h5ad`, 132 real
panel targets in HCT116 with 18,400 controls) is what should settle it.

## What is banked regardless

`emit_counts(mode="quasi", spread=...)` is a **magnitude-preserving** control channel: it
moves `reach` and `fid` while provably leaving the pseudobulk — and therefore `pds`, `mse` and
`nmae` — alone. That is the first degree of freedom found in this project that one metric can
see and the others cannot.
