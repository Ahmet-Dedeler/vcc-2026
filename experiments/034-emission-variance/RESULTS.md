# 034 — Global emission variance, and a rule for which harness to believe

**Date:** 2026-08-22 · `mode="quasi"` with a **global** `spread` (every gene gets the same
fraction of Poisson noise), on top of the shipped v6 configuration. Run on both harnesses,
which turned out to be the point.

## The two harnesses disagree, sharply

| spread | **somatic HCT116** (132 VCC-panel targets, **150** reference cells) | **h1-fast** (40 targets, 2025 panel, **400** reference cells) |
|---|---|---|
| | fid / avg6 | fid / avg6 |
| 0 (even) | 0.522 / 0.1042 | 0.5135 / 0.1823 |
| 0.35 | 0.536 / 0.1120 | **0.5171 / 0.1892** |
| 0.7 | **0.554** / 0.1213 | 0.5032 / 0.1814 |
| 1.0 | 0.554 / 0.1142 | 0.4970 / 0.1615 |
| alpha 0.5, spread 0.7 | **0.582** / 0.1229 | **0.3716** / 0.0861 |

On the somatic harness `fid` climbs to 0.595 and looks like the biggest unclaimed pool on the
board — **no team has positive `fid`**, best is −0.005 scaled. On `h1-fast` the same
configuration **collapses** `fid` to 0.372.

## Which one is right, and why

`fid = k / max(n_pred, n_real)`. Adding emission variance reduces `n_pred`, because genes with
no real effect stop clearing significance. Whether that helps or hurts depends entirely on
`n_real`, **which is set by the reference's cell count**:

- **150 reference cells** → the reference's own Wilcoxon is underpowered, `n_real` is small,
  a reduced `n_pred` still clears the yield cap, and `fid` becomes pure precision — which
  improves, because the surviving calls are the ones with real signal.
- **400 reference cells** → `n_real` is large, `n_pred` falls below it, and the yield term
  `min(1, n_pred/n_real)` bites. At alpha 0.5 with spread 0.7 it bites hard.

**The competition uses 400.** So `h1-fast` is the trustworthy harness on this question and the
somatic `fid` result is an artefact of Orion's ~154 cells per perturbation.

## The rule this establishes

Neither harness is "the good one". They are wrong about different things, and the split is
predictable from how they were built:

| question | trust | because |
|---|---|---|
| `fid` yield, `reach` depth, anything through `n_real` | **h1-fast** | 400 reference cells, matching the competition |
| which genes carry signal — `pds`, `jac`, the CPM gate | **somatic HCT116** | actual VCC panel in a somatic line |

That also resolves the CPM-gate disagreement (exp 031 §5): gate 20 beat gate 5 on the somatic
harness and lost on `h1-fast`. It is a gene-population question, so the somatic answer stands.

## What actually survives — and it is cross-validated

At **spread 0.35**, both harnesses agree, despite disagreeing on everything else:

| | somatic | h1-fast |
|---|---|---|
| avg6, spread 0 → 0.35 | 0.1042 → 0.1120 | 0.1823 → 0.1892 |
| **delta** | **+0.0078** | **+0.0069** |

Same size, from different biology and different cell counts. On `h1-fast` the gain is carried
by `reach` (0.1669 → 0.1967, **+18%**) with `fid` up slightly (+0.004), `pds` up (+0.007) and
`jac` down (−0.005).

Pushing further is where they diverge: spread 0.7 and 1.0 keep buying `reach` (up to 0.2198)
and start losing `fid` and `jac`. **0.35 is the joint optimum, and it is the only setting both
harnesses endorse.**

Estimated board value: `reach` 0.1584 → ~0.187 (+0.038 scaled), `fid` +0.013 scaled, `jac`
−0.014 scaled → **≈ +0.007 overall**.


---

# 034b — the exact-total emitter, and the pinned optimum

## A bug in the first version, and why it mattered

The first `quasi` scattered the random part with a plain Poisson draw, which preserves each
gene's panel total only **in expectation**. The residual lands directly in the pseudobulk —
which is exactly what `pds`, `mse` and `nmae` read. Measured on `h1-fast` at spread 0.25, that
alone moved `pds` between **0.8404 and 0.8583** across two seeds: it was putting our best
metric (best on the whole board) at risk to buy `reach`.

Repairing each gene's total exactly — scatter, then add or remove the residual count-by-count
— makes the emitted pseudobulk **bit-identical** to the even emitter's. Verified directly:
max absolute gene-total error 0 at every spread, while the variance-to-mean ratio rises 0.284
→ 0.468 → 0.777 for spread 0 → 0.25 → 0.7.

## The pinned sweep, `h1-fast`, multiple seeds

| spread | pds | nmae | reach | fid | jac | **avg6** |
|---|---|---|---|---|---|---|
| 0.00 | 0.8506 / 0.8513 | 0.9450 | 0.1669 | 0.5135 | 0.0730 | 0.1824 |
| 0.15 | 0.8506 | 0.9447 | 0.1829 | 0.5169 | 0.0710 | 0.1863 |
| 0.20 | 0.8506 / 0.8513 | 0.9451 | 0.1874 | 0.5182 | 0.0704 | **0.1875** |
| 0.25 | 0.8494–0.8506 | 0.9449 | 0.1901 | 0.5180 | 0.0694 | **0.1873** |
| 0.30 | 0.8506 / 0.8519 | 0.9452 | 0.1969 | 0.5171 | 0.0684 | **0.1879** |
| 0.40 | 0.8500 | 0.9449 | 0.1942 | 0.5151 | 0.0663 | 0.1849 |
| 0.60 | 0.8513 | 0.9449 | 0.1951 | 0.5003 | 0.0614 | 0.1752 |

**`pds` is 0.8494–0.8519 across every spread and seed — the same band the even emitter
occupies (0.8506–0.8513). `nmae` moves in the fourth decimal.** The construction holds.

Broad plateau at **spread 0.20–0.30**, avg6 0.1873–0.1879 against 0.1824 at spread 0.
`reach` +14% to +18%, `fid` +0.004, `jac` −0.004. Beyond 0.4, `reach` saturates while `fid`
and `jac` keep falling.

**Chosen: spread 0.25**, the middle of the plateau.

Board translation: `reach` 0.1584 → ~0.180 (+0.031 scaled), `fid` 0.5016 → ~0.506
(+0.012 scaled), `jac` −0.001 → **≈ +0.007 overall**, with `pds`, `mse` and `nmae` provably
unchanged. The downside is bounded to `jac`, which is worth 0.002 scaled to us in total.

---

# 034c — confidence-targeted spread adds nothing over a uniform one

Experiment 013 proposed the variance channel as a *targeted* device: trusted genes at the
variance floor, untrusted ones noisier, with trust coming from cross-donor agreement. The
exact-total emitter makes that testable against a plain uniform spread.

| config | reach | fid | jac | avg6 |
|---|---|---|---|---|
| no spread | 0.1669 | 0.5135 | 0.0730 | 0.1824 |
| **uniform spread 0.25** | 0.1901 | 0.5180 | 0.0694 | **0.1873** |
| targeted, k = 500 trusted | 0.1864 | 0.5183 | 0.0695 | 0.1874 |
| targeted, k = 1500 | 0.1919 | 0.5176 | 0.0694 | 0.1878 |
| targeted, k = 4000 | 0.1893 | 0.5176 | 0.0697 | 0.1874 |
| targeted, k = 1500, spread 0.5 | 0.2041 | 0.5123 | 0.0635 | 0.1843 |

Every targeted variant lands inside the noise of the uniform one. **The cross-donor confidence
ranking contributes nothing.**

That is worth stating because it relocates the mechanism. The gain is not "we know which genes
to trust" — it is that **a uniform variance level acts as a signal-to-noise filter on our own
DE calls**. Genes with a real predicted effect keep clearing significance; genes riding on the
even emitter's mean-versus-median artefact stop. Our confidence signal is no better at that
than the predicted effect size already is.

Consistent with experiment 022: per-gene confidence schemes keep failing to beat uniform
treatment, now including the one place a two-channel design seemed to have room.

**Ships the simpler thing:** a single global `spread`, no confidence computation.

---

# 034d — the CPM gate, resolved at the shipped operating point

Experiment 031 §5 left this open: gate 20 beat gate 5 on the somatic harness (+0.0032 avg6)
and lost on `h1-fast` (−0.0115), and both were measured without the emitter change. Re-run at
the actual operating point (spread 0.25), two seeds:

| | pds | nmae | reach | fid | jac | avg6 |
|---|---|---|---|---|---|---|
| **gate 5** (v6/v8) | 0.8506 | **0.9449** | 0.1911 | 0.5181 | 0.0695 | **0.1877** |
| gate 20 | 0.8526 | 0.9613 | 0.1703 | 0.5127 | 0.0695 | 0.1770 |

**Gate 5 wins by 0.0107**, and the whole difference is `nmae` (0.9449 vs 0.9613) plus `reach`.

## Why the two harnesses split, and which one governs

Gate 20 *is* slightly better for `pds` on both harnesses (+0.002 here, +0.010 on the somatic
one) — that part is a genuine gene-population effect and the somatic harness is right about it.

But the `nmae` penalty is a **cell-count** effect, and it dwarfs the `pds` gain. `nmae` is
scored over the genes the reference calls DE. With 400 reference cells the Wilcoxon has enough
power to call many 5–20 CPM genes; gating at 20 means predicting *no change* on every one of
them, and each contributes a full-magnitude miss. With the somatic harness's 150 cells the
reference finds fewer of those genes, so the penalty barely registers — `nmae` differs by
0.0025 there against 0.0164 here.

**The competition has 400.** Per the rule in §034: anything routed through `n_real` belongs to
`h1-fast`. **Gate 5 stands**, which is what v6 shipped and what v8 keeps.

## gamma, and why this harness cannot settle it

| gamma | 0.65 | 0.85 | 1.05 |
|---|---|---|---|
| avg6 | 0.1900 | 0.1877 | 0.1824 |
| pds | 0.8577 | 0.8506 | 0.8468 |

`h1-fast` prefers 0.65, and the entire difference is `pds` — the one metric this harness is
documented to be unreliable on, because `pds` saturates on a 40-perturbation panel. LODO, which
computes `pds` exactly on 259 real VCC-panel perturbations, puts the optimum at 0.85 (0.8078 vs
0.8043). **LODO governs `pds`**; the somatic sweep at 132 panel perturbations is running as a
third opinion.

---

# 034e — what the emitted cells actually look like

The mechanism claims a minimum-variance emitter makes nearly every gene significant because 400
near-identical CPM values sit at the control distribution's mean, which for skewed counts is
above its median. Measured directly on the submission artefacts, for three genes spanning the
tested range (context A, 400 emitted cells against 2,000 reference control cells):

| gene | source | mean | median | **sd** | **distinct values** |
|---|---|---|---|---|---|
| 60 CPM | reference | 62.78 | 50.14 | 67.84 | 1,308 |
| | v6 (even) | 61.02 | 48.34 | **21.27** | **3** |
| | v8 (quasi 0.25) | 61.02 | 48.40 | 29.72 | 267 |
| 400 CPM | reference | 399.43 | 383.36 | 183.61 | 1,939 |
| | v6 (even) | 374.24 | 386.68 | **21.14** | **3** |
| | v8 (quasi 0.25) | 374.24 | 383.83 | 70.96 | 350 |

**The even emitter produces three distinct CPM values per gene across 400 cells**, against the
reference's ~1,900. A rank-sum test between a 3-point distribution and a 1,900-point one will
find a difference almost regardless of the means — which is exactly why `n_pred` covers the
whole gene gate and `fid` collapses to chance.

Two things confirm the construction rather than just the story. The **means are identical to
two decimals** (374.24 for both v6 and v8), so the pseudobulk identity holds cell-by-cell, not
merely in aggregate. And v6's **median sits above** v8's and the reference's (386.68 vs 383.83
vs 383.36), which is the mean-versus-median offset the mechanism predicts.

## The optimum is deliberately *less* realistic than the data

v8's sd is 71 against the reference's 184 — 3.4x closer than v6 but still well under-dispersed.
Matching the reference exactly would mean `spread` ≈ 1.0, and the sweep is clear that this is
worse (avg6 0.1752 against 0.1877 at 0.25).

So the scoring optimum is **not** the realistic value. Full realism costs statistical power, and
the DE metrics reward calling genes; minimum variance costs precision, and they punish calling
everything. The optimum sits at about a quarter of the way from one to the other, and there is
no reason it should have landed on realism.
