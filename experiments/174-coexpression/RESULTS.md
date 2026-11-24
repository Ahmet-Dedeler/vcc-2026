# 174 — recipient target↔gene co-expression as a sign predictor: **NULL**, and the mechanism is the shared programme

Date: **2026-09-02**, coordinator-run, ~1 min compute. h1-test, 18,400 control cells, log-CPM,
Pearson `r(t, g)` across cells for the 85 targets × 10,779 DE-universe genes; truth = the reference
DE table (`164/de164.a15.npz`). Scripts `coex174.py`, `coex174b.py`. No build, no submission.

## Hypothesis
The 70% of response energy that is line-specific (161) is not predictable from per-gene *marginal*
control features (161: R² 1.7%). A (target, gene) **pairwise** feature from the recipient's own
cells — if `t` and `g` covary cell-to-cell in *this* line, `t→g` regulation is active in *this*
line — is information no donor carries and nobody had tested.

## Result

| key | sign acc on the reference pool | on the weak-donor half (`\|cons\|` < median) |
|---|---:|---:|
| `H` (shipped) | 0.5373 | 0.5076 |
| `cons` (un-whitened) | 0.5518 | 0.5131 |
| **`−sign(r)`** | **0.5158** | **0.5193** |
| `+sign(r)` | 0.4842 | — |

Weak but real, in the biologically expected direction (activation → knockdown lowers the target's
partners), and it *beats* `H` on the half where the donors are uninformative. `H`'s accuracy is
0.552 where `−sign(r)` agrees with it and 0.522 where it disagrees.

**Then condition on cell state by removing top principal components of the control cells:**

| PCs removed | `−sign(r)` acc | on weak half | `\|r\|` median |
|---:|---:|---:|---:|
| 0 | 0.5158 | 0.5193 | 0.0256 |
| **2** | **0.4885** | **0.4855** | 0.0067 |
| 5 | 0.4884 | 0.4861 | 0.0056 |
| 10 | 0.4965 | 0.4981 | 0.0051 |
| 40 | 0.4970 | 0.4993 | 0.0049 |

**Two principal components carry the entire signal.** Once the dominant cell-state axes are removed
the residual co-expression is noise (`|r|` median 0.005) and marginally *anti*-informative.

## Mechanism
In control cells, target and gene covary through the dominant shared programme (proliferation /
stress). A proliferation-associated target and a proliferation-associated gene are positively
correlated; knocking the target down slows the programme and the gene falls. `−sign(r)` is right
for exactly those pairs — **because both ride the same axis, not because `t` regulates `g`.** This is
the hubness signal that spectral whitening removes for `pds` (AGENTS.md §4.2), re-derived from the
recipient's own cells. It carries no target-specific information, and adding it to the consensus
(`sign(cons_std − λ·r_std)`, λ = 0.5) reads **0.5428 — worse than `cons` alone at 0.5518.**

## What this closes
Recipient-cell pairwise co-expression as a source of line-specific perturbation signal, at
Perturb-seq control-cell noise levels. With 161 (marginals) this closes the recipient's control
cells as a source of the line-specific component at both the per-gene and the pairwise level. The
only remaining recipient signal not yet tested at the pair level would need the *perturbed* cells,
which the challenge does not give us.

`reach(key=|r|)` = 0.1737 (below production); `|H|·(1+agree)` = 0.2368 vs `|H|` 0.2353 — noise.
