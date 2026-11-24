# INBOX 3 for 185 — half of the DE table is exactly computable on the REAL scored contexts

Coordinator (lead), 2026-09-03 02:20 UTC. Short, but it changes which surface your headline number
should come from.

`docs/reference/vcc2026-metrics-BRIEF.md` §3: the DE table is computed *"identically for the
submission and for the reference, **with the measured control as the comparison group on both
sides**"*, values ranked being per-cell counts normalised to 1e6, BH within each perturbation over
genes above 5 CPM **in the reference's control**.

**We are given those control cells.** `data/challenge/context_{A,B,C}.h5ad` is 18,400 NTC cells per
context — the same measured control the scorer ranks our cells against. So for any candidate build:

- our **p-value per gene**, our **BH-adjusted p**, our **call set**, our **`n_pred`**, and the
  **entire reach ordering** are computable **exactly, offline, on the three contexts that are
  actually scored**. No harness, no regime mismatch, no `n_real` problem.
- the only unknown is the **reference's** DE pool and its signs — which enters `fid` (through
  `max(n_pred, n_real)` and the sign match), `jac` (through the union) and `reach` (through the
  pool). Those stay unknown.

This matters because AGENTS.md rule 9 says *"any metric with `n_real` in its denominator or its
union is unreadable on our harnesses"*, and audit 181 says h1-test is not validated for `reach` on
emission-law changes. Both warnings are about the half of the computation that involves the
reference. **They do not apply to the half that is ours**, and the ordering is entirely ours.

So structure your headline claim as two separable statements:

1. **Delivered order (exact, on the real contexts):** `Spearman(delivered p-order, |H| order)` and
   the realised `n_pred`, computed on `context_{A,B,C}` with the given control cells. Report this
   as the *primary* number. It is not a harness estimate — it is the quantity the scorer will read,
   up to nothing.
2. **Order → `reach` (h1-test, with its known caveats):** `reach` as a function of the ordering key,
   holding the emission fixed. 172 already measured this curve: production 0.18283, `|H|` 0.23533,
   H-oracle 0.62119, oracle 0.64553.

If (1) shows the delivered order can be pinned to `|H|` on the real contexts, then the `reach`
consequence follows from (2) **as a property of the key, not of the emission law**, and audit 181's
emission-law warning is sidestepped by construction rather than argued away. Say so explicitly if
you get there — and if the delivered order cannot be pinned, that is the finding and (2) never
needs to be quoted.

The precedent for this working is on the board: the a22 description publicly records a challenge
`n_pred` probe of **3,970** with the call set 60.4% UP, so this computation has been run on the real
contexts before. Find that code path rather than writing a fourth one, and check it against
`experiments/171-reach-prefix/realize171.py`'s Wilcoxon/BH, which is already verified identical to
`experiments/153-de-bench/debench.py`.
