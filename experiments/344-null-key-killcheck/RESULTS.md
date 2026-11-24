# 344 — No construct-level null exists on the challenge controls; the empirical-null key is dead at reading 1

September 9, 2026, 22:58Z, 30 seconds total. **STOP_no_null**, the first registered stop rule. Truth-free;
challenge control files only. Registration `f54fcc46`; outputs `work/results.json`, `work/nullkey-{A,B,C}.npz`.

| context | 5-CPM universe | constructs with any BH-significant gene | median significant genes per construct | Spearman(key, CPM) | LOCO recall, key vs CPM (46-construct bootstrap 95%) |
|---|---:|---:|---:|---:|---|
| A | 9,929 | 16 / 46 | **0** | 0.097 | 0.386 vs 0.478, Δ −0.092 [−0.290, +0.066] |
| B | 9,626 | 8 / 46 | **0** | 0.045 | 0.248 vs 0.391, Δ −0.143 [−0.533, +0.181] |
| C | 10,124 | 10 / 46 | **0** | 0.033 | 0.133 vs 0.630, Δ −0.497 [−0.810, −0.180] |

Each non-targeting construct's 400 cells, tested by the scorer's own Wilcoxon against the pooled 18,400
controls exactly as a perturbation is tested, produce essentially no significant genes. The reference's null
constructs behave as clean replicates of the pool. So the premise of construction A in
`docs/sprint-2026-09-09/claude-emission-fable.jsonl` §A — that a gene-level "DE propensity" estimable from the
controls exists and ranks the 2,200-gene block better than CPM — is false on the actual challenge controls.
Where any significant set exists, ranking by the other constructs' mean |z| recalls it worse than CPM does,
resolved in C. The BH-frequency variant is worse still.

**Consequences.** No H1 four-arm test, no preregistration, no build. 192's "unbuildable" same-line LOTO DE-rate
key stays unbuildable: its information comes from real perturbations' DE calls, which the null constructs do
not carry. The negative is scoped to construct-vs-pool null structure; it says nothing about keys derived from
other information sources.

**Side finding, no verdict.** The scorer's per-perturbation DE calls are not inflated by construct-level
technical structure among the NTCs, so a real perturbation's significant set is signal plus sampling noise, not
signal plus clonal artefact.
