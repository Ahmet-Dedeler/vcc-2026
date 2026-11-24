# 185 — Rank-controlled emission: the p-value order is a free variable, not a spread consequence

You are a research subagent. Work ONLY in `experiments/185-rank-controlled-emission/`.
No `vcc`, no `cell-eval2`, no submission build, no writes to `submissions/` or outside your
directory. Peak RSS ≤ 12 GB. `uv run python` always.

## READ FIRST
`AGENTS.md` §4 (scoring), §6.1 (format contract), §11.1 rules 1–13 **including the audit-181
amendments**. `docs/reference/vcc2026-metrics-BRIEF.md` §3–§6 (the exact DE table, fidelity, reach,
jaccard definitions — read the actual text, not a summary). Then:
- `experiments/171-reach-prefix/RESULTS.md` **in full** — it closed the spread-schedule family and
  its §12c mechanism is the premise of this experiment.
- `experiments/172-sign-confidence-key/RESULTS.md` — the key-quality half.
- `experiments/164-callset-artefact/`, `157-reach-order/`, `177-resample-emitter/`,
  `179-hy-kt8000-build/`, `181-audit-0902/RESULTS.md`.

## The premise 171 established, and the one it did not test

171 proved:
1. The reach ordering rule "our significant calls first" is a **mathematical no-op** (BH-adjusted p
   is monotone in p). Reach with the call set replaced by nothing or by everything is **0.18283
   either way**. Only the **p-value order** can move `reach`.
2. Under the shipped zero-inflated gamma-Poisson emitter, `|z|` is **U-shaped in the emitted spread
   `s`** with its minimum at the shipped 0.25, and **both branches are dominated by control
   expression rather than by our predicted effect**. So no setting of the per-gene spread vector
   delivers the ordering key we hold. The spread-schedule family is closed.

**What 171 did not test is the general case.** Its knob was one scalar per gene inside a fixed
parametric law. The emitted count matrix is a far larger object, and the metric constrains it much
less than the emitter does:

- `pds` and `mse` read the **group pseudobulk profile** = the per-gene total over the 400 cells.
- `nmae` and every `lfc` read the **arithmetic mean of per-cell CPM values**, which — when every
  emitted cell carries the **same library total `T`** — is exactly `1e6 · (gene total) / (400 · T)`.
- `fid`, `jac` and `reach` read only the **Wilcoxon rank-sum statistic** of our 400 per-cell CPM
  values against the reference's real control cells, gene by gene.

So at fixed per-cell total `T` and fixed per-gene total `M_g`, **the pseudobulk-reading members are
exactly pinned while the rank-reading members are free**, because a rank-sum statistic depends on
the *shape* of the per-cell count histogram and the mean does not. Concretely, for a gene with
`M_g` counts to distribute over 400 cells, putting one count in `M_g` cells and zero in the rest, or
`M_g` counts in one cell and zero in 399, give the **same mean** and **opposite rank shifts** against
a control group that is mostly zeros.

**Hypothesis under test: the emitted p-value order is a designable object.** If so, the ~0.05 raw
`reach` that 172 measured as the difference between our delivered order and our own `|H|` key is
deliverable, and the delivered order stops being a side effect of expression.

## Known-answer rows — all mandatory, abort if any fails (rule 10)
| gate | expected |
|---|---|
| a15 h1-test `reach`, native BH, fresh emission | **0.18283** |
| a15 median `n_pred` | **6035** |
| a8-flat h1-test `reach` / `n_pred` | **0.20291** / **3581** |
| `|H|` key reach on the same pool | **0.23533** |
| oracle (our own signs perfectly ordered) | **0.64553** |
| your Wilcoxon/BH implementation vs `experiments/153-de-bench/debench.py` | identical |

Reuse `experiments/171-reach-prefix/realize171.py`'s code path — same donors (095's h1-test
artifact: `k562_gwps` + `orion_hek293t`), same gate convention, so 157's board-validated `reach`
calibration still applies. **`data/donors/h1_2025` covers 85/85 h1-test targets and must never
enter.**

## Work

### Part A — controllability (the structural question)
For one perturbation on h1-test, take the shipped per-gene totals `M_g` and the real control cells.
Ask: **over legal allocations of `M_g` counts to 400 cells at fixed per-cell total `T`, what is the
achievable range of the Wilcoxon z for gene `g`?** Compute the reachable `z` interval per gene
(the extremes are the concentrated and the maximally-spread allocations; you can compute both
exactly). Report the distribution of that interval's width across genes, split by control expression
decile. **This is the experiment's core measurement**: if the reachable interval is wide for most
genes, the order is designable; if it collapses for the low-expression genes that dominate the pool,
the mechanism is dead and you should say so and stop.

Note the two coupling constraints and handle them honestly: counts are integers, and the per-cell
totals must stay equal to `T` (or at least be held fixed, since CPM normalisation is per cell). A
per-gene allocation changes the per-cell totals unless you repair them. Say exactly how you repair,
and prove the repair leaves the per-gene totals bit-identical (this is rule 12's non-vacuous form —
hash the delivered per-gene panel totals per target and require the digest to match a15's).

### Part B — realise a prescribed order
Build an emission law that targets a prescribed per-gene rank shift `z*_g = sign(H_g) · φ(|H_g|)`
for a monotone `φ`, subject to Part A's reachable interval (clip where infeasible). Emit, run the
real Wilcoxon/BH, and measure:
- **Spearman(delivered p-order, `|H|` order)** — the delivery fidelity. Incumbent order's
  correlation with `|H|` is what 171 measured as expression-dominated; report both.
- `reach`, `fid`, `jac`, `n_pred`, and the per-gene panel-total digest (must be a15's).
- The mean `P(1)` and the `P(k)` curve (a15's mean `P(1)` is **0.7317**; the best 171 arm 0.8171;
  the bar is `P₀ = 0.90`).

### Part C — the call-set policy, now that it is free
With `n_pred` a design variable and the called genes choosable, compute on h1-test the achievable
`fid` and `jac` as functions of `n_pred` and of the key used to choose. Remember
`F_p = k / max(n_pred, n_real)` and `J_p = |R∩R̂| / |R∪R̂|`, and that **`n_real` is a distribution,
not a scalar** — audit 181 puts the board's median at **90–110** with a heavy tail, roughly half
h1-test's 222. Report the optimum as a function of assumed `n_real`, not at one value.

### Part D — BOARD RETRODICTION (mandatory; audit 181's standing rule)
This family has board points. **You must retrodict them before your h1-test numbers may be read.**
- `a21` (2026-09-01, `resample_x 0.14`): board `reach` **−0.0185**, `nmae` −0.064, `fid` −0.029,
  `jac` −0.017, `n_pred` 4,420 → 1,141, `fid` 0.5034 → 0.4945.
- `a22`: the second emission-law board point.
- Rule 13 charges emission-law arms a **−0.0146 scaled `reach` toll**, which is exactly
  `mean(a21, a22)`.
- The h1-test → board `reach` delta shrink factor is **0.4277** (164 §6).
State plainly whether your arm is in a21/a22's class. It differs from them in one respect that may
matter: **a21 changed the delivered profile** (resampling from control cells moves the lfc), whereas
this arm holds the per-gene totals bit-identical. If you argue the toll does not apply, you must
support it with a measurement, not with the argument alone.

### Part E — six-member price
Price all six (rule 12). `pds` and `nmae` should be provably invariant — prove it with the digest,
do not assert it. `mse` is **not** invariant: the jackknife term `C_hat` reads our emitted noise, so
price it through 096's `claim` route (`claim(a15 as shipped)` = **0.24581**; a15's raw `mse`
**0.88228**; a15 sits **+0.01313 of claim headroom above the `mse` knee**, so `mse` is flat for any
arm above it — check whether yours still is). Then apply the composed slot gate: SHIP only if
predicted > **+0.0019** (shape-only) or > **+0.0089** (amplitude-changing).

## Reporting
`RESULTS.md`: known-answer table, Part A's reachable-interval distribution, Part B's delivery
fidelity and members, Part C's frontier, Part D's retrodiction, Part E's price, and a verdict.
A negative with a mechanism is a full result.

## Territory
184 owns Tahoe and cross-line transfer. 186 owns competitive intel. 187 owns amplitude / `mse` /
`nmae`. You own the emission law and the DE call set. Do not modify production code.
