# 196 — The emitted library depth is a free parameter nobody has ever swept, and it quantises the profile we deliver

You are a research subagent. Work ONLY in `experiments/196-emitted-depth/`.
No `vcc`, no `cell-eval2`, no submission build, no writes outside your directory.
Peak RSS ≤ 12 GB. `uv run python`.

## READ FIRST
`AGENTS.md` §2 (the hard numbers and the caps), §4, §6.1 (the format contract), §11.1 rules 1–13.
`docs/reference/vcc2026-metrics-BRIEF.md` §1–§3. Then:
- `src/vcc2026/donor_transfer.py:321` `emit_counts`, and `scripts/make_submission_pm_cd4.py`
  around lines 2327–2370 and 2836 (the call-site contract on `lib_sizes`).
- `experiments/043-template-emitter/RESULTS.md` (Hamilton rounding — the only place integer
  rounding has been thought about), `096-mse-floor/`, `171-reach-prefix/`,
  `185-rank-controlled-emission/RESULTS.md`, `191-rc-hybrid-build/RESULTS.md`.

## The observation

The quasi emitter converts a relative-expression profile into counts as

```
L      = int(median(lib_sizes))          # the CONTROL cells' own median library, ~20,000 UMI
totals = floor(n * L * profile_p)        # n = 400 cells
```

**`L` is inherited from the data. It has never been treated as a parameter, and nothing in the
scoring contract requires it to equal the reference's depth.** The scorer normalises every cell to
1e6 CPM before ranking (BRIEF §3) and normalises the group profile to `TS = 5e4` (§1), so absolute
depth cancels in every metric *except* through integer rounding and through our own emitted noise.

**And the rounding is coarse where it matters most.** At `L = 20,000` and `n = 400`, a gene at
1 CPM receives `floor(400·20000·1e-6) = 8` counts across 400 cells; at 5 CPM, 40. The DE gate admits
genes above **5 CPM in the control**, so a large share of the ~10,700 gated genes — the exact
population `nmae`, `fid`, `jac` and `reach` are scored on — is delivered at a few tens of counts.
**The profile the scorer reads is a coarsely rounded version of the profile we computed**, and the
rounding error is largest, in relative terms, precisely on the low-expression genes.

## The headroom, and the two constraints

- **Counts per cell ≤ 1,000,000** (§2). At ~20,000 we are **50× under**.
- **Stored entries in `.X` ≤ 4,750,000,000** (§2). Measured on the real artifacts by 191:
  `a15` **50.06%** of the cap, `rc_hybrid` **34.68%**. So roughly **2×** headroom for a15 and
  **2.9×** for rc before the format rejects the file. Deeper cells mean more nonzeros, so this is
  the binding constraint, not the per-cell cap.

## The competing effects, which is why this needs measuring rather than reasoning

- **Profile fidelity (up).** Larger `L` shrinks the relative rounding error of every gene's
  delivered total, so the delivered pseudobulk moves toward the intended one. `pds`, `mse` and
  `nmae` all read that pseudobulk.
- **`mse`'s jackknife (down, possibly).** `expr_mse_unbiased_capped_norm` subtracts
  `ρ · min(Ĉ, C_p)`, where `Ĉ` is the delete-one jackknife of **our own** emitted noise. Deeper
  cells mean less emitted sampling noise, so `Ĉ` falls and `claim` with it. 191 measured
  `Ĉ/C_p` = **1.094** for a15 and **2.786** for rc_hybrid, with rc at the cap on 48/48 rows. So
  a15 has little margin and rc has a lot — **the sign of the `mse` effect may differ between the
  two emitters, and that is itself worth knowing.**
- **DE power (up).** More counts per gene means more power, so `n_pred` rises. That is controllable
  — 185/191 established that the call set and the p-order are designable at fixed per-gene totals —
  but it must be held fixed or explicitly re-optimised, not left to drift.

## Work

### Phase 0 — known-answer rows (mandatory, abort on failure)
- Reproduce a15's h1-test row through the existing path: `reach` **0.18283**, `n_pred` **6035**,
  `jac` **0.06824**, `nmae` **0.96659**, and 096's `claim(a15 as shipped)` = **0.24581**, raw `mse`
  **0.88228**. Donors `k562_gwps` + `orion_hek293t` only on h1-test (`h1_2025` covers 85/85 and is
  pure leakage).
- Quantify the defect before trying to fix it: at production's `L`, what is the **relative error
  between the intended `profile_p` and the delivered `totals/(n·L)`**, per gene, by control-CPM
  decile? And what share of the gated pool sits below 100 delivered counts? If that error is
  negligible the experiment is over — report it and stop.

### Phase 1 — sweep `L`
Scale `L` by 1× (identity row), 2×, 4×, 8×, and at whatever multiple hits 90% of the nnz cap. For
each: delivered-vs-intended profile error, nnz as a share of the 4.75e9 cap, max counts per cell
against 1e6, `pds` (**PLOLO**, 162's instrument — `debench.Surface.pds()` must never size a `pds`
arm, 188 Phase 3), `mse` via 096's `claim` route on the challenge panel, and
`nmae`/`reach`/`fid`/`jac` on h1-test with a **significance-selected pool** (188: `|lfc|·√cpm` at
the reference's BH depth, median ~89 genes/row; a top-N-by-magnitude pool inverts `nmae`'s sign).

Do this for **both** emitters — the shipped quasi law and 191's `rc_hybrid` — because their `Ĉ/C_p`
margins differ 2.5× and the `mse` term may move in opposite directions.

### Phase 2 — the second, cheaper question in the same file
`fid` counts a called gene with **no direction** as a miss (`BRIEF` §4: *"A gene the submission
calls but leaves without a direction stays in `n_pred` and counts as a miss"*), and a gene is
adjudicable only when the **reference** gives it a defined non-zero lfc. Measure, on the challenge
contexts, **how many genes we call for which our own predicted lfc is exactly zero** (delivered
total equal to the control's after rounding). Every one is a guaranteed miss we are paying for. If
the count is material, say what it is worth and how to stop calling them.

### Phase 3 — falsifiers
1. **Format contract**: any arm exceeding the nnz cap or 1e6 counts/cell is dead on arrival —
   report the margin, do not assume it.
2. **Rule 12**: price all six. `pds`/`nmae` are *not* automatically pinned here — unlike an
   ordering arm, changing `L` changes the delivered totals themselves.
3. **The `mse` sign**: if `claim` falls out of saturation, `mse` degrades; check the branch
   explicitly rather than assuming the cap holds.

## Deliverable
`RESULTS.md`: the rounding-error table by expression decile, the `L` sweep on both emitters with all
six members, the zero-lfc call count, and a verdict. "The rounding is already negligible" is a full
result and closes a parameter that has sat unexamined since 2026-08-21.

## Territory
192 owns the call set. 193 owns the operator's eigenvalue map. 194 owns the donor tables. 195 owns
fitting the transform. You own the **emitted depth and integer rounding**.
