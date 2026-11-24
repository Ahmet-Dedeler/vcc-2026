# 213 — The transport law is the one load-bearing assumption nobody has ever tested

You are a research subagent. Work ONLY in `experiments/213-transport-law/`.
No `vcc`, no `cell-eval2`, no submission build, no writes outside your directory.
Peak RSS ≤ 12 GB. `uv run python`.

## READ FIRST
`AGENTS.md` §4.2, **§8.2**, §11.1 rules 1–13. Then:
1. `docs/CURRENT_STATE.md`, the **2026-09-04 00:00 UTC** entry — the table of everything now closed
   and why this is what is left.
2. `experiments/212-gamma-cosine/RESULTS.md` — **the gap is not inside the operator** (its entire
   cost is 1.3575× against a 1.655× gap; abolishing it still leaves 1.219× short).
3. `experiments/197-cosine-inversion/RESULTS.md` — our cosine **0.079249**, leaders **≥0.1312**,
   **+0.00365 overall per +0.01 of cosine**.
4. `experiments/206-cluster-method/RESULTS.md` §5 — they are on **public data with a method**.
5. `experiments/156-delivered-geometry/`, `026-lfc-space/`, `141-pan-coordinate-calibration/`,
   `104-context-heldout-credit-adapter/` — the nearest prior work. **Know exactly what each tested.**
6. `src/vcc2026/donor_transfer.py` and `scripts/make_submission_pm_cd4.py` —
   `profile_from_effect` and the staging chain.

## The argument

Every route to a better direction is now closed: donor lines (166, 205), cross-line adaptation (184,
on 50 Tahoe lines), learned transforms (195, at the cell-line boundary), the operator's eigenvalue
map (193), **the operator's exponent (212)**, donor re-estimation (194), donor weighting (204),
donor–donor correlation (209), the hierarchical estimator (203), and **global amplitude (211)**.
Only the CD4 dose remains, and it is built.

Meanwhile 206 established the leaders hold **1.66–1.76× our cosine on public data with a method**,
and 212 proved it is **not** in the operator.

**One assumption has been load-bearing since v1 and has never been varied.** The donor effect is

```
delta_donor = log1p(TS · p_pert) − log1p(TS · p_ctrl)          (in the DONOR's frame)
```

and it is applied to the recipient by **addition in that same space**, then rendered by
`profile_from_effect`. Every experiment in this project has varied what goes *into* that law — which
donors, what weights, what operator, what magnitude, what dose — and **none has varied the law
itself.** 200's ablation removed ten *components*; the transport was in all ten arms.

`log1p(TS·p)` is **linear below ~20 CPM and logarithmic above** (169 §6.4 measured exactly this, and
used it to explain why a power transform of `h` concentrates on high-*expression* rather than
large-*fold-change* genes). So the shipped law is neither a pure log-fold-change nor a pure
difference: it is a hybrid whose behaviour changes across the expression range, and **the recipient's
basal enters differently under each alternative.**

That matters because transport is exactly where the recipient's identity should enter, and five
experiments have failed to make the recipient's basal profile useful **as a feature** (160, 161,
174, 175, 104). None of them asked whether it is already entering **wrongly**.

## Work

### Phase 0 — known-answer rows (mandatory, abort on failure)
- Reproduce PLOLO's a15 `pds` **0.795784** and a23 **0.789121**, and the a15 → a23 mass-fix price
  **+0.004014** scaled (162; 204 and 209 both reproduce these exactly).
- Reproduce our shipped transfer cosine **0.079249** by 197's route, and 212's forward measurement
  that the **CD4 dose at f = 0.34 raises PLOLO cosine +13.83%**.
- Reproduce 169 §6.4's characterisation of the effect space (linear below ~20 CPM, logarithmic
  above; `|h|` ~2.4× larger on a 100-CPM gene than a 10-CPM one for the same fold change).
- **Rule 11:** path and sha256 for every operator. Production is **29,394 × 17,942**, sha
  **`e2f65ecf`**; nine near-identically-named operators exist.

### Phase 1 — enumerate the law, and measure what each one assumes
Write down the shipped law explicitly as a map
`(donor p_ctrl, donor p_pert, recipient p_ctrl) → recipient p_pert`, then enumerate the
alternatives. At minimum:
1. **shipped**: additive in `log1p(TS·p)` (the identity row);
2. **multiplicative / pure log-fold-change**: `p_rec_pert = p_rec_ctrl · (p_don_pert / p_don_ctrl)`;
3. **additive in raw CPM**;
4. **a variance-stabilising space** (`sqrt`, Anscombe, or `log(p + c)` with `c` fitted once);
5. **rank/quantile transport**: map the donor's within-gene quantile shift onto the recipient;
6. a **one-parameter family interpolating (1) and (2)** — e.g. `log(p + c)` with `c` sweeping from
   `1/TS` (log-fold-change) to large (additive) — so the shipped law is an interior point and the
   result is a curve, not a bake-off.

**Family (6) is the important one.** A bake-off of five named laws invites a lucky winner; a
one-parameter family with the incumbent inside it gives a maximiser and a curvature, and the
falsifier is whether the optimum is distinguishable from the incumbent.

### Phase 2 — score, with COSINE primary
For each law, on **PLOLO** (162's four folds, challenge frame):
- **Δcos — primary.** 204 measured this instrument as resolvable to **±0.20% relative
  ≈ ±0.00006 overall**; 212 used it for the same purpose. Report the interval.
- `pds` — **`debench.Surface.pds()` must never size a `pds` arm** (188 Phase 3).
- **|t|-weighted `q_w`** (187 §5).
- Report the cosine and `pds` maximisers **separately**. That divergence is now a documented pattern
  in this project — it is what 062 got wrong on the CD4 dose (207 recovered +0.0044) and what 170
  and 195 got wrong on gamma (212 confirmed it did not pay there).

**The operator interacts and you must handle it.** The whitener was fitted on effects in the shipped
space; a different transport changes what the operator is whitening. Report each law **both** with
the shipped operator and with an operator re-fitted in that law's own space, and say which
comparison you consider fair. 212's finding — the operator trades cosine for retrieval at 96.2% of
the board's rate — should be re-checked under any law that wins, because that ratio is what makes
gamma unpayable and it need not be law-invariant.

### Phase 3 — six members and the gate
Rule 12. `mse` through 096's `claim` route; note 211's standing rule: **`claim` is at its ceiling in
the one context that reads it (context A, `q_A` = 1.000000), so no arm may be re-priced on amplitude
grounds.** `nmae`/`reach`/`fid`/`jac` on h1-test with a **significance-selected pool** (188:
`|lfc|·√cpm` at the reference's BH depth, median ~89 genes/row; a top-N-by-magnitude pool inverts
`nmae`'s sign). Donors on h1-test: `k562_gwps` + `orion_hek293t` only.
Rule 13: profile-class, **no `reach` toll**. Gate **+0.0019** shape-only, **+0.0089** if amplitude
moves — a transport change will move `A`, so expect the higher bar and re-solve α (211: α\* = 0.9808
on `rc_callset`, and 0.9879 on `cd4f034`).

### Phase 4 — falsifiers
1. **The incumbent must be an interior point of family (6)** and reproduce the shipped numbers
   exactly at its own parameter. If it does not, your family does not contain the shipped law.
2. **Permutation control**: transport a *deranged* recipient basal (pair each donor effect with the
   wrong context's control). Any law whose gain survives that is measuring amplitude, not transport.
   209's random control was the maximiser over all 36 of its arms — that is the outcome to watch for.
3. **Leak-free fold**: report all four PLOLO folds including `kolf2`; never quote an absolute from
   the Orion pair, whose truth carries the long-gene artefact.
4. **The `mse` identity**: any law that changes the delivered profile must be priced through 096,
   not asserted exempt. `rc_hybrid` and `rc_callset` both reproduced `mse` to sixteen significant
   figures *because* their pseudobulk was identical; yours will not be.

## Deliverable
`RESULTS.md`: known-answer rows, the family and what each member assumes, the sweep with **cosine
primary** and both maximisers, the operator-interaction table, the six-member price at the re-solved
α, all four falsifiers, and a verdict in overall units.

**"The shipped transport is at the optimum of the family" is a full result** — it would close the
last load-bearing assumption in the stack, and it is worth knowing on a curve rather than on the
fact that nobody has looked.

## Territory
210 owns the CD4 dose (built, submitting at the reset). 211 closed amplitude, 212 closed gamma —
**do not reopen either.** You own the transport law and the render.
