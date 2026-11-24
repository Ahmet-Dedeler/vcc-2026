# 423: rank-preserving carrier magnitudes (the mechanism behind cqawesome's pure-nmae steps). Mandate from the PI, 2026-09-23 01:30Z

**Current state.** Read `docs/FINDINGS.md` §0–§1, `experiments/416-de-carrier-menu/RESULTS.md`,
`experiments/417-bulk-sharpen/RESULTS.md` + `BOARD-RESULT-417c.md` and `experiments/422-carrier-sign-only/RESULTS.md`.

Our best entry, 417a (0.2435334), is:
- a bulk (pds/mse) set by ballast cells;
- DE-channel carriers (fid/reach/jac/nmae) from the rc_callset emitter (`src/vcc2026/rc_emit_k.py`), with the clean4q80
  sign prior.

Our nmae (0.9366 raw) is the worst member relative to the field (cqawesome 0.836).

On the public board, several of cqawesome's steps moved **only nmae**, with pds/mse/fid/jac identical to 4 decimals
(e.g. 09-18 02:42 → 06:27: nmae 0.8629 → 0.8423). That is only possible if per-cell VALUES change while every gene's
Wilcoxon rank sum against the reference controls, and each block's summed counts, stay fixed.
- The Wilcoxon statistic depends only on how many control cells each of our cells exceeds.
- A carrier cell already above every control value can be pushed higher at zero rank cost, which raises that gene's
  arithmetic-mean CPM and hence its LFC magnitude.
- Moves inside a cell's interval between adjacent control values are rank-free in both directions.

**What we know on H1 (PI, 416's base_s0, analytic).**
- Our carriers' |LFC| on the reference's DE genes is small: median |ℓ|/|r| ≈ 0.16–0.19.
- A global LFC scale k gives nmae 0.9521 → 0.9479 at k ≈ 1.3 (−0.004); k 2 is worse.
- `sign × const` is worse (0.9726), because |ℓ| carries sign reliability: the top quintile of |ℓ| has sign accuracy 0.61,
  the rest 0.54.
- For pure magnitude changes the board read the sign right and was about 2.3× more sensitive than H1 (410's weld: H1 nmae
  +0.024, board +0.054).
- For carrier DESIGN changes H1 mis-signed nmae twice (392, 417c). Your arm must be a pure magnitude change.

**Territory.** Nobody else is on the DE channel right now. 417 is building bulk arms. Do not touch the builder, vcc,
cell-eval2 or the bulk. Work in the H1 six-member harness, copying 416's `run416.py` stack into your directory (its
known answers: base_s0 reproduces 416's `metrics-base_s0.json` to ≤ 1e-12).
- (a) Implement a **rank-preserving magnitude post-pass** on the emitted carrier block. For each gene, scale the carriers'
  arithmetic-mean CPM LFC toward a target, using only moves that leave the gene's Wilcoxon U against the REFERENCE
  control cells unchanged: push top cells above the control max, and move within control-value intervals. The ballast
  must keep the block's per-gene totals exactly: design-mode for touched genes, i.e. the ballast gives back or absorbs
  exactly the counts moved.
  - Account for library-size side effects: adding counts to a cell lowers every other gene's CPM in it. Measure and
    bound the resulting rank changes on untouched genes.
- (b) Verify on H1 that pds, mse, fid, jac and reach are unchanged to realisation noise, and that per-gene U is identical
  for the touched genes. Then measure nmae for:
  - global magnitude factor k ∈ {0.8, 1.2, 1.4, 1.6, 2.0};
  - k applied only to genes in the top |ℓ| quantiles (q50, q80);
  - up-genes only vs down-genes only (feasibility differs).
- (c) Report the feasible fraction (what share of the target LFC change each gene reaches under the rank constraint),
  the nmae prize at 3 seeds with sd, and an emitter-level spec for the builder (where in `ballast_rows` / the writer the
  pass would go, and its known-answer test).

**Deliverable:** `experiments/423-rank-preserving-magnitude/RESULTS.md`:
- the known answers;
- the table (nmae mean ± sd over 3 seeds, the other five members' deltas, feasible fraction);
- the recommended rule, with its predicted board nmae range. Use H1 × 1…2.3 for pure-magnitude changes, and say why.

**Constraints.**
- H1 arms strictly serial, ≤ 8 GB RSS.
- Memory gate before each run: `memory_pressure` ≥ 40 % free and free+inactive ≥ 20 GiB. ~15 GB builds share this Mac.
- Write only in experiments/423-rank-preserving-magnitude/.
