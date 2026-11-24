# 227 — `nmae` is a SEPARATE AXIS from cosine, and the rank-2 team is 0.045 raw ahead of us on it. Find the mechanism.

You are a research subagent. Work in `experiments/227-nmae-axis/`. **No `vcc`, no submission build,
no production change. No `cell-eval2` while another agent is running one — check first.** Peak RSS
≤ 12 GB. `uv run python`. Four other agents are live (221, 224, 225, 226); write only inside your
own directory.

## The observation, from the 2026-09-04 04:33Z board (657 entries)

We are **rank 3 at 0.21183**. `Aginglab.com` is **rank 2 at 0.21288**, ahead by **+0.00105**.
Per member:

| member | us | Aginglab | their edge |
|---|---:|---:|---:|
| `pds` | +0.74537 | +0.71512 | −0.03025 |
| `mse` | +0.17957 | +0.17825 | −0.00132 |
| **`nmae`** | **+0.11074** | **+0.18398** | **+0.07323** |
| `fid` | −0.01648 | +0.00323 | +0.01970 |
| `reach` | +0.22294 | +0.20081 | −0.02213 |
| `jac` | +0.02884 | −0.00411 | −0.03295 |

**Their entire advantage is `nmae`.** Raw: **ours 0.93389, theirs 0.88902** — a gap of **0.04487
raw = +0.07323 scaled = +0.0122 overall**, which is essentially **the whole gap to first place**
(`cqawesome` leads us by +0.01248).

## Why this is not the direction problem we assumed

**Their implied transfer cosine is LOWER than ours: 0.09887 vs 0.09947** (197's inversion on raw
`mse`; their raw `mse` 0.81851 vs ours 0.81726). **A team with worse direction than us is 0.045 raw
better on `nmae`.** Direction cannot explain that.

**And across the board it is worse than uncorrelated:** over the 55 teams above overall 0.15,
`corr(implied cosine, scaled nmae)` = **−0.394**. The top of the `nmae` table is full of teams with
*catastrophic* cosine — `CellWM` +0.19847 at cosine **−1.648**, `Liaoyi Xu` +0.17673 at **−1.330**,
`AgingLab_BI` +0.16967 at **−2.640**. **`Aginglab.com` is the only team in the top eight by `nmae`
with a positive cosine at all.** So there is a way to score `nmae` that costs direction — most of
that table found it — and Aginglab appears to have found one that does **not**.

## What is already closed, and why it does NOT close this

- **187** closed **global and per-context amplitude**: the joint `mse`+`nmae` optimum is exactly
  `VCC_ALPHA = 0.982`, and within the `mse`-feasible window the collectable `nmae` gain is **0.0012
  raw**. **187 itself measured the gap to Aginglab at 0.054 raw** and concluded *"`nmae` is a
  direction problem, not an amplitude problem."* **That inference is now falsified by the board:**
  Aginglab's direction is worse than ours. 187's *measurement* stands — amplitude cannot do it — but
  its *diagnosis* was a false dichotomy. **There is a third thing.**
- **158** closed the **per-target scale** (re-priced 2026-09-04: reason 2 stands — the oracle
  per-target scale costs `nmae` −0.0245 overall, so a better key makes it worse).
- **211** closed global amplitude for the second time; **216** put the dose turnover at `f* = 1.018`.

**Do not re-open amplitude in any form.** Your variable is not a scale.

## Where the third thing probably lives — read §4.1 of `AGENTS.md` closely

`nmae` reads **log2FC of per-cell CPM *arithmetic* means**; `mse` reads the **log-pseudobulk**.
`AGENTS.md` §4.1 states it plainly: *"a different functional of the same counts than the
log-pseudobulk `mse` reads, so the α that optimises one is not the α that optimises the other."*

**Two functionals of the same counts means the emission distribution moves them independently.**
That is a *mechanism* lever, not an information lever — and mechanism is where this project has won
twice (`rc_hybrid` +0.0093 and `rc_callset` +0.0084, both emission-side, both at bit-identical or
near-identical pseudobulk).

Candidates, in the order I would try them:

1. **Emission spread / dispersion.** `emit_counts(mode="quasi", spread=...)` changes cell-to-cell
   variance while preserving gene totals (034). Jensen's inequality makes the arithmetic mean of CPM
   a function of that variance **at fixed pseudobulk** — so this can move `nmae` with `mse` and `pds`
   pinned. 034 found it useful for `fid`/`reach`; **nobody has swept it against `nmae`.**
2. **The DE pool / call set.** `nmae` is evaluated on a pool. 188 established the pool must be
   **significance-selected** (`|lfc|·√cpm` at the reference's BH depth, median ~89 genes/row) and
   that a top-N-by-magnitude pool **inverts `nmae`'s sign**. We changed the call set with
   `rc_callset` (`n_pred` 4,480 → 2,599, +0.0084 overall). **Check what that did to `nmae`** — it
   may have bought ordering and paid here.
3. **The log2FC estimator itself** — pseudocount, CPM normalisation, zero handling. A different but
   equally valid estimator of the same quantity can shift a ratio-of-means statistic materially.

## Work

### Phase 0 — known-answer rows, mandatory before anything
Reproduce **188's exact `nmae` identity** (`nmae = q|k−1| + (1−q)(k+1)`, **|t|-weighted** `q`, slopes
−(2q−1)≈−0.074 below and **+1** above `k=1`), **187's** joint optimum at α 0.982, and our own board
`nmae` (raw **0.93389**, scaled **+0.11074**) from a local instrument. If you cannot retrodict our
own board `nmae`, stop — you have no instrument and nothing after it counts.

### Phase 1 — is the axis real and is it reachable at fixed cosine?
Sweep each candidate above on PLOLO (three cancer folds, `kolf2` **separately** — it is pluripotent
and inflates), reporting **all six members every time**. The gate is:
**`nmae` up by ≥ 0.01 raw with `pds` and aggregate cosine pinned to within ±0.2% relative.**
An arm that buys `nmae` by losing cosine is the failure mode 90% of that leaderboard table is in —
**report it as a failure, not as a trade.**

### Phase 2 — retrodict Aginglab
Can any setting reproduce **their** six-tuple — raw `nmae` 0.88902 at implied cosine 0.09887, `fid`
+0.00323, `jac` −0.00411? Their **`jac` is negative and their `fid` is positive**, the opposite of
our profile (`jac` **world 1st** at +0.02884, `fid` −0.01648). **That shape is a clue about their
call set**: they appear to call *more, looser* than we do. Say what construction produces it.

### Phase 3 — falsifiers, all mandatory
1. **Cosine-pinned control** — every reported `nmae` gain must come with the aggregate cosine and
   `pds` deltas beside it, on the same run. No exceptions.
2. **Shuffle null** for any pool/call-set change.
3. **Leave-top-3-out retention** on `nmae` (222 §9: on `h1union`, `pds` retained 14.6% and cosine
   72.7% — the channels differ, so measure, do not assume).
4. **October**: the panel is **replaced 2026-10-22 with 300 targets disjoint from today's**, and only
   the last entry after that counts. State whether your mechanism is panel-dependent. An
   emission-side lever should be **panel-independent by construction** — say so explicitly if it is,
   because that is worth a great deal relative to the source ladder.

## Reporting
**Do not convert PLOLO to board points with a scalar** — the factor is **not identified** (1.384 on
`cd4f034`, 1.000 on `h1union`, matched functional). Report PLOLO units, name the fold set and the
functional. 221 owns the reconciliation.

`RESULTS.md`: the retrodiction rows, the swept axes with all six members, the falsifiers, a Phase-2
answer, and a one-line verdict. **"`nmae` is unreachable at fixed cosine, and here is the mechanism
that makes Aginglab's number possible" is a complete result.** So is finding the lever.

## Territory
221 owns calibration · 224 learned transfer · 225 sources-as-donors · 226 the operator's fitting
sample. **You own `nmae` and the emission's effect on it.** Do not price sources, do not fit transfer
models, do not touch the builder or `submissions/`.
