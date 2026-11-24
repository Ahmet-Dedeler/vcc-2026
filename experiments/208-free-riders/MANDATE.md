# 208 — Compose the free riders, and adjudicate 200's `jac` convention with the board point we now hold

You are a research subagent. Work ONLY in `experiments/208-free-riders/`.
No `vcc`, no `cell-eval2`, no submission build, no writes outside your directory
(you may read `src/` and `scripts/`; you may NOT modify them). Peak RSS ≤ 12 GB. `uv run python`.

## READ FIRST
`AGENTS.md` §4, §6.1, §8.1, §8.2, §11.1 rules 1–13. Then:
1. `experiments/196-emitted-depth/RESULTS.md` — **the `floor` → Hamilton finding and its own
   verdict: "free rider, not a slot."**
2. `experiments/200-production-ablation/RESULTS.md` — the ten-component ablation, and specifically
   **the tier-vs-flat-1.00 trade it could not adjudicate**, handed to 192 with three seeds emitted.
3. `experiments/202-callset-rider/BOARD-RESULT.md` and
   `experiments/185-rank-controlled-emission/BOARD-RESULT-rc_hybrid.md` — **the two board points
   that make this experiment possible**, both landed today.
4. `experiments/192-predicted-callset/RESULTS.md` §4.2 — the board's own `n_real` distribution.

## Why now

Today's two submissions took us 0.1740865 → **0.1918638** (rank 11 of 643) entirely inside the
ORDERING channel, and both proved the same identity: **`pds` and `mse` reproduced to sixteen
significant figures** while the order and then the call set changed underneath them. That identity
is what makes a *composed* emission arm cheap: every component below is free on the three profile
members by construction.

**Two riders are sitting unclaimed, and one blocked question can now be answered.**

### Rider A — Hamilton rounding. Measured, +0.00120, vetoed alone.
`emit_counts` at `donor_transfer.py:380` does `totals = floor(n·L·p)`. Every gene is short ~0.5
counts; because the metrics read *ratios* and 0.5 is a *fixed count*, that is a **size tilt of
−0.01234 log2** across challenge A/B/C. Largest-remainder (Hamilton) rounding at the **same `L`**
removes **97%** of it (+0.00039 residual) and 196 measured:

| | `floor` | `hamilton` |
|---|---|---|
| PLOLO `pds`, 4 folds | 0.739016 | **0.742224** = **+0.001182 overall** — *above* the `m=∞` ceiling of 0.742137 |
| `mse` `claim` | 0.24648 saturated | 0.24653 **still saturated**, Δ +0.00000 |
| DE (native / wald89 / matched) | — | **+0.00237 / +0.00011 / +0.00237** (quasi); **+0.00337 / +0.00260 / +0.00338** (rc) |
| nnz | 50.06% of cap | 50.13% |

Mechanism: `pds` is hurt by the **correlated** part of the rounding error, not its energy — floor's
error adds 1.4–2.5× its own energy to the squared error, Hamilton's 0.88–1.25×. That is why halving
the energy captures 103% of the entire depth prize.

### Rider B — the tier-vs-flat-1.00 trade, which 200 could not call
200 measured that **`mse` does not need the tier**: flat spread 1.00 at φ4 holds `claim` **0.24653**
and raw `mse` **0.88228** *bit-equal to a15's*. So tier-vs-flat is a **free four-member trade**,
worth **−0.005580** under a fixed-`n_pred` `jac` convention or **+0.002585** under a scaled one.
200 could not adjudicate: *"the one board point (a21) favours fixed-NP but tests the wrong
direction."*

**We now hold a board point that tests the right direction.** `rc_callset` moved board `jac` raw
**0.03031 → 0.04043** (+0.01012, scaled −0.000253 → **+0.027335**) at `n_pred` 4,480 → **2,599**,
with `pds`/`mse`/`nmae` provably unmoved. That is a large, clean, single-variable `jac` move with a
known `n_pred` change — **exactly the observation the two conventions disagree about.** Use it.

## Work

### Phase 0 — known-answer rows (mandatory, abort on failure)
- Reproduce 196's Hamilton table: the challenge tilt **−0.01234 → +0.00039**, PLOLO `pds`
  **0.739016 → 0.742224**, and `claim` staying saturated at 0.24653.
- Reproduce 200's flat-1.00 row: `claim` **0.24653** and raw `mse` **0.88228**, bit-equal to a15's.
- Reproduce both of today's board rows from the artifacts on disk: `rc_hybrid`
  (`submissions/191-rc_hybrid.h5ad`, sha `e0b0677f…`) and `rc_callset`
  (`submissions/202-rc_callset.h5ad`, sha `2a894aa5…`) — their h1-test DE tables and their challenge
  `n_pred` (4,480 and 2,599).
- **Rule 11:** path and sha256 for every operator; the shipped one is `29,394 × 17,942`, sha
  `e2f65ecf`, and nine near-identically-named operators exist.

### Phase 1 — adjudicate the `jac` convention against the board
Write both conventions down explicitly as formulae, predict `rc_callset`'s board `jac` under each
from `rc_hybrid`'s artifact plus the known `n_pred` change, and compare to the realised
**+0.01012 raw**. Report which convention retrodicts and by how much. Then re-price 200's
tier-vs-flat trade under the surviving one. **If neither retrodicts, say so** — that is a finding
about the instrument, and it means the trade stays unadjudicated rather than being called on a
coin-flip.

### Phase 2 — compose, on top of `rc_callset` (the current incumbent)
Price, marginally to `rc_callset`, each rider alone and then composed:
- A: Hamilton rounding;
- B: flat 1.00 spread (only if Phase 1 says it is positive);
- A ∘ B.

**Do not assume additivity** — 178 closed the rider family precisely because *"gains compose
sub-additively, costs do not."* Measure the composition, and if you can only afford one composition
run, measure A ∘ B rather than inferring it.

Endpoints: PLOLO `pds`; **Δcos on PLOLO, which 204 measured as resolvable to ±0.20% relative
(≈±0.00006 overall)** and is the sharper direction instrument; `mse` through 096's `claim` route,
reporting whether `claim` stays saturated (a15 has only **5.6%** margin, so this is the binding
check); and `nmae`/`reach`/`fid`/`jac` on h1-test with a **significance-selected pool** (188:
`|lfc|·√cpm` at the reference's BH depth, median ~89 genes/row — a top-N-by-magnitude pool inverts
`nmae`'s sign). Donors on h1-test: `k562_gwps` + `orion_hek293t` only; `data/donors/h1_2025` covers
85/85 and is pure leakage.

### Phase 3 — the gate, and the honest board conversion
Rule 13's `reach` toll does **not** apply to a profile-preserving arm — refuted on the board by
`rc_hybrid` (see the amendment in §11.1 rule 13). Use the transfer factors **today's two board
points fix**, not the stale ones: `reach` **×1.290** measured on `rc_hybrid`, and for `jac` the
realised-vs-predicted ratio from `rc_callset` (**1.10×** on its own prescribed retrodiction).
Gate: **+0.0019** shape-only. Report the composed arm's predicted board Δ against `rc_callset`.

### Phase 4 — falsifiers
1. **Per-gene panel-total digest identity** against `submissions/202-rc_callset.h5ad` on 900/900
   blocks. Hamilton changes the *totals themselves* — so unlike the previous two arms this is
   **NOT** a pseudobulk-preserving change, and `pds`/`mse`/`nmae` are **NOT** pinned. Say so
   loudly, and price all six members without the exemption.
2. A **known-answer emission**: `floor` at the shipped policy must reproduce a15's / `rc_callset`'s
   DE table bit-identically through your code path, or your path is not theirs.
3. **nnz against the 4.75e9 cap** (`rc_callset` sat at 33.22%) and the **1e6 counts/cell** cap,
   which 196 found binds first at `m ≈ 45`.

## Deliverable
`RESULTS.md`: known-answer rows, the `jac` convention adjudication against the board point, the
composed six-member price marginal to `rc_callset`, the gate verdict, and — if it clears — the exact
builder invocation plus a preregistered per-member prediction. **Do not build. Do not submit.**

## Territory
207 owns the CD4 residual's dose and form. 199 owns D/E/F readiness. **You own emission-side riders
and rounding.** Do not touch the operator (193 closed it), the donor tables (194, 204), or the call
set's composition (192/202 — you may change the *spread policy* but not the union rule).
