# 233 — How much of our score SURVIVES a disjoint panel? The biggest unmeasured question for the prize.

You are a research subagent. Work in `experiments/233-october-retention/`. **No `vcc submit`, no
model name, no description, no production change.** You **may** re-fit the operator and re-extract
banks **into your own `work/`** — never over the shipped artifacts. Peak RSS ≤ 12 GB; **check
`memory_pressure` and `sysctl vm.swapusage` before any job over ~4 GB and hold if swap free < 2 GB**
— 229 is running 20 GiB builds and a job that exceeds memory here is killed **without a traceback**
(219 recorded this). `uv run python`. Three other agents are live (221, 226, 229, 232).

## The question

**On 2026-10-22 the 300-target panel is replaced by one DISJOINT from today's, and only the last
entry after that date is prize-eligible.** Everything we have measured — every arm, every gate, our
whole 0.2118310 — is measured on today's panel.

**Nobody has ever measured how much of it survives.** 199 rehearsed the transition **mechanically**
(timings, bit-identity, that the pipeline runs) and says so explicitly: *"No slot spent, no
`vcc submit`, no `cell-eval2`."* **Its synthetic-panel artifacts were built and never scored.**

**If a large fraction of our performance is panel-fitted, then every +0.003 we are chasing today is
misallocated, and we should be spending the next 48 days on transition robustness instead.** That is
the decision this experiment informs, and it is worth more than any arm currently in flight.

## Why it is feasible now
`DEF-RUNBOOK` §5 / 199 already provide **two synthetic 300-target panels, both asserted disjoint from
A/B/C**: `def-bundle-r` (realistic: K562 ∩ Orion ∩ CD4) and `def-bundle` (pessimistic: any donor).
**Reuse them. Do not invent a sampling model.** 221's cloud `cell-eval2` worker has just passed its
known-answer row (nine of ten metrics bit-exact, `nmae` to 3.3e-16) if you need heavy scoring.

## Work

### Phase 0 — known-answer rows, mandatory
1. **Re-fitting the operator for the A/B/C bundle must reproduce `e2f65ecf…`** (`DEF-RUNBOOK` §4).
2. **Re-extracting the CD4 banks for the A/B/C panel must reproduce the shipped banks byte-for-byte**
   on all 36 data files (done 2026-09-03; reproduce it).
3. **Your PLOLO baseline on today's panel must match the shipped stack at Δ = 0.0e+00** on the three
   cancer folds.
**If any fail, stop.** Everything downstream is a ratio against these.

### Phase 1 — the retention measurement
For **each** synthetic panel (realistic and pessimistic), and for today's panel as the control:
1. **Re-fit the operator for that panel.** This is not optional — `DEF-RUNBOOK` §0.1: the operator is
   a **function of the panel** and fails **silently** on a new one, with **98.8% of coordinates
   paired with the wrong gene while every check passes**. 212 priced that trap at **−0.0429 overall
   on `pds` alone**. **Assert the axis matches the panel, and demonstrate the failure mode
   deliberately once so your assertion is proven to fire.**
2. Re-extract the CD4 banks for that panel.
3. Run PLOLO on the three cancer folds, in the challenge frame, and report **all six members** plus
   aggregate cosine.
4. **Report the RETENTION RATIO per member**: new panel / today's panel.

### Phase 2 — attribute the loss
Decompose which components lose the most. `DEF-RUNBOOK` §5.-2's coverage table is the prior:
`k562_gwps` **269 → 150**, our shipped H1 sources **64 → 34** and **20 → 3**, `kolf2` and `feng_ipsc`
**→ 0**, `dld1_lowmoi` **67 → 73 (improves)**. **Say which of those coverage changes actually carries
the score change**, and whether the loss is in DIRECTION (`pds`/`mse`/`nmae`) or ORDERING
(`fid`/`reach`/`jac`) — §8.1's split, which now also carries a precision asymmetry (§8.1b).

### Phase 3 — the four sources that only exist in October
`nadig_hepg2`, `nadig_jurkat`, `rpe1`, `k562_essential` carry **0/300 of today's panel and ~30 each
of October's**. **On a synthetic panel they become priceable for the first time.** Price them. This
is the only opportunity to do so before the day itself, and 225 flagged their October value as
*"the single largest untested claim in this file."*

## Falsifiers
1. **Panel disjointness** — assert it, do not assume it.
2. **Operator-axis identity** — per Phase 1.1, with the deliberate failure demonstration.
3. **Two panels, not one.** If realistic and pessimistic disagree in direction, report both and draw
   no single conclusion.
4. **Coverage honesty** — report carrying rows per source per panel, and use the **carrying-row**
   median for any residual sizing. The harness mirror at
   `experiments/162-retrodicting-validator/plolo.py:363-369` silently returns `f_per_unit = 0.0`
   below 50% coverage, so a live source reads exactly zero and looks dead (223). **Low-coverage
   sources are exactly what a disjoint panel produces.**

## Deliverable
`RESULTS.md`: the Phase-0 rows, a **retention table per member per panel**, the attribution, the
four October-only sources priced, and a one-line verdict of the form *"we retain X% of our
direction and Y% of our ordering on a disjoint panel, and the loss is carried by Z."*

**Either answer reallocates the next 48 days.** High retention means keep chasing arms. Low retention
means the transition is the whole game and we should be hardening it now.

**Propose no slot. Build no submission.**
