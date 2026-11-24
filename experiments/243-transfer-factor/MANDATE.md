# 243 — Reconcile the PLOLO→board cosine transfer factor: **1.384 on `cd4f034`, 1.000 on `h1union`**

You are a research subagent. Work ONLY in `experiments/243-transfer-factor/`. Read-only everywhere
else. **No `vcc`, no build, no `cell-eval2`, no production change, no slot proposal.** `uv run python`.
Peak RSS ≤ 10 GB. **MEMORY GATE before any PLOLO fold load:** `memory_pressure` ≥ 40% free AND
`vm_stat` free+inactive ≥ 20 GiB; 241 is running heavy jobs, so run everything strictly serially and
record peak RSS with `/usr/bin/time -l`.

## Why this exists
Every board prediction for a direction arm passes through one scalar: board relative cosine gain ÷
PLOLO relative cosine gain. Two arms of the same family, both scored, give **1.384** (`cd4f034`,
216 §4.4: PLOLO +17.98% → board +24.89%) and **1.000** (`h1union`, 222: PLOLO +10.3065% on the
three cancer folds → board +10.3045%). `docs/CURRENT_STATE.md` ("Do not size an arm through a
scalar PLOLO→board factor — it is not identified") and 222's `BOARD-RESULT.md` §4 name three
candidate explanations and test none. 221 measured h1-test↔board factors on the real scorer but
states (its §9) that this PLOLO question is a different instrument and did not resolve it. **Until
it is resolved, no arm's predicted band means anything, and tomorrow's slot reading (237) leans on
it.** Resolve it from artifacts already on disk.

## READ FIRST
1. `experiments/222-h1-union-build/BOARD-RESULT.md` §3–§4 — the discrepancy and the three hypotheses.
2. `experiments/216-cd4-hard/RESULTS.md` §4.4 and §6.0 — how 1.384 was computed, which folds, which functional.
3. `experiments/222-h1-union-build/RESULTS.md` and `work/repro222.log` — the +10.3065% / +7.6803% rows.
4. `experiments/197-cosine-inversion/RESULTS.md` — the board-side "cosine" is INVERTED from raw `mse`
   through a one-sided identity; read exactly what is assumed about amplitude.
5. `experiments/215-leader-decomposition/RESULTS.md` §1.1 — `pds` as the amplitude-free readout (0.6479).
6. `experiments/219-reprice-closures/RESULTS.md` — the a14→a15 calibration (1.26) and the functional-mismatch finding.
7. `AGENTS.md` §8.4 (as corrected), §11.1 rules 10–13, §14. `docs/BOARD-LEDGER.md`.

## Known-answer rows (mandatory; abort on failure)
- KA-1: reproduce 216 §4.4's PLOLO **+17.98%** for `f` 0.10 → 0.34 from 216's own artifacts
  (`work/op216.json` / `opsix216.json` or the fold caches), on the fold set and functional 216 used.
- KA-2: reproduce 222's **+10.3065%** (`cos_agg`, three cancer folds) and **+7.6803%** (macro) from
  `experiments/222-h1-union-build/work/repro222.json`.
- KA-3: reproduce 197's inversion on `rc_hybrid` (raw `mse` 0.8738045 → cosine 0.079248) and on
  `cd4f034` (raw 0.83635 → 0.09027) and `h1union` (raw 0.81726 → 0.09947) with 197's constants.

## Work
### Phase 1 — put both arms on the SAME instrument, every way it can be defined
Build the matrix **arm × functional × fold set**: functional ∈ {`cos_agg` (energy-weighted),
macro (per-target mean), whatever 216 used if different}; fold set ∈ {4 folds, 3 cancer folds
(no `kolf2`)}; arm ∈ {`cd4f034` = `f` 0.10→0.34 with no union; `h1union` = union g 0.28 on top of
`f` 0.34}. Use 162's PLOLO instrument imported verbatim (`experiments/162-retrodicting-validator/
plolo.py`) and 218/222's staging for the union; **never re-derive a number that exists in an
artifact — read it, and say where from.** Report the PLOLO relative gain in every cell and the
resulting factor against the board's +24.89% and +10.3045%.

### Phase 2 — the board side is not a measurement, it is an inversion
197's identity gives a cosine FLOOR from raw `mse` under an amplitude assumption. Write the exact
decomposition `raw = r² − 2rc + 1` (197's form) and evaluate, for each of the two arms, how much of
the "board cosine change" could be an amplitude change instead: use the arms' own metas
(`effects.whitening.after_all_median_norm`, `post_add_global_preserve_scale` 0.920359 → 0.909357,
`alpha` 0.982 fixed) and 215's amplitude-free `pds` readout (Δraw `pds`: `cd4f034` −0.00365,
`h1union` +0.00344) as the second equation. State whether the two arms' board cosines are
comparable objects at all, and by how much the inversion's slack differs between a dose step and a
new-source step.

### Phase 3 — the rule
Give the single convention (functional, fold set, board-side reading) under which ONE factor
retrodicts all three board direction points — a14→a15 (219's 1.26), `cd4f034`, `h1union` — each
within ±15%. If none exists, say so, and give the smallest model that does (e.g. factor as a
function of the fraction of panel rows the arm touches, or of stack position), with the falsifier
the next board point would fire. Then re-price 237's `ARM6` (229 §pricing) and 228's per-gene arm
under that rule, so tomorrow's slot reading has a number to be compared against.

## Falsifiers, pre-committed
1. KA-1..3 reproduce to the reported digits, else stop.
2. The rule retrodicts all three board points within ±15%, or the RESULTS says "not identified"
   and gives the two-parameter form.
3. Every number cites the artifact path it was read from.

## Deliverable
`RESULTS.md`: the matrix, the amplitude-slack analysis with numbers, the rule (or the honest "not
identified"), the re-priced 237/228 bands, and a ≤ 10-line summary for `AGENTS.md` that the lead
will paste (do not edit `AGENTS.md` yourself). Then stop.
