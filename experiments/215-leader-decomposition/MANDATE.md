# 215 — Decompose the leader's 27-submission history into CHANNELS, with the instruments we now hold

You are an intel/analysis subagent. Work ONLY in `experiments/215-leader-decomposition/`.
Read-only on the repo. No `vcc submit`, no builds, no `cell-eval2`.
**Post nothing anywhere, contact nobody, create no accounts, join nothing, and put none of our
method into any query, form, issue or public field** (`AGENTS.md` §14). Read-only research only.
Peak RSS ≤ 8 GB. `uv run python`.

## READ FIRST
`AGENTS.md` §8, **§8.2**, §14, and the rule-13 table in §11.1. Then:
1. `experiments/206-cluster-method/RESULTS.md` — the previous pass. **Do not repeat it.** It already
   established that the four `pm-` accounts are one team (Prima Mente), that the board's cosine floor
   is one-sided and vacuous above raw `mse` 1.024, that `cqawesome` is untraceable after two sweeps,
   and that the field's implied `k_eff` of 10–14 exceeds the ~11 CRISPR-perturbed human lines that
   exist in total.
2. `experiments/197-cosine-inversion/RESULTS.md` — **the inversion, and its exact identifiability
   condition.** This is your main instrument.
3. `experiments/212-gamma-cosine/RESULTS.md` — the measured exchange rates
   `d(overall mse*)/d(cos) = +0.3655` and `d(overall pds)/d(cos) = −0.3517` on the gamma axis.
4. `docs/BOARD-LEDGER.md` and `experiments/186-intel-0903/board/trajectories.md` — our own points
   and the saved per-team histories.
5. `experiments/185-…/BOARD-RESULT-rc_hybrid.md`, `202-…/BOARD-RESULT.md` — our two ordering points,
   which calibrate what an ordering-only move looks like on the board.

## Why this pass is different from 206's

206 asked *who* and *whether it is data or method*, and answered both. It did **not** ask **what each
individual arm did**, and at the time we lacked the instruments to.

**We now hold measured exchange rates for every channel:**

| channel | signature on the board | rate |
|---|---|---|
| **direction** (cosine) | `mse` and `nmae` move together, `pds` moves with them | `mse` +0.3655/cos, and 197's +0.00365 overall per +0.01 cos |
| **ordering** (p-order, call set) | `reach`/`fid`/`jac` move; **`pds`, `mse`, `nmae` EXACTLY unchanged** | our own `rc_hybrid` (+0.0549 `reach`, `pds`/`mse` identical to 16 s.f.) and `rc_callset` (+0.0276 `jac`) |
| **amplitude** | `mse` moves against `nmae`; `pds` scale-invariant, so **`pds` unchanged** | 187/211: `pds` exactly invariant to a global α |
| **geometry** (operator) | `pds` up, cosine **down**, at 96.2% cancellation | 212 |

**Those four signatures are nearly orthogonal on the six-tuple.** So a team's arm-to-arm delta can
be *decomposed* rather than guessed at — and `cqawesome` has **27 submissions** with a
`pds` trajectory running **0.627 → 0.8703**, the largest climb on the board.

## Work

### Phase 0 — known-answer rows (mandatory, abort on failure)
Your decomposition must reproduce **our own** arms, whose channel we know for certain:
- `rc_hybrid` (ordering-only): must decompose to **Δcosine ≈ 0, Δamplitude ≈ 0**, all move in
  ordering. `pds` and `mse` moved by **exactly 0** to sixteen significant figures.
- `rc_callset` (ordering-only, call set): same, with the move in `jac`/`reach`.
- `a11` (α 1.0796 → 0.838, pure amplitude): must decompose to **Δcosine ≈ 0**.
- `a15` → the CD4 dose step (a14 → a15, direction): must show a **positive Δcosine**; 197 measured
  **+6.8%** for that pair and 212 measured **+13.83%** forward at `f = 0.34`.
- `a19` (operator ridge, geometry): `pds` and `mse` should move against each other.
If your decomposition misclassifies any of these five, it is not calibrated and nothing downstream
is readable. **Report the confusion matrix over these five before anything else.**

### Phase 1 — decompose every leader trajectory, arm by arm
For `cqawesome` (27), `Aginglab.com` (11), the four Prima Mente accounts, `Garrrbage`, `Slurm Shadi`
and `Haystack`, take every consecutive pair in the saved snapshots and classify the move into the
four channels with a magnitude and an uncertainty. Use 186's saved snapshots plus a fresh board
pull; state how many snapshots each pair is resolved by, and **do not classify a pair that only one
snapshot separates** unless the move exceeds the seed noise.

**Report, per team, a channel budget**: how much of their total climb came from direction, ordering,
amplitude and geometry. That table is the deliverable.

### Phase 2 — isolate the direction arms
For every arm anywhere on the board that decomposes to a **large positive Δcosine with ordering and
amplitude near zero**, report: the team, the arm, the date, the model-name pair, the implied Δcosine
in absolute units, and **what our own history has at that date** (there may be a public event —
a dataset release, a paper, a repo commit — that coincides).

197's identifiability condition binds and you must respect it: the inversion is one-sided,
`cos ≥ (1.023611 − raw)/2.074495`, **vacuous above raw `mse` 1.02361**, so **558 of 612 teams carried
nothing** in 197's sweep. For arms in the vacuous region report a **bound**, never a point, and say
which of your classifications rest on bounds.

### Phase 3 — the question that decides how we spend October
206 concluded **method/public data** on the grounds that three mutually independent pipelines sit
above cosine 0.11, two of them lone individuals, and that `k_eff` 10–14 exceeds the ~11 lines that
exist. **Test that conclusion against the trajectories.**

- If the leaders' cosine arrived in **one or two large steps**, it is an acquisition or a method
  switch, and it has a date we can look at.
- If it accumulated in **many small steps**, it is tuning of something they already had, and the
  channel budget tells us which knob.
- If their cosine has been **flat and high since their earliest resolvable arm**, they started with
  something we do not have.

`Haystack` is the cleanest natural experiment on the board and 206 flagged it: `reach` **0.2050 at
its third submission**, when `pds` was 0.686 — worse than our first entry — and it barely moved while
`pds` climbed 0.12 and `mse` fell 7×. Decompose that trajectory too; a team whose ordering was
structural from day one and whose direction climbed steadily is the shape that separates the two
hypotheses.

## Deliverable
`RESULTS.md`: the five-arm calibration confusion matrix, the per-team channel budgets, the isolated
direction arms with dates and bounds, and a verdict on Phase 3 with the evidence that decides it.
Plus an explicit list of what you could not resolve and why.

**A clean "their cosine was there from the start and never moved" is a full result** — it would mean
the remaining 48 days should go to D/E/F robustness and to the measured gains we already hold, not
to chasing a direction we have now closed eleven ways.

## Territory
213 owns the transport law. 214 owns component removal and rounding. **You own the board's history
and its decomposition.** Do not re-run 206's attribution work; extend it.
