# 216 — The leaders DOUBLED their cosine in the days after the CD4 atlas was published. Push that lever to its limit.

You are a research subagent. Work ONLY in `experiments/216-cd4-hard/`
(downloads to `data/external/cd4_cell2026/` — that is the ONE exception).
No `vcc`, no `cell-eval2`, no submission build, no production changes.
Peak RSS ≤ 12 GB. Download budget ≤ 60 GB; report the inventory before exceeding 10 GB.
`uv run python`. **Post nothing anywhere and put none of our method into any query** (`AGENTS.md` §14).

## READ FIRST
1. `experiments/215-leader-decomposition/RESULTS.md` — **the finding this mandate executes.
   Verify its Phase-2 windows yourself before building on them.**
2. `experiments/207-cd4-dose/RESULTS.md` — the 25-point sweep, its operating point, and
   **its Phase 2, which found the library form loses on all 57 settings.**
3. `experiments/212-gamma-cosine/RESULTS.md` — the forward measurement, **+13.83% PLOLO cosine at
   `f = 0.34`**, and the exchange rates.
4. `experiments/197-cosine-inversion/RESULTS.md` — the inversion, its one-sided identifiability, and
   **+0.00365 overall per +0.01 of cosine**.
5. `experiments/062-cd4-post-whiten/`, `047`, `050`, `051`, `054`, `060`, `085` — the CD4 family, and
   **how our banks were built and from what source**.
6. `experiments/210-cd4-build/PREREGISTRATION.md` — `cd4f034`, entry `QSe3JlNUYZ2BzvlcsL14`,
   submitted 2026-09-04 00:00 UTC. Its board result may land during your run; check for
   `experiments/210-cd4-build/BOARD-RESULT.md`.

## The finding

215 decomposed 707 board arm-pairs into four channels with a calibrated instrument (5/5 on our own
arms; out-of-sample ρ median **+0.944** on 13 direction-labelled arms and **+0.010** on 7
amplitude-labelled ones). It found:

- **`cqawesome`'s climb is 93% DIRECTION**, and their cosine floor went **0.06427 → 0.13211 (×2.06)**
  between `m0830b` and `m0903c`, with amplitude excluded at **12.4×** the conservative leak bound.
- **`jw_13288` did the same: 0.06333 → 0.13117 (×2.07)**, excluded at 14.9×.
- **Five days ago the two best directions on the board were BELOW ours** (0.0643 and 0.0633 against
  our 0.07221). Over 424 teams, the highest cosine anyone has *ever* carried on their earliest
  observed row is 0.08616. **Nobody started with this.**
- `Aginglab.com` is the control: their apparent direction gain is **93% amplitude** (ρ 0.07).
- **The only public event in the window** is the *Cell* publication of the **Marson/Zhu genome-scale
  CD4⁺ T-cell Perturb-seq atlas** — Crossref record created **2026-08-28**,
  **10.1016/j.cell.2026.08.002**, CC BY-NC. Both doubling windows open **1 and 2 days later**.
  Our own a15 shipped a CD4 residual on 08-30 at **5%** weight for **+6.8%** cosine.

215 is careful that this is coincidence of dates, not causation — the underlying data has been
public since 2025-12-23. **Your job is to find out whether it is causal, and either way to take the
lever to its limit**, because it is the only lever in this project that has ever moved cosine.

## Work

### Phase 0 — known-answer rows (mandatory, abort on failure)
- Reproduce PLOLO's a15 `pds` **0.795784**, a23 **0.789121**, mass-fix price **+0.004014**.
- Reproduce 207's sweep at its shipped points: `f = 0.10` reproducing a15's cosine **0.079249**, and
  212's forward **+13.83%** at `f = 0.34`.
- Reproduce our CD4 bank build byte-for-byte from its recorded source (199 verified the
  re-extraction reproduces the shipped banks on **36 of 36** files).
- **Rule 11:** path and sha256 for every operator; production is **29,394 × 17,942**, sha
  **`e2f65ecf`**, and nine near-identical operators exist.

### Phase 1 — is the published atlas the same object as our banks?
Resolve **10.1016/j.cell.2026.08.002** and its data-availability statement. Then answer, concretely:
does the *Cell* release provide a **processed** version — different QC, different guide assignment,
different normalisation, more cells, more targets, or per-guide rather than per-target tables —
that differs from what our banks were built from? Compare **target coverage against the 300 panel**
(ours is 276) and, on the shared targets, the **effect-vector cosine between their table and ours**.

**If the published version is materially different and better, that is the finding**, and it would
explain a doubling in a way a dose change cannot. **If it is the same object, say so** — then the
date coincidence is a coincidence and Phase 2 is the whole experiment.

Licence: CC BY-NC. §7 permits public data to **train or condition** a model and forbids looking up
an effect to revise output. We already ship this atlas as a learned residual. **Record the licence
position; do not resolve it — that is Ahmet's call**, exactly as 182 and 205 handled it.

### Phase 2 — take the dose to its limit, and re-open the form
207 stopped at `f = 0.80` (its cosine maximiser, **+31.5%**) and shipped 0.34 to protect `pds`.
**A ×2.06 move is +106% cosine — more than three times 207's maximiser.** So:

1. **Extend the sweep well past 0.80** — to 1.0, 1.5, 2.0 and beyond if the curve is still rising —
   with **Δcos on PLOLO primary** (204: resolvable to ±0.20% relative ≈ ±0.00006 overall), `pds` and
   `q_w` secondary. Report where cosine turns over, if it does.
2. **Re-open the library form against 207's Phase 2.** It lost on all 57 settings with `nmae`
   exploding because delivered amplitude collapsed to 0.52× — **but that was at the old operator and
   without re-solving α.** Re-test with α re-solved per arm (211: α\* = 0.9808 / 0.9879) and, as a
   separate arm, **with the operator RE-FIT to include CD4 in its basis**. 207 found CD4 must enter
   **un-whitened**; an operator fitted with CD4 in it is a different question from whitening a CD4
   residual through an operator that never saw it.
3. **CD4 as a fourth equal-weight donor** (the 1/k route). Note the arithmetic that makes this
   interesting: 3 → 4 lines gives only **√(4/3) = 1.15×** under 159's law, so **the 1/k route cannot
   explain ×2.06 on its own** — if the sweep does produce a large gain, the mechanism is not simply
   "a fourth line", and identifying what it is matters more than the points.

### Phase 3 — six members and the gate
Rule 12, all six. `mse` through 096's `claim` route; **211's standing rule: `claim` is at 100.0000%
of its cap in the one context that reads it, so no arm may be re-priced on amplitude grounds.**
`nmae`/`reach`/`fid`/`jac` on h1-test with a **significance-selected pool** (188: `|lfc|·√cpm` at the
reference's BH depth, median ~89 genes/row). **h1-test cannot see this arm** — CD4 shares 4 of its 85
targets and the median sizing rule zeroes the residual silently (200 found it, 207 verified it).
**Assert non-vacuity on every h1-test row or report the member as unmeasured with the mechanism.**
Rule 13: profile-class, no `reach` toll. Gate **+0.0019** shape-only, **+0.0089** if amplitude moves.

### Phase 4 — falsifiers
1. **Target-shuffled CD4** at every dose. 207 measured it strictly negative at all 12 doses with no
   permutation tie — reproduce that, and extend it to your new doses.
2. **The `pds` guard**: 207 chose 0.34 as the largest dose where `|Δpds|` stays under 162's Gate-1
   threshold of 1.645σ. Report σ at every dose; if the cosine optimum is far past the `pds` guard,
   **give both operating points and the six-member price of each**, and let the coordinator choose.
3. **Leak-free fold**: all four PLOLO folds including `kolf2`.
4. **215's own caveat**: `cqawesome` had 5 submissions inside a 2.5-day snapshot gap, so its
   doubling may be one arm or five. Do not build an argument that requires it to be one.

## Deliverable
`RESULTS.md`: known-answer rows, the Phase-1 comparison of the published atlas against our banks,
the extended dose curve with its turnover, the three forms priced, the six-member price at each
candidate operating point, all four falsifiers, and a verdict in overall units. If anything clears,
write the exact builder invocation and a preregistered per-member prediction and **stop** — the
coordinator builds and submits.

**This is the highest-value open question in the project.** The board says the only channel that
matters moved by ×2 last week, on a lever we already hold and have under-dosed since August.

## Territory
213 owns the transport law. 214 owns component removal and rounding — **read its result before you
start**, because if PerturbMatch and the CD4 dose turn out to be substitutes, your sweep must be run
with PerturbMatch off. **You own the CD4 atlas, its dose, and its form.**
