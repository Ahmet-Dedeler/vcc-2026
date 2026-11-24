# 214 — Two measured, one-flag positives nobody has composed: drop PerturbMatch, and Hamilton the rounding

You are a research subagent. Work ONLY in `experiments/214-cheap-positives/`.
No `vcc`, no `cell-eval2`, no submission build, no writes outside your directory
(read `src/` and `scripts/`; do not modify them). Peak RSS ≤ 12 GB. `uv run python`.

## READ FIRST
`AGENTS.md` §8.1, §8.2, §11.1 rules 1–13 **with the rule-13 table**. Then:
1. `experiments/200-production-ablation/RESULTS.md` — **the ten-component ablation. Its
   PerturbMatch row is half of this mandate.**
2. `experiments/208-free-riders/RESULTS.md` — **the Hamilton rider and, more importantly, its
   lesson that a rider's value is not a property of the rider.**
3. `experiments/196-emitted-depth/RESULTS.md` — where Hamilton was found and measured.
4. `experiments/210-cd4-build/PREREGISTRATION.md` + `RESULTS.md` — `cd4f034`, submitted at the
   2026-09-04 00:00 UTC reset. **It is the incumbent you must price against**, not `rc_callset`.
5. `experiments/058-perturbmatch/`, `061-perturbmatch-hct-count/`, `064-pm-cd4-composition/`,
   `069-perturbmatch-aligned-whitener/`, `071-production-pm-cd4/` — the PerturbMatch family, so you
   know what it was originally worth and on what instrument.

## The two arms

### Arm A — DROP PerturbMatch
200 ablated all ten of a15's components on all six members. **Nine are positive or neutral. One is
negative:**

| rank | component | six-member sum |
|---:|---|---:|
| 1–6 | union whitening, equal donor weights, `--overdispersion 4.0`, tiered spread, `VCC_ALPHA`, CPM gate 10 | +0.0184 to +0.0440, all **load-bearing** |
| 7 | CD4 residual | +0.002635 |
| 8 | mass fix | +0.000628 |
| 9 | pruned basis | +0.000373 |
| **10** | **PerturbMatch** | **−0.001615** |

200's own verdict was *"none net negative outside its gate; nothing should be dropped"* — i.e. it
declined to act on a −0.0016 that sits inside the ±0.0019 noise floor **as a standalone arm**. That
is correct for a standalone arm and **wrong as a reason not to compose it.** Removing a component
that costs us 0.0016 is the cheapest possible change in the entire stack: one flag, no new code, no
new artifact, and 200 already measured it on all six members.

**But verify before you believe it.** PerturbMatch was promoted on real evidence (058: +0.00517
effect-level PDS across both Orion holdouts; 064: +0.01074 / +0.00640 composed with CD4). If it has
gone from +0.005 to −0.0016, **something changed underneath it** — most likely the operator, the CD4
dose, or the emission. Find out *which*, because that mechanism matters more than the flag.

### Arm B — Hamilton rounding
`emit_counts` at `donor_transfer.py:380` uses `floor`, short-changing every gene ~0.5 counts — a
**−0.01234 log2 size tilt**. Largest-remainder rounding removes 97% of it. 196 measured PLOLO `pds`
0.739016 → **0.742224** (above the `m=∞` depth ceiling) with `claim` still saturated. 208 measured it
marginal to `rc_callset` at **+0.002899**, composition with its other rider **additive to 9e−6**.

**208's lesson is binding here:** its `reach` half was **+0.00607** raw against `rc_hybrid` and
**+0.00006** against `rc_callset`, because `rc_callset` had already driven the delivered order to
ρ +0.9993. **So re-price Arm B against `cd4f034`, not against anything earlier.**

Also note 208's mechanism finding: Hamilton's `pds` gain is **not** direction (Δcos +0.02% relative,
inside 204's ±0.20% floor) — it removes a shared systematic tilt, the way whitening does.

## Work

### Phase 0 — known-answer rows (mandatory, abort on failure)
- Reproduce 200's PerturbMatch ablation row (**−0.001615** six-member) and at least two of its
  load-bearing rows, on its own instrument.
- Reproduce 196's Hamilton row: challenge tilt **−0.01234 → +0.00039**, PLOLO `pds`
  **0.739016 → 0.742224**, `claim` saturated at 0.24653.
- Reproduce 208's marginal-to-`rc_callset` figure **+0.002899** and its `floor_tier` gate
  (`rc_callset`'s h1-test DE table bit-identical, max |diff| 0 on 85 × 10,779).
- Reproduce PLOLO's a15 `pds` **0.795784**, a23 **0.789121**, mass-fix price **+0.004014**.
- **Rule 11:** path and sha256 for every operator. Production is **29,394 × 17,942**, sha
  **`e2f65ecf`**; nine near-identically-named operators exist.

### Phase 1 — why did PerturbMatch turn negative?
Price it at each stage of the stack's history — against a15, `rc_hybrid`, `rc_callset` and
`cd4f034` — and identify what flipped it. Candidate interactions, in the order I'd check them:
the **CD4 dose** (207 tripled it; 064 composed PerturbMatch *with* CD4 at the old dose and found
them synergistic, so a dose change is the obvious suspect), the **operator** (069 fitted a
PerturbMatch-aligned whitener; production does not ship it), and the **emission** (185/191/202).
**If PerturbMatch and the CD4 dose are substitutes, that is a more valuable finding than either
flag**, because it says the two components carry the same information and 207's sweep was run with
PerturbMatch still in.

### Phase 2 — compose, on `cd4f034`
Price A, B and A∘B marginally to `cd4f034`. **Do not assume additivity** — 178 closed the rider
family precisely because *"gains compose sub-additively, costs do not"*, and 208 found additivity
held for its pair only after measuring it.

Endpoints: **Δcos on PLOLO** (204: resolvable to ±0.20% relative ≈ ±0.00006 overall — the sharp
direction instrument); PLOLO `pds` (**`debench.Surface.pds()` must never size a `pds` arm**, 188
Phase 3); `mse` through 096's `claim` route; `nmae`/`reach`/`fid`/`jac` on h1-test with a
**significance-selected pool** (188: `|lfc|·√cpm` at the reference's BH depth, median ~89
genes/row — a top-N-by-magnitude pool inverts `nmae`'s sign). Donors on h1-test:
`k562_gwps` + `orion_hek293t` only; `data/donors/h1_2025` covers 85/85 and is pure leakage.

### Phase 3 — the gate
Rule 13's table: Arm A is **profile-class** (no `reach` toll); Arm B is **emission-law,
profile-CHANGING** — Hamilton alters the per-gene panel totals themselves (208 measured the digest
differing on **900/900** blocks, max |ΔT| = 1, ~6,955 genes/block, all UP), **so it does carry the
−0.0146 toll and claims no exemption.** Price all six for both.

**211's standing rule applies:** `claim` is at its ceiling in the one context that reads it
(context A, `q_A` = 1.000000), so **no arm here may be re-priced on amplitude grounds.** If an arm
moves `A`, re-solve α (211: α\* = 0.9808 on `rc_callset`, **0.9879 on `cd4f034`**) and say what it
cost. Gate **+0.0019** shape-only, **+0.0089** if amplitude moves.

### Phase 4 — falsifiers
1. **Identity rows**: with both flags off, your pipeline must reproduce `cd4f034` exactly — assert
   the per-gene panel-total digest, not just the scalars.
2. **Arm A's sign must not depend on the instrument.** 200 read −0.001615 on its own surface; report
   PLOLO and h1-test separately and say if they disagree. PerturbMatch was promoted on Orion
   holdouts, whose truth carries the long-gene artefact — never quote an absolute from that pair.
3. **Leak-free fold**: report all four PLOLO folds including `kolf2`.
4. **The composition must be measured, not inferred**, and if A∘B is worse than A alone, say so.

## Deliverable
`RESULTS.md`: known-answer rows, the PerturbMatch history and the mechanism that flipped it, the
composed six-member price against `cd4f034`, the gate verdict, and — if it clears — the exact
builder invocation plus a preregistered per-member prediction. **Do not build. Do not submit.**

**"PerturbMatch is a substitute for the CD4 dose" would be the most valuable outcome here**, more
than either flag's points, because it would mean 207's dose sweep was run with a redundant component
still installed and the true optimum is elsewhere.

## Territory
213 owns the transport law. 211 closed amplitude and 212 closed gamma — **do not reopen either.**
You own component removal and emission rounding.
