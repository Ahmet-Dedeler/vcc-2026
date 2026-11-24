# 236 — How do the four October-only sources COMPOSE? 233 priced them singly; the stack curve is what we actually ship.

You are a research subagent. Work in `experiments/236-october-stack/`. **No `vcc`, no submission
build, no production change. No `cell-eval2` while another agent runs one — check first.**
`uv run python`. Two other agents are live (221, 235 — **235 is editing
`scripts/make_submission_pm_cd4.py` and `scripts/build_union_whitener.py`; if you read either, take
a byte copy into your own `work/` and pin its sha**).

**MEMORY GATE (corrected — an earlier rule of mine was wrong and stalled work).** `vm.swapusage
free` is **not** a pressure reading: macOS swap is a dynamically sized file whose free tail sits at
0.8–1.7 GB while `memory_pressure` reports 86% free. **Proceed when `memory_pressure ≥ 40% free`
AND `vm_stat` free+inactive **≥ 20 GiB** AND swap free **≥ 400 MiB**.** Heavy jobs **strictly
serially**; report peak RSS with `/usr/bin/time -l`. A job over the limit is killed **without a
traceback** (219, 226).

## What 233 established, and the one thing it did not

**Read `experiments/233-october-retention/RESULTS.md` §5 first, and reuse its instrument verbatim.**
It priced the four sources that carry **0/300 of today's panel** — `nadig_hepg2`, `nadig_jurkat`,
`rpe1`, `k562_essential` — on both `DEF-RUNBOOK` §5 synthetic disjoint panels, **individually**:

| panel | best → worst, relative `cos_agg` @ g = 0.34 |
|---|---|
| pessimistic | `rpe1` **+17.66%**, `nadig_jurkat` **+17.27%**, `nadig_hepg2` **+14.82%**, `k562_essential` **+11.71%** |
| realistic | `rpe1` **+10.49%**, `k562_essential` **+8.85%**, `nadig_jurkat` **+6.44%**, `nadig_hepg2` **+5.51%** |

**Dose-monotone at 0.05/0.10/0.20/0.34 in all eight rows**, shuffled control within ±0.7% and
negative at the top dose in 5 of 8. **All eight clear §8.4's shape-only bar; three clear the
amplitude bar alone.** These are single sources carrying only **8–16%** of the panel.

**233 never composed them.** That is your whole job, and it is the number the October plan turns on.

## Why the answer is not obvious in either direction
225 measured today's sources **saturating at 1.84× the best single at five sources, 1.77× at
thirteen** — *"the fifth source buys 10% of what the first bought"*. **If the October four saturate
the same way, the stack is ~1.8× the best single (~+30% relative cosine on the pessimistic panel) —
which would be transformative, since our board cosine is 0.099 and the field's ceiling looks like
0.13–0.14.** If they saturate harder — plausible, since all four are *low-coverage* sources drawing
on overlapping biology (two Nadig lines, plus `k562_essential`, an essential-gene screen) — the
stack is worth little more than `rpe1` alone.

## Work

### Phase 0 — known-answer rows (mandatory)
Reproduce **233's eight single-source rows** to the last reported digit on both panels, through
233's own instrument, plus its coverage rows (`k562_gwps` **269 → 150** pessimistic, Orion
**291/293**). **If they do not reproduce, stop.**

### Phase 1 — the stack curve, on BOTH panels
Add the four **cumulatively** in 233's own dose form (un-whitened post-whitening residual, sized on
**carrying rows** — the all-rows median is `0.0000` for all four and would deliver identically
zero), on top of the **shipped stack as it will exist in October** (three donors + CD4 + the H1
union, all re-coveraged for that panel). Report **1 → 2 → 3 → 4** and the **saturation ratio vs the
best single**, on **both** panels. Sweep the greedy order; report whether order matters.

### Phase 2 — the six members, not just cosine
Price the best stack on **all six**. **`nmae` decides this, not cosine** (§8.1d/8.1f): 228's arm
survived read-back on cosine and `pds` and died at **−0.00702** on the composite; 229's `ARM6`
looked like **+0.00268** on a rate card and measured **−0.00219**. **Run 188's exact identity and
report `e = ΔM/ΔU`** — a source addition should read as SIGNS (h1union's precedent), a rescale as
MAGNITUDES (`e` ≈ 0.07 or below is the death signature).

### Falsifiers, all mandatory
1. **Whole-row derangement of each source at matched dose** (233's null — it holds magnitude fixed,
   so a real gain must survive it).
2. **`k562_essential` must never be used on the fold that holds `k562_gwps` out** — same cell line,
   same lab (233's rule; keep it).
3. **Both panels.** If pessimistic and realistic disagree in direction, report both and conclude
   nothing.
4. **Leave-top-3-out retention** on cosine and `pds` separately.
5. **Operator re-fit per panel**, with the axis asserted — the operator is a function of the panel
   and fails silently, 98.8% of coordinates mis-paired, priced at **−0.021 on `pds` alone** (233).

## Deliverable
`RESULTS.md`: the reproduced Phase-0 rows, the **stack curve on both panels with the saturation
ratio**, the six-member price of the best stack with `e`, the falsifiers, and a one-line verdict of
the form *"the October four compose to X% relative cosine, saturating at Y× the best single, and the
six-member composite is Z."*

**Propose no slot and build nothing** — none of this is submittable until the panel exists on
2026-10-22. **You are sizing an asset, not building an arm.**
