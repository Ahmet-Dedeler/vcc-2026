# 241 — Re-run 200's ten-component ablation on the renderer we ACTUALLY ship. 219's #1 recommendation, never harvested.

You are a research subagent. Work in `experiments/241-ablation-reprice/`. **No `vcc`, no submission
build, no production change. No `cell-eval2` while another agent runs one — check first.**
`uv run python`. Two other agents are live (221, **240 — it owns `scripts/def_pipeline.py` and
`scripts/def_postbuild.py`; do not open either**).

**MEMORY GATE:** proceed when `memory_pressure ≥ 40% free` AND `vm_stat` free+inactive **≥ 20 GiB**
AND swap free **≥ 400 MiB**. Heavy jobs strictly serially; `/usr/bin/time -l` for peak RSS. A job
over the limit is killed **without a traceback** (219, 226).

## Why this exists

**219 ranked this #1 of five, confidence HIGH, cost LOW — *"the harness exists, no build"* — and its
own mandate then forbade it from running it. Nobody picked it up in the month since.** 239's audit
found it still unharvested **and still being quoted**.

**The defect.** 200's ten-component leave-one-out ablation — the table the **October component
ranking** rests on — was measured with **`VCC_MASS_FIX=minnorm`**, a renderer **only a23 ever
shipped**, and **a23 scored −0.00073 and was abandoned.** 219 §4: *"200's 'six load-bearing, four
inside their gate' ordering is **not safe to quote** until it is re-run with `VCC_MASS_FIX` unset."*

**It is still being quoted.** `docs/DEF-RUNBOOK.md:176` carries 200's PerturbMatch row —
**−0.001615 six-member, "the only NET-NEGATIVE component of the ten"** — as live evidence for an
Oct-22 ship decision, and §3.1 says *"this row and the CD4 row must both be re-decided before
Oct 22."*

**Three of 200's ten rows were priced on only TWO members** (`pds` + `mse`), which is §8.4's exact
error — 219 measured the missed DE propagation at **2.51×** for that class. The CD4 residual row is
one of them, and the board has since measured that arm directly at **+0.0130557**.

## Work

### Phase 0 — known-answer rows (mandatory)
Reproduce **200's own published rows on its own renderer** (`VCC_MASS_FIX=minnorm`) before changing
anything — at minimum the PerturbMatch **−0.001615**, the CD4 residual **+0.002635**, `VCC_ALPHA`
**+0.018542**, and the pruned-basis **+0.000373**. **If they do not reproduce, stop.**

### Phase 1 — re-price all ten on the shipped renderer, all six members
`VCC_MASS_FIX` **unset**, against the **currently shipped** stack (`h1union`, board **0.2118310**),
three cancer folds, **`kolf2` separately**. **Never carry a member at zero** — that is §8.4's whole
lesson and the reason 200's ordering is unsafe. Report the **full ten-row table, before and after**,
and say which rows **change rank**.

### Phase 2 — the two rows the runbook says must be re-decided
1. **PerturbMatch.** Is it still net-negative on the shipped renderer? 200 declined to drop it
   because −0.001615 sits inside the ±0.0019 floor. **Note 231 has since refuted 214's Arm A (drop
   PerturbMatch) on arithmetic** — its +0.0088 required a board raw `reach` of **0.31223** against a
   world best of **0.290835** — so a "drop it" conclusion here must reconcile with that. **If they
   disagree, say so and do not resolve it by preference.**
2. **The CD4 residual row**, which 200 priced on two members at **+0.002635** and the board has
   since measured at **+0.0130557** — **≥ 4.95× understated**, exactly §8.4's failure. Confirm the
   direction of that correction on the shipped renderer.

## Falsifiers
1. **Reproduce 200's own numbers first** (Phase 0) — a re-price that cannot retrodict the original
   is measuring something else.
2. **All six members on every row.** Report `nmae` explicitly with 188's exact identity and
   `e = ΔM/ΔU` — **§8.1d: a component whose contribution is magnitude reads `e` ≈ 0.07 or below.**
3. **`kolf2` separately** — pluripotent, gains 3–5×, and it has already forced one band reissue.
4. **Do not convert PLOLO to board with a scalar** — the factor is **not identified** (1.384 vs
   1.000; 231 found a 7.84× spread over four pairs). Quote PLOLO units and name fold set, functional,
   gate and tag.

## Deliverable
`RESULTS.md`: the Phase-0 retrodiction, the **ten-row before/after table**, which rows change rank,
the two Phase-2 answers, the falsifiers, and a one-line verdict. Then say **exactly which lines of
`docs/DEF-RUNBOOK.md` §3.1 must change** — do not edit it yourself; 240 is editing the runbook.

**Propose no slot. Build nothing.** This re-prices the October component ranking, which is the
decision it serves.
