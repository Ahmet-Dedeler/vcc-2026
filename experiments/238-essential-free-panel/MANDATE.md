# 238 — Re-price the October asset on a panel that EXCLUDES essential genes, as Arc's real panel does.

You are a research subagent. Work in `experiments/238-essential-free-panel/`. **No `vcc`, no
submission build, no production change.** `uv run python`. Two other agents are live (221, 235 —
**235 owns `scripts/make_submission_pm_cd4.py` and `scripts/build_union_whitener.py`; do not open
either**).

**MEMORY GATE:** proceed when `memory_pressure ≥ 40% free` AND `vm_stat` free+inactive **≥ 20 GiB**
AND swap free **≥ 400 MiB**. Heavy jobs **strictly serially**; `/usr/bin/time -l` for peak RSS. A
job over the limit is killed **without a traceback**.

## The problem, and it is mine

233 and 236 measured the four October-only sources (`nadig_hepg2`, `nadig_jurkat`, `rpe1`,
`k562_essential`) at **+0.0066…+0.0104 six-member composite** on a disjoint panel. **Both are
internally sound** — 236's Phase 0 reproduces 233 at max |Δ| **0.0e+00** on 157/157 rows. **The
flaw is in the panels they were given, and I wrote both mandates.**

- **236 §1:** *"all four are common-essential libraries … `rpe1`, whose panel-overlapping rows are
  the essential ones."*
- **`DEF-RUNBOOK` §5 / 199 §3:** `panel ∩ essential list` = **0 (real A/B/C)** vs **34 / 30** on the
  two synthetic draws.
- The synthetic panels were drawn from **donor coverage** (`def-bundle-r` = K562 ∩ Orion ∩ CD4;
  `def-bundle` = any donor). **Donor libraries screen essentials. Arc's panel does not contain
  them.**

**So the rows that make these four sources valuable exist only because the sampling model put
essential genes in the panel, and Arc's one observed panel has none.** There is a strong assay
reason for that — an essential knockdown kills the cell, so there is no transcriptional response to
score — which suggests Arc will exclude them again.

## Work

### Phase 0 — establish the fact, precisely (mandatory)
1. **Reproduce the `0 / 34 / 30` row** from `DEF-RUNBOOK` §5 with your own code, and **state which
   essential-gene list is being used and where it comes from.**
2. **Reproduce 236's headline** for the four sources on the existing pessimistic panel
   (+18.9% relative `cos_agg`, composite +0.0066…+0.0104). If it does not reproduce, stop.
3. **Quantify the dependence directly**: on the existing synthetic panels, what fraction of each
   source's carrying rows is an essential gene? If it is ~100%, the conclusion follows immediately
   and the rest is confirmation.

### Phase 1 — the honest panel
Build a **third** synthetic panel: **300 targets disjoint from A/B/C AND carrying ~0 essential
genes**, matching Arc's one observed panel on both properties. Reuse 199's bundle machinery; change
only the exclusion. **State the draw explicitly and assert both properties.**

### Phase 2 — re-run 233's and 236's measurements on it
- The four October sources, individually and composed, on the essential-free panel.
- **The overall retention number** (233's 65% / 0.137 projection) on the essential-free panel — it
  may also move, and in either direction.
- Report **all six members**, `kolf2` separately, with 233's whole-row-derangement null.

### Falsifiers
1. **Operator re-fit per panel with the axis asserted** — it fails silently otherwise (98.8%
   mis-paired, −0.021 on `pds`).
2. **Carrying-row sizing** — the all-rows median is `0.0000` for all four; the old rule delivers
   identically zero.
3. **`k562_essential` never used on the fold holding `k562_gwps` out.**
4. **Report the essential fraction of each source's carrying rows on all three panels.**

## Deliverable
`RESULTS.md`: the reproduced `0/34/30` row with its list provenance, the essential fraction per
source per panel, the four sources re-priced on the essential-free panel, the retention number, and
a one-line verdict.

**"The October asset is worth ≈0 on an essential-free panel" is the expected answer and a complete
result** — it would retract a strategy claim I made publicly, which is exactly why it must be
measured rather than argued. **If instead the sources survive, say what carries them.**

**Propose no slot. Build nothing.**
