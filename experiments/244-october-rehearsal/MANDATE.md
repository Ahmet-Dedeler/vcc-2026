# 244 — The Oct-22 rehearsal, end to end, on the FIXED driver, on a panel that is disjoint AND essential-free at the pessimistic coverage draw

You are a build-and-verify subagent. Work in `experiments/244-october-rehearsal/`. **You are
authorised to build packages ONLY into `experiments/244-october-rehearsal/work/` — never into
`submissions/`, and NEVER run `vcc submit`.** `uv run python`. **Builds and `vcc prep` are 20–27 GiB
each: run them strictly serially, never alongside another build or prep (`ps aux | grep -E
"make_submission|vcc prep|def_pipeline"` first), and gate every heavy stage on `memory_pressure`
≥ 40% free AND `vm_stat` free+inactive ≥ 30 GiB.** Record peak RSS and wall time per stage.

## Why this exists
240 fixed the Oct-22 driver (`scripts/def_pipeline.py`) and proved its defaults rebuild the shipped
`h1union` package byte-for-byte on A/B/C. 235 added guards (the operator-panel mismatch now raises;
the carrying-row rule covers every registered source). **But the full eleven-stage sequence has not
been run end to end since those changes, and never on a low-coverage panel that is essential-free
the way Arc's is** (`AGENTS.md` §8.1g; 238). `docs/DEF-RUNBOOK.md` says do not improvise on the
day. This is the rehearsal that earns that sentence. 240 §"next" also names one open hole: the
`whitener` stage's skip branch is unguarded.

## READ FIRST
`docs/READ-FIRST.md`; `experiments/247-oct22-guards/RESULTS.md` (the guard fixes and the CD4 carrying-row sizing
flag for D/E/F — your run must exercise them); `docs/DEF-RUNBOOK.md` (all of it; §2 is the sequence, §3 the decision table, §4 the known-answer
rows, §5 coverage responses, §5.0 the measured over-dose factor); `experiments/240-driver-fix/
RESULTS.md`; `experiments/235-october-hardening/RESULTS.md` (its third panel: disjoint AND
essential-free; `work/carry235.json`); `experiments/238-essential-free-panel/RESULTS.md` (its
panel/bundle scripts); `experiments/180-def-readiness/RESULTS.md` §10 (the original rehearsal);
`scripts/def_pipeline.py`, `scripts/def_postbuild.py`; `AGENTS.md` §6, §11.1 rule 10, §14.

## Work
1. **Panel.** Use an essential-free panel disjoint from A/B/C at the PESSIMISTIC coverage draw
   (§5's table: CD4 union ≈ 174/300, `panel_support` ≈ 141/167). Reuse 235's or 238's bundle if it
   meets both properties; say which and cite its census. If none does, draw one with 238's
   scripts and record the census (0 essential genes, disjoint from the 300 A/B/C targets).
2. **Run the driver as the runbook says to run it on the day** — defaults, nothing exported, all
   stages, `--work` in your directory, an opaque `--name`. Capture every stage's stdout.
3. **Verify, with numbers:** the operator was re-fit (sha ≠ `e2f65ecf…`; 235's axis guard passes);
   the CD4 banks re-extracted and the coverage table reproduces 235's `carry235.json` rows for
   this panel; the H1-union frame at low coverage (`carrying_rows`, `g_per_unit_fraction`,
   `solved_fraction`, and `shipped_median_rule_would_solve` = 0); **CD4's over-dose factor on this
   panel** (`residual_median_row_norm_over_carrying_rows` ÷ `…over_all_rows` — compute it from the
   staged residual; §5.0 predicts ~1.24× at this draw); the built artifact's `effects.cd4.fraction`
   and `target_f_requested` 0.34, `h1_union_g` 0.28; `def_postbuild.py` **20/20**; `vcc prep
   --dry-run --contexts D,E,F` clean with nnz vs the 4.75e9 cap and max counts/cell vs 1e6.
4. **The skip-branch hole (240 §next).** Experiment 247 (`experiments/247-oct22-guards/RESULTS.md`) may already have
   closed this and the other 242 items on the driver and post-build check — read its RESULTS first and VERIFY its
   guards fire in your end-to-end run rather than re-implementing them. If 247 did not close it, demonstrate it: run the `whitener` stage against a stale
   operator file and show whether the driver silently reuses it. If it does, add the smallest
   possible guard to `scripts/def_pipeline.py` (no deletions; a raise keyed on the operator's
   recorded panel/axis sha vs the bundle's), show it FIRING on the stale operator and NOT firing on
   the fresh one, and then **re-run 240's row 1** (defaults on A/B/C rebuild `h1union`
   byte-for-byte: `.h5ad` `8df4dadd…`, `.vcc` `f13ce40b…`) so the change is proven inert on the
   shipped path. If the hole is already closed, say so with the line that closes it.
5. **Timing.** A stage-by-stage table (wall, peak RSS) against the runbook's ≈ 29 min figure, and
   the disk used.

## Known-answer rows (abort on failure)
- 235's guard fires on a wrong-panel operator (reproduce its 4-way demonstration once).
- The coverage table reproduces 235's rows for the same bundle.
- If you touched `scripts/def_pipeline.py`: 240's row 1 byte-for-byte on A/B/C.

## Deliverable
`RESULTS.md`: the panel and its census, the stage table, every verification number, the over-dose
factor measured on this panel, the skip-branch verdict with the diff (if any) and its firing
demonstration, and one line: **RUNBOOK EXECUTABLE AS WRITTEN / NOT**, with the exact runbook lines
that must change. Do not edit `docs/DEF-RUNBOOK.md`; put the proposed edits in your RESULTS. Stop.
