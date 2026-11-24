# 254 — CD4's partial-condition rows as their own residual (246 §3.7), inert by default; then BUILD the composed arm `cd4svd5+partial` for a later slot. DO NOT SUBMIT.

You are a build subagent. Work in `experiments/254-partial-rows/`. **The only file you may edit outside it is
`scripts/make_submission_pm_cd4.py`** (changed today by 222, 247, 249, 248, 251 — read their RESULTS; the current
builder is git HEAD; do not undo anything). You may build into `submissions/` under the `254-` prefix only.
**NEVER `vcc submit`.** `uv run python`. Heavy stages strictly serial and alone (`ps aux | grep -E
"make_submission|vcc prep|def_pipeline"`), gated on `memory_pressure` ≥ 40% free AND `vm_stat` free+inactive
≥ 30 GiB; record peak RSS. 252 (harness) and 253 (docs) are running; neither touches the builder. Commit the
builder with its known-answer hashes.

## Why
246 §3.7: production's `complete mean3 only` rule zeroes CD4's partial-condition rows (targets carried by 1–2 of
the three CD4 conditions; 15–24 rows per fold). Carrying-row sized at 0.34 as their own residual (`C2_partial@0.34`)
they read macro +0.587% / `pds` +0.00178 (3/3, 3/3) — NOT build-worthy alone (+0.0009…+0.0018, concentrated) —
but **composed with the k = 5 deflation (`E_cd4_svd5+partial@0.34`) the pair reads `pds` +0.01405 (3/3), macro
−0.699%, leave-top-3-out 61.5%, price +0.0032 / +0.0047 / +0.0100** (246 §4), additive to within 0.0002 of the
parts. If tomorrow's `cd4svd5` board point lands, this composition is the next one-variable step on top of it;
the capability and the package must exist before that decision, not after.

## READ FIRST
`docs/READ-FIRST.md`; 246 RESULTS §3.7, §4 and `lib246.py` (family C: how the partial rows are extracted —
`target_condition_count` in {1, 2} — and `merged_many`); 249 RESULTS §1 (the deflation flag; the composed arm
must run deflation on the COMPLETE rows exactly as `cd4svd5` does, then add the partial-row residual sized on
its carrying rows), 249's `build249.sh`/`stage249.py`/falsifier scripts (reuse); 247 §3 (sizing rules);
`experiments/222-h1-union-build/` gate scripts; `AGENTS.md` §11.1 rule 10, §14.

## The change
`VCC_CD4_PARTIAL_G=<g>` (unset = inert). When set, the CD4 rows with `target_condition_count` in {1, 2} are
assembled as their own residual (mean over the conditions that carry them, `gene_norm`, the same rescale as the
complete rows), sized by the carrying-row rule at `g` against the same base denominator, and merged into the
one add + one preserve with the complete-row CD4 residual (deflated or not) and the union(s). Meta:
`effects.cd4_partial` (rows, carrying rows, solved fraction, `shipped_median_rule_would_solve`) and
`config.cd4_partial_g`, written only when on. The complete-row path, the deflation and the unions must be
byte-identical with the flag off.

## Known-answer rows (abort on failure)
1. Flag unset: staging `abeb0749…` reproduces, and the A/B/C package `8df4dadd…` on 12,648,807,800 B (251
   rebuilt it two hours ago through HEAD's builder; rebuild once if its verified copy is gone).
2. `VCC_CD4_PARTIAL_G=0.34` alone (no deflation) reproduces 246's `C2_partial@0.34` on the three cancer folds
   through 162's instrument on the builder's staged effects: macro +0.587%, `pds` +0.00178 — to the digit.
3. `VCC_CD4_DEFLATE_K=5` + `VCC_CD4_PARTIAL_G=0.34` reproduces `E_cd4_svd5+partial@0.34`: `pds` +0.01405, macro
   −0.699%, and its leave-top-3-out 61.5% (246's per-target artifacts).

## Build (conditional on the three rows)
`submissions/254-cd4svd5-partial-f034.{h5ad,vcc}` with `VCC_CD4_DEFLATE_K=5 VCC_CD4_PARTIAL_G=0.34` and everything
else the shipped configuration; 222/249's full gate set against `submissions/249-cd4svd5-f034.h5ad.meta.json`
as the incumbent (rule-10 diff = exactly `cd4_partial_g` plus its consequences; scale-corrected profile split
by partial-row; non-vacuity; `n_pred` ±2%; delivered amplitude within 0.5% via 243's stream; `vcc prep
--dry-run`). Record sha256 and bytes. Write a preregistration draft (`PREREGISTRATION-DRAFT.md`, six-member
prediction vs `cd4svd5` from 246 §4, five falsifiers) for the lead to commit before any slot.

## Deliverable
`RESULTS.md`: the three rows, every gate's number, the package, the builder commit, the draft preregistration,
and one line **BUILDABLE / NOT BUILDABLE**. No model name, no description, no submission. Then stop.
