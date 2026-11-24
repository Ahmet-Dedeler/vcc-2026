# 247 — Close the Oct-22-path guards 242 found, each proven FIRING, and make CD4's dose panel-invariant for D/E/F

You are a fix-and-prove subagent. Work in `experiments/247-oct22-guards/`. **The ONLY files you may
edit outside it:** `scripts/def_postbuild.py`, `scripts/def_pipeline.py`, `scripts/def_census.py`,
`scripts/check_submission_novelty.py`, `scripts/make_submission_pm_cd4.py` (one flag, default off,
see §C), `.gitignore` (one line, §D), `submissions/SUBMITTED_LEDGER.json` (backfill only, §D).
**Do NOT edit `scripts/arm_submit.sh`** — a process is executing it (bash reads scripts
incrementally); report its one-line patch in RESULTS instead. Never `vcc submit`. You may build a
package ONLY into `experiments/247-oct22-guards/work/`, and only for the byte-identity row in §C.
`uv run python`. **Memory:** 241, 243, 245, 246 are running. Anything ≥ 15 GiB (a staging run or a
build) waits for `memory_pressure` ≥ 40% free AND `vm_stat` free+inactive ≥ 30 GiB, and runs alone
(`ps aux | grep -E "make_submission|vcc prep|def_pipeline|stage2"` first). Record peak RSS.
Commit each file you change with its known-answer row hashes in the message (`git add <file>`,
never `git add -A`).

## READ FIRST
`experiments/242-guard-audit/RESULTS.md` (all of §1 and §2 — this mandate is its fix list),
`experiments/240-driver-fix/RESULTS.md` (row 1: the byte-identity gate, `.h5ad` `8df4dadd…` on
12,648,807,800 B, `.vcc` `f13ce40b…`), `experiments/235-october-hardening/RESULTS.md` (the axis
guard; `axis_sha256` `edf6ba88…` for A/B/C), `experiments/222-h1-union-build/work/stage-on028.json`
(the CD4 frame: `residual_median_row_norm_at_unit_fraction` 4.464717, `whitened_base_median_row_norm`
4.229296, `solved_fraction` 0.32207201067578145), `docs/DEF-RUNBOOK.md` §5 and §5.0,
`experiments/223-silent-zero-audit/work/overdose223.json` (the 1.4222× at 58.7%), `AGENTS.md` §11.1
rule 10 (**every guard must be demonstrated firing, not asserted**), §14.

## A — `def_postbuild.py` operator block (242 #1, #7)
Assert the field that carries the answer: `meta["operator"]["axis_sha256"]` == the axis sha
recomputed from THIS bundle (235's construction), and `meta["operator"]["sha256"]` == the sha256 of
the operator file the driver pinned for this panel (read the pin, hash the file). Fix #7 so the row
count is actually compared (`rows == n_constructs × 400 × n_contexts`). **Demonstrate:** on
`submissions/222-h1union-g028.h5ad.meta.json` the new checks PASS (20/20 → 22/22 or whatever the
count becomes); on a COPY of that meta in `work/` with the operator block replaced by
`data/donors/whitener_union.npz`'s (sha `5608e807…`, same `n_cols`) the operator check FAILS; on a
copy with a wrong row count the row check FAILS. Quote the verbatim FAIL lines.

## B — `def_pipeline.py` and `def_census.py` (242 #8, #10, #18)
- #8: `other_builders()` must catch all six launch forms 242 replayed (interpreter flags, `runpy`,
  the `uv` shim). Demonstrate on each form with a sleeping dummy process.
- #10: the census stage must be fail-closed on D/E/F: add minimum-coverage thresholds that the
  PESSIMISTIC draw still passes (§5's table: CD4 union 174, `panel_support` pm 141 / cd4 167,
  0-carrier targets ≤ 5) — e.g. fire if any production donor covers < 100 targets, CD4 union < 120,
  or 0-carrier targets > 5 — and demonstrate it firing on a synthetic census with collapsed coverage
  and passing on 235's pessimistic-draw census.
- #18: the `whitener` stage's skip branch: when an existing operator is passed, assert its recorded
  `axis_sha256`/panel matches the bundle, else raise. Demonstrate with a stale operator (the A/B/C
  file against a synthetic D/E/F bundle from 235/238) and with the fresh one.

## C — CD4's dose must not silently inflate with coverage (242 #2; 222/223's 1.24–1.42×)
`solve_cd4_fraction` sizes by the median over ALL rows; below ~92% coverage that over-doses by
`median_carrying / median_all` (1.012–1.017× on A/B/C; **1.4222× measured at 58.7%**). Add ONE flag,
default off, e.g. `VCC_CD4_SIZING=carrying` with `VCC_CD4_TARGET_G=<float>`, that sizes CD4 on its
carrying rows exactly as `solve_carrying_row_fraction` does for the union. **Flag off must be
byte-inert:** the staging-stage `effects_array_sha256` `87a4b6d1…` (222's `stage-off.json`) must
reproduce, and if memory allows within your run, the full A/B/C package sha (240 row 1). The D/E/F
target is the A/B/C DELIVERED carrying-row dose: `g_cd4 = 0.34 × (median_carrying / median_all)`
measured on A/B/C — compute it, report it to full precision, and show that with the flag ON at that
`g_cd4` the A/B/C `effects.cd4.fraction` reproduces **0.32207201067578145** to ≤ 1e-9 (numeric, not
byte, identity — say so). Then on the pessimistic panel (235/238's bundle) report the delivered
per-row dose under both rules relative to A/B/C's. Wire the driver: `SHIPPED_ENV` keeps the legacy
path on A/B/C (so 240's row 1 stays byte-for-byte) and sets the carrying rule with that `g_cd4` on
any non-A/B/C bundle; `def_postbuild`'s shipped-arms check accepts the new fields. If you cannot
make the flag inert, stop and report — do not ship a builder change that moves a shipped byte.

## D — the novelty ledger (242 #3)
`submissions/SUBMITTED_LEDGER.json` is gitignored (`submissions/*`) and five entries stale, so
`check_submission_novelty.py check` passes already-scored packages (demonstrated on `abc-knee.vcc`).
Un-ignore that one file (`!submissions/SUBMITTED_LEDGER.json`), backfill every scored package since
2026-09-02 from `docs/BOARD-LEDGER.md` and the `BOARD-RESULT*.md` files (entry id, package path,
sha256 where recorded — hash the file if it is on disk), and demonstrate: `check submissions/abc-knee.vcc`
→ ALREADY SCORED, non-zero exit; `check submissions/229-stack6-g035.vcc` → OK novel. Report the
one-line patch that makes `scripts/arm_submit.sh` call the check before submitting (do not apply it).

## Deliverable
`RESULTS.md`: for every item A–D the diff summary, the known-answer rows with the verbatim PASS and
FAIL lines, the measured `g_cd4` and the two delivered-dose numbers, the commit hashes, and one line:
**OCT-22 PATH GUARDED / NOT**, listing anything from 242's table you did not close and why. Then stop.
