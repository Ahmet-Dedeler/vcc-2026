# 531 — Final-phase reader, per-row commands, and the Oct 22 → Nov 5 tabletop

2026-09-29. No `vcc`, no contact, nothing written outside this directory. This VM has no gitignored challenge matrices, so the tabletop was walked against the runbook and the code, and the reader was tested on synthetic captures plus one real entry-id lookup. `uv run python experiments/531-cloud-final-phase-tooling/test_readdef.py` exited 0: 74 row-verdicts, including both-lose, the q-curve turn, the DE veto at exactly 0.009, prior-only reads under the 0.009 floor, read468's boundaries, and members hidden.

The operational conclusion: R0 is still the first upload, and the last upload is still an already-scored package started by Nov 4 20:00 UTC. Three feedback rows cannot be built until knobs land (ρ\* scale and the CD4 ×1.5 stack, both named as experiment 530; the called-gene scale s, which 530 was not named as owning). A 7.5 h upload makes two slots a day the exception, so the plan's full arm list does not fit unless uploads move.

## 1. Reader

```bash
uv run python experiments/531-cloud-final-phase-tooling/readdef.py ROW ENTRY BASE [ALT]
```

`ROW` is a section-2 id: `1`, `2`, `3`, `4`, `5`, `6`, `7`, `8-9`, `10-12`, `13-14`. The aliases `8`, `9`, `10`, `11`, `12`, `13`, `14` map onto the grouped rows. `ENTRY`, `BASE` and optional `ALT` are capture files (the `capture_board_entry.py` wrapper or a bare entry) or ids under `data/board-entries/`. Exit 0 when a verdict was printed, including a decision not to adopt. Exit 2 when a reading row was given hidden or absent members. Exit 1 on a missing file or an unknown row.

It prints every raw delta and every scaled delta, the scaled delta's contribution to overall (`scaled / 6`, the same convention as `experiments/456-zs-riders/read0928.py`), then `BIN` and `ADOPTION`.

How the plan's words were encoded:

| Row | What is read | Adopt |
|---|---|---|
| 1 | nothing decides | Baseline. No adoption. |
| 2 | bulk = `(Δscore_pds + Δscore_mse) / 6`, exactly | The better side if it is strictly positive. Both ≤ 0: R0's ρ\* stands. A tie adopts neither. One side alone is not a both-lose read. |
| 3 | `Δscore_avg`, exactly. First arm is the lower q, second is the higher q | Each q must beat the day-2 best. If the higher q is not strictly better than the lower q, the curve has turned: adopt the lower q and stop. If the higher q is better, adopt it. Neither beats the base: the day-2 best stands. |
| 4 | `Δscore_avg`, exactly, each arm on its own | Adopt every arm that is strictly positive. The LODO cap and CD4 ×1.5 are different switches. |
| 5 | `Δscore_avg` exactly, plus bulk and DE | Adopt if overall is strictly positive, unless bulk is up and the DE contribution `(Δfid+Δreach+Δjac+Δnmae)/6` is `< -0.009`. A cost of exactly 0.009 does not veto. Two w that both clear: the better overall. |
| 6 | `(Δscore_fid + Δscore_nmae + Δscore_reach) / 6`, threshold `≥ +0.0015` | Prior-only. The 0.009 floor does not apply, and a large overall that comes from other members does not adopt. |
| 7 | read468 on this pair | Falsifier if raw pds or mse exceeds 1e-6, or raw jac exceeds 0.0010: do not adopt. Else `≥ +0.0015` transfers, `≥ +0.0005` weak-adopt, `≤ -0.0005` HARM closes the gate, and the open interval keeps le1. Prior-only. |
| 8-9 | raw nmae, raw pds, raw mse, `Δscore_avg` | Adopt only if nmae `≤ -0.002`, pds and mse are pinned (`|raw| ≤ 1e-6`), and overall `> 0.009`. Exactly 0.009 does not clear. Under 0.009 asks for a second seed. Not prior-only. Two scales that both clear: the more negative nmae, not both. |
| 10-12 | `Δscore_avg`, exactly, and the members must be present | Scored and strictly better: adopt the combination. Scored and not better: keep the base. The ρ\* re-check is a separate row-2 call against the combination. |
| 13-14 | nothing | Re-upload the best scored package. No new arm. |

A capture whose twelve score fields or `score_avg` are missing prints `MEMBERS HIDDEN` and adopts nothing. Rows 1 and 13-14 do not read members, so a dark record still returns baseline or reserve.

## 2. Commands for rows 2–12

R0 is whatever `uv run python scripts/october_recipe.py command` prints with exit 0 after the last A/B/C read. On 2026-09-29 that command exits 2. With `--allow-pending` it prints the line below and treats `lodo_hi_pct`, `sharpen_pt`, `programme_source`, `prior_gate`, `mde_sign` as off. Do not pass `--allow-pending` on Oct 22.

```text
DEF_DAY_LODO=1 DEF_DAY_DE_MAG=dn_k250_zs bash scripts/def_day.sh \
  data/challenge <work> <stem> 1.75 0.5 noprog 4
```

Every arm below is that command, in a fresh work directory, with the one change named in the row. `<work>` is a new directory every time (`experiments/531` is not a build dir). The positional args stay `1.75 0.5 noprog 4` unless the row changes one of them. Env that the 10-05 reads turn on (percentile cap, per-target sharpening, programme, prior gate, MDE sign) stays at R0's value except where the row says otherwise.

`def_day.sh` itself diffs the artifact against 446b's meta (`scripts/def_day.sh:314-328`) and must log `postbuild PASS`. That allow-list is the October-versus-A/B/C list. It is not the one-variable gate. The one-variable gate is a second postbuild against the D/E/F base this row is compared with:

```bash
uv run python scripts/def_postbuild.py \
  --meta <arm-work>/<stem>_final.h5ad.meta.json \
  --bundle data/challenge \
  --cd4-root <arm-work>/cd4_compact \
  --incumbent <base-work>/<base-stem>_final.h5ad.meta.json \
  --allow-config-diff <keys>
```

`def_postbuild.py` passes when the moved keys are a subset of the allowed set (`scripts/def_postbuild.py:465-469`). The feedback gate is stronger: the moved set must contain every must-move key and nothing outside must-move ∪ companion. A knob that did not fire still exits 0. Read the `moved=` list and fail the arm if it is wrong. Do not submit that arm.

| Row | Command, once the named knob exists | Compared with | Must move | May also move | Gate |
|---|---|---|---|---|---|
| 2 | `DEF_DAY_RHO_SCALE=0.98` and `=1.02` on the R0 command. Two fresh dirs. | R0 | `bulk_ctx_scale`, `rho_scale` | none | Row 2 of the reader. Also: every shipped context's ntclevel ρ is within 1% of `0.9793434758 × scale` and `≤ 1`, and the log says `CONVERGED`, not `nearest safe build`. |
| 3 | Positional q `1.5` and `2.0`. `DEF_DAY_SHARPEN_PT` must be unset, which is only a one-variable arm if R0 left it unset. | Day-2 best | `bulk_sharpen_q`, `bulk_ctx_scale` | none | Row 3. Stop submitting the higher q once the reader says the curve turned. |
| 4a LODO | If R0 shipped `DEF_DAY_LODO_HI_PCT=96.64697522769474`, the arm is fixed cap 6: that variable unset, `DEF_DAY_LODO_HI=6`. If R0 shipped neither (cap 3, the 495-fails branch), skip this half. | R0 | `bulk_load_sha256` | `bulk_load`, `bulk_ctx_scale` | The dump key's json `hi` is 6, not the percentile and not 3. Reader row 4 on this pair. |
| 4b CD4 | `DEF_DAY_CD4_STACK=x1.5` on the R0 command. Dump subshell only: `VCC_CD4_TARGET_G=0.5165482917974462` (1.5 × `0.3443655278649641`). The level builds keep the carrying g. | R0 | `bulk_load_sha256` | `bulk_load`, `bulk_ctx_scale` | Level-meta `cd4_target_g` is still `0.3443655278649641`. Reader row 4. |
| 5 | Positional w `0` and `1`, only after `def_day_helpers.py nctx` on R0's sameline artifact is `> 0`. | R0 | `sameline` | `bulk_ctx_scale` | Skip the row when nctx is 0 (the bytes would match R0). Reader row 5. |
| 6 | `DEF_DAY_PROGRAMME_SOURCE=h1_2025_full` and positional flags `""` instead of `noprog`. Skip if the 09-30 read closed the source, or if R0 already shipped it (then the arm is the off command: the variable unset and flags `noprog`). | The off package | `sameline` | `bulk_ctx_scale` | Reader row 6. Prior-only, so one entry is enough. |
| 7 | R0's programme command plus `DEF_DAY_PRIOR_GATE=le1`. Run only if row 6 adopted the programme (or R0 already shipped it). | Programme on, gate off | `de_prior_gate` | none | Reader row 7. A moved `bulk_ctx_scale` is unexpected (469 retraced the factors). |
| 8-9 | `DEF_DAY_DE_CALLED_SCALE=1.25` and `=1.6`. This is not `DEF_DAY_DE_MAG` (that is the zs rider and stays at R0's value). | Current best | `de_called_scale` | none | Reader row 8-9. pds and mse raw pinned, so a moved bulk key fails the gate even if the reader would also reject it. |
| 10-11 | The R0 command with every switch the reader adopted, all at once. One fresh dir. | The current best single package | the union of the adopted must-move keys | the union of their companions | Reader row 10-12. Rule 5 is satisfied only because each switch was already scored alone. |
| 12 | Row 2's two `DEF_DAY_RHO_SCALE` values on the combination command. | The combination (scale 1) | `bulk_ctx_scale`, `rho_scale` | none | Reader row 2, not row 10. Same convergence check as row 2. |

### Knobs that do not exist yet

**ρ\* scale (row 2 and row 12).** `rho_star` in `scripts/def_day_helpers.py:28-32` and `experiments/448-def-rehearsal-q175/relevel_iter.py:28-29` takes `sd` and the fixed `q = 0.31 / (0.31 + 1.78)` and has no multiplier. Unset, the target is `0.9793434758`. The ×0.98 target is `0.9597566062`. The ×1.02 target is `0.9989303453`. Experiment 530 is the named owner. Once it lands, `DEF_DAY_RHO_SCALE` multiplies that target in both places, defaults to unset = 1 (bit-identical re-level), and the level meta records `rho_scale`. The commands above are the commands after that lands. Until then the row cannot be built, and a hand-edited factor string is not a substitute: the secant loop would walk it back to the unscaled target.

**CD4 stack (row 4b).** `scripts/def_day.sh:194-197` hardcodes the dump as `VCC_H1_UNION_G=0.56` and the carrying CD4 default. `scripts/october_recipe.py:76-77` refuses `cd4_stack` other than `x1`. Experiment 530 is the named owner. Once it lands, `DEF_DAY_CD4_STACK=x1.5` is legal, acts only inside the dump subshell, and does not change the level build's `cd4_target_g` (the postbuild arm check would fail if it did). Do not `october_recipe.py decide cd4_stack x1.5`: that would change R0. The feedback arm passes the env on that one command. The command above is the command after 530 lands.

**Called-gene scale (rows 8-9).** No builder variable multiplies the delivered log-fold change of called genes. `experiments/513-reach-walk-anatomy/scale513.py` is analytic. `DEF_DAY_DE_MAG=dn_k250_zs` is a different rider and is already in the recipe. This knob was not named as 530's work. The command above is the command once a pass-through exists; the row is skipped if it does not exist by the day it would be built.

**Fixed LODO cap 6 (row 4a).** `scripts/def_lodo_key.py:61-63` has `--hi` (default 3) and `--hi-pct` (overrides `--hi`). `scripts/def_day.sh:179-183` can pass `--hi-pct` only. Unset, the cap is 3, not 6. `DEF_DAY_LODO_HI` does not exist. It is not part of the named 530 scope. The command above is the command once a pass-through exists that passes `--hi` and refuses to be combined with `--hi-pct`.

## 3. Tabletop, Oct 22 00:00 UTC → Nov 5 23:59 UTC

Clock times below assume the bundle is downloadable at **2026-10-22 12:00 UTC**. If it is later, every Oct 22 time shifts by the same amount and the slip rule in the stress section decides whether Oct 23's slot is an upload completion rather than a new arm. Durations used: `def_day.sh` about 1 h 45 min (452's rehearsal was 1 h 47 min, 465 was 1 h 43 min), real `vcc prep` about 20 min, upload **planned at 7.5 h** (the 2026-09-29 observation; the 15 min figure is the old line). Scoring time is unknown until the first D/E/F status returns.

### Oct 21, the entry gate (runbook section 1, done before the 22nd)

Confirm the pinned shas in runbook section 1, `def_preflight.py --verify`, `def_identity.py --verify`, `def_census.py --verify`, `vcc --version` with `prep --help` still showing `--contexts`, and ≥ 60 GB free. Do the `vcc` upgrade here if it is going to happen at all, then re-run the section 4 prep row on a rehearsal artifact before anyone goes to sleep. Ask for the two downloads (the bundle, and the CLI only if it was not upgraded today).

### Oct 22, R0 only

| UTC | What | Gate, fail closed |
|---|---|---|
| 00:00–11:00 | Slots reset at 00:00. Nothing to build. Check `vcc datasets list` on the hour from 08:00. Do not run the day-page `uv tool upgrade vcc-cli`. | No D/E/F id yet: keep waiting. An id of `controls` only is still A/B/C (`docs/DEF-RUNBOOK.md:34-35`). |
| 12:00 | `cp -R data/challenge data/challenge-abc`, then `vcc datasets download controls` into `data/challenge`. | Manifest exists. Contexts are three identifier names. Copy of A/B/C is intact. |
| 12:15 | `october_recipe.py show`, then `october_recipe.py command` **without** `--allow-pending`. | Exit 0. Any `PENDING` or `# REFUSED` stops the day: an undecided switch is not "off". |
| 12:20 | `DEF_DAY_KEEP_PINS=1` plus the resolved command, fresh work dir `experiments/NNN-def-r0/work`. Claim the next model name with `next_model_name.py --claim r0` before the build, and do not print it anywhere public except the eventual `-m`. | Setup log reaches `repin-panel-support rc=0`. Identity JSON written into the round's RESULTS. |
| 12:25 | Preflight, inside that run. | 66/66. `ntc_id` column present and its category count equal to the manifest (absent is finding 5; do not continue on a silent pass). Gene header `gene_name`. No duplicate targets. `cells_per_pert` is 400. Context labels match the filenames. |
| 12:40 | CD4 extract and census. | If the log says the CD4 extraction gate fired: stop, run `experiments/244-october-rehearsal/mapping244.py`, then rerun the whole day in a **new** work dir with `DEF_DAY_CD4_THRESHOLDS=relaxed`. Census: every production donor ≥ 100, CD4 union ≥ 120, CD4 complete ≥ 100, zero-carrier targets ≤ 5. More than 5: re-check symbols before believing it; there is no override flag (finding 16). |
| 13:10 | Whitener. | New sha, not `e2f65ecf…`. Axis sha not `edf6ba88…`. A stale npz in the work dir must refuse (`verify_existing_operator`), not skip. |
| 13:15 | Same-line. | Log the three calls and whether each was unanimous. nctx 0 with no programme source is a valid R0. nctx 0 with the programme source set is a dead day (finding 9): turn the source off and rerun in a fresh dir rather than miss R0. |
| 13:20–15:10 | Dump, lvl0, re-level loop, audits, prep dry-run. | `CONVERGED` or an explicit nearest-safe ship with every ρ ≤ 1. Zero-gene audit says 0. Prep dry-run: 360,000 cells, the manifest's three context names, not a hardcoded `D,E,F`. |
| 15:10 | `postbuild PASS` against 446b. Record shas of the h5ad and, after the real prep, the vcc. | Moved keys are only the October allow-list. `cd4_target_g` is the carrying value. |
| 15:20 | Real `vcc prep` with `--contexts` taken from the manifest. Commit the pins (`DEF_DAY_KEEP_PINS=1` left the builder pinned). Then restore the builder from `pin-artifacts.json`. | `git diff` of the builder is the pin you just reviewed. After restore, the builder sha matches the pre-pin sha. |
| 15:30 | `vcc submit` of R0. `-m` the claimed queue name. `-d` empty. One in flight. | Novelty check passes. Ledger row written. |
| 16:00, 18:00, 20:00, 22:00 | `Content-Range: bytes */<size>` against the pending upload URL in `~/.config/vcc/state.json`. | The range is moving. A range unchanged for 30 min is a stall: restart that upload. A slow range that is still advancing is the 09-29 pattern (148 MB after 3.6 h, then recovery): leave it. |
| 15:30–23:00 | While the upload runs, prebuild nothing that pins the builder. The next arm's work dir can be created and its command written down. Do not start a second submit. | Oct 22 slot 2 expires unused if the upload is still in flight at 16:00. That is the intended outcome under a 7.5 h upload. |

Decision at the first terminal status, before any second package is submitted:

- `vcc status <id>` and `vcc --json status <id>`, then `scripts/capture_board_entry.py`.
- Members present (all twelve score fields and `score_avg` numeric): the feedback plan is on. Row 1 of the reader records R0. No adoption.
- Status returns nothing, or a terminal record whose members are null: the feedback plan is off. DEF-RUNBOOK section 6 stands. Do not build rows 2–12. The Nov 4 upload is this same package, under the reserved final name.
- Capture exits 1 on a published record with null scores: that is not "still scoring". Today the script cannot write that record (finding 8). Treat a repeated reject of a terminal id as members-hidden and switch to the dark plan.

### Oct 23 → Nov 3, one new arm a day unless the upload was fast

The feedback order is the plan's order. Skip a row whose switch cannot apply (no same-line match, programme closed on 09-30, knob not landed, 495 closed the percentile so the fixed-6 half has nothing to compare with).

| Day | Arms, in order, each only if the previous gate passed | Reader |
|---|---|---|
| Oct 23 | ρ\* ×0.98, and ×1.02 only if the first upload has published | row 2, both sides together |
| Oct 24 | q 1.5, then q 2.0 only if 1.5 has not yet shown a turn and the slot is free | row 3 |
| Oct 25 | the LODO alternative R0 did not ship, then CD4 ×1.5 if the slot is free | row 4, each arm alone |
| Oct 26 | w 0, then w 1 | row 5. Skip the day if nctx was 0. |
| Oct 27 | programme on (or off, if R0 shipped it) | row 6 |
| Oct 28 | prior gate le1, only if the programme is on | row 7 |
| Oct 29–30 | called-gene s 1.25, then 1.6 | row 8-9 |
| Oct 31–Nov 2 | the combination of every adopted switch, then its two ρ\* re-checks | row 10-12, then row 2 against the combination |
| Nov 3 | buffer. A slipped arm, or nothing. | — |

Daily shape, slot reset 00:00 UTC:

| UTC | What | Gate |
|---|---|---|
| 08:00 | Yesterday's status is terminal and captured. Run the reader. Write the bin in the round's RESULTS before any new submit. | A non-terminal status: do not submit a package whose decision needs that read. Prebuilding the next command is fine. |
| 08:30–11:00 | Build today's arm in a fresh work dir. `KEEP_PINS=1`. | Postbuild vs 446b passes, and the sibling diff's moved set equals the row's keys. Audits clean. Prep dry-run clean. |
| 11:00 | Submit slot 1 if the upload window is the constraint. `-m` next queue name, claimed when the command was written. `-d` empty. | Not a re-submit of a scored package. Novelty check. |
| 16:00 | Slot-2 decision. | Start slot 2 only if slot 1 has **published** and the second package is already built and gated. Otherwise slot 2 expires. A 7.5 h upload started at 11:00 finishes at 18:30 and does not earn slot 2. |
| overnight | Leave the upload. Check the byte range at 20:00 and 23:00. | Stall rule as on Oct 22. |

Build the next day's arm while today's upload is in flight only when the builder has been restored and no pin is held. Never run `cell-eval2` or `vcc prep` beside a build.

### Nov 4 and Nov 5, the last entry

| UTC | What | Gate |
|---|---|---|
| Nov 4 12:00 | Freeze the candidate. It is the best package the reader has adopted, and if that package is a combination it must itself have a capture with members. | No members yet on the combination: the candidate is the best package that does have members. Do not assemble a new combination today. |
| Nov 4 16:00 | The `.vcc` exists, shas are written, pins are committed, postbuild of that exact file passed. `-m` is the reserved name `Be extremely selfish. Be extremely generous.` Claim nothing from the queue for this upload. `-d` empty. | The bytes are a package that has already been scored, or the scored combination unchanged. An unscored edit stops the upload. |
| Nov 4 20:00 | `vcc submit` has started. This is the deadline in the 09-29 amendment, not a target. | One in flight. Byte-range check at 21:00, 23:00, then Nov 5 08:00 and 12:00. |
| Nov 5 12:00 | Retry only if the Nov 4 upload died or stalled and was restarted. Same bytes, same reserved name. | A 7.5 h upload started at 12:00 finishes at 19:30, with margin before 23:59. Do not start a new arm. |
| Nov 5 23:59 | Deadline. The last upload that completed is the prize entry. | Nothing is submitted after a start time that cannot finish. |

HCASMC stays sealed until this candidate is frozen. If it is opened, it is one confirmation run and it may only choose among packages that already have members. It does not authorise a new build on Nov 4 (finding 15).

### Model names

Measured 2026-09-29 with `next_model_name.py --list`: queue entries 1–22 are used, 23–56 are free (34 names). The reserved final name is not in the queue and the helper will not hand it out (`docs/MODEL-NAMES.md:33`). R0 and every feedback arm take the next free name at the moment the command is written, via `--claim`, so two arms never share one. The Nov 4 upload uses the reserved sentence and no other. Names are not invented and not reused. `-d` stays empty on every submit in the window.

## 4. Stress

**A 2.8 GB upload that takes 7.5 h.** The 09-29 upload of 2.77 GB ran 03:08Z to about 10:40Z. Under that rate, R0 started at 15:30 on Oct 22 finishes about 23:00 and Oct 22's second slot is gone. A slot-1 start at 11:00 finishes at 18:30, so slot 2 of the same day is also gone. The sustainable rate is one new package a day, with the next arm prebuilt during the crawl. The plan's "two slots a day" (`docs/FINAL-PHASE-FEEDBACK-PLAN.md:46`) and the runbook's "15 min upload, do not start after 21:00" (`docs/DEF-RUNBOOK.md:707-710`) both describe a different network. The amendment at `docs/DEF-RUNBOOK.md:719-724` is the one this timeline follows. Counted arms if every row fires: R0, two ρ\*, two q, two row-4, two w, programme, gate, two scales, combination, two ρ\* re-checks: 16 queue names, plus the reserved re-upload. Oct 23 through Nov 3 is 12 days. At one slot a day, three arms do not fit. Drop from the bottom of the order (the ρ\* re-check first, then the second scale) and keep R0, the ρ\* pair, and the combination. The combination has to be scored before Nov 4 12:00 or it cannot be the last entry.

**A reboot mid-build.** AGENTS section 6.2 order, applied to this day: (1) compare the builder sha with `git show HEAD:scripts/make_submission_pm_cd4.py`. If `pin-artifacts.json` and its before-copy exist, `--restore --record` that file, not `pin-lodo.json` (restoring the LODO record leaves the first pin). If the machine died after `tmp.replace` and before the JSON, there is no record; `git checkout -- scripts/make_submission_pm_cd4.py` still works because the pin was not committed. (2) Re-arm waiters only if one was left; this timeline does not depend on a waiter to submit. (3) Board snapshot, so a published R0 is captured before the next submit replaces it. (4) Move the interrupted work dir aside (`build_X.aborted-…`) and rebuild in a fresh dir. Do not resume in place: `def_day.sh` will treat a pre-existing h5ad incorrectly if the dir is reused, and a fresh run whose on-disk builder is already the pin will restore the pin and call it success (526 finding 3). (5) Relaunch anything that was watching, with a resume note. Then re-check the builder sha before the next `def_day`.

**A `vcc` CLI upgrade in the window.** The day page's first command is `uv tool upgrade vcc-cli` (`docs/DEF-RUNBOOK.md:33`). Section 1 says not to upgrade in the window without re-running the rehearsal (`docs/DEF-RUNBOOK.md:251`), and section 6.7 says a mid-window upgrade re-runs the prep row first (`docs/DEF-RUNBOOK.md:731`). This timeline upgrades on Oct 21 or not at all. If the version changes anyway: no submit until `vcc prep --help` still has `--contexts` and a dry-run on the already-built rehearsal artifact exits 0. An upgrade during an in-flight upload is left alone until that upload finishes; the CLI is not restarted under it.

**Each bundle surprise in 526's question-1 table** (`experiments/526-cloud-october-review2/RESULTS.md:81-99`):

| 526 row | On the day |
|---|---|
| Context labels other than D/E/F, three identifier names, files and manifest agree | `def_day.sh` proceeds (`:137-141`). The hand `vcc prep` / `vcc submit` lines in the runbook still say `--contexts D,E,F`. Use the manifest's names or the prep rejects the file. Finding 7. |
| Context count ≠ 3 | Preflight fails (`scripts/def_preflight.py:189`). Stop. |
| Header missing or the wrong column name | Preflight fails (`:204-205`, `:219`). Stop. |
| `.var` index ≠ `gene_names.csv` | Preflight fails, and the builder fails again. Stop. |
| Duplicate targets | Preflight fails (`:223`). Stop. |
| `obs['context']` ≠ the filename | Preflight fails (`:126`). Stop. |
| `ntc_id` present, category count ≠ manifest | Preflight fails (`:137`). Stop. The count is categories, not codes used. |
| `ntc_id` absent | Preflight passes and `C_p` becomes the random-draw proxy (`scripts/make_submission_pm_cd4.py:5256-5261`). The package looks normal. Do not ship it. Finding 5. |
| `cells_per_pert` follows the manifest, including a value other than 400 | Preflight requires the key (`:187`) and does not require 400. `vcc prep` fails after the build. Finding 14. |
| Panel length ≠ 300 | LODO's `assert len(panel) == 300` (`experiments/438-october-stack-repricing/lib438.py:205`) stops R0, because LODO is on. A different 300 symbols is allowed by that assert. |
| Gene axis reordered, bundle path `data/challenge` | LODO's `assert official == OFFICIAL` holds, because both sides read the same path (`lib438.py:203`). The CD4 `V5_abc` basis then raises (`scripts/make_submission_pm_cd4.py:3202-3210`). The runbook says proceed (`docs/DEF-RUNBOOK.md:696-699`). Do not proceed. Finding 6. |
| `pert_counts.csv` columns swapped, header intact | Preflight and the LODO bind read `target_gene` by name. Continue. |
| No same-line bank | Valid. Same-line stays off. Continue. If the programme source is on, unmatched contexts are still switched on at w 0 (`scripts/def_sameline.py:250-259`) provided the source clears 150 rows and coverage 0.5. Below that floor the source raises and, if it was part of R0, the day stops. Finding 9. |
| A context matched to a bank already in the donor mix | Continues at w 0.5, which is the decided value, not `auto`. Row 5 is then in play. |
| A registered line whose call is not unanimous | That context stays off. Not a crash. Row 5 runs only if some other context did switch on. |
| More than 5 targets in no production donor | Census stops the day (`scripts/def_census.py:69-70`, `:237-239`). The runbook says re-check symbols and record it. There is no second command once the symbols are confirmed. Finding 16. |
| Stale whitener left in the work dir | `verify_existing_operator` refuses (`scripts/def_pipeline.py:425-450`). The runbook's warning that this skip is silent (`docs/DEF-RUNBOOK.md:308-317`) is stale. Delete the two files and refit. Do not edit the driver on the day to "fix" a closed hole. Finding 13. |

**Members hidden versus returned.** The first terminal D/E/F status is the branch. Returned members: capture, then the reader, one row at a time, and never a DE adoption under 0.009 from a single entry unless the row is prior-only (6 and 7). Hidden members: section 6 of the runbook, and the reader exits 2 if it is pointed at that capture. A published entry with null scores is a third state the capture script currently rejects (finding 8); it is the dark branch, not a retry loop. Mixed across our own entries (R0 returns members, a later entry does not, or the reverse) : the dark entry is not a read. Do not let it replace a package whose members you have. The last upload is chosen only from packages with members.

**The queue and the reserved name.** 34 free names as of this writing. The final upload does not consume one. Validation slots from Sep 29 through Oct 21 can consume the rest before R0 exists (finding 12). If `next_model_name.py` prints `queue exhausted`, append screened names; do not invent one in the submit command, and do not fall back to the reserved sentence for a mid-window arm.

**The last-entry rule under the stress cases.** The Nov 4 20:00 upload is the best package that has members at 16:00 that day. A 7.5 h upload started then finishes Nov 5 about 03:30. A reboot during that upload resumes the same bytes; it does not start a new build. A CLI upgrade does not change the bytes. A bundle surprise that invalidates R0 is discovered on Oct 22, not on Nov 4. An unscored combination, an unscored ρ\* tweak, and an HCASMC-inspired edit are all out.

## 5. Findings

1. **High.** Two slots a day does not survive a 7.5 h upload, and the full feedback list is about 16 arms against 12 days before the final upload has to be chosen. Evidence: `docs/DEF-RUNBOOK.md:719-724` (the 03:08Z–10:40Z upload and "one in flight"), `docs/FINAL-PHASE-FEEDBACK-PLAN.md:46` (two rows can share a day), `docs/DEF-RUNBOOK.md:707-710` (the 15 min budget those two rows were written against). Fix: one new submit a day unless the previous upload has already published by 16:00 UTC; prebuild the next arm during the crawl; if by Nov 1 the unread list is longer than the days left before Nov 4 12:00, drop from the bottom of section 2 (ρ\* re-check, then the second called-gene scale) and keep the combination scored.

2. **High.** A reboot after the pin replaces the builder and before the record is written leaves a pinned builder and nothing for the runbook's restore command to read. Evidence: `scripts/def_pin_panel_artifacts.py:207-225` (replace, then write the JSON and the before-copy), `docs/DEF-RUNBOOK.md:157` (restore needs the record), `experiments/526-cloud-october-review2/RESULTS.md:49-54`. Fix: apply 526's `patches/pin_record_before_replace.patch` before Oct 22 (write the copy and the JSON, then replace).

3. **High.** A fresh `def_day` started on that pinned builder records the pinned sha as `BSHA0` and the exit trap restores the pin, exit 0. Evidence: `scripts/def_day.sh:67` and `:100-120`, `experiments/526-cloud-october-review2/RESULTS.md:54-58`. Fix: apply 526's `patches/def_day_refuse_pinned.patch` (refuse at start when the on-disk builder is not `git show HEAD:scripts/make_submission_pm_cd4.py`). Until that lands, the Oct 22 recovery checklist checks the sha by hand before every `def_day`.

4. **High.** The default exit trap restores the builder before the commit the runbook requires, and the success log tells the operator to re-run the pin stages to get the diff back. Evidence: `scripts/def_day.sh:354-357`, `docs/DEF-RUNBOOK.md:328` (commit the pins before submit). Fix: the submit path exports `DEF_DAY_KEEP_PINS=1`, commits the reviewed diff, then `--restore --record`. The timeline above does that. The script's default should say so in the day page, not only in a comment at `:27-28`.

5. **High.** A control file with no `ntc_id` column passes preflight and the build takes the random-draw `C_p`. Evidence: `scripts/def_preflight.py:133-137` (the check runs only when the column exists), `scripts/make_submission_pm_cd4.py:5256-5261`, `experiments/526-cloud-october-review2/RESULTS.md:90`. Fix: when `per_context[c].n_ntc_ids` is set and the column is absent, preflight fails. Until that patch, Oct 22 treats a missing column as a stop even if preflight prints PASS.

6. **High.** The runbook says a moved `gene_names.csv` proceeds. The pinned CD4 basis raises on a gene-order mismatch, and LODO's same-path assert does not catch it. Evidence: `docs/DEF-RUNBOOK.md:696-699`, `scripts/make_submission_pm_cd4.py:3202-3210`, `experiments/438-october-stack-repricing/lib438.py:203-205`, `experiments/526-cloud-october-review2/RESULTS.md:92`. Fix: rewrite that runbook bullet to stop. A changed axis is a different panel for `V5_abc`; do not ship and do not refit the basis on the day.

7. **High.** Hand prep and submit hardcode `--contexts D,E,F`. The driver reads the manifest and accepts any three identifier names. Evidence: `docs/DEF-RUNBOOK.md:201`, `:327-329`, `:729`; `scripts/def_day.sh:137-141`; `experiments/526-cloud-october-review2/RESULTS.md:83`. Fix: every hand command takes the three names from `manifest.json`. A label set other than D/E/F is a proceed, not a stop, once those commands match the driver.

8. **High.** `capture_board_entry.py` will not write a published entry whose scores are null, so "members hidden" and "capture failed, retry" are the same exit code. Evidence: `scripts/capture_board_entry.py:79-84` (published and superseded require all twelve score fields), `docs/DEF-RUNBOOK.md:732-741` (the branch is "members come back" versus "nothing comes back"), `docs/FINAL-PHASE-FEEDBACK-PLAN.md:5-6`. Fix: a terminal status with null scores writes a capture marked members-hidden and exits 0. The reader already treats that file as exit 2 and adopts nothing. Retry only while the status is non-terminal.

9. **High, conditional on the programme being in R0.** A programme source with fewer than 150 rows left or coverage under 0.5 raises, and `def_day` dies when the artifact switches on zero contexts. That blocks the first upload, which the feedback plan is not allowed to delay. Evidence: `scripts/def_sameline.py:245-249`, `scripts/def_day.sh:168-169`, `docs/FINAL-PHASE-FEEDBACK-PLAN.md:15-16`. Fix: if R0's resolved command includes `DEF_DAY_PROGRAMME_SOURCE` and the sameline stage dies for that reason, rerun immediately in a fresh dir with the source unset and flags `noprog`. Record it. Do not debug the bank on Oct 22.

10. **High, conditional on per-target sharpening being in R0.** Row 3's positional q is ignored when `DEF_DAY_SHARPEN_PT` is set. Evidence: `scripts/def_day.sh:93-98` ("scalar q not used"). Fix: row 3 runs only if R0 left sharpening off. If R0 shipped `stercile:…`, changing the positional q is a no-op and unsetting the variable is a second change. Skip row 3 in that branch rather than submit a two-variable arm. The command table says the same.

11. **High.** Rows 2, 4b and 12 have no knob. ρ\* is a fixed formula. The dump is hardcoded H1 ×2 / CD4 ×1, and the recipe refuses any other CD4 stack. Evidence: `scripts/def_day_helpers.py:28-32`, `scripts/def_day.sh:194-197`, `scripts/october_recipe.py:76-77`. Fix: land experiment 530's `DEF_DAY_RHO_SCALE` and `DEF_DAY_CD4_STACK` before Oct 22, with unset bit-identical to today, and a meta key so the sibling diff can see the change. The commands in section 2 are written as they will be after that. Do not hand-edit factors or `VCC_CD4_TARGET_G` on a level build in the meantime.

12. **Medium.** The called-gene scale (rows 8-9) and the fixed LODO cap of 6 (row 4a) are also missing, and they were not named as 530's work. Evidence: no called-gene scale in `scripts/make_submission_pm_cd4.py` (the zs rider is the only `VCC_DE_MAG`); `scripts/def_day.sh:179-183` passes `--hi-pct` only; `scripts/def_lodo_key.py:61` defaults `--hi` to 3. Fix: a `DEF_DAY_DE_CALLED_SCALE` pass-through and a `DEF_DAY_LODO_HI` pass-through that cannot be combined with `--hi-pct`. If either is absent on its day, skip the row. Unset hi-pct is cap 3, not the "fixed 6" the plan is asking to compare.

13. **Medium.** The runbook still tells the operator the whitener skip is silent, and section 2 step 10 still diffs against the h1union meta. `def_day.sh` diffs against 446b and the skip is closed. Evidence: `docs/DEF-RUNBOOK.md:305` versus `scripts/def_day.sh:314`; `docs/DEF-RUNBOOK.md:308-317` versus `scripts/def_pipeline.py:425-450`. The one-page says it supersedes section 2's command block where they differ (`docs/DEF-RUNBOOK.md:13`), but those two paragraphs are the ones an operator will follow when something looks wrong. Fix: strike the skip warning and point step 10 at 446b's meta and at `def_day.sh`'s allow-list.

14. **Medium.** Preflight does not require `cells_per_pert == 400`. A manifest that says something else fails at `vcc prep`, after the build. Evidence: `scripts/def_preflight.py:187` (the key is required; the value is not pinned to 400), `experiments/526-cloud-october-review2/RESULTS.md:91`. Fix: one preflight row, `cells_per_pert == 400`, before the whitener.

15. **Medium.** Opening HCASMC "while there is still time to act" conflicts with the last-entry rule. Evidence: `docs/DEF-RUNBOOK.md:464-469`, `docs/FINAL-PHASE-FEEDBACK-PLAN.md:17-19`. Fix: one sentence in the runbook. A disagreement can only select an already-scored package. It cannot start a new build after the candidate is frozen.

16. **Medium.** Census stops the day above 5 zero-carrier targets, and the runbook's response is "re-check symbols" with no command for the case where the symbols are right. That can push R0 past the upload window. Evidence: `scripts/def_census.py:67-70` and `:237-239`, `docs/DEF-RUNBOOK.md:690-692`. Fix: after the symbol check is written down, a documented one-shot override in a fresh dir, used only for R0, with the count recorded. Do not re-weight donors.

17. **Medium.** `def_postbuild.py` cannot express "these keys moved and no others" ; it expresses "nothing outside this set moved". A silent no-op passes. Evidence: `scripts/def_postbuild.py:465-469`. Fix: the feedback gate in section 2, run by hand until the script grows a `--require-config-diff`. The equality check is the gate, not the exit code alone.

18. **Medium.** ρ\* ×1.02 aims at `0.9989303453`. The loop deletes any build with a context past ρ 1 (`scripts/def_day.sh:238-241`) and will then ship the nearest safe build, which can be R0's target rather than the arm. Evidence: the target computed from `relevel_iter.py:28-29`; the 1% band around 0.99893 crosses 1. Fix: the row-2 gate in the command table. Do not submit a `nearest safe` artifact as the ×1.02 arm.

19. **Medium.** 34 queue names are free and nothing reserves a block for the final window. Evidence: `docs/MODEL-NAMES.md:33` and the queue that follows; `scripts/next_model_name.py` (measured 2026-09-29: 1–22 used, 23–56 free). The reserved sentence is protected. The other 34 are not. Fix: append screened names before Oct 21, and do not let validation claims consume the last 20. A mid-window `queue exhausted` is a stop, not a prompt to type a name.

20. **Low.** `def_day.sh` cds to `~/Code/virtual-cell-challenge-2026` (`scripts/def_day.sh:58`). That is the production clone. A worktree or this VM dies on the first line. Fix: derive the repo from the script path before Oct 22 if there is any chance the day runs anywhere else. On the production Mac, no change is required.

21. **Low.** The stall response is not written down. The amendment says to check the byte range and does not say what a stuck range means. Evidence: `docs/DEF-RUNBOOK.md:722-724`. The 09-29 upload was slow and then recovered, so a slow-but-moving range must be left alone. Fix: the 30 min unchanged-range rule in the timeline. Restart only then, and only the same bytes.
