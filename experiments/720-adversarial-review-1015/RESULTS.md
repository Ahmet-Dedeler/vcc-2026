# 720 — Adversarial review of the Oct 22 → Nov 5 plan (2026-10-15, session [020377], fifth lane of the 10-15 council)

Class: day-of procedure review. No recipe value changes; nothing here moves a score. Reviewed as of HEAD d8edf4f8 plus the
readiness lane's uncommitted `experiments/690-def-readiness-1015/` and `691-release-watcher/`. Targets, as the coordinator
set them: the October recipe as `scripts/october_recipe.py command` resolves it, `docs/DEF-RUNBOOK.md`, 690's
`OCT22-CHECKLIST.md`, `docs/FINAL-PHASE-FEEDBACK-PLAN.md`, and our public footprint. Primary sources for the external facts
are in `docs/NEWS-2026-10-15.md` (site JS, rules, FAQ, Arc's 2025 and 2026 emails, Discord).

## 0. What resolves today
`uv run python scripts/october_recipe.py command --bundle data/challenge --work <w> --name <n>` (no `--allow-pending`) now
resolves: every switch is DECIDED (the coordinator decided rho_scale / lodo_strength / lodo_hi / callset_block off at
22:33Z). It prints:
`DEF_DAY_LODO=1 DEF_DAY_LODO_HI_PCT=96.64697522769474 DEF_DAY_SHARPEN_PT=tercilec:2.5,2.0,1.75 DEF_DAY_DE_MAG=dn_k250_zs
DEF_DAY_MDE_SIGN=clean4:0.05 bash scripts/def_day.sh <bundle> <work> <name> 1.75 0.5 noprog 4`; `check` rc 0.

## 1. Findings, severity-ordered (S1 could lose the round; S2 could cost a day; S3 housekeeping)

**F1 (S1) — The planned last move is refused by our own tooling.** FINAL-PHASE plan §1 rule 2 and DEF-RUNBOOK §6 (09-29
amendment) make the Nov 4/5 entry "the best scored package re-uploaded unchanged". `scripts/arm_submit.sh:20` runs
`check_submission_novelty.py check`, which exits 1 for any package whose `.vcc` sha OR prediction digest is in
`submissions/SUBMITTED_LEDGER.json`, and its docstring says bypassing is "deliberately impossible" (only
`--allow-unidentified`, for non-builder packages). So on Nov 4 the armed path says "novelty check refused -- NOT submitting".
Fix options (coordinator's call): (a) a documented manual sequence in DEF-RUNBOOK §6 — `vcc submit <vcc> -m "<reserved
name>"`, then `check_submission_novelty.py record … --note final-reupload`, then `capture_board_entry.py`; or (b) an explicit
`FINAL_REUPLOAD=1` branch in `arm_submit.sh` that skips `check`, logs the override and still records. Either way, rehearse it
once in the final phase on an entry we would re-upload anyway (see F2).

**F2 (S1) — Nobody knows whether the server accepts a byte-identical, already-scored file.** Our rule "never re-submit a
scored package" means it has never been tried. The CLI sends `model_name` + `file_name` and a resumable upload; no client-side
dedupe. If the server dedupes by content or file name, the last-entry plan fails on the evening it matters. Cheap insurance:
re-upload R0's `.vcc` as the second entry of a quiet final-phase day under a fresh queue name and record what happens. If
refused, the fallback is a re-prep (same prediction, different `.vcc` bytes), which `vcc prep` produces deterministically.

**F3 (S1) — Final-phase feedback: the prior should be DARK.** In 2025 Arc wrote (10-31-2025) "once the final test set is
released, your submissions will receive a confirmation email and a success notification in the app, but you will not see a
score", and (11-10-2025) "your score will remain hidden until winners are announced". The 2026 rules keep "no leaderboard
displayed for final entries"; the FAQ is silent; the gated `/api/leaderboard/team/<id>?get_final=true` may only serve the
post-announcement view. The feedback plan is correctly conditional and `capture_board_entry.py` already writes
`members_hidden` on null members (exit 0). Action: keep the dark plan (DEF-RUNBOOK §6) as the default mental model; do not
schedule probe arms before the first `vcc status` read.

**F4 (S2) — The release watcher gives up too early and cannot upgrade the CLI.** `watch691.sh:73-78` retries the download
6 × 60 s, then `STOP … exit 1` and AUTO_R0 never runs. Cold Turkey's morning block (~03:30–07:40Z on 10-05 and 10-06)
overlaps the 07:00Z release and lasts longer than 6 min. And 2025's release-day email required a tool upgrade ("you'll need to
upgrade to the latest version of the cell-eval Python package"); a `vcc-cli` release around Oct 22 is likely, and a stale CLI
may not list or download the new dataset. Fix: retry with backoff for ≥ 3 h; after N failures run the standing-approved
`uv tool upgrade vcc-cli && vcc skill install` once and retry; notify on every failure, not only the last.

**F5 (S2) — The watcher picks the FIRST new dataset id** (`watch691.sh:70`). If Arc publishes two datasets together (a
final bundle and, say, the A/B/C ground truth, as a 2025-style release might), the wrong archive is downloaded; `release691.py`
unpack refuses only on structural grounds (one manifest, no traversal). Prefer the id whose filename or description matches
`final|test|DEF`; with two candidates and no match, stop and leave NEXT.txt for a human.

**F6 (S2) — Every submission after the `final_week` gate flips is a FINAL entry.** The CLI has no flag (`create_submission`
sends only model_name/file_name; `is_final` is server-set), and the gate may flip at 00:00Z while the bundle lands at 07:00Z.
An A/B/C package submitted on Oct 22 would become our standing final entry and spend a final slot. Rule: nothing is submitted
on Oct 22 before R0. (The coordinator's 10-21 slots are the last validation slots.)

**F7 (S2) — 690a re-opens PENDING switches.** Today's recipe has none; applying 690a adds `de_uprn`, `mde_sign_keygate`,
`de_dnmag` as PENDING, and `command` refuses again until each is decided (off is a decision). Decide them at apply time, not at
the 10-21 freeze.

**F8 (S3) — Submission cadence assumptions.** 2025 allowed one final entry per day; the 2026 rules say two. FINAL-PHASE §5's
"one new arm a day" is still the right plan given uploads of up to 7.5 h and 12.5 h scoring queues (10-06; Arc's #challenge-
status says capacity was expanded 10-07).

**F9 (S3) — Public footprint.** All 50 of our board model names since 09-02 are quotes or jokes with no method word (list in
`work/`); the only leaks remain the known 08-22 "Attention is not all you needed" and 08-24 "Alpha 0.40: probing the mse
floor". `-d` has been empty since 08-27. The Arc blog mini-profile (10-15) is Ahmet's text; the coordinator reviews it when
live. The Jason Levine trade offer (10-03) is unanswered; answering it would be a §14 and rules risk.

**F10 (S3) — Timing facts the checklist already carries** (verified against the site's JS): release 07:00Z Oct 22; the app's
final-week countdown targets 23:59Z Nov 5 while the static widget says 08:00Z; "Safe deadline Nov 5 08:00Z" is the right
conservative choice. Team members can be added only until the release.

## 2. Things checked and found sound
- `capture_board_entry.py` reads `vcc --json status` (token auth), not the public leaderboard, so final entries are
  capturable; `members_hidden` and the final-phase `rank` exemption are handled (690d).
- `def_day.sh` carries pass-throughs for every feedback-plan row that is still live (RHO_SCALE, CD4_STACK, LODO_HI,
  LODO_STRENGTH, CALLSET_BLOCK, RC_CALLSET_CPM, DE_UPRN, MDE_SIGN_KEYGATE); rows 8–9 (called-gene scale) have no knob and are
  already marked "skip".
- The watcher never touches `data/challenge` (A/B/C stays for the other lanes); the work dir rule (under `experiments/`) is
  recorded in 691.
- The model-name queue has ~55 unused names plus the reserved final name; enough for the final phase.
- Scorer: no public change; `anchor_version` r4 on every snapshot; cell-eval2 0.18.0 moves no number.

## 3. Not reviewed yet
`def_day.sh`'s own stage order against a SAME-PANEL D/E/F (690's R-SAME covers it); `def_postbuild.py`'s allow-list for the
final re-upload (no build, so none needed); the 2025 "final results leaderboard" timing vs. the 2026 "late November".

## 4. Second pass (2026-10-15 ~05:10Z) on 691 v2's RESUME path, the 720 fixes (1ba3851e, f9024c8d, 0758daa3) and arm1015f
Verified sound: `vcc datasets download -d DIR` is a real 0.2.2 flag; `release691.py diff` catches an in-place replacement of
`controls` (filename/description change, or unavailable → available); the pin check reads `builder_sha256_after`, which
`def_pin_panel_artifacts.py` does write; `def_day.sh` logs "builder after restore:" from its EXIT trap on every path, and 692
turns a STOP into NOT READY, so a stopped R0 cannot be prepped; arm1015f's FINAL_REUPLOAD branch parses the novelty refusal's
"  entry <id>" line and MODEL-NAMES' reserved-name line as they are written today; the branch fails closed without a final-phase
capture of the ledgered entry.

Findings (S2 = could cost the day's hours, S3 = housekeeping):
- **R1 (S2) No "already running" guard.** `watch691.sh:45` overwrites `watch691.pid` without checking it. Two live watchers (a
  human re-launch while the first is alive) both pass the R0-dir checks during the gate wait (no def_day yet), and the gate's
  "day-path chain" test does not see a watcher, so both can enter `def_day.sh` on the same `$R0W` seconds apart. Add at start:
  if the recorded pid is alive and its args contain `watch691.sh`, STOP.
- **R2 (S2) Early release defeats the baseline.** The baseline is the catalog at arm time (10-21 ~21:00Z). If Arc lists the
  D/E/F dataset as `available: true` before the arm, `diff` never fires and the watcher polls to CAP. Checklist §B should say:
  read the baseline line in `watch691.log`; if it already lists more than `controls`, do not arm, download by hand.
- **R3 (S2) CAP on a late re-launch.** The RESUME branch skips the CAP check (fine), but the download loop breaks at CAP after one
  failure (line 119) and the gate wait exits 3 at CAP (line 199). A re-launch after 10-23 18:00Z therefore dies on the first
  hiccup. Checklist §B: any re-launch after CAP passes `CAP691=<new>`; `AUTO_R0=1` must be repeated on every re-launch (env,
  not persisted).
- **R4 (S2) Automatic CLI upgrade needs a follow-up step.** `upgrade_once` fires on any error text matching
  `upgrade|version|unsupported|426` (line 74) — "version" appears in unrelated TLS/Python errors. Harmless once, but the day's
  real `vcc prep` then runs on a new CLI: DEF-RUNBOOK §6.7's prep known-answer must be re-run when `watch691.log` has an
  UPGRADE line. Tighten the regex to `upgrade|unsupported|426|client version`.
- **R5 (S3) A stopped R0 is silent.** `watch691.sh:162` logs "already holds a finished def_day … nothing to do" and exits 0 for
  a STOP'd run too (same trap marker). 692 will write NOT READY, but no notification wakes anyone. Grep `STOP:` there and
  `notify "R0 STOPPED"`.
- **R6 (S2) The identical-bytes test can burn hours on a deterministic refusal.** If the server refuses a re-upload before an
  entry exists, arm1015f's retry loop (lines 97–111) treats it as "failed before any entry" and retries with backoff until CAP.
  For `FINAL_REUPLOAD_TEST=1`, cap the retries at 2 so the refusal is read in minutes and the fallback (re-prep → new `.vcc`
  bytes, same prediction digest, which the branch also accepts) can be scheduled the same day.
- **R7 (S3) Pick regex.** `test` in `FINAL_RE` also matches a hypothetical "validation ground truth (test …)" release; the
  ambiguity path then STOPs correctly, but a single such dataset with no final bundle would be downloaded (WARNING only) and
  stopped only at summarize's A/B/C guard. Acceptable; noting it.
Everything else in the RESUME path (interrupted-dir move-aside before the gate, pin-held STOP with the restore command, finished
→ defer to 692, PICK691 resume, fresh `$DL/$TS` reuse for the CLI's own download resume) reads correctly.

**Closure (05:50Z):** [17adab] fixed R1–R6 in bde82e21 (R7 noted). Verified by reading the scripts, not the message: watch691.sh
STOPs when the recorded pid is alive and its args contain `watch691.sh` (pid-specific, not a `pgrep -f` pattern; `trap - EXIT`
avoids a false END line); r0ready692.sh has the same guard on `<stem>.r0ready692.pid`; the list-refusal regex is
`upgrade|unsupported|426|client version`; an UPGRADE logs and notifies that the §6.7 prep KA must be re-run; a STOPPED R0
notifies and exits 1 on both the resume and the own-build path; arm1015f caps the submit at 2 tries under
`FINAL_REUPLOAD_TEST=1`. Both production processes re-armed on the fixed scripts (pids 62641 / 62642, 23:20Z) with a
controls-only baseline. R2/R3 live in the checklist §B. 720 is closed; anything new goes in a 72x.

## 5. 694 R-ABC read-readability pass (10-15 21:45Z, asked by the coordinator for the 10-16 s1 slot)
Verdict: READABLE as preregistered (a37cbcae). Evidence: base 524dr2 captured at full precision (dVHLws6KdMxudq0Yn5ap);
on the 446bj2z path both riders moved pds/mse/fid/jac by exactly 0.0 raw (630arm−446bj2z; DN−630arm) and reach by
0.0 / +9.6e-8, so F1 (1e-9), F2 (1e-6), F3 (1e-5) are backed by precedent (the hooks' ~1e-7 jackknife drift does not apply:
the riders keep per-cell libraries exact); ntclevel417's rho reads only per-block pseudobulk totals (decompose) and the NTC
constructs (cp46), so the re-level retraces 524dr2's factors and gate694's bulk_ctx_scale / bulkid658 checks fail closed
before a slot; prior-only noise ~1e-5 vs a NULL edge of +0.0004; prediction +0.0013 = ledger's +0.0009297 + 0.0003736.
Caveats: addends read on the 446bj2z path (SUBADDITIVE live; bins handle it); READY ~02:10Z and the read needs a normal
scoring queue; chain CAP 10-16 17:00Z.

## 6. Lab charter boundary review (10-15 ~23:20Z, docs/LAB-2026-10-16.md, for the coordinator)
Verified: repo private; gate602 checks memory/pressure/swap/disk/FIFO and `other_builders()` only, with stale Oct 4–7 windows;
.gitignore does not cover experiments/**/work/, *.parquet, *.ckpt; aws_reap only STOPS instances. Holes sent to the coordinator:
H1 git index (branch/stash/reset/rebase/merge-pull and `git add <dir>` in the shared checkout), H2 machine (gate blind to
def_day setup, `vcc prep/submit`, day-path chains; no Oct 22 window; gate-pattern names), H3 data dirs (downloads must not touch
data/*), H4 network (no > 1 GB transfers in slot/upload/bundle windows), H5 cloud cost (budget alarm, per-launch/daily caps,
`Project=lab` tags, reap volumes/snapshots, ring-fence production credits); housekeeping: reserve wear (R7 lines), Anthropic usage
concurrency, coordinator re-runs dossier KAs.

**§5 closure (10-16 06:50Z):** the 694 read landed as predicted, COMPOSES (+0.0012429 vs 524dr2; predicted +0.0013; bulk exact,
re-level factors identical, F1 IDENTICAL, postbuild 42/42). `gate694`'s strict diff nevertheless FAILED on `bulk_load_sha256` and
`sameline`, which are rebuilt per run (the day's LODO dump and same-line record carry run-specific shas even when their content and
every total are identical), so the gate produced a false NOT READY and the coordinator read the package by hand. Lesson for strict
config diffs (also for Oct 22's R0 gates): whitelist per-run artifact keys or compare content digests, never file shas, and keep the
behavioural checks (F1 bulk identity, bulk_ctx_scale, postbuild) as the real gates. The readability verdict was right; the gate's
key list was not. Also consistent with 721: K-UP (10-16 s2) read NULL on fid + jac (+0.00012); K stays 2,600.
**§5 resolved by 8fa198b7 (10-16):** gate694 v2 is content-aware and `scripts/def_metadiff.py` (tests 8/8) now serves the strict
diff in OCT22 §C 6a and rider690: per-run artifact keys pass when their arrays are equal after the work-dir token is replaced; a
changed array, provenance text or same-line w still FAILs. R0's work dir is kept through the final as the diff's reference.
