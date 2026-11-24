# 499 — adversarial review of the October day-of code changes of 2026-09-27 (cloud agent, code only)

You are reviewing code in the repo `Ahmet-Dedeler/virtual-cell-challenge-2026` (branch `master`). Read `AGENTS.md` first
(stable rules), then `docs/DEF-RUNBOOK.md` (what the October D/E/F day runs). You have NO gitignored data on this VM
(no challenge controls, donors, harness files): work from the code and synthetic inputs only. Never run `vcc`, never
touch credentials. Commit ONLY to the branch `cloud/499-october-code-review`, ONLY inside
`experiments/499-october-code-review/`. Do not open a PR.

## Why this matters
On Oct 22 the organisers release three new cell lines (D/E/F); the leaderboard goes dark, and `scripts/def_day.sh`
builds the final submission unattended. A defect there ships a wrong or stale build that no score will reveal. The
production Mac runs **macOS `/bin/bash` 3.2** (no associative arrays; `"${arr[@]}"` on an empty array under `set -u` is an
"unbound variable" error in bash < 4.4; no `${var,,}`, no `mapfile`), **BSD grep/sed** (`sed -i ''`, no `\|` alternation
after a `$` anchor), and zsh for interactive shells. This VM has GNU bash 5.x: flag bash-3.2-only issues as PLAUSIBLE
unless you can reproduce them (e.g. with a bash-3.2 container if available).

## What changed today (review all of it)
1. **476: 485's review fixes merged onto 474.** The merged script is `experiments/476-def-day-merge/work/def_day.resolved.sh`
   (to replace `scripts/def_day.sh`, which is 474's version at HEAD). Inputs: 474's script (`scripts/def_day.sh`), the
   pre-474 base (`experiments/476-def-day-merge/work/def_day.base.sh`), 485's patched script
   (`experiments/476-def-day-merge/work/def_day.theirs.sh`), 485's findings and patches
   (`experiments/485-cloud-october-review/RESULTS.md`, `patches/`). The merge notes are in
   `experiments/476-def-day-merge/RESULTS.md`. 485's other three patches (helpers, `scripts/def_pipeline.py`,
   `experiments/448-def-rehearsal-q175/relevel_iter.py`) apply unchanged; the patched copies are in
   `experiments/476-def-day-merge/repro/*.new`, and `repro/repro_fixes.py` passes on them.
2. **497: two October pass-throughs** on top of 476: `experiments/497-october-switches/patch_defday497.py` (applied copies:
   `work/def_day.497.sh`, `work/def_day_helpers.497.py`; unit test `test_helpers497.py`). `DEF_DAY_LODO_HI_PCT` passes
   `--hi-pct` to `scripts/def_lodo_key.py` (committed today) and pins the key as `lodo_lo0.5_pct_def`;
   `DEF_DAY_SHARPEN_PT=<tercile|stercile>:q,q,q` exports the builder's per-target sharpening form
   (`VCC_BULK_SHARPEN_Q`, parsed in `scripts/make_submission_pm_cd4.py`, search `BULK_SHARPEN_PT`) for the level builds,
   and `meta_matches` / the postbuild allow list are extended.
3. **496: the emitter/builder patches for an optional `mde_add_sign`** (not applied yet):
   `experiments/496-uncalled-sign-hook/patch_emitter496.py` (on `src/vcc2026/rc_emit_k.py`) and `patch_builder496.py`
   (on `scripts/make_submission_pm_cd4.py`; env `VCC_DE_MDE_SIGN=clean4:<tau>`).

## What to look for
- (a) any path where the day ships a build other than the one it measured (stale h5ad, wrong factors, wrong q or
  sharpening form, wrong key), or links `_final` to an artifact whose meta disagrees with `def_day.final.env`;
- (b) any path that stops the day when a safe build exists, or continues when it should stop;
- (c) bash 3.2 / BSD-tool incompatibilities, `set -u` traps, quoting and word-splitting errors;
- (d) whether every new knob is inert when unset (same exports, argv, checks, outputs);
- (e) the merge itself: anything of 474's rehearsed behaviour or 485's fixes lost or duplicated in the resolution;
- (f) 497: the percentile-cap pin id and key file naming, the per-target form through meta_matches, final.env, the 5b
  allowed diff; interactions with `DEF_DAY_PRIOR_GATE`, `DEF_DAY_PROGRAMME_SOURCE`, same-line mode;
- (g) 496: is the emitter byte-for-byte unchanged in behaviour when `mde_add_sign` is absent; is the builder knob's gene
  axis (`rc_control.keep`) right; are the pins and the config key right.

## Deliverable
`experiments/499-october-code-review/RESULTS.md`: one section per finding with severity (critical / high / medium / low),
file:line, what goes wrong on the day, CONFIRMED (reproduced; include the command and output) / PLAUSIBLE / REFUTED, and a
minimal fix as a patch under `experiments/499-october-code-review/patches/` (against the file it names, `diff -u` from the
repo root). Put reproduction scripts under `repro/`; they must run without the gitignored data (synthetic inputs, stub
drivers). Lead with a table of findings ranked by severity. Say explicitly which areas you checked and found clean. Be
concrete and terse; no speculation presented as fact.
