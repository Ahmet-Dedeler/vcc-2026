# 691 — the Oct 22 release watcher (2026-10-15, session [17adab]; coordinator-approved AUTO_R0, to be armed on 10-21)

Class: day-of procedure. Nothing fitted; no score moves.

## What it does (`watch691.sh`, helpers `release691.py`)
From START (default 10-21 22:00Z; Arc's site widget puts the final test set release at 10-22 07:00Z, [020377]) it polls
`vcc --json datasets list` every 5 min against a baseline taken at arm time. A new or changed available dataset is the
release. Then, unattended: download with retries (the CLI resumes and verifies the server's checksum) → unpack as
untrusted data into a fresh `data/challenge-def` (data/challenge stays A/B/C) → the driver's builder-free stages (a tripped
preflight guard does not stop the identity read) → summary (SAME-PANEL vs NEW-PANEL against A/B/C, gene axis identical,
identity calls) → `NEXT.txt` with the exact command for that branch. With AUTO_R0=1 (approved by [400b6e] 10-15, armed
on 10-21 after the final rehearsal) it also builds R0 with def_day.sh, but only if preflight, the light stages, the summary
and a plain `october_recipe.py command` (no `--allow-pending`) are all clean and the gates pass (no HOLD, builder == HEAD,
no builder or day-path chain, free+inactive ≥ 26 GiB, disk ≥ 60 GiB). On SAME-PANEL it pre-places 453's operator 66610be7.
It never runs the real `vcc prep` or `vcc submit`.

## Tests (stub `vcc_stub691.sh`, no network; records in `work/testrecords/`)
| test | bundle released | result |
|---|---|---|
| new | 690's foreign controls, ef-bundle-r panel | detected on poll 3; unpacked 6 files; preflight rc 1 (641 depth row, expected); identity/cd4/census rc 0; branch NEW-PANEL; NEXT.txt = the recipe command |
| same | 690's foreign controls, A/B/C panel | branch SAME-PANEL; NEXT.txt adds the operator pre-place line |
| defb | 180's def-bundle (A/B/C relabelled D/E/F) | preflight rc 0; identity Jurkat/HeLa/CAL33; SAME-PANEL; AUTO_R0=1 reached the gate, waited on "builder != HEAD" (the DN build held it), STOP at CAP, no work dir created |
| unpack refusals | crafted zips | `../` traversal, absolute path, symlink member, a `.sh` member, two manifests, existing destination: all REFUSED, nothing extracted outside the staging dir |

Finding for the day: `def_pipeline.py` records stage logs with `relative_to(ROOT)`, so a work dir outside the repo
crashes the driver after the stage ran. The day's work dirs must live under `experiments/` (they always have).

## End-to-end test
690's R-SAME (v2) runs this watcher with AUTO_R0=1 on a stub release of build690b's depth-matched SAME-PANEL bundle: the
whole automated path with no guard override, through def_day. See 690's RESULTS.

## v2 (2026-10-14 23:10Z): resume after a restart + review 720's fixes (records in `work/testrecords/v2_*`)
- **Resume** (AGENTS 6.2 step 6 re-launch): a release already detected skips the poll; a bundle matching its unpack record
  (`release691.py verify`) skips the download; the cheap stages re-run. An R0 dir with a finished def_day is left for
  692; a running one is left alone; an interrupted one is moved aside and rebuilt, unless the builder on disk is that
  dir's own pin (`builder_sha256_after` of its pin records), in which case it STOPs and logs the exact restore command.
  Another session's live build is not mistaken for our pin. Tested: resume with a verified bundle; (a) interrupted,
  not holding the builder → moved aside, gate waits; (b) its own pin → STOP + restore command; (c) finished → nothing to do.
- **720 F4**: the download retries with backoff (60 s doubling to 15 min) for 4 h, notifies on each failure, and after the
  3rd failure runs ONE standing-approved `uv tool upgrade vcc-cli` + `vcc skill install`, logging the version before and
  after (a uv tool, not the project venv: no numerics KA). A list refusal that names the version also triggers it. Tested:
  4 failures → upgrade after the 3rd → success on the 5th.
- **720 F5**: `release691.py pick` chooses among several new ids by final/test/D-E-F in the id, filename or description.
  Ambiguous → nothing is downloaded, NEXT.txt lists the candidates, and a re-launch with `PICK691=<id>` resumes with that
  choice. Tested: two final-looking ids → STOP; one final + one A/B/C refresh → the final one; PICK691 → resumes.
- **A/B/C guard**: a bundle whose contexts are A/B/C makes the summary rc 3, so AUTO_R0 never builds it. Tested with
  data/challenge zipped.

## LaunchAgent draft (2026-10-14 23:30Z; NOT installed, waits on Ahmet's yes via the coordinator)
`launchagent/`: `com.ahmet.vcc.watch691.plist` (RunAtLoad, StartInterval 1800, AbandonProcessGroup so the re-launched
processes survive the job's exit, PATH with ~/.local/bin and /opt/homebrew/bin) runs `relaunch691.sh`. Between 10-21 18:00Z
and 10-23 18:00Z it re-launches the watcher and the 692 adapter only if no copy is running AND the log does not end in END,
so a process killed by a reboot resumes and one that ended on purpose waits for a person. Dry-tested: outside the window it
does nothing; with both running it does nothing; killed mid-run → it would re-launch; ended or never armed → it does not.
Install and uninstall commands are in `launchagent/README.md`.

## Armed (production) 2026-10-14 23:16Z
watch691 pid 54923 (AUTO_R0=1; baseline: `controls` only; polls from 10-21 22:00Z, cap 10-23 18:00Z) and r0ready692 pid
54924 (cap 10-23 18:00Z). Commands in docs/PLAN-2026-10-15.md. Any edit of either script = kill + re-arm (they resume).

## v3 (2026-10-14 23:30Z): pre-decided overrides, frozen copies
- `release691.py overrides` + the watcher: with `AUTO_ALLOW_DEPTH691=1`, a uniform depth shift (ratios within 2 %, in
  [0.6, 1.8], depth rows the only failures) gets `DEF_DAY_ALLOW_DEPTH=1`. With `AUTO_ALLOW_WIDE691=1`, spreads ≤ 3.5 get
  `DEF_DAY_ALLOW_WIDE_LIBS=1`. Preflight re-runs with the override and must pass, and def_day gets the same env. Nine
  fixture cases plus an integration run on a uniformly thinned (×0.80) def-bundle (`work/testrecords/v3_depth`). Both
  flags default OFF.
- Both scripts copy themselves to `<work>/.frozen-<pid>-<name>.sh` and exec it, so editing the source never corrupts a
  live run. The single-instance guard and the relauncher match the basename. Tested.
