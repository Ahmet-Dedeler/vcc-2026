# Oct 22 checklist (v3.2, 690; 2026-10-15). One page. The guard-by-guard actions stay in 644's DECISION-CARD.md.

Times are UTC. **Release = Arc's email or the catalog flip (`vcc --json datasets list`), whichever comes first.** Arc's
10-15 entrants email says only "next week" with no hour; 10-22 07:00Z is the site widget's guess (midnight Pacific). The
watcher polls the catalog every 5 min from 10-19 00:00Z to 10-27 00:00Z (both moved 10-16, coordinator's YES); an AUTO_R0 before
the 10-20 freeze can use only DECIDED switches (plain `command` refuses while one is PENDING) and a later entry may replace it. **Safe deadline: Nov 5 08:00Z.** The
widget says submissions close then; the rules say 23:59Z, so never rely on the tail. Only the LAST final entry counts.

## A. Freeze, by 10-21 20:00Z (each line is a gate; a failed line moves the freeze, not the day)
1. **No PENDING switch:** `uv run python scripts/october_recipe.py show` lists none, and a plain `command` (no
   `--allow-pending`) prints the def_day line. 690a's de_dnmag gets decided by DN671's read. Off is a decision.
2. **Rehearsed on the final tree, judged from outside:** R-NEW (the final recipe on foreign controls + a new panel) DONE,
   postbuild PASS; and **R-KA PASS**: the whole automated chain (watcher → def_day → 692), pinned to 524dr2's decisions
   snapshot, rebuilds 524dr2's board-scored prediction (novelty: "ALREADY BEEN SCORED", entry dVHLws6K = 0.2730889;
   690 RESULTS §9). Any decide or patch after that reruns both (~4 h, in a builder gap).
3. **Known answers on HEAD**: `def_preflight --verify`, `def_identity --verify`, `def_census --verify`,
   `scripts/test_october_recipe.py`, `644/verify644.py` (restore its json after), `410/run_builder_tests.py` (56/0/1).
   Builder on disk == HEAD (`git status --short scripts/` empty), with no uncommitted patch. With DNMAG v2 pinned (693p):
   `693/work/ka693b.log` reads KAU2 PASS and DNV2 PASS (the v2 tree rebuilds the scored DN entry wzaVLFTm…).
4. **Machine**: on AC, lid open, `caffeinate -is` running; disk ≥ 100 GiB free (the day's floor is 60; PIE assets would
   take ~75); no other session's downloads scheduled for 10-22.
5. **Site reachable**: `curl -s -o /dev/null -w '%{http_code}' https://virtualcellchallenge.org/api/leaderboard` = 200.
   Cold Turkey has blocked the site in morning windows (10-05/10-06, ~03:30–07:40Z), which overlaps the 07:00Z release.
   **Ahmet allows virtualcellchallenge.org + github.com for 10-22 → 11-05** (his app; one ask through the coordinator).
6. `vcc whoami` → can_submit true; `vcc --version` == PyPI's latest (`uv tool upgrade vcc-cli && vcc skill install` if newer).
7. Names: claim R0's `-m` with `scripts/next_model_name.py --claim R0`; the reserved final name is in MODEL-NAMES.md.
   `-d` stays EMPTY. **The server runs a content filter on `-m`** (10-19 00:01Z: 1019s1's name containing "shut up" was refused
   at entry creation, before any upload). arm1015 then retries the SAME name until its cap, so a refused R0 would loop for days.
   Screen every claimed name for insults, profanity or violence before arming, and keep a `name_guard.sh` beside each arm1015.
   The coordinator armed `experiments/670-coord-1015/name_guard.sh` (df102a0c) for R0, 1019s2 and 1020s1. On a refusal with no
   pending upload, it claims the next queue name and re-launches arm1015, at most 3 times.
7b. **The real final dataset id is `controls-final`** (`vcc_2026_controls_final.zip`). Arc LISTS it before release:
   since ≥ 10-20 18:58Z, `vcc --json datasets list` shows it with `"available": false, "unavailable_reason":
   "final_week_not_active"`. watch691 fires only when an entry is new AND available (release691.py diff; tested 10-20 on today's list
   with the flag flipped: diff → [controls-final], pick → controls-final). Any check that counts listed ids without `available`
   false-fires. The coordinator's guard_release did exactly that at 18:58:20Z: it wrote NOGO for every A/B/C slot and cancelled W1's
   10-21 s1. Never treat "listed" as "released"; the release is `available: true` (or data/challenge-def/manifest.json).
8. **The last A/B/C slot is 10-21** (review 720 F6). When the site's final_week gate flips (possibly 10-22 00:00Z, before
   the 07:00Z bundle), every submit is server-set `is_final`: no A/B/C waiter may be armed for 10-22.

## B. Arming (watcher + 692 ARMED since 10-14 23:16Z; current pids in PLAN Notes; arm1015 R0 on 10-21)
Before any (re-)arm: the watcher's `baseline:` log line must list only `controls`. If the D/E/F dataset is already listed,
the diff can never fire: do not arm, download by hand (720 R2). A re-launch after 10-27 00:00Z needs `CAP691=<new>`
(and `AUTO_R0=1`, as always) or it stops at once (720 R3). Any edit of watch691.sh / r0ready692.sh = kill + re-arm; both
resume, and a second copy on the same work dir refuses to start (720 R1). Both run from frozen per-launch copies, so an
edit of the source never corrupts a live run. The watcher is armed with the pre-decided overrides (uniform depth shift,
spreads ≤ 3.5; 690 CONTINGENCY-DEPTH.md); when one fires, NEXT.txt and watch691.log say so.
Launch every long-lived job through `scripts/detach.sh` (setsid). On 10-15 a Claude app quit/restart SIGKILLed every plain
`nohup` waiter, caffeinate and snapshot loop, with no reboot and no END lines (cb081c82).
```bash
AUTO_R0=1 AUTO_ALLOW_DEPTH691=1 AUTO_ALLOW_WIDE691=1 nohup bash scripts/detach.sh bash experiments/691-release-watcher/watch691.sh >> experiments/691-release-watcher/work/watch691.nohup 2>&1 &
nohup bash scripts/detach.sh bash experiments/692-def-r0/r0ready692.sh data/challenge-def experiments/692-def-r0/work/r0 r0def 202610270000 \
  >> experiments/692-def-r0/work/r0ready692.nohup 2>&1 &
# PRIOR = the log of our last A/B/C arm that day; if none: PRIOR=experiments/692-def-r0/work/prior-stub.log; echo "$(date -u +%FT%TZ) END" > $PRIOR
nohup bash scripts/detach.sh bash experiments/670-coord-1015/arm1015.sh R0 r0def experiments/692-def-r0/work/r0def.ready.log \
  experiments/692-def-r0/work/r0def.run.log submissions/r0def.vcc "<name from next_model_name.py --claim R0>" "R0 D/E/F" \
  "$PRIOR" 202610190000 202610270000 > experiments/692-def-r0/work/arm-R0.nohup 2>&1 &
```
It polls from 10-19 00:00Z. On release it downloads, unpacks into `data/challenge-def` (A/B/C stays in data/challenge),
runs the builder-free stages, writes `work/summary.json` + `work/NEXT.txt`, and if everything is clean builds **R0**
with def_day.sh into `experiments/692-def-r0/work/r0`. It never preps for real and never submits. After a reboot,
re-launch the same three commands (AGENTS §6.2 step 6). A HOLD file in the watcher's work dir pauses R0.
Also armed (coordinator, 10-16; pids in PLAN): the R0 submit waiter (arm1015 R0, floor 10-19 00:00Z, cap 10-27 00:00Z, prior
`experiments/692-def-r0/work/prior-stub.log`) and the **release guard** `experiments/670-coord-1015/guard_release.sh`.
On a release (data/challenge-def/manifest.json, a `RELEASE DETECTED` line in watch691.log, which watch691 writes on
detection, or an AVAILABLE non-controls id in `vcc --json datasets list`; a listed but unavailable one, such as controls-final before 10-22, is not a release) the guard writes NOGO for every A/B/C slot 1016s1..1021s2, because every
submit after the release is FINAL. After a reboot, re-launch it too (`nohup bash scripts/detach.sh bash …guard_release.sh`).
**Reboot recovery, readiness lane (10-19; macOS has staged a 27.2 RC restart and auto-install is on until Ahmet turns it off).**
In this order, after `uptime` shows a fresh boot. Before every launch, `ps -axo comm=,args= | grep <script>` must show no live copy
(the scripts' single-instance guards refuse a second copy anyway).
1. Builder: `shasum -a 256 scripts/make_submission_pm_cd4.py` must be e99bbd4a… (HEAD). If it differs and
   `experiments/692-def-r0/work/r0/pin-artifacts.json` exists (R0 was killed mid-build), restore it:
   `uv run python scripts/def_pin_panel_artifacts.py --restore --record experiments/692-def-r0/work/r0/pin-artifacts.json`.
   A relaunched watcher prints this same command and stops if you skip it.
2. watch691: the same command as above (`AUTO_R0=1 AUTO_ALLOW_DEPTH691=1 AUTO_ALLOW_WIDE691=1 …`). Before a release it polls.
   After one (work/RELEASE-DETECTED) it RESUMES: it verifies the bundle and skips the download, keeps a FINISHED R0, moves an
   INTERRUPTED R0 dir aside and rebuilds (~4 h), and stops on a STOPPED R0 (then §C2).
3. r0ready692: the same command as above (it re-reads its work dir).
4. wait708: `nohup bash scripts/detach.sh bash experiments/690-def-readiness-1015/wait708.sh > experiments/690-def-readiness-1015/work/wait708.nohup 2>&1 &`.
   If the reboot hit DURING day708 (day708-d708def.log exists without a "DECISIONS WRITTEN" line), wait708 refuses to run twice.
   First move that log and its .nohup aside (`mv … ….aborted-<HHMMSS>`), then relaunch.
5. Coordinator lane: `bash experiments/670-coord-1015/rearm_coord.sh` (safe to run twice: each of arm1015 R0, name_guard R0,
   guard_release, board_snapshot, pisces_watch and caffeinate starts only if no copy runs). A pending R0 upload resumes from
   `~/.config/vcc/state.json` through arm1015's retry. With steps 1–5, §B restores everything if no session is around.
R-KA's rehearsal watcher runs in its own work dir with a stub vcc, so it cannot trip the guard. After the release, stop or HOLD
every A/B/C build chain: AUTO_R0's gate waits for any running builder or day-path chain (and ≥ 26 GiB, ≥ 60 GiB disk).
692 waits for def_day, checks its verdicts (postbuild PASS, DONE, builder restored, UPRN fail-safe 0, not past the cliff),
runs the real prep and the novelty check, writes READY, and commits the round's records (690c's builder.pinned.py
included). arm1015 then submits through the slot path (536 retry + resume, capture, `work/NOGO-R0` vetoes).

## C. The day, by hand (whoever is awake; expected ≈ 07:00 → 10:00Z)
1. `cat experiments/691-release-watcher/work/NEXT.txt`, then `summary.json`. Record into the round's RESULTS.md: contexts,
   **branch** (SAME-PANEL / NEW-PANEL), gene axis identical (if not: **STOP**, 531 F6), per-context depth and spread, and
   each context's identity call, n/8, min_top1 and top 5 (the identity reading rule).
   Also record (a record only: it gates nothing and never delays R0) how many D/E/F targets are RPE1 ∪ Jurkat essential-screen
   targets: `comm -12 <(tail -n +2 data/challenge-def/pert_counts.csv | cut -d, -f1 | sort -u) <(cat data/donors/rpe1/perts.txt
   data/donors/nadig_jurkat/perts.txt | sort -u) | wc -l` (A/B/C: 0; our self-drawn R-NEW panel: 9; the lab predicts ≤ 3, FINDINGS §1
   4b597910). ≤ 3 = A/B/C reads and LABC are the right instruments for D/E/F; a large count = D/E/F was drawn differently.
2. If R0 did not start (watch691.log says why): fix that one cause, then run NEXT.txt's commands (its work dir is the one
   692 waits on, `experiments/692-def-r0/work/r0`). If R0 started and STOPPED, that dir is spent: 692 and arm1015 then
   end NOT READY / NOT SUBMITTING; rerun in a fresh dir and re-arm both with it. SAME-PANEL means 453's operator 66610be7
   is pre-placed first (NEXT.txt carries the line). Guard trips: DECISION-CARD §1–§4.
3. R0 DONE: `experiments/692-def-r0/work/r0def.run.log` says READY (or NOT READY and why). Read the 554 dose row in
   postbuild.json before the waiter fires; touch `experiments/670-coord-1015/work/NOGO-R0` to veto. If R0 was built by hand
   from NEXT.txt, run 692 on that work dir. If watch691.log has an `UPGRADE` line, the vcc CLI changed on the day: re-run
   DEF-RUNBOOK §6.7's prep KA under NOGO-R0 before the first submit (720 R4). An AUTO OVERRIDE line (depth or spread)
   goes into the round's RESULTS.md with the three medians and spreads.
4. After the submit: `data/board-entries/<id>.json` exists (arm1015 captures it; 690d accepts a final entry without rank),
   and the entry is in BOARD-LEDGER with what it reads.
4b. **708 (science; prereg 69bd6105), AFTER R0 is submitted:** re-price the panel-dependent switches on D/E/F's REAL panel
   with held-out public lines as truth (K562-GWPS / Orion HCT116 / HEK293T PLOLO folds). Part A = `dayof438.py` verbatim
   (H1 × CD4 stack weights, 438 §5's rule); Part B = the 705 rider and the 707 swap on the same folds. Inputs from R0's
   def_day work dir `experiments/692-def-r0/work/r0`: `whitener_union_pruned_def.npz` and `cd4_compact/`, so do NOT clean
   that dir before 708 has run. Cost ≈ 1 h 15 m: 3 held-out fold operators (~15 min, 13 GB each, serial), then dayof438
   (~30 min, ≤ 8.8 GB) and dump708/score708 (~30 min). No build runs beside the fits. A move ships only as a LATER final
   entry: its own def_day build with `DEF_DAY_H1_STACK` / `DEF_DAY_CD4_STACK` (690i: H1 x1..x3, CD4 x0.75..x2), the
   usual gates and arm1015. R0 stays the decided recipe; a NULL 708 changes nothing.
   **AUTOMATED (10-19 01:36Z): `experiments/690-def-readiness-1015/wait708.sh` (pid 66950 at arming, log work/wait708.log)
   runs the command below ONCE, after a real release AND R0 is done (SUBMITTED ok, or R0's terminal END with no arm1015 R0 alive).
   Do NOT also run it by hand: wait708 and a hand run would collide, and wait708 refuses if day708-d708def.log already exists.
   If wait708 died (check its pid), re-arm it rather than running the command.** **If W1 (673, H1 ×2.5 on R0's recipe) reads BELOW
   (10-21 s1; SIXA Amendment 2), re-arm wait708 BEFORE the release with `WAIT708_R710=HARM` in front of the same command, so decide708
   itself blocks any H1 > ×2 move.** wait708's log is no one's ready-log, so stopping it is safe. The decision file lands at
   `experiments/708-def-panel-repricing/work/decision708-d708def.txt`, and wait708.log copies it.
   **The day's command, with the three A/B/C board reads decide708 takes (set 10-19 by readiness on the coordinator's 710 relay):**
   `R705=CLEARED R707=WEAK R710=POSITIVE nohup bash scripts/detach.sh bash experiments/708-def-panel-repricing/day708.sh
   data/challenge-def experiments/692-def-r0/work/r0 d708def > experiments/708-def-panel-repricing/work/day708-d708def.nohup 2>&1 &`
   - R705 CLEARED: lstr705 o7ObvpGu05zvmCwsFtiC (10-17 s1), primary −0.000240 > −0.0005.
   - R707 WEAK: sw707 AYbAiHRjuHdil2zCjlsk (10-17 s2), primary +0.0008417, inside +0.0005…+0.0010.
   - R710 POSITIVE: UW4LsxTUkeLumDCBotfU (10-19 s1, 0.2765675), read710 Δoverall +0.0004824 vs ALLSW ≥ +0.0003. decide708 blocks an
     H1 > ×2 move only on R710=HARM; POSITIVE reads "does not block". Tested 10-19 on the r708abcdef rehearsal, HARM vs POSITIVE
     (the documented values were HARM|NOT_HARM|NONE).
   - The env passes through detach.sh to day708 and its decide708 child. day708's arming line logs R705/R707 only, so the R710
     value shows only in decision708-d708def.txt, and only if Part A's winner raises H1.
5. First `vcc status <id>` after scoring. The prior is DARK (2025 hid final scores until the winners), so DEF-RUNBOOK §6 is
   the default branch; members returned → FINAL-PHASE-FEEDBACK-PLAN (maximum-information use, documented).
6. **The last entry** (submitted AND scored before Nov 5 08:00Z; start the upload by Nov 4 ~20:00Z). A new package goes
   through arm1015.sh with the reserved final name. Re-uploading an already-scored final package goes through
   `experiments/670-coord-1015/arm1015f.sh` with `FINAL_REUPLOAD=1` (690f: reserved name, ledgered package whose capture is
   a final entry, ledger key `<sha>#reupload-<id>`). The server's handling of identical bytes is untested: run ONE early
   harmless test (`FINAL_REUPLOAD=1 FINAL_REUPLOAD_TEST=1`, R0 re-uploaded on a day nothing else is ready), never first on
   Nov 4.
6a. **A post-R0 final entry REPLACES R0** (only the last entry counts), so a buggy one costs far more than 708's
   +0.001–0.002. It is submitted only if EVERY line holds; otherwise nothing is submitted and R0 stays final:
   - (a) built by the same def_day path with ONLY the recipe changes 708 decided (`october_recipe.py decide …`, evidence
     = decision708-<tag>.txt); its final meta's config diff vs R0's (`experiments/692-def-r0/work/r0/r0def_final.h5ad.meta.json`,
     sameline.path ignored) equals EXACTLY the keys those switches act on (run_arm_generic's EXPECT rule; e.g. H1/CD4 dose →
     bulk_load(+sha), bulk_ctx_scale **and bulk_sharpen_pt** (coordinator, 10-20: W1 673 false-FAILed without it; a dose changes the stack
     dump that tercilec's key ratio is normalised on, so `edges` / `key_ratio_c` / `key_ratio_maxdev` move while every behavioural field
     — mode, q, counts, carried key + sha, moved_vs_tercile, n_carried — must stay equal: check with `experiments/670-coord-1015/regate673.py`'s
     rule); lodo_strength → + pert_scale(+sha)). Use the CONTENT-aware diff, never a plain one:
     `uv run python scripts/def_metadiff.py <new>_final.h5ad.meta.json <new work dir> experiments/692-def-r0/work/r0/r0def_final.h5ad.meta.json
     experiments/692-def-r0/work/r0 --expect <keys>`. Every work dir writes its own path into the dump's `env` and the same-line
     `provenance`, so `bulk_load_sha256` / `sameline.sha256` differ with identical content. A plain diff then FAILS falsely
     (R-ABC vs 524dr2, 10-16) and would keep a worse R0. Keep R0's work dir (its npz are the reference) until the final closes;
   - (b) R0's gates pass again: postbuild all PASS, zero-gene 0, CLIFF OK, every rho ≤ 1 and in band, novelty rc 0;
   - (c) where it applies (only DE riders changed), bulk identity vs R0 (`bulkid658.py`: per-(context, target, gene) sums
     identical); a weights move changes the bulk by design, so (b)'s CLIFF / rho band is then its only bulk guard;
   - (e) the same tree as R0: the meta's `hashes.implementation` and `de_sign_emitter_sha256` equal R0's, and def_day.sh is
     the sha in run690.log's last "R-KA tree:" line (decide708 checks that line; 690i APPLIED for any dose beyond CD4 x1.5);
   - (d) one attempt = ONE SUBMISSION (clarified 10-19 with the coordinator). WHICH candidate goes in (X0 = R0, X1 = R0 + 708's
     decided changes, X2 = X1 + 811) is fixed in advance by `experiments/670-coord-1015/SIXA-DECISION.md` incl. Amendments 1–3 (HEAD 88200fd1+; 47249e76 = the base rule): build X1,
     then X2 if eligible; ship the most complete one that passed every line here; if none passed, nothing is submitted and R0 stays
     final. The choice follows that file's pre-written order, never a score. **No fix-and-retry:** each candidate is built at most
     once.
     - A build interrupted by a reboot, OOM or crash before its gates may be rebuilt with unchanged inputs.
     - A gate failure shown to be a judge false-FAIL (as 697 and 698) is re-gated on the finished build, not rebuilt.
     - A real gate failure retires that candidate.
     No new replacement attempt after 11-04 12:00Z (the upload must start by Nov 4 ~20:00Z for the 08:00Z widget close).
   - (f) the name: claim it from the queue (`next_model_name.py --claim <tag>`) and screen it for insults, profanity or
     violence before arming (§A 7). Arm `name_guard.sh` beside its arm1015, with spec = the arm1015 args minus the name
     (one per line: SL TAG READY-LOG RUN-LOG PKG NOTE PRIOR-LOG FLOOR CAP): `nohup bash scripts/detach.sh bash
     experiments/670-coord-1015/name_guard.sh <spec> <guard cap> > experiments/670-coord-1015/work/name_guard-<SL>.nohup 2>&1 &`.
     Its PRIOR is R0's slot, never `submit-R0.log` itself: arm1015 `kill -0`s the first pid its prior log names, and pids wrap here
     (10-20), so R0's 3409 will be someone else's by then. Arm `nohup bash scripts/detach.sh bash experiments/670-coord-1015/prior_relay.sh
     R0 <stub> <cap> > experiments/670-coord-1015/work/prior_relay-R0.nohup 2>&1 &` and pass `<stub>` as arm1015's prior-log. It writes the
     stub (R0's entry id + END, no pid) only once submit-R0.log is terminal and no `bash arm1015.sh R0` runs. R0 is that log's only
     writer, so an old END cannot trigger it.
   - Amendments in force (SIXA-DECISION): a swap-OFF X2 needs C0's (699) board read to be NOT HARM; a W1 (673) BELOW blocks H1 > ×2 via
     `WAIT708_R710=HARM` (4b); a W2 (674) BELOW means X1 may not carry CD4 > ×1, applied by the coordinator at decision time.
     Amendment 4/4b (coordinator, 10-20): on a swap-ON day with X2 eligible, X2 keeps the swap only if C1:C0 on DEFo has ci95[0] > 0; else
     X2 = C0 (R0 + 811) if C0 − B0 passes its own bars, else C1 (one X2 build). The W1 BELOW route is automated by
     `experiments/670-coord-1015/w1route.sh` (log work/w1route.log): do not also re-arm wait708 by hand for W1.
     Amendment 5: X2 is only ever the priced construction (C1 = R0 + swap + 811, or C0 = R0 + 811), never on a 708 weights move or rider veto;
     ship order X2, else X1, else R0. The whole rule, mechanically: `experiments/670-coord-1015/decide6a.py` (16/16 tests).

## D. Timings (measured; plan with the right column)
| step | typical | worst seen |
|---|---|---|
| poll → download → unpack → light stages | ~1 min (stub test; real download ~630 MB) | Cold Turkey block: until it lifts |
| def_day R0, the October recipe as shipped (UPRN + DNMAG v2 ON; setup, dump, 3 levelled builds, audits) | **4 h 10 m** (R-ABC 10-16: dump+lvl0 81 min, lvl1 84, lvl2 84; 524dr2 without UPRN/DNMAG: 1 h 41 m, ~34 min per level) | +~80 m if a 4th level or a relaxed CD4 rerun is needed |
| where the extra ~50 min per level build goes | DNMAG v2 ≈ 15 min (693); the rest is UPRN's rank-neutral UP pass (630), the next speed target | — |
| memory on the day | every R0 gate: free+inactive ≥ 20 GiB AND memory_pressure free ≥ 40 % (AGENTS §9) = def_day's stage floor `DEF_DAY_MIN_FREE_GIB` 20; no swap file; quit heavy apps if fi < 20 | — |
| real `vcc prep` | 15 s, 20.5 GiB RSS (R-NEW 10-16) | — |
| upload 2.8 GB | 11–19 min | 7.5 h (09-29) |
| scoring queue | 33–62 min | 12.5 h (10-06: cancel recipe in memory) |
So R0 is READY ≈ release + 4 h 15 m (A/B/C-like re-level) to ≈ release + 5 h 45 m (R-NEW 10-16: a new panel whose re-level
needed 3 levelled builds after one overshoot; 91 + 85 + 82 + 81 min), and scored ≈ 1–2 h later. A second build the same day costs
the same again. The real `vcc prep` takes 15 s but 20.5 GiB RSS.

## E. Late recipe changes (Oct 15 → 21): what each kind touches
| kind | example | touch | lead time |
|---|---|---|---|
| a switch def_day already carries | q, sharpen_pt, de_mag, mde_sign, de_uprn, de_dnmag (690a), rc_callset_cpm, cd4_stack | `october_recipe.py decide <sw> <v> --evidence …` | minutes, plus R-NEW + R-KA reruns (~4 h) |
| a builder `VCC_*` knob with no pass-through | 16 today, e.g. VCC_CARRIER_GAMMA, VCC_DONOR_WEIGHTS, VCC_WHITEN_LAM | a 690a-style patch: validate, export to the level builds, final-meta check, final.env, postbuild ALLOW, recipe switch, tests | ~2 h + rehearsal |
| a new day-of input computed from D/E/F's controls | the PIE blend ([bf003b]'s `700/dumphook700.sh`) | already plumbed: 690g's `DEF_DAY_DUMP_HOOK` (validated dump swap, optional 672 per-context npz); decide `dump_hook <path>`, then rehearse WITH the hook (offline, HF_HUB_OFFLINE=1, ~75 GB assets on disk) | the decide + a ~2.5 h rehearsal |
| an emission-law change | — | builder + the all-unset KA (== ec3a5ae6) + rehearsal | ≥ 1 day |
def_day.sh unsets every `VCC_*` at start, so nothing reaches a build except through a pass-through.

**Pending late switches (10-15 21:35Z):** 705 `lodo_strength` 0.25 (read 10-17 s1; CLEARED unless HARM ≤ −0.0005) and 707's
day-of dump hook `experiments/707-strong-tier-direction/dumphook707.sh` (read 10-17 s2 / 10-18 s1; ships only on TRANSFERS ≥ +0.0010).
For each that clears: `october_recipe.py decide …`, then ONE rehearsal on the NEW-panel bundle, both riders together if both clear:
`RIDER_LSTR=1 RIDER_HOOK=1 nohup bash scripts/detach.sh bash experiments/690-def-readiness-1015/rider690.sh >> experiments/690-def-readiness-1015/work/rider690.nohup 2>&1 &`
(drop the flag of a switch that did not clear). It judges against R-NEW's own final meta (R-NEW PASSED 10-16) and waits for
R-ABC2's driver to exit (RIDER_AFTER). If 711 adopts K 3,400, rider690-K is the same command with RIDER_K=1 (DEF_DAY_RC_CALLSET_CPM=3000)
plus the adopted θ / hook flags and `RIDER_AFTER="rabc2_695.sh rabc3_696.sh"` (coordinator's order: R-ABC2 → rider690 → R-ABC3 → rider690-K).
Its gates are 522 §8 (a) KA, the def_day rider/hook lines, and a content-aware config diff only in bulk_load(+sha)/bulk_ctx_scale
(+ bulk_sharpen_pt with LSTR: θ renames the dump `_s0.25` and changes the tercilec key sha; the level builds keep their pert_scale). DRY 10-15: 522 KA 3/3 PASS, env composes. Not rehearsed by 10-20 12:00Z (its cap) = the switch ships OFF.

## F. What the evidence behind the shipped recipe covers (690 RESULTS §10)
Board: A/B/C only. **R-ABC 0.2743318 (10-16) is the exact October recipe** (same resolved env, tree and def_day path as
R0): +0.0012429 vs 524dr2, UPRN + DNMAG composing at 95 % of their separate reads; pds/mse/fid/jac identical to 524dr2 (694). Sealed and fresh-panel: the bulk F = M6+P2g (R6 C10b, R7a)
and the CLEAN4 hook (R7b). Unsealed: UPRN, DNMAG and PIE (board read + local prices only). There is no honest offline
estimate of the D/E/F score: no leak-free, powered truth for a new line exists on disk. Say so in every summary.

## G. After Nov 5 (the final is closed)
- **Re-enable sleep** (Ahmet, sudo; set 10-15 ~22:00Z so a closed lid no longer sleeps the Mac): `sudo pmset -a disablesleep 0`,
  then `pmset -g | grep SleepDisabled` must show 0. Until then the Mac stays on AC: SleepDisabled does not stop a battery
  running out.
- Stop the long-lived jobs (caffeinate, board_snapshot, pisces_watch, any waiter) and remove an installed 691 LaunchAgent if one was
  ever installed (`launchctl bootout gui/$(id -u)/com.ahmet.vcc.watch691`).
