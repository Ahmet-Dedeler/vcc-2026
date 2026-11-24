# 690 — D/E/F day-of readiness on the current tree (2026-10-15, session [17adab]; in progress)

Class: day-of procedure (scorer mechanics / delivery). No recipe value changes; nothing here moves a score.
Role from the 10-15 five-session split ([400b6e] coordinator, [020377] news, [bf003b] science, [17adab] readiness).

## 1. Known answers on HEAD (2026-10-14 22:2xZ; builder fe4678c3 == HEAD, before 631's patch)
- `def_preflight --verify` PASS, `def_identity --verify` PASS (A/B/C = Jurkat / HeLa / CAL33, 8/8), `def_census --verify`
  GATE PASS (CD4 union 291, complete 276, zero-carrier 0). `work/verify1.log`.
- `scripts/test_october_recipe.py` 6/6; `644/verify644.py` all_ok (33 rows; its rewritten json was restored, the run's
  copy is `work/verify644-1015.json`); `410-pertscale-revisit/run_builder_tests.py` 56 passed / 0 failed / 1 skipped
  (the canonical runner; `scripts/test_make_submission_pm_cd4.py` run directly fails one stale test, not a gate).
- Environment unchanged since 09-26: numpy 2.5.2, scipy 1.18.0, h5py 3.16.0, anndata 0.13.2.
- `vcc` 0.2.2 is still the newest on PyPI (0.2.2 published 09-23). `cell-eval2` 0.18.0 (10-01) says "Scoring is
  unchanged: rule_version stays 3, rule digest bit-identical to 0.16.0's"; we run 0.16.0, so no upgrade is needed.

## 2. Day-of hazards found
1. **Same panel.** Arc's 10-01 post says the same sgRNAs were used in all six lines, so D/E/F may carry A/B/C's 300
   targets. `build_union_whitener.py` reads only `gene_names.csv` and `pert_counts.csv` from the bundle, so on the same
   panel the day's fit reproduces A/B/C's operator `e2f65ecf` byte for byte. `def_pipeline.py` then aborts by design
   ("the re-fit operator is byte-identical to the A/B/C one: the panel did not change"), and a plain copy then trips the
   plolo registry at the key stage (453 attempts 1 and 2). 524dr / 524dr2 ran this path by pre-placing 453's re-saved
   operator (`66610be7`, array-equal to e2f65ecf, its own identity) in the work dir. It works, but it was not a day-of
   branch in the runbook. Rehearsal R-SAME runs it on non-A/B/C controls.
2. **PENDING switches.** `october_recipe.py command` (without `--allow-pending`) refuses: rho_scale, lodo_strength, lodo_hi,
   callset_block are PENDING. Each must be DECIDED (off is a decision) before Oct 22, or the day uses `--allow-pending`,
   which the decision card forbids.
3. **Late riders unreachable.** def_day.sh passes DEF_DAY_DE_UPRN (630dd) and DEF_DAY_MDE_SIGN_KEYGATE (612dd) through,
   but the recipe has no switch for them; 631's VCC_DE_DNMAG has no def_day pass-through. Fixed by 690a (below).
4. **def_day.sh refuses an uncommitted builder** (builder on disk must equal `git show HEAD:`). Any builder patch must be
   committed before a day-path build.

## 3. 690a (UNAPPLIED patch, the coordinator applies)
`patch690a.py`: DEF_DAY_DE_DNMAG pass-through in def_day.sh (validation with the builder's grammar, export with the DE
riders, final-meta check of variant + dn_mag module 0a8cc560, final.env, postbuild ALLOW), recipe switches de_uprn /
mde_sign_keygate / de_dnmag (PENDING, dependency + grammar checks, installed-knob checks), and the recipe test fixture.
`test690a.py` 27/27 PASS (`work/test690a.log`), including the KA "resolved env and args unchanged while the switches are
PENDING". `--check` passes on the live tree.

## 4. Mock D/E/F bundles (`build690.py`, `work/bundles/build690.json`)
Real foreign control cells in the exact bundle format (595's readers): D = HCT116 (Orion NTC), E = H1 (Arc 2025 Flex
controls), F = HEK293T (Orion NTC); 18,400 cells, 46 ntc_id. `new/` carries ef-bundle-r's 300 targets, `same/` A/B/C's;
the context files are APFS clones (identical bytes). Build 1 min, peak RSS 11.8 GB.

| ctx | line | median library | q75/q25 |
|---|---|---:|---:|
| D | HCT116 | 16,804 (×0.839 of A/B/C) | 1.771 |
| E | H1 | 20,186 | 1.575 |
| F | HEK293T | 16,114 (×0.804) | 1.793 |

## 5. Light stages (driver `--stages preflight,identity,cd4-extract,census`; no builder pin)
- Preflight: the 641 depth row FAILs D and F, as designed (real D/E/F are released at ~20k, Arc). Rehearsals export
  `DEF_DAY_ALLOW_DEPTH=1`; recorded here.
- Identity: D = HCT116 (registry, same-line ON), E = TERA2 (H1 absent from CCLE, split; abstains), F = HAP1 (HEK293T's
  false call; blocked by 590's 0.65 floor). Same on both bundles.
- CD4 extract: union 300 (new) / 291 (same); census GATE PASS on both (no zero-carrier targets; panel_support predicted
  279 / 290 on new).

## 6. Rehearsals (armed: `run690.sh`, waits for dn671 and the coordinator's HOLD; R-NEW ~05:00Z, R-SAME ~07:15Z)
- 10-15 04:45Z: run690 re-armed (pid 48612). R-NEW now also waits for ka693b (KAu2 + DNv2) to EXIT. ka693b's gate sees only
  the python builder / vcc prep, not def_day.sh's setup stage, so after the HOLD lifted the two could have started in one poll.
  Order: KAu2 → DNv2 → R-NEW → R-KA → 694 R-ABC (the October recipe on A/B/C for 10-16 s1).
- R-NEW: def_day.sh on `work/bundles/new` with the resolved recipe (`command --allow-pending`) + DEF_DAY_ALLOW_DEPTH=1,
  then one real `vcc prep` (timed, deleted).
- R-SAME (v2): the Oct 22 automation itself. 691's watcher with AUTO_R0=1 on a stub release of build690b's
  depth-matched same-panel bundle (Orion HCT116 / HEK293T importance-subsampled to context A's library distribution;
  H1 as is). No guard override; the only rehearsal knob is R0_ALLOW_PENDING691=1, logged.

- **R-KA (10-16 07:27→09:12Z): PASS in substance; corrected verdict 690k on record.** The whole Oct 22 chain (691 watcher AUTO_R0
  on a stub release of 180's def-bundle, SAME-PANEL branch, then def_day, then 692 prep and novelty) with 524dr2's decisions
  snapshot rebuilt 524dr2's def_day output BYTE FOR BYTE (output.sha256 c5b6e359…). That is exactly the file the board-scored
  package was relabelled from (524dr2 meta relabel453.source_output_sha256; ledger entry dVHLws6K, h5ad 0c394160, 0.2730889), and
  def_metadiff vs 524dr2's final shows no real diff. Tree: def_day.sh 11ff102e (690i APPLIED, the AGENTS memory floor), builder e99bbd4a.
  run690 logged FAIL at 09:12:07Z because my check asked novelty whether the D/E/F package was the A/B/C relabel. Novelty hashes
  `context`, so that could never match. This was a design error in my check, not the chain. `rka690k.py` re-verified the six-link sha
  chain and appended the corrected PASS + tree line (decide708 reads PASS for every dose). run690's check now judges the def_day
  output sha + metadiff. Coverage: unset riders only (524dr2 leaves UPRN, DNMAG and 690i's doses unset); those keep their own KAs
  (R-ABC's board read, DNV2, test690i).

### 6b. 10-16: R-NEW found a real day-of failure (690j)
- R-NEW (new 300-target panel, foreign controls: the D/E/F shape) died at 07:27:51Z in stage 2b, `def_lodo_key.py` →
  lib438.residuals → sens415.cd4_frozen_and_partial → 062 prepare_cd4_mean(max_rss_gib=8.0): "MemoryError: cd4_rest centering:
  peak RSS 8.530 GiB reached 8.000 GiB". The guard tests the PROCESS's peak RSS (ru_maxrss), so it depends on the panel and on
  everything loaded before, not on free memory. 460's A/B/C KA already peaked at 8.0 GB. **R0 on a new panel would likely have
  died here on the day.** This is the find R-NEW was armed for.
- Fix (5e1aa480): def_lodo_key lifts 062's enforce_memory_limit to `DEF_DAY_LODO_MAX_RSS_GIB` (default 16) for its own process,
  covering sens415's 8.0 and lib438's 9.0 caps. Guard-only. sens415.py / 062 run.py are in R7a's pinned closure and are untouched.
- Proof pending: chain690j.sh runs def_lodo_key --ka / --ka-pct on fresh lib438 tags (must reproduce 8c82bef3 / dc374281), then
  reruns R-NEW alone (`ONLY690=rnew`), at the coordinator's GO-690j. Also the same morning: the stage-memory floor 23 → 20 (fc59a367)
  and the outer gates moved to the AGENTS §9 rule. The builder's own cap (def_pipeline --max-rss-gib 24 vs R-ABC's 20.2 GiB peak)
  is next in line; the rerun measures it on a new panel.
- **690j proof PASS (10-16 13:10Z):** def_lodo_key with the lifted guard, on fresh lib438 tags, reproduces A/B/C's 8c82bef3 (--ka)
  and dc374281 (--ka-pct: max|dw| 3.2e-14, max|dkey| 1.0e-15), each at 7.47 GiB peak; the lift was applied on 2 module objects.
  **The R-NEW rerun** (ONLY690=rnew, started 13:10:51Z) is past the step that killed it: the key fitted on the new panel (carried
  300/300, LODO cap 7.48) and the stack dump was written at 13:17:42Z.
- **R-NEW PASS (10-16 18:51Z): the October recipe as shipped, on a NEW 300-target panel with foreign controls, ran through def_day.**
  CONVERGED after 3 levelled builds at D 1.1425 / E 1.2837 / F 1.1880, final rho D 0.9755 / E 0.9807 / F 0.9717 (all ≤ 1, in band);
  postbuild PASS; 690h DNMAG LOW_CONTEXTS none; builder restored == HEAD; the real `vcc prep` took 15 s at 20.5 GiB RSS. The re-level was steeper
  than on A/B/C: lvl0 rho 0.61 / 0.61 / 0.54, local exponent ~3.9 to 4.5 vs ~2. lvl1 overshot D (1.54) and F (1.72) and was deleted as
  past the cliff. The secant then bracketed and converged. Wall time 5 h 40 m: setup + operator fit + LODO + dump + lvl0 91 min,
  lvl1 85, lvl2 82, lvl3 81. Builder peaks 17.4 / 19.4 / 12.4 / 16.7 GiB (cap 24).
- **rider690 PASS (10-17 11:06Z): θ 0.25 rehearsed on the NEW-panel bundle.** The 522 §8 (a) KA passed. def_day logged the rider; it
  CONVERGED after 3 levelled builds (D 1.1125 / E 1.2853 / F 1.1567); lvl1 overshot D 1.56 / F 1.82 and was deleted, as on R-NEW;
  327 min. Postbuild PASS, zero-gene 0. The content-aware diff vs R-NEW is exactly bulk_load (+sha), bulk_sharpen_pt and
  bulk_ctx_scale, θ's signature, the same as on A/B/C (695). θ now has its set-value rehearsal on both panels (A/B/C: R-ABC2).
  Then a deadlock was caught: 697's wait list named rider690.sh, which also matched the WAITING rider690-H (which waits for 697).
  Fixed (65ad66e1), and R-ABC2H started 11:22Z.
- **rider690-H PASS (10-17 22:00Z): θ 0.25 + the 707 dump hook rehearsed together on the NEW-panel bundle.** Started 16:23:45Z (pid
  85922; the gate saw builder e99bbd4a == HEAD); the 522 §8 (a) KA passed again. The hook (dumphook707.sh df21eb69) ran at 16:30:12Z
  on the θ dump `dump_h1x2_cd4x1_s0.25.npz` (60673a2f) and replaced it with 96e2d59a: 100 of 300 rows moved (top z_src tertile;
  donors k562_gwps / orion_hct116 / orion_hek293t cover 300/300), row norms preserved to 6e-15, change energy 2.1 % of the bulk,
  moved-row cosine to input 0.974–0.994. The final meta records the HOOK's dump sha, not the stack dump's. CONVERGED after 3 levelled
  builds at D 1.0582 / E 1.2490 / F 1.0809, final rho D 0.9788 / E 0.9710 / F 0.9774 (target 0.9793 ± 1 %). lvl0 rho 0.77 / 0.61 /
  0.69; lvl1 overshot D 1.34 / F 1.63 and was deleted as past the cliff, as in every new-panel run; lvl2 0.97 / 0.97 / 0.95. 336 min
  (θ alone: 327). Builder peaks 15.6 / 15.0 / 16.2 / 20.1 GiB (cap 24). 690h DNMAG LOW_CONTEXTS none; postbuild PASS; zero-gene 0;
  builder after restore e99bbd4a == HEAD; def_day.sh still 11ff102e (frozen).
  Two diffs, both content-aware (def_metadiff, rtol 1e-12):
  1. vs R-NEW (rider690's own check): bulk_ctx_scale, bulk_load, bulk_load_sha256, bulk_sharpen_pt. That is θ's signature, plus the
     hook's sha.
  2. vs the θ-only rider (`dr_rider.aside-162347`), which isolates the hook: exactly `bulk_load_sha256` (REAL: array bulk differs)
     and `bulk_ctx_scale`, the hook's signature, as on A/B/C (697). `work/dr_rider/metadiff_hook_vs_theta.json`.
  The hook's re-level factors sit 0.04–0.08 below θ alone's, so its bulk needs less lift to reach the NTC level. θ and the hook
  now each have a set-value rehearsal on both panels. The hook stays off in docs/october-decisions.json (707 WEAK) and enters only by
  a 708-style decision. K 3,400 gets its new-panel rider (rider690-K) only if read711 adopts it. The 12 GB final h5ad
  (`work/dr_rider/dr690r_final.h5ad`, no package) stays until 10-19 for any science read; it can be deleted after that.

### 6c. R-KA2 (preregistered 2026-10-18 03:25Z, before the run; requested by the coordinator [2d6f46])
**Question.** R-KA proved the Oct 22 chain on a frozen snapshot of an older recipe (524dr2's). Nothing has shown that R0 reproduces
the board-read R-ABC2 byte for byte when it resolves its recipe from the LIVE `docs/october-decisions.json` (691 → october_recipe →
def_day → 692). R-KA2 is that last link, between "what we read" and "what ships".
**Run.** `rka2_690.sh`, the R-KA section of run690.sh with three changes:
- `OCTOBER_DECISIONS` is UNSET (explicitly `env -u`), so october_recipe reads `docs/october-decisions.json`. At registration that file
  is sha256 04be2af3…, last commit 656dc69a; the run logs its sha at start and again after def_day.
- No `R0_ALLOW_PENDING691`: the plain `command`, as on the day (the live file has pending []).
- Fresh dirs `rka2_watch`, `def_rka2`, `dl_rka2`, `dr_ka2`, `r0ready_ka2` (R-KA's are kept). OUT692 / PKGDIR692 are in 690's work
  dir, COMMIT692=0, and the stub vcc means no notify. The real watch691 (12330), r0ready692 (96817) and arm1015 R0 (97739) never
  see any of it.
Same stub release as R-KA: `work/zips/final_defbundle.zip` (sha 6d7de771…, 180's def-bundle), the SAME-PANEL branch with 453's
operator 66610be7. Same gates as run690: no HOLD, builder == HEAD, no builder, no chain, AGENTS §9 memory, disk ≥ 35 GiB. Watch691
also requires ≥ 60 GiB. Latest start 18:30Z, so it ends before 23:00Z (R-ABC2's def_day took 4 h 05 m); otherwise DEFERRED until
after R0. Never submits; the trial package `dr_ka2/r0def.vcc` is deleted after 692.
**Tree at registration** (== R-ABC2's recorded tree, work/env.json):
- def_day.sh 11ff102e
- builder e99bbd4a (== HEAD)
- emitter rc_emit_k 0f56209d
- helpers 5354e449
The only scripts/ change since R-ABC2's build is def_metadiff's float tolerance (623d4373), and no build reads it.
All 24 switches have the same status/value as R-ABC2's at-resolve snapshot; the resolved env is identical (8 exports).
**PASS iff all of:**
1. the stub watcher exits rc 0 after AUTO_R0 ran def_day on the SAME-PANEL branch with def_day rc 0;
2. `dr_ka2/r0def_final.h5ad.meta.json` output.sha256 == **c0a636a00437120705a8d5d8ba035bebfc84a44e1c31fcdca812fa7f98aae2fb**
   (R-ABC2's def_day final, `experiments/695-rabc2/work/run/rabc2_final.h5ad.meta.json`; board 0.2739849, nnHQfcZHsT3as9EafA3i);
3. `def_metadiff` vs that final and run dir with `--expect ''`: PASS;
4. 692 writes READY (real prep rc 0, novelty rc 0) into `r0ready_ka2`;
5. the decisions file sha is unchanged across the run.
**Otherwise FAIL**, as a bug to find before Oct 22:
- a sha mismatch with an empty metadiff = nondeterminism or path dependence in the build;
- a non-empty metadiff = the live file or the chain resolves differently from what R-ABC2 built.
Judged on the def_day output (690k's like-stage rule), never on a relabelled package.

**R-KA2 PASS (2026-10-18 07:32:26Z): all five conditions hold.** Run 03:26:39 → 07:32:27Z (driver pid 18108, gate open at once).
1. The stub watcher detected the release at 03:26:45Z (zip 6d7de771), took the SAME-PANEL branch and pre-placed 453's operator.
   AUTO_R0 started def_day at 03:27:10Z from the plain `command` on the live file (the 8 exports, pending []); def_day DONE at
   07:31:45Z (4 h 05 m, as R-ABC2's 4 h 05 m); watcher rc 0 in 245 min.
2. def_day output.sha256 **c0a636a00437120705a8d5d8ba035bebfc84a44e1c31fcdca812fa7f98aae2fb == R-ABC2's**, byte for byte. The
   re-level retraced R-ABC2 exactly: lvl0 rho 0.8673 / 0.7994 / 0.7766; lvl1 at 1.0626 / 1.1069 / 1.1230 overshot F (1.0038) and
   was deleted; it converged at lvl2.
3. `def_metadiff --expect ''` PASS: diff []; bulk_load_sha256 and sameline.sha256 are provenance-only (work-dir strings).
   `work/rka2_metadiff.json`.
4. 692 READY r0def 7b014014… (trial package, deleted): real prep rc 0 in 15 s (peak 20.5 GB), package 2,781,921,280 bytes,
   novelty rc 0.
5. docs/october-decisions.json 04be2af3 → 04be2af3, unchanged.
Builder after restore e99bbd4a == HEAD. Isolation held:
- the real 692's run log was last written 10-16 03:50Z, and there is no real r0def.ready.log;
- the real watch691 (12330) shows no RELEASE line;
- the final h5ad was deleted after its sha was logged.
**What this closes.** The release-day path, resolving its recipe from the live decisions file, produces exactly the prediction the
board scored at 0.2739849 (R-ABC2, nnHQfcZHsT3as9EafA3i). From here, "what we read" and "what ships" differ only by the D/E/F
data. The proof covers this tree (def_day 11ff102e, builder e99bbd4a, emitter 0f56209d, helpers 5354e449, at 82c6c37f) and this
decisions file (04be2af3). Any later edit to either needs its own KA.

**R-KA3 PASS (2026-10-20 20:42:10Z): the same proof under macOS 27.2.** Ahmet installed macOS 27.2 (26B118) and the Mac rebooted at
~15:31Z. numpy 2.5.2 and scipy 1.18.0 link Apple Accelerate, which the OS update replaced, so the coordinator reran `rka2_690.sh`
unchanged (driver 7510, RKA2_LATEST 202610202000; prereg in PLAN 4123337c; reading rule: PASS → nothing; drift → record; structural
diff → investigate).
- def_day 15:40:49 → 20:41Z, the same tree (def_day 11ff102e, helpers 5354e449, builder e99bbd4a, emitter 0f56209d, at 4123337c).
- The re-level retraced R-KA2 / R-ABC2 digit for digit: lvl0 0.8673 / 0.7994 / 0.7766; lvl1 1.0626 / 1.1069 / 1.1230 → 0.9730 / 0.9689 /
  1.0038 (F over, deleted); converged lvl2.
- **def_day output.sha256 c0a636a0… == R-ABC2's**; `def_metadiff --expect ''` PASS; decisions 04be2af3 unchanged.
- 692 READY with the SAME trial package sha as R-KA2 (7b014014…): the prep and the packaging are byte-stable too.
- So the new Accelerate changes nothing in the build at the byte level. The frozen recipe and its proof stand on the OS that will
  run R0.
- Side fix found on the way (0581ff00): rka2_690.sh now also moves a stale `def_rka2.unzipped` aside. The first R-KA3 attempt stopped
  on the 10-18 one, and the coordinator moved it by hand.

### 6d. Judge fixes on 10-18 (not regressions)
- `scripts/def_metadiff.py` gained an absolute float floor, `--float-atol` (default 1e-12), beside rtol 1e-12 (f9b6e09d). 698's C1
  gate falsely failed on tercilec's round-off diagnostic key_ratio_maxdev (8.66e-15 vs 1.02e-14). Exact is now
  `--float-rtol 0 --float-atol 0`. The coordinator reviewed and accepted it. It also protects the §C 6a replacement protocol from
  the same false FAIL. test_def_metadiff 10/10; R-KA2's metadiff is still diff [].
- `experiments/694-rabc/test694.py` G6 was FLAKY by PYTHONHASHSEED (gate694's `_norm` replaced path tokens in set order), so a
  10/11 there was not a regression. Fixed (1cc0bada): longest token first; 11/11 for seeds 0–7. See 694 RESULTS. It could only
  cause false FAILs, so no verdict changes.
- (Repair: c2b702e0's §6c edit had dropped this file's `## 7. 690c` heading; restored here.)

## 7. 690c (UNAPPLIED; the coordinator applies it after 690a)
The final meta records the builder AFTER the day's pins (524dr2: d5cec53d), but def_day's trap restores HEAD and the pin
records keep only the pre-pin copy. 690c saves `<work>/builder.pinned.py` inside the trap before the restore, so the
round commits it and HEAD stays clean (the runbook's KEEP_PINS-and-commit route would leave HEAD pinned and break every
later def_day). test690c 7/7.

## 8. The one-pager
`OCT22-CHECKLIST.md`: freeze gates (10-21 20:00Z), arming the watcher, the day's hand steps, measured timings, and
the late-change surface (which kind of recipe change touches what, and its lead time).

## 9. The deeper pass (2026-10-14 23:30Z →): outside judges instead of our own PASS lines
Ahmet asked for a deeper, not bigger, version: what a top professional would add, from the rawest material, with judges
other than us. Every readiness check above was the pipeline grading itself (postbuild PASS, prep dry runs, stub tests).

**R-KA, the board-validated known answer of the whole Oct 22 chain (replaces R-SAME in run690 v3).**
- The chain: 691's watcher with AUTO_R0 on a stub release of 180's def-bundle (the SAME-PANEL branch, 453's operator
  pre-placed) → def_day → 692's real prep and novelty check.
- Recipe resolution is pinned to 524dr2's own decisions snapshot (`OCTOBER_DECISIONS=658/dr2/october-decisions.at-resolve.json`;
  same env as recipe524.json), on today's builder (631 + 690a/c/g, all inert unset).
- PASS iff the novelty checker recognises the prediction as 524dr2's board-scored one. Entry dVHLws6KdMxudq0Yn5ap, board
  0.2730889; ledger h5ad `0c394160…`, content `50ec8ff1…`. Cell names are `cell_{i}`, so the run name cannot spoil the
  content digest.
- The judge is the leaderboard plus the ledger, not us. Result: pending its window (after R-NEW).

**A day-of bug the design exposed (fixed, 305f494a).** r0ready692 ran the novelty check before placing
`<stem>.h5ad.meta.json` beside the package. The checker identifies a prediction from that meta and refuses an
unidentifiable one, so R0 would have ended NOT READY on Oct 22 and arm1015 would never have fired. The order now matches
run_arm_generic's. The production 692 was re-armed.

**R7 (sealed reserve): independent check of [bf003b]'s frozen claim** (f134a63d, cbf443c3). The claim encodes the decided
recipe correctly:
- R7a F = M6+P2g;
- R7b final `cpm+c4h`, veto pair `cpm+c4h:cpm`;
- 553's `+c4h` = CLEAN4 τ 0.05, known-answer-tied to 517's h1_c4-c4t0.05 = mde_sign clean4:0.05.

R7b's access had started (23:42Z) before the check was requested, so that half of the check is post hoc; R7a was cleared
before its access. R7 is then spent on F plus the hook.

## 10. What cannot be measured offline (stated, so nobody fakes it)
- **No honest mock final truth exists on disk.** A new line's own controls → its own truth, Arc-like panel, full gene axis,
  ~400 cells per target, absent from the build. Checked 10-14:
  - K562 and Orion HCT116 / HEK293T are production donors;
  - H1 (VCC 2025) feeds the union and the programme option;
  - Jurkat is context A's line;
  - Jiang 2025 is DE tables only, and Broad-16 pseudobulk signs;
  - Nadig HepG2 has 9,624 of 18,533 genes, a median of 45 cells per target, essential genes only (underpowered, not
    Arc-like);
  - Replogle RPE1 (scPerturb) is an essential screen on a reduced gene set;
  - Sunshine Calu-3 is R7b (sealed, now spent).

  So the October path's absolute D/E/F score cannot be estimated offline. The evidence for it remains the A/B/C board
  (524dr2 0.2731), the bulk transfer reads on fresh panels (R6: retention 0.70; C10b: M6+P2g TRANSFERS; R7a pending) and
  R7b (the hook only).
- **The late switches ship without sealed confirmation:** UPRN (630, board PARTIAL +0.00093), DNMAG (631), and PIE (702)
  if adopted. R7's grammar cannot express them, as the ledger states before the access.
- **PIE's held-out reads, corrected by [bf003b]'s audit700:** Replogle has weight 0 and K562 (CVCL_0004) is not among
  PIE's 71 training contexts (45 Tahoe CVCLs, 23 Jiang, orion.hct116 / hek293t, arc_vcc_25.ARC_H1). So the K562 GWPS read
  is clean, and an unseen-context read, as D/E/F will be. Orion and Jiang truths are contaminated. H1 is clean of h1-test's
  responses but is a SEEN context (190 other H1 perturbations and its control means in train / val), so the h1-test read is
  likely optimistic for PIE relative to D/E/F. My first note ("contaminated on K562 and H1") over-claimed and is withdrawn.
- **The local instruments map to the board loosely:** H1 → board ratios run 0.07× … 2.3× and invert on several members
  (docs/H1-BOARD-CALIBRATION.md). A late DE switch's offline price is weak evidence of its D/E/F value.

## 11. Day-of cost of DNMAG (observed 10-15 00:47Z; exact numbers when the DN build ends)
The coordinator's DN build (630arm + VCC_DE_DNMAG=d_m200, run_arm_generic) was still in its build stage after ≥ 139 min,
on its last context; 630arm, the same base without DNMAG, took 74 min (`build_630arm/timings.json`, build 4,451.6 s).
So the DOWN-magnitude pass at least doubles a build. If de_dnmag is adopted, every def_day level build carries it. R0's
def_day (524dr2: 1 h 41 m, with three ~33–35 min level builds) could become ~4–5 h, moving R0's READY from ~09:00Z to
~noon on a 07:00Z release; final-phase arms and rehearsals cost the same. Flagged to the coordinator.
- **Measured (693, 10-15 00:55Z):** DNMAG costs 6.04 s per block on real H1 blocks (zs → UP → DN, peak 3.9 GB), so ~90 min
  per 900-block build. That's the DN build's overrun. Rehearsal windows and R0 with de_dnmag ON at v1: plan ~4–5 h per
  def_day. 693 has a byte-identical 3–4× faster design (per-row heaps for the sub-gate payback, batched rank lookups), to be
  built only if DN reads TRANSFERS or PARTIAL.
- **Exact (10-15 01:45Z):** the DN build's build stage took **9,727.9 s = 162 min** against 630arm's 4,451.6 s = 74 min:
  +88 min for 900 blocks, matching 693's 6.04 s per block. DN was READY at 01:11:52Z and submitted for 10-15 s1 at 01:12Z.
