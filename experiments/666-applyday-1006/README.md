# 666: the "10-06 APPLY DAY v2" run (coordinator; starts 10-05 after the 18:00Z freeze)

The procedure is `docs/NEXT-SLOT-DECISION.md` "10-06 APPLY DAY v2". Pins: the LAST PINNED block in 606 RESULTS (10-02 14:56Z,
85d4be68). Preflight on the post-10-05 HEAD (05bb393d): test598 43/43, compose606 ALL PASS, final shas == pins.

- Steps 1–4b (by hand): 598 → 612 → 612dd → 630 → 630dd → 642 → 604 → 644 v3, each `--check` then apply, each with its suite.
- Step 5: `experiments/410-pertscale-revisit/run_builder_tests.py` → 56 passed / 0 failed / 1 skipped.
- Step 6: `NEWB=<post-apply builder sha> nohup bash experiments/666-applyday-1006/run_kas.sh`. It runs six KAs, serial (KAu, KAb,
  KAc, KAh, KAg, KAr), with commands copied verbatim from 598's and 612's KA plans, plus KA630 per step 6 and 606. It ends with
  `KAS ALL PASS` or `KAS FAIL <tags>` in `work/kas.log`.
- Step 7: run_arm_generic driver list += the post-598 def_pipeline sha (76314af6 per 650); `bash -n`; test_gate591.
- Step 8: commit with a pathspec (KA hashes in the message); send [537757] the final shas.
- Step 9: `october_recipe.py decide sharpen_pt tercilec:2.5,2.0,1.75` once KAb and KAc pass (C10b f6aa2ed9 + 599 7db7938c).

## Run log (2026-10-05, from the 18:00Z freeze)
- Freeze check: 606 RESULTS unchanged since cbaf7085 (timestamps only); last PINNED block 14:56Z; nothing building; tree clean.
- Step 1, 598: `--check` clean, test598 43/43 (pre-apply), applied, post598 ALL PASS (tercilec audit 98/104/98, the 4 moves);
  def_pipeline 76314af6.
- Steps 2 / 2b, 612 / 612dd: test612 23/23 before and after; 612 applied (3 builder edits), 612dd applied (6 def_day insertions).
- Steps 3 / 3b, 630 / 630dd: applied (rank_up.py 7e6b3e52 installed); test630 35/35; test630dd 54/54 (the checklist's "50/50" is
  now 54); test_rank_up630 seed 0 b_k100_q20r3p10 85 PASS bit for bit (265 s).
- Step 4, 642: applied cleanly. Its test cannot run on a tree where 642 is already live: setUpClass rebuilds a tree from the live
  files and re-applies 642, which refuses "already patched". P1's skip rule covers the 10-05 bundle patches, not 642's own apply.
  So it was run unchanged in a temp git worktree at the pre-apply HEAD (e630746b), with the venv, data/, submissions/ and the
  untracked inputs linked in (590's 28 profiles; 166's nadig_hepg2 donor). Result: 6/6 OK, worktree removed. **Checklist fix:**
  run test_hek293_642 BEFORE 642's apply.
- Step 4a, 604: the def_sameline WARNING is benign as documented; 35/35 and 5/5. Step 4b, 644 v3: OK. `bash -n` def_day.sh OK.
- Pins: 12/12 live shas == the PINNED 14:56Z block (def_day 54fcf040, def_sameline 5cd3ab07, def_identity 7073406f, def_day_helpers
  fb715e38, def_pipeline 76314af6, october_recipe 9dbcd914, builder fe4678c3, def_preflight 17f0836f, def_postbuild 557c2b52,
  gate524 b16e60c0, def_identity_prior 428f4a32, rank_up 7e6b3e52).
- Step 5: run_builder_tests 56 passed / 0 failed / 1 skipped.
- Step 7 (done early; it touches no build): run_arm_generic driver list += 76314af6; `bash -n` OK; test_gate591 ALL PASS.
- Step 6: run_kas.sh armed 18:18:35Z (pid 26113, NEWB fe4678c3). Steps 8 and 9 follow `KAS ALL PASS`.
