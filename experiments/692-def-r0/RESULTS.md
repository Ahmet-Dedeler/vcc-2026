# 692 — R0 on D/E/F: the readiness adapter (2026-10-15, session [17adab]; in progress)

`r0ready692.sh <bundle> <r0 work> <stem> [CAP]` turns a finished def_day R0 (691's AUTO_R0 build, or NEXT.txt run by hand)
into the evidence the coordinator's slot waiter `experiments/670-coord-1015/arm1015.sh` gates on, so R0 is submitted by
the same tested path as every slot (536 retry + resume, novelty record, capture, NOGO). It never submits.

Checks, each fail-closed (NOT READY + END): def_day finished (trap line + no process), postbuild PASS, DONE, no STOP, the
builder restored to its pre-pin sha; 630dd FAILSAFE_TOTAL 0 when present; the final build not past its cliff
(def_day_helpers past_cliff); the real `vcc prep` under the build gates; the novelty check. Provenance (warning only):
690c's builder.pinned.py == the final meta's builder sha. With COMMIT692=1 (default) it commits the round's small records
with an explicit pathspec once READY.

Day-of use (coordinator): arm1015.sh with ready-log `experiments/692-def-r0/work/r0def.ready.log`, run-log
`…/r0def.run.log`, package `submissions/r0def.vcc`, TAG r0def, FLOOR 202610220700.

Test so far: on 524dr2's finished run (its h5ad swept) the parse steps read `[524dr2] PASS` and `CLIFF OK final 524dr2-lvl2:
D 0.9735 E 0.9794 F 0.9885`, then NOT READY on the missing h5ad (fail-closed). The full path (prep, novelty, provenance)
runs inside 690's R-SAME.

## The day: D/E/F released 2026-10-22 00:04Z (readiness [401d9b]; OCT22-CHECKLIST §C1–§C2)

- **Release:** `controls-final` became available at 00:02Z (final phase open). watch691 detected it at 00:04:26Z and downloaded
  `vcc_2026_controls_final.zip` (643,356,568 bytes, sha256 d625f47f56a6b582…) to `data/def-download/20261022T0004Z/`. It unpacked
  6 files to `data/challenge-def`. vcc 0.3.0 (PyPI 0.3.0, no UPGRADE line).
- **Bundle:** panel `vcc2026-test-1`, contexts D/E/F, 300 targets × 400 cells, 18,400 controls per context, max library 53,865 / 51,808 /
  52,365. Preflight rc 0, no AUTO OVERRIDE.
- **Branch NEW-PANEL:** 129/300 targets shared with A/B/C; gene axis identical (18,533 genes, so no 531 F6 stop).
  RPE1 ∪ Jurkat essential-screen targets in the panel: **0** (A/B/C 0, lab prediction ≤ 3). A/B/C reads and LABC remain the right
  instruments.
- **Identity** (691 light, 8 methods; reading rule DEF-RUNBOOK top block):
  - D = **MOLM13** (AML, Myeloid), 8/8.
  - E = **U87MG** (glioblastoma, CNS), 8/8.
  - F = **SKLMS1** (leiomyosarcoma, Soft Tissue) by consensus, but **not unanimous**: top-1 votes CALU1 3, OS252 2, SKLMS1 1,
    CCLF_RARE_0020_T 1, A1207 1. F's lineage call is weak.
- **AUTO_R0 did not start.** The light stage `cd4-extract` stopped with "only 287 panel targets map across any condition; require
  290" (047's union gate). watch691 wrote NEXT.txt and stopped (00:05:09Z).
- **Mapping confirmed** (DEF-RUNBOOK 0.0 / §5, read-only on GWCD4i obs; `experiments/690-def-readiness-1015/work/def-day/mapping_def.json`):
  - A/B/C reproduces its known answer, 291/300.
  - D/E/F maps 287/300. Missing: ANG, C5orf22, DNAJC25, DPY30, GPATCH1, LENG1, LTBR, LTO1, RAB26, REV1, SURF4, TYRO3, VPS36.
  - No case-variant hits. The only substring near-hits are different genes (ANGEL1/2, ANGPT*, DNAJC25-GNG10 readthrough,
    GPATCH11).
  - So the gate firing is a real coverage fact, not a mapping bug. C5orf22, LENG1 and REV1 were missing on A/B/C too.
  - The agent-inspection session [2ed02f] found the same independently.
- **R0 launched by hand, 00:10:24Z** (def_day pid 17996) with the one fix, `DEF_DAY_CD4_THRESHOLDS=relaxed` (481 F1; rehearsed on
  efp, union 186):
  - Recipe: NEXT.txt's, re-resolved from the live decisions file with OCTOBER_DECISIONS unset (`october_recipe.py command --json`
    identical, problems []; `…/def-day/r0-recipe-relaxed.json`).
  - Work dir: the fresh `experiments/692-def-r0/work/r0`, which r0ready692 (5060) and arm1015 R0 (5181, floor 10-19) already wait on.
  - Gates at launch: builder e99bbd4a == HEAD; decisions 04be2af3, def_day 11ff102e, emitter 0f56209d, helpers 5354e449; no other
    builder or chain; free+inactive 26 GiB; memory_pressure 80 % free; disk 125 GiB.
  - The census floors (union ≥ 120, complete ≥ 100) still gate.
- **Re-level on D/E/F** (target rho* 0.9793, band ±1 %, cliff rho 1; ~75 min per build). E (U87MG, Cp46 sd 1.336 against D 0.500 and
  F 0.634) is far steeper than A/B/C's contexts:

  | build | factors D / E / F | rho D / E / F | outcome |
  |---|---|---|---|
  | lvl0 | 1 / 1 / 1 | 0.8247 / **1.5555** / 0.8777 | deleted: E past its cliff |
  | lvl1 | 1.0898 / 0.7935 / 1.0563 (ntclevel EV scales) | **1.0151** / 0.5822 / **1.0800** | deleted: D and F past their cliffs |
  | lvl2 | 1.0737 / 0.8968 / 1.0294 (secant; E's local exponent 4.25) | 0.9762 / 0.9319 / 0.9749 | **SAFE**; D and F in band |

  The first safe build exists, so R0 ships (best_safe at the worst). lvl3 holds D and F and steps E to 0.9063 (local exponent
  4.70). Plan B (`experiments/690-def-readiness-1015/planb_r0.sh`, 42dea53f: fresh dir, DEF_DAY_RHO_SCALE 0.98) is not needed.
- **lvl3 CONVERGED** (05:17:54Z): D 1.0737 / E 0.9063 / F 1.0294 → rho D 0.9762 / **E 0.9764** / F 0.9749, every context in band.
  Final output `fa83298813978cd8…` (r0def_final.h5ad, kept). 690h DNMAG: no LOW context. 630dd uprn FAILSAFE_TOTAL 0.
- **Postbuild STOP (05:17:56Z), 41/42 rows pass.** The one failure is 554's emitter steering-dose row: D's median `prior_lambda`
  0.8975 against the 0.9 floor (margin −0.0025). D's other numbers are healthy (MDE 0.954, missed calls 1); E and F read
  prior 1.0 / MDE 1.0, missed 24 / 39.
  - The guard was built to catch collapses (dose 0–0.44: 429a, 458*, 465efr). A/B/C's prior always read 1.0.
  - So D (MOLM13) gets ~10 % less CLEAN4-prior steering on its median block, a capacity effect of D's library classes.
  - 554 RESULTS prescribes exactly this: 0.85–0.9 is a ~10–15 % dose loss; "stop, look, and override with
    DEF_DAY_ALLOW_LOW_DOSE=1 if nothing else is wrong".
  - Not yet run on this build: the zero-gene audit and the prep dry run (they come after postbuild in def_day).
- Consequences: r0ready692 wrote NOT READY (05:18:13Z) and arm1015 R0 logged NOT SUBMITTING + END (05:18:14Z). wait708 read R0
  as done and started day708 at 05:21:08Z, as the coordinator ordered.
- **Prepared, not launched:** `planb_r0.sh` gains `PLANB_MODE=lowdose`.
  - It reruns the same recipe in the fresh r0b with `DEF_DAY_ALLOW_LOW_DOSE=1` (rho-scale off, max_levels 4).
  - It refuses unless postbuild's only failing row is 554's, every lambda is ≥ 0.85, and day708 is no longer running.
  - `ka_r0b.sh` touches NOGO-R0 unless r0b's final output is fa832988… byte for byte. The chain is deterministic: R-KA2 rebuilt
    R-ABC2 in another dir.
  - Dry run clean. ETA ≈ 5 h of builds after day708.
- **The coordinator's call (05:2xZ):** apply 554's preregistered override by **re-gating the finished build** (the 697/698/673
  precedent), not a rebuild. `experiments/690-def-readiness-1015/regate_r0.sh` started 05:27:15Z with evidence in
  `experiments/692-def-r0/work/r0g-ready/`. All PASS by 05:27:37Z:
  - **KA:** def_postbuild re-run exactly as def_day's step 5b, with `--builder r0/builder.pinned.py`. a397abb1 == the meta's
    builder, so the production builder is untouched. It reproduces def_day's 42 rows identically (41 pass).
  - **Override** (`DEF_DAY_ALLOW_LOW_DOSE=1`): exactly one row changes, 554's, now "ALLOWED"; 42/42 pass.
  - **def_day's step 6:** the zero-gene audit reads 0 in 0 blocks (lfc < −5: 17, > +3: 11). The driver's prep stage is only
    `vcc prep --dry-run`, which the real prep below covers.
  - FAILSAFE_TOTAL 0; DNMAG LOW none; CLIFF OK lvl3.
  - Then it runs r0ready692's real `vcc prep` behind the same gate (waits out day708's fold fits), novelty, and READY into a fresh
    ready log. The coordinator re-arms arm1015 R0 on it.
  - The lowdose rebuild mode of planb_r0.sh stays unused.
- **R0 SUBMITTED, 2026-10-22 06:05:02Z: entry `k2OKr0zPvJq48BCtdUx8`** (final status published, partition test, panel
  vcc2026-test-1, anchors testD/E/F-r1, `is_final` true, members hidden in the final phase).
  - Package `submissions/r0def.vcc` sha256 5c82ea9e…, 2,671,319,040 bytes. vcc prep rc 0 in 24 s (13.7 GB RSS).
  - Prediction identity fa83298813978cd8… (r0def-lvl3; re-gate READY 05:37:26Z, records 31fe4efd).
  - Upload 05:38:15Z → 06:05:02Z (27 min). vcc 0.3.0's `--wait` ended on "published", with no hang.
  - Model name from the queue ('I am persistent but open-minded.'), `-d` empty. Captured to
    `data/board-entries/k2OKr0zPvJq48BCtdUx8.json`, novelty-ledgered by arm1015.
  - **What R0 is:** the frozen October recipe, decisions 04be2af3 (H1 ×2 / CD4 ×1 stack, LODO key strength 0.25,
    tercilec sharpen, DE riders dn_k250_zs / b_k100_q20r3p10 / d_m200 / clean4:0.05). Plus two day-of, preregistered
    contingencies: DEF_DAY_CD4_THRESHOLDS=relaxed (CD4 union 287 < 290) and 554's DEF_DAY_ALLOW_LOW_DOSE=1 (D prior dose
    0.8975), applied by re-gate.
  - Re-level converged at lvl3: D 1.0737 / E 0.9063 / F 1.0294.
- **day708 (wait708, 05:21 → 06:09:50Z, rc 0, 48 min)** re-priced the panel-dependent switches on D/E/F's real panel with
  held-out public folds. Decisions (`experiments/708-def-panel-repricing/work/decision708-d708def.txt`):
  - WEIGHTS move to h3c1.5 (H1 ×3, CD4 ×1.5; production-key confirmation +0.00177, floor-kept +0.00172).
  - Rider 705 ON (+0.00032).
  - Swap 707 ON (+0.00182).
  - These feed the coordinator's 6a replacement decision (X1/X2). R0 stays the final entry until a replacement passes every
    6a line.
- **X1** (06:13Z, experiments/692-def-r0/work/x1, driver x1_690.sh 2e3778cf) = R0 recipe + H1 ×3 + CD4 ×1.5 + hook707. Its log
  confirms the stack dump at VCC_H1_UNION_G 0.84 / VCC_CD4_TARGET_G 0.5165… and the 707 hook replacing it (3754b49a).
  `chain_x1x2.sh` (cee075de) settles it, including the 6a (d) 554 re-gate, then runs gate6a.
- **X2 = C0** (SIXA 1ebd3a1f) = R0 + dumphook811c0.
  - **Blocker found 07:0xZ:** hook811 refused D/E/F because 136 targets were missing from 811's banks (pin 37e34335 covered
    ABC + L1–L6 only).
  - The lab extended the banks to DEF (manifest 4cc6d82a) and re-pinned hook811.py (c7cff26d, eb79f7f5). The ruling (SIXA
    344c4c31): data inputs only.
  - **LZ PASS: C0 − B0 +0.00999 [+0.00659, +0.01273], 6/6** (CONFIRMATION-LEDGER 79868244).
  - **Preflight PASS 08:10:40Z:** the hook on r0's stack dump gives output 3b4f0c05… (carriers / perts / genes / env equal;
    bulk changed on 299/300 rows, median relative change 0.183).
  - Open: banks/manifest.json is uncommitted (disk 4cc6d82a, HEAD 37e34335). X2 builds after X1 once the coordinator writes GO-X2.
- **X1 READY + 6A PASS (10:1xZ).**
  - def_day CONVERGED at lvl2: D 1.0835 / E 0.9557 / F 1.0681 → rho D 0.9783 / E 0.9715 / F 0.9703. lvl0 deleted on E 1.1730,
    lvl1 deleted on F 1.0745. Output d0d4bd9c….
  - Postbuild failed only 554's row, with D prior_lambda 0.8975 again (identical to R0's; MDE 0.954, missed 1). So chain_x1x2 ran
    the 6a (d)(2) re-gate: KA reproduced the 42 rows, the override flipped only 554, zero-gene 0 in 0 blocks, prep rc 0, novelty.
  - **READY x1def 3528fb51fcc8f35a… in experiments/692-def-r0/work/x1g-ready.**
  - gate6a (rerun `work/gate6a_x1.rerun.out`, json `x1g-ready/gate6a_x1.json`), all PASS:
    - the diff vs R0 is {bulk_ctx_scale, bulk_load, bulk_load_sha256, bulk_sharpen_pt};
    - bulk_sharpen_pt moved only edges / key_ratio_c / key_ratio_maxdev;
    - the stack dump is at H1 0.84 / CD4 0.5165…; bulk_load_sha256 is the 707 hook's dump 3754b49a, not the stack dump 7fd7f98e;
    - in band; implementation and emitter == R0's.
  - **Three judge fixes in gate6a_x1.py before that PASS** (none touches the build). The chain's first run crashed on (1):
    1. relative paths are resolved against the repo;
    2. the meta records bulk_load as a file NAME, so the sha decides;
    3. (e): a pinned builder can never byte-equal R0's across work dirs, since the pin embeds the work dir and each
       cd4_compact/*/meta.json sha, which records a measured `peak_rss_gib` (0.3128 vs 0.3096). It is now judged by content:
       identical after normalising those, the cd4 metas differ only in peak_rss_gib, and every other implementation entry is
       byte-equal.
    The R0-vs-R0 negative control still FAILs every X1-specific row.
- **X2 launched 10:11:55Z** by chain_x1x2 in experiments/692-def-r0/work/x2 (x2def): GO-X2 08:13:31Z; the chain's own preflight
  reproduced 3b4f0c05 (ref `work/def-day/hook811c0_ref-101149`). r0ready692 → x2-ready.
- **Disk sweep (11:3xZ, the coordinator's yes): freed 21 GiB (77 → 98 free).**
  - Removed both names of each final prediction h5ad: r0 (r0def_final.h5ad = r0def-lvl3.h5ad, inode 262099155, 11,516,333,144 B)
    and x1 (x1def_final.h5ad = x1def-lvl2.h5ad, inode 262174855, 11,524,721,752 B).
  - Before the rm: each package matched its logged READY sha (r0def.vcc 5c82ea9e…, x1def.vcc 3528fb51…), each inode had exactly
    2 links, and lsof showed no reader.
  - Kept: every meta, npz, cd4_compact, ntclevel, hook dump, log and both packages. x2 untouched.
  - Restore (content, not the original byte layout):
    `uv run python scripts/restore_prediction_h5ad.py submissions/r0def.vcc <out.h5ad>`, and the same for x1def.vcc.
- **X2 = C0 (R0 + dumphook811c0): READY + 6A PASS, and SHIPPED as the replacement final** (coordinator's 6a decision, SIXA 63803087).
  - Re-level:

    | build | factors D / E / F | rho D / E / F | outcome |
    |---|---|---|---|
    | lvl0 | 1 / 1 / 1 | 0.8248 / **1.6792** / 0.9247 | deleted |
    | lvl1 | 1.0897 / 0.7637 / 1.0291 | **1.0096** / 0.5381 / **1.0329** | deleted |
    | lvl2 | 1.0757 / 0.8801 / 1.0150 | 0.9767 / 0.9172 / 0.9781 | safe; E out of band |
    | lvl3 | 1.0757 / 0.8924 / 1.0150 | 0.9767 / 0.9757 / 0.9781 | **CONVERGED** (15:28:38Z) |

  - Postbuild failed only 554's row, with D prior_lambda **0.8994** (≥ 0.85). chain_x1x2 ran the 6a (d)(2) re-gate:
    **READY x2def e39d2ba6ee716442… in experiments/692-def-r0/work/x2g-ready** (15:30:47Z). Prediction identity 94ae91b6….
  - gate6a --candidate x2, by the frozen judge (gate6a_x1.py sha 360f3125, 777b6db6), json `x2g-ready/gate6a_x2.json`, **PASS**:
    - the diff vs R0 is {bulk_ctx_scale, bulk_load_sha256};
    - X2's stack dump == R0's ada1baa2 by content; the used dump (bulk_load_sha256 aeb99967) == the 811 hook on R0's stack dump
      by content (ref `work/def-day/hook811c0_ref-101149`);
    - doses 0.56 / 0.3443…; in band; same tree (pinned builder 7a756472 ≡ R0's a397abb1 after the accepted normalisation).
    - NOTE: bulk_sharpen_pt's data fields moved (edges / key_ratio_c / key_ratio_maxdev, ~1e-15 per the coordinator).
      Behavioural fields are identical, so it was accepted under the W1 rule.
  - **Submitting since 17:49:09Z** (arm1015 R1, model name from the queue, `-d` empty).
  - **X1 retired, never submitted** (x1g-ready, 3528fb51).
  - Keep x2, r0 and x1's metas and npz until the final closes.
- **X2 SUBMITTED, 2026-10-22 18:20:14Z: entry `S5RxOrm4ZbqtnPRIX8NX`** (published, partition test, panel vcc2026-test-1, `is_final`
  true, members hidden). Package x2def.vcc e39d2ba6…; upload 17:49:09Z → 18:20:14Z (31 min). Ledgered and captured by arm1015 R1.
  **This is now our last final entry** (it replaces R0 k2OKr0zPvJq48BCtdUx8; only the last counts).
