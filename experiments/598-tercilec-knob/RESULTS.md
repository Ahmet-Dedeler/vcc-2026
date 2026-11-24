# 598 — `tercilec`, 493h's carried-only key tercile: an UNAPPLIED patch, 43/43 tests

On A/B/C, tercilec moves exactly the 4 targets 492/493h predicted:
- 98 / 104 / 98;
- E2F4 and SLC25A13 go HIGH → MID; EPHB2 and TMEM209 go LOW → MID.

`tercile` and `stercile` stay byte-identical, and unset stays verbatim.

2026-10-02, session "Explore contribution areas [57577c]", requested by the coordinator. No builds, no scorer. Peak RSS
0.30 GB. Outside the coordinator's windows.

**Class (AGENTS 11.1c):**
- This is a transfer-rule fix to the per-target sharpening (`sharpen_pt`, DECIDED 10-02 `tercile:2.5,2.0,1.75`).
- Held-out read (493h): efp, 3 October-like folds. The carried-only rule gives +0.84 % rel cos (3/3), against the plain
  rule's +0.48 % (2/3).
- Tries for this switch: 5 rules on 438's panels (plain, tieMID, tieLOW, carried-only, S_P2), 0 on the board.
- On A/B/C it moves 4 of 300 targets. The switch is about panels where non-carried targets are many: an October panel
  with thin donor coverage.

## What the patch does (`patch598.py`, 538/554 style)
Six files, applied AFTER the 10-05 bundle (538 → 554 → 590 → 598). It refuses:
- a tree whose def_day.sh lacks 554's mark;
- a second apply;
- a partial apply;
- any anchor that does not match exactly once.

### 1. `scripts/make_submission_pm_cd4.py`
- `VCC_BULK_SHARPEN_Q=tercilec:<low>,<mid>,<high>` parses like the two older forms: lowercased, three qs in (1, 3.5],
  needs `VCC_PERT_SCALE_BALLAST`.
- **`sharpen_q_by_tercile(..., carried=None)`.** In tercilec, the floor, the cut points and the tie blocks (511's floor-cut
  and cap rules) are read on **carried targets only**. Every non-carried target is then put in MID.
  - With every target carried it IS `tercile`, bit for bit: q and audit, 6 random keys plus the tie cases.
  - For tercile and stercile the carried mask is all-True, so every array the function builds is today's: q bytes and
    audits are identical on every synthetic case and on the real dump.
- **Audit (meta `bulk_sharpen_pt`).** It keeps today's keys and, for tercilec, adds:
  - `n_carried` and `noncarried_to_mid`;
  - `moved_vs_tercile`: {target: [plain tier, tercilec tier]}, which makes the KA checkable from the meta;
  - `carried_key`, `carried_key_sha256`, `key_ratio_c`, `key_ratio_maxdev`.
- **New `dump_key_carried(dump_env, …)`.** `carried` is read from the key that scaled **the loaded dump**: the dump env's
  `VCC_PERT_SCALE`, resolved through `PERT_SCALE_ARTIFACTS`, with the sha and the target order checked.
  - It is not this build's own `VCC_PERT_SCALE`. On A/B/C, 479p2g's level build ran `wshare_lo0.3_hi3` while its dump
    ran `lodo415_lo0.5_hi3`. On D/E/F the level builds run DEF_ENV's `wshare_lo0.3_hi3_def` while the dump runs the
    LODO id (def_day.sh 2b). That id is in the registry from def_day's pin until its exit trap.
  - **The artifact must reproduce the dump's own ratio:** ‖bulk_t‖/‖carriers_t‖ = c·w_t, constant to 1e-6 relative on
    every target with a carrier, or the build refuses. A/B/C reads 9.1e-15.
  - A wrong key (the build's wshare, the cap-6 LODO, wshare lo 0.5) misses by 0.5–2.7 and refuses. An unknown or
    absent id refuses. Without this check a wrong key would mis-assign the terciles silently.
  - The 49-target cap block (w = 3) and the w = 1 block cannot be told apart in r-space (493h). That is why the flags
    come from the artifact rather than being inferred from r.
- The LOAD block passes `carried` for tercilec only.

### 2–4. Day-of and driver files
- **`scripts/def_day.sh`** (the post-bundle 5851f84b): DEF_DAY_SHARPEN_PT accepts `tercilec:*,*,*`. A tercilec value
  needs `def dump_key_carried` in the builder, or def_day dies at start. Doc lines updated.
- **`scripts/def_day_helpers.py` meta_matches** and **`scripts/def_pipeline.py` september_rows**: both select the
  per-target forms with an explicit `("tercile", "stercile")` tuple. Without `tercilec` in it:
  - the helper float()s the value → rc 2, and the re-level loop stops;
  - the driver raises ValueError after a good h5ad exists (499's failure mode).
  Both are shown by T9/T10 on the pre-598 files.

### 5–6. Recipe and gate
- **`scripts/october_recipe.py`** checks sharpen_pt's form for all three modes, as the builder parses it. The check is
  stricter only on malformed values; the DECIDED `tercile:2.5,2.0,1.75` passes. A tercilec value also needs:
  - the builder's dump_key_carried;
  - def_day.sh's tercilec case.
- **`experiments/524-october-recipe-dr/gate524.py`**: a tercilec run must show:
  - `n_carried`, and `noncarried_to_mid` = total − n_carried;
  - `carried_key` = the LODO id def_day logged ("LODO key ON for the stack dump: <id>"), else the build's own key;
  - `key_ratio_maxdev` ≤ 1e-6.

### File shas
Canonical = HEAD 54b8886d + 538 → 554 → 590, rebuilt in a temp tree:

| file | before | after 598 |
|---|---|---|
| scripts/make_submission_pm_cd4.py | 4daae237c5fe4fd2 | 5cf0e0124452e1b1 |
| scripts/def_day.sh | 5851f84b9f238c93 | 2f18ce7aec0ab397 |
| scripts/def_day_helpers.py | e3ff45bcefc1a6d3 | fb715e3833f3f89d |
| scripts/def_pipeline.py | 9a718b9eaa24f10c | 76314af6b24959db |
| scripts/october_recipe.py | 503c45a00bf74164 | 9dbcd9149020f1e9 |
| experiments/524-october-recipe-dr/gate524.py | 7a43d2f5ae0f03c0 | b16e60c08e48e5d3 |

Full shas are in `work/test598.json`. They move if HEAD's copies move before the apply. post598 rebuilds them from HEAD,
so nothing here is a frozen table.

## Tests (`test598.py`, 43/43 PASS, `work/test598.log`)
The suite runs on temp trees only (the live tree is checked untouched). Its canonical tree skips any bundle step whose
`--check` says "already". All three do on an applied tree, so the suite keeps working after 10-05 commits the bundle.

- **Patch mechanics:**
  - `--check` writes nothing;
  - a second apply, a partial apply and a pre-bundle def_day.sh (b425c8bc's 0b53ba89) are refused;
  - py_compile ×5 and `bash -n` pass;
  - the day-of pin (`def_pin_panel_artifacts.rewrite`) commutes with the patch.
- **Synthetic arrays** (r = c·w exactly as `apply_pert_scale` scales, with 1e-14 noise):
  - floor block straddling the lower cut → all LOW;
  - 120-target cap block on the upper cut → all HIGH, also with 30 non-carried mixed in, and equal to 493h's rule;
  - efp-like 150/300 non-carried:
    - every non-carried target is MID;
    - carried targets split 50/50/50 with the floor LOW;
    - equal to 493h's w-space rule;
    - insensitive to the noise draw;
  - zero-carrier targets → MID (tercile and tercilec) or LOW (stercile, whose value 0 is finite: today's rule);
  - refusals.
- **The real A/B/C dump** (479p2g's `dump_h1x2_cd4x15_lodo05.npz`, 00ccfcba):
  - the old function reproduces 479p2g's recorded audit exactly (instrument KA);
  - new tercile and stercile equal the old;
  - tercilec gives 98 / 104 / 98, moves exactly the 4 targets, and equals 493h's w-space rule on the 8c82bef3
    artifact's w for all 300 targets.
  - Edges, tercile → tercilec: [0.30808, 0.90538] → [0.29975, 0.92510].
  - The moved targets:
    - EPHB2: r 0.3026, w 0.524;
    - TMEM209: r 0.3077, w 0.533;
    - E2F4: r 0.9117, w 1.580;
    - SLC25A13: r 0.9102, w 1.577.
- **The other files:**
  - env parse;
  - def_day.sh case lines run under bash;
  - meta_matches;
  - september_rows;
  - october_recipe `check`: clean with the patched pair; one problem with either unpatched file; malformed forms flagged;
  - gate524's row on 538's fixture: PASS; FAIL on a wrong key, maxdev, a count, a missing n_carried; no row for tercile;
  - `ka598.py` on synthetic metas.
- **`scripts/test_october_recipe.py` is stale on its own.** This needs the coordinator's attention.
  - On the live tree 1 test fails since 10-02's sharpen_pt decision (`test_live_state…` expects the old command).
  - On the post-bundle tree 5 fail: 538's new PENDING switches make rc 2 where the tests expect 0 or 3.
  - 598 changes nothing in it: the same 5 fail on the canonical and the patched tree (T11). It is not 598's to fix, and a
    checklist step that runs it would fail.

## On real day-of dumps (`analyse598.py`, `work/analyse598.json`)
This is the patched builder's own `dump_key_carried` and `sharpen_q_by_tercile` on def_day rehearsal dumps. Each
registry is rebuilt from the run's pin records, which is what the level builds saw.

| run | dump env key | carried | maxdev | tercile | tercilec | moved |
|---|---|---|---|---|---|---|
| 497dr (A/B/C relabelled, LODO pct = cap 6) | lodo_lo0.5_pct_def (c998b18c) | 294 | 7.9e-15 | 100/100/100 | 98/104/98 | the same 4 |
| 526 efp_inert (October-like, 410 key) | wshare_lo0.3_hi3_def (c3e429dd) | 295 | 8.0e-15 | 100/100/100 | 100/101/99 | 1 (SEPTIN7 HIGH→MID) |
| 510 efr (realistic coverage, LODO pct) | lodo_lo0.5_pct_def (dc7cbc83) | 300 | 9.1e-15 | 117/83/100 | identical | 0 |

- With all three donors carrying, tercilec moves 0–4 targets. Its effect grows with the number of non-carried targets.
  493h's efp gain came from folds with 153–156 non-carried. A D/E/F panel with thin donor coverage is where it matters.
- The day-of id resolves, and the ratio holds to 1e-14 on every real dump. On D/E/F the knob will not refuse for a
  benign reason.

## KA plan for the coordinator
Run after applying 598 to the live tree, before committing it. Three builds of ~40 min, serial, through 527's KA runner,
like 535's KAs.

```bash
REPO=~/Code/virtual-cell-challenge-2026; cd "$REPO"
uv run python experiments/598-tercilec-knob/post598.py                    # ALL PASS (light; also run it first)
NEWB=$(shasum -a 256 scripts/make_submission_pm_cd4.py | cut -d' ' -f1)   # the applied builder (5cf0e012... on today's HEAD)
DL=$REPO/experiments/417-bulk-sharpen/work/dump_h1x2_cd4x15_lodo05.npz; DLS=00ccfcba1b67fa24ff743ee0c64067e01158147ffe49c5b8ca37cef67600b667
REFZ=$REPO/submissions/446bj2z-cap15-partial50-clean4q80-pertscale3-ballast-sharpen175-stacklodo05load-ctxlevelm-mde010coarse-zr5-zs250.h5ad.meta.json
REFP=$REPO/submissions/479p2g-cap15-partial50-clean4q80-pertscale3-ballast-sharpenp2g-stacklodo05load-ctxlevelr-mde010coarse-zr5.h5ad.meta.json
clean(){ unset VCC_SAMELINE VCC_DE_PRIOR_GATE VCC_UNION_SOURCES VCC_RC_AUDIT_DUMP VCC_DE_MDE_SIGN VCC_RC_CALLSET_CPM VCC_RC_CALLSET_BLOCK VCC_DE_MAG VCC_BULK_CTX_SCALE; }
# KA598a: unset == 446bj2z (ec3a5ae6), 535's base exactly
( clean; export VCC_BULK_LOAD="$DL:$DLS" VCC_DE_MDE=0.1,0.1,0.3,coarse VCC_DE_ZERO_REPAIR=-5 VCC_BULK_CTX_SCALE="A=1.0962,B=1.1461,C=1.1681" VCC_DE_MAG=dn_k250_zs
  KA_Q=1.75 bash experiments/527-arc-hook-arm/run_ka_527.sh KA598a "$REFZ" "$NEWB" )
# KA598b: tercile == 479p2g (dd82740d); 479p2g's ARMED env (run479p2g.log line 1), no DE_MAG
( clean; export VCC_BULK_LOAD="$DL:$DLS" VCC_DE_MDE=0.1,0.1,0.3,coarse VCC_DE_ZERO_REPAIR=-5 VCC_BULK_CTX_SCALE="A=1.0957,B=1.1505,C=1.1807"
  KA_Q=tercile:2.5,2.0,1.75 bash experiments/527-arc-hook-arm/run_ka_527.sh KA598b "$REFP" "$NEWB" )
# KA598c: tercilec on the same env; run_ka_527 prints "KA FAIL" (the artifact must differ), the checker decides
( clean; export VCC_BULK_LOAD="$DL:$DLS" VCC_DE_MDE=0.1,0.1,0.3,coarse VCC_DE_ZERO_REPAIR=-5 VCC_BULK_CTX_SCALE="A=1.0957,B=1.1505,C=1.1807"
  KA_Q=tercilec:2.5,2.0,1.75 bash experiments/527-arc-hook-arm/run_ka_527.sh KA598c "$REFP" "$NEWB" )
uv run python experiments/598-tercilec-knob/ka598.py experiments/527-arc-hook-arm/work/build_KA598c/ka527-KA598c.h5ad.meta.json
```

**Pass conditions:**
- KA598a and KA598b: "KA PASS (identical)" in `experiments/527-arc-hook-arm/work/build_KA598{a,b}/runKA598{a,b}.log`,
  with driver rc 0.
- KA598c: `ka598.py` prints KA598c PASS. It checks three things:
  - `bulk_sharpen_pt` equals `work/ka598c_expected.json`: 98/104/98, the 4 moves, n_carried 294, carried_key
    lodo415_lo0.5_hi3 / 8c82bef3, maxdev ≤ 1e-12;
  - the config diff vs 479p2g is exactly `bulk_sharpen_pt`, plus `de_sign_emitter_sha256`, since the emitter moved from
    edb5b8dc to 0f56209d after 479p2g;
  - the output differs from dd82740d.
- KA598b is the first rebuild of 479p2g's artifact on builder 4daae237 / emitter 0f56209d. 535 showed byte identity only
  for 446bj2z / 446bj2zmz.

**Then:** commit the six files, with the KA hashes in the message. Then `october_recipe.py decide sharpen_pt
tercilec:2.5,2.0,1.75` if the coordinator adopts it.
- The 10-05 dress rehearsal resolves with the value decided at that time.
- A rehearsal with tercilec would show `moved_vs_tercile` in its meta: 4 targets on relabelled A/B/C, as 497dr.

## Files
| file | what |
|---|---|
| `patch598.py` | the patch (`--root <repo> [--check]`) |
| `test598.py` | 43 tests on temp trees (`work/test598.log`, `work/test598.json`) |
| `post598.py` | the post-apply check (live six files == HEAD + 598; identities on the real dump); validated ALL PASS on an applied scratch worktree, and FAIL on a one-byte change |
| `ka598.py`, `work/ka598c_expected.json` | the KA598c checker and its expected record |
| `analyse598.py`, `work/analyse598.json` | tercilec on three real day-of dumps |
| `work/builder.post598.py` | the patched builder text (5cf0e012), for reading |
