# 499 — adversarial review of the October day-of code (2026-09-27)

Code only, on `cloud/499-october-code-review` at `40c164e0`. No challenge controls, no `vcc`.
Repro: `bash experiments/499-october-code-review/repro/run_repro.sh` (log: `repro/last-repro.log`).
This VM is GNU bash 5.2.21. No docker/podman, so nothing below is a reproduced bash-3.2 failure.

No critical finding. Nothing here links `_final` to a build whose factors or q disagree with `def_day.final.env`. The two high findings stop the day; they do not publish a different build.

| id | severity | status | where | what happens on the day |
|---|---|---|---|---|
| F1 | high | CONFIRMED | `scripts/def_pipeline.py:206`; `experiments/497-october-switches/work/def_day.497.sh:214-218` | `DEF_DAY_SHARPEN_PT` makes every level build's arm assert raise `ValueError`. The day deletes that h5ad and exits. `meta_matches` and the link never run. |
| F2 | high | CONFIRMED | `experiments/496-uncalled-sign-hook/patch_builder496.py:17-19`; `def_day.497.sh:289` then `:302`; `scripts/def_postbuild.py:456-460,524-525` | Applying 496 rewrites `de_sign_emitter_sha256` on every build, knob unset included. The day links `_final`, then rule-10 exits 1. No `DONE`. |
| F3 | low | CONFIRMED | `def_day.497.sh:321` | On a per-target build the closing log prints `q None` and omits `bulk_sharpen_pt`. `final.env` still has the form. |

## F1 — per-target sharpening crashes the arm assert, then the day deletes the build

**Severity:** high. **Status:** CONFIRMED.

`september_rows` always does `scalar("VCC_BULK_SHARPEN_Q", "bulk_sharpen_q", float)` (`scripts/def_pipeline.py:206`). `assert_build_carries_arms` calls it after every non-dump build (`scripts/def_pipeline.py:704`; 485's copy still does, at `experiments/476-def-day-merge/repro/def_pipeline.py.new:718`). The day exports the per-target string before those builds (`def_day.497.sh:198`, `SQ` from `:82-87`).

`float("stercile:2.5,2.0,1.75")` raises `ValueError`. That is not a failed row: it is an uncaught exception, driver exit 1. `build()` then deletes the h5ad and meta it just required to be non-empty and `die`s (`:216-218`). `die` is `exit 1`, so `best_safe` does not run.

Until 479 is in the builder the day dies earlier, at the `sharpen_q_by_tercile` grep (`:85`). HEAD has no such function (`make_submission_pm_cd4.py:1056-1058` is `float()` and `(1, 2]`). Once 479 is applied, the builder can write a good h5ad whose meta has `bulk_sharpen_pt` and no scalar `bulk_sharpen_q`, and this assert still crashes. `meta_matches` would accept that meta.

Unset, `SQ="$Q"` (`1.75`) and the float row stays `('bulk_sharpen_q', 1.75, 1.75, True)`.

**Repro** (`repro/last-repro.log`):

```
ValueError: could not convert string to float: 'stercile:2.5,2.0,1.75'
scalar 1.75 row: ('bulk_sharpen_q', 1.75, 1.75, True)
after driver rc=1 h5=nonempty
deleted h5=gone meta=gone
meta_matches ... stercile:2.5,2.0,1.75
exit 0
```

After `patches/def_pipeline.py.patch` (applies to HEAD and to the 485 copy): a matching `bulk_sharpen_pt` row is `ok=True`, scalar `1.75` is unchanged, a tercile meta against a stercile request is `ok=False`, and an empty env still expects `bulk_sharpen_q` absent.

## F2 — 496's emitter pin fails rule-10 with the knob off

**Severity:** high. **Status:** CONFIRMED. The patches are not on HEAD. This fires only if they are applied before the day.

`patch_builder496.py` edit 1 replaces `DE_SIGN_EMITTER_SHA256` with the sha of the patched emitter on every apply. The config key is top-level (`make_submission_pm_cd4.py:6752`, inside the `DE_SIGN_PRIOR` spread). `de_mde_sign` is written only when `DE_MDE_SIGN is not None` (edit 6). The day wipes every `VCC_*` at `def_day.497.sh:64` and never sets `VCC_DE_MDE_SIGN`, so the knob cannot leak in. The pin still moves.

Measured on temp copies: HEAD emitter `edb5b8dcbee519683b9d50c3556d4470ca6290594d4d1611be3d55c4fb7d41c7` (same value as `experiments/455-driver-q175/work/ka446b/455-ka-446b.h5ad.meta.json` and `submissions/460dr-october-path-lodo-abc-q1.75.h5ad.meta.json`). Patched emitter and the new pin are both `1658ae2c400f75ee1659a8ffcaebca569bf3c8afd55ffd1c3006414344ef2791`.

Rule-10 (`def_postbuild.py:456-460`) treats that one key as moved. Today's `ALLOW` (`def_day.497.sh:298-301`) does not name it. `unexpected` is non-empty, so postbuild `SystemExit(1)` (`:524-525`). The hard link of `_final` is at `:289`, before that check. The script then dies, so there is no `DONE` line. Counts are the pre-496 counts: with `mde_add_sign` absent, `move_m` and `mde_rec` match the unpatched block (see clean list).

The default incumbent path `submissions/446b-cap15-….meta.json` is not in this checkout. The two metas above record the same sha as the current pin, which is what rule-10 compares.

**Fix:** `patches/def_day.497.sh.patch` adds `de_sign_emitter_sha256` to `ALLOW` only when the builder source pin differs from the incumbent's recorded sha. On HEAD they match, so the line adds nothing:

```
src=edb5b8dcbee51968 inc=edb5b8dcbee51968 ALLOW=whitener_id,pert_scale,pert_scale_sha256,bulk_load,bulk_load_sha256,bulk_ctx_scale
src=1658ae2c400f75ee inc=edb5b8dcbee51968 ALLOW=...,de_sign_emitter_sha256
```

A config sha that moves while the source pin still equals the incumbent still fails the gate. `bash -n` on the patched script exits 0. Both the quote-escaped pattern and the `.` pattern extract the pin on this Python; the patch uses `.` so the shell line has no escaped quote.

## F3 — closing log drops the per-target form

**Severity:** low. **Status:** CONFIRMED.

`def_day.497.sh:321` prints `bulk_sharpen_q` only. On a PT meta that key is absent, so the log says `q None`. `def_day.final.env` is written from `$SQ` at `:291`, so the file is right. The same day-script patch adds `| pt`.

```
OLD meta: bulk_ctx_scale {'D': 1.0, 'E': 1.0, 'F': 1.0} | q None | sameline null
NEW meta: bulk_ctx_scale {'D': 1.0, 'E': 1.0, 'F': 1.0} | q None | pt {'mode': 'stercile', 'q': [2.5, 2.0, 1.75]} | sameline null
```

## Checked and clean

- **No wrong-ship path.** `meta_matches` (`def_day.497.sh:286-287`) runs before `ln -f` (`:289`). On `stercile:2.5,2.0,1.75` with matching `bulk_sharpen_pt`, no scalar q, and `D=E=F=1`, it exits 0. F1 deletes the h5ad before that. F2 links a count-identical build and then stops before `DONE`.
- **476 merge.** `def_day.resolved.sh` has no conflict markers. `bash -n` exits 0. Still present: `UNIT` from the manifest (`:54`), `--cd4-thresholds "$CD4T"` (`:94`), delete-before-build, driver rc check, `past_cliff` rc 2 = `die`, `[bulk-dump] wrote`, `best_safe`, `meta_matches` before the link, rule-10 before the zero-gene audit (`ZERO block total` gate). 485's driver hunk does not touch line 206; F1's patch applies to `repro/def_pipeline.py.new`.
- **497 knobs unset.** `SQ` stays `1.75`. `${DEF_DAY_LODO_HI_PCT:+--hi-pct …}` expands to zero words when unset and when empty (`hi_argc=0`). `LID` stays `lodo_lo0.5_hi3_def`. `[ "$SQ" != "$Q" ]` is false, so the allow list does not gain sharpen keys. The dump is a subshell (`:185-189`) and the sharpen export is after it (`:198`).
- **497 knobs set.** `--hi-pct` is two words (`--hi-pct` `96.64697522769474`). Pin id `lodo_lo0.5_pct_def` matches `scripts/def_pin_panel_artifacts.py:182` (`[a-z0-9_.-]+`). File name is `pert_scale_lodo_def_lo0.5_pct.npz` (`:170`). `def_lodo_key.py` puts `hi_pct` in the json sidecar only (`:124`); the npz is `targets,w,key,carried,n_donors` (`:120-121`). The pin hashes the npz. The shell passes P through; it does not snap 96.647 to cap 6.
- **Factor grep** (`:243`) returns `D=1.01,E=0.99,F=1.1` for uppercase labels. A lowercase label matches only the empty tail (`lower []`), which takes the "no factors" branch and ships the nearest safe build. October labels are D, E, F.
- **PRIOR_GATE / programme / same-line.** Independent of the new knobs. Programme plus `noprog` dies (`:136`). Prior gate without a programme dies (`:140`). Programme with 0 contexts dies (`:157-158`). Prior-gate export is after the dump and requires `SL` (`:200-202`). 5b adds `de_prior_gate` only when set, and the sharpen keys only when `SQ!=Q`. Both can be on together.
- **496 emitter, key absent.** Executing the `mde_add` block from HEAD `rc_emit_k.py` and from the patched file on the same synthetic arrays: `move_m equal True`, `mde_rec equal True`, keys unchanged (no `mde_n_sign_*`). With `mde_add_sign` present, `move_m` differs and the diagnostic reports `mde_n_sign_cover` 5, `mde_n_sign_flip` 4.
- **496 gene axis.** `patch_builder496.py:56` indexes `P_down[rc_control.keep]`. `keep` is the boolean full-gene mask; `ref_lin` is `cpm[keep]` (`src/vcc2026/rc_emit.py:752-767`). On a synthetic axis, `P_down[keep]` equals `P_down[flatnonzero(keep)]` with `equal_nan=True` (naive `array_equal` is false only because of NaNs) and the shape equals `ref_lin`. `clean4_prior_challenge.npz`: genes and `P_down` are `(18533,)`, 6346 finite, and `psign` agrees with `where(P_down>0.5, -1, 1)` on all 6346. The hook's non-finite → 0 is the 496 patch, not `psign`.
- **bash 3.2.** No `declare -A`, `mapfile`, or `${var,,}` in the merged day script. The new shell is `:+` and `[ -n ]`, same form as the existing `--force` / `--programme-source` lines. `POINTS` is appended inside `build()` before any `"${POINTS[@]}"`; a failed build `exit`s first. Not reproduced on bash 3.2.

## Coordinator's disposition (2026-09-27 06:45Z)
All three findings accepted; F1 was also live on the production path, not only the day script.
- **F1 (high) — confirmed on our own logs, and it reached further than the day script.** 479's KA-T driver run logged
  `ValueError: could not convert string to float: 'tercile:1.75,1.75,1.75'` (rc 1); the KA passed only because the KA
  runner reads the meta written before the assert. `run_arm_generic.sh` runs build + postbuild in one driver call, so the
  10-02 pair (chain479b), 498's combined arms, 497dr and 497e would all have stopped at "UNEXPECTED postbuild failure count".
  - Fix: this review's `patches/def_pipeline.py.patch`, applied as driver **44714b9a** (cb9e432a). Checked on KA-T's and
    KA-S's real metas: the tercile row is ok, a stercile request against a tercile meta fails, the scalar row is unchanged.
  - `run_arm_generic.sh` pins 44714b9a (inserted before 62b5110d, so chain476's anchor still matches once).
  - 476's and 497's `repro/def_pipeline.py.new` are now the fixed driver + 485's patch (**9b2c7eca**). 485's
    `repro_fixes.py` passes on both.
  - chain479b was killed before it started and re-armed behind chain495, via `relaunch479b.sh`: chain495's gate regex names
    chain479b.sh, so a waiting chain479b would have deadlocked it.
- **F2 (high) — accepted into 497's day-script patch.** 5b adds `de_sign_emitter_sha256` only when the builder's emitter pin
  differs from the incumbent's. Tested under /bin/bash 3.2 (edb5b8dc vs edb5b8dc: nothing added; 1658ae2c vs edb5b8dc: the
  key is added). 497's day script is now **939e32ec**, and chain497e waits for it.
- **F3 (low) — accepted** (the closing meta line prints `bulk_sharpen_pt`).
- A second defect found while applying these: `commit_after_chain.sh` put its file paths in argv. The chains' "no other
  build" regexes match any bash process whose args contain `run_arm_generic.sh` or `def_day.sh`, so the 476 and 497
  committers would have deadlocked chain495 → chain476. The committer now reads a spec file, and all five were relaunched.
