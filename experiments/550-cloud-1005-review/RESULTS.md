# 550 — review of the 10-05 October apply

Static review plus temp-copy and synthetic runs. The live nine files were not modified. `vcc` was not run. No challenge matrices, donor tables, or the callblock artifact were on this VM.

The first pass of this review did not have the snapshot directory. Follow-up 1 checked out `experiments/550-cloud-1005-review/` from `b185ac12` and nothing else from master. 537's v4 text is `snapshots/537-RESULTS-v4-20260930T2031Z.md`. The conflicting sentence in the task ("apply 530's bundle; re-test `patch_defday528.py` on the bundled def_day.sh") is that snapshot's §5.3, lines 666–667, and it is the sequence that was run. Swap and steer for `04fb24ec` are in `snapshots/446bj2zcp.h5ad.meta.json` and match the builder defaults 549 uses.

Follow-up 1 uses the snapshot added in `b185ac12` (`git checkout origin/master -- experiments/550-cloud-1005-review/` only). That snapshot's §5.3 (`snapshots/537-RESULTS-v4-20260930T2031Z.md:666-667`) is the "apply 530, then re-test `patch_defday528.py`" sequence already shown to exit 1 and leave the block out.

## What was run

| run | result |
|---|---|
| `uv run python experiments/538-cloud-bundle-extras/test_extras538.py` | 13 tests, OK, 29.954s. Live nine files still at the pre-test hashes (re-checked at the end of `sequences550.py`). |
| `apply_october538.py --root <temp> --check` | rc 0, `patches apply cleanly (not written)`. Temp root byte-identical. `def_day.sh` before `0b53ba89…ba459c`, after `cb82c34f…d01f5c` (same pair 538 recorded). |
| `repro/sequences550.py` | rc 0. Log: `repro/last-sequences.txt`. |
| `repro/test_bulk550.py` | identical h5ads PASS; one gene +7 in group B/T2 FAIL; `.vcc` (tar + `pred.h5ad.zst`) PASS. |
| `repro/busy_sigstop550.sh` | PASS. |
| `repro/patch_census_preseed550.py --self-test` | PASS. Not applied to any live file. |
| `repro/procs_clear550.sh` | `CLEAR`, exit 0. |
| `repro/test_followup550.py` | rc 0. 549's EXPECT set equals the draft's three keys. A drifted swap passes against 446bj2z and fails against the cp snapshot. |
| `repro/patch_hook_conditional550.py --self-test` | rc 0. Not applied to the live recipe. |

`git log --since='2026-09-30 11:40:00 +0000'` on the nine files and on the callees they import or exec (`def_pin_panel_artifacts.py`, `def_lodo_key.py`, `def_census.py`, `make_submission_pm_cd4.py`, `rc_emit_k.py`, `callblock528.py`, `feat528.py`, `def_preflight.py`, `build_submission.py`) is empty. `git diff abcb541d^..HEAD` on that set is empty. No anchor drift.

## Findings

### HIGH — do not re-run `chain502.sh`, `fix526.sh`, or `chain476.sh` as the known-answer driver

`experiments/502-october-de-riders/chain502.sh:22-28` runs **before** the wait. If `work/pre502.def_day.sh` and `work/pre502.def_day_helpers.py` exist and the sha pair is not `e29e793a…` / `60b2a777…` or `610db782…` / `20142142…`, it copies those backups over the live scripts. The wait then requires def_day prefix `e29e793a` (`:33`). On 10-05 the cap `202610031200` is already past, so the next line exits 1. The revert has already happened.

Repro (temp tree, clock forced to `202610050400`): `repro/sequences550.py` section "chain502 clobber". Output: `reverted …` then `cap hit, exit 1`, and the day file's contents were `PRE502-DAY`. The backups are not in this checkout (gitignored); they are on the Mac if chain502 finished there.

`experiments/526-cloud-october-review2/fix526.sh:21-22` caps at `202610021200`. On 10-05 the first failed prerequisite check exits 1 and does not run the KA. It does not pin `def_day.sh`. `experiments/476-def-day-merge/chain476.sh:34-36` restores five pre-476 files on any failure.

The KA **command** those scripts used is still the right command. The scripts are not.

### HIGH — nothing local checks that a DE-only arm left the bulk identical to 446bj2z

`experiments/428-zero-gene-audit/run_arm_generic.sh:90-97` compares the bulk to **410b** inside a band (0.85–1.15 when `VCC_BULK_CTX_SCALE` is set, which 549 sets; else 0.80–1.03). A 10% bulk move still prints `CLIFF OK`. `scripts/def_postbuild.py` has no count-sum digest. The meta's `output.sha256` is the whole h5ad; `effects.final_sha256` includes the DE. There is no pre-slot gate for "pds and mse equal 446bj2z to ≤ 1e-6".

`repro/bulk_identity550.py` streams CSR (or a `.vcc`: tar member `pred.h5ad.zst`, decompressed with the `zstd` binary, same as `scripts/restore_prediction_h5ad.py:91-108`) and compares per-`(context, target_gene)` summed counts. It does not densify the matrix. Synthetic test: `repro/test_bulk550.py`.

Run it after `experiments/428-zero-gene-audit/work/build_446bj2zcph/cph549.h5ad` exists and **before the slot is armed**, against the 446bj2z prediction h5ad or its `.vcc`. `build549.sh:62` leaves the h5ad on disk after prep, so the check can run when `build549.sh` exits 0. Do not use it to compare an October day build to 446bj2z: October refits the operator and re-levels; the bulk is not 446bj2z's pinned dump `00ccfcba…` or scales `A=1.0962,B=1.1461,C=1.1681` (`build549.sh:47-53`).

### HIGH — `build549.sh` can leave 537's driver stopped, and the wait has no cap

`experiments/549-cph-arm/build549.sh:36-43` sends SIGSTOP to the pid in `experiments/537-arm-pipeline/work/drive537.pid` when its args contain `drive537`, then waits with no `202610061200` cap (the other loops at `:29` and `:46` have one). The EXIT trap (`:18-20`) sends SIGCONT.

`repro/busy_sigstop550.sh`: after SIGTERM the child went from `T` to `S` and the log contained `RESUMED`. After SIGKILL of the pauser the child stayed `T` until this review sent SIGCONT. A reboot clears the process table (the driver is dead, not left in `T`). A sleep does not: a stopped driver stays stopped. `drive537.sh` is not in the repo, so the comment that the job is a direct child named `time|uv|python` was not checked against the driver. `pgrep -P` only sees direct children; a job that is not a direct child is invisible to the wait, and a direct `uv`/`python` that never exits holds the stop past 10-06.

### HIGH — 530 then 528 does not install the call-propensity block. The one applier does.

`experiments/528-callprop-block/patch_defday528.py:70` anchors on the MDE `grep` immediately followed by `fi`. 530 moved that grep. `:84-89` writes only after every anchor matches once, and `:23-24` exits `already patched` if `DEF_DAY_CALLSET_BLOCK` is already present.

On a temp copy, `patch_bundle530.py` then `patch_defday528.py`: 528 rc 1, `anchor found 0 times`, def_day bytes unchanged, `DEF_DAY_CALLSET_BLOCK` absent. `apply_october538.py` (530, then that miss, then 538): rc 0, the block is in the file, `bash -n` clean, a second apply exits `already patched`, and 528-after-538 exits `already patched` without writing. Log: `repro/last-sequences.txt`.

The committed 537 RESULTS does not contain the "re-test `patch_defday528.py`" sentence. That sentence, run as written, stops with exit 1 and leaves the block out. 538's RESULTS ("Commands the coordinator runs on 10-05") is the sequence that installs it. Do not also run `patch_bundle530.py`.

### MEDIUM — `DEF_DAY_CENSUS_MAX_ZERO` dies on the directory `recipe_dr524.sh` pre-seeds

`experiments/524-october-recipe-dr/recipe_dr524.sh:44-45` copies the whitener npz and json into the work dir, then execs `def_day.sh`. The patched census check (inserted by `patch_extras538.py:172-173`, after `mkdir -p "$W"`) dies when `ls -A` is non-empty. Unset, the block is skipped.

Repro in `sequences550.py`: pre-seeded dir with the knob set exits 1 `not empty`; the same dir with the knob unset prints `CENSUS_OK`; an empty dir with the knob set prints `CENSUS_OK`.

The knob is not a resolver switch, so `decide` cannot turn it on, and the 02:00 rehearsal does not set it. It kills the day only if someone exports it on a pre-seeded dir (the rehearsal, or an October work dir that already holds the whitener). `repro/patch_census_preseed550.py` is the unapplied fix: allow only `whitener_union_pruned_def.npz` and `.json`; any other file still dies. `--self-test` passed. It refuses to edit a path inside this repo.

### MEDIUM — `build549.sh` `busy()` is the pattern `recipe_dr524.sh` already replaced, and the two scripts do not see each other

`build549.sh:25` matches any bash line that **contains** `run_arm_generic.sh` or `def_day.sh`. `recipe_dr524.sh:18-22` matches an argv word that **ends** in the script name. `repro/busy_sigstop550.sh`: a canned `build549.sh` line does not match; a bash whose argv mentions `run_arm_generic.sh` does; a live `bash -c 'sleep 30; true' run_arm_generic.sh` made the same awk return busy. `other_builders()` caught a `uv`/`python` sleeper whose argv contained `make_submission_pm_cd4.py`, and a failed python check is fail-closed (`[ "" != 0 ]`). `build549.sh` does not match itself.

`recipe_dr524.sh`'s `builders()` does not name `build549.sh`. A 549 wait loop is invisible to the rehearsal; once 549's child is `run_arm_generic.sh`, the rehearsal does see it. A false busy parks 549 until `202610061200`.

### MEDIUM — the 02:00 rehearsal cannot be a test of the call-propensity block, and `command` without `--allow-pending` still exits 2 after the two decides

`recipe_dr524.sh` resolves `october_recipe.py command --allow-pending` on the **unpatched** tree (`:33-34`), then runs `def_day.sh` (`:47`), then `gate524.py --write` (`:62`). `callset_block` is not a switch until the applier. Unpatched `decide callset_block` on a temp decisions file: rc 1 `unknown switch`, file unchanged (`sequences550.py`).

`chain1006.sh:17-20` waits until 20:00Z and until `recipe_dr524.sh`'s pid is dead. It greps `gates524.log` for `^[tag] PASS`, `^[tag] CLIFF OK`, `novelty rc=0`, and `package sha256:`. It does not execute the nine files. Do not re-run `gate524.py --write` after the patch: `:170-174` overwrites the log only when every gate passes. `arm1005.sh` / `arm1001.sh` / `chain1005.sh` / `chain1001.sh` do not execute the nine files.

A morning `decide` of an **existing** switch survives the JSON edit. Temp decisions with `mde_sign` already `DECIDED clean4:0.05` (and a sentinel `de_mag` evidence string), then the applier: status and value and the sentinel were intact; `lodo_hi` and `callset_block` landed `PENDING` immediately before `cd4_thresholds` (`patch_bundle530.py:662-707`, `patch_extras538.py:453-485`).

After `decide callset_block=prop` and `decide mde_sign=clean4:0.05` on the **temp** tree, `command` without `--allow-pending` exited 2 (`PENDING: lodo_hi_pct, sharpen_pt, rho_scale, rc_callset_cpm, hook_estimator, lodo_strength, lodo_hi`). With `--allow-pending`, rc 0 and the env was `DEF_DAY_LODO=1 DEF_DAY_DE_MAG=dn_k250_zs DEF_DAY_MDE_SIGN=clean4:0.05 DEF_DAY_CALLSET_BLOCK=prop`. The live `docs/october-decisions.json` hash was unchanged.

`eq_head` in `recipe_dr524.sh:27` compares `def_day.sh` to HEAD **during the wait**. Patching while that process is still in the until-loop makes the rehearsal wait until 20:00Z and exit 1 (disk ≠ HEAD). Committing the patch during the wait makes the rehearsal run the patched script. The process that must have exited is `recipe_dr524.sh` itself, not merely "no `def_day.sh` child".

### LOW — a context whose control libraries fall in one quantile bin makes `callblock528.py` die

`experiments/528-callprop-block/feat528.py:67-70`: `np.quantile` on a constant library, then `L[b == k].mean()` for empty bins. Synthetic constant library of 20 cells: 9 of 10 bin means were NaN (`sequences550.py`). `callblock528.py:98-99` exits if the key is non-finite, and the patched day treats that as `die`. Real A/B/C libraries span thousands of counts; this is a degenerate-control failure, not the A/B/C artifact path.

### LOW — a kill during install is not one transaction

`apply_october538.py:144-153` `os.replace`s each of the nine files in turn. A kill in the middle trips `partial_patched` (`:81-94`). Repro: copy only the patched `def_day.sh` onto an otherwise unpatched tree; `--check` exits 1 `partial extras (have def_day LODO_HI; missing …); restore and re-apply`. Restore the nine from git and re-apply. Do not match the bare prefix `DEF_DAY_LODO_HI`.

## D — can the patched day express 549's DE arm?

549's builder env (`build549.sh:53-55`) is `VCC_RC_CALLSET_BLOCK=prop:<abc npz>:04fb24ec…` with no swap/steer field, and `VCC_DE_MDE_SIGN=clean4:0.05`. `build535.sh:80` is the same block string. The builder (`scripts/make_submission_pm_cd4.py:2021-2030`) reads a missing swap as 1.0 and a missing steer as `kept`, and `block_keep = 2200 - round(1.0 * 2200) = 0`. It passes that steer into the emitter at `:6021`. The emitter's own default `block_steer="key"` (`src/vcc2026/rc_emit_k.py:806`) is not what the builder uses. The patched FINALCHECK uses the same defaults (swap 1.0, steer `kept`).

On the patched day, `DEF_DAY_CALLSET_BLOCK=prop` (no suffix) builds `prop:$W/callblock_def.npz:<sha>` and does **not** append `:swap:steer` (`patch_extras538.py:54-63`, confirmed present in the patched temp file). `DEF_DAY_MDE_SIGN=clean4:0.05` is exported as `VCC_DE_MDE_SIGN="$MDS"` (`patch_bundle530.py:220-221`); the `hookbank:` rewrite does not touch `clean4:`.

`def_day.sh:66` wipes every `VCC_*`, then re-exports only what it sets. `def_pipeline.py` `resolve_env` (`:164-179`) fills `SHIPPED_ENV` keys that are absent. The build stage (`:717-718`) is not `bare_env`. Child env is `{**os.environ, **resolved}` (`:511-513`). Synthetic parent with `VCC_DE_MDE` absent, `VCC_DE_MDE_SIGN=clean4:0.05`, and a `prop:` block: the child had `VCC_DE_MDE=0.1,0.1,0.3,coarse`, `VCC_DE_ZERO_REPAIR=-5`, `VCC_DE_SIGN_PRIOR=clean4q80`, and both parent exports. That is what `clean4` needs (`make_submission_pm_cd4.py:822-825`). `run_arm_generic.sh:34` sets the same sign prior for 549. The pin stages use `bare_env=True` (`:692`, `:697`); they do not need these exports.

`callblock528.py` `main` is the function both paths call. The day argv is `--bundle "$B" --contexts "$CTXS" --out "$W/callblock_def.npz"`. The artifact command (`callblock528.py:15-16`) is the same flags on `data/challenge` and `A,B,C`. Monkeypatching `features` and calling `main` twice on one synthetic bundle wrote identical npz bytes (sha `4ab61ebd…`, len 3388), and `key_A` matched `key_from_features` at `BETA=22`, `MEDIAN_X_H1=0.00877448551593013`.

Not the same arm on the bulk. 549 pins `VCC_BULK_LOAD` to `dump_h1x2_cd4x15_lodo05.npz` and fixed context scales. October builds the dump that day and re-levels. Also not the same if `rc_callset_cpm` is decided to something other than 2200: at swap 1.0 `block_keep` stays 0, but the block length changes. 549 does not set `VCC_RC_CALLSET_CPM`.

## E — 549 vs the build gates that do exist

`run_arm_generic.sh:81-82` compares `EXPECT` as a set (`set(d) == want`). Order does not matter. Repro printed `set equal True` and `string equal False` for `rc_callset_block,de_mde_sign` vs the reverse. `build549.sh:51` builds that string and appends `de_sign_emitter_sha256` only when the emitter sha differs from the reference meta. Artifact sha is checked at `:34` against `04fb24ec…`. Builder/emitter prefixes `4daae237c5fe4fd2` / `0f56209d57e8385a` are checked at `:33`. `busy()` does not match `build549.sh`'s own argv; it does match a bash mention of `run_arm_generic.sh` or `def_day.sh`, and it does not miss a python builder (`other_builders`). SIGSTOP is the finding above.

## G — 549 against 537's CPH draft

Draft: `snapshots/537-PREREGISTRATION-draft-CPH-20260930T1635Z.md:6-17`. Build path: `experiments/549-cph-arm/build549.sh:47-56`. The env matches the draft and §5.2 item 1 of the v4 snapshot (`:597-607`): same bulk load, MDE, zero repair, context scales, `VCC_DE_MAG`, the two unsets, then `VCC_RC_CALLSET_BLOCK=prop:…:04fb24ec…` and `VCC_DE_MDE_SIGN=clean4:0.05`. The cp snapshot's block is `source prop`, `swap 1.0`, `steer kept`, `block_keep 0`, sha `04fb24ec…` (`snapshots/446bj2zcp.h5ad.meta.json:94-118`). 549 does not pass swap or steer, and the builder defaults are those values.

| falsifier | enforced before a slot? |
|---|---|
| Config diff vs 446bj2z is exactly `rc_callset_block`, `de_mde_sign`, `de_sign_emitter_sha256` | Yes, as a **key set**. `build549.sh:50-51` adds `de_sign_emitter_sha256` because the emitter file is `0f56209d…` and 446bj2z's meta records `edb5b8dc…`. `run_arm_generic.sh:76-86` stops before prep unless `set(diff) == set(EXPECT)`. Order does not matter. The draft's string and 549's string are the same set (`test_followup550.py`). |
| The block's swap, steer, and sha equal 446bj2zcp's | No. The gate compares keys, not values. A swap of 0.5 still differs from 446bj2z's absent key, so the same EXPECT passes (`test_followup550.py`, drifted vs z, rc 0). |
| Config diff vs 446bj2zcp is exactly `de_mde_sign` | No. `run_arm_generic.sh` takes one reference meta, and 549 passes 446bj2z's. The cp snapshot has no `de_mde_sign`. Against 446bj2z it differs in exactly `de_sign_emitter_sha256` and `rc_callset_block`. A built meta that is that config plus `de_mde_sign` differs from the cp snapshot in exactly `de_mde_sign`. A drifted swap differs in `de_mde_sign` and `rc_callset_block` and would be caught only by that second comparison. `repro/config_keys550.py` is that comparison. Run it on the built meta against `snapshots/446bj2zcp.h5ad.meta.json` with expected key `de_mde_sign`, after the h5ad meta exists and before the slot. |
| pds and mse equal 446bj2z's (and 446bj2zcp's) to ≤ 1e-6 raw | Not locally. Nothing in the build runs `cell-eval2`. The board read is what reports those members. The bulk-leak half of that falsifier is `repro/bulk_identity550.py` (item E): same per-(context, target) count sums means the bulk channel those two members read did not move. `|Δjac|` and `Δfid` in the draft's later falsifiers are board-only too. |

## H — hook conditional on the block, CLEAN4 as the default source

v4 §5.4 (`snapshots/537-RESULTS-v4-20260930T2031Z.md:680-688`): `DEF_DAY_MDE_SIGN=clean4:0.05` should be on only with `DEF_DAY_CALLSET_BLOCK=prop`, unless 496a and a D/E/F entry both favour the hook alone; CLEAN4 stays the default source, and Arc moves only on an explicit later read.

The patched resolver can store both switches as `DECIDED` (`mde_sign=clean4:0.05`, `callset_block=prop`) and emit both. It cannot express the conditional. `mde_sign` and `callset_block` are independent. `_problems_530` accepts `clean4:`, `programme:`, and `hookbank:` with no mention of the block (`patch_bundle530.py`, the check that becomes the patched recipe's lines that start `if mde and not (mde.startswith("clean4:")`). A temp recipe after 538, with `mde_sign` decided `clean4:0.05` and `callset_block` still `PENDING`, `check --allow-pending` exits 0 and exports the hook alone.

`repro/patch_hook_conditional550.py` is not applied. On a copy of the 538-patched `october_recipe.py` it refuses:

- the hook without `callset_block=prop`, unless the switch evidence contains `hook-alone-explicit`
- a source other than `clean4` (`programme`, `hookbank`, `file`), unless the evidence contains `arc-source-explicit`

`clean4:0.05` with `callset_block=prop` passes with neither token. `decide` does not write either token. Self-test: unpatched hook-alone rc 0; patched hook-alone rc 3; clean4+prop rc 0; `hookbank:` without the token rc 3; the same value with `arc-source-explicit` rc 0; hook-alone with `hook-alone-explicit` rc 0. The live nine files were unchanged. It refuses to edit a path inside this repo.

## What this VM could not run

- The inert KA against `ac49d319…`. No `ef-bundle`. The patched script still contains `DONE. Next` (`scripts/def_day.sh:361` and `:363` are unchanged; the patch inserts the callset env line in front of `log "final:"`, not in front of `DONE`).
- 528's D3 byte-compare to `04fb24ec…`. The npz is not here. Equivalence above is the same function and the same CLI shape on synthetic features, not the A/B/C bytes.
- A real `callblock528.py` on challenge or D/E/F h5ads. `cell-eval2`. `vcc`. The Mac path `~/Code/virtual-cell-challenge-2026` (the patched script derives `REPO` from `dirname`).
- The Mac's scored `446bj2zcp` package, beyond the snapshot meta in this directory. Swap, steer, `block_keep`, and the artifact sha were read from that snapshot (section G).
- `drive537.sh`. Whether a re-run of the patched gate would pass a real `524dr` meta. 453dr's `bulk_load` is `dump_h1x2_cd4x1.npz`, which is what the new else branch (`patch_bundle530.py:794-795`) requires when stack and strength are off.

## 10-05 command sequence

1. Wait until `recipe_dr524.sh`'s **process** has exited (the wait loop included, not "it is 04:00"). Then:

```bash
bash experiments/550-cloud-1005-review/repro/procs_clear550.sh
```

Exit 0 and the line `CLEAR` means no bash/python/uv argv word is `recipe_dr524.sh` or one of the nine files. Exit 2 means do not patch. `chain1006.sh` may already be sleeping until 20:00Z; it does not read those files. Do not use `pgrep -f` on a bare name.

2. Do not run `chain502.sh`, `fix526.sh`, `chain476.sh`, or `patch_bundle530.py`, or `patch_defday528.py` on its own. Do not re-run `gate524.py --write`.

3. Apply, on the repo root:

```bash
uv run python experiments/538-cloud-bundle-extras/apply_october538.py --root . --check
uv run python experiments/538-cloud-bundle-extras/apply_october538.py --root .
uv run python experiments/538-cloud-bundle-extras/test_extras538.py
bash -n scripts/def_day.sh
uv run python scripts/october_recipe.py check --allow-pending
```

`--check` must print `patches apply cleanly (not written)` and `def_day.sh` after-hash `cb82c34f808d49700f2a8e369978811ad30eb07e7d79d27ceed7f1be76d01f5c`. A second apply must exit `already patched`.

4. Inert known answer. Fresh work dir. This is **not** `october_recipe.py command` (that already exports `DEF_DAY_LODO=1` and `DEF_DAY_DE_MAG=dn_k250_zs`). From the repo root, so the derived `REPO` is this checkout:

```bash
rm -rf experiments/502-october-de-riders/work/ka1005-inert
env -u DEF_DAY_LODO -u DEF_DAY_LODO_HI -u DEF_DAY_LODO_HI_PCT -u DEF_DAY_LODO_STRENGTH \
    -u DEF_DAY_DE_MAG -u DEF_DAY_MDE_SIGN -u DEF_DAY_CALLSET_BLOCK -u DEF_DAY_CENSUS_MAX_ZERO \
    -u DEF_DAY_RHO_SCALE -u DEF_DAY_RC_CALLSET_CPM -u DEF_DAY_CD4_STACK -u DEF_DAY_HOOK_ESTIMATOR \
    -u DEF_DAY_PROGRAMME_SOURCE -u DEF_DAY_PRIOR_GATE -u DEF_DAY_SHARPEN_PT -u DEF_DAY_FORCE \
    -u DEF_DAY_KEEP_PINS -u DEF_DAY_ALLOW_NO_NTC_ID \
    DEF_DAY_CD4_THRESHOLDS=relaxed \
    bash scripts/def_day.sh experiments/238-essential-free-panel/work/ef-bundle \
    experiments/502-october-de-riders/work/ka1005-inert 474efp 1.75 auto pt:noprog 4
```

That directory is only a fresh work dir. It does not invoke `chain502.sh`. Require `DONE. Next` in `ka1005-inert/def_day.log` and `output.sha256` `ac49d319f5982507003b1f2070660d2d03ee351fd13c9fdd30df8f5fb08ec0b3` in `474efp_final.h5ad.meta.json`. `def_day.sh` wipes `VCC_*` itself; it does not wipe `DEF_DAY_*`, which is why the `env -u` list is the whole knob set.

5. Commit through `scripts/commit_after_chain.sh <spec>`. The spec's argv must not contain `def_day.sh` (the FILE lines live in the spec). `EXPECT` can pin `scripts/def_day.sh=cb82c34f808d4970`. `REQ` must match a line only the KA log writes, including the sha. Refuse the commit if the log contains `RESTORED`.

6. Only then `decide` the new switches, on the patched resolver. `callset_block=prop` (swap 1.0, steer `kept`, the 549/535 default). `mde_sign=clean4:0.05` only if that read is actually in. A decide of `mde_sign` made **before** 02:00 on the unpatched file is preserved by the applier; `callset_block` cannot be decided before the applier. The applier does not couple them: a decided hook with the block still pending is accepted (section H). To get the refusal, after the applier and before `decide mde_sign`, copy `scripts/october_recipe.py` out of the repo and run `repro/patch_hook_conditional550.py --apply-to` on that copy, then review the diff. The script refuses to edit the live file.

7. Rider-on rehearsal: `october_recipe.py command --allow-pending` (without the flag it exits 2 while `lodo_hi` and the 530 switches are still PENDING). That command will not reproduce `446bj2zcph`. Its DE exports for this form are `DEF_DAY_CALLSET_BLOCK=prop` and `DEF_DAY_MDE_SIGN=clean4:0.05`, on top of the already-decided `DEF_DAY_LODO=1` and `DEF_DAY_DE_MAG=dn_k250_zs`. The 10-06 read of `524dr` is the pre-patch recipe. It is not a test of the block.

## Local verification on the Mac (session "Explore contribution areas", 2026-09-30 ~22:20Z)
- **chain502 clobber: CONFIRMED.** `work/pre502.def_day.sh` / `pre502.def_day_helpers.py` exist (Sep 28). Live pair is
  `0b53ba89` / `20142142`, which is neither allowed pair at `chain502.sh:22-28`, so a re-run copies the pre-502 files over
  the live scripts before its expired cap exits. Never re-run chain502/fix526/chain476 as the KA driver.
- **Inert KA command: CONFIRMED** equal to `chain502.sh:52` (`DEF_DAY_CD4_THRESHOLDS=relaxed bash scripts/def_day.sh <efp>
  <fresh dir> 474efp 1.75 auto pt:noprog 4`); the `env -u` list is needed because def_day wipes `VCC_*` but not `DEF_DAY_*`.
- **Loose bulk gate: CONFIRMED** (`run_arm_generic.sh:90-97`: 0.85–1.15 vs 410b when `VCC_BULK_CTX_SCALE` is set).
- **recipe_dr524's wait compares def_day/helpers/builder/emitter to HEAD every loop: CONFIRMED** (`recipe_dr524.sh:27`).
- `procs_clear550.sh` now: BUSY (pid 99304 `recipe_dr524.sh`, rc 2), as it should be before 10-05.
- `config_keys550.py` on the real metas: 446bj2zcp vs 446bj2z → exactly `de_sign_emitter_sha256`, `rc_callset_block`: PASS.
- **Caveat on `bulk_identity550.py`:** its `.vcc` path reads the whole package and zstd's output into memory (≈ 12+ GB for a
  real 2.77 GB package). Use it on two `.h5ad` files only (restore a `.vcc` with `scripts/restore_prediction_h5ad.py`
  first); its h5ad path streams row by row (slow, ~0.2 GB RAM).
- **First real use (coordinator, 2026-09-30 ~22:30Z, 549's `cph549.h5ad` vs 446bj2zcp's `callprop535.h5ad`):** 895/900
  (context, target) groups identical; in 5 groups exactly one count moved between two low-count genes (+1/−1, group totals
  equal; e.g. A MAPK12: VPS9D1 30→29, PLD2 50→51): the CLEAN4 hook's whole-count moves landing on a neighbouring gene,
  pds/mse impact far below 1e-6 raw (549 RESULTS, 79d01405). **Read the check's per-gene diffs, not only its exit code:**
  a +1/−1 swap between genes with equal group totals is not a bulk leak; a leak moves group totals or many genes.
