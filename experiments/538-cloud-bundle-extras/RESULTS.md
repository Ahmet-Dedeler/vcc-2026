# 538 — October bundle extras (NOT applied)

Scorer-mechanics / October recipe plumbing. Default-off. No score moved, and no prediction changes until a switch is `DECIDED`. The live tree is untouched: `scripts/def_day.sh` is still `0b53ba89bd4774896d38b9f75c37145f060bd807b34073096f8d1fea09ba459c`.

Line numbers below are the **patched copies** after `apply_october538.py` (530's bundle, then 538's extras). 528's own patch does not write on that tree; 538 re-inserts its snippets. Do not read these lines against the live files.

## What each knob does

1. **`DEF_DAY_LODO_HI=<h>`** (531 #12). Needs `DEF_DAY_LODO=1`. `1 < h ≤ 20`, plain decimal (`[1-9][0-9]*(?:\.[0-9]+)?`); refuses `1`, `21`, `1e1`, `06`, `0.5`, and any value set together with `DEF_DAY_LODO_HI_PCT`. Unset keeps today's cap 3 or the percentile form.
   - Header: `scripts/def_day.sh:22-24`. Validation: `:119-131` (needs LODO at 121, mutual exclusion at 122, regex at 123-129, `grep -q 'add_argument("--hi"'` at 130 — `--hi-pct` contains `--hi`).
   - Key name and argv: `:318-325`. Fixed cap writes `pert_scale_lodo_def_lo0.5_hi${h}.npz` / id `lodo_lo0.5_hi${h}_def` and passes `--hi "$DEF_DAY_LODO_HI"`. Percentile stays the `elif`. Unset passes neither `--hi` nor `--hi-pct` (script default cap 3, filename `…_hi3.npz`).
   - `def_day.final.env`: printf at `:513`.
   - `riders_match` fourth argument gains `,lodo_hi=<h>` only when the knob is set (`:476`, `:481`). Optional key in `scripts/def_day_helpers.py:94`; a recorded `config.lodo_hi` must match (`:134-142`); absence passes, because the level meta does not carry the cap. The pin's json does (`def_lodo_key.py` already writes `hi`).
   - Resolver switch `lodo_hi` PENDING, `docs/october-decisions.json:177-185`, before `cd4_thresholds`. Mapping `scripts/october_recipe.py:45`; range/deps `_problems_538` at `:168-183`; export marker `'--hi "$DEF_DAY_LODO_HI"'` at `:201`. PENDING exports nothing.
   - gate524: `SWITCH_KEYS["DEF_DAY_LODO_HI"] = {bulk_load, bulk_load_sha256}` at `:51`. Value row "LODO fixed cap" at `:210-223`: sidecar `pert_scale_lodo_def_lo0.5_hi<h>.json` field `hi` and log id `lodo_lo0.5_hi<h>_def`.

2. **Call-propensity block** (528, re-anchored). `DEF_DAY_CALLSET_BLOCK=prop[:<swap>[:<steer>]]`, `0 < swap ≤ 1`, steer `kept|key|cpm`. Unset: the `if` bodies are skipped; `CSB` stays empty; no export.
   - VALIDATE `scripts/def_day.sh:191-202`, immediately before `SQ="$Q"` at `:203`. COMPUTE `:302-310` (`callblock528.py` at 305). EXPORT `:378` `export VCC_RC_CALLSET_BLOCK="$CSB"`. FINALCHECK `:490-500` (source `prop`, sha256, swap default 1.0, steer default `kept`). Allow-list key `rc_callset_block` at `:527`. Env record printf at `:515`.
   - Those four snippets are the constants from `experiments/528-callprop-block/patch_defday528.py`, byte-identical (`test_528_anchor_miss_leaves_bytes_and_snippets_match`).
   - Switch `callset_block` PENDING, `docs/october-decisions.json:186-194`. Mapping `october_recipe.py:46`. Builder need `VCC_RC_CALLSET_BLOCK`. Emitter/callblock paths via `OCTOBER_EMITTER` / `OCTOBER_CALLBLOCK` (`:54-55`); installed check `:214-219`.
   - gate524 `SWITCH_KEYS["DEF_DAY_CALLSET_BLOCK"] = {rc_callset_block}` at `:52`. Value row "call-propensity block" at `:224-246` (`config.rc_callset_block`: source, 64-hex sha256, swap, steer).
   - `riders_match` optional key `callset` (`helpers.py:75-89`, `:143-158`). This one **is** in the level meta and is checked. A 4-key 530 spec (no `lodo_hi`, no `callset`) still parses. Unknown keys or a missing required key return None (exit 2). When both extras are unset and a 530 knob is on, the fourth argv word is exactly `rc=${_rc},stack=${_st},rho=${_rh},strength=${_ls}`.

3. **`def_postbuild.py --require-config-diff <k1,k2,…>`** (531 #17). Argparse `:212-216`. Check `:475-483`, after the existing unexpected-diff `chk` (`:473`). Fails unless every listed key is in `moved` **and** nothing outside the allow-list moved. Required keys are **not** auto-allowed; `--allow-config-diff` keeps its meaning. Empty (the default) adds no `chk` and does not change today's text or exit code.

4. **`DEF_DAY_CENSUS_MAX_ZERO=<n>`** (531 #16). Integer `6 ≤ n ≤ 30`, no leading zero (`06`, `6.0`, `5`, `31` refused). Header `:25`. Range and fresh-dir check `:86-94`, after `die()` (`:85`) and before `log "start"` (`:97`), so `def_day.log` does not exist yet. Reused means the work dir exists and `ls -A` is non-empty. Setup forwards `${DEF_DAY_CENSUS_MAX_ZERO:+--census-max-zero "$DEF_DAY_CENSUS_MAX_ZERO"}` at `:235`. After setup, reads `census.json` → `consensus.targets_with_0_donors` and logs `census override: max zero-carrier <n> (count <c>)` at `:242-245`. Recorded in `def_day.final.env` at `:514`.
   - `scripts/def_pipeline.py` argparse `--census-max-zero` at `:546-548`. Census stage `:649-659` appends `--max-zero-carrier` only when the flag is passed, and re-checks the range. `def_census.py` is not edited; its existing `--max-zero-carrier` (default floor 5) is what receives the override. Unset: the census argv has no extra word.

5. **REPO from the script** (531 #20). `scripts/def_day.sh:74-78`:

   ```
   REPO="$(CDPATH= cd "$(dirname "$0")/.." && pwd)"
   [ -f "$REPO/scripts/def_day.sh" ] || { echo "STOP: derived REPO=$REPO has no scripts/def_day.sh" >&2; exit 1; }
   cd "$REPO" || exit 1
   ```

   `CDPATH=` keeps `cd` from printing. On the production Mac, `dirname` of `…/virtual-cell-challenge-2026/scripts/def_day.sh` is that checkout. In a worktree it is the worktree. A copy that does not live in `scripts/` is refused. This VM did not execute the Mac path.

## Applier

`apply_october538.py [--check] [--root <repo>]` copies the nine files to a stage, runs `experiments/530-cloud-october-bundle/patch_bundle530.py`, then `experiments/528-callprop-block/patch_defday528.py`, then `patch_extras538.py`. It installs with a sibling temp and `os.replace` only if the whole sequence succeeded. `--check` does not copy back and prints `patches apply cleanly (not written)`. A second apply (markers already present) exits `already patched` before writing. An anchor miss in 530 or in 538 is non-zero and leaves the destination untouched. 528's anchor miss after 530 is expected and is printed, not a failure.

Files (9): `scripts/def_day.sh`, `scripts/def_day_helpers.py`, `experiments/448-def-rehearsal-q175/relevel_iter.py`, `scripts/def_sameline.py`, `scripts/october_recipe.py`, `docs/october-decisions.json`, `experiments/524-october-recipe-dr/gate524.py`, `scripts/def_postbuild.py`, `scripts/def_pipeline.py`.

On a throwaway copy of this tree the applier printed, for `def_day.sh`, before `0b53ba89…ba459c` and after `cb82c34f808d49700f2a8e369978811ad30eb07e7d79d27ceed7f1be76d01f5c`. The live file stayed at the before hash.

## Tests

`uv run python experiments/538-cloud-bundle-extras/test_extras538.py` → **13 tests, OK, 25.224s**. Temp copies only. `tearDownClass` asserts the live nine files are byte-identical to the start of the run.

| test | what it shows |
|---|---|
| `test_check_clean_second_apply_refused_anchor_miss` | `--check` on an unpatched tree exits 0 and writes nothing; a second apply is refused; a corrupted REPO anchor exits non-zero and leaves the dest untouched |
| `test_528_anchor_miss_leaves_bytes_and_snippets_match` | 528 after 530 exits 1 with the anchor found 0 times and does not write; VALIDATE/COMPUTE/EXPORT/FINALCHECK in the 538 copy equal 528's constants, and VALIDATE sits immediately before `SQ="$Q"` |
| `test_bash_n_and_repo_derivation` | `bash -n` clean; `bash -x` shows `REPO=<resolved worktree>`; a script not under `scripts/` is refused |
| `test_unset_matches_530_plus_callset` | knobs unset: 15 uv calls, the LODO=1+DE_MAG day (17), and a rho-only day (`DEF_DAY_RHO_SCALE=1.02`, 38) match the 530+reanchored-callset baseline, including `final.env`. Rho-only riders argv ends `rc=-,stack=x1,rho=1.02,strength=-`. LODO unset passes neither `--hi` nor `--hi-pct`; the file is `…_hi3.npz` |
| `test_lodo_hi_set_and_refusals` | `HI=6` and `HI=20` pass `--hi`, the `hi<h>` filename, the log id, `final.env`, and `lodo_hi=` on riders. Refusals: no `LODO=1`; both cap forms; `1`, `21`, `1e1`, `06`, `0.5`. `HI+PCT+strength` dies on "two cap forms". `HI+strength` without PCT dies on 530's "percentile cap" |
| `test_callset_set_and_refusals` | `prop` and `prop:0.5:key` export `VCC_RC_CALLSET_BLOCK`. `banana`, `prop:2`, `prop:0`, `prop:0.5:nope` die before `callblock528.py` |
| `test_census_override_and_refusals` | `8` logs `census override: max zero-carrier 8 (count 2)` and writes `final.env`. A reused dir dies "fresh work dir". `5`, `31`, `6.0`, `06`, `0`, `abc` die with no uv records. `30` succeeds |
| `test_pipeline_forwards_max_zero_carrier` | patched census dry-run with the flag unset matches the unpatched stdout and has no `--max-zero-carrier`. `--census-max-zero 8` adds `--max-zero-carrier 8`. `5`, `31`, `06` are refused |
| `test_postbuild_require_config_diff` | synthetic new-panel world. Flag absent: patched stdout/stderr/rc identical to unpatched, rc 0. `k1,k2` both moved and allowed → pass. Only `k1` moved → fail `missing=['k2']`. An extra `k3` outside the allow-list → fail mentioning `k3` |
| `test_recipe_pending_exports_nothing_new` | `command --allow-pending` is exactly `DEF_DAY_LODO=1 DEF_DAY_DE_MAG=dn_k250_zs bash scripts/def_day.sh B W N 1.75 0.5 noprog 4`. Both new switches are PENDING and sit before `cd4_thresholds`. `cd4_stack` stays DECIDED `x1` |
| `test_recipe_decided_and_refusals` | decided `lodo_hi=6` and `callset_block=prop:0.5:key` export; bad ranges, missing `LODO=1`, both cap forms, and bad grammar exit 3 |
| `test_riders_match_extension` | a 4-key 530 spec passes; absent `lodo_hi` passes; recorded `6` matches; callset `prop` and `prop:0.5:key` match; bad hi / missing block / bad swap / bad steer fail; an unknown key and a short spec exit 2 |
| `test_gate_rows` | empty env PASS. `LODO_HI=6` with sidecar `hi` 6 PASS; sidecar `hi` 3 FAIL. Callset `prop` with a matching meta PASS. `prop:0.5:key` against swap 1.0 / steer `kept` FAIL |

"Byte-identical when unset" is the day-stub comparison in `test_unset_matches_530_plus_callset`: the baseline is 530 applied plus `callset_edits()` only (528's script cannot apply), with the hard-coded REPO rewritten to the temp repo; the full 538 copy derives REPO. Both are run until the first real data access. A stub `uv` records argv. Recorded paths are normalised. Incoming `VCC_*` and `DEF_DAY_*` are scrubbed before the test env is applied.

## Commands the coordinator runs on 10-05

Only after `recipe_dr524.sh`'s **process** has exited. Do not edit `def_day.sh` while that bash process is running. This applier **includes** 530. Do not also run `patch_bundle530.py` on the same tree.

```bash
uv run python experiments/538-cloud-bundle-extras/apply_october538.py --root . --check
uv run python experiments/538-cloud-bundle-extras/apply_october538.py --root .
uv run python experiments/538-cloud-bundle-extras/test_extras538.py
bash -n scripts/def_day.sh
uv run python scripts/october_recipe.py check --allow-pending
```

Then, still the coordinator's and not run here: the helper suites that already exist (497, 502, 517), the inert known answer (efp relaxed through `def_day.sh` with every new knob unset, byte-for-byte `ac49d319f5982507003b1f2070660d2d03ee351fd13c9fdd30df8f5fb08ec0b3`, as 530's RESULTS requires), and one rider-on rehearsal for whatever `october_recipe.py decide` has turned on. The commit after the inert KA goes through `scripts/commit_after_chain.sh`.

`--allow-pending` stays today's command until someone runs `decide`. Both new switches are PENDING.

## What is wrong or risky in 530 and 528

- **528 does not apply after 530.** Its VALIDATE anchor is the MDE `grep` followed by `fi`. 530 moved that grep, so the anchor occurs 0 times. 528 writes only after every anchor matches, so the file is left byte-identical and the process exits 1. Evidence: `test_528_anchor_miss_leaves_bytes_and_snippets_match`, and the applier's line `528 did not apply (anchor miss after 530); extras re-anchors the call-propensity block`. 538 therefore owns the insertion. If a later edit makes 528 apply cleanly, 538 skips the callset edits when `DEF_DAY_CALLSET_BLOCK` is already in the file.
- **Fixed cap and 530's strength rider cannot both be on.** 530's strength check (`def_day.sh:174` on the patched copy) still requires `DEF_DAY_LODO_HI_PCT`. 538 refuses `LODO_HI` together with `LODO_HI_PCT`, and that check runs first. `HI+PCT+strength` dies on "two cap forms". `HI+strength` without the percentile dies on 530's "needs the percentile cap". This was not relaxed. A strength arm still has to use the percentile cap.
- **`DEF_DAY_LODO_HI` is a prefix of `DEF_DAY_LODO_HI_PCT`.** The first `--check` on an unpatched tree false-positived as "partial extras" because the marker was the bare prefix. Markers are now `'--hi "$DEF_DAY_LODO_HI"'`, `'"lodo_hi": {'`, `'"DEF_DAY_LODO_HI":'`, and `'{"lodo_hi", "callset"}'`. After that change `--check` is clean. Do not match the bare prefix.
- **Census is not a `def_day.sh` child of `def_census.py`.** The override is one argv word on the setup `def_pipeline.py` invocation, and the pipeline forwards `--max-zero-carrier` only when that word is present. A future edit that calls `def_census.py` directly will drop the override.

## What this VM could not test

No gitignored data (`data/challenge`, donors, the A/B/C callblock artifact). Not run:

- 528's D3 byte-compare against the real callblock artifact.
- Real `def_lodo_key.py` (the sidecar `hi` field and the npz). The day stub writes a json `{"hi": <from --hi or 3>}` and the argv is what was checked. gate524 reads that stub sidecar.
- Real `def_census.py` and a real zero-carrier count. The log line uses a stub `census.json` with `consensus.targets_with_0_donors = 2`. The pipeline dry-run checks argv forwarding only.
- The inert known answer, `cell-eval2`, and `vcc`. Nothing was submitted.
- The literal path `~/Code/virtual-cell-challenge-2026`. Derivation was checked by resolving `dirname` of the temp script; equality with the Mac checkout is by construction, not by execution.

## Coordinator's local check (2026-09-30 11:40Z, on the Mac, before merge)
- Branch `cloud/538-bundle-extras` (3d7b4bbd): only these four files were taken onto master (`git checkout <branch> -- …`).
- `test_extras538.py`: **13 tests OK in 24.7 s** on this checkout (real repo files copied to temp trees).
- `apply_october538.py --root . --check`: `patches apply cleanly (not written)`; the nine live files hash identical before
  and after the tests and the check (`shasum -a 256`, compared byte for byte).
- Applied on 10-05 only, after `recipe_dr524.sh`'s process exits (the commands above), then the inert KA and the committer.
