# 530 — October pass-through bundle (patches only; not applied)

The five knobs from 529's plan are one string-anchored patch, `patch_bundle530.py`. It was not applied to this
checkout. `scripts/make_submission_pm_cd4.py` and the emitter were not edited. 529's `DEF_DAY_CALLSET_BLOCK` is not
in this bundle. 526's six-line "builder is not HEAD" insertion is not part of the patch either: none of the anchors
touch the gap after `BSHA0=`, and the same patch applies whether or not those lines are already present (HEAD does
not have them; a copy with them inserted also applies).

Line numbers below are the files **after** the patch, measured on a temp copy. They are not the live scripts.

## What each knob does

Unset means the variable is absent. A set value is validated and then changes one export or one argument.

1. **`DEF_DAY_RHO_SCALE=<s>`** (`0.9 ≤ s ≤ 1.1`). The re-level target is ρ* × s instead of ρ*.
   - `scripts/def_day.sh:119` refuses a scale outside the interval, and only when the variable is set.
   - `scripts/def_day.sh:344-345` logs `re-level rho-scale=<s>` and appends `--rho-scale <s>` to `relevel_iter.py`.
     The append is `${DEF_DAY_RHO_SCALE:+...}`, so an unset variable adds no word.
   - `scripts/def_day.sh:365` passes the same scale as an extra positional to `inband`; `:388` passes `--scale` to
     `best_safe`. Both expansions are empty when unset.
   - `experiments/448-def-rehearsal-q175/relevel_iter.py:29-35` adds `--rho-scale` (default 1). The target is
     `base if scale == 1.0 else base * scale`, with `base = exp(2 * sd * Φ⁻¹(q))`, `q = 0.31/(0.31+1.78)`, `sd = 0.01`.
   - `scripts/def_day_helpers.py:30-37` is the same branch for `rho_star`. `inband` (`:134`) takes an optional second
     positional; `best_safe` (`:148`) takes a leading `--scale`. Omitted, both call `rho_star()` with scale 1.

2. **`DEF_DAY_CD4_STACK=x1|x1.5`**. x1 is today's dump. x1.5 is 457's one-off, inside the dump subshell only.
   - `scripts/def_day.sh:128` accepts only `x1` or `x1.5` (unset counts as `x1`; the case runs no external command).
   - `scripts/def_day.sh:277-291`: x1 writes `dump_h1x2_cd4x1.npz` and does not export `VCC_CD4_TARGET_G`. x1.5 writes
     `dump_h1x2_cd4x15.npz` and exports `VCC_CD4_TARGET_G=0.5165482917974462` (`1.5 * 0.3443655278649641`, the literal
     in `experiments/457-cd4-half-def/def_day_cd4x15.sh`) in that subshell only. Level builds keep DEF_ENV's
     `0.3443655278649641`. A level meta whose `cd4_target_g` moved is a leak and fails `riders_match` and gate524.

3. **`DEF_DAY_RC_CALLSET_CPM=<n>`** (`100 ≤ n ≤ 5000`). `scripts/def_day.sh:123-126` checks the integer and that the
   builder contains `VCC_RC_CALLSET_CPM`. `scripts/def_day.sh:306` exports it **after** the dump subshell, with the
   other level-build riders. Unset leaves the builder's own 2200 (the key is always written; this script does not
   export the variable). `scripts/def_day.sh:436` adds `rc_callset_cpm` to the postbuild allow-list only when the
   variable was exported. The day script still wipes every incoming `VCC_*` at start, so a parent
   `VCC_RC_CALLSET_CPM` does not survive an unset knob.

4. **Arc hook** (517's `patch_defday517.py` and `patch_sameline517.py`, folded in).
   - `DEF_DAY_MDE_SIGN` is `clean4:<tau>`, `programme:<tau>[:<tau_c4>]`, or
     `hookbank:<bank>:<tau>[:<tau_c4>[:le<N>]]` (`scripts/def_day.sh:98-115`).
   - `hookbank` builds a second same-line artifact at `$W/hookprog.npz` (`scripts/def_day.sh:232-242`) and exports
     `VCC_DE_MDE_SIGN=file:<npz>:<sha>:<rest>` (`:305`). The meta records source `file`, not `hookbank`.
   - `DEF_DAY_HOOK_ESTIMATOR=meanz_ne` (`:143-149`) is passed as `--programme-estimator` only when set (`:239`).
     Default pdown does not enter `programme_meanz_ne` and adds no provenance key
     (`scripts/def_sameline.py:182`, `:272-274`, `:293`).

5. **`DEF_DAY_LODO_STRENGTH=<theta>`**. A non-zero plain decimal in `[-1, 1]` (scientific notation is refused, so the
   filename stays a plain token). Requires `DEF_DAY_LODO=1` and `DEF_DAY_LODO_HI_PCT` (`scripts/def_day.sh:132-141`).
   - `scripts/def_day.sh:255` saves `LJSON` from the LODO key **before** the rider. The cap's `hi` is read from that
     file (`:269`). 522's sidecar has no `hi`.
   - `scripts/def_day.sh:257-262` runs `experiments/522-october-bulk-key/def_strength_key522.py` and points the pin's
     `--pert-scale` at `pert_scale_lodo_strength.npz`. Unset, `LK` is still the key just written.
   - `scripts/def_day.sh:283-285` renames the dump to `dump_h1x2_cd4x1_s<theta>.npz` (or `cd4x15` when both are on;
     `-` becomes `m`). The level meta's `pert_scale_sha256` stays the 410 key, because `VCC_PERT_SCALE` is still
     exported only inside the dump subshell. The theta is carried by the dump filename, the log line
     `LODO strength rider: theta <theta>`, and `work/run/pert_scale_lodo_strength.json`.

`riders_match` (`scripts/def_day_helpers.py:171`) stays a 3-argument call when none of rc / non-x1 stack / rho /
strength / hook-estimator is set (`scripts/def_day.sh:397-411`). The 3-argument path uses 517's grammar
(`_mde_sign_request`), so `clean4`, `programme` and the exported `file:` string are checked, including sha, tau,
tau_c4 and `le`. A 4th argument `rc=<-|N>,stack=<x1|x1.5>,rho=<-|s>,strength=<-|theta>` (`_october530_spec`,
`:74`) checks the dump filename (`_dump_name530`, `:66`), `cd4_target_g == 0.3443655278649641`, `rc_callset_cpm`
(2200 when rc is `-`), and `rho_scale` only if the meta actually has the key. The builder does not record the loop
target; a wrong recorded value fails, an absent key passes.

## Resolver, decisions, gate

`docs/october-decisions.json`: `rho_scale` (`:105`), `rc_callset_cpm` (`:113`), `hook_estimator` (`:121`) and
`lodo_strength` (`:129`) are **PENDING**. `cd4_stack` stays **DECIDED `x1`** (`:70-74` evidence mentions 530). x1 is
not exported, so today's `--allow-pending` command is unchanged. `mde_sign` stays PENDING and gains `if_hookbank`
(`:103`).

`scripts/october_recipe.py:41-44` maps the four new switches. `cd4_stack` is handled beside `ENV_OF` (`:82-95`):
DECIDED `x1.5` exports `DEF_DAY_CD4_STACK=x1.5`; DECIDED `x1` exports nothing; anything else is a problem. A
DECIDED theta of `0` is a problem (`:93-95`) even though `is_on` treats `"0"` as off, because theta 0 still rewrites
the npz. `_problems_530` (`:107`) checks ranges and the hook/strength dependencies. `installed_problems` matches the
export line, not a comment (`:149-168`): MDE is `export VCC_DE_MDE_SIGN="$MDS"`. An unpatched `def_day.sh` still
refuses `x1.5` (`cd4_stack` is in the message).

`experiments/524-october-recipe-dr/gate524.py:47-50` allows the meta keys each knob can move. `_mde_ok530` (`:101`)
replaces the one-colon partition: `clean4`, `programme` (tau_c4 defaults to 0.05) and `hookbank` (source must be
`file`, sha length 64, `le` present exactly when the spec has `le<N>`). `:173-185` checks the composed dump name
and that `cd4_target_g` stayed at DEF_ENV; with both knobs off, `bulk_load` must still be `dump_h1x2_cd4x1.npz`
(that name is in `BASE_ALLOWED`, so a silent rename would otherwise pass the diff). Rho (`:186`) matches
`rho-scale=<s>` with a non-digit boundary, so `1.02` does not match `1.020`. Strength (`:189`) reads the sidecar
theta and the log line. The hook estimator (`:203`) must be `meanz_ne`, the log must contain `estimator meanz_ne`,
and `mde_sign` must start with `hookbank:`.

## Bit-identity when unset

| Knob | What stays today's |
|---|---|
| rho | No `--rho-scale`, no extra `inband`/`best_safe` word, no log line. Scale `1.0` **inside** `relevel_iter.py` / `rho_star` does not multiply (`base if scale == 1.0`), and omitted stdout equals `--rho-scale 1` byte for byte. Setting `DEF_DAY_RHO_SCALE=1` is a real argv change; the identity path is the unset variable. |
| CD4 | Unset and explicit `x1` both write `dump_h1x2_cd4x1.npz` and do not export `VCC_CD4_TARGET_G`. The case statement does not spawn a process. |
| RC | No `VCC_RC_CALLSET_CPM` export, no `rc_callset_cpm` on the allow-list. A parent export is wiped with the other `VCC_*`. |
| hook | The `hookbank` case does not match. `clean4` still exports the request (`MDS` is the request). `meanz_ne` is not entered; the provenance splat is `{}` when the estimator is pdown. |
| strength | `LJSON` is the same path the old log read (`${LK%.npz}.json`). The pin's `--pert-scale` is the LODO key. The dump name is not suffixed. |
| riders | The `else` branch is today's two-argument `riders_match` (mag, sign). A renamed `bulk_load` does not fail that call. |

Proved by running the patched and unpatched scripts against the same stubs (real helpers and `relevel_iter.py`, stub
`uv` for the driver, same-line, LODO, ntclevel, postbuild and audit). With no new knob, and again with today's
October env (`DEF_DAY_LODO=1`, `DEF_DAY_DE_MAG=dn_k250_zs`), the `uv` argv list and the `VCC_*` environment of every
call match, and so does `def_day.final.env`. A pre-set `VCC_RC_CALLSET_CPM=9999` is absent from both. Explicit
`DEF_DAY_CD4_STACK=x1` matches the bare patched run. `DEF_DAY_RHO_SCALE=1.2` exits before any `uv` call.

## Tests

`uv run python experiments/530-cloud-october-bundle/test_bundle530.py` on copies only (synthetic inputs; the live
seven files' sha256 are compared at the end).

```
test_applies_with_and_without_526_insertion ... ok
test_bad_scale_stops_before_any_build ... ok
test_bash_n_and_strength_regex ... ok
test_check_does_not_write ... ok
test_decisions_new_switches_are_pending_and_cd4_stays_x1 ... ok
test_gate_rows ... ok
test_knobs_on_change_only_their_exports ... ok
test_partial_bundle_refuses ... ok
test_relevel_and_helper_scale ... ok
test_resolver_every_combination_of_the_new_switches ... resolver combinations 32 ok
test_resolver_rejects_bad_values_and_unpatched_cd4 ... ok
test_resolver_stock_allow_pending_is_todays_command ... ok
test_riders_match_accepts_and_rejects ... ok
test_second_apply_refuses_and_keeps_bytes ... ok
test_unset_matches_todays_exports_and_arguments ... unset identity: 15 uv calls match; october 17 match
test_zz_live_tree_untouched ... ok

Ran 16 tests in 19.564s
OK
```

Covered: `bash -n`; `--check` does not write; a second apply exits `already patched` and does not change bytes; a
partial marker set is refused; apply with and without 526's six lines; the resolver's 32 on/off combinations of the
new switches, plus boundaries (`0.9`, `1`, `1.1`, `100`, `2200`, `5000`) and the refusals (scale, cpm, `x2`, hook
without `hookbank`, theta `0` / `1e-1` / out of range, strength without LODO or the cap, a raw `file:` mde, and
`x1.5` against an unpatched `def_day.sh`); stock `--allow-pending` prints exactly
`DEF_DAY_LODO=1 DEF_DAY_DE_MAG=dn_k250_zs bash scripts/def_day.sh B W N 1.75 0.5 noprog 4`; `riders_match` accepts
the matching meta and rejects each mismatch, including 526's tau `0.40` false pass staying a fail; gate524 passes
an all-knobs fixture and fails the rho prefix, the CD4 leak, a silent x1.5 rename, the wrong cpm, a `hookbank`
source in the meta, a missing or extra `le`, and a missing strength sidecar. A dropped same-line record still fails.

## How the coordinator applies it

Do this only after `recipe_dr524.sh`'s **process** has exited. Do not edit `def_day.sh` while that process is
running (bash reads the script as it goes).

```bash
uv run python experiments/530-cloud-october-bundle/patch_bundle530.py --root . 
bash -n scripts/def_day.sh
# helper suites that already exist: 497, 502, 517
# then the inert known answer: efp relaxed through def_day.sh with none of the new knobs set
# must reproduce ac49d319f5982507003b1f2070660d2d03ee351fd13c9fdd30df8f5fb08ec0b3 byte for byte, as 497 and 502 did
# then one rider-on rehearsal on relabelled A/B/C, with every knob that october_recipe.py decide has turned on
uv run python scripts/october_recipe.py check --allow-pending
```

A second run of the patch exits `already patched` before writing. `--check` parses and syntax-checks and does not
replace the targets. New switches stay PENDING until `october_recipe.py decide`; until then the `--allow-pending`
command is today's. The commit after the inert KA goes through a spec file (`scripts/commit_after_chain.sh`), not
from this experiment.

## Coordinator's adaptations on merge (2026-09-29 14:15Z)
- `patch_decisions` is now a JSON-level edit: the string anchors broke after `docs/october-decisions.json` gained 534's
  `held_out` lines (b8681568). It sets `mde_sign.needs` / `if_hookbank` and the `cd4_stack` evidence, and inserts the four
  PENDING switches (each with a `held_out` line) before `cd4_thresholds`. Output keeps the file's format (indent 1, UTF-8).
- `test_bundle530.py`: ROOT is the checkout (was `/workspace`); the stubbed day runs commit their stub tree in a throwaway
  git repo because the live `def_day.sh` now carries 526's refuse-not-HEAD check (600adea8); the "with/without 526"
  test tests whichever variant the live tree is not.
- Local result on this tree: `patch_bundle530.py --root . --check` → "patches apply cleanly (not written)";
  `test_bundle530.py` → Ran 16 tests, OK (resolver combinations 32; unset identity 15 / October 17 uv calls match).
- NOT applied. Apply after `recipe_dr524.sh`'s process exits (10-05), then the inert KA (efp relaxed == ac49d319) and one
  rider-on rehearsal, then commit via a spec file.
