# 526 — October-path review, round 2 (2026-09-29)

Reviewed the October (D/E/F) path since `2d1fc214`: `8960c102` (MDE sign hook), `39d38196` (tercile tie rules), `b1162b13` (DE mag / MDE sign pass-throughs), `3eb21c57` (`VCC_RC_CALLSET_CPM`), `22cbca95` (524). Also re-read the Oct 22 dependencies that have no commit in that range: `def_pipeline.py`, `def_sameline.py`, `def_lodo_key.py`, `def_preflight.py`, `def_pin_panel_artifacts.py`, `lib438.py`, and `def_day.sh`'s re-level caller. `experiments/448-def-rehearsal-q175/relevel_iter.py` has no commit since `2d1fc214`. No challenge h5ad, donor table, or harness matrix was read. Nothing here was applied to `scripts/` or to another experiment. `vcc` was not run.

Repros (each exits 0 while the defect, or the classification, is still in the tree):

```text
uv run python experiments/526-cloud-october-review2/repro/repro_q1.py
uv run python experiments/526-cloud-october-review2/repro/repro_recipe.py
uv run python experiments/526-cloud-october-review2/repro/repro_bitid.py
uv run python experiments/526-cloud-october-review2/repro/repro_sharpen.py
uv run python experiments/526-cloud-october-review2/repro/repro_gate524.py
uv run python experiments/526-cloud-october-review2/repro/repro_pin.py
uv run python experiments/526-cloud-october-review2/repro/repro_524wait.py
```

All seven exited 0 on this checkout. `git apply --check` passed for every patch below. The gate patch, applied only in a temp tree, turns the two false passes into fails and still passes the matching artifact.

## Findings

### 1. `gate524.py` passes an artifact that does not carry the recipe — would silently degrade the 10-06 read

`experiments/524-october-recipe-dr/gate524.py:112-115` and `:134-135`.

The docstring (`gate524.py:8-9`) requires every ON switch to show up with its resolved value, and the same-line record to agree with the recipe. The checks do less than that.

| What the recipe asked | What the meta had | Gate |
|---|---|---|
| `w 0.5`, flags `noprog` | no `sameline` key at all | `ALL GATES PASS`. The key is in `BASE_ALLOWED` (`:36`), so the diff against 453dr is allowed, and the `else` only appends a line |
| `DEF_DAY_MDE_SIGN=clean4:0.05` | `de_mde_sign.tau = 0.40` | `ALL GATES PASS`. The check is `bool(de_mde_sign)` |
| `DEF_DAY_PRIOR_GATE=le1` | `de_prior_gate.le = 99` | `ALL GATES PASS`. The check is `bool(de_prior_gate)` |
| `DEF_DAY_PROGRAMME_SOURCE=h1_2025_full` | log never names that bank | `ALL GATES PASS`. The meta stores `prog`, not the bank (`make_submission_pm_cd4.py:5656-5661`) |

`repro_gate524.py` builds those two trees, runs the real `gate524.py`, and gets `ALL GATES PASS` on both. A third tree whose values match the recipe also passes, so the harness is not failing everything. `def_day_helpers.py riders_match` (`:100-109`) rejects the tau 0.40 meta (rc 1). `def_pipeline.py september_rows` (`:244-247`) checks `de_mag.variant` and does not check `de_mde_sign` or `de_prior_gate`. A dropped same-line record is a run `def_day.sh` will finish: no registered bank logs "same-line mode off" and continues (`def_day.sh:167-171`). On the relabelled rehearsal that is an integration change versus 453dr, and the falsifiers in `524/PREREGISTRATION.md` are there to catch it. The gate currently blesses it.

A per-target sharpening assignment is also invisible. The sharpen check (`gate524.py:104-108`) compares `mode` and the q list. `repro_sharpen.py` shows 511 moving a lower-cut tie from q 2.0 to q 2.5 while `mode` and the q list stay `stercile / [2.5, 2.0, 1.75]`. Both metas would pass. There is no oracle for the right assignment in the gate; the patch does not invent one.

Patch: `patches/gate524_values.patch`. Compare tau, compare `le` or `q`, fail when the record is absent, require `programme-source <bank>` in `def_day.log`. Not applied.

### 2. A control file with no `ntc_id` column passes preflight and changes the dispersion — would silently degrade

`scripts/def_preflight.py:133-137` checks `n_ntc_ids` only when the column exists. `repro_q1.py` builds a 3-context, 18,533-gene fixture with every other check passing and the column omitted: 0 failed checks.

`scripts/make_submission_pm_cd4.py:5220-5225` then takes the other branch and sets `C_p` from `random_draw_cp` (`:5191`), the 096 random proxy. The submission is a normal count matrix. A column that is present with the wrong category count fails closed (the same repro: 2 categories against a manifest of 46).

Patch: none. The day should fail closed in preflight when `per_context[c].n_ntc_ids` is set and the column is absent, instead of falling through to the random proxy. Left unpatched because the fallback is an intentional builder branch; the missing piece is the preflight row, and I did not want a half-edit of `def_preflight.py` proposed without the manifest field's optionality spelled out. The classification is the deliverable.

### 3. A reboot mid-pin leaves the builder pinned, and a fresh `def_day` then treats the pin as the original — would break Oct 22

`scripts/def_pin_panel_artifacts.py:207-225`. `tmp.replace(BUILDER)` runs, and only then does the script write `pin-artifacts.json` and `pin-artifacts.builder-before.py`. `repro_pin.py` replays that order: a kill after the replace leaves the pinned bytes and `restored_to_head is None`. The runbook command (`docs/DEF-RUNBOOK.md:157`) needs the record.

`def_day.sh:100-120` restores only from `$W/pin-artifacts.json`, and only from the EXIT trap. A reboot does not run the trap. The trap does not read `pin-lodo.json`. Restoring `pin-artifacts.json` is the right record when the before-copy exists: the second pin writes a different sibling (`pin-lodo.builder-before.py`), so the first before-image is still HEAD.

A fresh work dir does not notice. `BSHA0` (`def_day.sh:67`) is the sha of whatever is on disk. The new pin's before-image is that already-pinned file, the trap writes it back, the sha matches `BSHA0`, and the trap exits 0. `repro_pin.py` shows `match=True` on the pinned bytes. `recipe_dr524.sh:17-25` will not restore either: `eq_head` waits until `202610052000` and exits 1.

The replace itself is atomic (`os.replace`). The builder file is not half-written. What is missing is the record.

Patches: `patches/pin_record_before_replace.patch` (write the copy and the JSON, then replace) and `patches/def_day_refuse_pinned.patch` (refuse at start when the on-disk builder is not `git show HEAD:scripts/make_submission_pm_cd4.py`). Not applied.

### 4. `recipe_dr524.sh` can start beside a live Python build, and a failed gate exits 0 after `vcc prep` has overwritten the package — would break the 10-06 candidate

`experiments/524-october-recipe-dr/recipe_dr524.sh:18` and `:50-56`.

`builders()` keeps lines whose comm is `bash` and whose arguments contain `def_day.sh` (or `chain4xx` / `chain5xx` / `run_arm_generic.sh` / …). `repro_524wait.py` feeds it eight lines. It hits `bash scripts/def_day.sh …` and `bash -c echo scripts/def_day.sh is mentioned`. It does not hit `python scripts/make_submission_pm_cd4.py` or `uv run python scripts/make_submission_pm_cd4.py` or `python scripts/build_union_whitener.py`. `def_pipeline.other_builders` (`def_pipeline.py:360-365`) is the check that sees those. It also does not hit its own argv or `chain1006.sh` (the §6.2 self-match). The wait caps are timestamps (`FLOOR` / `CAP`), not iteration counts.

The script is `set -uo pipefail` (`:9`), not `set -e`. After `gate524 rc=1` the next command is `log "END recipe_dr524 …"`, so the process exits 0. The repro's copy of that tail exits 0. `chain1006.sh` does not read that rc. It reads `gates524.log`, which `--write` does not replace on failure (`gate524.py:149-150`). Two consequences:

- A re-run that dies before prep leaves the previous package and the previous log. `ready()` (`chain1006.sh:36-41`) still blesses that older package. That part is the previous success, and it is the right package.
- Prep (`recipe_dr524.sh:50`) overwrites `$STEM.vcc` before the gate. A re-run that dies at the gate leaves the new bytes and the old log. `ready()` compares shas and refuses. The good package is gone, and 10-06 slot 1 falls through to `446bj2zk16`. The sha check stops a wrong submit. It does not save the candidate.

`chain1006.sh:17-19` can also `break` at 23:30Z while `recipe_dr524.sh` is still alive. `arm1006.sh:43-44` re-checks the sha at submit time, so an overwrite between the decision and the submit is refused. The collision with a Python builder is the one that is not refused.

The READY greps match the lines that are actually written. `gate524.py:70` writes `[524dr] CLIFF OK` and `:85` writes `[524dr] PASS`. `run_arm_generic.sh:83` and `:96` write `[446bj2zk16] PASS` and `[446bj2zk16] CLIFF OK` with no timestamp; `novelty rc=0` and `package sha256:` go through `log()` and carry one. `chain1006.sh:37-39` and `arm1006.sh:16-18,43` use the same three patterns, and `awk '{print $NF}'` on `package sha256:` is the digest. `repro_524wait.py` greps both log shapes and gets a match, including the K-arm tag. `arm1006.sh:31` uses `grep -qE "NOT SUBMITTING| END$"`. No `declare -A` in the three scripts.

Patch: `patches/recipe_dr524_fail_closed.patch`. `builders()` also calls `other_builders()`, a failed gate restores `$STEM.vcc.prev` and exits 1. Not applied.

## Question 1 — a D/E/F bundle that is not A/B/C

`repro_q1.py` drives `def_preflight.preflight` on synthetic bundles (4 control cells, the real 18,533-gene width so the density check is in range) and `def_census.gate` on a fake census dict. Same-line is the real `def_sameline.py` on a two-gene bundle.

| Difference | What happens |
|---|---|
| Context labels other than D/E/F, files and manifest agree, three identifier names | Proceeds. `def_day.sh:61` and `:137-141` read the manifest and require `len == 3` and `str.isidentifier`. Not hardcoded |
| Context count ≠ 3 | Fail closed. `def_preflight.py:189`. The repro's two-context bundle fails "manifest declares 3 contexts" |
| `gene_names.csv` / `pert_counts.csv` header missing or the wrong column name | Fail closed. `def_preflight.py:204-205` requires `columns[0] == "gene_name"`; `:219` requires `manifest["pert_col"]`. The repro's headerless gene file fails |
| Gene order: `.var` index ≠ `gene_names.csv` | Fail closed at preflight (`:113`) and again in the builder (`make_submission_pm_cd4.py:5215-5216`) |
| Duplicate targets | Fail closed (`def_preflight.py:223`) |
| `obs['context']` label ≠ the filename's context | Fail closed (`:126`). Repro writes `NOT_D` into `context_D.h5ad` |
| `ntc_id` present, category count ≠ `per_context[c].n_ntc_ids` | Fail closed (`:137`). The count is `len(categories)`, not how many codes are used |
| `ntc_id` absent | Silent. Finding 2 |
| `cells_per_pert` follows the manifest | Preflight requires the key (`:187`) and, when an `n_cells` column exists, that it equals the manifest (`:229`). It does not hardcode 400. `build_submission.py:77` uses `manifest.get("cells_per_pert", 400)`. The 400-cell rule is `vcc prep`, which this review did not run. `def_day.sh:350` dry-runs prep and dies if that stage exits non-zero |
| Panel length ≠ 300 | Fail closed once LODO runs. `lib438.py:205` (`assert len(panel) == 300`) inside `bind()`, which `def_lodo_key.py` uses. LODO is DECIDED on. A length-300 panel of different symbols passes that assert |
| Gene axis reordered, bundle path is `data/challenge` | LODO's `assert official == OFFICIAL` (`lib438.py:203`) compares the file to the tuple `plolo.py:102-106` read from the same path at import. Same path, same bytes, the assert holds. The CD4 check does not: `DEF_ENV` sets `VCC_CD4_DEFLATE_BASIS=V5_abc` off A/B/C (`def_pipeline.py:104-107`), and `_verify_cd4_deflation_basis_axis` (`make_submission_pm_cd4.py:3166-3174`) raises unless the basis axis equals the bundle's symbols in order. The repro execs that function on a swapped axis and gets `ValueError` |
| `pert_counts.csv` columns swapped (`n_cells,target_gene`) | Preflight and `lib438.bind` read `target_gene` by name. `plolo.py:97-99` reads column 0 at import (`"400"` in the repro) and then `bind()` replaces `plolo.PANEL` (`lib438.py:216`) before the key is fit |
| No same-line bank | Valid artifact, same-line off. `def_sameline.py:187-191` logs "no registered bank" and leaves the context out. The repro's `--force D=NOT_A_LINE,E=ALSO_NO,F=NOPE` exits 0 with `nctx 0`. `def_day.sh:171` logs it and continues. A programme source that yields 0 contexts dies (`:168-169`). An unknown programme-source name dies (`def_sameline.py:240-241`, repro rc 1). A source with fewer than 150 rows left or coverage under 0.5 dies (`:245-249`); that floor was not executed here because the donor tables are not on this VM |
| Context matched to a bank already in the donor mix | Proceeds at the recipe's w. `REGISTRY` (`def_sameline.py:56-57`) maps K562 to `k562_gwps` and HCT116 to `orion_hct116`. Those are `PRODUCTION_DONORS` (`def_census.py:55`). `def_sameline.py` does not import that tuple. The decided w is 0.5 (`october-decisions.json`), not `auto` (which would use `w_strong` 2 for those banks, `def_sameline.py:75-76,152`) |
| A registered line whose call is not unanimous | Context unchanged (`def_sameline.py:180-185`). Not a crash |
| Targets absent from every production donor | 0–5 such targets: census passes (`def_census.py:69-70,237-239`). The repro's `gate()` passes 5 and fails 6. Those rows are `MultiDonorModel` `fallback="mean"` (`def_census.py:14-16`): a real prediction, the panel-mean direction, not a zero row. More than 5, or a donor under 100 targets, or CD4 union under 120, or CD4 complete under 100: the gate fails and the driver stops |
| Stale whitener left in the work dir | Fail closed. `verify_existing_operator` (`def_pipeline.py:425-460`) refuses to skip the fit when `axis_sha256` differs or is missing. `n_cols` is explicitly not the check |

## Question 2 — resolver and pass-throughs

`repro_recipe.py` re-executes `october_recipe.py` against copies of `docs/october-decisions.json`.

Live file, `--allow-pending`: rc 0, env exactly `{DEF_DAY_LODO: 1, DEF_DAY_DE_MAG: dn_k250_zs}`, args `… 1.75 0.5 noprog 4`. That is q, w, flags `noprog`, max_levels 4, LODO, zs. Pending sharpen (`if_s_p2` = stercile), percentile cap, programme, prior gate, and MDE sign are off.

Live file, no `--allow-pending`: rc 2, stdout contains `# REFUSED ` plus the day-of CD4 comment (`october_recipe.py:123-127`). The refused line still names `DEF_DAY_DE_MAG`. Honoring the rc does not run it.

| Combination | Result |
|---|---|
| Programme `h1_2025_full` on, prior still pending, `--allow-pending` | rc 0. Prior absent. Flags are `""` (`with_programme_source`), not `pt:noprog` |
| Prior `le1` on, programme still pending | rc 3. Problem: prior without programme. `def_day.sh:150-151` dies on the same pair. The builder dies too if `VCC_DE_PRIOR_GATE` is set and sameline has `noprog` or is unset (`make_submission_pm_cd4.py:698-699`) |
| Programme, prior `le1`, MDE `clean4:0.05`, pct `96.64697522769474`, sharpen `stercile:2.5,2.0,1.75` all DECIDED | rc 0. All five exports present. Flags `""` |
| Programme on and `with_programme_source` contains `noprog` | rc 3. `def_day.sh:147` dies as well. The builder would keep CLEAN4 (`:674`, `prog = "noprog" not in flags`) |
| `lodo_hi_pct` on, `lodo` off | rc 3. `def_day.sh:78-79` dies as well |
| `cd4_stack` `x1.5` | rc 3. The dumper is H1×2 / CD4×1 |
| Bash `${6-pt:noprog}` | An empty 6th argument stays `""`. An omitted one becomes `pt:noprog`. The resolver always passes the argument (`october_recipe.py:83`) |

`def_day.sh:66` unsets every `VCC_*`, then `:210-216` exports a rider only when the `DEF_DAY_*` variable is non-empty. The repro's wipe leaves `VCC_DE_MAG=dn_k250_zs` and `VCC_DE_MDE_SIGN`, `VCC_RC_CALLSET_CPM`, `VCC_DE_PRIOR_GATE` unset. `VCC_RC_CALLSET_CPM` is not in `ENV_OF` (`october_recipe.py:33-41`) and not in `def_day.sh`. A shell value cannot ride into the October build. `gate524.py` has a `DEF_DAY_RC_CALLSET_CPM` switch the resolver cannot emit; an unknown export makes the gate exit 1 (`:51-53`).

The shell's `DEF_DAY_DE_MAG` check allows any `[a-z0-9_]+` (`def_day.sh:86`). The builder's regex refuses anything that is not `<sel>_k<k>_zs` with k in range (`make_submission_pm_cd4.py:756-758`). `clean4:abc` passes the shell (`:90`) and `SystemExit`s in the builder (`:789-791`). Neither is ignored.

`is_on` treats `off`, `0`, `none`, and `""` as off (`october_recipe.py:45,57-58`). The string `false` is ON. The builder then refuses `VCC_DE_MDE_SIGN=false` (`:789`). Fail closed one layer later. Cosmetic at the resolver.

## Question 3 — unset means bit-identical

`repro_bitid.py` and `repro_sharpen.py`.

| Knob | Unset path | Config key |
|---|---|---|
| `VCC_DE_MDE_SIGN` | The block at `make_submission_pm_cd4.py:784` is not entered. `DE_MDE_SIGN` stays `None`. The config splat at `:6893` adds nothing. The emitter reads `mde_add_sign` only inside `if "mde_add_sign" in lfc_target` (`rc_emit_k.py:999`). No RNG on that branch | Key absent. `riders_match` requires it absent (`def_day_helpers.py:101-104`). `def_day.sh:320` adds `de_mde_sign` to the postbuild allow-list only when the variable is set |
| `VCC_DE_MAG` | Same shape (`:748`, splat `:6907`). The scorer-control block is behind `DE_MAG is not None` | Key absent unless set (`def_day.sh:319`) |
| `VCC_RC_CALLSET_CPM` | `RC_CALLSET_CPM = 2200` (`:1961`) and the `if` at `:1964` requires a non-empty string. `"2200"` would assign the same int. `isdigit()` rejects a sign and a decimal | `rc_callset_cpm` is written either way (`:6971`). 453dr already has 2200, so the diff is empty |
| `VCC_BULK_SHARPEN_Q` scalar (the knob `DEF_DAY_SHARPEN_PT` off) | `sharpen_q_by_tercile` is not called. The carrier array is built only `if BULK_SHARPEN_PT is not None` | No `bulk_sharpen_pt` key. 511 does not run |
| `VCC_BULK_SHARPEN_Q=stercile:…` | This is the ON path. A lower-cut tie block of size ≥ 2 goes LOW (`make_submission_pm_cd4.py:4276-4280`). The previous function left that edge to `np.digitize`, which put it in MID. The repro: old q `[2.0, 2.0, 2.0, 2.0, 1.75, 1.75]`, new q `[2.5, 2.5, 2.5, 2.0, 1.75, 1.75]`. A no-tie vector matches. A NaN is excluded from the quantile and assigned MID; the previous function filled it with `nanmedian` and included it, and the cuts move (old edges 2.6667/3.3333, new 2.3333/3.6667 on the repro's vector). Audit keys are unchanged: `mode, q, edges, counts, floor_to_low, tie_to_mid, nonfinite` | The key records mode and q, not the per-target assignment. Finding 1 |

`de_sign_emitter_sha256` moves on every build, MDE sign unset included. The constant is `1658ae2c400f75ee1659a8ffcaebca569bf3c8afd55ffd1c3006414344ef2791` (`make_submission_pm_cd4.py:652`). `sha256(src/vcc2026/rc_emit_k.py)` equals it. 446b's meta is `edb5b8dc…`. The value is written inside the `DE_SIGN_PRIOR` block (`:6879`), and that prior is in `SHIPPED_ENV`. `def_day.sh:321-326` adds the key to the allow-list exactly when the source pin differs from the incumbent's recorded sha, which it now does. `def_postbuild.py:465-468` then allows any difference under that key name. A finished build cannot record a third sha: `de_sign_prior_full` (`:1188-1190`) raises unless the file hashes to the constant, and the config writes the constant. Output counts do not change, because the new emitter code is the `mde_add_sign` branch and that branch is not entered.

`def_postbuild.py:319-323` still pins `abc_targets` to `ABC_TARGET_SET_SHA` (the 519 fix). It does not re-read the live `pert_counts.csv`.

## Question 4 — the 524 chain

Findings 1 and 4. The READY line formats match; see finding 4.

`recipe_dr524.sh` uses `--allow-pending` on purpose (`524/PREREGISTRATION.md:9-16`): a switch still PENDING at 10-05 02:00Z runs off and is listed in `recipe524.json`. The gate exits 1 if that JSON has `problems` or an export outside `SWITCH_KEYS` (`gate524.py:51-53`).

## Question 5 — reboot mid-`def_day`

Finding 3. State left behind:

| When the machine dies | What is on disk |
|---|---|
| During a stage, before any pin | Work dir partial (`driver-*.log`, maybe a half stage). Builder is HEAD. Re-run wants a fresh work dir (`def_day.sh:55-56`): it deletes a tag's h5ad before rebuilding and will not treat a pre-existing file as the build it just asked for |
| After `tmp.replace`, before the JSON and the before-copy | Builder pinned. No record. Runbook restore cannot run. `git checkout -- scripts/make_submission_pm_cd4.py` still can, because the pin is not committed (`def_pin_panel_artifacts.py:217` tells the operator to commit it afterwards) |
| After both pins, trap not run | Builder at the LODO pin. `pin-artifacts.builder-before.py` is HEAD. `pin-lodo.builder-before.py` is the post-first-pin file. Restoring `pin-artifacts.json` returns HEAD. Restoring `pin-lodo.json` leaves the first pin |
| Fresh `def_day` without that restore | Trap succeeds and the builder stays pinned. Finding 3 |

`recipe_dr524.sh:39` moves `work/run` aside on retry. It does not move the `.vcc`. Finding 4.

## Checked and found correct

- 519's five defects stay fixed in the files this round re-read. `installed_problems` matches the export line for `DEF_DAY_DE_MAG` and `DEF_DAY_MDE_SIGN` (`october_recipe.py:92-96`), not a comment. A nonzero `command` prints `# REFUSED`.
- Resolver refuses prior-without-programme, percentile-without-LODO, programme-plus-`noprog`, and `cd4_stack` other than `x1`. The live `--allow-pending` command is the two DECIDED exports plus `1.75 0.5 noprog 4`.
- `def_day.sh` wipes `VCC_*` before the rider exports. The dump subshell (`:196-200`) ends before those exports, so the riders are on the level builds only.
- Unset MDE sign, unset DE mag, unset call-set CPM, and scalar q do not enter their new branches and do not add a config key (the call-set key was already there, value 2200).
- Emitter file hash equals the new pin. The allow-list movement of `de_sign_emitter_sha256` is the 499 F2 case, now actually different from 446b, and it is gated on that difference.
- `chain1006` / `arm1006` greps match `gate524.py` and `run_arm_generic.sh`. Timestamp caps, process-exit wait (`chain1006.sh:17-18` checks the pid, with a 23:30 break), `grep -qE "NOT SUBMITTING| END$"`, no `declare -A`, and `recipe_dr524.sh` is not named `chain5NN.sh`.
- Preflight fails closed on headers, duplicate targets, context count, obs label, var order, and a present-but-wrong `ntc_id` count. Census fails closed above 5 zero-donor targets and below the donor/CD4 floors. Programme-source with an unknown bank dies. CD4 `V5_abc` raises on a gene-axis mismatch, which is the check LODO's same-file assert does not provide.
- Same-line off when nothing matches, and w 0.5 on an in-mix bank, are the decided behaviour (`sameline_w` 0.5, evidence in `october-decisions.json`). They are not silent knobs.

## Patches (not applied)

| File | What it changes |
|---|---|
| `patches/gate524_values.patch` | Tau, `le`/`q`, a missing same-line record, and the programme bank in the log |
| `patches/pin_record_before_replace.patch` | Write the before-image and the JSON, then replace the builder |
| `patches/def_day_refuse_pinned.patch` | Refuse a fresh run whose builder is not HEAD |
| `patches/recipe_dr524_fail_closed.patch` | `other_builders()` in the wait; a failed gate restores the previous `.vcc` and exits 1 |
