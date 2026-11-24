# 607b delta review

HEAD `7ed5e52c122c337ac56824b4e8f5460120a699bb` (pinned commit; `git rev-parse HEAD` matched before this review).

Scope is the unapplied 10-06 delta only: 630's module `d8c39c3f` → `7e6b3e52`, 630dd, and 604's octile flip. 607's review of `7dbb1007` is not repeated.

## Findings

### MAJOR — the fail-safe report never sees a `def_day` build

`experiments/630-rank-neutral-up/patch630dd.py:77` reads `$FINAL.rc_audit.jsonl`. Nothing in the day script creates that file.

- `scripts/def_day.sh:66` unsets every `VCC_*` variable on entry. The script never exports `VCC_RC_AUDIT_DUMP` (`grep` of `scripts/def_day.sh` is empty). `build()` (`scripts/def_day.sh:227`) does not set it either.
- The builder writes the jsonl only when that variable is set, and it appends (`scripts/make_submission_pm_cd4.py:942-944`, `:6086-6089`).
- The per-target `ur_*` counters do land in the meta the C check already opens, under `output.rc_callset` (`scripts/make_submission_pm_cd4.py:6190`, emitter `rc_callset` from `scripts/def_day.sh:76`). The report does not read them.

Executed under `def_day`'s real flags (`set -uo pipefail`, `scripts/def_day.sh:57`, no `set -e`), with `def_day`'s `log` and the exact report pipeline cut from `C_NEW`:

- Missing file: exit 0, the log is a `FileNotFoundError` for `…/final.rc_audit.jsonl`, and `FAILSAFE_TOTAL` is absent (0 matching lines).
- Malformed jsonl: exit 0, a `JSONDecodeError` traceback, and again no `FAILSAFE_TOTAL` line.
- A one-row dump with `ur_guard_failed: 1`: exit 0 and `630dd: uprn FAILSAFE_TOTAL 1`.

Failure scenario: every target hits `ur_invariant_failed` or `ur_guard_failed` and keeps its zs block. The meta still records `config.de_uprn.variant`, so the build looks like the pass ran. The coordinator's KA rule (grep `630dd: uprn FAILSAFE_TOTAL`, `n > 0` is a KA FAIL) does not fire, because the line is never printed. `test630dd.py` locks this in: "a missing dump is not fatal" requires exit 0 and does not require the line.

The report code itself is fine when the file exists, and it cannot abort the day (the pipeline ends in `|| true`, and there is no `set -e`). The hole is the missing producer.

Fix, inside `build()` in the 630dd patch, before the driver:

- `rm -f "$W/$N-$tag.rc_audit.jsonl"` (the builder opens with `"a"`, so a retry would otherwise double-count).
- `export VCC_RC_AUDIT_DUMP="$W/$N-$tag.rc_audit.jsonl"`.

`FINAL` is `$W/$N-$tag` with no `.h5ad` suffix, so that path is `$FINAL.rc_audit.jsonl`. Set it inside `build()`, not once next to the `VCC_DE_UPRN` export: a single path would concatenate lvl0, lvl1, … into one file. The driver's non-bare stages keep parent `VCC_*` (`scripts/def_pipeline.py:511-513`), so no driver change is required. Pin stages use `bare_env` and strip `VCC_*`; leave them alone.

Also change the KA token to "the line is present and the number is 0". `n > 0` does not catch a missing line. On `FileNotFoundError` / `JSONDecodeError`, print `630dd: uprn FAILSAFE_TOTAL -1` before returning, so a grep for `FAILSAFE_TOTAL 0` fails closed. Keep `|| true`.

### MINOR — `apply()` labels every `AssertionError` as an invariant failure

`experiments/630-rank-neutral-up/src_rank_up.py:270-275` catches `AssertionError` with no message check and returns the input block with `ur_invariant_failed = 1`.

The post-conditions in `up_pass2` (`:241-251`: library, per-gene total, negative count, U record) and the `rank_mag` helpers the comment names (`cell_bounds` convergence, `_take`) are the intended fail-safe. The same handler also catches the precondition at `:119-120` ("base lfc does not match the block") and any future `assert` in the call stack.

Executed counterexample: `de_stats` patched to return `lfc + 1` on the synthetic deep-ballast block (`test_neggain630.make_case`, 396 carriers, 4 ballast at 995k), then `apply(..., "b_k100")`. `up_pass2` raises `AssertionError: base lfc does not match the block` before it commits a count. `apply` does not raise. It returns a fresh `int64` copy of the input (same shape, same row sums, input array unmodified) and a 19-key dict with `ur_invariant_failed == 1.0`. A bad variant still raises `ValueError` (`parse_spec` is outside the `AssertionError` net).

Failure scenario: a systematic disagreement between `de_stats` and `up_pass2` skips the pass on every target, the h5ad stays the zs block, and the meta says the pass was on. That is a safe matrix. It is the wrong one for the knob, and finding 1 means the day log will not say so.

Fix: catch only an allow-list (the four post-conditions, the `kept gene` U message, and the `rank_mag` helper messages named in the docstring). Re-raise every other `AssertionError` so a programming `assert` still stops the build.

### NOTE — `flips_vs_clean4` counts flips against the already-flipped prior

`experiments/604-identity-prior/def_identity_prior.py:186-189`. At λ 0 with mask `10010001`, `prov["flips_vs_clean4"]` is 0 while 2381 signs differ from the original CLEAN4 psign (`flipped_genes` 2381). The code comment says this. The npz sidecar's `flips_vs_clean4_block` is against the original psign (205 on the top-2200 block for mask `10000000`). The mask string is stored, so this does not collapse artifact shas (confirmed below). Read `flipped_genes` and `flip_octiles`, not `flips_vs_clean4`, to see whether a mask fired.

## Confirmed

### 630 module `7e6b3e52` (`c30d64ce`; file unchanged from that commit to HEAD)

File sha `7e6b3e52e2ad374b7b6e21e92327fcf65771433693eb3846c09c0189f3f2959f`. `git diff 69972aa1 HEAD -- experiments/630-rank-neutral-up/src_rank_up.py` is the wrapper, `DIAG_KEYS`, and `invariant_failed` on the guard and success returns. `up_pass2` is unchanged. `git diff c30d64ce HEAD` on that file is empty.

Bit-identical to `d8c39c3f` (`git show 69972aa1:…`, loaded beside HEAD) when nothing fires, for the count block:

- `apply(..., "b_k100")` on the deep case and the mixed-depth case: blocks equal. Diagnostics differ by exactly one key, `ur_invariant_failed = 0.0`; every other value matches.
- `up_pass2` on both cases with an explicit selection: blocks equal and the stats dicts equal (`n_moved == 1`).

All three paths return the same 19 `DIAG_KEYS` (`src_rank_up.py:257-260`):

- `test_neggain630.py` 14/14. N4 guard, N5 invariant (stubbed `AssertionError`), N6 key-set equality, N7 a real `apply`.
- Independent probe: normal `apply` on the deep case, the guard stub, and the lfc-mismatch path. Each dict's keys equal `DIAG_KEYS`. Guard and mismatch return the input block, `dtype int64`, same shape, same row sums, and a copy (the input object is not aliased and is not mutated).

`test_rank_up630.py` was not run. It calls `evaluate323.load_reference()`, which needs `experiments/323-current-family-validation/work/reference.npz`. That file is not in this checkout (`FileNotFoundError`). `experiments/630-rank-neutral-up/work/blocks/` is absent too.

### 630dd

`uv run python experiments/630-rank-neutral-up/test630dd.py`: 45/45. Live `scripts/def_day.sh` unchanged. `test598.canon_tree` needs no h5ad; it rebuilds the post-bundle tree from `git show HEAD:` plus the 641 → 598 → 612 → 630 patches.

- Commutes with 612dd: both orders produce def_day.sh sha `0c4c97e021c7…` (T2). `bash -n` passes on 630dd alone and on 612dd+630dd (T3).
- Unset `DEF_DAY_DE_UPRN`: no export, no meta check, no `final.env` line, `ALLOW` unchanged (T4). Each of V/E/C/F/A is inside `[ -n "${DEF_DAY_DE_UPRN:-}" ]` or `[ -n "${VCC_DE_UPRN:-}" ]`. DOC is comments. `def_day.sh:66` has already cleared inherited `VCC_*`, so a pre-set `VCC_DE_UPRN` does not reach the C check.
- Regex. The shell pattern in `patch630dd.py:56`, the builder pattern in `patch_builder630.py:59`, and `VARIANT_RE` in `src_rank_up.py:35` are the same string `(s|b)_k(\d{2,3})(?:_q(\d+)r(\d+)p(\d+))?`. For every integer `n` in `0..999`, `0.5 <= n/100 <= 3` matches `50 <= n <= 300` (no mismatches). Sampled variants (`b_k100`, `b_k50`, `b_k300`, `b_k100_q20r3p10`, `b_k100_q0r0p0`, `b_k050`, and the refusals `b_k49`, `b_k301`, `b_k100_q20`, `b_k100_faz`, `B_k100`) agree across shell, builder, and `parse_spec`. Surrounding whitespace is the one difference: the builder strips (`patch_builder630.py:53`) and `parse_spec` strips; the shell does not, so ` b_k100` is refused by def_day and would be accepted by the builder. That fails closed.
- Meta keys. The 630 meta insertion is the `de_mag` anchor inside `"config"` (`scripts/make_submission_pm_cd4.py:6946`, `:7050-7053`; `patch_builder630.py:95-96` writes `de_uprn.variant` and `de_uprn.module_sha256`). The C check reads those two keys (`patch630dd.py:71-72`). T4 passes the matching meta and refuses a wrong variant, a wrong module sha, and a missing `de_uprn`. The pinned sha equals the module file sha above.
- The report cannot abort def_day. Flags are `set -uo pipefail` only. Executed: missing file, malformed file, and a real fail-safe row all return exit 0. When the file has `ur_*` rows, the per-context json and `FAILSAFE_TOTAL` are printed (T5, and the shell probe above).

### 604 octile flip (`0b87b2ba`)

`uv run python experiments/604-identity-prior/test604flip.py`: 21/21. F1 committed v3 shas `f65d12fc` / `6cd84dc2` / `3b1e22f0` (λ 0 / 0.5 / 1) rebuild with no mask and with `--flip 00000000`.

Recomputed from the tracked npz files (`clean4_prior_challenge.npz` sha matches `7e532f84…`, `modules625_v8.npz` sha matches `c24e3304…`, `flex_refs604.npz` via its sidecar):

- 18,533 genes. CLEAN4 covered: 6,346. None have `P == 1/2`.
- Octile values on the file are exactly `{-1, 0, 1, 2, 3, 4, 5, 6, 7}`. Every challenge gene is in the file by name.
- Covered-gene partition: octile 0 = 794, 1 = 793, 2 = 793, 3 = 793, 4 = 793, 5 = 793, 6 = 793, 7 = 794. Sum 6,346. No covered gene has octile `< 0` or `> 7`.
- Uncovered genes: 12,187, all octile `-1`.

`flip_vector` (`def_identity_prior.py:111-113`) requires finite `P` and `octile >= 0`, so `np.clip` never assigns an uncovered gene to octile 0. Executed, mask `10010001`: `fv` selects 2,381 genes, all covered, equal to octiles `{0, 3, 7}`. λ 0 returns `1 - P` and negated psign on those genes, the original prior on the other covered genes, and the original psign with NaN `pdown` on uncovered genes (A, B, and C). The caller's `P` and `psign` arrays are not mutated (`max |ΔP| = 0`). A loop over λ 0, 0.5, 1 on one pair of arrays matches three fresh calls, so `build604.py:32-33` does not flip twice.

None, omitted, and `00000000` take the early returns (`def_identity_prior.py:105-110`) and add no flip keys. F1 shows the artifacts are byte-identical to the committed no-mask v3 files.

Artifact sha includes the mask. `build604.py --lam 0 --flip 10000000` twice: sha `c317f35f7d205cdb…` both times. `--flip 01000000`: sha `f7c2522d1523fef9…`. `pdown_A` and `programme_A` differ. The npz `provenance` array records `flip_octiles` and `flipped_genes` (794 vs 793). Two different masks cannot share an sha: the eight octiles are non-empty and disjoint, so the numeric arrays differ, and the provenance string differs even if an octile were empty.
