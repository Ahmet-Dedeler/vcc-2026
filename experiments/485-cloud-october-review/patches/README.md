# 485 patches

Apply from the repo root with `patch -p1`. Dry-run checked on this branch: all four apply clean.

**`def_pipeline.py.patch` and `def_day.sh.patch` land together.** The shell starts requiring the dump stage to exit 0. The driver on master exits 1 after every successful dump, because `assert_build_carries_arms` looks for a meta file the dump never writes (`make_submission_pm_cd4.py` saves the npz and `SystemExit(0)`s before `write_submission`). Applying only the shell patch makes every dump die. Applying only the driver patch leaves the shell ignoring that exit code.

`def_day_helpers.py.patch` is required by the shell patch (`inband`, `best_safe`, `meta_matches`, `past_cliff` exit 2).

`relevel_iter.py.patch` is what makes a flat or repeated-factor secant exit 2 with no `VCC_BULK_CTX_SCALE` line, so the shell's existing fallback can run. The 452 rehearsal's exponents (1.89–2.73) are unchanged: the patched and unpatched scripts print the same stdout on that input.

The files under `repro/*.new` are the patched sources these diffs were taken from. Regenerating:

```bash
diff -u --label a/scripts/def_day.sh --label b/scripts/def_day.sh \
  scripts/def_day.sh experiments/485-cloud-october-review/repro/def_day.sh.new \
  > experiments/485-cloud-october-review/patches/def_day.sh.patch
```

Same pattern for `scripts/def_day_helpers.py`, `scripts/def_pipeline.py`, and `experiments/448-def-rehearsal-q175/relevel_iter.py`.
