# README — 373-blend-rows-patch

Pre-staged exact-anchor extension adding the third `VCC_CTX_A_PRIOR` mode,
`blend_rows`, on top of 366's patch. Status: chain applies cleanly
(366 +155 lines, 373 +116 lines), all 100 tests pass (45 inherited from 366,
55 new). Nothing in `scripts/` was touched; the PI applies per §1 when a
blend_rows artifact clears its gate.

## 1. How the PI applies the patch to the real builder

```bash
cd ~/Code/virtual-cell-challenge-2026
# a. verify every anchor still exists verbatim in current HEAD + record hashes
uv run python experiments/373-blend-rows-patch/patch373.py --check
# b. rebuild the combined patched copy from the live working tree (refuses if
#    already patched, refuses if any anchor occurs != once)
uv run python experiments/373-blend-rows-patch/patch373.py
# c. run the suite against the combined copy
uv run python experiments/373-blend-rows-patch/test373.py   # expect 100 passed, 0 failed
# d. apply to the real builder and inspect the diff (THIS is the only step that writes scripts/)
cp experiments/373-blend-rows-patch/work/patched/make_submission_pm_cd4.py scripts/make_submission_pm_cd4.py
git diff --stat          # expect: scripts/make_submission_pm_cd4.py | ~271 +++++---... (1 file)
git diff scripts/make_submission_pm_cd4.py   # expect 4 hunks marked 366/CTX_A_PRIOR + 7 marked 373/blend_rows
```

Rules for step (d): stage the path explicitly
(`git add scripts/make_submission_pm_cd4.py`, never `git add -A` while workers
run); put the pre/mid/post sha256 printed by `--check` in the commit message
(254 §1 pattern); the flag unset must reproduce the incumbent staging hash
before any flagged build (254 KA-1 pattern).

Known hashes (builder HEAD `1318cd3c…`, 2026-09-10 UTC):
366-only `95913dc0…` (matches 366's record exactly — the chain reproduces 366),
366+373 `c6f366ca…`.

## 2. Building a `blend_rows` artifact `.npz`

Loader contract (structure validated at import; axis/panel checks at apply):
required `targets` (1-D unicode, non-empty, duplicate-free),
`genes` (1-D official symbols, non-empty, duplicate-free),
`matrix` (`len(targets)`×`len(genes)` float64, all finite — NaN refused),
`beta` (scalar in (0,1] — 0 refused), `mode` (`'blend_rows'`); optional
`context` (default `'A'`, A–F), `norm_rule` (default `'fixed_norm_on_covered'`,
only value accepted). `p_down`/`vector`/`tau`/`alpha`/`scale_rule` are refused
in this mode, and `targets`/`matrix`/`beta`/`norm_rule` are refused in 366's two
modes (both directions tested). Targets outside the bundle panel are refused at
apply (the panel is unknowable at import); genes hitting panel-coordinate
positions are refused at apply; a zero artifact row (`‖m‖==0`) raises at apply.

```bash
uv run python - <<'EOF'
import numpy as np, pandas as pd
official = pd.read_csv("data/challenge/gene_names.csv")["gene_name"].astype(str).to_numpy()
panel = set(pd.read_csv("data/challenge/pert_counts.csv")["target_gene"].astype(str))
opos = {g: i for i, g in enumerate(official)}
# example: 40 covered genes, 12 panel targets, beta doses the blend
targets = np.asarray(sorted(["JUN", "FOS", "EGR1", "MYC", "TP53", "JUND", "FOSB",
                             "ATF3", "EGR2", "KLF6", "IER2", "DUSP1"]))
assert set(targets.tolist()) <= panel, "off-panel target — the build will refuse it"
genes = np.asarray(official[1000:1040])
assert not (set(genes.tolist()) & panel), "panel-coordinate support — refused"
rng = np.random.default_rng(0)
mat = rng.normal(0, 1, (len(targets), len(genes)))   # placeholder: real rows from the PI's source
assert np.all(np.isfinite(mat))
np.savez("experiments/373-blend-rows-patch/work/ctx_a_prior_blend.npz",
         targets=targets, genes=genes, matrix=mat, beta=np.asarray(0.35),
         mode=np.asarray("blend_rows"), context=np.asarray("A"),
         norm_rule=np.asarray("fixed_norm_on_covered"))
print("wrote blend artifact", mat.shape)
EOF
```

Build with `VCC_CTX_A_PRIOR=<that file>` plus the 319p base env
(`VCC_NORMREL_EXP=-1.5 VCC_NORMREL_MODE=cap VCC_CD4_PARTIAL_G=0.34 VCC_ALPHA=0.982`).

## 3. One-key config-diff expectation (rule 10)

Flag off: config block byte-identical to the incumbent (no `ctx_a_prior*` keys).
Flag on with a blend_rows artifact: exactly five added keys, nothing
changed/removed:

```
config.ctx_a_prior_mode / ctx_a_prior_sha256 / ctx_a_prior_context / ctx_a_prior_beta / ctx_a_prior_targets
```

(sign_block/add_vector keep 366's `alpha`/`tau` pair instead of `beta`/`targets`;
evaluated from the shipped line in `test373.py::test_config_expression`.)
Postbuild allowance via the existing `--allow-config-diff` flag (254 precedent):

```bash
uv run python scripts/def_postbuild.py --meta <pkg>.h5ad.meta.json --bundle data/challenge \
  --allow-config-diff ctx_a_prior_mode,ctx_a_prior_sha256,ctx_a_prior_context,ctx_a_prior_beta,ctx_a_prior_targets
```

The `[ctx-a-prior]` stdout line prints `rows_touched` for this mode
(`touched_entries` for 366's modes); the `effects.ctx_a` audit carries
`rows_touched/covered_genes/mean_cosine/beta/artifact_sha256`.

## Discrepancies / design notes (all deliberate, all tested)

1. Targets-vs-panel validation is split by necessity: duplicates/NaN/shape/beta
   at import, panel-membership at apply (the loader cannot see `spec.perts`).
2. `beta=0` artifacts are refused at load; the "β=0 ⇒ identity" test uses a
   hand-built prior dict to prove the arithmetic path is the identity there.
3. A design row exactly zero on the covered genes is skipped (the fixed-norm
   blend is the identity there) and excluded from `mean_cosine`; it still counts
   in `artifact_targets` but not in `rows_touched`.
4. `build_ctx_a_context_effects` gains `panel_targets=None`; old callers are
   unaffected (`kw.pop`, old path never sees the kwarg — proven by the
   re-run 45/45).
