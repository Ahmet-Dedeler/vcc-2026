# README — 366-ctx-a-prior-patch

Pre-staged exact-anchor patch adding `VCC_CTX_A_PRIOR` to a copy of the builder.
Status: patch applies cleanly (+155 lines), all 45 tests pass. Nothing in `scripts/`
was touched; the PI applies the patch to the real builder when the prior clears its gate.

## 1. How the PI applies the patch to the real builder

```bash
cd ~/Code/virtual-cell-challenge-2026
# a. verify the anchors still exist verbatim in current HEAD + record pre/post hashes
uv run python experiments/366-ctx-a-prior-patch/patch366.py --check
# b. rebuild the patched copy from the live working tree (refuses if already patched,
#    refuses if any anchor occurs != once)
uv run python experiments/366-ctx-a-prior-patch/patch366.py
# c. run the suite against the patched copy
uv run python experiments/366-ctx-a-prior-patch/test366.py   # expect 45 passed, 0 failed
# d. apply to the real builder and inspect the diff (THIS is the only step that writes scripts/)
cp experiments/366-ctx-a-prior-patch/work/patched/make_submission_pm_cd4.py scripts/make_submission_pm_cd4.py
git diff --stat          # expect: scripts/make_submission_pm_cd4.py | ~155 +++++---... (1 file)
git diff scripts/make_submission_pm_cd4.py   # expect exactly 4 hunks, all marked 366/CTX_A_PRIOR
```

Rules for step (d): stage the path explicitly (`git add scripts/make_submission_pm_cd4.py`,
never `git add -A` while workers run); put the pre/post sha256 printed by `--check` in the
commit message (254 §1 pattern); the flag unset must reproduce the incumbent staging hash
before any flagged build (254 KA-1 pattern).

## 2. Building the artifact `.npz`

Loader contract (validated at import; axis checks at apply): keys `genes` (str, subset of
the official 18,533 axis by NAME), `mode` (`sign_block`|`add_vector`), `vector` (float,
same shape, finite); optional `p_down` (finite ∈ [0,1], sign_block only),
`tau` (default 0.0 ∈ [0,0.5]), `alpha` (default 1.0, finite ≠ 0),
`scale_rule` (default `'gate_median'`, only value accepted), `context` (default `'A'`).
Zero support on the 300 panel-coordinate columns or the build raises.

### sign_block from 370 (`jurkat_block_prior.npz`: bank axis, `P_down`)

```bash
uv run python - <<'EOF'
import numpy as np, pandas as pd
z = np.load("experiments/370-jurkat-block-sign-prior/work/jurkat_block_prior.npz", allow_pickle=True)
official = pd.read_csv("data/challenge/gene_names.csv")["gene_name"].astype(str).to_numpy()
bank_genes = z["genes"].astype(str)
keep = np.array([g in set(official) for g in bank_genes])
assert keep.sum() == 8284, keep.sum()   # 365 KA-1: bank ∩ official axis
np.savez("experiments/366-ctx-a-prior-patch/work/ctx_a_prior_sign.npz",
         genes=bank_genes[keep], mode=np.asarray("sign_block"),
         vector=np.zeros(keep.sum()),            # required but unused when p_down present
         p_down=np.asarray(z["P_down"], dtype=np.float64)[keep],
         tau=np.asarray(0.0),                    # PI doses (363 owns the value)
         alpha=np.asarray(1.0),
         scale_rule=np.asarray("gate_median"), context=np.asarray("A"))
print("wrote ctx_a_prior_sign.npz", keep.sum(), "genes")
EOF
```

### add_vector from 350 (`jurkat_common.npz`: already official axis)

```bash
uv run python - <<'EOF'
import numpy as np, pandas as pd
z = np.load("experiments/350-same-line-priors-data/work/jurkat_common.npz", allow_pickle=True)
official = pd.read_csv("data/challenge/gene_names.csv")["gene_name"].astype(str).to_numpy()
panel = set(pd.read_csv("data/challenge/pert_counts.csv")["target_gene"].astype(str))
assert list(z["genes"].astype(str)) == list(official), "350 axis drift — stop, do not ship"
vec = np.asarray(z["vector"], dtype=np.float64).copy()
vec[[i for i, g in enumerate(official) if g in panel]] = 0.0   # loader refuses panel support
assert np.median(np.abs(vec[vec != 0.0])) > 0
np.savez("experiments/366-ctx-a-prior-patch/work/ctx_a_prior_add.npz",
         genes=official, mode=np.asarray("add_vector"), vector=vec,
         alpha=np.asarray(1.0),                   # PI doses (362 owns the value)
         scale_rule=np.asarray("gate_median"), context=np.asarray("A"))
print("wrote ctx_a_prior_add.npz, support", int((vec != 0).sum()))
EOF
```

Build with `VCC_CTX_A_PRIOR=<that file>` plus the 319p base env
(`VCC_NORMREL_EXP=-1.5 VCC_NORMREL_MODE=cap VCC_CD4_PARTIAL_G=0.34 VCC_ALPHA=0.982`).

## 3. One-key config-diff expectation (rule 10)

Flag off: config block byte-identical to the incumbent (no `ctx_a_prior*` keys anywhere).
Flag on: exactly five added keys, nothing changed/removed:

```
config.ctx_a_prior_mode / ctx_a_prior_sha256 / ctx_a_prior_alpha / ctx_a_prior_tau / ctx_a_prior_context
```

plus their consequence `effects.ctx_a.*`. Postbuild and the driver need NO source change;
pass the keys explicitly (254 precedent):

```bash
uv run python scripts/def_postbuild.py --meta <pkg>.h5ad.meta.json --bundle data/challenge \
  --allow-config-diff ctx_a_prior_mode,ctx_a_prior_sha256,ctx_a_prior_alpha,ctx_a_prior_tau,ctx_a_prior_context
```

Note for the driver (`scripts/def_pipeline.py --allow-arm-drift` is for dropped shipped arms,
not this): `assert_build_carries_arms` has no `ctx_a_prior` row, so a flagged build passes
it silently — the PI should add a row (239 §2.1 class) when this graduates from pre-staged
patch to shipped arm. Until then the `[ctx-a-prior]` stdout line + `effects.ctx_a` audit are
the visibility.

## Discrepancies vs DESIGN.md 365 (all deliberate, all tested)

1. Insertion anchor is the knee tail (`del frames`/`gc.collect()` → `del donors, ...`),
   not L5325 `context_effects = None`: inserting there runs before the knee control frames
   load and would be overwritten by the knee rebuild. Knee+prior is refused at import, so
   both orders coincide when off (proven by the inert test).
2. Artifact format follows the 366 mandate (`genes/mode/vector[/p_down]/tau/alpha/
   scale_rule/context`), not 365 §2's `targets/matrix/sign/measured` layout.
3. Config keys follow the mandate (`ctx_a_prior_{mode,sha256,alpha,tau,context}`), not 365
   §2's (`ctx_a_prior/mode/dose/targets`); `replace` mode is not implemented (365 C1:
   vacuous on A/B/C — the loader's mode gate enforces this).
4. Axis-by-NAME mapping and panel-support checks run at apply time (official axis unknown
   at import); everything axis-independent is validated at import. `sign` for vector-only
   artifacts is `sign(vector)`, 0 = untouched (365 §2a); the else-`down` default applies to
   the `p_down` tau-rule exactly as mandated. Artifact-absent genes are untouched.
5. `add_vector` scale `s` is per-row (mandate: "the row's universe genes"), not 365 §2's
   pooled median; universe = `b_ctrl > 0` genes. No `def_postbuild.py`/`def_pipeline.py`
   patch: the existing `--allow-config-diff` mechanism covers the five keys.
