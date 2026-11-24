"""431: per-block rc diagnostics from a build run with VCC_RC_AUDIT_DUMP (one JSON line per block).
Answers where the A/B/C MDE dose binds (426a: `mde_lambda` median 0.84 over 900 blocks, H1 1.0): the per-context
distribution of mde_lambda / prior_lambda / mde_overfill_frac and what correlates with a low dose.
usage: uv run python analyse431.py <rc_audit.jsonl>"""
import json
import sys

import numpy as np

rows = [json.loads(l) for l in open(sys.argv[1]) if l.strip()]
print(f"blocks {len(rows)}")
keys = ["mde_lambda", "prior_lambda", "mde_overfill_frac", "mde_fd", "mde_fu", "mde_n_up", "mde_n_dn", "ballast_K",
        "prior_n_flip_target", "clip_frac", "forced", "call_set_unsuppressible", "zero_repaired", "zs_n_moved"]
for ctx in sorted({r["context"] for r in rows}):
    R = [r for r in rows if r["context"] == ctx]
    print(f"\ncontext {ctx}: {len(R)} blocks")
    for k in keys:
        v = np.array([r[k] for r in R if k in r], dtype=float)
        if v.size:
            print(f"  {k:26s} median {np.median(v):10.4f}  q10 {np.quantile(v, .1):10.4f}  min {v.min():10.4f}  "
                  f"frac<1 {np.mean(v < 1):.3f}" if k.endswith("lambda") else
                  f"  {k:26s} median {np.median(v):10.4f}  q10 {np.quantile(v, .1):10.4f}  q90 {np.quantile(v, .9):10.4f}")
lam = np.array([r.get("mde_lambda", np.nan) for r in rows], dtype=float)
ok = np.isfinite(lam)
print("\nSpearman of mde_lambda with block features (all contexts):")
from scipy.stats import spearmanr  # noqa: E402
for k in ["prior_lambda", "mde_overfill_frac", "ballast_K", "prior_n_flip_target", "mde_n_dn", "mde_n_up", "forced"]:
    v = np.array([r.get(k, np.nan) for r in rows], dtype=float)
    m = ok & np.isfinite(v)
    if m.sum() > 10:
        print(f"  {k:24s} {spearmanr(lam[m], v[m]).correlation:+.3f}  (n {m.sum()})")
