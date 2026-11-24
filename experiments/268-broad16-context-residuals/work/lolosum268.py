"""Pooled LOLO table with paired target bootstraps from lolo268_K<K>.csv (targets are the resampling unit; each target carries all 12 folds)."""
import sys, numpy as np, pandas as pd
K = sys.argv[1]; df = pd.read_csv(f"lolo268_K{K}.csv"); rng = np.random.default_rng(268)
arms = [a for a in dict.fromkeys(c.split("|")[0] for c in df.columns if "|" in c)]
tg = df.groupby("target").mean(numeric_only=True); fold = df.groupby("held").mean(numeric_only=True)
def boot(v, n=4000):
    v = np.asarray(v); idx = rng.integers(0, len(v), (n, len(v))); return v[idx].mean(1)
print(f"K={K}: targets {df.target.nunique()}, folds {df.held.nunique()}; ceiling cos(rA,rB) median {df.ceiling_cos_AB.median():.3f} mean {df.ceiling_cos_AB.mean():.3f}; "
      f"target-only cos(full) mean {df.base_cos_full.mean():.3f}; target-only retrieval mean {df.base_retrieval.mean():.3f}; |residual|^2/|response|^2 median {df.res_norm_share.median():.2f}")
print("\n| arm | cos(pred residual, held residual) mean [95% CI] | folds > 0 | Δcos(full response) mean [95% CI] | folds > 0 | Δretrieval mean [95% CI] | folds > 0 |\n|---|---|---|---|---|---|---|")
for a in arms:
    cells = []
    for met in ("cos_res", "dcos_full", "dretrieval"):
        col = f"{a}|{met}"; b = boot(tg[col]); cells.append(f"{tg[col].mean():+.4f} [{np.quantile(b, .025):+.4f}, {np.quantile(b, .975):+.4f}]"); cells.append(f"{int((fold[col] > 0).sum())}/12")
    print(f"| `{a}` | " + " | ".join(cells) + " |")
print("\nPaired differences (per-target means, bootstrap over targets):\n| comparison | metric | Δ mean [95% CI] | P(Δ ≤ 0) |\n|---|---|---|---|")
main = ["b_recipient_offset", "c1_target_x_basal_ridge", "c2_target_x_basal_programmes", "nn_basal_nearest_line"]
for a in main:
    for ctrl in (f"{a}__WRONG_recipient", f"{a}__SHUFFLED_labels", "a_target_only"):
        if ctrl not in arms: continue
        for met in ("cos_res", "dcos_full", "dretrieval"):
            d = tg[f"{a}|{met}"] - tg[f"{ctrl}|{met}"]; b = boot(d)
            print(f"| `{a}` − `{ctrl}` | {met} | {d.mean():+.4f} [{np.quantile(b, .025):+.4f}, {np.quantile(b, .975):+.4f}] | {(b <= 0).mean():.3f} |")
