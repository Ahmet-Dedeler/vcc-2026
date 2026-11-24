"""434b analytic (H1): direct L1 regression of the truth on oracle-free features, on the gate, out-of-fold by target.
The L1-optimal prediction of a gate cell is the conditional MEDIAN of l given features and significance, so fit it with
the quantile (0.5) loss, oriented by our sign (y = s * l), weights 1/sum_p|l| so the loss is nmae's per-target mean.
This bounds what our existing signals (our delivered lfc/q, six donors, expression, target strength) can buy on nmae."""
from pathlib import Path

import numpy as np
from sklearn.ensemble import HistGradientBoostingRegressor
from sklearn.model_selection import GroupKFold

ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
z = np.load(HERE / "work/allcells434.npz", allow_pickle=True)
L, Q, D, T = z["L"].astype(float), z["Q"].astype(float), z["D"].astype(float), z["T"].astype(float)
r = np.load(ROOT / "experiments/323-current-family-validation/work/reference.npz", allow_pickle=True)
keep = r["keep"].astype(bool); ref_lin = np.asarray(r["ref_lin"], float)
ref_lin = ref_lin[keep] if ref_lin.shape[0] == keep.shape[0] else ref_lin
cnt = ref_lin * np.median(np.asarray(r["libraries"], float)) / 1e6
G = np.isfinite(T); K, I = np.nonzero(G)
l = T[K, I]; lh = L[K, I]; s = np.sign(lh)
with np.errstate(all="ignore"):
    dm = np.nanmean(D[:, K, I], 0)
    s = np.where(s == 0, np.sign(np.nan_to_num(dm)), s); s = np.where(s == 0, 1.0, s)
    nd = np.isfinite(D[:, K, I]).sum(0)
    agree = np.nansum(np.sign(D[:, K, I]) * s[None], 0)
    supp = np.nan_to_num(np.nanmean(D[:, K, I] * s[None], 0))
    smin = np.nan_to_num(np.nanmin(D[:, K, I] * s[None], 0)); smax = np.nan_to_num(np.nanmax(D[:, K, I] * s[None], 0))
    mabs = np.nan_to_num(np.nanmean(np.abs(D[:, K, I]), 0))
    rown = np.abs(L).sum(1); ncall = (Q < 0.05).sum(1)
    drown = np.nan_to_num(np.nanmean(np.nansum(np.abs(D), 2), 0))
    gvar = np.nanstd(np.nan_to_num(D), axis=(0, 1))
X = np.column_stack([np.abs(lh), np.clip(-np.log10(Q[K, I] + 1e-300), 0, 50), np.log10(cnt[I] + 1e-4), nd, agree, supp, smin,
                     smax, mabs, np.log10(rown[K] + 1e-3), np.log10(ncall[K] + 1), np.log10(drown[K] + 1e-3), np.log10(gvar[I] + 1e-5)])
y = s * l
den = np.array([np.abs(l[K == k]).sum() for k in range(K.max() + 1)] if True else [])
w = 1.0 / np.maximum(den[K], 1e-9)

def nmae(pred):
    return float(np.mean([np.abs(pred[K == k] - l[K == k]).sum() / np.abs(l[K == k]).sum() for k in np.unique(K)]))

base = nmae(lh)
mde = lh + s * np.where(cnt[I] >= 0.3, np.minimum(0.1 / np.sqrt(np.maximum(cnt[I], 1e-3)), 0.1), 0.0)
print(f"base {base:.5f}  MDE {nmae(mde) - base:+.5f}")
for lr, leaves, msl in ((0.05, 15, 200), (0.05, 31, 100), (0.03, 15, 500)):
    for seed in (0,):
        pred = np.zeros(len(y))
        for tr, te in GroupKFold(5).split(X, y, K):
            m = HistGradientBoostingRegressor(loss="quantile", quantile=0.5, max_iter=400, learning_rate=lr, max_leaf_nodes=leaves,
                                              min_samples_leaf=msl, l2_regularization=1.0, random_state=seed)
            m.fit(X[tr], y[tr], sample_weight=w[tr])
            pred[te] = m.predict(X[te])
        lhat = s * pred
        print(f"L1-GBM lr {lr} leaves {leaves} msl {msl}: nmae {nmae(lhat):.5f} d {nmae(lhat) - base:+.5f}; "
              f"sign of pred vs ours flipped on {np.mean(pred < 0):.3f}; median |pred| {np.median(np.abs(pred)):.3f} (|l| {np.median(np.abs(l)):.3f})")
# per-target breakdown for the first config

# ---- decomposition for the chosen config (lr 0.03, leaves 15, msl 500): magnitude vs flips, feature subsets, per target
names = ["|ours|", "-log10 q", "log cnt", "n donors", "agree", "support", "smin", "smax", "mean|d|", "log rownorm",
         "log ncalls", "log donor rownorm", "log gene var"]

def oof(cols, mono=None, seed=0):
    pred = np.zeros(len(y))
    for tr, te in GroupKFold(5).split(X, y, K):
        m = HistGradientBoostingRegressor(loss="quantile", quantile=0.5, max_iter=400, learning_rate=0.03, max_leaf_nodes=15,
                                          min_samples_leaf=500, l2_regularization=1.0, random_state=seed,
                                          monotonic_cst=mono)
        m.fit(X[tr][:, cols], y[tr], sample_weight=w[tr])
        pred[te] = m.predict(X[te][:, cols])
    return pred

allc = list(range(X.shape[1]))
p_all = oof(allc)
print(f"\nfull: {nmae(s * p_all) - base:+.5f}; our signs only (no flips, max(pred,0)): {nmae(s * np.maximum(p_all, 0)) - base:+.5f}; "
      f"flips only at our |l|+MDE magnitude: {nmae(np.sign(p_all + 1e-12) * np.abs(mde)) - base:+.5f}")
for drop in allc:
    cols = [c for c in allc if c != drop]
    print(f"  drop {names[drop]:18s}: {nmae(s * oof(cols)) - base:+.5f}")
sub = {"ours only (|ours|, q, cnt)": [0, 1, 2], "donors only (+cnt)": [2, 3, 4, 5, 6, 7, 8, 12],
       "ours + donors, no target-level": [0, 1, 2, 3, 4, 5, 6, 7, 8, 12]}
for k_, cols in sub.items():
    print(f"  subset {k_:32s}: {nmae(s * oof(cols)) - base:+.5f}")
pt = np.array([np.abs(s[K == k] * p_all[K == k] - l[K == k]).sum() / np.abs(l[K == k]).sum()
               - np.abs(lh[K == k] - l[K == k]).sum() / np.abs(l[K == k]).sum() for k in np.unique(K)])
print(f"per target: improved {np.mean(pt < 0):.2f} of {len(pt)}; median d {np.median(pt):+.4f}; worst {pt.max():+.4f}; best {pt.min():+.4f}")
