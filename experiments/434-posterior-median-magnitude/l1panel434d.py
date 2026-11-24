"""434d: does 434b's gain survive PANEL-AVAILABLE donors only (k562_gwps, orion_hct116, orion_hek293t; 430 section 0.3)? The A/B/C panel has no k562_essential coverage (0/300), and
430 traced H1's donor-sign advantage to essential-gene donors. Here: drop k562_essential from the features entirely,
then report the out-of-fold nmae change on (i) all 58 gated targets and (ii) the targets k562_essential does NOT cover
(the panel-like subset), for the full feature set and the donors-only subset. Same model as 434b's chosen config."""
from pathlib import Path

import numpy as np
from sklearn.ensemble import HistGradientBoostingRegressor
from sklearn.model_selection import GroupKFold

ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
z = np.load(HERE / "work/allcells434.npz", allow_pickle=True)
L, Q, D, T = z["L"].astype(float), z["Q"].astype(float), z["D"].astype(float), z["T"].astype(float)
donors = [str(x) for x in z["donors"]]
ESS = [donors.index(d) for d in ("k562_essential", "rpe1", "nadig_jurkat")]   # essential-gene screens: 0/300 panel coverage (430)
ess_cov = np.isfinite(D[ESS]).any(2).any(0)              # targets any essential-screen donor covers
D = np.delete(D, ESS, axis=0)                           # panel-available donors only: k562_gwps + two Orion lines
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
den = np.array([np.abs(l[K == k]).sum() if (K == k).any() else 1.0 for k in range(K.max() + 1)])
w = 1.0 / den[K]
mde = lh + s * np.where(cnt[I] >= 0.3, np.minimum(0.1 / np.sqrt(np.maximum(cnt[I], 1e-3)), 0.1), 0.0)

def nmae(pred, sel=None):
    ks = [k for k in np.unique(K) if sel is None or sel[k]]
    return float(np.mean([np.abs(pred[K == k] - l[K == k]).sum() / np.abs(l[K == k]).sum() for k in ks]))

def oof(cols, train_sel=None):
    pred = np.zeros(len(y))
    for tr, te in GroupKFold(5).split(X, y, K):
        if train_sel is not None:
            tr = tr[train_sel[K[tr]]]
        m = HistGradientBoostingRegressor(loss="quantile", quantile=0.5, max_iter=400, learning_rate=0.03, max_leaf_nodes=15,
                                          min_samples_leaf=500, l2_regularization=1.0, random_state=0)
        m.fit(X[tr][:, cols], y[tr], sample_weight=w[tr])
        pred[te] = m.predict(X[te][:, cols])
    return pred

panel_like = ~ess_cov
print(f"gated targets {len(np.unique(K))}; covered by k562_essential {int(ess_cov[np.unique(K)].sum())}; panel-like {int(panel_like[np.unique(K)].sum())}")
for nm, cols in (("full (panel donors)", list(range(X.shape[1]))), ("donors only (+cnt)", [2, 3, 4, 5, 6, 7, 8, 12])):
    for tn, tsel in (("train all", None), ("train panel-like only", panel_like)):
        p = oof(cols, tsel)
        print(f"{nm:26s} {tn:22s}: all {nmae(s * p) - nmae(lh):+.5f} | panel-like {nmae(s * p, panel_like) - nmae(lh, panel_like):+.5f} "
              f"| ess-covered {nmae(s * p, ess_cov) - nmae(lh, ess_cov):+.5f}")
print(f"MDE: all {nmae(mde) - nmae(lh):+.5f} | panel-like {nmae(mde, panel_like) - nmae(lh, panel_like):+.5f} | ess-covered {nmae(mde, ess_cov) - nmae(lh, ess_cov):+.5f}")
print(f"base nmae: all {nmae(lh):.5f} | panel-like {nmae(lh, panel_like):.5f} | ess-covered {nmae(lh, ess_cov):.5f}")
