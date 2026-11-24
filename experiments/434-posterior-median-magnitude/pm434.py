"""434 analytic (H1): calibrated per-cell sign confidence -> posterior-median magnitude for nmae.
nmae reads only the reference's significant genes, so the L1-optimal prediction for a cell is the median of l given
the features AND significance. With a = P(our sign is right | features, significant) and F the distribution of |l| on
the gate (per expression bin), that median is s * F^-1(1 - 1/(2a)) when a > 1/2 (0 otherwise; the sign flips when
a < 1/2). a is fitted out-of-fold over targets (GroupKFold 5) with gradient boosting on oracle-free features.
Compares: base (426 base_s0_bps), MDE add (426's c 0.1 cap 0.1 min 0.3), the rule (replace / raise-only / add),
and the bounds (our signs with oracle magnitudes; in-sample a)."""
from pathlib import Path

import numpy as np
from sklearn.ensemble import HistGradientBoostingClassifier
from sklearn.model_selection import GroupKFold

ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
z = np.load(HERE / "work/allcells434.npz", allow_pickle=True)
L, Q, D, T = z["L"].astype(float), z["Q"].astype(float), z["D"].astype(float), z["T"].astype(float)
r = np.load(ROOT / "experiments/323-current-family-validation/work/reference.npz", allow_pickle=True)
keep = r["keep"].astype(bool)
ref_lin = np.asarray(r["ref_lin"], float)
ref_lin = ref_lin[keep] if ref_lin.shape[0] == keep.shape[0] else ref_lin
cnt = ref_lin * np.median(np.asarray(r["libraries"], float)) / 1e6
print("cnt", cnt.shape, np.quantile(cnt, [.01, .5, .99]))
G = np.isfinite(T)
K, I = np.nonzero(G)
l = T[K, I]; lh = L[K, I]; s = np.sign(lh)
# cells where we deliver exactly 0: take the donor-mean sign (they are rare)
dm = np.nanmean(D[:, K, I], 0)
s = np.where(s == 0, np.sign(np.nan_to_num(dm)), s); s = np.where(s == 0, 1.0, s)
y = (np.sign(l) == s).astype(int)
nd = np.isfinite(D[:, K, I]).sum(0)
agree = np.nansum(np.sign(D[:, K, I]) * s[None], 0)
supp = np.nanmean(D[:, K, I] * s[None], 0); supp = np.nan_to_num(supp)
mabs = np.nan_to_num(np.nanmean(np.abs(D[:, K, I]), 0))
rown = np.abs(L).sum(1); ncall = (Q < 0.05).sum(1)
drown = np.nan_to_num(np.nanmean(np.nansum(np.abs(D), 2), 0))
gvar = np.nanstd(np.nan_to_num(D), axis=(0, 1))          # gene responsiveness across donor rows
X = np.column_stack([np.log10(np.abs(lh) + 1e-4), np.clip(-np.log10(Q[K, I] + 1e-300), 0, 50), np.log10(cnt[I] + 1e-4),
                     nd, agree, supp, mabs, np.log10(rown[K] + 1e-3), np.log10(ncall[K] + 1), np.log10(drown[K] + 1e-3),
                     np.log10(gvar[I] + 1e-5)])
names = ["log|ours|", "-log10 q", "log cnt", "n donors", "agree", "support", "mean|d|", "log rownorm", "log ncalls",
         "log donor rownorm", "log gene var"]
print("gate cells", len(y), "our sign acc", y.mean().round(4), "targets", len(np.unique(K)))
ahat = np.zeros(len(y)); gkf = GroupKFold(5)
for tr, te in gkf.split(X, y, K):
    m = HistGradientBoostingClassifier(max_iter=300, learning_rate=0.05, max_leaf_nodes=15, min_samples_leaf=200,
                                       l2_regularization=1.0, random_state=0).fit(X[tr], y[tr])
    ahat[te] = m.predict_proba(X[te])[:, 1]
# calibration
qs = np.quantile(ahat, np.linspace(0, 1, 11))
print("\ncalibration (deciles of a-hat): mean a-hat | observed acc | n | share of sum|l|")
tot = np.abs(l).sum()
for a0, a1 in zip(qs[:-1], qs[1:]):
    b = (ahat >= a0) & (ahat <= a1)
    print(f"  {ahat[b].mean():.3f} | {y[b].mean():.3f} | {b.sum():6d} | {np.abs(l[b]).sum() / tot:.3f}")
# |l| distribution on the gate per expression bin (fitted out of fold too)
bins = np.quantile(np.log10(cnt[I] + 1e-4), np.linspace(0, 1, 9)); bins[-1] += 1e-9
eb = np.clip(np.digitize(np.log10(cnt[I] + 1e-4), bins) - 1, 0, 7)

def quant_rule(a, fold_tr):
    out = np.zeros(len(a))
    for b in range(8):
        tr = fold_tr & (eb == b)
        absl = np.sort(np.abs(l[tr]))
        sel = eb == b
        p = np.clip(1 - 0.5 / np.maximum(a[sel], 0.5), 0, 1)       # quantile level for a >= 1/2
        out[sel] = np.where(a[sel] > 0.5, np.quantile(absl, p) if absl.size else 0.0, 0.0)
    return out

m_rule = np.zeros(len(y)); m_flip = np.zeros(len(y))
for tr_idx, te_idx in gkf.split(X, y, K):
    trm = np.zeros(len(y), bool); trm[tr_idx] = True
    mr = quant_rule(ahat, trm); mf = quant_rule(1 - ahat, trm)
    m_rule[te_idx] = mr[te_idx]; m_flip[te_idx] = mf[te_idx]

def nmae(pred):
    return float(np.mean([np.abs(pred[K == k] - l[K == k]).sum() / np.abs(l[K == k]).sum() for k in np.unique(K)]))

base = nmae(lh)
mde = lh + s * np.where(cnt[I] >= 0.3, np.minimum(0.1 / np.sqrt(np.maximum(cnt[I], 1e-3)), 0.1), 0.0)
rows = [("base (ours)", lh), ("MDE add c0.1 cap0.1 min0.3", mde),
        ("rule: replace, s*m", s * m_rule),
        ("rule: replace with flips", np.where(ahat >= 0.5, s * m_rule, -s * m_flip)),
        ("rule: raise-only max(|ours|, m) on our sign", s * np.maximum(np.abs(lh) * (np.sign(lh) == s), m_rule)),
        ("rule: MDE add, then raise to m", s * np.maximum(np.abs(mde), m_rule)),
        ("BOUND our signs, oracle |l|", s * np.abs(l)),
        ("BOUND oracle signs, our |l| + MDE", np.sign(l) * np.abs(mde))]
print(f"\nnmae on h1-test gates (58 targets): base {base:.5f}")
for name, p in rows:
    print(f"  {name:48s} {nmae(p):.5f}  d {nmae(p) - base:+.5f}")
m = HistGradientBoostingClassifier(max_iter=300, learning_rate=0.05, max_leaf_nodes=15, min_samples_leaf=200,
                                   l2_regularization=1.0, random_state=0).fit(X, y)
from sklearn.inspection import permutation_importance  # noqa: E402
pi = permutation_importance(m, X, y, n_repeats=3, random_state=0, scoring="roc_auc")
print("\npermutation importance (AUC drop, in-sample):", {n: round(v, 4) for n, v in zip(names, pi.importances_mean)})
from sklearn.metrics import roc_auc_score  # noqa: E402
print("out-of-fold AUC", round(roc_auc_score(y, ahat), 4))
np.savez_compressed(HERE / "work/pm434_oof.npz", ahat=ahat, m_rule=m_rule, K=K, I=I)
