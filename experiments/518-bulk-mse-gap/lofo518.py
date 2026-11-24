"""518 — how much of the per-target magnitude pot is capturable with leak-free features: a leave-one-fold-out fit.

Multiplier on the production key: m_p = exp(theta . z_p), z = centred percentile ranks (in [-0.5, 0.5]; missing -> 0) of
    kv1      the LODO key vector itself (its information beyond the (0.5, 3) clip)
    cd4      CD4 T-cell DE count (log1p nDE, mean over the conditions carried; Marson GWCD4i DESeq2, adj p < 0.05)
    pre      our pre-key row energy (the stack composition's non-panel norm)
    srcde    the fold sources' own DE-count proxy (feat518 logde_mean)
    srccos   the fold sources' mean pairwise row cosine (feat518 src_cos)
theta fitted to maximise the mean relative cos_agg gain on the TRAINING folds (headroom540's first-order route on arms540's
saved per-target arrays; ridge 1e-3), scored on the held-out fold. Cancer folds rotate; kolf2 is scored with theta fitted on
the three cancer folds. Also: each feature alone (fitted the same way) and the in-sample optimum (an upper bound).

    uv run python experiments/518-bulk-mse-gap/lofo518.py
"""
from __future__ import annotations

import sys

sys.dont_write_bytecode = True
import json  # noqa: E402
import warnings  # noqa: E402
from pathlib import Path  # noqa: E402

import numpy as np  # noqa: E402
from scipy.optimize import minimize  # noqa: E402

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from cd4feat518 import W540, WORK, BASE, FOLDS, CANCER, cd4_features  # noqa: E402

warnings.filterwarnings("ignore")
FEATS = ("kv1", "cd4", "pre", "srcde", "srccos")


def crank(v):
    v = np.asarray(v, float)
    ok = np.isfinite(v)
    out = np.zeros(v.size)
    if ok.sum() > 1:
        r = np.argsort(np.argsort(v[ok]))
        out[ok] = r / (ok.sum() - 1) - 0.5
    return out


def load():
    fa = np.load(WORK / "feat518_arrays.npz")
    data = {}
    for f in FOLDS:
        z = np.load(W540 / f"arrays540-{f}.npz")
        x, a, d = z[f"{BASE}__x"], z[f"{BASE}__a"], z["d"]
        dz = np.load(WORK / f"dump518-{f}.npz")
        targets = [str(t) for t in dz["targets"]]
        cd4 = cd4_features(targets)["lognde_mean"]
        pre = np.linalg.norm(dz["fST"][:, np.setdiff1d(np.arange(dz["fST"].shape[1]), dz["panel_positions"])], axis=1)
        Z = np.column_stack([crank(np.log(np.maximum(dz["kv1"], 1e-6))), crank(cd4), crank(pre),
                             crank(fa[f"{f}__logde_mean"]), crank(fa[f"{f}__src_cos"])])
        base = float(np.mean([x[j].sum() / np.sqrt(a[j].sum() * d.sum()) for j in range(3)]))
        data[f] = (x, a, d, Z, base)
    return data


def gain(theta, D, cols):
    x, a, d, Z, base = D
    m = np.exp(Z[:, cols] @ theta)
    c = np.mean([(m * x[j]).sum() / np.sqrt((m * m * a[j]).sum() * d.sum()) for j in range(3)])
    return c / base - 1


def fit(data, train, cols):
    obj = lambda th: -np.mean([gain(th, data[f], cols) for f in train]) + 1e-3 * float(th @ th)
    best = None
    for s in range(4):
        th0 = np.random.default_rng(s).normal(0, 0.2, len(cols))
        r = minimize(obj, th0, method="Nelder-Mead", options=dict(maxiter=4000, xatol=1e-5, fatol=1e-9))
        if best is None or r.fun < best.fun:
            best = r
    return best.x


def main():
    data = load()
    out = dict(features=FEATS, sets={})
    sets = {"all": list(range(len(FEATS)))} | {f"only_{n}": [i] for i, n in enumerate(FEATS)} | {"kv1+cd4": [0, 1],
                                                                                               "cd4+pre": [1, 2], "kv1+cd4+pre": [0, 1, 2]}
    for nm, cols in sets.items():
        rec = dict(heldout={}, theta={})
        for h in CANCER:
            tr = [f for f in CANCER if f != h]
            th = fit(data, tr, cols)
            rec["heldout"][h] = gain(th, data[h], cols)
            rec["theta"][h] = th.tolist()
        th = fit(data, list(CANCER), cols)
        rec["heldout"]["kolf2 (fit on 3 cancer)"] = gain(th, data["kolf2"], cols)
        rec["theta"]["cancer3"] = th.tolist()
        rec["insample_cancer3"] = float(np.mean([gain(th, data[f], cols) for f in CANCER]))
        rec["lofo_cancer_mean"] = float(np.mean([rec["heldout"][h] for h in CANCER]))
        out["sets"][nm] = rec
        print(f"{nm:14s} LOFO cancer mean {100 * rec['lofo_cancer_mean']:+.2f} % (" +
              " ".join(f"{100 * rec['heldout'][h]:+.2f}" for h in CANCER) + f")  kolf2 {100 * rec['heldout']['kolf2 (fit on 3 cancer)']:+.2f} %"
              f"  in-sample {100 * rec['insample_cancer3']:+.2f} %  theta(cancer3) {np.round(th, 3).tolist()}")
    (WORK / "lofo518.json").write_text(json.dumps(out, indent=1, default=float))


if __name__ == "__main__":
    main()
