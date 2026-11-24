"""518 — a per-target correction learned from PLOLO itself: the residual ideal multiplier u_p of OTHER folds as a feature.

u_p(f) = (sum_frames x_p / sum_frames a_p) / (sum x / sum a) on fold f's 446bj2 bulk (cap 3) or cap 6 (arms540's arrays).
Feature for a held-out fold h: z_p = centred rank of the mean over the training folds of log(max(u_p, 0.05)) (targets absent
from every training fold -> 0). Multiplier exp(theta z) on the key; theta fitted on the training folds (each training fold
scored with the feature built from the OTHER training folds, so the fit never sees a fold's own u), scored held out.
Leakage note: a cancer fold's u is built with the other cancer lines as SOURCES, so the held-out cancer line's data enters
the feature through the prediction (a bias against the feature, not for it). kolf2 is the clean test: its data never
enters a cancer fold; its feature is built from the three cancer folds.

    uv run python experiments/518-bulk-mse-gap/xfold518.py
"""
from __future__ import annotations

import sys

sys.dont_write_bytecode = True
import json  # noqa: E402
from pathlib import Path  # noqa: E402

import numpy as np  # noqa: E402
from scipy.optimize import minimize_scalar  # noqa: E402

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from lofo518 import crank  # noqa: E402

ROOT = HERE.parents[1]
WORK = HERE / "work"
W540 = ROOT / "experiments/540-lodo-clip-scan/work"
FOLDS = ("k562_gwps", "orion_hct116", "orion_hek293t", "kolf2")
CANCER = FOLDS[:3]


def load(base):
    D = {}
    for f in FOLDS:
        z = np.load(W540 / f"arrays540-{f}.npz")
        x, a, d = z[f"{base}__x"], z[f"{base}__a"], z["d"]
        t = [str(s) for s in np.load(WORK / f"dump518-{f}.npz")["targets"]]
        u = (x.sum(0) / a.sum(0)) / (x.sum() / a.sum())
        D[f] = dict(x=x, a=a, d=d, t=t, lu=dict(zip(t, np.log(np.maximum(u, 0.05)))),
                    base=float(np.mean([x[j].sum() / np.sqrt(a[j].sum() * d.sum()) for j in range(3)])))
    return D


def feat(D, h, train):
    t = D[h]["t"]
    v = np.array([np.mean([D[f]["lu"][k] for f in train if k in D[f]["lu"]]) if any(k in D[f]["lu"] for f in train) else np.nan
                  for k in t])
    r = crank(v)
    r[~np.isfinite(v)] = 0.0
    return r


def gain(D, h, z, th):
    x, a, d = D[h]["x"], D[h]["a"], D[h]["d"]
    m = np.exp(th * z)
    return float(np.mean([(m * x[j]).sum() / np.sqrt((m * m * a[j]).sum() * d.sum()) for j in range(3)])) / D[h]["base"] - 1


def main():
    out = {}
    for base in ("V1_0.5_3@R16", "V1_0.5_6@R16"):
        D = load(base)
        res = {}
        for h in FOLDS:
            train = [f for f in CANCER if f != h]
            zs = {f: feat(D, f, [g for g in train if g != f]) for f in train} if len(train) > 1 else {}
            if len(train) >= 2:
                obj = lambda th: -np.mean([gain(D, f, zs[f], th) for f in train]) if all(len([g for g in train if g != f]) for f in train) else 0
            th = minimize_scalar(lambda th: -np.mean([gain(D, f, feat(D, f, [g for g in train if g != f]), th) for f in train]),
                                 bounds=(-2, 3), method="bounded").x
            zh = feat(D, h, train)
            res[h] = dict(theta=float(th), heldout=gain(D, h, zh, th), at_theta1=gain(D, h, zh, 1.0))
        out[base] = res
        print(f"{base}: " + " | ".join(f"{h[:10]} held-out {100 * v['heldout']:+.2f} % (theta {v['theta']:+.2f}; at theta 1 {100 * v['at_theta1']:+.2f} %)"
                                        for h, v in res.items()))
    (WORK / "xfold518.json").write_text(json.dumps(out, indent=1, default=float))


if __name__ == "__main__":
    main()
