"""410 robustness: wshare key with strictly positive floors, bootstrap CI over targets, key distribution."""
import json, sys
from pathlib import Path
import numpy as np
ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "experiments/410-pertscale-revisit"))
sys.path.insert(0, str(ROOT / "experiments/162-retrodicting-validator"))
import plolo
from size410 import donor_rows, GAMMA, LAM, OFFICIAL
rng = np.random.default_rng(20260922)
res = {}
for held, fold in plolo.FOLDS.items():
    z = np.load(ROOT / f"experiments/162-retrodicting-validator/work/fold-{held}.npz", allow_pickle=True)
    pre = np.asarray(z["pre"], float); truth = np.asarray(z["truth"], float)
    scored = np.asarray(z["scored"]); wcols = np.asarray(z["wcols"])
    targets = [str(t) for t in z["targets"]]
    op = plolo.Operator(fold["op"])
    P = op.apply(pre, wcols, GAMMA, LAM)[:, scored]; T = truth[:, scored]
    live = (np.abs(T).sum(1) > 0) & (np.abs(P).sum(1) > 0)
    P, T = P[live], T[live]; tg = [t for t, k in zip(targets, live) if k]; n = len(tg)
    rows, H = [], []
    for s in fold["sources"]:
        r, h = donor_rows(s, tg); rows.append(r); H.append(h)
    H = np.vstack(H)
    W = []
    for r in rows:
        W.append(op.apply(r, wcols, GAMMA, LAM)[:, scored])
    R = [r[:, scored] for r in rows]
    key = np.full(n, np.nan)
    for p in range(n):
        have = [k for k in range(len(R)) if H[k, p]]
        if len(have) < 2: continue
        ips = [0.5*(W[a][p] @ R[b][p] + W[b][p] @ R[a][p]) for i, a in enumerate(have) for b in have[i+1:]]
        m = np.mean([W[k][p] for k in have], axis=0)
        key[p] = np.mean(ips) / max(m @ m, 1e-12)
    # uncarried -> median key (neutral)
    med = np.nanmedian(key); keyf = np.where(np.isnan(key), med, key)
    pt = (P*T).sum(1); pp = (P*P).sum(1); d = (T*T).sum()
    def agg(w, idx=None):
        idx = np.arange(n) if idx is None else idx
        return (w[idx]*pt[idx]).sum() / np.sqrt((w[idx]**2*pp[idx]).sum() * (T[idx]**2).sum())
    base = agg(np.ones(n))
    q = np.quantile(keyf, [0, .05, .1, .25, .5, .75, .9, .95, 1])
    out = dict(n=n, key_quantiles=q.tolist(), frac_le0=float((keyf <= 0).mean()), rules={})
    # normalise key by its median so floors are relative
    kn = keyf / med
    for lo in (0.02, 0.05, 0.1, 0.2, 0.3):
        for hi in (1.0, 2.0, 3.0, 1e9):
            for g in (0.5, 1.0):
                w = np.clip(kn, lo, hi) ** g
                gain = agg(w)/base - 1
                bs = []
                for _ in range(400):
                    idx = rng.integers(0, n, n)
                    bs.append(agg(w, idx)/agg(np.ones(n), idx) - 1)
                out["rules"][f"lo{lo}_hi{hi}_g{g}"] = [float(gain), float(np.quantile(bs, .025)), float(np.quantile(bs, .975)), float((w <= lo**g + 1e-12).mean())]
    res[held] = out
    best = sorted(out["rules"].items(), key=lambda kv: -kv[1][0])[:6]
    print(held, "n", n, "base", round(base, 5), "frac<=0", round(out["frac_le0"], 3), "q", np.round(q, 3).tolist())
    for k, v in best: print("   ", k, "gain %+.3f CI [%+.3f, %+.3f] at-floor %.2f" % tuple(v))
    for k in ("lo0.1_hi2.0_g1.0", "lo0.2_hi2.0_g1.0", "lo0.1_hi1000000000.0_g1.0", "lo0.2_hi3.0_g0.5"):
        v = out["rules"][k]; print("    [ref]", k, "gain %+.3f CI [%+.3f, %+.3f] at-floor %.2f" % tuple(v))
(ROOT / "experiments/410-pertscale-revisit/work/robust410.json").write_text(json.dumps(res, indent=1))
