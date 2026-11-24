"""410 check: PLOLO folds with the clean key (key410), energy-normalised weights, and the final per-target scale
distribution (what a welded arm would do to the median target's DE channel)."""
import json, sys
from pathlib import Path
import numpy as np
ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "experiments/410-pertscale-revisit"))
sys.path.insert(0, str(ROOT / "experiments/162-retrodicting-validator"))
import plolo
from key410 import shared_energy_key, weights, GAMMA
rng = np.random.default_rng(20260922)
res = {}
RULES = [(0.1, 1e9), (0.2, 1e9), (0.3, 1e9), (0.2, 3.0), (0.3, 3.0), (0.3, 2.0), (0.5, 2.0), (0.5, 3.0), (0.5, 1e9)]
for held, fold in plolo.FOLDS.items():
    z = np.load(ROOT / f"experiments/162-retrodicting-validator/work/fold-{held}.npz", allow_pickle=True)
    pre = np.asarray(z["pre"], float); truth = np.asarray(z["truth"], float)
    scored = np.asarray(z["scored"]); wcols = np.asarray(z["wcols"])
    targets = [str(t) for t in z["targets"]]
    op = plolo.Operator(fold["op"])
    P = op.apply(pre, wcols, GAMMA, 0.0)[:, scored]; T = truth[:, scored]
    live = (np.abs(T).sum(1) > 0) & (np.abs(P).sum(1) > 0)
    P, T = P[live], T[live]; tg = [t for t, k in zip(targets, live) if k]; n = len(tg)
    key, carried, nd = shared_energy_key(tg, fold["sources"], op, wcols, cols=scored)
    pt = (P*T).sum(1); pp = (P*P).sum(1)
    def agg(w, idx=None):
        idx = np.arange(n) if idx is None else idx
        return (w[idx]*pt[idx]).sum() / np.sqrt((w[idx]**2*pp[idx]).sum() * (T[idx]**2).sum())
    base = agg(np.ones(n))
    sstar = pt/pp
    rho = float(np.corrcoef(np.argsort(np.argsort(key[carried])), np.argsort(np.argsort(sstar[carried])))[0,1])
    out = dict(n=n, carried=int(carried.sum()), spearman_sstar=rho, rules={})
    print(f"{held}: n {n} carried {int(carried.sum())} base agg {base:.5f} rho(key,s*) {rho:+.3f}")
    for lo, hi in RULES:
        w, med = weights(key, carried, lo, hi)
        norm = np.sqrt(pp.sum() / (w**2 * pp).sum())       # energy-preserving
        wf = w * norm
        gain = agg(w)/base - 1
        bs = [agg(w, i)/agg(np.ones(n), i) - 1 for i in (rng.integers(0, n, n) for _ in range(400))]
        q = np.quantile(wf, [0.1, 0.25, 0.5, 0.75, 0.9])
        out["rules"][f"lo{lo}_hi{hi}"] = dict(gain=float(gain), ci=[float(np.quantile(bs,.025)), float(np.quantile(bs,.975))],
                                             final_scale_q=q.tolist(), frac_below_half=float((wf < 0.5).mean()))
        print(f"   lo {lo:<4} hi {hi:<6g} gain {gain:+.3f} CI [{np.quantile(bs,.025):+.3f},{np.quantile(bs,.975):+.3f}]  "
              f"final scale q10/25/50/75/90 {np.round(q,2).tolist()}  frac<0.5 {(wf<0.5).mean():.2f}")
    res[held] = out
(ROOT / "experiments/410-pertscale-revisit/work/check410.json").write_text(json.dumps(res, indent=1))
