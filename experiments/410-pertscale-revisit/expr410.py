"""410 side question: does the held line's OWN expression of the target predict the transferable signal (s*, cosine)?"""
import sys, json
from pathlib import Path
import numpy as np
ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "experiments/410-pertscale-revisit")); sys.path.insert(0, str(ROOT / "experiments/162-retrodicting-validator"))
import plolo
from key410 import shared_energy_key, weights, GAMMA
TS = 50000.
OFF = {g: i for i, g in enumerate(plolo.OFFICIAL)}
for held, fold in plolo.FOLDS.items():
    z = np.load(ROOT / f"experiments/162-retrodicting-validator/work/fold-{held}.npz", allow_pickle=True)
    pre = np.asarray(z["pre"], float); truth = np.asarray(z["truth"], float); bctrl = np.asarray(z["b_ctrl"], float)
    scored = np.asarray(z["scored"]); wcols = np.asarray(z["wcols"]); targets = [str(t) for t in z["targets"]]
    op = plolo.Operator(fold["op"])
    P = op.apply(pre, wcols, GAMMA, 0.0)[:, scored]; T = truth[:, scored]
    live = (np.abs(T).sum(1) > 0) & (np.abs(P).sum(1) > 0)
    P, T = P[live], T[live]; tg = [t for t, k in zip(targets, live) if k]
    cpm = np.array([np.expm1(bctrl[OFF[t]]) / TS * 1e6 if t in OFF else np.nan for t in tg])
    pt = (P*T).sum(1); pp = (P*P).sum(1); tt = (T*T).sum(1)
    cos = pt/np.sqrt(pp*tt); ss = pt/pp
    key, carried, _ = shared_energy_key(tg, fold["sources"], op, wcols, cols=scored)
    kn = key/np.median(key[carried])
    bins = [(-1, 1), (1, 5), (5, 10), (10, 30), (30, 100), (100, 1e9)]
    print(f"{held}: n {len(tg)}  held-line target CPM quantiles {np.round(np.nanquantile(cpm,[.05,.25,.5,.75,.95]),1).tolist()}")
    for lo, hi in bins:
        m = (cpm >= lo) & (cpm < hi)
        if m.sum() == 0: continue
        print(f"   CPM [{lo:>4},{hi:>6}) n {m.sum():3d}  mean cos {cos[m].mean():+.4f}  median s* {np.median(ss[m]):+.3f}  "
              f"frac s*<=0 {(ss[m]<=0).mean():.2f}  truth-norm median {np.median(np.sqrt(tt[m])):.2f}  key/med median {np.nanmedian(kn[m]):+.2f}")
    # does a context gate add to the key? agg cosine with key weights vs key*gate
    def agg(w): return (w*pt).sum()/np.sqrt((w*w*pp).sum()*tt.sum())
    w, _ = weights(key, carried, 0.3, 3.0); base = agg(np.ones(len(tg)))
    for thr, f in ((1, 0.1), (5, 0.1), (5, 0.3), (10, 0.3), (10, 0.1)):
        g = np.where(np.nan_to_num(cpm, nan=1e9) < thr, f, 1.0)
        print(f"   gate CPM<{thr} x{f}: alone {agg(g)/base-1:+.3f}  key {agg(w)/base-1:+.3f}  key*gate {agg(w*g)/base-1:+.3f}  (n gated {(g<1).sum()})")
