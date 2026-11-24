"""410 x 228: do the per-TARGET key and the per-GENE gain stack in aggregate transfer cosine (PLOLO, cancer folds)?"""
import sys
from pathlib import Path
import numpy as np
ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "experiments/410-pertscale-revisit")); sys.path.insert(0, str(ROOT / "experiments/162-retrodicting-validator"))
import plolo
from key410 import shared_energy_key, weights, GAMMA
rng = np.random.default_rng(5)
for held in ("k562_gwps", "orion_hct116", "orion_hek293t"):
    fold = plolo.FOLDS[held]
    z = np.load(ROOT / f"experiments/162-retrodicting-validator/work/fold-{held}.npz", allow_pickle=True)
    pre = np.asarray(z["pre"], float); truth = np.asarray(z["truth"], float)
    scored = np.asarray(z["scored"]); wcols = np.asarray(z["wcols"]); targets = [str(t) for t in z["targets"]]
    op = plolo.Operator(fold["op"])
    W = op.apply(pre, wcols, GAMMA, 0.0)
    g = np.load(ROOT / f"experiments/228-theta-gain-plolo/work/gain-{held}.npz", allow_pickle=True)
    live = np.asarray(g["live"], bool); gr = np.asarray(g["gain_reg"], float)
    med = float(np.median(gr[live])); m = np.where(live, 1.0 + 0.25 * (gr / med - 1.0), 1.0)
    variants = {"base": W, "gain": W * m[None, :]}
    T = truth[:, scored]
    live_t = np.abs(T).sum(1) > 0
    tg = [t for t, k in zip(targets, live_t) if k]
    key, carried, _ = shared_energy_key(tg, fold["sources"], op, wcols, cols=scored)
    out = []
    for vname, V in variants.items():
        P = V[:, scored][live_t]; TT = T[live_t]
        pt = (P * TT).sum(1); pp = (P * P).sum(1)
        def agg(w): return (w * pt).sum() / np.sqrt((w * w * pp).sum() * (TT * TT).sum())
        for lo, hi in ((None, None), (0.3, 3.0), (0.3, 1e9)):
            w = np.ones(len(tg)) if lo is None else weights(key, carried, lo, hi)[0]
            out.append((vname, "none" if lo is None else f"key {lo}/{hi:g}", agg(w)))
    base = out[0][2]
    print(held, " ".join(f"[{a} + {b}: {c/base-1:+.3f}]" for a, b, c in out))
