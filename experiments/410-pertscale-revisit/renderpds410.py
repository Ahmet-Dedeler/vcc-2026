"""Does the RENDER (exp, normalise, log read-back) make per-target scaling cost pds? PLOLO cancer folds, challenge b_ctrl."""
import sys
from pathlib import Path
import numpy as np
ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "experiments/410-pertscale-revisit")); sys.path.insert(0, str(ROOT / "experiments/162-retrodicting-validator"))
import plolo
from key410 import shared_energy_key, weights, GAMMA
for held in ("k562_gwps", "orion_hct116", "orion_hek293t"):
    fold = plolo.FOLDS[held]
    z = np.load(ROOT / f"experiments/162-retrodicting-validator/work/fold-{held}.npz", allow_pickle=True)
    pre = np.asarray(z["pre"], float); truth = np.asarray(z["truth"], float)
    scored = np.asarray(z["scored"]); wcols = np.asarray(z["wcols"]); targets = [str(t) for t in z["targets"]]
    gate = np.asarray(z["gate"], bool)
    op = plolo.Operator(fold["op"])
    W = op.apply(pre * 0.982, wcols, GAMMA, 0.0)
    key, carried, _ = shared_energy_key(targets, fold["sources"], op, wcols, cols=scored)
    T = truth[:, scored]
    res = []
    for ctx in ("A", "B", "C"):
        b = plolo.challenge_bctrl(ctx)
        cg = gate & (np.expm1(b) / np.expm1(b).sum() * 1e6 > 10)
        base = plolo.deliver(W, b, cg)[:, scored]
        row = [plolo.endpoints(base, T)["pds"]]
        for lo, hi in ((0.3, 3.0), (0.3, 1e9), (1.0, 1e9)):
            w, _ = weights(np.nan_to_num(key, nan=np.nanmedian(key)), carried, lo, hi)
            e2 = (W[:, scored] ** 2).sum(1); c = np.sqrt(e2.sum() / (w * w * e2).sum())
            d = plolo.deliver(W * (w * c)[:, None], b, cg)[:, scored]
            row.append(plolo.endpoints(d, T)["pds"])
        res.append(row)
    res = np.array(res)
    print(f"{held}: pds rendered base {res[:,0].mean():.5f}  Δ clip(0.3,3) {res[:,1].mean()-res[:,0].mean():+.5f}  "
          f"Δ clip(0.3,inf) {res[:,2].mean()-res[:,0].mean():+.5f}  Δ clip(1,inf) {res[:,3].mean()-res[:,0].mean():+.5f}  (mean over A/B/C b_ctrl)")

