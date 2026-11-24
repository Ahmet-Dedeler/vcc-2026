"""410: do extra lines (KOLF2 iPSC, CD4 T cells x3 conditions) sharpen the shared-energy key? Cancer PLOLO folds.
Every donor's rows are divided by that donor's median row norm over the fold targets (scale-free mixing)."""
import sys
from pathlib import Path
import numpy as np
ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "experiments/410-pertscale-revisit")); sys.path.insert(0, str(ROOT / "experiments/162-retrodicting-validator"))
import plolo
from key410 import raw_whiten, weights, GAMMA
import size410
CD4 = ROOT / "experiments/395-protected-signs-arm/work/build/cd4_compact"
plolo.DONOR_PATHS.update({"cd4_rest": CD4 / "cd4_rest", "cd4_stim8hr": CD4 / "cd4_stim8hr", "cd4_stim48hr": CD4 / "cd4_stim48hr"})
rng = np.random.default_rng(1)
def key_from(rows, H, op, wcols, cols, normalize):
    n = rows[0].shape[0]
    if normalize:
        rows = [r / np.median(np.linalg.norm(r[h], axis=1)) if h.any() else r for r, h in zip(rows, H)]
    Wr = [raw_whiten(op, r, wcols)[:, cols] for r in rows]; R = [r[:, cols] for r in rows]
    key = np.full(n, np.nan)
    for p in range(n):
        idx = [k for k in range(len(R)) if H[k][p]]
        if len(idx) < 2: continue
        ips = [0.5 * (Wr[i][p] @ R[j][p] + Wr[j][p] @ R[i][p]) for a, i in enumerate(idx) for j in idx[a + 1:]]
        m = np.mean([Wr[k][p] for k in idx], axis=0); key[p] = np.mean(ips) / max(m @ m, 1e-300)
    return key
for held in ("k562_gwps", "orion_hct116", "orion_hek293t"):
    fold = plolo.FOLDS[held]
    z = np.load(ROOT / f"experiments/162-retrodicting-validator/work/fold-{held}.npz", allow_pickle=True)
    pre = np.asarray(z["pre"], float); truth = np.asarray(z["truth"], float)
    scored = np.asarray(z["scored"]); wcols = np.asarray(z["wcols"]); targets = [str(t) for t in z["targets"]]
    op = plolo.Operator(fold["op"])
    P = op.apply(pre, wcols, GAMMA, 0.0)[:, scored]; T = truth[:, scored]
    live = (np.abs(T).sum(1) > 0) & (np.abs(P).sum(1) > 0)
    P, T = P[live], T[live]; tg = [t for t, k in zip(targets, live) if k]; n = len(tg)
    pt = (P * T).sum(1); pp = (P * P).sum(1); tt = (T * T).sum(1); ss = pt / pp
    def agg(w): return (w * pt).sum() / np.sqrt((w * w * pp).sum() * tt.sum())
    base = agg(np.ones(n))
    cache = {}
    def rows_of(name):
        if name not in cache: cache[name] = size410.donor_rows(name, tg)
        return cache[name]
    variants = {"V0 raw": (list(fold["sources"]), False), "V0 norm": (list(fold["sources"]), True),
                "+kolf2": (list(fold["sources"]) + ["kolf2"], True),
                "+cd4x3": (list(fold["sources"]) + ["cd4_rest", "cd4_stim8hr", "cd4_stim48hr"], True),
                "+kolf2+cd4x3": (list(fold["sources"]) + ["kolf2", "cd4_rest", "cd4_stim8hr", "cd4_stim48hr"], True)}
    print(f"{held}: n {n} base {base:.5f}")
    for vname, (donors, norm) in variants.items():
        rr = [rows_of(d) for d in donors]
        key = key_from([r for r, h in rr], [h for r, h in rr], op, wcols, scored, norm)
        carried = np.isfinite(key)
        rho = np.corrcoef(np.argsort(np.argsort(key[carried])), np.argsort(np.argsort(ss[carried])))[0, 1]
        out = []
        for lo, hi in ((0.3, 3.0), (0.3, 1e9)):
            w, _ = weights(np.nan_to_num(key, nan=np.nanmedian(key)), carried, lo, hi)
            out.append(agg(w) / base - 1)
        print(f"   {vname:14s} carried {carried.sum():3d}  rho(key,s*) {rho:+.3f}  gain hi3 {out[0]:+.3f}  uncapped {out[1]:+.3f}")
