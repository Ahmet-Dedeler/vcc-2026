"""Per-gene decomposition of a variant's nmae change (H1, from stats npz): by direction, sign agreement, |lfc| quantile."""
import sys, json, numpy as np
from pathlib import Path
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import run423 as RN
E = RN.E
z = np.load(RN.ROOT / "experiments/323-current-family-validation/work/reference.npz")
genes_k = z["genes"][z["keep"]]; targets = z["targets"]; RQ = z["q"]; RL = z["lfc"]
base = np.load(HERE / f"work/h1/stats-base_{sys.argv[2]}.npz")
for v in sys.argv[1].split(","):
    st = np.load(HERE / f"work/h1/stats-{v}_{sys.argv[2]}.npz")
    rows = []
    for i in range(len(targets)):
        S = (RQ[i] < .05) & np.isfinite(RL[i]) & (RL[i] != 0) & (genes_k != targets[i])
        if S.sum() < 10: continue
        l0 = np.nan_to_num(base["lfc"][i]); l1 = np.nan_to_num(st["lfc"][i]); r = RL[i]
        den = np.abs(r[S]).sum()
        c = np.where(S, (np.abs(l1 - r) - np.abs(l0 - r)) / den, 0.0)
        moved = S & (l1 != l0)
        a = np.abs(l0); thr = np.quantile(a[np.isfinite(base["lfc"][i])], 0.8)
        rows.append(dict(tot=c.sum(), n_moved=int(moved.sum()), nS=int(S.sum()),
                         right=c[moved & (np.sign(l0) == np.sign(r))].sum(), wrong=c[moved & (np.sign(l0) != np.sign(r))].sum(),
                         up=c[moved & (l0 > 0)].sum(), dn=c[moved & (l0 < 0)].sum(),
                         top20=c[moved & (a >= thr)].sum(), rest=c[moved & (a < thr)].sum(),
                         acc_moved=float((np.sign(l0) == np.sign(r))[moved].mean()) if moved.any() else np.nan,
                         acc_moved_w=float((np.abs(l0[moved]) * (np.sign(l0) == np.sign(r))[moved]).sum() / max(np.abs(l0[moved]).sum(), 1e-12)) if moved.any() else np.nan,
                         ratio_moved=float(np.median(np.abs(l0[moved]) / np.abs(r[moved]))) if moved.any() else np.nan,
                         overshoot=int((moved & (np.sign(l1) == np.sign(r)) & (np.abs(l1) > np.abs(r))).sum())))
    agg = {k: float(np.nanmean([rw[k] for rw in rows])) for k in rows[0]}
    print(v, json.dumps({k: round(x, 5) for k, x in agg.items()}))
    print("   targets with gain<-0.005:", sum(rw["tot"] < -0.005 for rw in rows), " loss>+0.005:", sum(rw["tot"] > 0.005 for rw in rows), " total:", len(rows))

# direction-split sign accuracy of moved genes
if len(sys.argv) > 3:
    for v in sys.argv[1].split(","):
        st = np.load(HERE / f"work/h1/stats-{v}_{sys.argv[2]}.npz")
        acc = {"up": [], "dn": []}; rat = {"up": [], "dn": []}; n = {"up": 0, "dn": 0}
        for i in range(len(targets)):
            S = (RQ[i] < .05) & np.isfinite(RL[i]) & (RL[i] != 0) & (genes_k != targets[i])
            if S.sum() < 10: continue
            l0 = np.nan_to_num(base["lfc"][i]); l1 = np.nan_to_num(st["lfc"][i]); r = RL[i]
            moved = S & (l1 != l0)
            for d, m in (("up", moved & (l0 > 0)), ("dn", moved & (l0 < 0))):
                acc[d] += list((np.sign(l0) == np.sign(r))[m]); rat[d] += list((np.abs(l0) / np.abs(r))[m]); n[d] += int(m.sum())
        print(v, {d: (n[d], round(float(np.mean(acc[d])), 3), round(float(np.median(rat[d])), 3)) for d in acc})
