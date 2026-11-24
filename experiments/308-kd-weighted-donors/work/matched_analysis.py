"""308 -- the target-set-matched control analysis for the only 3/3 cell (arithmetic on matched308's saved arrays)."""
import json, sys
import numpy as np
from pathlib import Path
W = Path(__file__).resolve().parent
CAN = ("k562_gwps", "orion_hct116", "orion_hek293t")
rows, pt = {}, {}
for f in CAN:
    for r in json.loads((W / f"matched308-{f}.json").read_text())["rows"]:
        rows[(f, r["arm"])] = r
    z = np.load(W / f"matchedpt308-{f}.npz", allow_pickle=True)
    pt[f] = {k: np.asarray(z[k]) for k in z.files}
rng = np.random.default_rng(20260907); B = 20000
def bp(diffs):
    acc = np.zeros(B)
    for d in diffs:
        d = d[np.isfinite(d)]; acc += d[rng.integers(0, d.size, size=(B, d.size))].mean(1)
    m = acc / len(diffs); lo, hi = np.percentile(m, [2.5, 97.5])
    return float(np.mean([float(np.nanmean(x)) for x in diffs])), float(lo), float(hi)
def nm(a):
    U = np.concatenate([np.load(W / f"matchedid308-{f}.npz")[f"{a}|U"] for f in CAN])
    M = np.concatenate([np.load(W / f"matchedid308-{f}.npz")[f"{a}|M"] for f in CAN])
    D = np.concatenate([np.load(W / f"matchedid308-{f}.npz")[f"{a}|D"] for f in CAN])
    return float(np.mean(1 + (U - 2 * M) / D))
arms = [r["arm"] for r in json.loads((W / "matched308-k562_gwps.json").read_text())["rows"]]
BASE, ARM = "cd4svd5", "soft:abs:p1.0:dis"
L, out = [], {}
def P(s=""):
    print(s, flush=True); L.append(s)
P(f"{'row':34s} {'pooled dpds [95% CI]':34s} folds>0  per-fold                nmae d")
for a in arms:
    if a == BASE: continue
    diffs = [pt[f][f"{a}|pds"] - pt[f][f"{BASE}|pds"] for f in CAN]
    p, lo, hi = bp(diffs); pf = [float(np.nanmean(d)) for d in diffs]
    out[a] = dict(point=p, lo=lo, hi=hi, per_fold=pf, d_nmae=nm(a) - nm(BASE))
    P(f"{a:34s} {p:+.5f} [{lo:+.5f}, {hi:+.5f}]   {sum(1 for x in pf if x>0)}/3   "
      + " ".join(f"{x:+.4f}" for x in pf) + f"   {out[a]['d_nmae']:+.5f}")
P("")
P("ARM minus each MATCHED control (paired over targets, pooled over the three folds):")
for c in [a for a in arms if a not in (BASE, ARM)]:
    diffs = [pt[f][f"{ARM}|pds"] - pt[f][f"{c}|pds"] for f in CAN]
    p, lo, hi = bp(diffs); pf = [float(np.nanmean(d)) for d in diffs]
    out[f"{ARM} - {c}"] = dict(point=p, lo=lo, hi=hi, per_fold=pf)
    P(f"  {ARM} - {c:36s} {p:+.5f} [{lo:+.5f}, {hi:+.5f}]  folds>0 {sum(1 for x in pf if x>0)}/3")
(W / "matched_analysis.json").write_text(json.dumps(out, indent=2))
(W / "matched_analysis.log").write_text("\n".join(L) + "\n")
