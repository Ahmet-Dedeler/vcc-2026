"""Reviewer check: Testh's disclosed 'cross-target geometry' lever priced on h1-test (analytic, bulk only).

D' = D - beta * (D V_k^T) V_k, V_k = top-k right-singular vectors of the stacked row-normalised donor deltas of the
evaluated targets (other lines only: K562-GWPS, Orion HCT116, Orion HEK293T; never H1 data). Applied to our predicted
bulk effects (435's harness sums), pds recomputed from the spec (vcc2026 brief section 1). Known answer: 435's metrics
json pds 0.8952380952 (base, seed 0). Controls: a random orthonormal basis of the same rank; a basis from our own
predictions (self-deflation).
"""
import glob, json, sys
from pathlib import Path
import numpy as np

R = Path("~/Code/virtual-cell-challenge-2026")
TS = 5e4
ref = np.load(glob.glob(str(R / "experiments/323-*/work/reference.npz"))[0], allow_pickle=True)
genes = ref["genes"].astype(str); targets = list(ref["targets"].astype(str)); perts = list(ref["perts"].astype(str))
S_true = np.asarray(ref["sums"], float)
st = np.load(glob.glob(str(R / "experiments/435-zs-on-coarse/work/h1/stats-base_s0_bps_prodem_mde0.1cap0.1min0.3coarse.npz"))[0])
assert list(st["targets"].astype(str)) == targets
S_pred = np.asarray(st["sums"], float)
ctrl_i = [i for i, p in enumerate(perts) if p not in set(targets)]
assert len(ctrl_i) == 1, perts[:3]
def logb(S): return np.log1p(TS * S / S.sum(1, keepdims=True))
b_ctrl = logb(S_true[ctrl_i])[0]
tix = [perts.index(t) for t in targets]
D_true = logb(S_true[tix]) - b_ctrl
D_pred = logb(S_pred) - b_ctrl

def pds(Dp, Dt, drop):
    keep = ~drop
    P = Dp[:, keep]; T = Dt[:, keep]
    nP = np.linalg.norm(P, axis=1); nT = np.linalg.norm(T, axis=1)
    C = (P @ T.T) / np.outer(np.where(nP > 0, nP, 1), np.where(nT > 0, nT, 1))
    d = 1 - C
    d[(nP == 0)[:, None] | (nT == 0)[None, :]] = 1.0
    n = len(d); out = []
    for p in range(n):
        own = d[p, p]
        less = (d[p] < own).sum(); ties = (d[p] == own).sum()
        k = less + 0.5 * (ties - 1)
        out.append(1 - k / (n - 1))
    return float(np.mean(out)), np.array(out)

gpos = {g: i for i, g in enumerate(genes)}
drop85 = np.zeros(len(genes), bool)
for t in targets:
    if t in gpos: drop85[gpos[t]] = True
want = json.load(open(glob.glob(str(R / "experiments/435-zs-on-coarse/work/h1/metrics-base_s0_bps_prodem_mde0.1cap0.1min0.3coarse.json"))[0]))["raw"]["pds"]
base, base_pt = pds(D_pred, D_true, drop85)
print(f"KA pds: mine {base:.10f} vs harness {want:.10f} -> {'PASS' if abs(base - want) < 1e-9 else 'FAIL'}  (dropped {drop85.sum()} target genes)")
if abs(base - want) >= 1e-9:
    sys.exit(1)

# donor rows for the evaluated targets, mapped onto the reference axis (memory-mapped; only these rows are read)
rows = []
for bank in ("k562_gwps", "orion_hct116", "orion_hek293t"):
    d = R / "data/donors" / bank
    dl = np.load(d / "delta.npy", mmap_mode="r")
    bp = [l.strip() for l in open(d / "perts.txt")]; bg = [l.strip() for l in open(d / "genes.txt")]
    col = np.array([gpos.get(g, -1) for g in bg]); ok = col >= 0
    pidx = {p: i for i, p in enumerate(bp)}
    n = 0
    for t in targets:
        if t in pidx:
            v = np.zeros(len(genes)); v[col[ok]] = np.asarray(dl[pidx[t]], float)[ok]
            v[drop85] = 0.0
            nv = np.linalg.norm(v)
            if nv > 0: rows.append(v / nv); n += 1
    print(f"{bank}: {n} of {len(targets)} targets")
M = np.vstack(rows)
_, _, Vt = np.linalg.svd(M, full_matrices=False)
rng = np.random.default_rng(0)
res = {"ka_pds": base, "grid": []}
Dself = D_pred.copy(); Dself[:, drop85] = 0
_, _, Vself = np.linalg.svd(Dself / np.linalg.norm(Dself, axis=1, keepdims=True), full_matrices=False)
for k in (1, 3, 5, 10, 14, 20, 30):
    Vk = Vt[:k]
    Q, _ = np.linalg.qr(rng.standard_normal((len(genes), k))); Vr = Q.T
    for beta in (0.25, 0.5, 0.75, 1.0):
        row = {"k": k, "beta": beta}
        for name, V in (("donor", Vk), ("random", Vr), ("self", Vself[:k])):
            Dp = D_pred - beta * (D_pred @ V.T) @ V
            v, pt = pds(Dp, D_true, drop85)
            # squared-error proxy for mse's numerator: sum ||b'-b||^2 relative to the base
            se = float(((Dp - D_true)[:, ~drop85] ** 2).sum()) / float(((D_pred - D_true)[:, ~drop85] ** 2).sum())
            row[name] = {"dpds": round(v - base, 5), "se_ratio": round(se, 4), "wins": int((pt > base_pt).sum()), "losses": int((pt < base_pt).sum())}
        res["grid"].append(row)
        print(f"k {k:>2} beta {beta:.2f} | donor dpds {row['donor']['dpds']:+.5f} se {row['donor']['se_ratio']:.4f} (+{row['donor']['wins']}/-{row['donor']['losses']}) "
              f"| random {row['random']['dpds']:+.5f} se {row['random']['se_ratio']:.4f} | self {row['self']['dpds']:+.5f} se {row['self']['se_ratio']:.4f}")
Path(__file__).with_suffix(".json").write_text(json.dumps(res, indent=1))
