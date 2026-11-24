"""450k: does a NEAR line of the same tissue carry information? Truth = Orion HCT116 (colon, MSI-high). Near line = DLD-1
low-MOI screen (colon, MSI-high; different lab). Cross-line = K562-GWPS, Orion HEK293T. Per-target centred cosine on common
targets (sample up to 1,500, seed 20260925) and genes, each target's own gene zeroed; and the programme (P_down over rows,
the evaluated targets excluded) as a sign predictor of the truth rows' signs on their top-100 genes."""
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
B = {"truth_hct116": ROOT / "data/donors/orion_hct116", "dld1 (near, colon)": ROOT / "experiments/223-silent-zero-audit/work/tables/dld1_lowmoi",
     "k562 (cross)": ROOT / "data/donors/k562_gwps", "hek293t (cross)": ROOT / "data/donors/orion_hek293t"}
L = {k: (np.load(d / "delta.npy", mmap_mode="r"), [l.strip() for l in open(d / "perts.txt")], [l.strip() for l in open(d / "genes.txt")]) for k, d in B.items()}
ct = set(L["truth_hct116"][1])
for k in B:
    ct &= set(L[k][1])
cg = set(L["truth_hct116"][2])
for k in B:
    cg &= set(L[k][2])
cg = sorted(cg)
ess = set(l.strip() for l in open(ROOT / "data/donors/k562_essential/perts.txt"))
rng = np.random.default_rng(20260925)
print(f"common targets {len(ct)}, common genes {len(cg)}")
for sname, pool in (("non-essential", sorted(t for t in ct if t not in ess)), ("essential", sorted(t for t in ct if t in ess))):
    ts = list(rng.choice(pool, size=min(1500, len(pool)), replace=False))
    M = {}
    for k, (D, P, G) in L.items():
        pi = {p: i for i, p in enumerate(P)}; gi = {g: i for i, g in enumerate(G)}
        rows = np.array([pi[t] for t in ts]); cols = np.array([gi[g] for g in cg]); o = np.argsort(rows)
        X = np.empty((len(ts), len(cols))); X[o] = np.asarray(D[rows[o]], dtype=np.float64)[:, cols]
        own = {g: j for j, g in enumerate(cg)}
        for r, t in enumerate(ts):
            if t in own:
                X[r, own[t]] = 0.0
        M[k] = X
    T = M["truth_hct116"]; Tc = T - T.mean(0)
    u = lambda A: A / np.maximum(np.linalg.norm(A, axis=1, keepdims=True), 1e-12)
    top = np.argsort(-np.abs(T), axis=1)[:, :100]
    print(f"\n== {sname}: {len(ts)} targets")
    for k in ("dld1 (near, colon)", "k562 (cross)", "hek293t (cross)"):
        P = M[k]; Pc = P - P.mean(0)
        cc = np.sum(u(Pc) * u(Tc), 1).mean()
        sg = np.mean(np.sign(np.take_along_axis(P, top, 1)) == np.sign(np.take_along_axis(T, top, 1)))
        # programme: P_down of the bank over its rows NOT in the evaluated set, as a per-gene sign for every target
        D, Pn, G = L[k]; gi = {g: i for i, g in enumerate(G)}; cols = np.array([gi[g] for g in cg])
        rest = np.array([i for i, p in enumerate(Pn) if p not in set(ts)])
        sub = rest if len(rest) <= 3000 else rng.choice(rest, 3000, replace=False)
        R = np.asarray(D[np.sort(sub)], dtype=np.float64)[:, cols]
        pd = (R < 0).sum(0) / np.maximum((R != 0).sum(0), 1)
        prog = np.where(pd > 0.5, -1.0, 1.0)
        ps = np.mean(np.take_along_axis(np.broadcast_to(prog, T.shape), top, 1) == np.sign(np.take_along_axis(T, top, 1)))
        print(f"  {k:20s} centred cos {cc:.4f}   row sign@top100 {sg:.4f}   programme sign@top100 {ps:.4f}")
