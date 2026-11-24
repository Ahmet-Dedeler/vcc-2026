"""549: which genes differ, and by how much, in the 5 (context, target) groups where 446bj2zcph's summed counts differ from
446bj2zcp's (bulk_identity550.py found max |diff| 1 count). Reads only those groups' rows (backed CSR)."""
import sys
import h5py
import numpy as np

L, R = sys.argv[1], sys.argv[2]
GROUPS = [("A", "MAPK12"), ("B", "AKT2"), ("B", "GRAMD1A"), ("C", "STK3"), ("C", "ZBTB48")]

def col(obs, name):
    n = obs[name]
    if isinstance(n, h5py.Group):
        cats = [c.decode() if isinstance(c, bytes) else str(c) for c in n["categories"][:]]
        return np.array([cats[i] for i in n["codes"][:]])
    return np.array([x.decode() if isinstance(x, bytes) else str(x) for x in n[:]])

def sums(path):
    out = {}
    with h5py.File(path, "r") as f:
        ctx, tgt = col(f["obs"], "context"), col(f["obs"], "target_gene")
        var = f["var"]; idx = var.attrs.get("_index", "_index")
        genes = np.array([g.decode() if isinstance(g, bytes) else str(g) for g in var[idx][:]])
        ip, ind, dat = f["X/indptr"], f["X/indices"], f["X/data"]
        for c, t in GROUPS:
            rows = np.where((ctx == c) & (tgt == t))[0]
            v = np.zeros(len(genes))
            for r in rows:
                s, e = int(ip[r]), int(ip[r + 1])
                np.add.at(v, ind[s:e], dat[s:e])
            out[(c, t)] = (v, len(rows))
    return out, genes

a, genes = sums(L); b, _ = sums(R)
for k in GROUPS:
    va, n = a[k]; vb, _ = b[k]
    d = np.nonzero(va != vb)[0]
    tot = va.sum()
    for g in d:
        print(f"{k[0]} {k[1]:8s} cells {n}  gene {genes[g]:12s} cph {va[g]:.0f}  cp {vb[g]:.0f}  (group total {tot:.0f}; "
              f"rel {abs(va[g] - vb[g]) / max(vb[g], 1):.2e})  totals equal: {va.sum() == vb.sum()}")
