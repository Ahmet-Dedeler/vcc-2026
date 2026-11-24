"""462h — the lab-signature test in CRISPRi, with cells: Replogle 2022 K562-essential and RPE1-essential (Weissman lab, the
same essential library, two lines; scPerturb single-cell counts). 462b did this on the Broad-16 (Cas9 KO): with disjoint
controls, lines shared about a third of their honest programme. Here, per line:
  shared    P_down(T1 vs full C) vs P_down(T2 vs full C)        (target cells split in half; one control)
  disjoint  P_down(T1 vs C1)     vs P_down(T2 vs C2)            (no shared noise: the honest within-line reliability)
and across lines (independent by construction): K562 half-A vs RPE1 half-B, and full vs full.
Programmes over the targets both screens carry (control out), on the genes both measure; also on genes expressed >= 5
CPM in both lines' controls. usage: uv run python experiments/462-broad16-lab-signature/weissman462.py
"""
import itertools
import json
import time
from pathlib import Path

import h5py
import numpy as np

ROOT = Path(__file__).resolve().parents[2]; HERE = Path(__file__).resolve().parent
FILES = {"K562": ROOT / "data/external/scperturb/ReplogleWeissman2022_K562_essential.h5ad",
         "RPE1": ROOT / "data/external/scperturb/ReplogleWeissman2022_rpe1.h5ad"}
TS, CHUNK = 5e4, 4000
t0 = time.time()
log = lambda *x: print("[462h]", *x, f"({time.time() - t0:.0f}s)", flush=True)


def cat(n):
    return n["categories"][:].astype(str)[n["codes"][:]] if isinstance(n, h5py.Group) else n[:].astype(str)


def sums(path, seed):
    rng = np.random.default_rng(seed)
    with h5py.File(path, "r") as f:
        pert = cat(f["obs"]["perturbation"])
        genes = cat(f["var"]["gene_name"]) if "gene_name" in f["var"] else f["var"]["_index"][:].astype(str)
        X = f["X"]
        keys, kidx = np.unique(pert, return_inverse=True)
        half = rng.integers(0, 2, pert.size)
        S = np.zeros((len(keys), 3, X.shape[1]))
        for s in range(0, pert.size, CHUNK):
            e = min(pert.size, s + CHUNK)
            blk = np.asarray(X[s:e], dtype=np.float64)
            for hv in (0, 1):
                sel = half[s:e] == hv
                np.add.at(S[:, 1 + hv], kidx[s:e][sel], blk[sel])
        S[:, 0] = S[:, 1] + S[:, 2]
        n = np.bincount(kidx, minlength=len(keys))
    return list(keys), genes, S, n


def lp(c):
    return np.log1p(TS * c / c.sum())


def pdn(M):
    nz = (M != 0).sum(0); dn = (M < 0).sum(0)
    with np.errstate(invalid="ignore", divide="ignore"):
        return np.where(nz > 0, dn / np.maximum(nz, 1), np.nan)


def r(a, b, m):
    m = m & np.isfinite(a) & np.isfinite(b)
    return float(np.corrcoef(a[m] - .5, b[m] - .5)[0, 1])


def main():
    D = {}
    for i, (line, path) in enumerate(FILES.items()):
        keys, genes, S, n = sums(path, 462 + i)
        D[line] = dict(keys=keys, genes=genes, S=S, n=n)
        log(f"{line}: {len(keys) - 1} perturbations, {int(n[keys.index('control')])} controls, {len(genes)} genes")
    kg, rg = D["K562"]["genes"], D["RPE1"]["genes"]
    common_g = sorted(set(kg) & set(rg))
    T = sorted((set(D["K562"]["keys"]) & set(D["RPE1"]["keys"])) - {"control"})
    log(f"common targets {len(T)}, common genes {len(common_g)}")
    P, expressed = {}, np.ones(len(common_g), bool)
    for line, d in D.items():
        gi = {g: j for j, g in enumerate(d["genes"])}
        cols = np.array([gi[g] for g in common_g])
        S = d["S"][:, :, cols]
        c = d["keys"].index("control"); C, C1, C2 = S[c, 0], S[c, 1], S[c, 2]
        cpm = 1e6 * C / C.sum(); expressed &= cpm >= 5
        rows = [d["keys"].index(t) for t in T]
        lf = lambda k, ctrl: np.stack([lp(S[i, k]) - lp(ctrl) for i in rows])
        P[line] = {"full": pdn(lf(0, C)), "T1C": pdn(lf(1, C)), "T2C": pdn(lf(2, C)),
                   "T1C1": pdn(lf(1, C1)), "T2C2": pdn(lf(2, C2))}
    res = {"targets": len(T), "genes": len(common_g), "expressed_both": int(expressed.sum())}
    for tag, m in (("all_common_genes", np.ones(len(common_g), bool)), ("expressed_both", expressed)):
        rec = {}
        for line in P:
            rec[f"{line}_shared"] = r(P[line]["T1C"], P[line]["T2C"], m)
            rec[f"{line}_disjoint"] = r(P[line]["T1C1"], P[line]["T2C2"], m)
        rec["cross_halves"] = float(np.mean([r(P["K562"][a], P["RPE1"][b], m)
                                             for a, b in (("T1C1", "T2C2"), ("T2C2", "T1C1"))]))
        rec["cross_full"] = r(P["K562"]["full"], P["RPE1"]["full"], m)
        within = np.mean([rec["K562_disjoint"], rec["RPE1_disjoint"]])
        rec["cross_over_within_disjoint"] = float(rec["cross_halves"] / within)
        res[tag] = rec
    print(json.dumps(res, indent=1))
    (HERE / "work/weissman462.json").write_text(json.dumps(res, indent=1))


if __name__ == "__main__":
    main()
