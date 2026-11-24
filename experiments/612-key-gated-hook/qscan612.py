"""612: the gate quantile q on 480's folds (same machinery as kgh612.py): hybrid(q) = CLEAN4 sign where |key| rank < q,
key sign elsewhere; MDE-law magnitude. Reports nmae per surface and the min-max regret vs the best single source per q."""
import sys
from pathlib import Path
import numpy as np
ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "experiments/480-nmae-instrument")); sys.path.insert(0, str(ROOT / "experiments/537-arm-pipeline"))
import fold480 as F
from common480 import load_donor_rows, nmae_rows
import foldsign537 as FS
C4ALL = ["k562_gwps", "k562_essential", "nadig_jurkat", "nadig_hepg2"]
C4B = {"h1": C4ALL, "orion_hct116": C4ALL, "orion_hek293t": C4ALL, "k562_gwps": ["nadig_jurkat", "nadig_hepg2"]}
QS = [0.0, 0.3, 0.5, 0.65, 0.8, 0.9, 0.95, 1.0]
res = {}
for name in C4B:
    S = F.load_surface(name)
    T, cnt, genes, targets = S["T"], S["cnt"], S["genes"], S["targets"]
    cons0 = np.nan_to_num(np.nanmean(np.stack([load_donor_rows(d, targets, genes, pc_cpm=1.0)[0] for d in S["donors"]]), 0))
    P4 = FS.clean4_P(genes, C4B[name]); s4 = np.where(np.isfinite(P4), np.where(P4 > 0.5, -1.0, 1.0), 0.0)
    S4 = np.broadcast_to(s4, T.shape); sc = np.sign(cons0); add = F.mde_add(cnt)[None, :]
    rank = np.full(T.shape, np.nan)
    for i in range(T.shape[0]):
        a = np.abs(cons0[i]); nz = a > 0
        if nz.sum() > 10: rank[i, nz] = (np.argsort(np.argsort(a[nz])) + 0.5) / nz.sum()
    res[name] = {}
    for q in QS:
        low = np.where(np.isfinite(rank), rank < q, True)
        hyb = np.where(low & (S4 != 0), S4, np.where(sc != 0, sc, S4))
        res[name][q] = float(np.nanmean(nmae_rows(hyb * add, T)))
    res[name]["key"] = float(np.nanmean(nmae_rows(sc * add, T)))
    res[name]["c4"] = float(np.nanmean(nmae_rows(S4 * add, T)))
print("nmae by gate quantile q (q 0 = key everywhere except key-less genes; q 1 = CLEAN4 wherever it has a sign)")
print("%-14s" % "surface" + "".join("%9s" % f"q{q}" for q in QS) + "%9s%9s" % ("key", "C4"))
for n, r in res.items():
    print("%-14s" % n + "".join("%9.5f" % r[q] for q in QS) + "%9.5f%9.5f" % (r["key"], r["c4"]))
print("\nregret = nmae(q) - min over (key, C4, all q) per surface; max over surfaces:")
for q in QS:
    reg = [r[q] - min([r["key"], r["c4"]] + [r[x] for x in QS]) for r in res.values()]
    print("  q %.2f: max regret %.5f  mean regret %.5f" % (q, max(reg), np.mean(reg)))
