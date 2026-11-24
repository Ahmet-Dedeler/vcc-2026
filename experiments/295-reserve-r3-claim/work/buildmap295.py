"""295 item 2 — the MAP side.  Reads only the 28 conservation-map banks (279's sources); no reserve
file is opened here.

For every one of the 174 distinct reserve genes it stores that gene's UNIT residual row in each
bank that carries it, on 279's frozen 2,739-gene axis, built by `build279.residual_rows` /
`load_unit` with the FULL 17,106-gene held list so that the centring vector is the map's own
(DEFINITION-276 §2), then subset to the reserve genes.

KA-4 (the construction check): the pairwise cosines of those stored rows must reproduce
`held279.npz` `Ccross` on the reserve genes.

Writes work/map295.npz, work/ka_map295.json.
"""
from __future__ import annotations

import csv
import json
import resource
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path("~/Code/virtual-cell-challenge-2026")
HERE = ROOT / "experiments/295-reserve-r3-claim"
WORK = HERE / "work"
sys.path.insert(0, str(ROOT / "experiments/276-conservation-map"))
sys.path.insert(0, str(ROOT / "experiments/279-znf-state-dependence"))
import build279 as B  # noqa: E402
import inv279 as I  # noqa: E402

HELD279 = ROOT / "experiments/279-znf-state-dependence/work/held279.npz"


def main() -> None:
    t0 = time.time()
    z = np.load(HELD279, allow_pickle=True)
    genes_held = [str(g) for g in z["genes"]]
    axis = [str(g) for g in z["axis"]]
    keys = [str(k) for k in z["keys"]]
    lines = [str(l) for l in z["lines"]]
    Ccross = z["Ccross"]
    pairs_cross = z["pairs_cross"]

    with open(WORK / "reserve_genes295.csv") as fh:
        rows = list(csv.DictReader(fh))
    rgenes = sorted({r["gene"] for r in rows})
    gpos = {g: i for i, g in enumerate(rgenes)}
    held_pos = {g: i for i, g in enumerate(genes_held)}
    rheld = np.asarray([held_pos[g] for g in rgenes], dtype=np.int64)   # index into the 17,106
    print(f"{len(rgenes)} distinct reserve genes; axis {len(axis)}; {len(keys)} banks", flush=True)

    U = np.zeros((len(keys), len(rgenes), len(axis)), np.float32)
    A = np.zeros((len(keys), len(rgenes)), bool)
    for q, k in enumerate(keys):
        t1 = time.time()
        rel = I.CANDIDATES[k][0]
        idx, u, _ = B.load_unit(rel, axis, genes_held)      # idx = positions into genes_held
        pos = {int(i): j for j, i in enumerate(idx)}
        for g in rgenes:
            j = pos.get(held_pos[g])
            if j is not None:
                U[q, gpos[g]] = u[j]
                A[q, gpos[g]] = True
        del idx, u
        print(f"  {k:22s} {int(A[q].sum()):4d}/{len(rgenes)} reserve rows "
              f"({time.time()-t1:.1f}s)", flush=True)

    # ---- KA-4: reproduce Ccross on the reserve genes from the stored rows -------------------
    worst, n_chk = 0.0, 0
    for p, (i, j) in enumerate(pairs_cross):
        both = A[i] & A[j]
        if not both.any():
            continue
        c = (U[i][both] * U[j][both]).sum(axis=1)
        ref = Ccross[p][rheld[both]]
        ok = np.isfinite(ref)
        if ok.any():
            worst = max(worst, float(np.max(np.abs(c[ok] - ref[ok]))))
            n_chk += int(ok.sum())
    print(f"KA-4  {n_chk} (pair, reserve gene) cosines vs held279 Ccross: max|diff| {worst:.2e} "
          f"-> {'PASS' if worst < 1e-5 else 'FAIL'}", flush=True)

    # also confirm every reserve gene's Ccross availability matches A
    nb = A.sum(axis=0)
    np.savez_compressed(WORK / "map295.npz", U=U, avail=A,
                        genes=np.asarray(rgenes, dtype=object),
                        keys=np.asarray(keys, dtype=object),
                        lines=np.asarray(lines, dtype=object),
                        axis=np.asarray(axis, dtype=object),
                        held_pos=rheld)
    out = dict(n_genes=len(rgenes), n_banks=len(keys), n_axis=len(axis),
               ka4_n_checked=n_chk, ka4_max_abs_diff=worst, ka4_pass=bool(worst < 1e-5),
               banks_per_gene=dict(min=int(nb.min()), median=float(np.median(nb)),
                                   max=int(nb.max())),
               seconds=time.time() - t0,
               peak_rss_gib=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 2**30)
    (WORK / "ka_map295.json").write_text(json.dumps(out, indent=1))
    print(json.dumps(out, indent=1))
    if not out["ka4_pass"]:
        raise SystemExit("KA-4 FAILED — abort")


if __name__ == "__main__":
    main()
