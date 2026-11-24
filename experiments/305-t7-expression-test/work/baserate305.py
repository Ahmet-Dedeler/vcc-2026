"""305 step 5 -- POST-HOC negative comparator (declared in RESULTS.md as a deviation added AFTER the
outcomes were seen): the base rate of the frozen four-criterion conjunction over ARBITRARY genes.

The frozen test asks whether a candidate gene is unusually repressed under its construct.  It cannot be
read without knowing how often an arbitrary gene of the same expression class satisfies the identical
conjunction under the identical construct.  This script applies C1-C4 verbatim to EVERY testable gene
(control CP50k >= 0.5) of every tested construct, and to the expression-matched band around each
candidate.  It can only weaken a positive reading, never create one.

    PYTHONDONTWRITEBYTECODE=1 uv run python baserate305.py
"""
from __future__ import annotations

import csv
import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import lib305 as X  # noqa: E402
from expr305 import MAG, MIN_CP50K, P_RAW, Split  # noqa: E402

BAND = 2.0        # expression-matched band: control CP50k within a factor of 2 of the candidate's


def main() -> None:
    P = X.L.Logger(X.HERE / "baserate305.log")
    SP = {s: Split(s) for s in X.SPLITS}
    pairs = list(csv.DictReader(open(X.HERE / "candidates305.csv")))

    out = []
    for r in pairs:
        c, g, grp = r["construct"], r["candidate_gene"], r["group"]
        s = r["split"] if grp != "ntc" else "train"
        sp = SP[s]
        cA, cB, cAB, loo = sp.refs(c)
        pAB, pA, pB = X.L.prof(cAB), X.L.prof(cA), X.L.prof(cB)
        e_pool = sp.Pmat[sp.cidx[c]] - pAB
        e0 = X.L.prof(sp.half[(c, 0)]) - pA
        e1 = X.L.prof(sp.half[(c, 1)]) - pB
        E = sp.Pmat - pAB[None, :]                                   # every construct, every gene
        keep = [sp.cidx[c2] for c2 in sp.targeting if c2 != c]       # the same null population as C4
        N = len(keep)
        rank = (E[keep] <= e_pool[None, :]).sum(0) + 1               # 1 = most negative of all
        p_raw = (1.0 + (rank - 1)) / (1.0 + N)
        cp = X.L.TS * cAB / cAB.sum()
        testable = cp >= MIN_CP50K
        conj = (e_pool < 0) & (e_pool <= MAG) & (e0 < 0) & (e1 < 0) & (p_raw <= P_RAW) & testable
        j = sp.gpos[g]
        band = testable & (cp >= cp[j] / BAND) & (cp <= cp[j] * BAND)
        row = dict(pair_id=r["pair_id"], group=grp, construct=c, candidate_gene=g, split=s,
                   candidate_cp50k=float(cp[j]), candidate_meets=bool(conj[j]),
                   n_testable_genes=int(testable.sum()), n_meeting=int(conj.sum()),
                   base_rate=float(conj.sum() / max(1, testable.sum())),
                   n_band=int(band.sum()), n_meeting_band=int((conj & band).sum()),
                   base_rate_band=float((conj & band).sum() / max(1, band.sum())),
                   n_rank1=int(((rank == 1) & testable).sum()),
                   n_rank1_mag=int(((rank == 1) & testable & (e_pool <= MAG)).sum()),
                   candidate_e=float(e_pool[j]),
                   e_quantile_among_meeting=(float((np.sort(e_pool[conj]) <= e_pool[j]).mean())
                                             if conj.any() else float("nan")),
                   strongest_meeting=[[str(sp.genes[k]), float(e_pool[k])]
                                      for k in np.argsort(np.where(conj, e_pool, np.inf))[:5]
                                      if conj[k]])
        out.append(row)
        P(f"  {grp:8s} {c[:26]:26s} {g:9s} CP50k {cp[j]:7.3f} | candidate meets conjunction: "
          f"{'YES' if conj[j] else 'no ':3s} | base rate over {int(testable.sum()):5d} testable genes "
          f"{conj.sum():4d} ({conj.sum() / max(1, testable.sum()):.3%}); in the x{BAND:g} expression band "
          f"{int((conj & band).sum())}/{int(band.sum())} ({(conj & band).sum() / max(1, band.sum()):.3%})")

    P("\nper-construct summary (one row per tested construct):")
    seen = {}
    for r in out:
        seen.setdefault((r["construct"], r["split"]), r)
    for (c, s), r in seen.items():
        P(f"  {c[:34]:34s} {s:11s} {r['n_meeting']:4d}/{r['n_testable_genes']} genes meet C1-C4 "
          f"({r['base_rate']:.3%}); rank-1 genes {r['n_rank1']}, of them |e|>=0.10: {r['n_rank1_mag']}")

    X.dump(dict(band=BAND, rows=out), X.HERE / "baserate305.json")
    P.save()


if __name__ == "__main__":
    main()
