"""493g summary: direction members of each arm vs the q 1.75 base under the October LODO key, per panel/fold (realised 16).
rel cos = cos_agg(arm)/cos_agg(base) - 1; dpds_bd = board-emulated pds (lam 0.5); dpds_lit = literal pds.

    uv run python experiments/493-pertarget-q-october/analyse493g.py
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
PANELS = ("efp", "efr", "real")
FOLDS = ("k562_gwps", "orion_hct116", "orion_hek293t")
ARMS = ("P2g", "P2 (no guard)", "S_P2", "L30g", "CTRL P2g permuted")


def main():
    out = {}
    for p in PANELS:
        rows = {}
        for f in FOLDS:
            fn = HERE / "work" / f"comp493g-{p}-{f}.json"
            if not fn.exists():
                continue
            d = json.loads(fn.read_text())
            b = d["arms"]["base q1.75"]["realised"]["16"]
            rows[f] = dict(frac_floor=d["frac_floor"], terciles=d["tercile_counts"])
            for a in ARMS:
                x = d["arms"][a]["realised"]["16"]
                rows[f][a] = dict(r=x["cos_agg"] / b["cos_agg"] - 1, dpds_bd=x["pds_lam0.5"] - b["pds_lam0.5"],
                                  dpds_lit=x["pds"] - b["pds"], dK=x["K_mean"] - b["K_mean"])
        if not rows:
            continue
        summ = {}
        for a in ARMS:
            v = [rows[f][a] for f in rows]
            summ[a] = dict(r=float(np.mean([x["r"] for x in v])), pos=int(sum(x["r"] > 0 for x in v)), n=len(v),
                           dpds_bd=float(np.mean([x["dpds_bd"] for x in v])), dpds_lit=float(np.mean([x["dpds_lit"] for x in v])),
                           dK=float(np.mean([x["dK"] for x in v])))
        out[p] = dict(folds=rows, mean=summ)
        print(f"{p}: floor shares {[round(rows[f]['frac_floor'], 3) for f in rows]}; P2g terciles {[rows[f]['terciles']['P2g'] for f in rows]}")
        for a, s in summ.items():
            print(f"  {a:18s} rel cos {100 * s['r']:+.2f} % ({s['pos']}/{s['n']}) dpds bd {s['dpds_bd']:+.4f} lit {s['dpds_lit']:+.4f} "
                  f"dK {s['dK']:+.2f}")
    (HERE / "work" / "analyse493g.json").write_text(json.dumps(out, indent=1, default=float))


if __name__ == "__main__":
    main()
