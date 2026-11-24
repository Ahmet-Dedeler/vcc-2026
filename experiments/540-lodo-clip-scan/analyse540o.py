"""540o summary: the LODO key's cap on 438's October-like panels (and cap 6 x the per-target arms), direction members vs the
October default (clip 0.5, 3; q 1.75), per panel, mean of the three cancer folds (realised 16).

    uv run python experiments/540-lodo-clip-scan/analyse540o.py
"""
import json
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
PANELS = ("efp", "efr", "real")
FOLDS = ("k562_gwps", "orion_hct116", "orion_hek293t")
ARMS = ("V1_0.5_4", "V1_0.5_6", "V1_0.5_8", "V1_0.5_1e+09", "V1_0.5_3 + P2g", "V1_0.5_3 + S_P2", "V1_0.5_6 + P2g", "V1_0.5_6 + S_P2")


def main():
    out = {}
    for p in PANELS:
        rows = {}
        for f in FOLDS:
            fn = HERE / "work" / f"comp540o-{p}-{f}.json"
            if not fn.exists():
                continue
            d = json.loads(fn.read_text())
            b = d["arms"]["base q1.75"]["realised"]["16"]
            rows[f] = {a: dict(r=d["arms"][a]["realised"]["16"]["cos_agg"] / b["cos_agg"] - 1,
                               dpds_lit=d["arms"][a]["realised"]["16"]["pds"] - b["pds"],
                               dpds_bd=d["arms"][a]["realised"]["16"]["pds_lam0.5"] - b["pds_lam0.5"]) for a in ARMS}
            rows[f]["keys"] = d["keys"]
        if not rows:
            continue
        mean = {a: dict(r=float(np.mean([rows[f][a]["r"] for f in rows])), pos=int(sum(rows[f][a]["r"] > 0 for f in rows)),
                        dpds_lit=float(np.mean([rows[f][a]["dpds_lit"] for f in rows])),
                        dpds_bd=float(np.mean([rows[f][a]["dpds_bd"] for f in rows]))) for a in ARMS}
        # additivity: (cap6 + X) vs cap6 + (cap3 + X), X in P2g / S_P2
        for x in ("P2g", "S_P2"):
            mean[f"additivity {x}"] = dict(combined=mean[f"V1_0.5_6 + {x}"]["r"],
                                           sum=float(np.mean([(1 + rows[f]["V1_0.5_6"]["r"]) * (1 + rows[f][f"V1_0.5_3 + {x}"]["r"]) - 1 for f in rows])))
        out[p] = dict(folds=rows, mean=mean)
        print(f"{p}: cap shares at 3 -> (4/6/8) " + str({f: {k: round(v['frac_cap'], 3) for k, v in rows[f]['keys'].items()} for f in rows}))
        for a in ARMS:
            m = mean[a]
            print(f"  {a:18s} rel cos {100 * m['r']:+.2f} % ({m['pos']}/{len(rows)}) dpds lit {m['dpds_lit']:+.4f} bd {m['dpds_bd']:+.4f}")
        for x in ("P2g", "S_P2"):
            m = mean[f"additivity {x}"]
            print(f"  additivity {x}: combined {100 * m['combined']:+.2f} % vs product of the parts {100 * m['sum']:+.2f} %")
    (HERE / "work" / "analyse540o.json").write_text(json.dumps(out, indent=1, default=float))


if __name__ == "__main__":
    main()
