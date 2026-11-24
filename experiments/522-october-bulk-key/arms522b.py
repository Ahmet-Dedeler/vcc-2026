"""522b (EXPLORATORY, written after arms522's read) — the rider's feature: src only vs src + CD4 (the panel's own CD4
extraction; October-portable). On A/B/C the CD4-augmented strength was the best stand-alone form (+3.24 %), on the
October-like panels the worst of the strength forms. Arms on the percentile cap (P_PROD) and on cap 6:
    pct x Scd4 th (th 0.25, 0.5), cap 6 x Scd4 0.5;  reference rows pct, pct x S0.25 (== arms522's, asserted |d| < 1e-12)

    uv run python experiments/522-october-bulk-key/arms522b.py
"""
from __future__ import annotations

import sys
from pathlib import Path

sys.dont_write_bytecode = True
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import json  # noqa: E402
import time  # noqa: E402

import numpy as np  # noqa: E402

import fold522 as FD  # noqa: E402

WORK = HERE / "work"


def main():
    out = {}
    for p, f in FD.surfaces(("abc", "efp", "efr", "real")):
        t0 = time.time()
        S = FD.load(p, f)
        z, zc = FD.zsrc(S), FD.zsrc(S, extra=("cd4_robust",))
        kv, car = S.kv, S.carried
        w3 = FD.wclip(kv, car, 0.5, 3.0)
        w6 = FD.wclip(kv, car, 0.5, 6.0)
        wp = FD.wclip(kv, car, 0.5, FD.pct_hi(kv, car, FD.P_PROD))
        ref = json.loads((WORK / f"arms522-{p}-{f}.json").read_text())["arms"]
        A = {"cap3": w3, "pctPROD": wp, "pctPRODxS0.25": wp * np.exp(0.25 * z), "pctPRODxScd4_0.25": wp * np.exp(0.25 * zc),
             "pctPRODxScd4_0.5": wp * np.exp(0.5 * zc), "cap6xScd4_0.5": w6 * np.exp(0.5 * zc)}
        rec = {}
        base = None
        for nm, w in A.items():
            s, _ = S.score(w)
            base = base or s
            rec[nm] = dict(cos_agg=s["cos_agg"], pds=s["pds"], r=s["cos_agg"] / base["cos_agg"] - 1, dpds=s["pds"] - base["pds"])
            if nm in ref:
                d = abs(s["cos_agg"] - ref[nm]["cos_agg"])
                assert d < 1e-12, f"KA {p}/{f} {nm}: {d}"
        out[f"{p}/{f}"] = rec
        print(f"[{p}/{f}] " + " ".join(f"{nm} {100 * v['r']:+.2f}" for nm, v in rec.items() if nm != "cap3")
              + f" ({time.time() - t0:.0f}s)", flush=True)
        (WORK / "arms522b.json").write_text(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()
