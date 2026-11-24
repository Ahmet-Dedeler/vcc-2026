"""308 -- the retrodiction row: `STACK` (h1union's stack = the UN-deflated CD4 at f = 0.34 + the H1 union at g = 0.28)
on each cancer fold, so 261's pricing recipe can retrodict the BR249 board step inside THIS experiment.

    PYTHONDONTWRITEBYTECODE=1 uv run python stack308.py <fold>
"""
from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
WORK = HERE / "work"
sys.path.insert(0, str(HERE))
import lib308 as Z  # noqa: E402

Y, L, X = Z.Y, Z.Y.L, Z.Y.X


def main(held: str) -> None:
    t0 = time.time()
    st = Y.stage(held)
    ident = X.Identity(st)
    rows = {}
    for name in ("cd4svd5", "STACK"):
        if name == "cd4svd5":
            parts, aud = Y.arm_parts(st, st["whitened"])
        else:
            fr, pu = Y.plolo.solve_fraction(st["whitened"], st["cd4"], Y.CD4_F)
            from compose218 import solve_fraction_carrying
            fx, px = solve_fraction_carrying(st["whitened"], st["union"], Y.H1_G)
            parts = [(st["cd4"], float(fr)), (st["union"], float(fx))]
            aud = dict(cd4_fraction=float(fr), cd4_per_unit=float(pu), union_fraction=float(fx), union_per_unit=float(px),
                       note="STACK: un-deflated CD4 at the staged fraction")
        st_arm = dict(st, whitened=st["whitened"])
        final = L.compose(st_arm, parts)
        e = L.score(st_arm, final)
        m = ident.measure(final)
        rows[name] = dict(fold=held, arm=name, **{k: v for k, v in e.items() if not k.startswith("_")},
                          nmae_wald89=m["nmae_wald89"], top89_precision=m["top89_precision"],
                          delivered_norm=m["delivered_norm"], audit=aud)
        np.savez_compressed(WORK / f"stackid308-{held}-{name}.npz", U=m["U"], M=m["M"], D=m["D"])
        print(f"[{held}/{name}] cos_agg {e['cos_agg']:.6f} pds {e['pds']:.6f} nmae89 {m['nmae_wald89']:.5f}", flush=True)
    (WORK / f"stack308-{held}.json").write_text(json.dumps(dict(rows=list(rows.values())), indent=2))
    print(f"[{held}] {time.time()-t0:.0f}s, peak {Y.peak_gib():.2f} GiB", flush=True)


if __name__ == "__main__":
    main(sys.argv[1])
