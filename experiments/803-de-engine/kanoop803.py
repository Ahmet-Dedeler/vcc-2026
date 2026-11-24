#!/usr/bin/env python3
"""803 follow-up KA: the sharpening wrapper with q None (arms shnoop*) must equal the plain E0 runs (g1pre / zs, collected
in outmerged/) array for array in members.npz, every surface and seed. Writes work/out/ka_noop803.json; exit 1 on a miss.
    uv run python experiments/803-de-engine/kanoop803.py <new out dir> <outmerged dir>"""
import json
import sys
from pathlib import Path

import numpy as np

new, old = Path(sys.argv[1]), Path(sys.argv[2])
PAIRS = (tuple(tuple(x.split(":")) for x in sys.argv[3].split(",")) if len(sys.argv) > 3
         else (("shnooppre", "g1pre"), ("shnoopzs", "zs")))
OUTNAME = sys.argv[4] if len(sys.argv) > 4 else "ka_noop803.json"
rec, ok = {}, True
for sf in ("h1test", "h1tv"):
    for a, b in PAIRS:
        for s in (0, 1, 2):
            fa, fb = new / sf / a / f"s{s}/members.npz", old / sf / b / f"s{s}/members.npz"
            if not (fa.exists() and fb.exists()):
                rec[f"{sf}|{a}|s{s}"] = "missing"
                ok = False
                continue
            za, zb = np.load(fa), np.load(fb)
            bad, last_bit = [], {}
            for k in zb.files:
                x, y = za[k], zb[k]
                if x.shape == y.shape and (np.array_equal(x, y, equal_nan=True) if x.dtype.kind in "fc" else np.array_equal(x, y)):
                    continue
                # cell-eval2's own means are not bit-deterministic on identical inputs (599 hashseed599): tolerate <= 1e-12 there
                if k == "ce_raw":
                    dx = max(abs(v - json.loads(str(y))[m]) for m, v in json.loads(str(x)).items())
                elif k.startswith("ce_") and x.shape == y.shape:
                    dx = float(np.nanmax(np.abs(x - y)))
                else:
                    bad.append(k)
                    continue
                if dx <= 1e-12:
                    last_bit[k] = dx
                else:
                    bad.append(k)
            rec[f"{sf}|{a}|s{s}"] = ("bit-equal" if not last_bit else {"bit-equal except cell-eval2 last bits": last_bit}) if not bad else {"differs": bad}
            ok &= not bad
rec["pass"] = ok
Path(new / OUTNAME).write_text(json.dumps(rec, indent=1))
print(json.dumps(rec, indent=1))
sys.exit(0 if ok else 1)
