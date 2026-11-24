"""810 -- the LABC instrument read's known-answer gate and summary (EC2 only, after price810 ran on the three LABC folds; the Mac
was under a local-job hold). For each LABC fold: price810's own KA keys (parts sum == F.fin, sgu0 == inc x5) must all be present
and < 1e-12, and inc / 707 must equal 801's ka-labc arrays (ka804 unmodified, `experiments/801-lab-folds/work/ka-labc/bulk/`) on
x, a, pds, pdsl (< 1e-12): 14 keys per fold, a missing key fails. Writes <run>/ka-labc810.json and aborts on any failure.

    uv run python experiments/810-strong-direction/labc810.py experiments/810-strong-direction/work/labc
"""
import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
REF = ROOT / "experiments/801-lab-folds/work/ka-labc/bulk"
NEED = ["parts sum == F.fin"] + [f"sgu0 == inc: {k}" for k in ("x", "a", "pds", "pdsl", "K")]
run = Path(sys.argv[1])
out, bad = {}, []
for held in ("k562_gwps", "orion_hct116", "orion_hek293t"):
    stem = f"LABC__{held}"
    ka = {k: json.loads((run / "meta" / f"{stem}.json").read_text())["ka"].get(k) for k in NEED}
    for a in ("inc", "707"):
        r, z = np.load(REF / stem / f"{a}.npz"), np.load(run / "bulk" / stem / f"{a}.npz")
        assert list(r["targets"]) == list(z["targets"]), stem
        for k in ("x", "a", "pds", "pdsl"):
            ka[f"{a} == 801 ka-labc {k}"] = float(np.max(np.abs(np.asarray(z[k], float) - np.asarray(r[k], float))))
    out[stem] = ka
    bad += [(stem, k, v) for k, v in ka.items() if v is None or not v < 1e-12]
    print(stem, len(ka), "keys, worst", max(v for v in ka.values() if v is not None), flush=True)
(run / "ka-labc810.json").write_text(json.dumps(out, indent=1))
assert not bad, f"KA FAILED: {bad}"
print(f"KA gate PASS: {sum(len(v) for v in out.values())} keys")
