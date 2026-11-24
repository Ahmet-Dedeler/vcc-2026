"""810 -- the A-only instrument read's known-answer gate (light; MBs, seconds). For each of the 18 A folds in <run>/bulk:
price810's own KA keys (parts sum == F.fin, sgu0 == inc x5) present and < 1e-12, and inc / 707 == 801's ka-aonly arrays
(ka804 unmodified, experiments/801-lab-folds/work/ka-aonly/bulk/) on x, a, pds, pdsl (< 1e-12): 14 keys per fold; a missing key
fails. Also decomp810's compose KAs from <run>/parts/<fold>-comp.json. Writes <run>/ka-aonly810.json; aborts on any failure.

    uv run python experiments/810-strong-direction/aonly810.py experiments/810-strong-direction/work/aonly
"""
import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
REF = ROOT / "experiments/801-lab-folds/work/ka-aonly/bulk"
NEED = ["parts sum == F.fin"] + [f"sgu0 == inc: {k}" for k in ("x", "a", "pds", "pdsl", "K")]
run = Path(sys.argv[1])
stems = [f"A{k}__{h}" for k in range(1, 7) for h in ("k562_gwps", "orion_hct116", "orion_hek293t")]
out, bad = {}, []
for stem in stems:
    meta = run / "meta" / f"{stem}.json"
    if not meta.exists():
        bad.append((stem, "meta missing", None))
        continue
    ka = {k: json.loads(meta.read_text())["ka"].get(k) for k in NEED}
    comp = json.loads((run / "parts" / f"{stem}-comp.json").read_text())["ka"]
    for k in ("parts sum == reference fin", "h1c1 == reference fPP"):
        ka[f"decomp810 {k}"] = comp.get(k)
    for a in ("inc", "707"):
        r, z = np.load(REF / stem / f"{a}.npz"), np.load(run / "bulk" / stem / f"{a}.npz")
        assert list(r["targets"]) == list(z["targets"]), stem
        for k in ("x", "a", "pds", "pdsl"):
            ka[f"{a} == 801 ka-aonly {k}"] = float(np.max(np.abs(np.asarray(z[k], float) - np.asarray(r[k], float))))
    out[stem] = ka
    bad += [(stem, k, v) for k, v in ka.items() if v is None or not v < 1e-12]
(run / "ka-aonly810.json").write_text(json.dumps(out, indent=1))
n = sum(len(v) for v in out.values())
worst = max((v for d in out.values() for v in d.values() if v is not None), default=None)
print(f"{len(out)} folds, {n} keys, worst {worst}")
assert not bad, f"KA FAILED: {bad}"
print("KA gate PASS")
