"""435 known answers and the reading. KA1: prod emitter + coarse MDE (no zs) == 426's k_emit426 coarse run, stats
bit-equal. KA2: prod emitter + plain MDE + zs == 429's MDE + zs run, bit-equal. Then per seed: coarse+zs vs coarse
alone (426's runs), all six raw members from the metrics files."""
import json
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
W = ROOT / "experiments/435-zs-on-coarse/work/h1"; W426 = ROOT / "experiments/426-mde-magnitude/work/h1"
W429 = ROOT / "experiments/429-zs-production/work/h1"

def eq(a, b):
    A, B = np.load(a), np.load(b)
    bad = [k for k in sorted(set(A.files) & set(B.files)) if not np.array_equal(A[k], B[k])]
    return (not bad) and set(A.files) == set(B.files), bad

ok = True
for nm, a, b in (("KA1 prod coarse == k426 coarse", W / "stats-base_s0_bps_prodem_mde0.1cap0.1min0.3coarse.npz",
                  W426 / "stats-base_s0_bps_mde0.1cap0.1min0.3coarse.npz"),
                 ("KA2 prod MDE+zs == 429 MDE+zs", W / "stats-base_s0_bps_prodem_mde0.1cap0.1min0.3_dnk250zs.npz",
                  W429 / "stats-base_s0_bps_mde0.1cap0.1min0.3_dnk250zs.npz")):
    if not a.exists():
        print(nm, "PENDING", a.name); ok = False; continue
    e, bad = eq(a, b); ok &= e
    print(nm, "PASS (bit-equal)" if e else f"FAIL {bad}")
MEM = ["pds", "mse", "nmae", "fid", "reach", "jac"]
SLOPE = {"pds": 2.2138, "mse": -1.051, "nmae": -1.6532, "fid": 3.3645, "reach": 1.1292, "jac": 2.7263}

def raw(p):
    m = json.loads(p.read_text())
    r = m.get("raw", m)
    return {k: r[k] for k in MEM if k in r}

rows = []
for s in (0, 1, 2):
    a = W / f"metrics-base_s{s}_bps_prodem_mde0.1cap0.1min0.3coarse_dnk250zs.json"; b = W426 / f"metrics-base_s{s}_bps_mde0.1cap0.1min0.3coarse.json"
    if a.exists() and b.exists():
        ra, rb = raw(a), raw(b)
        d = {k: ra[k] - rb[k] for k in MEM if k in ra and k in rb}
        ov = sum(SLOPE[k] * d[k] for k in d) / 6
        rows.append((s, d, ov))
        print(f"seed {s}: " + " ".join(f"{k} {d[k]:+.5f}" for k in d) + f"  overall {ov:+.5f}")
if rows:
    dn = np.array([r[1]["nmae"] for r in rows]); ov = np.array([r[2] for r in rows])
    print(f"zs on coarse: nmae {dn.mean():+.5f} ± {dn.std(ddof=1) if len(dn) > 1 else 0:.5f}; overall {ov.mean():+.5f} (429: zs on plain MDE nmae -0.00179 ± 0.00007)")
