"""496z: does the hook stack on the zs rider? Per seed: (hook + zs) vs (zs alone, 435's coarse + dn_k250_zs runs), all six
raw members, overall-equiv with 435's board slopes; and the determinism check (the re-run pre-pass stats equal the first
run's, moved to work/h1/bak/). Additivity: the hook's delta on the zs base vs its delta on the plain base (496's runs).
usage: uv run python experiments/496-uncalled-sign-hook/stack496z.py"""
import json
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
W = ROOT / "experiments/496-uncalled-sign-hook/work/h1"; W435 = ROOT / "experiments/435-zs-on-coarse/work/h1"
MEM = ["pds", "mse", "nmae", "fid", "reach", "jac"]
SLOPE = {"pds": 2.2138, "mse": -1.051, "nmae": -1.6532, "fid": 3.3645, "reach": 1.1292, "jac": 2.7263}
raw = lambda p: json.loads(p.read_text())["raw"]
rows = []
for s in (0, 1, 2):
    hz = W / f"metrics-base_s{s}_bps_hookem_mde0.1cap0.1min0.3coarseas0.05_dnk250zs.json"
    if not hz.exists():
        print(f"seed {s}: pending"); continue
    pre, bak = W / f"stats-base_s{s}_bps_hookem_mde0.1cap0.1min0.3coarseas0.05.npz", W / f"bak/stats-base_s{s}_bps_hookem_mde0.1cap0.1min0.3coarseas0.05.npz"
    A, B = np.load(pre), np.load(bak)
    det = set(A.files) == set(B.files) and all(np.array_equal(A[k], B[k]) for k in A.files)
    z = raw(W435 / f"metrics-base_s{s}_bps_prodem_mde0.1cap0.1min0.3coarse_dnk250zs.json")
    b = raw(W435 / f"metrics-base_s{s}_bps_prodem_mde0.1cap0.1min0.3coarse.json")
    h = raw(W / f"bak/metrics-base_s{s}_bps_hookem_mde0.1cap0.1min0.3coarseas0.05.json")
    hzr = raw(hz)
    d_on_zs = {m: hzr[m] - z[m] for m in MEM}
    d_on_base = {m: h[m] - b[m] for m in MEM}
    ov = lambda d: sum(SLOPE[m] * d[m] for m in MEM) / 6
    rows.append((s, d_on_zs, d_on_base))
    print(f"seed {s}: determinism {'PASS' if det else 'FAIL'} | hook on zs: " + " ".join(f"{m} {d_on_zs[m]:+.5f}" for m in MEM)
          + f" ov {ov(d_on_zs):+.5f} | hook on base: nmae {d_on_base['nmae']:+.5f} reach {d_on_base['reach']:+.5f} ov {ov(d_on_base):+.5f}")
if rows:
    nz = np.array([r[1]["nmae"] for r in rows]); nb = np.array([r[2]["nmae"] for r in rows])
    print(f"nmae: hook on zs {nz.mean():+.5f} (sd {nz.std(ddof=1) if len(nz) > 1 else 0:.5f}) vs hook on base {nb.mean():+.5f}: "
          f"{'ADDITIVE (within 20 %)' if abs(nz.mean() - nb.mean()) <= 0.2 * abs(nb.mean()) else 'NOT additive'}")
