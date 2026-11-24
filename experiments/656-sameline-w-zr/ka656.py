"""656 structural KA: each zero-repair seed-0 arm == 450's stored non-repaired seed-0 arm (q, p, lfc, sums) on every target
without a non-own zero total (lfc < -30 off the own-target diagonal). usage: ka656.py <w> <p> -> KA PASS / KA FAIL (rc 0/1)"""
import glob, sys
from pathlib import Path
import numpy as np
ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "experiments/521-strength-nreal"))
import lib521 as L  # noqa: E402
w, p = sys.argv[1], sys.argv[2]
d = f"h1c_w{w}_clean4_f{p}s0"
mine = glob.glob(str(ROOT / f"experiments/656-sameline-w-zr/work/{d}/stats-base_s0_bps*zr-5*.npz"))
ref = ROOT / f"experiments/450-same-line-programme/work/{d}/stats-base_s0_bps_prodem_mde0.1cap0.1min0.3coarse.npz"
assert len(mine) == 1 and ref.exists(), (mine, ref)
a, b = np.load(mine[0]), np.load(ref)
own = L.load_state(0)["own"]
hit = ((np.asarray(b["lfc"]) < -30) & ~own).any(1)
ok = True
for k in ("q", "p", "lfc", "sums"):
    x, y = np.asarray(a[k]), np.asarray(b[k])
    eq = np.array([np.array_equal(x[i], y[i], equal_nan=True) if x.dtype.kind == "f" else np.array_equal(x[i], y[i]) for i in range(x.shape[0])])
    bad = int((~eq & ~hit).sum())
    print(f"  {d} {k}: rows equal {int(eq.sum())}/{len(eq)}; targets with a non-own zero total {int(hit.sum())}; unequal outside them {bad}")
    ok &= bad == 0
print("KA PASS" if ok else "KA FAIL")
sys.exit(0 if ok else 1)
