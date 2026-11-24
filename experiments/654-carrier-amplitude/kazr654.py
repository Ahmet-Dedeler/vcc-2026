"""654 A1 structural KA: base_zr s<k> == base s<k> (q, p, lfc, sums) on every target without a zero-total gated gene in base
(lfc < -30 marks one), and differs only inside those targets. usage: kazr654.py <seed> -> KA PASS / KA FAIL (rc 0 / 1)"""
import glob, sys
from pathlib import Path
import numpy as np
D = Path(__file__).resolve().parent / "work/h1_cpm"
s = int(sys.argv[1])
b = np.load(D / f"stats-base_s{s}_bps_prodem_mde0.1cap0.1min0.3coarse_dnk250zs.npz") if (D / f"stats-base_s{s}_bps_prodem_mde0.1cap0.1min0.3coarse_dnk250zs.npz").exists() else np.load(Path(__file__).resolve().parents[2] / f"experiments/630-rank-neutral-up/work/h1/stats-base_s{s}_bps_prodem_mde0.1cap0.1min0.3coarse_dnk250zs.npz")
zf = [f for f in glob.glob(str(D / f"stats-base_s{s}_*zr-5*_dnk250zs.npz"))]
assert len(zf) == 1, zf
z = np.load(zf[0])
zt = (np.asarray(b["lfc"]) < -30).any(1)
ok = True
for k in ("q", "p", "lfc", "sums"):
    x, y = np.asarray(z[k]), np.asarray(b[k])
    same_rows = np.array([np.array_equal(x[i], y[i], equal_nan=True) if x.dtype.kind == "f" else np.array_equal(x[i], y[i]) for i in range(x.shape[0])])
    print(f"  {k}: rows equal {int(same_rows.sum())}/{len(same_rows)}; zero-total targets {int(zt.sum())}; unequal outside them {int((~same_rows & ~zt).sum())}")
    ok &= bool((~same_rows & ~zt).sum() == 0)
print(f"  file {Path(zf[0]).name}")
print("KA PASS" if ok else "KA FAIL")
sys.exit(0 if ok else 1)
