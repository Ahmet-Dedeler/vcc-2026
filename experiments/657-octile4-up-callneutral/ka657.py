"""657 known answers. usage: ka657.py E <seed>  -> run657e tvc4h (c4h) == 618's hk_c4h stats (q, p, lfc, sums; pre-zs and zs)
                             ka657.py R <seed>  -> the copy (rnup657, unsuffixed 630) == the ORIGINAL rnup630 (rnup657orig) on the same c4h dumps"""
import sys
from pathlib import Path
import numpy as np
ROOT = Path(__file__).resolve().parents[2]; HERE = Path(__file__).resolve().parent
X = ROOT / "experiments/618-h1-trainval-surface/work"
kind, s = sys.argv[1], int(sys.argv[2])
tag = f"base_s{s}_bps_zr-5_prodem_mde0.1cap0.1min0.3coarse"
pairs = ([(HERE / f"work/c4h/stats-{tag}{suf}.npz", X / f"hk_c4h/stats-{tag}{suf}.npz") for suf in ("", "_dnk250zs")] if kind == "E"
         else [(HERE / f"work/stacks/stats-rnup_c4h_s{s}_b_k100_q20r3p10.npz", HERE / f"work/stacks_orig/stats-rnup_c4h_s{s}_b_k100_q20r3p10.npz")])
ok = True
for a, b in pairs:
    za, zb = np.load(a), np.load(b)
    for k in ("q", "p", "lfc", "sums"):
        x, y = np.asarray(za[k]), np.asarray(zb[k])
        e = x.shape == y.shape and (np.array_equal(x, y, equal_nan=True) if x.dtype.kind == "f" else np.array_equal(x, y))
        print(f"  KA-{kind} s{s} {a.name} {k}: {'equal' if e else 'DIFFERS'}"); ok &= bool(e)
print("KA PASS" if ok else "KA FAIL"); sys.exit(0 if ok else 1)
