"""654 known answers: run654's `cpm base <seed>` stats == the stored production-construction stats, bit for bit
(q, p, lfc, sums; pre-zs and zs). seed 0 vs 537's work/h1_cpm, seed 1 vs 630's work/h1 (the base 631 decomposed).
usage: uv run python experiments/654-carrier-amplitude/ka654.py <seed>   -> prints KA PASS / KA FAIL, rc 0 / 1"""
import sys
from pathlib import Path
import numpy as np
ROOT = Path(__file__).resolve().parents[2]
seed = int(sys.argv[1])
ref_dir = ROOT / ("experiments/537-arm-pipeline/work/h1_cpm" if seed == 0 else "experiments/630-rank-neutral-up/work/h1")
mine = Path(__file__).resolve().parent / "work/h1_cpm"
ok = True
for suf in ("", "_dnk250zs"):
    f = f"stats-base_s{seed}_bps_prodem_mde0.1cap0.1min0.3coarse{suf}.npz"
    a, b = np.load(mine / f, allow_pickle=True), np.load(ref_dir / f, allow_pickle=True)
    for k in ("q", "p", "lfc", "sums"):
        x, y = np.asarray(a[k]), np.asarray(b[k])
        same = x.shape == y.shape and np.array_equal(x, y, equal_nan=True) if x.dtype.kind == "f" else np.array_equal(x, y)
        print(f"  {f} {k}: {'equal' if same else 'DIFFERS'}")
        ok &= bool(same)
print("KA PASS" if ok else "KA FAIL")
sys.exit(0 if ok else 1)
