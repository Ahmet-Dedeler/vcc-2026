#!/usr/bin/env python3
"""654 helper: count zero-block-total gated genes in run496-lineage stats (stats-*.npz; key `lfc`, gated axis).
A zero total delivers lfc = log2((0 + 1e-9) / (m_ctrl + 1e-9)) <= about -30 for any gated gene (control CPM > 5), while one
count anywhere keeps lfc above about -20, so lfc < -30 <=> zero total. Production's VCC_DE_ZERO_REPAIR=-5 floors exactly these;
a harness run without --zero-repair charges them ~+30 error each.
usage: uv run python experiments/654-carrier-amplitude/zerocount654.py <stats.npz> [<stats.npz> ...]  (optional --pools f to
       also count those inside lib521's H1 pools at regime f; H1 only)"""
import sys
from pathlib import Path

import numpy as np

args = [a for a in sys.argv[1:] if not a.startswith("--")]
f = None
if "--pools" in sys.argv:
    f = float(sys.argv[sys.argv.index("--pools") + 1])
    args = [a for a in args if a != sys.argv[sys.argv.index("--pools") + 1]]
    sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "experiments/521-strength-nreal"))
    import lib521 as L  # noqa: E402
    POOL = L.pools(L.load_state(0), f)
for p in args:
    lfc = np.asarray(np.load(p)["lfc"], float)
    z = lfc < -30
    line = f"{Path(p).name}: zero-total gated cells {int(z.sum())} in {int(z.any(1).sum())} targets"
    if f is not None and z.shape == POOL.shape:
        line += f"; inside f {f:g} pools {int((z & POOL).sum())}"
    print(line)
