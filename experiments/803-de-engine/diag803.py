#!/usr/bin/env python3
"""803 diagnosis of the Linux-vs-Mac UPRN gap: on the MAC's own dumped h1-test s0 blocks (630/work/blocks, 435's
construction), compute the pre-pass and zs-pass stats on Linux (evaluate323.block_stats + vcc2026.rank_mag.postpass, as
rnup630 does), write them as the src zs stats for rnup630 (so its own zs KA compares Linux with Linux), and record how they
differ from the Mac's 435 stats (same counts in, so any difference is the stats arithmetic itself).
    uv run python experiments/803-de-engine/diag803.py <src dir with blocks/> <mac zs stats.npz> <out json>
"""
import json
import re
import sys
from pathlib import Path

import numpy as np
import scipy.sparse as sp

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "experiments/400-uncalled-knn"))
import run400 as H  # noqa: E402
from vcc2026 import rank_mag as ZS  # noqa: E402

E = H.E
TAG = "base_s0_bps_prodem_mde0.1cap0.1min0.3coarse"
src, mac, outj = Path(sys.argv[1]), np.load(sys.argv[2]), Path(sys.argv[3])
r, c, _ = E.load_reference()
files = sorted((src / "blocks").glob(f"{TAG}_*.npz"), key=lambda f: int(re.search(r"_(\d+)\.npz$", f.name).group(1)))
assert len(files) == 85
spec = ZS.parse_spec("dn_k250_zs")
acc = {k: [] for k in ("sums", "jk", "p", "q", "lfc")}
for f in files:
    z = np.load(f)
    X = sp.csr_matrix((z["data"], z["indices"], z["indptr"]), shape=tuple(z["shape"])).toarray().astype(np.int64)
    s, j, p, q, l = E.block_stats(X, c)
    Xz, _ = ZS.postpass(X, int(z["n_car"]), c, l, spec)
    s, j, p, q, l = E.block_stats(Xz, c)
    for k, v in zip(acc, (s, j, p, q, l)):
        acc[k].append(v)
st = {k: np.asarray(v) for k, v in acc.items()}
np.savez(src / f"stats-{TAG}_dnk250zs.npz", targets=r["targets"], **st)
rec = {}
for k in ("sums", "jk", "p", "q", "lfc"):
    a, b = np.asarray(st[k], float), np.asarray(mac[k], float)
    fin = np.isfinite(a) & np.isfinite(b)
    rec[k] = dict(equal=bool(np.array_equal(a, b, equal_nan=True)), n_diff=int((a[fin] != b[fin]).sum()),
                  max_abs=float(np.max(np.abs(a[fin] - b[fin]), initial=0.0)),
                  max_rel=float(np.max(np.abs(a[fin] - b[fin]) / np.maximum(np.abs(b[fin]), 1e-300), initial=0.0)),
                  nan_pattern_equal=bool(np.array_equal(np.isnan(a), np.isnan(b))))
rec["calls_diff"] = int(((st["q"] < 0.05) != (np.asarray(mac["q"]) < 0.05)).sum())
outj.parent.mkdir(parents=True, exist_ok=True)
outj.write_text(json.dumps(rec, indent=1))
print(json.dumps(rec))
