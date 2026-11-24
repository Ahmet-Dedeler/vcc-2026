#!/usr/bin/env python3
"""625 Amendment 2 — the 8-module design (CLEAN4 octiles), before any read.

Same CLEAN4 file and octile edges as modules625.npz (the expression half is dropped). Hadamard column j <-> octile COLMAP[j];
column 0 (the all-flip row's column, decoded sd 1.20σ vs 0.66σ for the others; coordinator's arithmetic) holds octile 3,
the weakest prior (median |P - .5| 0.005). Arms: 8 Sylvester rows (row 0 flips all 8; rows 1-7 flip the 4 columns with
H = -1) + 2 single-octile additivity checks (octiles 0 and 7). Writes work/modules625_v8.npz and prints its sha256.
"""
import hashlib
import json
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
z = np.load(HERE / "work/modules625.npz", allow_pickle=True)
module16 = z["module"].astype(int)
octile = np.where(module16 >= 0, module16 // 2, -1).astype(np.int16)
COLMAP = np.array([3, 0, 1, 2, 4, 5, 6, 7])           # column j holds octile COLMAP[j]
col_of_octile = np.argsort(COLMAP)                    # octile k sits in column col_of_octile[k]
column = np.where(octile >= 0, col_of_octile[np.maximum(octile, 0)], -1).astype(np.int16)
H2 = np.array([[1, 1], [1, -1]]); H = np.kron(np.kron(H2, H2), H2)
assert np.allclose(H @ H.T, 8 * np.eye(8)) and (H[:, 0] == 1).all()
flips = np.zeros((8, 8), bool); flips[0] = True
for k in range(1, 8):
    flips[k] = H[k] == -1
checks = np.zeros((2, 8), bool); checks[0, col_of_octile[0]] = True; checks[1, col_of_octile[7]] = True
np.savez_compressed(HERE / "work/modules625_v8.npz", genes=z["genes"], module=column, octile=octile, colmap=COLMAP,
                    octile_edges=z["octile_edges"], hadamard=H, flips=flips, check_flips=checks)
sha = hashlib.sha256((HERE / "work/modules625_v8.npz").read_bytes()).hexdigest()
info = {"modules625_v8_sha256": sha, "colmap_column_to_octile": COLMAP.tolist(),
        "column_sizes": np.bincount(column[column >= 0], minlength=8).tolist(), "flips_per_arm": flips.sum(1).tolist(),
        "check_columns": [int(np.flatnonzero(c)[0]) for c in checks], "check_octiles": [0, 7]}
(HERE / "work/modules625_v8.json").write_text(json.dumps(info, indent=1))
print(json.dumps(info, indent=1))
