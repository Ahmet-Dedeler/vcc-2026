#!/usr/bin/env python3
"""625 — FIX the 16 readout modules now, before any read: gene sets on the challenge axis (gene_names.csv order).

module(g) = 2 * octile(CLEAN4 P_down) + half(A/B/C-pooled control log-CPM), over genes where production's CLEAN4 P_down is
defined (391's clean4_prior_challenge.npz, the block prior's source); -1 elsewhere. Octile edges = quantiles of P_down over
the defined genes; the expression half = the median of mean_c log1p(ctrl CPM_c) (c in A, B, C; 622's abc_cpm622.npz) over the
defined genes. Also writes the 16x16 Sylvester-Hadamard flip design (arm 1 flips all 16; arm k >= 2 flips the 8 modules
where H[k-1] = -1), the decode matrix, and two single-module ADDITIVITY check arms (modules 0 and 15). Prints the sha256s.
"""
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
genes = pd.read_csv(ROOT / "data/challenge/gene_names.csv")["gene_name"].astype(str).to_numpy()
c4 = np.load(ROOT / "experiments/391-consensus-block-signs/work/clean4_prior_challenge.npz", allow_pickle=True)
assert list(c4["genes"].astype(str)) == list(genes), "CLEAN4 file is not on the challenge axis"
P = c4["P_down"].astype(float)
abc = np.load(ROOT / "experiments/622-identity-track/work/abc_cpm622.npz", allow_pickle=True)
assert list(abc["genes"].astype(str)) == list(genes)
lexp = np.mean([np.log1p(abc[c].astype(float)) for c in ("ctxA", "ctxB", "ctxC")], 0)
ok = np.isfinite(P) & np.isfinite(lexp)
edges = np.quantile(P[ok], np.linspace(0.125, 0.875, 7))
octile = np.digitize(P, edges)
half = (lexp > np.median(lexp[ok])).astype(int)
module = np.where(ok, octile * 2 + half, -1).astype(np.int16)
H2 = np.array([[1, 1], [1, -1]])
H = H2
for _ in range(3):
    H = np.kron(H, H2)
assert H.shape == (16, 16) and np.allclose(H @ H.T, 16 * np.eye(16))
flips = np.zeros((16, 16), bool)
flips[0] = True                                   # arm 1: all modules flipped
for k in range(1, 16):
    flips[k] = H[k] == -1                           # arm k+1: 8 modules
checks = np.zeros((2, 16), bool); checks[0, 0] = True; checks[1, 15] = True
np.savez_compressed(HERE / "work/modules625.npz", genes=genes, module=module, octile_edges=edges,
                    expr_median=float(np.median(lexp[ok])), hadamard=H, flips=flips, check_flips=checks)
sha = hashlib.sha256((HERE / "work/modules625.npz").read_bytes()).hexdigest()
sizes = np.bincount(module[module >= 0], minlength=16)
info = {"modules_npz_sha256": sha, "n_defined": int(ok.sum()), "module_sizes": sizes.tolist(),
        "octile_edges": edges.tolist(), "flips_per_arm": flips.sum(1).tolist()}
(HERE / "work/modules625.json").write_text(json.dumps(info, indent=1))
print(json.dumps(info, indent=1))
