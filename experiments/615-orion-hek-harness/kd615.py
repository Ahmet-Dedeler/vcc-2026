#!/usr/bin/env python3
"""615 QC: CRISPRi knockdown of each harness target's OWN gene in its cells vs the control pool (labels sanity).
Streams the truth in row blocks and reads only the 100 target columns. -> work/kd615.json"""
import json
from pathlib import Path
import h5py
import numpy as np
W = Path(__file__).resolve().parent / "work"
with h5py.File(W / "hek-truth.h5", "r") as h:
    lab = h["obs/target_gene"][:].astype(str)
    genes = h["var/_index"][:].astype(str)
    ip = h["X/indptr"][:]
    targets = sorted(set(lab) - {"non-targeting"})
    gi = {g: i for i, g in enumerate(genes)}
    cols = np.array([gi.get(t, -1) for t in targets])
    colpos = {c: k for k, c in enumerate(cols) if c >= 0}
    n = len(lab)
    V = np.zeros((n, len(targets)))
    lib = np.zeros(n)
    for r0 in range(0, n, 5000):
        r1 = min(n, r0 + 5000)
        a0, a1 = int(ip[r0]), int(ip[r1])
        d = h["X/data"][a0:a1].astype(np.float64)
        ix = h["X/indices"][a0:a1]
        rows = np.repeat(np.arange(r0, r1), np.diff(ip[r0:r1 + 1]))
        np.add.at(lib, rows, d)
        m = np.isin(ix, cols[cols >= 0])
        for c in np.unique(ix[m]):
            k = colpos[int(c)]
            sel = m & (ix == c)
            np.add.at(V[:, k], rows[sel], d[sel])
cpm = 1e6 * V / lib[:, None]
ctrl = lab == "non-targeting"
out = {}
for k, t in enumerate(targets):
    if cols[k] < 0:
        out[t] = None
        continue
    own, base = float(cpm[lab == t, k].mean()), float(cpm[ctrl, k].mean())
    out[t] = {"own_cpm": own, "ctrl_cpm": base, "remaining": own / base if base > 0 else None}
rem = np.array([v["remaining"] for v in out.values() if v and v["remaining"] is not None])
summ = {"targets": len(targets), "with_ctrl_expression": int(len(rem)), "remaining_median": float(np.median(rem)),
        "remaining_q75": float(np.quantile(rem, 0.75)), "n_remaining_gt_0.7": int((rem > 0.7).sum())}
print(json.dumps(summ))
(W / "kd615.json").write_text(json.dumps({"summary": summ, "per_target": out}, indent=1))
