"""812 — call-set composition by head size: K fixed (head + block = 2,600 prescribed), head = top-H |key| (design sign),
block = top-(2,600 - H) control CPM among the rest (CLEAN4 sign where covered, else the design sign; the q80 protection
is ignored, so both arms carry the same simplification). fid and jac analytic at f 1 on the realized truth; the base row
is the same emulation at H 400 (today), so deltas compare like with like. Does a key-ranked head buy precision x a1
(the pack's possible fid route) or lose a0?

usage: uv run python experiments/812-de-signs/head812.py   -> work/head812.json
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import surf812 as SF  # noqa: E402

INV = dict(fid=3.32482, jac=2.72081)
KTOT = SF.N_HEAD + SF.N_CPM
out = {}
for surf in ("h1test", "h1tv"):
    S = SF.load(surf, 0)
    adj, sig, call = SF.masks(S)
    rl, key = S["rl"], S["key"]
    T, G = rl.shape
    c4 = S["c4"]
    c4s = np.where(np.isfinite(c4), np.where(c4 > 0.5, -1.0, 1.0), np.nan)
    rows = {}
    for H in (0, 200, 400, 800, 1200, 1600, 2600):
        callm = np.zeros((T, G), bool)
        sg = np.zeros((T, G))
        for i in range(T):
            hk = np.argsort(-np.abs(key[i]), kind="stable")[:H]
            hk = hk[np.abs(key[i, hk]) > 0]
            hd = np.zeros(G, bool)
            hd[hk] = True
            bl = np.argsort(-np.where(hd, -np.inf, S["ref_lin"]), kind="stable")[:KTOT - hd.sum()]
            callm[i, hd] = True
            callm[i, bl] = True
            sg[i] = np.where(hd, np.sign(key[i]), np.where(np.isfinite(c4s), c4s, np.sign(key[i])))
        callm &= ~S["own"]
        fr = SF.fid_rows(callm, sg, sig, adj, rl)
        inter = (callm & sig).sum(1)
        uni = (callm | sig).sum(1)
        jr = np.where(uni == 0, 1.0, inter / np.maximum(uni, 1))
        M0 = callm & adj & ~sig
        M1 = callm & sig
        rows[H] = dict(fid=float(np.nanmean(fr)), jac=float(jr.mean()),
                       a0=round(float((np.sign(sg[M0]) == np.sign(rl[M0])).mean()), 4),
                       a1=round(float((np.sign(sg[M1]) == np.sign(rl[M1])).mean()), 4),
                       precision=round(float(M1.sum() / max((callm & adj).sum(), 1)), 4))
    b = rows[400]
    for H, r in rows.items():
        r["dfid"] = round(r["fid"] - b["fid"], 5)
        r["djac"] = round(r["jac"] - b["jac"], 5)
        r["d_overall_fidjac"] = round((r["dfid"] * INV["fid"] + r["djac"] * INV["jac"]) / 6, 5)
        print(surf, H, r, flush=True)
    out[surf] = {str(k): v for k, v in rows.items()}
json.dump(out, open(SF.HERE / "work/head812.json", "w"), indent=1)
