"""812 step 1b — is sign predictability a property of WHICH genes we call? CLEAN4-sign accuracy on every non-significant
adjudicable gene (called or not, head excluded) by control CPM x CLEAN4 confidence, on both DE-H1 surfaces; and an
analytic emulation of re-composing the 2,200-gene block by CLEAN4 confidence instead of control CPM (same K, head
unchanged, every block gene delivered with the CLEAN4 sign; base = today's block under the same delivery assumption).

usage: uv run python experiments/812-de-signs/pool812.py     -> work/pool812.json
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import surf812 as SF  # noqa: E402

INV = dict(fid=3.32482, jac=2.72081)
out = {}
for surf in ("h1test", "h1tv"):
    S = SF.load(surf, 0)
    adj, sig, call = SF.masks(S)
    rl = S["rl"]
    T, G = rl.shape
    c4 = S["c4"]
    cov = np.isfinite(c4)
    conf = np.where(cov, np.abs(c4 - 0.5), 0.0)
    c4s = np.where(c4 > 0.5, -1.0, 1.0)
    cpm = S["cpm"]
    ok = (np.sign(rl) == c4s[None, :])
    M = adj & ~sig & ~S["head"] & cov[None, :]
    tab = {}
    for lo, hi in ((5, 10), (10, 20), (20, 50), (50, 100), (100, 300), (300, 1e9)):
        for clo, chi in ((0, 0.05), (0.05, 0.1), (0.1, 0.2), (0.2, 0.51)):
            gm = (cpm >= lo) & (cpm < hi) & (conf >= clo) & (conf < chi)
            m = M & gm[None, :]
            n = int(m.sum())
            tab[f"cpm{lo}-{hi:g}|conf{clo}-{chi}"] = dict(n=n, genes=int(gm.sum()), acc=round(float(ok[m].mean()), 4) if n else None,
                                                         inblock=round(float((S["block"] & gm[None, :]).sum() / max(T * gm.sum(), 1)), 3))
    R = dict(table=tab)
    # analytic block re-composition (sign = CLEAN4 where covered, else the emitted sign of today's arm; head = emitted)
    base_sign = np.where(S["head"], S["ol"], np.where(cov[None, :], c4s[None, :], S["ol"]))
    call0 = (S["head"] | S["block"]) & ~S["own"]

    def score(callm, sgn):
        fr = SF.fid_rows(callm, sgn, sig, adj, rl)
        inter = (callm & sig).sum(1)
        uni = (callm | sig).sum(1)
        jr = np.where(uni == 0, 1.0, inter / np.maximum(uni, 1))
        M0 = callm & adj & ~sig
        a0 = float((np.sign(sgn[M0]) == np.sign(rl[M0])).mean())
        return float(np.nanmean(fr)), float(jr.mean()), a0
    f0, j0, a00 = score(call0, base_sign)
    R["base"] = dict(fid=f0, jac=j0, a0=a00)
    rows = {}
    for floor in (10, 20, 50):
        for mix in (0.25, 0.5, 1.0):   # share of the block re-chosen by confidence
            blk = np.zeros((T, G), bool)
            nsw = int(round(mix * SF.N_CPM))
            for i in range(T):
                cand = ~S["head"][i] & (cpm >= floor) & ~S["own"][i]
                keepn = SF.N_CPM - nsw
                order_cpm = np.argsort(-np.where(cand, S["ref_lin"], -np.inf), kind="stable")
                kept = order_cpm[:keepn]
                blk[i, kept] = True
                sc = np.where(cand & ~blk[i], conf + 1e-9 * np.log1p(cpm), -np.inf)
                blk[i, np.argsort(-sc, kind="stable")[:nsw]] = True
            callm = (S["head"] | blk) & ~S["own"]
            f1, j1, a01 = score(callm, base_sign)
            rows[f"floor{floor}|mix{mix}"] = dict(dfid=round(f1 - f0, 5), djac=round(j1 - j0, 5), a0=round(a01, 4),
                                                  d_overall_fidjac=round(((f1 - f0) * INV["fid"] + (j1 - j0) * INV["jac"]) / 6, 5))
    R["recompose"] = rows
    out[surf] = R
    print(surf, "base", R["base"])
    for k, v in rows.items():
        print("  ", k, v)
json.dump(out, open(SF.HERE / "work/pool812.json", "w"), indent=1)
