"""480 resume-2 audit, folds: recompute RESULTS.md §4.1's fold rows with fresh code (no common480 / fold480 imports).

For each cancer fold, the truth gate comes from gates480.py's table (scorer rules; its own known answers are in §1),
donors exclude the held-out line, and the constructions follow §4.1's definitions:
  cons      mean of the available donors' log2((P_pert + 1 CPM) / (P_ctrl + 1 CPM)), 0 where none
  mde_cons  sign(cons) x min(0.1 / sqrt(cnt), 0.1) for cnt >= 0.3 (cnt = control CPM x median library / 1e6)
  a_analog  mde_cons + 0.575 x cons

    uv run python audit480b.py     # -> work/audit480b.json
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
TS = 5e4
DONORS = ["k562_gwps", "orion_hct116", "orion_hek293t"]

QUOTED = {  # RESULTS.md §4.1 (nmae) and §4.1's calibration column (a_analog − cons@1)
    "orion_hct116": dict(scored=127, mde_cons=0.9844, cons1=0.9356, a_analog=0.9417, d=0.0061),
    "orion_hek293t": dict(scored=125, mde_cons=0.9971, cons1=0.9905, a_analog=0.9851, d=-0.0054),
    "k562_gwps": dict(scored=71, mde_cons=0.9753, cons1=0.8731, a_analog=0.8967, d=0.0236),
}


def donor_log2(name, targets, genes, pc=1.0):
    base = ROOT / "data/donors" / name
    prow = {p.strip(): i for i, p in enumerate(open(base / "perts.txt"))}
    gcol: dict[str, int] = {}
    for i, g in enumerate(open(base / "genes.txt")):
        gcol.setdefault(g.strip(), i)
    col = np.array([gcol.get(g, -1) for g in genes])
    have = col >= 0
    delta = np.load(base / "delta.npy", mmap_mode="r")
    b = np.load(base / "b_ctrl.npy").astype(float)[col[have]]
    c = np.expm1(b) / TS * 1e6
    out = np.full((len(targets), len(genes)), np.nan)
    for k, t in enumerate(targets):
        if t in prow:
            p = np.maximum(np.expm1(b + np.asarray(delta[prow[t]], float)[col[have]]) / TS * 1e6, 0.0)
            with np.errstate(divide="ignore", invalid="ignore"):
                v = np.log2((p + pc) / (c + pc))
            out[k, have] = np.where(np.isfinite(v), v, np.nan)
    return out


res, bad = {}, 0
for fold, want in QUOTED.items():
    z = np.load(HERE / f"work/gate-{fold}.npz", allow_pickle=True)
    targets, genes = z["targets"].astype(str), z["genes"].astype(str)
    q, lfc = z["q"].astype(float), z["lfc"].astype(float)
    gate = (q < 0.05) & np.isfinite(lfc) & (lfc != 0) & (genes[None, :] != targets[:, None])
    scored = gate.sum(1) >= 10
    gate &= scored[:, None]
    Tz = np.where(gate, lfc, 0.0)
    den = np.abs(Tz).sum(1)

    def rows(P):
        out = np.full(len(targets), np.nan)
        out[scored] = (np.abs(np.nan_to_num(P) - Tz) * gate).sum(1)[scored] / den[scored]
        return out

    with np.errstate(all="ignore"):
        cons = np.nan_to_num(np.nanmean(np.stack([donor_log2(d, targets, genes) for d in DONORS if d != fold]), 0))
    cnt = z["ctrl_cpm"].astype(float) * float(z["lib_med"]) / 1e6
    add = np.where(cnt >= 0.3, np.minimum(0.1 / np.sqrt(np.maximum(cnt, 1e-3)), 0.1), 0.0)[None, :]
    mde_cons = np.sign(cons) * add
    got = dict(scored=int(scored.sum()), zero=float(np.nanmean(rows(np.zeros_like(Tz)))),
               mde_cons=float(np.nanmean(rows(mde_cons))), cons1=float(np.nanmean(rows(cons))),
               a_analog=float(np.nanmean(rows(mde_cons + 0.575 * cons))))
    dd = rows(mde_cons + 0.575 * cons) - rows(cons)
    dd = dd[np.isfinite(dd)]
    bs = dd[np.random.default_rng(0).integers(0, dd.size, (4000, dd.size))].mean(1)
    got["d"] = float(dd.mean())
    got["d_ci"] = [float(np.quantile(bs, 0.025)), float(np.quantile(bs, 0.975))]
    for k, v in want.items():
        tol = 0 if k == "scored" else 6e-5
        ok = abs(got[k] - v) <= tol
        bad += not ok
        print(f"{'ok ' if ok else 'BAD'} {fold:14s} {k:9s} quoted {v:+.4f} recomputed {got[k]:+.5f}")
    print(f"    {fold:14s} zero {got['zero']:.6f}  Δ CI [{got['d_ci'][0]:+.4f}, {got['d_ci'][1]:+.4f}]")
    res[fold] = got
json.dump(res, open(HERE / "work/audit480b.json", "w"), indent=1)
print("fold audit:", "all rows reproduce" if not bad else f"{bad} rows differ")
sys.exit(1 if bad else 0)
