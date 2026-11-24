#!/usr/bin/env python3
"""706b: is the strong targets' remaining gap DIRECTION that the raw donor rows already carry?

budget706 showed per-target amplitude is near its predictable optimum on held-out folds, so the strong-tier gap is
direction. This test asks whether, for the targets the donors themselves call strong, the plain (unwhitened) donor
direction transfers better than our whitened stack composition.

Surfaces: 522's 9 October-like folds (fold522.Oct, known answers asserted on load). The truth is the held-out line's own
donor rows (checked: cos 0.99 with data/donors/<line>/delta.npy on the scored genes). So the fold's two SOURCE donors'
raw delta rows are in the truth's space.
Arm: on the selected rows, the stack composition `fin` is replaced by a norm-preserving blend toward the source donors'
equal-weight unit mean d_p = Σ_j r_j/|r_j| (panel columns zeroed):
    fin'_p = |fin_p| · normalise((1 − β) fin_p/|fin_p| + β d_p/|d_p|).
Then production's key (percentile cap), q 1.75 row sharpening, render and ballast run exactly as the instrument does
(S.score). Selections are tertiles of z_src (the donors' strength rank, predictable on the day): top third, top two
thirds, all. β in 0.25, 0.5, 1.
Reports ΔX (cos_agg) and Δpds vs the base, with a paired target bootstrap of ΔX (2,000 draws, seed 706; X recomputed per
frame from the per-target x, a, d arrays, frame-averaged).
Variant (--centre): U's across-target mean (over targets with source rows) is removed before the blend, so only the
target-specific part of the donor direction enters. Diagnostics: the share of the truth's, base rows' and U's energy in
their across-target mean (the board's A/B/C truth has ~1 %: the context-mean baseline scores mse 0.986-0.992).
usage: uv run python experiments/706-strong-tier-budget/dirswap706.py [--centre]
"""
from __future__ import annotations

import csv
import json
import sys
import time
from pathlib import Path

sys.dont_write_bytecode = True
HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(ROOT / "experiments/522-october-bulk-key"))
sys.path.insert(0, str(ROOT / "experiments/438-october-stack-repricing"))
import numpy as np  # noqa: E402

import fold522 as FD  # noqa: E402

OUT = HERE / "work"
CENTRE = "--centre" in sys.argv
NONESS = "--noness" in sys.argv   # 522's essential-like (DepMap mean Chronos < -0.5 or K562-essential) removed: never
                                  # swapped, and every X / bootstrap is read on the essential-free targets only
GENES = [r["gene_name"] for r in csv.DictReader(open(ROOT / "data/challenge/gene_names.csv"))]
rng = np.random.default_rng(706)


def donor_rows(src, targets):
    D = ROOT / "data/donors" / src
    dg = (D / "genes.txt").read_text().split()
    dp = (D / "perts.txt").read_text().split()
    delta = np.load(D / "delta.npy", mmap_mode="r")
    gpos = {g: i for i, g in enumerate(dg)}
    ppos = {p: i for i, p in enumerate(dp)}
    cols = np.array([gpos.get(g, -1) for g in GENES])
    ok = cols >= 0
    R = np.zeros((len(targets), len(GENES)))
    have = np.zeros(len(targets), bool)
    for i, t in enumerate(targets):
        if t in ppos:
            R[i, ok] = np.asarray(delta[ppos[t]], float)[cols[ok]]
            have[i] = bool(np.any(R[i]))
    return R, have


def X_of(x, a, d, idx):
    return float(np.mean([x[c, idx].sum() / np.sqrt(a[c, idx].sum() * d[idx].sum()) for c in range(x.shape[0])]))


def boot(arr0, arr1, B=2000, keep=None):
    x0, a0, x1, a1, d = (np.asarray(arr0["x"], float), np.asarray(arr0["a"], float), np.asarray(arr1["x"], float),
                         np.asarray(arr1["a"], float), np.asarray(arr0["d"], float))
    live = np.flatnonzero((d > 0) & (np.ones(len(d), bool) if keep is None else keep))
    point = X_of(x1, a1, d, live) / X_of(x0, a0, d, live) - 1
    draws = []
    for _ in range(B):
        s = rng.choice(live, len(live))
        draws.append(X_of(x1, a1, d, s) / X_of(x0, a0, d, s) - 1)
    return 100 * point, [100 * float(np.quantile(draws, .025)), 100 * float(np.quantile(draws, .975))]


def main():
    res = {}
    panels = tuple(a for a in sys.argv[1:] if not a.startswith("--")) or ("efp", "efr", "real")
    for panel, held in FD.surfaces(panels):
        t0 = time.time()
        S = FD.load(panel, held, check=True)
        wp = FD.wclip(S.kv, S.carried, 0.5, FD.pct_hi(S.kv, S.carried, FD.P_PROD))
        z = FD.zsrc(S)
        s0, a0 = S.score(wp)
        HOST, ATTR = (S.F, "fST") if panel == "abc" else (S, "fin")
        NONPANEL, SCORED = HOST.nonpanel, HOST.scored
        TRUTH_FULL = HOST.truth_full
        rows, haves = zip(*(donor_rows(src, S.targets) for src in S.sources))
        donor_axis = np.zeros(len(GENES), bool)
        for R in rows:
            donor_axis |= np.any(R != 0, axis=0)
        U = np.zeros_like(getattr(HOST, ATTR))
        for R, h in zip(rows, haves):
            nr = np.linalg.norm(R, axis=1, keepdims=True)
            U += np.where(nr > 0, R / np.maximum(nr, 1e-300), 0.0)
        U[:, ~NONPANEL] = 0.0
        have = np.linalg.norm(U, axis=1) > 0
        def shared(M, rows):
            M = M[rows][:, SCORED]
            mu = M.mean(0)
            return float(len(M) * (mu @ mu) / (M * M).sum())
        diag = {"shared_truth": shared(TRUTH_FULL, have & (S.d > 0)), "shared_U": shared(U, have),
                "shared_fin": shared(getattr(HOST, ATTR), have)}
        if CENTRE:
            U[have] = U[have] - U[have].mean(0, keepdims=True)
            U[:, ~NONPANEL] = 0.0
        ESSL = np.zeros(S.n, bool)
        if NONESS:
            st = np.load(ROOT / "experiments/522-october-bulk-key/work/strength522.npz")
            pt = [str(t) for t in st[f"{panel}__targets"]]
            ii = np.asarray([pt.index(t) for t in S.targets])
            ESSL = (np.asarray(st[f"{panel}__dep_mean"], float)[ii] < -0.5) | np.asarray(st[f"{panel}__k562_essential"])[ii]
            have = have & ~ESSL
        tz = np.full(S.n, -1)
        tz[have] = FD_tiers(z[have])
        fin0 = getattr(HOST, ATTR).copy()
        rec = {"base": {"cos_agg": s0["cos_agg"], "pds": s0["pds"]}, "n_have": int(have.sum()), "diag": diag,
               "n_essential_like": int(ESSL.sum())}
        for sel_name, sel in (("top", tz == 2), ("top2", tz >= 1), ("all", have)):
            for beta in (0.25, 0.5, 1.0):
                F = fin0.copy()
                nf = np.linalg.norm(F, axis=1)
                nu = np.linalg.norm(U, axis=1)
                for i in np.flatnonzero(sel & (nf > 0) & (nu > 0)):
                    m = (1 - beta) * F[i] / nf[i] + beta * U[i] / nu[i]
                    F[i] = nf[i] * m / np.linalg.norm(m)
                setattr(HOST, ATTR, F)
                s1, a1 = S.score(wp)
                ch = F - fin0
                ce = float((ch * ch).sum())
                notscored = np.ones(len(GENES), bool)
                notscored[SCORED] = False
                axis = {"change_outside_scored": float((ch[:, notscored] ** 2).sum() / ce) if ce > 0 else 0.0,
                        "base_energy_outside_donor_axis_sel": float((fin0[sel][:, ~donor_axis] ** 2).sum() / max((fin0[sel] ** 2).sum(), 1e-300))}
                dX, ci = boot(a0, a1, keep=(~ESSL) if NONESS else None)
                rec[f"{sel_name}_b{beta:g}"] = {"dX_pct": dX, "ci95": ci, "dpds": s1["pds"] - s0["pds"], "n_sel": int(sel.sum()), **axis}
        setattr(HOST, ATTR, fin0)
        rec["seconds"] = round(time.time() - t0, 1)
        res[f"{panel}/{held}"] = rec
        line = " | ".join(f"{k} {v['dX_pct']:+.2f}% [{v['ci95'][0]:+.2f},{v['ci95'][1]:+.2f}] pds {v['dpds']:+.4f}"
                          for k, v in rec.items() if isinstance(v, dict) and "dX_pct" in v and k.endswith(("b0.5", "b1")))
        print(f"{panel}/{held} [ess-like {int(ESSL.sum())}] (X {s0['cos_agg']:.4f}; shared truth {diag['shared_truth']:.3f} U {diag['shared_U']:.3f} "
              f"fin {diag['shared_fin']:.3f}): {line} | {rec['seconds']} s", flush=True)
        tag = ("c" if CENTRE else "") + ("n" if NONESS else "") + ("-" + "-".join(panels) if panels != ("efp", "efr", "real") else "")
        (OUT / f"dirswap706{tag}.json").write_text(json.dumps(res, indent=1))


def FD_tiers(v, n=3):
    r = np.argsort(np.argsort(v, kind="stable"), kind="stable")
    return np.minimum((r * n) // len(v), n - 1)


if __name__ == "__main__":
    main()
