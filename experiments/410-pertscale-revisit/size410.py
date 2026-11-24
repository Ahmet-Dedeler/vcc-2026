"""410 — re-size the per-target magnitude lever on the four PLOLO folds (deep pseudobulk truth).

158 sized it on the thinned Orion harnesses (aggregate cosine 0.026 / 0.014, noisy truth) and closed it
on an nmae cost measured through the welded DE channel. This script re-measures only the BULK side
(the quantity `expr_mse` pools): the aggregate transfer cosine of the whitened prediction, the oracle
per-target scale ceiling, and what build-time keys capture of it. pds is exactly invariant to a
strictly positive per-row scale, so nothing here can move it.

Per fold (held library = truth, the other libraries = sources, fold-held whitening operator):
    P   = whiten(pre)                        the production-shaped prediction (a15-era 'pre' cached by 162)
    T   = truth                              the held library's pseudobulk effect
    agg = <P,T>_F / (|P|_F |T|_F) on the fold's scored columns
    oracle scale s_p = max(<P_p,T_p>, 0)/|P_p|^2, sse_oracle = 1 - sum s_p <P_p,T_p> / |T|^2
Keys (all from the SOURCE donors only):
    cos   mean pairwise cosine of source donor rows
    share mean pairwise <d_i,d_j> / |mean_i d_i|^2       (shared-energy shrinkage factor)
    wshare the same with the fold operator applied to the donor rows
Weights w_p = clip(key, lo, hi)^g or rank-based as in 097; global scale re-optimised (free for mse).
Output: work/size410.json
"""
from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "experiments/162-retrodicting-validator"))
import plolo  # noqa: E402

OUT = ROOT / "experiments/410-pertscale-revisit/work"
OUT.mkdir(parents=True, exist_ok=True)
GAMMA, LAM = 0.65, 0.0

OFFICIAL = plolo.OFFICIAL
OFF_INDEX = {g: i for i, g in enumerate(OFFICIAL)}


def donor_rows(name: str, targets) -> tuple[np.ndarray, np.ndarray]:
    """Raw donor delta rows for `targets` on the official axis (zeros off the donor's axis)."""
    path = plolo.DONOR_PATHS[name]
    genes = [l.strip() for l in (path / "genes.txt").read_text().splitlines() if l.strip()]
    perts = [l.strip() for l in (path / "perts.txt").read_text().splitlines() if l.strip()]
    pidx = {p: i for i, p in enumerate(perts)}
    delta = np.load(path / "delta.npy", mmap_mode="r")
    cols_src, cols_off = [], []
    for j, g in enumerate(genes):
        k = OFF_INDEX.get(g)
        if k is not None:
            cols_src.append(j)
            cols_off.append(k)
    cols_src = np.asarray(cols_src)
    cols_off = np.asarray(cols_off)
    out = np.zeros((len(targets), len(OFFICIAL)), dtype=np.float64)
    have = np.zeros(len(targets), dtype=bool)
    for r, t in enumerate(targets):
        i = pidx.get(str(t))
        if i is None:
            continue
        out[r, cols_off] = np.asarray(delta[i], dtype=np.float64)[cols_src]
        have[r] = True
    return out, have


def agg_cos(P, T):
    return float((P * T).sum() / np.sqrt((P * P).sum() * (T * T).sum()))


def sse_weighted(P, T, w):
    x = (w * (P * T).sum(1)).sum()
    a = (w * w * (P * P).sum(1)).sum()
    d = (T * T).sum()
    return 1.0 - x * x / (a * d)


def rank01(x):
    return (np.argsort(np.argsort(x)) + 0.5) / len(x)


def main():
    res = {}
    for held, fold in plolo.FOLDS.items():
        t0 = time.time()
        z = np.load(ROOT / f"experiments/162-retrodicting-validator/work/fold-{held}.npz", allow_pickle=True)
        pre = np.asarray(z["pre"], dtype=np.float64)
        truth = np.asarray(z["truth"], dtype=np.float64)
        scored = np.asarray(z["scored"])
        wcols = np.asarray(z["wcols"])
        targets = [str(t) for t in z["targets"]]
        op = plolo.Operator(fold["op"])
        P_full = op.apply(pre, wcols, GAMMA, LAM)
        P = P_full[:, scored]
        T = truth[:, scored]
        live = (np.abs(T).sum(1) > 0) & (np.abs(P).sum(1) > 0)
        P, T = P[live], T[live]
        tg = [t for t, k in zip(targets, live) if k]
        n = len(tg)
        base_agg = agg_cos(P, T)
        base_sse = 1.0 - base_agg ** 2
        pt = (P * T).sum(1)
        pp = (P * P).sum(1)
        tt = (T * T).sum(1)
        s_or = np.maximum(pt, 0.0) / pp
        sse_or = 1.0 - (s_or * pt).sum() / tt.sum()
        cos_p = pt / np.sqrt(pp * tt)
        # keys from the SOURCE donors only
        rows, haves = [], []
        for s in fold["sources"]:
            r, h = donor_rows(s, tg)
            rows.append(r)
            haves.append(h)
        # restrict donor rows to the scored columns for the keys
        R = [r[:, scored] for r in rows]
        H = np.vstack(haves)
        nsrc = len(R)
        cos_key = np.zeros(n)
        share_key = np.zeros(n)
        ip_key = np.zeros(n)
        carried = H.sum(0) >= 2
        for p in range(n):
            have = [k for k in range(nsrc) if H[k, p]]
            if len(have) < 2:
                continue
            cs, ips = [], []
            for a_ in range(len(have)):
                for b_ in range(a_ + 1, len(have)):
                    u, v = R[have[a_]][p], R[have[b_]][p]
                    nu, nv = np.linalg.norm(u), np.linalg.norm(v)
                    cs.append(float(u @ v) / (nu * nv) if nu * nv > 0 else 0.0)
                    ips.append(float(u @ v))
            m = np.mean([R[k][p] for k in have], axis=0)
            cos_key[p] = np.mean(cs)
            ip_key[p] = np.mean(ips)
            share_key[p] = np.mean(ips) / max(float(m @ m), 1e-12)
        # whitened shared-energy key: <W d_i, d_j> / |W mean|^2 on the operator axis
        wshare_key = np.zeros(n)
        Wrows = []
        for k in range(nsrc):
            full = np.zeros((n, len(OFFICIAL)))
            full[:] = rows[k]
            Wk = op.apply(full, wcols, GAMMA, LAM) if (np.abs(full[:, wcols]).sum(1) > 0).sum() > n // 2 else None
            Wrows.append(None if Wk is None else Wk[:, scored])
        for p in range(n):
            have = [k for k in range(nsrc) if H[k, p] and Wrows[k] is not None]
            if len(have) < 2:
                continue
            ips = []
            for a_ in range(len(have)):
                for b_ in range(a_ + 1, len(have)):
                    ips.append(0.5 * (float(Wrows[have[a_]][p] @ R[have[b_]][p]) + float(Wrows[have[b_]][p] @ R[have[a_]][p])))
            m = np.mean([Wrows[k][p] for k in have], axis=0)
            wshare_key[p] = np.mean(ips) / max(float(m @ m), 1e-12)
        keys = {"cos": cos_key, "share": share_key, "wshare": wshare_key, "ip": ip_key,
                "rownorm": np.sqrt(pp), "truthnorm_ORACLE": np.sqrt(tt)}
        out = dict(n=n, carried=int(carried.sum()), base_agg=base_agg, base_sse=base_sse,
                   oracle_sse=float(sse_or), oracle_agg_equiv=float(np.sqrt(1 - sse_or)),
                   frac_sstar_le0=float((pt <= 0).mean()), keys={})
        for name, key in keys.items():
            ks = {}
            ks["spearman_sstar"] = float(np.corrcoef(rank01(key), rank01(pt / pp))[0, 1])
            ks["spearman_cos"] = float(np.corrcoef(rank01(key), rank01(cos_p))[0, 1])
            grid = {}
            for g in (0.5, 1.0, 2.0):
                for f in (0.0, 0.1, 0.3, 0.5):
                    w = np.maximum(rank01(key), f) ** g
                    grid[f"rank g{g} f{f}"] = float(np.sqrt(1 - sse_weighted(P, T, w)) / base_agg - 1)
            if name in ("share", "wshare", "cos"):
                for lo in (0.0, 0.05, 0.1, 0.2):
                    for g in (0.5, 1.0, 2.0):
                        w = np.clip(key, lo, 1.0) ** g
                        grid[f"clip lo{lo} g{g}"] = float(np.sqrt(1 - sse_weighted(P, T, w)) / base_agg - 1)
            ks["rel_agg_gain"] = grid
            best = max(grid, key=grid.get)
            ks["best"] = [best, grid[best]]
            out["keys"][name] = ks
        out["oracle_rel_agg_gain"] = float(np.sqrt(1 - sse_or) / base_agg - 1)
        out["seconds"] = time.time() - t0
        res[held] = out
        print(json.dumps({held: {k: v for k, v in out.items() if k != "keys"}}), flush=True)
        for name, ks in out["keys"].items():
            print(f"  {name:18s} rho(s*) {ks['spearman_sstar']:+.3f} rho(cos) {ks['spearman_cos']:+.3f} "
                  f"best {ks['best'][0]} {ks['best'][1]:+.4f}", flush=True)
        (OUT / "size410.json").write_text(json.dumps(res, indent=1))


if __name__ == "__main__":
    main()
