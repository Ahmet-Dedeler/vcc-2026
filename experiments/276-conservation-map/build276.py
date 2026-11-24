"""276 item 2 — build the per-source residual rows on the frozen common axis and the
per-target conservation statistic, with its pair count, bootstrap interval and shuffled null.

    PYTHONDONTWRITEBYTECODE=1 uv run python build276.py

Writes work/rows276.npz (residual rows per source), work/cons276.npz (the statistic and its
nulls) and work/build276.json.  Read-only outside work/.
"""
from __future__ import annotations

import json
import resource
import time
from pathlib import Path

import numpy as np

import lib276 as L

SEED = 20260906
B_PAIR = 2000
N_NULL = 200


def build_rows() -> tuple[list[str], list[str], np.ndarray, np.ndarray]:
    keys = list(L.SOURCES)
    axis = L.common_axis(keys)
    R = np.zeros((len(keys), len(L.PANEL), len(axis)), dtype=np.float32)
    A = np.zeros((len(keys), len(L.PANEL)), dtype=bool)
    for i, k in enumerate(keys):
        t0 = time.time()
        R[i], A[i] = L.load_source(k, axis)
        print(f"  {k:22s} {int(A[i].sum()):3d} panel rows  ({time.time()-t0:.1f}s)", flush=True)
    return keys, axis, R, A


def pair_cosines(keys, R, A, pairs):
    """(n_pairs, n_targets) cosine, NaN where either source lacks the target."""
    C = np.full((len(pairs), len(L.PANEL)), np.nan)
    for p, (i, j) in enumerate(pairs):
        both = A[i] & A[j]
        c = L.cosine_rows(R[i][both].astype(np.float64), R[j][both].astype(np.float64))
        C[p, both] = c
    return C


def shuffled_null(keys, R, A, pairs, n_draw=N_NULL, seed=SEED):
    """Target-shuffled pairings: source j's row replaced by a random OTHER panel target of j."""
    rng = np.random.default_rng(seed)
    nT = len(L.PANEL)
    out = np.full((n_draw, nT), np.nan)
    idx_by_src = {i: np.flatnonzero(A[i]) for i in range(len(keys))}
    for d in range(n_draw):
        acc = np.zeros(nT); cnt = np.zeros(nT)
        for (i, j) in pairs:
            both = np.flatnonzero(A[i] & A[j])
            if both.size == 0:
                continue
            pool = idx_by_src[j]
            pick = pool[rng.integers(0, pool.size, size=both.size)]
            clash = pick == both
            while clash.any():                     # never pair a target with itself
                pick[clash] = pool[rng.integers(0, pool.size, size=int(clash.sum()))]
                clash = pick == both
            c = L.cosine_rows(R[i][both].astype(np.float64), R[j][pick].astype(np.float64))
            ok = np.isfinite(c)
            acc[both[ok]] += c[ok]; cnt[both[ok]] += 1
        with np.errstate(invalid="ignore"):
            out[d] = np.where(cnt > 0, acc / np.maximum(cnt, 1), np.nan)
    return out


def boot_per_target(C, b=B_PAIR, seed=SEED):
    """Paired bootstrap over each target's own pairs."""
    rng = np.random.default_rng(seed)
    nT = C.shape[1]
    lo = np.full(nT, np.nan); hi = np.full(nT, np.nan)
    for t in range(nT):
        v = C[:, t][np.isfinite(C[:, t])]
        if v.size < 3:
            continue
        draws = v[rng.integers(0, v.size, size=(b, v.size))].mean(axis=1)
        lo[t], hi[t] = np.percentile(draws, [2.5, 97.5])
    return lo, hi


def main() -> None:
    t0 = time.time()
    keys, axis, R, A = build_rows()
    print(f"common axis {len(axis)} genes; {len(keys)} sources", flush=True)

    cross = L.pair_list(keys, cross_line_only=True)
    within = L.pair_list(keys, cross_line_only=False)
    print(f"{len(cross)} cross-line pairs, {len(within)} same-line pairs", flush=True)

    Ccross = pair_cosines(keys, R, A, cross)
    Cwithin = pair_cosines(keys, R, A, within)
    with np.errstate(invalid="ignore"):
        c_raw = np.nanmean(Ccross, axis=0)
        n_pairs = np.isfinite(Ccross).sum(axis=0)
        c_within = np.nanmean(Cwithin, axis=0)
        n_within = np.isfinite(Cwithin).sum(axis=0)
    lo, hi = boot_per_target(Ccross)
    null = shuffled_null(keys, R, A, cross)
    with np.errstate(invalid="ignore"):
        null_mean = np.nanmean(null, axis=0)
        null_lo = np.nanpercentile(null, 2.5, axis=0)
        null_hi = np.nanpercentile(null, 97.5, axis=0)

    # --- the same on the full shared axis of the three production donors: does the 2,994-gene
    #     common axis read the same direction as the donors' own 18,106-gene axis? -----------
    trio = ["k562_gwps", "orion_hct116", "orion_hek293t"]
    ax_trio = L.common_axis(trio)
    Rt, At = [], []
    for k in trio:
        r, a = L.load_source(k, ax_trio)
        Rt.append(r); At.append(a)
    diag = {}
    for (i, j) in [(0, 1), (0, 2), (1, 2)]:
        both = At[i] & At[j]
        c_big = L.cosine_rows(Rt[i][both].astype(np.float64), Rt[j][both].astype(np.float64))
        ii, jj = keys.index(trio[i]), keys.index(trio[j])
        c_small = L.cosine_rows(R[ii][both].astype(np.float64), R[jj][both].astype(np.float64))
        from scipy.stats import spearmanr
        diag[f"{trio[i]}|{trio[j]}"] = dict(
            n=int(both.sum()), n_genes_big=len(ax_trio), n_genes_common=len(axis),
            median_big=float(np.median(c_big)), median_common=float(np.median(c_small)),
            spearman=float(spearmanr(c_big, c_small).statistic),
            pearson=float(np.corrcoef(c_big, c_small)[0, 1]))
        print(f"  axis check {trio[i]}|{trio[j]}: median {diag[list(diag)[-1]]['median_big']:+.4f} "
              f"(18k) vs {diag[list(diag)[-1]]['median_common']:+.4f} (common), "
              f"rho {diag[list(diag)[-1]]['spearman']:+.3f}")

    np.savez_compressed(
        L.WORK / "rows276.npz", keys=np.asarray(keys, dtype=object),
        axis=np.asarray(axis, dtype=object), R=R, A=A,
        targets=np.asarray(L.PANEL, dtype=object))
    np.savez_compressed(
        L.WORK / "cons276.npz",
        targets=np.asarray(L.PANEL, dtype=object), keys=np.asarray(keys, dtype=object),
        pairs_cross=np.asarray(cross), pairs_within=np.asarray(within),
        Ccross=Ccross.astype(np.float32), Cwithin=Cwithin.astype(np.float32),
        c_raw=c_raw, n_pairs=n_pairs, c_within=c_within, n_within=n_within,
        boot_lo=lo, boot_hi=hi, null_mean=null_mean, null_lo=null_lo, null_hi=null_hi)

    ok = np.isfinite(c_raw)
    summary = dict(
        n_sources=len(keys), n_genes_common=len(axis),
        n_cross_pairs=len(cross), n_within_pairs=len(within),
        targets_with_score=int(ok.sum()),
        pairs_per_target=dict(min=int(n_pairs[ok].min()), median=float(np.median(n_pairs[ok])),
                              max=int(n_pairs.max())),
        c_raw=dict(mean=float(np.nanmean(c_raw)), median=float(np.nanmedian(c_raw)),
                   p10=float(np.nanpercentile(c_raw, 10)), p90=float(np.nanpercentile(c_raw, 90)),
                   share_le0=float(np.nanmean(c_raw[ok] <= 0))),
        c_within=dict(mean=float(np.nanmean(c_within)), median=float(np.nanmedian(c_within)),
                      n=int(np.isfinite(c_within).sum())),
        null=dict(mean=float(np.nanmean(null_mean)), median=float(np.nanmedian(null_mean)),
                  sd_over_targets=float(np.nanstd(null_mean)),
                  mean_band=[float(np.nanmean(null_lo)), float(np.nanmean(null_hi))]),
        axis_check=diag, seconds=time.time() - t0,
        peak_rss_gib=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 2**30)
    (L.WORK / "build276.json").write_text(json.dumps(summary, indent=1))
    print(json.dumps(summary, indent=1))


if __name__ == "__main__":
    main()
