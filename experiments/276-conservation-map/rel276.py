"""276 items 2b + 3 — the reliability-corrected conservation score, and the split-half
reliability OF THE SCORE ITSELF.

Sources with a measured per-target split-half reliability, all read read-only:
  * K562-GWPS / Orion HCT116 / Orion HEK293T — 271's batch-parity half rows
    (`experiments/271-.../work/half/rows-<line>.npz`: `A`, `B` on the official 18,533 axis,
    `r_half_raw` = 267's split-half residual cosine per target).
  * DLD-1 — its two low-MOI replicates as separate compact-donor banks
    (`experiments/223-.../work/tables/dld1_low{1,2}`); `dld1_lowmoi` is their mean.

Writes work/rel276.npz and work/rel276.json.

    PYTHONDONTWRITEBYTECODE=1 uv run python rel276.py
"""
from __future__ import annotations

import json
import resource
from pathlib import Path

import numpy as np
from scipy.stats import spearmanr

import lib276 as L

SEED = 20260906
R_FLOOR = 0.02          # DEFINITION-276 §3: pairs dropped when either reliability <= this
HALF271 = L.ROOT / "experiments/271-donor-reliability-estimation/work/half"

HALFED = {                      # key -> (line, how the halves were made)
    "k562_gwps": ("K562", "gemgroup parity"),
    "orion_hct116": ("HCT116", "batch parity"),
    "orion_hek293t": ("HEK293T", "batch parity"),
    "dld1_lowmoi": ("DLD-1", "two low-MOI replicates"),
}


def sb(r: np.ndarray) -> np.ndarray:
    """Spearman-Brown: half-depth reliability -> full-depth reliability."""
    with np.errstate(invalid="ignore", divide="ignore"):
        return 2.0 * r / (1.0 + r)


def load_halves(axis: list[str]):
    """Return per source: (rowsA, rowsB, availA, availB, r_half) on `axis`, residual-centred."""
    out = {}
    gi = {g: i for i, g in enumerate(L.OFFICIAL)}
    cols_off = np.asarray([gi[g] for g in axis], dtype=np.int64)
    for key in ("k562_gwps", "orion_hct116", "orion_hek293t"):
        z = np.load(HALF271 / f"rows-{key}.npz", allow_pickle=True)
        tg = [str(t) for t in z["targets"]]
        assert tg == L.PANEL, f"{key}: target order differs from the panel"
        A = np.asarray(z["A"], dtype=np.float64)[:, cols_off]
        B = np.asarray(z["B"], dtype=np.float64)[:, cols_off]
        aA = np.asarray(z["avail_A"], dtype=bool)
        aB = np.asarray(z["avail_B"], dtype=bool)
        A[aA] -= A[aA].mean(axis=0, keepdims=True)
        B[aB] -= B[aB].mean(axis=0, keepdims=True)
        A[~aA] = 0.0; B[~aB] = 0.0
        out[key] = dict(A=A, B=B, availA=aA, availB=aB,
                        r_half=np.asarray(z["r_half_raw"], dtype=np.float64),
                        real=np.asarray(z["real"], dtype=np.float64)[:, cols_off],
                        avail=np.asarray(z["avail"], dtype=bool))
    # DLD-1: the two replicates
    A, aA = L.load_source_path("experiments/223-silent-zero-audit/work/tables/dld1_low1", axis)
    B, aB = L.load_source_path("experiments/223-silent-zero-audit/work/tables/dld1_low2", axis)
    r = np.full(len(L.PANEL), np.nan)
    both = aA & aB
    r[both] = L.cosine_rows(A[both].astype(np.float64), B[both].astype(np.float64))
    out["dld1_lowmoi"] = dict(A=A.astype(np.float64), B=B.astype(np.float64),
                              availA=aA, availB=aB, r_half=r, real=None, avail=aA & aB)
    return out


def submap_score(rows: dict[str, np.ndarray], avail: dict[str, np.ndarray], keys: list[str]):
    """mean cross-line pair cosine over the given rows; returns (score, n_pairs, per-pair matrix)."""
    lines = [HALFED[k][0] if k in HALFED else L.SOURCES[k][1] for k in keys]
    pairs = [(i, j) for i in range(len(keys)) for j in range(i + 1, len(keys))
             if lines[i] != lines[j]]
    C = np.full((len(pairs), len(L.PANEL)), np.nan)
    for p, (i, j) in enumerate(pairs):
        ki, kj = keys[i], keys[j]
        both = avail[ki] & avail[kj]
        if not both.any():
            continue
        C[p, both] = L.cosine_rows(rows[ki][both], rows[kj][both])
    with np.errstate(invalid="ignore"):
        s = np.nanmean(C, axis=0)
    return s, np.isfinite(C).sum(axis=0), C, pairs


def boot_spearman(x, y, b=2000, seed=SEED):
    ok = np.isfinite(x) & np.isfinite(y)
    x, y = x[ok], y[ok]
    if x.size < 5:
        return np.nan, (np.nan, np.nan), int(x.size)
    rng = np.random.default_rng(seed)
    r0 = spearmanr(x, y).statistic
    draws = np.empty(b)
    for i in range(b):
        s = rng.integers(0, x.size, x.size)
        draws[i] = spearmanr(x[s], y[s]).statistic
    return float(r0), tuple(float(v) for v in np.percentile(draws, [2.5, 97.5])), int(x.size)


def main() -> None:
    z = np.load(L.WORK / "cons276.npz", allow_pickle=True)
    keys = [str(k) for k in z["keys"]]
    Ccross = np.asarray(z["Ccross"], dtype=np.float64)
    pairs = np.asarray(z["pairs_cross"])
    axis_full = [str(g) for g in np.load(L.WORK / "rows276.npz", allow_pickle=True)["axis"]]

    # ---- the half-map's own axis: the four halved sources only -------------------------------
    hkeys = ["k562_gwps", "orion_hct116", "orion_hek293t", "dld1_lowmoi"]
    axis_h = L.common_axis_paths([L.SOURCES[k][0] for k in hkeys])
    print(f"half-map axis: {len(axis_h)} genes")
    H = load_halves(axis_h)

    # KA-4: 271's reliability medians
    ka4 = {k: float(np.nanmedian(H[k]["r_half"][np.isfinite(H[k]["r_half"])]))
           for k in ("k562_gwps", "orion_hct116", "orion_hek293t")}
    print("KA-4 reliability medians:", {k: round(v, 4) for k, v in ka4.items()})
    ka_dld = float(np.nanmedian(H["dld1_lowmoi"]["r_half"][np.isfinite(H["dld1_lowmoi"]["r_half"])]))
    print(f"  DLD-1 replicate-pair reliability median {ka_dld:+.4f} "
          f"(n={int(np.isfinite(H['dld1_lowmoi']['r_half']).sum())})")

    # ---- item 3: STRICT split-half of the conservation score ---------------------------------
    rowsA = {k: H[k]["A"] for k in hkeys}; availA = {k: H[k]["availA"] for k in hkeys}
    rowsB = {k: H[k]["B"] for k in hkeys}; availB = {k: H[k]["availB"] for k in hkeys}
    sA, nA, CA, prs = submap_score(rowsA, availA, hkeys)
    sB, nB, CB, _ = submap_score(rowsB, availB, hkeys)
    rho_strict, ci_strict, n_strict = boot_spearman(sA, sB)
    print(f"STRICT split-half of the score: rho {rho_strict:+.3f} {ci_strict} on n={n_strict}")

    # a full-depth version of the same 4-source sub-map, for the prediction table
    rowsF, availF = {}, {}
    for k in hkeys:
        if k == "dld1_lowmoi":
            r, a = L.load_source_path("experiments/223-silent-zero-audit/work/tables/dld1_lowmoi", axis_h)
            rowsF[k] = r.astype(np.float64); availF[k] = a
        else:
            m = H[k]["avail"]
            R = H[k]["real"].copy()
            R[m] -= R[m].mean(axis=0, keepdims=True)
            R[~m] = 0.0
            rowsF[k] = R; availF[k] = m
    sFull, nFull, CFull, _ = submap_score(rowsF, availF, hkeys)

    # ---- item 2b: reliability-corrected score on the same four sources -----------------------
    rfull = {k: sb(H[k]["r_half"]) for k in hkeys}
    Crel = np.full_like(CFull, np.nan)
    dropped = 0
    for p, (i, j) in enumerate(prs):
        ki, kj = hkeys[i], hkeys[j]
        ri, rj = rfull[ki], rfull[kj]
        good = np.isfinite(CFull[p]) & (ri > R_FLOOR) & (rj > R_FLOOR)
        dropped += int((np.isfinite(CFull[p]) & ~good).sum())
        Crel[p, good] = CFull[p, good] / np.sqrt(ri[good] * rj[good])
    with np.errstate(invalid="ignore"):
        c_rel = np.nanmean(Crel, axis=0)
        n_rel = np.isfinite(Crel).sum(axis=0)
    print(f"reliability-corrected: n targets {int(np.isfinite(c_rel).sum())}, "
          f"median {np.nanmedian(c_rel):+.3f}, pairs dropped {dropped}")

    # ---- declared secondary: the FULL 26-source map's halves (shared sources -> upper bound) --
    #      half A = the four halved sources' A rows + every other source's full row
    zr = np.load(L.WORK / "rows276.npz", allow_pickle=True)
    Rall = np.asarray(zr["R"]); Aall = np.asarray(zr["A"])
    gi = {g: i for i, g in enumerate(axis_full)}
    colsh = np.asarray([gi[g] for g in axis_h if g in gi], dtype=np.int64)
    axis_hh = [g for g in axis_h if g in gi]
    rowsUA, availUA, rowsUB, availUB = {}, {}, {}, {}
    for i, k in enumerate(keys):
        if k in hkeys:
            gih = {g: q for q, g in enumerate(axis_h)}
            ch = np.asarray([gih[g] for g in axis_hh], dtype=np.int64)
            rowsUA[k] = H[k]["A"][:, ch]; availUA[k] = H[k]["availA"]
            rowsUB[k] = H[k]["B"][:, ch]; availUB[k] = H[k]["availB"]
        else:
            r = Rall[i][:, colsh].astype(np.float64)
            m = Aall[i]
            rc = r.copy(); rc[m] -= rc[m].mean(axis=0, keepdims=True); rc[~m] = 0.0
            rowsUA[k] = rc; availUA[k] = m
            rowsUB[k] = rc; availUB[k] = m
    uA, _, _, _ = submap_score(rowsUA, availUA, keys)
    uB, _, _, _ = submap_score(rowsUB, availUB, keys)
    rho_up, ci_up, n_up = boot_spearman(uA, uB)
    print(f"UPPER-BOUND split-half (shared full-depth sources included): rho {rho_up:+.3f} {ci_up} n={n_up}")

    np.savez_compressed(
        L.WORK / "rel276.npz", targets=np.asarray(L.PANEL, dtype=object),
        hkeys=np.asarray(hkeys, dtype=object), axis_h=np.asarray(axis_h, dtype=object),
        sA=sA, sB=sB, nA=nA, nB=nB, sFull=sFull, nFull=nFull,
        c_rel=c_rel, n_rel=n_rel, uA=uA, uB=uB,
        Crel=Crel.astype(np.float32), CFull=CFull.astype(np.float32),
        prs=np.asarray(prs), CA=CA.astype(np.float32), CB=CB.astype(np.float32),
        **{f"r_half|{k}": H[k]["r_half"] for k in hkeys},
        **{f"r_full|{k}": rfull[k] for k in hkeys})

    out = dict(
        axis_half_genes=len(axis_h),
        ka4_reliability_median=ka4, dld1_replicate_reliability_median=ka_dld,
        strict=dict(rho=rho_strict, ci=list(ci_strict), n=n_strict,
                    pairs=[[hkeys[i], hkeys[j]] for i, j in prs],
                    median_A=float(np.nanmedian(sA)), median_B=float(np.nanmedian(sB)),
                    n_pairs_median_A=float(np.median(nA[nA > 0])),
                    n_pairs_median_B=float(np.median(nB[nB > 0]))),
        upper_bound=dict(rho=rho_up, ci=list(ci_up), n=n_up),
        c_rel=dict(n=int(np.isfinite(c_rel).sum()), median=float(np.nanmedian(c_rel)),
                   mean=float(np.nanmean(c_rel)), p10=float(np.nanpercentile(c_rel, 10)),
                   p90=float(np.nanpercentile(c_rel, 90)), pairs_dropped=dropped,
                   share_le0=float(np.nanmean(c_rel[np.isfinite(c_rel)] <= 0))),
        c_full_submap=dict(median=float(np.nanmedian(sFull)),
                           n=int(np.isfinite(sFull).sum())),
        peak_rss_gib=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 2**30)
    (L.WORK / "rel276.json").write_text(json.dumps(out, indent=1))
    print(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()
