"""279 — split-half reliability for the HELD genes, the reliability-corrected conservation
`C_rel` (DEFINITION-276 §3, unchanged) and the noise-vs-state-dependence descriptor.

Four sources have independent halves for non-panel genes:
  * `k562_gwps`, `orion_hct116`, `orion_hek293t` — 267's batch-parity accumulators
    (`experiments/267-.../work/half/<line>/{sums.npz, acc_batch{A,B}.npy}`, read-only); the half
    rows are rebuilt here with `truth267.write_half_bank`'s arithmetic for the held genes only
    (267 deleted the half `delta.npy` files after use).  KA-3 checks the rebuild against 271's
    panel reliabilities.
  * `dld1_lowmoi` — its two low-MOI replicate banks, already on disk.

NOTE ON THE RESERVE: programme rule 2 reserves the **`orion_hek293t` PLOLO fold's split-half B
TRUTHS on the 300 panel targets**.  Nothing here touches a panel target — the held set excludes
all 300 by construction — and every statistic is also reported without `orion_hek293t`.

    PYTHONDONTWRITEBYTECODE=1 uv run python rel279.py

Writes work/rel279.npz, work/rel279.json.  Read-only outside work/.
"""
from __future__ import annotations

import json
import resource
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path("~/Code/virtual-cell-challenge-2026")
HERE = ROOT / "experiments/279-znf-state-dependence"
WORK = HERE / "work"
HALF267 = ROOT / "experiments/267-replicate-validated-oracle/work/half"
sys.path.insert(0, str(ROOT / "experiments/276-conservation-map"))
sys.path.insert(0, str(HERE))
import lib276 as L  # noqa: E402
import build279 as B  # noqa: E402

TS = 5e4
MIN_HALF = 5            # 267's rule
R_FLOOR = 0.02          # DEFINITION-276 §3
ROW_BLK = 1024

HALFED = {
    "k562_gwps":     ("K562",    "data/donors/k562_gwps"),
    "orion_hct116":  ("HCT116",  "data/donors/orion_hct116"),
    "orion_hek293t": ("HEK293T", "data/donors/orion_hek293t"),
    "dld1_lowmoi":   ("DLD-1",   "experiments/223-silent-zero-audit/work/tables/dld1_lowmoi"),
}
DLD_REPS = ("experiments/223-silent-zero-audit/work/tables/dld1_low1",
            "experiments/223-silent-zero-audit/work/tables/dld1_low2")


def sb(r: np.ndarray) -> np.ndarray:
    with np.errstate(invalid="ignore", divide="ignore"):
        return 2.0 * r / (1.0 + r)


def half_rows(line: str, side: str, axis: list[str], genes: list[str]):
    """The half's residual rows for `genes`, on `axis` — write_half_bank's arithmetic, in blocks."""
    d = HALF267 / line
    z = np.load(d / "sums.npz", allow_pickle=True)
    keep = np.asarray(z["keep"], dtype=np.int64)
    keep_names = [str(x) for x in z["keep_names"]]
    pert_names = [str(x) for x in z["pert_names"]]
    tot = np.asarray(z["tot"], dtype=np.int64)
    ncell = np.asarray(z["ncell"], dtype=np.int64)
    ctrl_sum = np.asarray(z["ctrl_sum"], dtype=np.float64)
    col = 0 if side == "A" else 1
    cs = ctrl_sum[col]
    b_all = np.log1p(TS * cs / cs.sum()).astype(np.float32)
    b_keep = b_all[keep]
    on_all_genes = line == "k562_gwps"

    ki = {g: i for i, g in enumerate(keep_names)}
    cols_axis = np.asarray([ki[g] for g in axis], dtype=np.int64)
    bank_perts = set((ROOT / f"data/donors/{line}/perts.txt").read_text().split("\n"))
    pi = {nm: i for i, nm in enumerate(pert_names)}
    tidx = {g: i for i, g in enumerate(genes)}
    have = [(g, pi[g]) for g in genes
            if g in pi and g in bank_perts and ncell[pi[g], col] >= MIN_HALF and tot[pi[g], col] > 0]
    if not have:
        return np.zeros(0, np.int64), np.zeros((0, len(axis)), np.float32)
    order = np.argsort([r for _, r in have])
    rows = np.asarray([have[i][1] for i in order], dtype=np.int64)
    out = np.empty((rows.size, len(axis)), np.float32)
    bk = b_keep[cols_axis]
    cols_raw = keep[cols_axis] if on_all_genes else cols_axis
    for lo in range(0, rows.size, ROW_BLK):
        r = rows[lo:lo + ROW_BLK]
        acc = np.load(d / f"acc_batch{side}.npy", mmap_mode="r")   # re-mapped per block so the
        blk = np.asarray(acc[r][:, cols_raw], dtype=np.float64)    # file's pages never accumulate
        acc._mmap.close(); del acc                                 # in this process's RSS
        prof = np.log1p(TS * blk / tot[r, col][:, None]).astype(np.float32)
        out[lo:lo + r.size] = prof - bk
    idx = np.asarray([tidx[have[i][0]] for i in order], dtype=np.int64)
    o = np.argsort(idx)
    idx, out = idx[o], out[o]
    out -= out.mean(axis=0, keepdims=True)          # residual centring, over the half's own rows
    return idx, out


def cos_pair(idx_a, A, idx_b, Bm, n_out):
    common = np.intersect1d(idx_a, idx_b, assume_unique=True)
    out = np.full(n_out, np.nan)
    if common.size == 0:
        return out
    pa = np.searchsorted(idx_a, common)
    pb = np.searchsorted(idx_b, common)
    for lo in range(0, common.size, 4096):
        hi = min(lo + 4096, common.size)
        x = A[pa[lo:hi]].astype(np.float64); y = Bm[pb[lo:hi]].astype(np.float64)
        nx = np.linalg.norm(x, axis=1); ny = np.linalg.norm(y, axis=1)
        d = nx * ny
        c = np.full(hi - lo, np.nan)
        ok = d > 0
        c[ok] = (x[ok] * y[ok]).sum(axis=1) / d[ok]
        out[common[lo:hi]] = c
    return out


def main() -> None:
    t0 = time.time()
    z = np.load(WORK / "held279.npz", allow_pickle=True)
    genes = [str(g) for g in z["genes"]]
    nG = len(genes)
    keys = list(HALFED)
    axis = L.common_axis_paths([HALFED[k][1] for k in keys])
    print(f"half-map axis: {len(axis)} genes (276's rel276 used 7543 on the panel)", flush=True)

    # ---- full-depth rows of the four sources on the half-map axis, and their pair cosines ------
    #      (computed and freed BEFORE the half rows are built, so the two big stages never
    #       overlap in memory)
    idxF, UF = {}, {}
    for k in keys:
        i, u, _ = B.load_unit(HALFED[k][1], axis, genes)
        idxF[k], UF[k] = i, u
        print(f"  full {k:16s} {i.size:6d} held rows ({time.time()-t0:.0f}s)", flush=True)
    lines0 = [HALFED[k][0] for k in keys]
    prs0 = [(i, j) for i in range(len(keys)) for j in range(i + 1, len(keys))
            if lines0[i] != lines0[j]]
    Cfull = np.full((len(prs0), nG), np.nan)
    for p, (i, j) in enumerate(prs0):
        g, cc = B.pair_cos(idxF[keys[i]], UF[keys[i]], idxF[keys[j]], UF[keys[j]])
        Cfull[p, g] = cc
    del UF, idxF

    # ---- per-source split-half reliability on the held genes ----------------------------------
    r_half = {}
    for k in keys:
        t1 = time.time()
        if k == "dld1_lowmoi":
            ia, A = B.residual_rows(DLD_REPS[0], axis, genes)
            ib, Bm = B.residual_rows(DLD_REPS[1], axis, genes)
        else:
            ia, A = half_rows(k, "A", axis, genes)
            ib, Bm = half_rows(k, "B", axis, genes)
        r_half[k] = cos_pair(ia, A, ib, Bm, nG)
        n = int(np.isfinite(r_half[k]).sum())
        print(f"  r_half {k:16s} n {n:6d}  median {np.nanmedian(r_half[k]):+.5f} "
              f"({time.time()-t1:.0f}s)", flush=True)
        del A, Bm

    # ---- KA-3 / KA-4: the same rebuild on the PANEL reproduces 271's medians ------------------
    ka = {}
    for k in ("k562_gwps", "orion_hct116", "orion_hek293t"):
        ia, A = half_rows(k, "A", axis, L.PANEL)
        ib, Bm = half_rows(k, "B", axis, L.PANEL)
        rp = cos_pair(ia, A, ib, Bm, len(L.PANEL))
        ka[k] = float(np.nanmedian(rp))
        del A, Bm
    ia, A = B.residual_rows(DLD_REPS[0], axis, L.PANEL)
    ib, Bm = B.residual_rows(DLD_REPS[1], axis, L.PANEL)
    ka["dld1_lowmoi"] = float(np.nanmedian(cos_pair(ia, A, ib, Bm, len(L.PANEL))))
    del A, Bm
    want = dict(k562_gwps=0.04706, orion_hct116=0.02756, orion_hek293t=0.02825,
                dld1_lowmoi=0.086321)
    print("KA-3/4 panel reliability medians (rebuilt on the 7,543-gene half axis):")
    for k, v in ka.items():
        print(f"    {k:16s} {v:+.5f}   (271/276 on their axes: {want[k]:+.5f})")

    # KA-4b: the rebuilt panel reliability against 271's OWN per-target array.  271's r_half_raw
    # lives on the 18,533-official axis, centred on the bank's every-perturbation mean with the
    # panel's own columns zeroed, so the medians are not the same number; the per-target ORDERING
    # must agree if the rebuild is measuring the same thing.
    from scipy.stats import spearmanr
    ka4b = {}
    for k in ("k562_gwps", "orion_hct116", "orion_hek293t"):
        z271 = np.load(ROOT / f"experiments/271-donor-reliability-estimation/work/half/rows-{k}.npz",
                       allow_pickle=True)
        assert [str(t) for t in z271["targets"]] == L.PANEL
        r271 = np.asarray(z271["r_half_raw"], dtype=np.float64)
        ia, A = half_rows(k, "A", axis, L.PANEL)
        ib, Bm = half_rows(k, "B", axis, L.PANEL)
        mine = cos_pair(ia, A, ib, Bm, len(L.PANEL))
        del A, Bm
        ok = np.isfinite(r271) & np.isfinite(mine)
        ka4b[k] = dict(spearman=float(spearmanr(r271[ok], mine[ok]).statistic),
                       pearson=float(np.corrcoef(r271[ok], mine[ok])[0, 1]),
                       n=int(ok.sum()), median_271=float(np.median(r271[ok])),
                       median_mine=float(np.median(mine[ok])))
        print(f"KA-4b {k:16s} rho(271's r_half_raw, rebuilt) {ka4b[k]['spearman']:+.3f} "
              f"(pearson {ka4b[k]['pearson']:+.3f}) on n={ka4b[k]['n']}")

    # ---- C_rel on the held genes (DEFINITION-276 §3) ------------------------------------------
    prs = prs0
    rfull = {k: sb(r_half[k]) for k in keys}
    Crel = np.full((len(prs), nG), np.nan)
    dropped = 0
    for p, (i, j) in enumerate(prs):
        ki, kj = keys[i], keys[j]
        ri, rj = rfull[ki], rfull[kj]
        good = np.isfinite(Cfull[p]) & (ri > R_FLOOR) & (rj > R_FLOOR)
        dropped += int((np.isfinite(Cfull[p]) & ~good).sum())
        Crel[p, good] = Cfull[p, good] / np.sqrt(ri[good] * rj[good])
    keep3 = [p for p, (i, j) in enumerate(prs) if "orion_hek293t" not in (keys[i], keys[j])]
    with np.errstate(invalid="ignore"):
        c_full4 = np.nanmean(Cfull, axis=0)
        c_rel = np.nanmean(Crel, axis=0)
        c_rel3 = np.nanmean(Crel[keep3], axis=0)
        n_rel = np.isfinite(Crel).sum(axis=0)
        rbar = np.nanmedian(np.vstack([r_half[k] for k in keys]), axis=0)
    print(f"C_rel: n {int(np.isfinite(c_rel).sum())}, median {np.nanmedian(c_rel):+.4f}, "
          f"pairs dropped at the r<=0.02 floor {dropped}", flush=True)

    np.savez_compressed(
        WORK / "rel279.npz", genes=np.asarray(genes, dtype=object),
        keys=np.asarray(keys, dtype=object), axis=np.asarray(axis, dtype=object),
        prs=np.asarray(prs), Cfull=Cfull.astype(np.float32), Crel=Crel.astype(np.float32),
        c_full4=c_full4, c_rel=c_rel, c_rel3=c_rel3, n_rel=n_rel, rbar=rbar,
        **{f"r_half|{k}": r_half[k] for k in keys})
    out = dict(axis_genes=len(axis), pairs=[[keys[i], keys[j]] for i, j in prs],
               ka_panel_reliability=ka, ka_reference=want, ka4b=ka4b,
               r_half_median={k: float(np.nanmedian(r_half[k])) for k in keys},
               r_half_n={k: int(np.isfinite(r_half[k]).sum()) for k in keys},
               c_rel=dict(n=int(np.isfinite(c_rel).sum()), median=float(np.nanmedian(c_rel)),
                          pairs_dropped=dropped),
               c_full4=dict(n=int(np.isfinite(c_full4).sum()), median=float(np.nanmedian(c_full4))),
               seconds=time.time() - t0,
               peak_rss_gib=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 2**30)
    (WORK / "rel279.json").write_text(json.dumps(out, indent=1))
    print(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()
