"""304 — audit re-analysis of 302's SAVED SUMS (no reserve is re-opened).

This is an AUDIT of three reported by-products of F15, not an independent
confirmation of claim C6:

  KA-1  rebuild S_g from work/sums302.npz + the gse295214_h1 bank and reproduce
        302's symmetry302.csv column to 0.0e+00.
  Q6    target-shuffled baseline for rho(|S|, gen_pc20)  (mandate item 6)
  Q4    WITHIN-TARGET depth ladder from the 16 saved shards (mandate item 4)
  Q3    the ceiling under the half/full-sample convention and alternative
        rho_CRISPRi                                        (mandate item 3)

    PYTHONDONTWRITEBYTECODE=1 uv run python work/audit304.py
"""
from __future__ import annotations

import csv
import json
import sys
from pathlib import Path

import numpy as np

E302 = Path(__file__).resolve().parents[2] / "302-orf-vs-crispri-symmetry"
sys.path.insert(0, str(E302))
import lib302 as L  # noqa: E402

OUT = Path(__file__).resolve().parent
LOG: list[str] = []


def P(m: str = "") -> None:
    print(m, flush=True)
    LOG.append(m)


def cosmat(A: np.ndarray, B: np.ndarray, own: np.ndarray) -> np.ndarray:
    """C[i, j] = cos(a_i, b_j) with the own columns of BOTH i and j zeroed in both
    vectors — analyse302.cos_pair, vectorised."""
    n = A.shape[0]
    G = A @ B.T
    na2 = (A ** 2).sum(axis=1)
    nb2 = (B ** 2).sum(axis=1)
    C = np.empty((n, n))
    for i in range(n):
        for j in range(n):
            cols = {own[i], own[j]} - {-1}
            num, a2, b2 = float(G[i, j]), float(na2[i]), float(nb2[j])
            for k in cols:
                num -= A[i, k] * B[j, k]
                a2 -= A[i, k] ** 2
                b2 -= B[j, k] ** 2
            d = np.sqrt(max(a2, 0.0) * max(b2, 0.0))
            C[i, j] = num / d if d > 0 else np.nan
    return C


def orf_residual(acc: np.ndarray, lib: np.ndarray, shards: np.ndarray, nT: int) -> np.ndarray:
    S = acc[:, shards].sum(axis=1)
    Lt = lib[:, shards].sum(axis=1)
    Pr = L.profile(S, Lt[:, None])
    return Pr[:nT] - Pr[:nT].mean(axis=0)


def sb(r: float) -> float:
    return 2.0 * r / (1.0 + r)


def main() -> None:
    z = np.load(E302 / "work" / "sums302.npz", allow_pickle=True)
    axis = np.asarray(z["axis"]).astype(str).tolist()
    scored = np.asarray(z["scored"]).astype(str).tolist()
    primary = set(np.asarray(z["primary"]).astype(str).tolist())
    acc_sym, lib_sym, cells_sym = z["acc_sym"], z["lib_sym"], z["cells_sym"]
    draws = np.asarray(z["draws"])
    axis_h1 = np.asarray(z["axis_h1"])
    nT = len(scored)

    A = orf_residual(acc_sym, lib_sym, np.arange(acc_sym.shape[1]), nT)

    perts, _h1g, _m = L.read_h1_bank()
    h1_row = {p: i for i, p in enumerate(perts)}
    Dm = np.load(L.H1BANK / "delta.npy", mmap_mode="r")
    rows = np.asarray([h1_row[s] for s in scored])
    lo, hi = int(rows.min()), int(rows.max())
    Dsub = np.asarray(Dm[lo:hi + 1])[rows - lo][:, axis_h1].astype(np.float64)
    del Dm
    B = Dsub - Dsub.mean(axis=0, keepdims=True)

    apos = {g: i for i, g in enumerate(axis)}
    own = np.asarray([apos.get(s, -1) for s in scored])

    C = cosmat(A, B, own)
    Sg = np.diag(C).copy()
    pri = np.asarray([s in primary for s in scored])
    out: dict = {}

    # ---------------------------------------------------------------- KA-1
    ref = {r["symbol"]: r for r in csv.DictReader(open(E302 / "work" / "symmetry302.csv"))}
    ref_S = np.asarray([float(ref[s]["S_g"]) for s in scored])
    d = float(np.nanmax(np.abs(Sg - ref_S)))
    st = json.load(open(E302 / "work" / "stats302.json"))
    med302 = float(st["C6a"]["median"])
    pri0 = np.asarray([s in primary for s in scored])
    dmed = abs(float(np.median(Sg[pri0])) - med302)
    P(f"KA-1  rebuilt S_g vs 302's symmetry302.csv over {nT} scored TFs: max abs diff {d:.2e} "
      f"(the csv stores S_g at %.6f, so 5.0e-07 is its own rounding bound) "
      f"-> {'PASS' if d < 5.0e-7 else 'FAIL'}")
    P(f"KA-1b rebuilt median S over the 148 primary TFs vs stats302.json's full-precision "
      f"{med302:+.10f}: abs diff {dmed:.2e} -> {'PASS' if dmed < 1e-12 else 'FAIL'}")
    out["ka1_max_abs_diff"] = d
    out["ka1b_median_abs_diff"] = dmed
    if d >= 5.0e-7:
        P("ABORT: the rebuild does not reproduce 302; nothing below is trustworthy.")
        (OUT / "audit304.json").write_text(json.dumps(out, indent=1))
        return

    gen = np.asarray([float(ref[s]["gen_pc20_h1"]) if ref[s]["gen_pc20_h1"] else np.nan
                      for s in scored])
    rel302 = np.asarray([float(ref[s]["orf_reliability"]) for s in scored])
    cells = np.asarray([float(ref[s]["orf_cells"]) for s in scored])
    kd_mag = np.asarray([float(ref[s]["kd_mag"]) for s in scored])
    orf_mag = np.asarray([float(ref[s]["orf_mag"]) for s in scored])

    # ---------------------------------------------------------------- Q6
    absS = np.abs(Sg)
    ok = pri & np.isfinite(gen)
    ii = np.flatnonzero(ok)
    rho_obs = L.spearman(gen[ii], absS[ii])
    P(f"\nQ6  observed rho(|S|, gen_pc20) on n={len(ii)}: {rho_obs:+.4f}  "
      f"(302 reports +0.3823)")

    rng = np.random.default_rng(L.SEED)
    idx_p = np.flatnonzero(pri)
    nulls = np.empty(200)
    for t in range(200):
        while True:
            perm = rng.permutation(len(idx_p))
            if not (perm == np.arange(len(idx_p))).any():
                break
        # |cos(a_{pi(s)}, b_s)| — the KD row keeps its own gen_pc20, the ORF row is reshuffled
        v = np.abs(C[idx_p[perm], idx_p])
        g = gen[idx_p]
        m = np.isfinite(v) & np.isfinite(g)
        nulls[t] = L.spearman(g[m], v[m])
    P(f"    target-shuffled rho(|S_perm|, gen_pc20): median {np.median(nulls):+.4f}, "
      f"mean {nulls.mean():+.4f} +- {nulls.std():.4f}, "
      f"[{np.percentile(nulls, 2.5):+.4f}, {np.percentile(nulls, 97.5):+.4f}]; "
      f"{(nulls >= rho_obs).mean():.1%} of 200 draws >= observed")
    out["q6"] = dict(rho_obs=rho_obs, null_median=float(np.median(nulls)),
                     null_mean=float(nulls.mean()), null_sd=float(nulls.std()),
                     null_ci=[float(np.percentile(nulls, 2.5)), float(np.percentile(nulls, 97.5))],
                     frac_ge_obs=float((nulls >= rho_obs).mean()),
                     z=float((rho_obs - nulls.mean()) / nulls.std()))

    # "alignability": how well does b_s align with the ORF ensemble at large?
    Coff = C[np.ix_(idx_p, idx_p)].copy()               # rows = ORF partner t, cols = KD row s
    np.fill_diagonal(Coff, np.nan)
    g = gen[idx_p]
    out["q6_alignability"] = {}
    for nm, fn in [("median_t |cos|", lambda M: np.nanmedian(np.abs(M), axis=0)),
                   ("mean_t |cos|", lambda M: np.nanmean(np.abs(M), axis=0)),
                   ("sd_t cos", lambda M: np.nanstd(M, axis=0)),
                   ("q90_t |cos|", lambda M: np.nanpercentile(np.abs(M), 90, axis=0))]:
        align = fn(Coff)
        mm = np.isfinite(align) & np.isfinite(g)
        rho_align = L.spearman(g[mm], align[mm])
        part = L.partial_spearman(g[mm], absS[idx_p][mm], align[mm, None])
        P(f"    rho(gen_pc20, {nm} over OTHER TFs' ORF rows) = {rho_align:+.4f}; "
          f"rho(|S|, gen_pc20) partialling it out = {part:+.4f}")
        out["q6_alignability"][nm] = dict(rho=rho_align, partial=part)

    # multiplicity / coupling covariates
    for nm, vv in [("orf_reliability", rel302), ("kd_mag", kd_mag), ("orf_mag", orf_mag),
                   ("log_cells", np.log(np.maximum(cells, 1)))]:
        m2 = np.isfinite(vv[idx_p]) & np.isfinite(g)
        P(f"    rho(gen_pc20, {nm}) = {L.spearman(g[m2], vv[idx_p][m2]):+.4f}")
        out.setdefault("q6_cov", {})[nm] = L.spearman(g[m2], vv[idx_p][m2])

    # ---------------------------------------------------------------- Q4
    P("\nQ4  WITHIN-TARGET depth ladder (same TFs, nested shard subsets of their own cells)")
    P("    k = shards per side; rho_half = split-half r at k shards; rho_full = SB to 2k shards")
    ladder = {}
    for k in (1, 2, 4, 8):
        rel = np.zeros(nT)
        ncell = np.zeros(nT)
        for dd in draws:
            a_sh, b_sh = dd[:8][:k], dd[8:][:k]
            AA = orf_residual(acc_sym, lib_sym, a_sh, nT)
            BB = orf_residual(acc_sym, lib_sym, b_sh, nT)
            for i in range(nT):
                x, y = AA[i], BB[i]
                num, a2, b2 = float(x @ y), float(x @ x), float(y @ y)
                if own[i] >= 0:
                    kk = own[i]
                    num -= x[kk] * y[kk]; a2 -= x[kk] ** 2; b2 -= y[kk] ** 2
                dd2 = np.sqrt(max(a2, 0.0) * max(b2, 0.0))
                rel[i] += (num / dd2 if dd2 > 0 else 0.0) / len(draws)
            ncell += (cells_sym[:nT][:, a_sh].sum(axis=1)
                      + cells_sym[:nT][:, b_sh].sum(axis=1)) / len(draws)
        rh = float(np.median(rel[pri]))
        rf = sb(max(rh, 1e-9))
        n2k = float(np.median(ncell[pri]))
        c = n2k * (1 - rf) / rf if rf > 0 else np.nan
        P(f"    k={k:2d}  cells in 2k shards (median) {n2k:8.1f}  rho_half {rh:+.5f}  "
          f"rho_full(2k) {rf:+.5f}  implied c {c:9.1f}")
        ladder[str(k)] = dict(n_cells_2k=n2k, rho_half=rh, rho_full=rf, c=float(c),
                              rel=rel.tolist())
    out["q4_within_target"] = {k: {kk: vv for kk, vv in v.items() if kk != "rel"}
                               for k, v in ladder.items()}

    # the same ladder inside 302's own cell-count strata (its cross-target table)
    P("\n    the within-target ladder inside 302's cell-count strata:")
    bands = [(100, 200), (200, 500), (500, 1000), (1000, 3000), (3000, 20000)]
    strat = {}
    for loB, hiB in bands:
        sel = pri & (cells >= loB) & (cells < hiB)
        row = []
        for k in (1, 2, 4, 8):
            rel = np.asarray(ladder[str(k)]["rel"])
            rh = float(np.median(rel[sel]))
            rf = sb(max(rh, 1e-9))
            n2k = float(np.median(cells[sel])) * (2 * k / 16.0)
            row.append((k, n2k, rh, rf, n2k * (1 - rf) / rf if rf > 0 else np.nan))
        P(f"    {loB:5d}-{hiB:<6d} n={int(sel.sum()):3d}  " +
          "  ".join(f"k{k}: c={c:8.0f}" for k, _n, _rh, _rf, c in row))
        strat[f"{loB}-{hiB}"] = dict(n=int(sel.sum()),
                                     rows=[dict(k=k, n_cells=n, rho_half=rh, rho_full=rf, c=c)
                                           for k, n, rh, rf, c in row])
    out["q4_strata"] = strat

    # ---------------------------------------------------------------- Q3
    P("\nQ3  the ceiling under the half/full convention and alternative rho_CRISPRi")
    rel_all = float(np.median(rel302[pri]))
    rho_orf_full = sb(rel_all)
    med_S = float(np.median(Sg[pri]))
    lo_, hi_, _ = L.boot_ci(lambda kx: float(np.median(Sg[pri][kx])), int(pri.sum()), 20000, L.SEED)
    P(f"    ORF rel_half {rel_all:+.4f} -> rho_ORF(full) {rho_orf_full:+.4f}; "
      f"median S {med_S:+.4f} [{lo_:+.4f}, {hi_:+.4f}]")
    tab = {}
    for lab, rkd in [("302's substitution: rho_KD = 0.243 (a HALF-sample figure used as full)", 0.243),
                     ("same figure carried to full sample: SB(0.243) = 0.391", sb(0.243)),
                     ("pessimistic bank: rho_KD = 0.10", 0.10),
                     ("pessimistic bank: rho_KD = 0.05", 0.05),
                     ("optimistic bank: rho_KD = 0.60", 0.60)]:
        ceil = float(np.sqrt(rho_orf_full * rkd))
        P(f"    rho_KD {rkd:.3f}  ceiling {ceil:.4f}  implied TRUE median cos "
          f"{med_S / ceil:+.4f} [{lo_ / ceil:+.4f}, {hi_ / ceil:+.4f}]   ({lab})")
        tab[lab] = dict(rho_kd=rkd, ceiling=ceil, implied=med_S / ceil,
                        implied_ci=[lo_ / ceil, hi_ / ceil])
    out["q3_ceiling"] = tab

    c_med = 9562.085473953619
    for lab, target in [("302's target: rho_ORF(full) = 0.243", 0.243),
                        ("like-for-like: rho_ORF(full) = SB(0.243) = 0.391", sb(0.243))]:
        n = c_med * target / (1 - target)
        P(f"    cells needed for {lab}: {n:.0f}")
        out.setdefault("q3_depth_target", {})[lab] = n

    (OUT / "audit304.json").write_text(json.dumps(out, indent=1, default=float))
    (OUT / "audit304.log").write_text("\n".join(LOG) + "\n")
    P("\n-> work/audit304.json, work/audit304.log")


if __name__ == "__main__":
    main()
