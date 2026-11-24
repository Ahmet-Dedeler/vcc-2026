"""290 step 1 — KA-1 (the row reader) and the response-derived feature block (PREREG §3, §5).

Loads each of the 28 banks once: computes that bank's per-gene response features, keeps its unit
rows, then rebuilds 279's `c_held` over all 364 cross-line pairs and asserts it at 0.0e+00.

    PYTHONDONTWRITEBYTECODE=1 uv run python build290.py

Writes work/feat290.csv, work/build290.json, work/build290.log, work/persrc290.npz.
"""
from __future__ import annotations

import json
import resource
import time

import numpy as np
import pandas as pd

import lib290 as L

T0 = time.time()
LOG: list[str] = []


def P(msg: str = "") -> None:
    line = f"[{time.time() - T0:7.1f}s] {msg}" if msg else ""
    print(line, flush=True)
    LOG.append(line)


def main() -> None:
    h = L.check_prereg()
    P(f"prereg {h} MATCH")

    per_src, genes_held, keys, axis = L.scan_sources()
    lines = [per_src[k]["line"] for k in keys]
    nG, nA = len(genes_held), len(axis)
    P(f"held genes {nG}; sources {len(keys)}; axis {nA} genes; {len(set(lines))} lines")
    assert (nG, len(keys), nA, len(set(lines))) == (17106, 28, 2739, 21), "source scan drifted"

    z = np.load(L.HELD279, allow_pickle=True)
    assert list(np.asarray(z["keys"]).astype(str)) == keys
    assert list(np.asarray(z["lines"]).astype(str)) == lines
    assert list(np.asarray(z["axis"]).astype(str)) == axis
    assert list(np.asarray(z["genes"]).astype(str)) == genes_held

    W, w_prolif, w_stress, psizes = L.programme_directions(axis)
    P("programme directions: " + "; ".join(
        f"{p.replace('HALLMARK_', '')} {psizes[p]['n_on_axis']}/{psizes[p]['n_set']}"
        for p in L.PROGRAMMES))
    P(f"  proliferation axis {int((w_prolif != 0).sum())} axis genes; "
      f"stress axis {int((w_stress != 0).sum())}")

    aidx = {g: i for i, g in enumerate(axis)}
    own_col = np.array([aidx.get(g, -1) for g in genes_held], dtype=np.int64)
    on_axis = own_col >= 0
    P(f"held genes on the 2,739-gene axis: {int(on_axis.sum())} ({on_axis.mean():.1%})")

    # ---- per-source features + unit rows ------------------------------------------------------
    QNAMES = (["mag", "breadth", "gen_cos", "pc5", "pc20"] + L.PE_COLS
              + ["proj_prolif", "proj_stress", "ontarget"])
    Q = {q: np.full((len(keys), nG), np.nan, np.float32) for q in QNAMES}
    idxs, Us = {}, {}
    src_info = {}

    for q, k in enumerate(keys):
        t1 = time.time()
        idx, U, nrm, blk, mu = L.load_unit(per_src[k]["rel"], axis, genes_held)
        idxs[k], Us[k] = idx, U
        n = idx.size

        # (a) magnitude — 279's scale-normalised row norm
        Q["mag"][q, idx] = (nrm / np.median(nrm)).astype(np.float32)
        # (b) breadth — count above the bank's 99th-percentile |residual|
        tau = float(np.percentile(np.abs(blk), L.BREADTH_Q))
        cnt = (np.abs(blk) > tau).sum(axis=1)
        Q["breadth"][q, idx] = np.log1p(cnt).astype(np.float32)
        # (c1) genericness — alignment with the bank's mean RAW row
        mun = float(np.linalg.norm(mu))
        Q["gen_cos"][q, idx] = (U @ (mu / max(mun, 1e-12)).astype(np.float32)).astype(np.float32)
        # (c3,c4) PC energy — exact eigendecomposition of the bank's Gram matrix
        if n >= L.PC_MIN_ROWS:
            G = blk.T @ blk
            ev, V = np.linalg.eigh(G)
            V = V[:, ::-1][:, :20]
            pr = (U.astype(np.float64) @ V) ** 2
            Q["pc5"][q, idx] = pr[:, :5].sum(axis=1).astype(np.float32)
            Q["pc20"][q, idx] = pr[:, :20].sum(axis=1).astype(np.float32)
            del G, V, pr
        # (d) programme engagement
        PJ = U @ W.T.astype(np.float32)
        for c, name in enumerate(L.PE_COLS):
            Q[name][q, idx] = PJ[:, c]
        Q["proj_prolif"][q, idx] = U @ w_prolif.astype(np.float32)
        Q["proj_stress"][q, idx] = U @ w_stress.astype(np.float32)
        # (f) on-target knockdown, where the gene is on the axis
        have_own = on_axis[idx]
        if have_own.any():
            rows = np.flatnonzero(have_own)
            Q["ontarget"][q, idx[rows]] = U[rows, own_col[idx[rows]]]

        src_info[k] = dict(line=per_src[k]["line"], dataset=per_src[k]["dataset"],
                           rel=per_src[k]["rel"], n_rows=int(n), tau=tau,
                           mean_breadth=float(np.expm1(np.log1p(cnt)).mean()),
                           median_breadth=float(np.median(cnt)),
                           mu_norm=mun, pcs=bool(n >= L.PC_MIN_ROWS))
        del blk
        P(f"  {k:22s} {per_src[k]['line']:16s} {n:6d} rows  tau {tau:.4f}  "
          f"median DE-proxy {np.median(cnt):6.0f}  ({time.time() - t1:.1f}s)")

    # ---- KA-1: 279's c_held, rebuilt over all cross-line pairs ---------------------------------
    cross = [(i, j) for i in range(len(keys)) for j in range(i + 1, len(keys))
             if lines[i] != lines[j]]
    assert np.array_equal(np.asarray(cross), np.asarray(z["pairs_cross"])), "pair list drifted"
    C = np.full((len(cross), nG), np.nan, np.float32)
    for p, (i, j) in enumerate(cross):
        g, c = L.pair_cos(idxs[keys[i]], Us[keys[i]], idxs[keys[j]], Us[keys[j]])
        C[p, g] = c
    with np.errstate(invalid="ignore"):
        c_held = np.nanmean(C, axis=0)
    ref = np.asarray(z["c_held"])
    dC = float(np.nanmax(np.abs(C - np.asarray(z["Ccross"]))))
    d = float(np.nanmax(np.abs(c_held - ref)))
    ka1 = dict(n_pairs=len(cross), max_abs_diff_c_held=d, max_abs_diff_Ccross=dC,
               n=int(np.isfinite(c_held).sum()), median=float(np.nanmedian(c_held)),
               pass_=bool(d == 0.0 and dC == 0.0 and len(cross) == 364))
    P(f"\nKA-1  c_held over {len(cross)} cross-line pairs: max|diff| vs held279.npz "
      f"{d:.1e} (Ccross {dC:.1e}), n {ka1['n']}, median {ka1['median']:+.6f} "
      f"-> {'PASS' if ka1['pass_'] else 'FAIL'}")
    if not ka1["pass_"]:
        raise SystemExit("KA-1 FAILED — abort (AGENTS.md:1375-1382)")
    del C, Us, idxs

    # ---- aggregate to the two halves ----------------------------------------------------------
    h1 = np.array([i for i, ln in enumerate(lines) if ln in L.HALF1_LINES])
    h2 = np.array([i for i, ln in enumerate(lines) if ln in L.HALF2_LINES])
    assert len(h1) + len(h2) == len(keys) and len(set(h1) & set(h2)) == 0
    P(f"\nhalf 1: {len(h1)} sources {[keys[i] for i in h1]}")
    P(f"half 2: {len(h2)} sources {[keys[i] for i in h2]}")

    out = pd.DataFrame({"gene": genes_held, "on_axis": on_axis.astype(float)})

    def agg(rows: np.ndarray, tag: str) -> None:
        with np.errstate(invalid="ignore"):
            mag = np.nanmean(Q["mag"][rows], axis=0)
            out[f"logmag_{tag}"] = np.log(np.where(mag > 0, mag, np.nan))
            out[f"breadth_{tag}"] = np.nanmean(Q["breadth"][rows], axis=0)
            gc = Q["gen_cos"][rows]
            out[f"gen_cos_{tag}"] = np.nanmean(gc, axis=0)
            out[f"gen_abscos_{tag}"] = np.nanmean(np.abs(gc), axis=0)
            out[f"gen_pc5_{tag}"] = np.nanmean(Q["pc5"][rows], axis=0)
            out[f"gen_pc20_{tag}"] = np.nanmean(Q["pc20"][rows], axis=0)
            for name in L.PE_COLS:
                out[f"{name}_{tag}"] = np.nanmean(Q[name][rows], axis=0)
            out[f"pe_prolif_abs_{tag}"] = np.nanmean(np.abs(Q["proj_prolif"][rows]), axis=0)
            out[f"pe_stress_abs_{tag}"] = np.nanmean(np.abs(Q["proj_stress"][rows]), axis=0)
            out[f"ontarget_{tag}"] = np.nanmean(Q["ontarget"][rows], axis=0)
            ns = np.isfinite(Q["mag"][rows]).sum(axis=0)
            out[f"n_src_{tag}"] = ns.astype(float)
            out[f"log_n_src_{tag}"] = np.log1p(ns)
        out[f"on_axis_{tag}"] = on_axis.astype(float)

    agg(h1, "h1")
    agg(h2, "h2")
    out.to_csv(L.WORK / "feat290.csv", index=False)
    np.savez_compressed(L.WORK / "persrc290.npz",
                        genes=np.asarray(genes_held, dtype=object),
                        keys=np.asarray(keys, dtype=object),
                        lines=np.asarray(lines, dtype=object),
                        half1=h1, half2=h2, on_axis=on_axis, **Q)

    cov = {}
    for tag in ("h1", "h2"):
        cov[tag] = {c: float(np.isfinite(out[f"{c}_{tag}"]).mean()) for c in L.FEATS290
                    if c != "on_axis"}
    P("\ncoverage over all 17,106 held genes (share finite):")
    for c in L.FEATS290:
        if c == "on_axis":
            continue
        P(f"  {c:34s} h1 {cov['h1'][c]:.3f}   h2 {cov['h2'][c]:.3f}")

    summary = dict(prereg_sha256=h, ka1=ka1, sources=src_info, programme_sizes=psizes,
                   half1=[keys[i] for i in h1], half2=[keys[i] for i in h2],
                   n_held=nG, n_axis=nA, n_on_axis=int(on_axis.sum()), coverage=cov,
                   seconds=time.time() - T0,
                   peak_rss_gib=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 2 ** 30)
    (L.WORK / "build290.json").write_text(json.dumps(summary, indent=1, default=float))
    (L.WORK / "build290.log").write_text("\n".join(LOG))
    P(f"peak RSS {summary['peak_rss_gib']:.2f} GiB; wrote work/feat290.csv")


if __name__ == "__main__":
    main()
