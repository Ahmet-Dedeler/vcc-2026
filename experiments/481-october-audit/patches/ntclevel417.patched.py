"""417: per-context bulk re-level for ANY bundle (A/B/C or D/E/F) from the artifact and the bundle's own control cells.

The mse cliff of context c sits at bulk scale s*_c = sqrt(claim_c / B_c):
  B_c      our artifact's across-perturbation spread (cell_eval2 `_across_pert_budget`, streamed from the h5ad),
  claim_c  = f_c * Cp46_c * sum(w_c)  (saturated: ballast's jackknife boost keeps our C_hat_p >= C_p),
  Cp46_c   the mean over the context's 46 NTC constructs of cell_eval2's `jackknife_correction` on each construct's
           400 cells (180's `per_construct_jackknife`), i.e. the reference's per-perturbation sampling noise proxy,
  f_c      1 for a context with no board (A/B/C's board fit: 1.00 / 1.00 / 0.9825).
187's fitted A/B/C claims equal f * Cp46 * sum(w) to 1.0000 / 0.9927 / 1.0000 (DEF-RUNBOOK, 2026-09-23).

The recommended scale is the EV-optimal point under cliff uncertainty sd (in ln scale) with the asymmetric local slopes
measured on A/B/C (raw mse per unit scale: -0.31 below the cliff, +1.78 above): the loss-weighted quantile
q = 0.31 / (0.31 + 1.78) = 0.148, s = s* exp(sd * Phi^-1(q)).  sd 0.01 for an unfitted context (f unknown), 0.002 on A/B/C.

KNOWN ANSWER (--ka): on an A/B/C artifact with f = (1, 1, 0.9825), Cp46 must reproduce 180's cp_all (14.4398 / 14.8057 /
16.4397) and the cliffs must equal ctxlevel417.py's (187's claims) to <= 0.5 %.

usage: uv run python experiments/417-bulk-sharpen/ntclevel417.py <artifact.h5ad> [--bundle data/challenge]
       [--contexts A,B,C] [--f A=1,B=1,C=0.9825] [--sd 0.01] [--ka]
"""
from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

import h5py
import numpy as np
from scipy.stats import norm

ROOT = Path(__file__).resolve().parents[2]
TS = 50000.0
SLOPE_BELOW, SLOPE_ABOVE = 0.31, 1.78
CP46_ABC = {"A": 14.439751080335085, "B": 14.805734135005, "C": 16.439650143279888}   # 180 knee180.json cp_all


def str_col(h, path):
    obj = h[path]
    if isinstance(obj, h5py.Group):
        if "categories" in obj:
            return np.asarray(obj["categories"][:]).astype(str)[np.asarray(obj["codes"][:])]
        return np.asarray(obj["values"][:]).astype(str)
    return np.asarray(obj[:]).astype(str)


def cp46(bundle: Path, ctx: str) -> tuple[float, float, int]:
    """Mean/sd of the delete-1 jackknife over the context's NTC constructs (180's per_construct_jackknife), each construct
    put on the reference's footing: the jackknife goes as 1/n (481: x2.003 at 200 cells, x0.5002 at 800), so a construct
    of n cells is multiplied by n / cells_per_pert. At n == cells_per_pert (A/B/C: 46 x 400) the factor is exactly 1.0 and
    every value is unchanged bit for bit (481 KA). A construct under 100 cells is refused."""
    import anndata as ad
    from scipy.sparse import csr_matrix
    from cell_eval2.moments import jackknife_correction
    cpp = int(json.loads((bundle / "manifest.json").read_text())["cells_per_pert"])
    handle = ad.read_h5ad(bundle / f"context_{ctx}.h5ad")
    matrix = handle.X.tocsr()
    ids = handle.obs["ntc_id"].astype(str).to_numpy()
    vals, sizes = [], []
    for name in dict.fromkeys(ids):
        rows = np.flatnonzero(ids == name)
        if rows.size < 100:
            raise SystemExit(f"{ctx}: NTC construct {name!r} has {rows.size} cells (< 100); its jackknife cannot proxy C_p")
        sizes.append(int(rows.size))
        vals.append(jackknife_correction(csr_matrix(matrix[rows]), np.zeros(rows.size, np.intp), 1, TS, chunk=64)[0]
                    * (rows.size / float(cpp)))
    if set(sizes) != {cpp}:
        print(f"  {ctx}: NTC construct sizes {min(sizes)}..{max(sizes)} (not all {cpp}): jackknife rescaled by n/{cpp}", flush=True)
    v = np.asarray(vals, dtype=np.float64)
    return float(v.mean()), float(v.std()), int(v.size)


def decompose(path: Path, bundle: Path, contexts):
    from cell_eval2.metrics.delta import _across_pert_budget, _row_weights
    with h5py.File(path, "r") as h:
        genes = str_col(h, "var/_index")
        lab = list(zip(str_col(h, "obs/context").tolist(), str_col(h, "obs/target_gene").tolist()))
        indptr = np.asarray(h["X"]["indptr"][:])
        blocks, start = [], 0
        for i in range(1, len(lab) + 1):
            if i == len(lab) or lab[i] != lab[start]:
                blocks.append((lab[start], start, i))
                start = i
        totals = np.zeros((len(blocks), genes.size), dtype=np.float64)
        for k, (_, lo_row, hi_row) in enumerate(blocks):
            lo, hi = int(indptr[lo_row]), int(indptr[hi_row])
            totals[k] = np.bincount(h["X"]["indices"][lo:hi].astype(np.int64),
                                    weights=h["X"]["data"][lo:hi].astype(np.float64), minlength=genes.size)
    keys = [k for k, _, _ in blocks]
    gpos = {g: i for i, g in enumerate(genes)}
    out = {}
    for ctx in contexts:
        with h5py.File(bundle / f"context_{ctx}.h5ad", "r") as hc:
            assert np.array_equal(str_col(hc, "var/_index"), genes), f"gene axis mismatch on {ctx}"
            pooled = np.bincount(np.asarray(hc["X"]["indices"][:], dtype=np.int64),
                                 weights=np.asarray(hc["X"]["data"][:], dtype=np.float64), minlength=genes.size)
        b_ctrl = np.log1p(TS * pooled / pooled.sum())
        rows = [i for i, k in enumerate(keys) if k[0] == ctx]
        tot = totals[rows]
        mu = np.log1p(TS * tot / tot.sum(1, keepdims=True))
        dh = mu - b_ctrl[None, :]
        cols = {j: gpos[keys[i][1]] for j, i in enumerate(rows) if keys[i][1] in gpos}
        w = _row_weights(len(rows), genes.size, cols)
        a = np.einsum("ij,ij->i", dh, dh)
        for j, col in cols.items():
            a[j] -= dh[j, col] ** 2
        out[ctx] = dict(A=float((a * w).sum()), B=float(_across_pert_budget(mu, cols, w)), wsum=float(w.sum()),
                        n_rows=len(rows), n_resolved=len(cols))
    return out


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("artifact")
    ap.add_argument("--bundle", default=str(ROOT / "data/challenge"))
    ap.add_argument("--contexts", default="A,B,C")
    ap.add_argument("--f", default="")
    ap.add_argument("--sd", type=float, default=0.01)
    ap.add_argument("--ka", action="store_true")
    a = ap.parse_args()
    bundle, ctxs = Path(a.bundle), [c.strip() for c in a.contexts.split(",")]
    f = {c: 1.0 for c in ctxs}
    if a.f:
        f.update({k: float(v) for k, v in (p.split("=") for p in a.f.split(","))})
    t0 = time.time()
    dec = decompose(Path(a.artifact), bundle, ctxs)
    q = SLOPE_BELOW / (SLOPE_BELOW + SLOPE_ABOVE)
    res, scales = {}, {}
    for c in ctxs:
        m, sd_c, n = cp46(bundle, c)
        if a.ka and c in CP46_ABC:
            assert abs(m / CP46_ABC[c] - 1) < 1e-9, f"KA Cp46 {c}: {m} vs {CP46_ABC[c]}"
        claim = f[c] * m * dec[c]["wsum"]
        cliff = float(np.sqrt(claim / dec[c]["B"]))
        s = cliff * float(np.exp(a.sd * norm.ppf(q)))
        scales[c] = s
        res[c] = dict(Cp46=m, Cp46_sd=sd_c, n_constructs=n, f=f[c], claim=claim, **dec[c], rho=dec[c]["B"] / claim, cliff=cliff, scale=s)
        print(f"  {c}: Cp46 {m:.4f} (sd {sd_c:.3f}, {n} constructs)  B {dec[c]['B']:.6f}  claim {claim:.6f}  rho {dec[c]['B'] / claim:.4f}"
              f"  cliff {cliff:.4f}  EV-scale (sd {a.sd}) {s:.4f}", flush=True)
    if a.ka:
        sys.path.insert(0, str(ROOT / "experiments/187-joint-amplitude")); sys.path.insert(0, str(ROOT / "experiments/096-mse-floor"))
        sys.path.insert(0, str(ROOT / "src"))
        import joint187
        mdl = joint187.Model(joint187.F_FIT)
        for c in ctxs:
            cl187 = float(np.sqrt(mdl.t[c]["claim"] / dec[c]["B"]))
            dev = res[c]["cliff"] / cl187 - 1
            print(f"  KA {c}: NTC-proxy cliff {res[c]['cliff']:.4f} vs 187-claim cliff {cl187:.4f} ({dev:+.4%})", flush=True)
            assert abs(dev) <= 0.005, f"KA cliff {c} off by {dev:+.4%}"
        print("  KA PASS (Cp46 == 180's cp_all; cliffs == 187's to <= 0.5 %)", flush=True)
    print("VCC_BULK_CTX_SCALE (relative to THIS artifact) = " + ",".join(f"{c}={scales[c]:.4f}" for c in ctxs))
    out = Path(a.artifact).with_suffix(".ntclevel.json")
    out.write_text(json.dumps(dict(artifact=str(a.artifact), sd=a.sd, q=q, per_context=res), indent=1))
    print(f"wrote {out} [{time.time() - t0:.0f}s]")


if __name__ == "__main__":
    main()
