"""414: what in the shipped bulk did NOT scale with the per-target key? (the board's pds -0.0036 that PLOLO/H1 miss)

Model, per context, on the scorer-space rows d_p = log1p(5e4 tot_p / sum) - b_ctrl (the 395 and 410b artifacts):
    d^395_p  = s_p + beta          d^410b_p = w_p s_p + beta          w_p = the shipped final scale (key x c, c 0.62410)
so beta, a pattern common to every row that the key did not scale, is the LS solution
    beta = sum_p (1-w_p)(d^410b_p - w_p d^395_p) / sum_p (1-w_p)^2
Reports |beta|^2 against row energies (floor rows w = 0.187), the fraction of the 395->410b change it explains, and
saves beta per context (work/beta414.npz) so PLOLO can emulate it (noise414.py).
Known answer: the 395 decomposition totals reproduce decomp414.json's A_c exactly.

    uv run python shared414.py
"""
from __future__ import annotations

import json
import os
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
os.environ.setdefault("NUMBA_CACHE_DIR", str(HERE / "work" / "numba-cache"))
import numpy as np  # noqa: E402

sys.path.insert(0, str(HERE))
from decomp414 import ARTS, controls, str_col  # noqa: E402
import h5py  # noqa: E402

ROOT = HERE.parents[1]
WORK = HERE / "work"
TS = 50000.0


def rows_of(path, b_cache):
    with h5py.File(path, "r") as h:
        genes = str_col(h, "var/_index")
        ctx_lab = str_col(h, "obs/context")
        tgt_lab = str_col(h, "obs/target_gene")
        indptr = np.asarray(h["X"]["indptr"][:])
        out = {}
        for ctx in "ABC":
            out[ctx] = {}
        lab = list(zip(ctx_lab.tolist(), tgt_lab.tolist()))
        start = 0
        for i in range(1, len(lab) + 1):
            if i == len(lab) or lab[i] != lab[start]:
                lo, hi = int(indptr[start]), int(indptr[i])
                tot = np.bincount(h["X"]["indices"][lo:hi].astype(np.int64),
                                  weights=h["X"]["data"][lo:hi].astype(np.float64), minlength=genes.size)
                out[lab[start][0]][lab[start][1]] = tot
                start = i
    if "b" not in b_cache:
        b_cache["b"] = controls(genes)
    return genes, out


def main() -> None:
    cache = {}
    genes, R0 = rows_of(ARTS["395"], cache)
    _, Rb = rows_of(ARTS["410b"], cache)
    _, Rw = rows_of(ARTS["410"], cache)
    z = np.load(ROOT / "experiments/410-pertscale-revisit/work/pert_scale_410_lo0.3_hi3.npz")
    targets = [str(t) for t in z["targets"]]
    w = np.asarray(z["w"], float) * 0.6240958524547175
    gpos = {g: i for i, g in enumerate(genes)}
    dec = json.loads((WORK / "decomp414.json").read_text())
    res, betas = {}, {}
    for ctx in "ABC":
        b = cache["b"][ctx]
        D = {}
        for tag, R in (("395", R0), ("410b", Rb), ("410", Rw)):
            T = np.vstack([R[ctx][t] for t in targets])
            D[tag] = np.log1p(TS * T / T.sum(1, keepdims=True)) - b[None, :]
        # zero each row's own target coordinate (the scorer drops it)
        for tag in D:
            for i, t in enumerate(targets):
                D[tag][i, gpos[t]] = 0.0
        # KA: A_c (unweighted mean row energy x 300/18532 weights == decomp's A)
        A_chk = float((np.einsum("ij,ij->i", D["395"], D["395"]) / 18532).sum())
        ka = abs(A_chk - dec["395"][ctx]["A"])
        out = dict(ka_A=ka)
        for tag in ("410b", "410"):
            num = ((1 - w)[:, None] * (D[tag] - w[:, None] * D["395"])).sum(0)
            beta = num / float(((1 - w) ** 2).sum())
            fit = w[:, None] * (D["395"] - beta[None, :]) + beta[None, :]
            change = D[tag] - w[:, None] * D["395"]          # what a pure per-row scale leaves
            e_change = float((change ** 2).sum())
            e_resid = float(((D[tag] - fit) ** 2).sum())
            floor = w < 0.19
            rowE = np.einsum("ij,ij->i", D[tag], D[tag])
            bb = float(beta @ beta)
            # cosine of each row with beta (how much of a row IS the shared pattern)
            cosb = (D[tag] @ beta) / np.maximum(np.linalg.norm(D[tag], axis=1) * np.sqrt(bb), 1e-300)
            cosb0 = (D["395"] @ beta) / np.maximum(np.linalg.norm(D["395"], axis=1) * np.sqrt(bb), 1e-300)
            # direct per-row scale fit a_p (no beta) and its cosine residual
            a = np.einsum("ij,ij->i", D[tag], D["395"]) / np.einsum("ij,ij->i", D["395"], D["395"])
            cos_same = np.einsum("ij,ij->i", D[tag], D["395"]) / (np.linalg.norm(D[tag], axis=1) * np.linalg.norm(D["395"], axis=1))
            out[tag] = dict(beta_energy=bb, beta_energy_over_floor_row=bb / float(np.median(rowE[floor])),
                            beta_energy_over_median_row=bb / float(np.median(rowE)),
                            frac_nonscaling_change_explained=1 - e_resid / e_change,
                            cos_row_beta_floor_median=float(np.median(cosb[floor])),
                            cos_row_beta_all_median=float(np.median(cosb)),
                            cos_row_beta_395_floor_median=float(np.median(cosb0[floor])),
                            a_fit_over_w_median_floor=float(np.median(a[floor] / w[floor])),
                            a_fit_over_w_median_rest=float(np.median(a[~floor] / w[~floor])),
                            cos_with_395_floor_median=float(np.median(cos_same[floor])),
                            cos_with_395_rest_median=float(np.median(cos_same[~floor])),
                            n_floor=int(floor.sum()))
            if tag == "410b":
                betas[ctx] = beta
            print(f"{ctx} {tag}: |beta|^2 {bb:.4f} = {bb / np.median(rowE[floor]):.2f} x floor-row energy "
                  f"({np.median(rowE[floor]):.3f}), {bb / np.median(rowE):.3f} x median row; explains "
                  f"{100 * (1 - e_resid / e_change):.1f}% of the non-scaling change; cos(row, beta) floor median "
                  f"{np.median(cosb[floor]):.3f} (395: {np.median(cosb0[floor]):.3f}); cos(410b row, 395 row) floor "
                  f"{np.median(cos_same[floor]):.4f} rest {np.median(cos_same[~floor]):.4f}; a/w floor "
                  f"{np.median(a[floor] / w[floor]):.3f} rest {np.median(a[~floor] / w[~floor]):.3f}  [KA A {ka:.1e}]",
                  flush=True)
        res[ctx] = out
    np.savez(WORK / "beta414.npz", genes=np.asarray(genes), **{f"beta_{c}": betas[c] for c in betas})
    (WORK / "shared414.json").write_text(json.dumps(res, indent=1))


if __name__ == "__main__":
    main()
