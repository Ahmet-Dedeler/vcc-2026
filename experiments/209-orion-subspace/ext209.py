"""209 — the dimension-sensitivity basis.

`subspace209.py` stores only the `r` directions its pre-declared criterion admits, so the
`subC_r5` / `subC_r10` arms in the first scoring pass silently re-used the rank-2 basis and
returned rank-2 numbers.  This rebuilds the SAME estimator (same rows, same columns, same
residualisation) and stores the top 20 eigenvectors, so dimension can be swept honestly.

`subCext_r2` must reproduce `subC_l1` exactly -- that is this script's known-answer row.
"""
from __future__ import annotations

import json
import os
import sys
from pathlib import Path

os.environ.setdefault("OPENBLAS_NUM_THREADS", "6")
os.environ.setdefault("OMP_NUM_THREADS", "6")

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
WORK = HERE / "work"
sys.path.insert(0, str(ROOT / "experiments/162-retrodicting-validator"))
sys.path.insert(0, str(ROOT / "experiments/204-repression-weighted-mixing"))
sys.path.insert(0, str(HERE))

import plolo  # noqa: E402

plolo.WORK = WORK / "plolo"
import prep204  # noqa: E402

prep204.WORK = WORK
from prep204 import prepared, shared_targets  # noqa: E402
from plolo import PANEL_SET  # noqa: E402
from subspace209 import excess_subspace  # noqa: E402

DONORS = ("k562_gwps", "orion_hct116", "orion_hek293t")
KEEP = 20


def main() -> None:
    tg = shared_targets()
    est = np.asarray([i for i, t in enumerate(tg) if t not in PANEL_SET])
    P = {d: np.asarray(prepared(d, tg, "shared")["effects"], dtype=np.float64)
         for d in DONORS}
    support = {d: (np.abs(P[d]).sum(0) > 0) for d in DONORS}
    C = np.flatnonzero(support["k562_gwps"] & support["orion_hct116"]
                       & support["orion_hek293t"])
    res = excess_subspace(P["k562_gwps"][np.ix_(est, C)],
                          P["orion_hct116"][np.ix_(est, C)],
                          P["orion_hek293t"][np.ix_(est, C)],
                          n_perm=1, kmax=KEEP)
    lam = res["eigvals"]
    # re-derive the top-KEEP eigenvectors regardless of the admitted r
    n = len(est)
    nk = (P["k562_gwps"][np.ix_(est, C)] ** 2).sum(1)
    Xk = P["k562_gwps"][np.ix_(est, C)]
    Eh = P["orion_hct116"][np.ix_(est, C)] - (
        (P["orion_hct116"][np.ix_(est, C)] * Xk).sum(1) / np.maximum(nk, 1e-300))[:, None] * Xk
    Ee = P["orion_hek293t"][np.ix_(est, C)] - (
        (P["orion_hek293t"][np.ix_(est, C)] * Xk).sum(1) / np.maximum(nk, 1e-300))[:, None] * Xk
    A = np.vstack([Eh, Ee])
    U, sig, Wt = np.linalg.svd(A, full_matrices=False)
    keep = sig > sig[0] * 1e-10
    U, sig, Wt = U[:, keep], sig[keep], Wt[keep]
    B = U * sig
    K = (B[:n].T @ B[n:]) / n
    K = 0.5 * (K + K.T)
    ev, Z = np.linalg.eigh(K)
    order = np.argsort(ev)[::-1]
    V = (Wt.T @ Z[:, order[:KEEP]]).T

    zc = np.load(WORK / "subspace-C.npz")
    ang = np.linalg.svd(V[:2] @ zc["basis"].T, compute_uv=False)
    print(f"eigenvalues top-{KEEP}: " + " ".join(f"{v:.4g}" for v in ev[order[:KEEP]]))
    print(f"principal cosines of the top-2 against subspace-C: {ang.round(8).tolist()}")
    np.savez_compressed(WORK / "subspace-Cext.npz", basis=V, cols=C,
                        eigvals=ev[order[:KEEP]], r=int(zc["r"]),
                        thresh=float(zc["thresh"]))
    (WORK / "ext209.json").write_text(json.dumps(dict(
        eigvals_top20=ev[order[:KEEP]].tolist(),
        principal_cosines_vs_subspaceC=ang.tolist(),
        thresh=float(zc["thresh"]), r_declared=int(zc["r"]),
        n_above_thresh=int((ev[order] > float(zc["thresh"])).sum()),
        n_above_null_median=int((ev[order] > json.loads(
            (WORK / "subspace209.json").read_text())["variants"]["C"]
            ["null_top_median"]).sum())), indent=2))
    print("wrote work/subspace-Cext.npz and work/ext209.json")


if __name__ == "__main__":
    main()
