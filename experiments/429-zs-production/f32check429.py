"""429: 423 section 5's float32 robustness check on the MDE carriers. Every DE-universe gene's rank sum recomputed with
per-cell values in float32 (x * 1e6 / lib in float32, log1p in float32, the h1-test controls the same way), before and
after the vendored dn_k250_zs pass, on the dumped seed-0 blocks (base and MDE carriers). Counts genes whose U moves.
usage: uv run python f32check429.py   (reads h1-test; peak RSS ~7 GB)"""
import sys
from pathlib import Path

import numpy as np
import scipy.sparse as sp

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
for _p in (str(ROOT / "experiments/323-current-family-validation"), str(ROOT / "src"), str(HERE)):
    sys.path.insert(0, _p)
import evaluate323 as E  # noqa: E402
import src_rank_mag as RM  # noqa: E402

r, cc, _ = E.load_reference()
del r
lab, genes, xr = E.read_real()
xc = xr[lab == "non-targeting"]
del xr
keep = np.asarray(cc.keep, dtype=bool)
lib_c = np.asarray(xc.sum(1)).ravel().astype(np.float32)
Xc = xc[:, np.flatnonzero(keep)].toarray().astype(np.float32)
del xc
Vc = np.log1p(Xc * (np.float32(1e6) / lib_c[:, None])).T      # (G, n) float32
del Xc
Vc.sort(axis=1)


def U32(block):
    x = block[:, keep].astype(np.float32)
    lib = block.sum(1).astype(np.float32)
    V = np.log1p(x * (np.float32(1e6) / lib[:, None])).T
    out = np.empty(V.shape[0])
    for g in range(V.shape[0]):
        lt = np.searchsorted(Vc[g], V[g], "left")
        rt = np.searchsorted(Vc[g], V[g], "right")
        out[g] = (2 * lt + (rt - lt)).sum()
    return out


tot = {}
for p in sorted((HERE / "work/blocks").glob("base_s0_bps_*.npz")):
    z = np.load(p)
    X = sp.csr_matrix((z["data"], z["indices"], z["indptr"]), shape=tuple(z["shape"])).toarray().astype(np.int64)
    X2, dg = RM.apply(X, int(z["n_car"]), cc, "dn_k250_zs")
    d = U32(X2) - U32(X)
    arm = "MDE" if "mde" in p.stem else "base"
    t = tot.setdefault(arm, [0, 0, 0])
    t[0] += int((d != 0).sum()); t[1] += int(dg["zs_n_moved"]); t[2] += 1
    print(f"{p.stem} ({z['target']}): float32 path, genes with U moved {int((d != 0).sum())} of {int(dg['zs_n_moved'])} moved; "
          f"max |dU| {float(np.abs(d).max())} half-ranks", flush=True)
for arm, (nu, nm, nb) in tot.items():
    print(f"TOTAL {arm}: {nb} blocks, {nm} moved genes, {nu} with a float32-path U change")
