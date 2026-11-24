"""195 — project every object of one LODO fold into that fold's eigenbasis, once.

`W_prod = V diag((ev+lam)^(-gamma/2)) V^T` and `V` is orthogonal, so a cosine between a
whitened query and a raw gallery vector is unchanged by writing both in `V`'s coordinates:

    Yq = Q V,  Yt = T V   =>   <Q W_prod, t> = <Yq * w(gamma), Yt>,  ||Q W_prod|| = ||Yq * w||

That turns a `gamma` sweep, a rank sweep, a learning curve and a permutation null into
elementwise rescalings of two cached matrices instead of a 7,760 x 7,760 matmul each time.
The top-k eigen-subspace `E_k` is then just the last k coordinates, so the adapter's reduced
statistics fall out with no extra projection.

Training queries are strictly cross-donor and strictly non-panel: for the fold that holds out
line H, one non-held donor's gene_norm'd effect is the query and the *other* non-held donor's
raw effect for the same perturbation is the answer, over the perturbations both measure minus
every VCC target gene.  The held line contributes nothing -- not to the queries, not to the
galleries, not to the operator, not to the basis.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from harness import (LODO, Bank, build_fold, common_perts, fold_basis_rows,  # noqa: E402
                     load_banks, sha256_arr)

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
OPS = HERE / "work" / "operators"
PREP = HERE / "work" / "prep"
BLOCK = 2048


def project(X: np.ndarray, V: np.ndarray, out: np.ndarray) -> None:
    for i in range(0, X.shape[0], BLOCK):
        out[i:i + BLOCK] = (np.asarray(X[i:i + BLOCK], np.float64) @ V).astype(np.float32)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--fold", required=True, choices=LODO)
    args = ap.parse_args()
    held = args.fold
    PREP.mkdir(parents=True, exist_ok=True)

    banks = load_banks()
    common = common_perts()
    brows = fold_basis_rows(banks)
    f = build_fold(held, banks, common, brows[held])
    del f["basis"]

    oppath = OPS / f"op.{held}.npz"
    z = np.load(oppath)
    ev, V = z["evals"], z["evecs"]
    op_sha = sha256_arr(ev, V)
    ref_sha = json.loads((HERE / "work" / "phase0.json").read_text())["operators"][held]["sha256"]
    assert op_sha == ref_sha, f"operator sha mismatch {op_sha} != {ref_sha}"
    print(f"operator {oppath.relative_to(ROOT)}  sha256 {op_sha}", flush=True)

    others = f["others"]
    vcc = set(pd.read_csv(ROOT / "data/challenge/pert_counts.csv")["target_gene"])

    # --- panel side (float64 throughout: this is where the identity row must be exact) ---
    panel = {n: f["parts"][n] @ V for n in others}
    panel_T = f["T"] @ V
    np.savez(PREP / f"{held}.panel.npz", T=panel_T, evals=ev,
             **{f"p_{n}": panel[n] for n in others})

    # --- training side: both cross-donor directions among the two non-held donors ---
    meta = {"held": held, "others": others, "operator": {"path": str(oppath.relative_to(ROOT)),
                                                         "sha256": op_sha},
            "ref": f["ref"], "sc": f["sc"], "directions": {}}
    for tgt in others:
        pred = [n for n in others if n != tgt][0]
        shared = [p for p in banks[tgt].perts if p in banks[pred].idx and p not in vcc]
        rt = np.fromiter((banks[tgt].idx[p] for p in shared), np.int64, len(shared))
        rp = np.fromiter((banks[pred].idx[p] for p in shared), np.int64, len(shared))
        n = len(shared)
        t0 = time.time()
        Yq = np.empty((n, V.shape[0]), np.float32)
        Yt = np.empty((n, V.shape[0]), np.float32)
        Q = (f["ref"] / f["sc"][pred]) * banks[pred].gn(rp)
        project(Q, V, Yq)
        del Q
        project(banks[tgt].raw[rt], V, Yt)
        tag = f"{held}.train.{tgt}_from_{pred}"
        np.save(PREP / f"{tag}.Yq.npy", Yq)
        np.save(PREP / f"{tag}.Yt.npy", Yt)
        (PREP / f"{tag}.perts.txt").write_text("\n".join(shared))
        meta["directions"][f"{tgt}_from_{pred}"] = {
            "n_perts": n, "tag": tag,
            "perts_sha256": hashlib.sha256("\n".join(shared).encode()).hexdigest()}
        print(f"  {tag}: {n} queries  ({time.time()-t0:.1f}s)", flush=True)
        del Yq, Yt

    (PREP / f"{held}.meta.json").write_text(json.dumps(meta, indent=2))
    print("wrote", PREP / f"{held}.meta.json")


if __name__ == "__main__":
    main()
