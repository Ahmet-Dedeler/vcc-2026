"""191 supplementary — `mse`'s DISPERSION half, measured on the two BUILT artifacts.

The digest identity (falsifier 1) pins `mse`'s PROFILE half and `pds` outright.  It says
nothing about the other term.  `expr_mse_unbiased_capped_norm`'s numerator is

    N_p = (1/G_p) ( ||b_hat_p - b_p||^2  -  rho * min(C_hat_p, C_p)  -  C_p )

and `C_hat_p` is the delete-one jackknife of OUR OWN emitted noise, which the rank-controlled
emitter changes on purpose: it replaces the quasi law's gamma-Poisson scatter with two-level
histograms that concentrate each gene's counts into far fewer cells.

185 §6a measured `claim = mean(min(Chat, K C_p)) * weight_sum` on the real challenge panel and
found a15 already SATURATED at `Chat/C_p = 1.085`, with every rank-controlled arm at 1.9-3.3,
so `min(.,.)` stays at the cap and raw `mse` is flat at 0.88228.  That was measured on a
RECONSTRUCTION of the emitted block.  This measures `Chat` on the shipped artifacts.

`C_p` is 096's proxy -- the delete-1 jackknife of each ntc construct's own 400 control cells,
read from the context's own `context_<id>.h5ad`, which is the same estimator the builder's
`VCC_ALPHA=knee` path uses.  No ground truth, no `cell-eval2` process; only the library's
`jackknife_correction` kernel, which the production builder already imports.

    uv run python claim191.py [n_targets_per_context]
"""
from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import anndata as ad
import h5py
import numpy as np
import scipy.sparse as sp

ROOT = Path("~/Code/virtual-cell-challenge-2026")
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT / "src"))
from vcc2026.donor_transfer import TS  # noqa: E402

BUNDLE = ROOT / "data/challenge"
CHUNK = 64
PACKAGES = {
    "a15_local": ROOT / "submissions/191-a15_local.h5ad",
    "rc_hybrid": ROOT / "submissions/191-rc_hybrid.h5ad",
}


def main() -> None:
    from cell_eval2.metrics.delta import PRED_TRACE_CAP_K
    from cell_eval2.moments import jackknife_correction

    n_targets = int(sys.argv[1]) if len(sys.argv) > 1 else 16
    genes = [line.strip() for line in
             (BUNDLE / "gene_names.csv").read_text().splitlines()[1:]]
    n_genes = len(genes)

    obs = {}
    handles = {}
    for label, path in PACKAGES.items():
        adata = ad.read_h5ad(path, backed="r")
        obs[label] = (np.asarray(adata.obs["target_gene"].astype(str).values),
                      np.asarray(adata.obs["context"].astype(str).values))
        adata.file.close()
        handles[label] = h5py.File(path, "r")

    out: dict[str, dict] = {}
    t0 = time.time()
    for context in ("A", "B", "C"):
        controls = ad.read_h5ad(BUNDLE / f"context_{context}.h5ad")
        matrix = controls.X.tocsr() if not sp.isspmatrix_csr(controls.X) else controls.X
        ids = np.asarray(controls.obs["ntc_id"].astype(str).values)
        cps = []
        for name in dict.fromkeys(ids.tolist()):
            rows = np.flatnonzero(ids == name)
            cps.append(float(jackknife_correction(
                sp.csr_matrix(matrix[rows]), np.zeros(rows.size, np.intp), 1, TS,
                chunk=CHUNK)[0]))
        c_p = float(np.mean(cps))
        del controls, matrix
        cap = float(PRED_TRACE_CAP_K) * c_p
        print(f"\nctx {context}: C_p {c_p:.4f} over {len(cps)} ntc constructs "
              f"(sd {np.std(cps):.4f}); cap K*C_p = {cap:.4f}  [{time.time() - t0:.0f}s]",
              flush=True)

        labels_all, contexts_all = obs["a15_local"]
        targets = sorted(set(labels_all[contexts_all == context]))
        chosen = [targets[i] for i in
                  np.linspace(0, len(targets) - 1, n_targets).astype(int)]
        for label in PACKAGES:
            labels, contexts = obs[label]
            handle = handles[label]
            data, indices, indptr = (handle["X"]["data"], handle["X"]["indices"],
                                     handle["X"]["indptr"])
            chats = []
            for target in chosen:
                sel = np.flatnonzero((contexts == context) & (labels == target))
                lo, hi = int(indptr[sel[0]]), int(indptr[sel[-1] + 1])
                block = sp.csr_matrix(
                    (data[lo:hi].astype(np.float64), indices[lo:hi].astype(np.int64),
                     np.asarray(indptr[sel[0]:sel[-1] + 2]) - lo),
                    shape=(sel.size, n_genes),
                )
                chats.append(float(jackknife_correction(
                    block, np.zeros(sel.size, np.intp), 1, TS, chunk=CHUNK)[0]))
            chats = np.asarray(chats)
            saturated = int(np.sum(chats >= cap))
            out.setdefault(label, {})[context] = dict(
                C_p=c_p, cap=cap, chat_mean=float(chats.mean()),
                chat_min=float(chats.min()), chat_max=float(chats.max()),
                ratio_mean=float(chats.mean() / c_p),
                ratio_min=float(chats.min() / c_p),
                claim_per_row=float(np.mean(np.minimum(chats, cap))),
                cap_saturated_rows=saturated, rows=int(chats.size),
            )
            row = out[label][context]
            print(f"  {label:<11} Chat mean {row['chat_mean']:8.3f}  "
                  f"Chat/C_p mean {row['ratio_mean']:6.3f} min {row['ratio_min']:6.3f}   "
                  f"min(Chat, cap) mean {row['claim_per_row']:8.4f}   "
                  f"rows at the cap {saturated}/{row['rows']}", flush=True)

    for handle in handles.values():
        handle.close()

    print("\nREADING")
    for label in PACKAGES:
        claim = np.mean([out[label][c]["claim_per_row"] for c in "ABC"])
        ratio = np.mean([out[label][c]["ratio_mean"] for c in "ABC"])
        sat = sum(out[label][c]["cap_saturated_rows"] for c in "ABC")
        rows = sum(out[label][c]["rows"] for c in "ABC")
        print(f"  {label:<11} mean Chat/C_p {ratio:6.3f}   mean min(Chat, K C_p) "
              f"{claim:8.4f}   rows at the cap {sat}/{rows}")
    same = all(
        abs(out["a15_local"][c]["claim_per_row"] - out["rc_hybrid"][c]["claim_per_row"])
        < 1e-9 for c in "ABC"
    )
    print(f"  min(Chat, K C_p) identical across the two arms on every context: {same}")
    with open(HERE / "work" / "claim191.json", "w") as handle:
        json.dump(out, handle, indent=2)


if __name__ == "__main__":
    main()
