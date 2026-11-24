"""Extract per-perturbation response profiles from a donor Perturb-seq dataset.

The 2026 task is zero-shot transfer: we never see a perturbation response in the target
context, so every prediction has to start from a response observed *somewhere else*.
Replogle K562-gwps is the only public screen covering the panel (272/300), so its
per-perturbation effect vectors are the core modelling asset.

What this writes, per donor:
  delta.npy       (n_perts, n_genes) float32 — log-pseudobulk effect vs control,
                  b_p - b_ctrl with b = log1p(TS * P / sum(P)), TS = 5e4.
                  Same space `expr_mse` and `pds_cosine` are computed in.
  b_ctrl.npy      (n_genes,) float32 — the donor's own control profile
  perts.txt       perturbation labels, row order of delta
  genes.txt       gene symbols, column order of delta
  meta.json       cell counts per perturbation, donor identity, parameters

Effects are stored on the DONOR's own gene axis. Aligning to the VCC axis is a separate
step, deliberately: the alignment is lossy (gwps carries 8,248 genes, 7,679 of the VCC
18,533) and which genes get dropped is a modelling decision, not a preprocessing one.

    uv run python scripts/build_donor_table.py --donor k562_gwps
    uv run python scripts/build_donor_table.py --donor all
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import anndata as ad
import h5py
import numpy as np
import scipy.sparse as sp

TS = 5e4

DONORS = {
    "k562_gwps": {
        "path": "data/external/scperturb/ReplogleWeissman2022_K562_gwps.h5ad",
        "control": "control",
        "label_col": "perturbation",
        "gene_col": None,
    },
    "k562_essential": {
        "path": "data/external/scperturb/ReplogleWeissman2022_K562_essential.h5ad",
        "control": "control",
        "label_col": "perturbation",
        "gene_col": None,
    },
    "rpe1": {
        "path": "data/external/scperturb/ReplogleWeissman2022_rpe1.h5ad",
        "control": "control",
        "label_col": "perturbation",
        "gene_col": None,
    },
    "nadig_jurkat": {
        "path": "data/external/nadig2025/jurkat_raw_singlecell.h5ad",
        "control": "non-targeting",
        "label_col": "gene",
        "gene_col": "gene_name",
    },
    "nadig_hepg2": {
        "path": "data/external/nadig2025/hepg2_raw_singlecell.h5ad",
        "control": "non-targeting",
        "label_col": "gene",
        "gene_col": "gene_name",
    },
}


def group_sums(path: Path, codes: np.ndarray, n_groups: int, n_genes: int,
               block: int = 16_384) -> np.ndarray:
    """Sum rows by group, streaming the matrix in row blocks.

    scPerturb stores K562-gwps as a DENSE gzipped array: 1,989,578 x 8,248 float32 is
    65.6 GB uncompressed, so it cannot be materialised on a 64 GB machine. We read row
    blocks and accumulate with a sparse indicator matmul, which keeps peak memory at
    one block (~0.5 GB) regardless of dataset size and stays BLAS-fast.
    """
    sums = np.zeros((n_groups, n_genes), dtype=np.float64)
    with h5py.File(path, "r") as f:
        X = f["X"]
        dense = isinstance(X, h5py.Dataset)
        n_rows = X.shape[0] if dense else X.attrs["shape"][0]
        for lo in range(0, n_rows, block):
            hi = min(lo + block, n_rows)
            blk = np.asarray(X[lo:hi]) if dense else _read_csr_block(X, lo, hi, n_genes)
            c = codes[lo:hi]
            keep = c >= 0
            if not keep.any():
                continue
            ind = sp.csr_matrix(
                (np.ones(int(keep.sum()), dtype=np.float32),
                 (c[keep], np.flatnonzero(keep))),
                shape=(n_groups, hi - lo),
            )
            sums += ind @ blk
            if (lo // block) % 20 == 0:
                print(f"    {hi:,}/{n_rows:,} rows", flush=True)
    return sums


def _read_csr_block(X, lo: int, hi: int, n_genes: int) -> np.ndarray:
    indptr = X["indptr"][lo:hi + 1]
    a, b = int(indptr[0]), int(indptr[-1])
    data, idx = X["data"][a:b], X["indices"][a:b]
    return sp.csr_matrix((data, idx, indptr - a), shape=(hi - lo, n_genes)).toarray()


def build(path: Path, control_label: str, out: Path, min_cells: int = 20,
          label_col: str = "perturbation", gene_col: str | None = None) -> dict:
    a = ad.read_h5ad(path, backed="r")
    labels = a.obs[label_col].astype(str).to_numpy()
    genes = np.asarray(a.var[gene_col] if gene_col else a.var_names, dtype=object)
    n_genes = a.n_vars
    print(f"  {a.n_obs:,} cells x {n_genes:,} genes", flush=True)

    uniq, counts = np.unique(labels, return_counts=True)
    if control_label not in set(uniq):
        raise SystemExit(f"no cells labelled {control_label!r}; saw {list(uniq[:8])}")

    # Keep control plus every perturbation with enough cells for a stable pseudobulk.
    keep_names = [n for n, c in zip(uniq, counts)
                  if n == control_label or c >= min_cells]
    skipped = len(uniq) - len(keep_names)
    index = {n: i for i, n in enumerate(keep_names)}
    codes = np.fromiter((index.get(l, -1) for l in labels), dtype=np.int64, count=len(labels))

    print(f"  {len(keep_names)} groups ({skipped} dropped for <{min_cells} cells)", flush=True)
    sums = group_sums(path, codes, len(keep_names), n_genes)

    # Some public files use Ensembl IDs as var_names and carry gene symbols in
    # a column; symbols can be duplicated (Nadig has two HSPA14 records). The
    # downstream model aligns by symbol, so collapse duplicates here instead
    # of silently keeping whichever duplicate a dict happens to retain.
    uniq_genes, inverse = np.unique(genes.astype(str), return_inverse=True)
    if len(uniq_genes) != len(genes):
        mapper = sp.csr_matrix(
            (np.ones(len(genes)), (np.arange(len(genes)), inverse)),
            shape=(len(genes), len(uniq_genes)),
        )
        sums = np.asarray(sums @ mapper)
        genes = uniq_genes.astype(object)
        n_genes = len(genes)
        print(f"  collapsed {len(inverse) - n_genes} duplicate gene symbols", flush=True)

    ci = index[control_label]
    Pc = sums[ci]
    b_ctrl = np.log1p(TS * Pc / Pc.sum()).astype(np.float32)

    tot = sums.sum(axis=1, keepdims=True)
    good = (tot.ravel() > 0)
    good[ci] = False
    prof = np.zeros_like(sums, dtype=np.float32)
    prof[good] = np.log1p(TS * sums[good] / tot[good]).astype(np.float32)
    delta = prof[good] - b_ctrl
    names = [n for i, n in enumerate(keep_names) if good[i]]
    cells = [int(counts[list(uniq).index(n)]) for n in names]

    out.mkdir(parents=True, exist_ok=True)
    np.save(out / "delta.npy", delta)
    np.save(out / "b_ctrl.npy", b_ctrl)
    (out / "perts.txt").write_text("\n".join(names))
    (out / "genes.txt").write_text("\n".join(map(str, genes)))

    norms = np.linalg.norm(delta, axis=1)
    meta = dict(
        source=str(path), control_label=control_label, target_sum=TS,
        n_perturbations=len(names), n_genes=int(n_genes),
        n_control_cells=int(counts[list(uniq).index(control_label)]),
        min_cells=min_cells, skipped_low_cells=int(skipped),
        cells_per_pert=dict(min=int(min(cells)), median=int(np.median(cells)),
                            max=int(max(cells))),
        effect_norm=dict(min=float(norms.min()), median=float(np.median(norms)),
                         max=float(norms.max())),
    )
    (out / "meta.json").write_text(json.dumps(meta, indent=2))
    print(f"  wrote {out}: delta {delta.shape}")
    print(f"  effect norm: median {np.median(norms):.3f}  "
          f"min {norms.min():.3f}  max {norms.max():.3f}")
    return meta


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--donor", default="k562_gwps", choices=[*DONORS, "all"])
    ap.add_argument("--outdir", default="data/donors")
    ap.add_argument("--min-cells", type=int, default=20,
                    help="drop perturbations with fewer cells; their pseudobulk is too noisy")
    args = ap.parse_args()

    todo = list(DONORS) if args.donor == "all" else [args.donor]
    for name in todo:
        cfg = DONORS[name]
        path, ctrl = cfg["path"], cfg["control"]
        if not Path(path).exists():
            print(f"skip {name}: {path} missing")
            continue
        print(f"\n=== {name} ===")
        build(
            Path(path),
            ctrl,
            Path(args.outdir) / name,
            args.min_cells,
            label_col=cfg["label_col"],
            gene_col=cfg["gene_col"],
        )


if __name__ == "__main__":
    main()
