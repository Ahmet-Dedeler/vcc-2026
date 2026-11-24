"""182 -- build ST-HVG inference inputs for H1 (never-seen line) and K562 (in-train line).

Preprocessing reproduces `state tx preprocess_train` exactly:
    sc.pp.normalize_total(adata)   # target_sum = median of per-cell total counts
    sc.pp.log1p(adata)
    obsm["X_hvg"] = adata[:, highly_variable].X
except that the HVG *selection* is replaced by the checkpoint's own frozen 2,000-gene
list (var_dims.pkl['gene_names']), in the checkpoint's order. Genes absent from a
source's axis are emitted as zero columns and recorded.

Writes work/{h1,k562}_in.h5ad with obs columns: gene / cell_line / gem_group.
"""

import json
import pickle
import sys

import anndata as ad
import h5py
import numpy as np
import pandas as pd
from scipy import sparse

SEED = 20260903
MODEL_DIR = "data/external/state-pretrained/ST-HVG-Replogle/zeroshot/hepg2"
OUT = "experiments/182-pretrained-state/work"
N_CTRL = 2048  # control cells used as the basal pool
N_PER_TARGET = 128  # placeholder rows per target (basal is resampled by the CLI)
GEM_GROUP = "1"  # an in-vocabulary gem_group so the batch encoder gets a trained row


def read_csr_rows(f, rows):
    """Read selected rows out of an AnnData CSR /X group without loading the matrix."""
    X = f["X"]
    indptr = X["indptr"][:]
    n_var = int(X.attrs["shape"][1])
    out = np.zeros((len(rows), n_var), dtype=np.float32)
    for i, r in enumerate(rows):
        a, b = int(indptr[r]), int(indptr[r + 1])
        if b > a:
            out[i, X["indices"][a:b]] = X["data"][a:b]
    return out


def cat_codes(f, col):
    g = f["obs"][col]
    cats = [c.decode() if isinstance(c, bytes) else str(c) for c in g["categories"][:]]
    return cats, g["codes"][:]


def load_h1_controls(rng, n):
    p = "data/external/vcc2025/test/adata_Test.h5ad"
    with h5py.File(p, "r") as f:
        cats, codes = cat_codes(f, "target_gene")
        ntc = cats.index("non-targeting")
        idx = np.where(codes == ntc)[0]
        pick = np.sort(rng.choice(idx, size=min(n, len(idx)), replace=False))
        Xc = read_csr_rows(f, pick)
        genes = [x.decode() if isinstance(x, bytes) else str(x) for x in f["var"]["_index"][:]]
    return Xc, genes, len(idx)


def load_k562_controls(rng, n):
    p = "data/external/scperturb/ReplogleWeissman2022_K562_essential.h5ad"
    with h5py.File(p, "r") as f:
        cats, codes = cat_codes(f, "perturbation")
        idx = np.where(codes == cats.index("control"))[0]
        pick = np.sort(rng.choice(idx, size=min(n, len(idx)), replace=False))
        Xc = f["X"][pick, :].astype(np.float32)
        genes = [x.decode() if isinstance(x, bytes) else str(x) for x in f["var"]["gene_name"][:]]
    return Xc, genes, len(idx)


def to_hvg(Xc, src_genes, model_genes, target_sum=None):
    """normalize_total(median) -> log1p -> reindex onto the checkpoint's gene list.

    `target_sum=None` reproduces scanpy's default (this source's own median library
    size). Passing an explicit value lets us match the training line's scale, which
    matters because H1's median is ~4x Replogle's and the model only ever saw the latter.
    """
    tot = Xc.sum(axis=1)
    assert (tot > 0).all(), "control cell with zero counts"
    if target_sum is None:
        target_sum = float(np.median(tot))
    Xn = np.log1p(Xc / tot[:, None] * target_sum)
    pos = {g: i for i, g in enumerate(src_genes)}
    cols = np.array([pos.get(g, -1) for g in model_genes])
    hit = cols >= 0
    H = np.zeros((Xc.shape[0], len(model_genes)), dtype=np.float32)
    H[:, hit] = Xn[:, cols[hit]]
    return H, target_sum, int(hit.sum())


def build(name, Xc, src_genes, model_genes, targets, line, rng, target_sum=None):
    H, target_sum, n_hit = to_hvg(Xc, src_genes, model_genes, target_sum)
    n_ctrl = H.shape[0]
    # placeholder rows: content is irrelevant (the CLI redraws basal from the control
    # pool), but keep them on-manifold by copying random controls.
    blocks = [H]
    labels = ["non-targeting"] * n_ctrl
    for t in targets:
        sel = rng.choice(n_ctrl, size=N_PER_TARGET, replace=True)
        blocks.append(H[sel])
        labels += [t] * N_PER_TARGET
    Hall = np.concatenate(blocks, axis=0)
    obs = pd.DataFrame(
        {
            "gene": pd.Categorical(labels),
            "cell_line": pd.Categorical([line] * len(labels)),
            "gem_group": pd.Categorical([GEM_GROUP] * len(labels)),
        },
        index=[f"c{i}" for i in range(len(labels))],
    )
    A = ad.AnnData(
        X=sparse.csr_matrix((len(labels), 1), dtype=np.float32),
        obs=obs,
        var=pd.DataFrame(index=["_dummy"]),
    )
    A.obsm["X_hvg"] = Hall
    path = f"{OUT}/{name}_in.h5ad"
    A.write_h5ad(path)
    return {
        "path": path,
        "n_control_cells_used": n_ctrl,
        "n_rows": int(len(labels)),
        "n_targets": len(targets),
        "n_per_target": N_PER_TARGET,
        "normalize_target_sum": target_sum,
        "model_genes_present_in_source": n_hit,
        "model_genes_total": len(model_genes),
    }


def main():
    rng = np.random.default_rng(SEED)
    vd = pickle.load(open(f"{MODEL_DIR}/var_dims.pkl", "rb"))
    model_genes = [str(g) for g in vd["gene_names"]]
    vocab = {str(p) for p in vd["pert_names"]}

    # h1-test panel targets that the checkpoint's one-hot vocabulary can express
    d = np.load("experiments/164-callset-artefact/de164.a15.npz", allow_pickle=True)
    h1_targets = sorted({str(p) for p in d["perts"]} & vocab)
    print(f"h1-test targets in ST vocabulary: {len(h1_targets)}/85")

    meta = {"seed": SEED, "model_dir": MODEL_DIR, "targets": h1_targets}

    # K562 first: its own median library size defines the training-line scale that the
    # scale-matched H1 arm is rebuilt against.
    Xk, kgenes, k_avail = load_k562_controls(rng, N_CTRL)
    print(f"K562 controls available {k_avail}, using {Xk.shape}")
    # same 20 targets, so H1 and K562 differ only in the recipient line
    meta["k562"] = build("k562", Xk, kgenes, model_genes, h1_targets, "k562", rng)
    meta["k562"]["n_control_cells_available"] = k_avail
    replogle_ts = meta["k562"]["normalize_target_sum"]
    del Xk

    Xc, genes, n_avail = load_h1_controls(rng, N_CTRL)
    print(f"H1  controls available {n_avail}, using {Xc.shape}")
    meta["h1"] = build("h1", Xc, genes, model_genes, h1_targets, "hepg2", rng)
    meta["h1"]["n_control_cells_available"] = n_avail
    # scale-matched arm: same cells, renormalised to the training line's library size
    meta["h1_scalematched"] = build(
        "h1_scalematched", Xc, genes, model_genes, h1_targets, "hepg2", rng, target_sum=replogle_ts
    )
    meta["h1_scalematched"]["n_control_cells_available"] = n_avail
    del Xc

    json.dump(meta, open(f"{OUT}/prep_meta.json", "w"), indent=2)
    print(json.dumps({k: v for k, v in meta.items() if k in ("h1", "h1_scalematched", "k562")}, indent=2))


if __name__ == "__main__":
    sys.exit(main())
