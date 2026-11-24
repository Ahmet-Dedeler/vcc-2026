"""Build a donor table from a scPerturb-format h5ad whose `.X` is stored CSC.

`scripts/build_donor_table.py` streams row blocks, which needs CSR (or dense). The scPerturb
mirrors of Frangieh / Xu / Tian are `csc_matrix`, where a row block would touch every column
chunk. This walks GENE blocks instead and accumulates the same group sums, so peak memory is
one gene block regardless of dataset size.

Everything downstream of the group sums is byte-identical in construction to
`scripts/build_donor_table.py`: same `TS = 5e4`, same `b = log1p(TS * P / sum(P))`, same
`delta = b_pert - b_ctrl`, same five output files. The table therefore satisfies the
experiment-161 invariant `sum(expm1(b_ctrl)) == sum(expm1(b_ctrl + delta)) == TS` and loads
through `src/vcc2026/donor_transfer.py:Donor.load`.

`--stratum-col/--stratum` restricts to one experimental arm. Frangieh 2021 pools three arms
(Control / IFN-gamma / Co-culture) in one file; mixing them would put the stimulation response
into the control profile and into every effect vector, so only `Control` is built.

    uv run python build_extra.py --source frangieh_melanoma
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import h5py
import numpy as np
import scipy.sparse as sp

REPO = Path(__file__).resolve().parents[2]
TS = 5e4

SOURCES = {
    "frangieh_melanoma": dict(
        path="data/external/scperturb/FrangiehIzar2021_RNA.h5ad",
        control="control", label_col="perturbation",
        stratum_col="perturbation_2", stratum="Control",
        tissue="melanoma (patient-derived, Perturb-CITE-seq)",
    ),
    "tian_ipsc_neuron": dict(
        path="data/external/scperturb/TianKampmann2021_CRISPRi.h5ad",
        control="control", label_col="perturbation",
        stratum_col=None, stratum=None,
        tissue="iPSC-derived neuron",
    ),
    "xucao_hek293": dict(
        path="data/external/scperturb/XuCao2023.h5ad",
        control="control", label_col="perturbation",
        stratum_col=None, stratum=None,
        tissue="kidney (HEK293)",
    ),
}


def _decode(a):
    if a.dtype.kind in "SO":
        return np.asarray([x.decode() if isinstance(x, bytes) else str(x) for x in a])
    return np.asarray(a)


def read_col(h, name):
    g = h["obs"][name]
    if isinstance(g, h5py.Group):
        cats = _decode(g["categories"][:])
        return cats[np.clip(g["codes"][:], 0, len(cats) - 1)]
    return _decode(g[:])


def var_names(h):
    v = h["var"]
    key = v.attrs.get("_index", "_index")
    if isinstance(key, bytes):
        key = key.decode()
    return _decode(v[key][:])


def group_sums_csc(h, codes: np.ndarray, n_groups: int, n_genes: int,
                   block: int = 2000) -> np.ndarray:
    """Sum rows by group from a CSC `.X`, one gene block at a time."""
    X = h["X"]
    indptr = X["indptr"][:]
    keep = codes >= 0
    ind = sp.csr_matrix(
        (np.ones(int(keep.sum()), dtype=np.float32),
         (codes[keep], np.flatnonzero(keep))),
        shape=(n_groups, len(codes)),
    )
    sums = np.zeros((n_groups, n_genes), dtype=np.float64)
    for g0 in range(0, n_genes, block):
        g1 = min(g0 + block, n_genes)
        a, b = int(indptr[g0]), int(indptr[g1])
        blk = sp.csc_matrix(
            (X["data"][a:b], X["indices"][a:b], indptr[g0:g1 + 1] - a),
            shape=(len(codes), g1 - g0),
        )
        sums[:, g0:g1] = (ind @ blk).toarray()
        print(f"    genes {g1:,}/{n_genes:,}", flush=True)
    return sums


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--source", required=True, choices=list(SOURCES))
    ap.add_argument("--outdir", default=str(Path(__file__).resolve().parent / "donors"))
    ap.add_argument("--min-cells", type=int, default=20)
    args = ap.parse_args()

    cfg = SOURCES[args.source]
    path = REPO / cfg["path"]
    out = Path(args.outdir) / args.source

    with h5py.File(path, "r") as h:
        labels = read_col(h, cfg["label_col"])
        genes = var_names(h)
        n_genes = len(genes)
        n_cells = len(labels)
        print(f"  {n_cells:,} cells x {n_genes:,} genes")

        in_stratum = np.ones(n_cells, dtype=bool)
        if cfg["stratum_col"]:
            s = read_col(h, cfg["stratum_col"])
            in_stratum = s == cfg["stratum"]
            print(f"  stratum {cfg['stratum_col']}=={cfg['stratum']!r}: "
                  f"{int(in_stratum.sum()):,} cells")

        lab = np.where(in_stratum, labels, "")
        uniq, counts = np.unique(lab[in_stratum], return_counts=True)
        cnt = dict(zip(uniq.tolist(), counts.tolist()))
        ctrl = cfg["control"]
        if ctrl not in cnt:
            raise SystemExit(f"no {ctrl!r} cells in stratum; saw {list(cnt)[:8]}")
        keep_names = [n for n in uniq.tolist()
                      if n == ctrl or cnt[n] >= args.min_cells]
        skipped = len(uniq) - len(keep_names)
        index = {n: i for i, n in enumerate(keep_names)}
        codes = np.fromiter((index.get(l, -1) for l in lab), dtype=np.int64, count=n_cells)
        print(f"  {len(keep_names)} groups ({skipped} dropped for <{args.min_cells} cells)")

        sums = group_sums_csc(h, codes, len(keep_names), n_genes)

    # collapse duplicate gene symbols exactly as the shipped builder does
    uniq_genes, inverse = np.unique(genes.astype(str), return_inverse=True)
    if len(uniq_genes) != len(genes):
        mapper = sp.csr_matrix((np.ones(len(genes)), (np.arange(len(genes)), inverse)),
                               shape=(len(genes), len(uniq_genes)))
        sums = np.asarray(sums @ mapper)
        genes = uniq_genes.astype(object)
        n_genes = len(genes)
        print(f"  collapsed {len(inverse) - n_genes} duplicate gene symbols")

    ci = index[ctrl]
    Pc = sums[ci]
    b_ctrl = np.log1p(TS * Pc / Pc.sum()).astype(np.float32)

    tot = sums.sum(axis=1, keepdims=True)
    good = tot.ravel() > 0
    good[ci] = False
    prof = np.zeros_like(sums, dtype=np.float32)
    prof[good] = np.log1p(TS * sums[good] / tot[good]).astype(np.float32)
    delta = prof[good] - b_ctrl
    names = [n for i, n in enumerate(keep_names) if good[i]]
    cells = [cnt[n] for n in names]

    out.mkdir(parents=True, exist_ok=True)
    np.save(out / "delta.npy", delta)
    np.save(out / "b_ctrl.npy", b_ctrl)
    (out / "perts.txt").write_text("\n".join(names))
    (out / "genes.txt").write_text("\n".join(map(str, genes)))
    norms = np.linalg.norm(delta, axis=1)
    (out / "meta.json").write_text(json.dumps(dict(
        source=str(cfg["path"]), tissue=cfg["tissue"], control_label=ctrl,
        stratum={cfg["stratum_col"]: cfg["stratum"]} if cfg["stratum_col"] else None,
        target_sum=TS, n_perturbations=len(names), n_genes=int(n_genes),
        n_control_cells=int(cnt[ctrl]), min_cells=args.min_cells,
        skipped_low_cells=int(skipped),
        cells_per_pert=dict(min=int(min(cells)), median=int(np.median(cells)),
                            max=int(max(cells))),
        effect_norm=dict(min=float(norms.min()), median=float(np.median(norms)),
                         max=float(norms.max())),
    ), indent=2))
    print(f"  wrote {out}: delta {delta.shape}")
    print(f"  effect norm median {np.median(norms):.3f}  "
          f"cells/pert median {np.median(cells):.0f}")

    import pandas as pd
    panel = set(pd.read_csv(REPO / "data/challenge/pert_counts.csv")["target_gene"])
    hit = panel & set(names)
    print(f"  VCC panel coverage: {len(hit)}/300")


if __name__ == "__main__":
    main()
