#!/usr/bin/env python3
"""383 — build STATE inference templates (log1p, support 18,080-gene order).

Reads challenge/harness/support files read-only; writes only to this experiment's
work/templates/. Inference itself uses the TSV-pad path (`--tsv`), so templates
contain CONTROL cells only; perturbations are cloned from controls at infer time.

Subcommands:
  identity  val-template controls subset (native log1p) + TSV for the 50 val perts
  h1        h1-test.h5ad controls log1p'd + TSV for the 85 H1 targets
  ctx       context_{A,B,C}.h5ad controls log1p'd + TSV for the 300 panel targets
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import h5py
import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
TDIR = HERE / "work" / "templates"
SUPPORT_GENES = HERE / "data" / "competition_support_set" / "gene_names.csv"
CHALLENGE_GENES = ROOT / "data" / "challenge" / "gene_names.csv"

sys.path.insert(0, str(ROOT / "src"))


def support_gene_list() -> list[str]:
    return [l.strip() for l in SUPPORT_GENES.open() if l.strip()]


def read_cats(h: h5py.File, path: str) -> np.ndarray:
    o = h[path]
    if isinstance(o, h5py.Group) and "categories" in o:
        return o["categories"].asstr()[:][o["codes"][:]].astype(str)
    return o.asstr()[:].astype(str)


def read_block(h: h5py.File, rows: np.ndarray, n_cols: int) -> np.ndarray:
    X = h["X"]
    if isinstance(X, h5py.Dataset):  # dense (e.g. val template, already log1p)
        rows = np.asarray(rows)
        order = np.argsort(rows, kind="stable")
        got = np.asarray(X[np.sort(rows)], dtype=np.float32)
        out = np.empty_like(got)
        out[order] = got
        return out
    out = np.zeros((len(rows), n_cols), dtype=np.float32)
    indptr = X["indptr"][:]
    data = X["data"]
    idx = X["indices"]
    for i, r in enumerate(rows):
        s, e = int(indptr[r]), int(indptr[r + 1])
        out[i, idx[s:e]] = data[s:e]
    return out


def read_sparse_block(h: h5py.File, rows: np.ndarray, n_cols: int) -> np.ndarray:
    return read_block(h, rows, n_cols)


def write_h5ad(path: Path, X: np.ndarray, obs: dict, var_names: list[str]) -> None:
    import anndata as ad

    import pandas as pd

    a = ad.AnnData(X=X, obs=pd.DataFrame(obs),
                   var=pd.DataFrame(index=pd.Index(var_names, name="_index")))
    a.write_h5ad(path)


def write_tsv(path: Path, perts: list[str], n: int) -> None:
    with path.open("w") as f:
        f.write("perturbation\tnum_cells\n")
        for p in perts:
            f.write(f"{p}\t{n}\n")


VAL_TEMPLATE = HERE / "data" / "competition_support_set-log1p" / "competition_val_template.h5ad"
H1_TEST = ROOT / "data" / "harness" / "h1-test.h5ad"


def cmd_identity(args) -> None:
    TDIR.mkdir(parents=True, exist_ok=True)
    genes = support_gene_list()
    with h5py.File(VAL_TEMPLATE, "r") as h:
        tg = read_cats(h, "obs/target_gene")
        ctl = np.where(tg == "non-targeting")[0]
        perts = sorted(set(tg) - {"non-targeting"})
        assert len(perts) == 50, perts
        rng = np.random.RandomState(383)
        take = np.sort(rng.choice(ctl, size=args.controls, replace=False))
        X = read_sparse_block(h, take, len(genes))
        print("val template X space check: max", float(X.max()))
    write_h5ad(TDIR / "identity_controls.h5ad", X,
               {"target_gene": ["non-targeting"] * len(take),
                "cell_type": ["H1"] * len(take)}, genes)
    write_tsv(TDIR / "identity.tsv", perts, args.cells_per_pert)
    print("identity:", len(take), "controls,", len(perts), "perts x", args.cells_per_pert)


def log1p_controls(h5_path: str, genes: list[str], var_names: list[str]) -> tuple[np.ndarray, np.ndarray]:
    """Return (log1p control matrix in support-gene order, control row idx)."""
    chal = [l.strip() for l in open(CHALLENGE_GENES)][1:]
    col_of = {g: j for j, g in enumerate(chal)}
    cols = np.array([col_of.get(g, -1) for g in genes], dtype=np.int64)
    n_missing = int((cols < 0).sum())
    print("support genes missing from challenge space:", n_missing,
          [g for g, c in zip(genes, cols) if c < 0][:10])
    with h5py.File(h5_path, "r") as h:
        var = read_cats(h, "var/_index")
        if not np.array_equal(var, np.array(chal)):
            raise ValueError("challenge h5ad var order != gene_names.csv")
        tg = read_cats(h, "obs/target_gene")
        ctl_rows = np.where(tg == "non-targeting")[0]
        print(Path(h5_path).name, "controls:", len(ctl_rows), "of", len(tg))
        import scipy.sparse as sp

        X = h["X"]
        indptr, data, idx = X["indptr"][:], X["data"], X["indices"]
        n = len(ctl_rows)
        out = np.zeros((n, len(genes)), dtype=np.float32)
        have = cols >= 0
        want = cols[have]
        # gather raw counts for wanted columns, chunked over rows
        for a in range(0, n, 2000):
            rows = ctl_rows[a:a + 2000]
            for i, r in enumerate(rows):
                s, e = int(indptr[r]), int(indptr[r + 1])
                row_idx, row_dat = idx[s:e], data[s:e]
                m = {c: v for c, v in zip(row_idx.tolist(), row_dat.tolist())}
                out[a + i, have] = [m.get(int(c), 0) for c in want]
            if a % 10000 == 0:
                print("  rows", a, "/", n, flush=True)
    return np.log1p(out), ctl_rows


def cmd_h1(args) -> None:
    TDIR.mkdir(parents=True, exist_ok=True)
    genes = support_gene_list()
    with h5py.File(H1_TEST, "r") as h:
        var = read_cats(h, "var/_index")
        print("h1 var == support gene order:", bool(np.array_equal(var, np.array(genes))))
        tg = read_cats(h, "obs/target_gene")
        perts = sorted(set(tg) - {"non-targeting"})
        assert len(perts) == 85, len(perts)
    Xlog, _ = log1p_h1_controls(H1_TEST, genes, args.controls)
    write_h5ad(TDIR / "h1_controls.h5ad", Xlog,
               {"target_gene": ["non-targeting"] * len(Xlog),
                "cell_type": ["H1"] * len(Xlog)}, genes)
    write_tsv(TDIR / "h1.tsv", perts, args.cells_per_pert)
    print("h1:", len(Xlog), "controls,", len(perts), "perts x", args.cells_per_pert)


def log1p_h1_controls(h5_path: str, genes: list[str], n_take: int) -> tuple[np.ndarray, list[str]]:
    with h5py.File(h5_path, "r") as h:
        tg = read_cats(h, "obs/target_gene")
        ctl = np.where(tg == "non-targeting")[0]
        rng = np.random.RandomState(383)
        take = np.sort(rng.choice(ctl, size=min(n_take, len(ctl)), replace=False))
        raw = read_sparse_block(h, take, len(genes)).astype(np.float64)
    return np.log1p(raw).astype(np.float32), take.tolist()


def cmd_ctx(args) -> None:
    TDIR.mkdir(parents=True, exist_ok=True)
    genes = support_gene_list()
    panel = [l.strip() for l in (ROOT / "data" / "challenge" / "pert_counts.csv").open()][1:]
    assert len(panel) == 300
    for ctx in ("A", "B", "C"):
        p = str(ROOT / f"data/challenge/context_{ctx}.h5ad")
        Xlog, _ = log1p_controls(p, genes, None)
        rng = np.random.RandomState(383)
        take = np.sort(rng.choice(len(Xlog), size=min(args.controls, len(Xlog)), replace=False))
        write_h5ad(TDIR / f"ctx{ctx}_controls.h5ad", Xlog[take],
                   {"target_gene": ["non-targeting"] * len(take),
                    "cell_type": [ctx] * len(take)}, genes)
        write_tsv(TDIR / f"ctx{ctx}.tsv", panel, args.cells_per_pert)
        print(ctx, ":", len(take), "controls, 300 perts x", args.cells_per_pert)


def main() -> None:
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(required=True)
    p = sub.add_parser("identity")
    p.add_argument("--controls", type=int, default=1000)
    p.add_argument("--cells-per-pert", type=int, default=100)
    p.set_defaults(f=cmd_identity)
    p = sub.add_parser("h1")
    p.add_argument("--controls", type=int, default=2000)
    p.add_argument("--cells-per-pert", type=int, default=100)
    p.set_defaults(f=cmd_h1)
    p = sub.add_parser("ctx")
    p.add_argument("--controls", type=int, default=2000)
    p.add_argument("--cells-per-pert", type=int, default=100)
    p.set_defaults(f=cmd_ctx)
    args = ap.parse_args()
    args.f(args)


if __name__ == "__main__":
    main()
