"""Can non-panel CRISPRi screens improve the whitening geometry?

Prediction still uses only the production donors. Auxiliary screens contribute
only rows to the signal-covariance estimate, so target-panel overlap is not
required. This first-stage screen uses the donor-intersection operator axis for
speed; any winner must be repeated on the production union axis.
"""
from __future__ import annotations

import argparse
import gc
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))
from vcc2026.whiten import Whitener  # noqa: E402

MAIN = {
    "k562_gwps": ROOT / "data/donors/k562_gwps",
    "orion_hct116": ROOT / "data/donors/orion_hct116",
    "orion_hek293t": ROOT / "data/donors/orion_hek293t",
}
AUX = {
    "h1_2025": ROOT / "data/donors/h1_2025",
    "k562_essential": ROOT / "data/donors/k562_essential",
    "rpe1": ROOT / "data/donors/rpe1",
    "nadig_jurkat": ROOT / "data/donors/nadig_jurkat",
}
PANEL = set(pd.read_csv(ROOT / "data/challenge/pert_counts.csv")["target_gene"].astype(str))
OFFICIAL = set(pd.read_csv(ROOT / "data/challenge/gene_names.csv")["gene_name"].astype(str))


def load_bank(path: Path) -> dict:
    genes = np.asarray((path / "genes.txt").read_text().splitlines(), dtype=str)
    perts = np.asarray((path / "perts.txt").read_text().splitlines(), dtype=str)
    delta = np.load(path / "delta.npy", mmap_mode="r")
    if delta.shape != (len(perts), len(genes)):
        raise ValueError(
            f"{path}: delta shape {delta.shape} != ({len(perts)}, {len(genes)})"
        )
    return {
        "delta_path": path / "delta.npy",
        "n_rows": int(delta.shape[0]),
        "genes": genes,
        "pos": {g: i for i, g in enumerate(genes)},
        "gset": set(genes),
        "idx": {str(p): i for i, p in enumerate(perts)},
    }


def pds(pred: np.ndarray, truth: np.ndarray) -> float:
    p = pred / (np.linalg.norm(pred, axis=1, keepdims=True) + 1e-12)
    t = truth / (np.linalg.norm(truth, axis=1, keepdims=True) + 1e-12)
    dist = 1.0 - p @ t.T
    own = np.diag(dist)[:, None]
    n = dist.shape[1]
    rank = (dist < own).sum(1) + 0.5 * ((dist == own).sum(1) - 1)
    return float(np.mean(1.0 - rank / (n - 1)))


def column_mean(bank: dict, src: np.ndarray, chunk_rows: int = 512) -> np.ndarray:
    delta = np.load(bank["delta_path"], mmap_mode="r")
    total = np.zeros(len(src), dtype=np.float64)
    for lo in range(0, bank["n_rows"], chunk_rows):
        hi = min(lo + chunk_rows, bank["n_rows"])
        total += np.asarray(delta[lo:hi, src], dtype=np.float64).sum(axis=0)
    return total / bank["n_rows"]


def centered_rows(bank: dict, rows: np.ndarray, src: np.ndarray) -> np.ndarray:
    delta = np.load(bank["delta_path"], mmap_mode="r")
    mean = column_mean(bank, src)
    return np.asarray(delta[np.ix_(rows, src)], dtype=np.float64) - mean


def aligned_rows(bank: dict, feat: list[str], ref_norm: float) -> np.ndarray:
    out = np.zeros((bank["n_rows"], len(feat)), dtype=np.float32)
    dst, src = [], []
    for j, gene in enumerate(feat):
        k = bank["pos"].get(gene)
        if k is not None:
            dst.append(j)
            src.append(k)
    if src:
        src_arr = np.asarray(src, dtype=np.int64)
        mean = column_mean(bank, src_arr)
        delta = np.load(bank["delta_path"], mmap_mode="r")
        for lo in range(0, bank["n_rows"], 512):
            hi = min(lo + 512, bank["n_rows"])
            out[lo:hi, np.asarray(dst)] = (
                np.asarray(delta[lo:hi, src_arr], dtype=np.float32)
                - mean.astype(np.float32)
            )
    med = float(np.median(np.linalg.norm(out, axis=1))) or 1.0
    return out.astype(np.float64) * (ref_norm / med)


def combine_basis(parts: dict[str, np.ndarray], names: list[str], equal_source: bool) -> np.ndarray:
    rows = []
    if equal_source:
        target_n = float(np.mean([len(parts[n]) for n in names]))
        for n in names:
            rows.append(parts[n] * np.sqrt(target_n / len(parts[n])))
    else:
        rows = [parts[n] for n in names]
    return np.vstack(rows)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--axis", choices=["inter", "union"], default="inter")
    ap.add_argument("--gammas", nargs="+", type=float, default=[0.45, 0.55, 0.65, 0.75, 0.85])
    ap.add_argument(
        "--max-genes",
        type=int,
        default=0,
        help="deterministically subsample the operator axis for a low-memory screen",
    )
    ap.add_argument("--gene-seed", type=int, default=20260822)
    ap.add_argument(
        "--held",
        nargs="+",
        choices=["orion_hct116", "orion_hek293t"],
        default=["orion_hct116", "orion_hek293t"],
        help="folds to run; use one per process to release LAPACK workspace",
    )
    ap.add_argument(
        "--variants",
        nargs="+",
        default=None,
        help="optional subset of variant names to fit",
    )
    args = ap.parse_args()

    main = {name: load_bank(path) for name, path in MAIN.items()}
    aux = {name: load_bank(path) for name, path in AUX.items() if (path / "delta.npy").exists()}
    common = [
        p for p in pd.read_csv(ROOT / "data/challenge/pert_counts.csv")["target_gene"].astype(str)
        if all(p in bank["idx"] for bank in main.values())
    ]
    variants = {
        "base_rows": ([], False),
        "base_equal_source": ([], True),
        "plus_h1_rows": (["h1_2025"], False),
        "plus_essential_rows": (["k562_essential", "rpe1"], False),
        "plus_jurkat_rows": (["nadig_jurkat"], False),
        "plus_all_no_h1_rows": (
            ["k562_essential", "rpe1", "nadig_jurkat"],
            False,
        ),
        "plus_all_rows": (
            ["h1_2025", "k562_essential", "rpe1", "nadig_jurkat"],
            False,
        ),
        "all_equal_source": (
            ["h1_2025", "k562_essential", "rpe1", "nadig_jurkat"],
            True,
        ),
    }
    if args.variants:
        unknown = sorted(set(args.variants) - set(variants))
        if unknown:
            raise SystemExit(f"unknown variants: {unknown}")
        variants = {name: variants[name] for name in args.variants}
    results: dict[str, dict[str, float]] = {}

    for held in args.held:
        others = [n for n in main if n != held]
        feat = sorted((main[held]["gset"] & OFFICIAL) - PANEL)
        truth_cols = np.array([main[held]["pos"][g] for g in feat])
        truth_rows = np.asarray([main[held]["idx"][p] for p in common])
        truth = centered_rows(main[held], truth_rows, truth_cols)

        pred_parts = {}
        donor_norms = {}
        for name in others:
            cols = np.array([main[name]["pos"].get(g, -1) for g in feat])
            good = cols >= 0
            arr = np.zeros((len(common), len(feat)), dtype=np.float64)
            rows = np.asarray([main[name]["idx"][p] for p in common])
            arr[:, good] = centered_rows(main[name], rows, cols[good])
            donor_norms[name] = float(np.median(np.linalg.norm(arr, axis=1))) or 1.0
            pred_parts[name] = arr
        ref_norm = float(np.median(list(donor_norms.values())))
        pred = sum(pred_parts[n] * (ref_norm / donor_norms[n]) for n in others) / len(others)

        if args.axis == "inter":
            op_genes = sorted(set.intersection(*[main[n]["gset"] for n in others]) & set(feat))
        else:
            op_genes = sorted(set.union(*[main[n]["gset"] for n in others]) & set(feat))
        if args.max_genes and len(op_genes) > args.max_genes:
            rng = np.random.default_rng(args.gene_seed)
            keep = np.sort(
                rng.choice(len(op_genes), size=args.max_genes, replace=False)
            )
            op_genes = [op_genes[i] for i in keep]
        fpos = {g: i for i, g in enumerate(feat)}
        wcols = np.array([fpos[g] for g in op_genes])
        parts = {
            n: aligned_rows(main[n], op_genes, ref_norm) for n in others
        }
        parts.update({
            n: aligned_rows(bank, op_genes, ref_norm) for n, bank in aux.items()
        })
        print(
            f"fold {held}: panel={len(common)} score_genes={len(feat)} "
            f"operator_genes={len(op_genes)} aux={list(aux)}",
            flush=True,
        )

        for variant, (extra, equal_source) in variants.items():
            names = others + [n for n in extra if n in parts]
            basis = combine_basis(parts, names, equal_source)
            print(f"  fitting {variant}: {basis.shape}", flush=True)
            whitener = Whitener.fit(basis)
            for gamma in args.gammas:
                transformed = pred.copy()
                before = float(np.median(np.linalg.norm(pred[:, wcols], axis=1)))
                transformed[:, wcols] = whitener.apply(
                    pred[:, wcols], gamma, renorm=False
                )
                after = float(np.median(np.linalg.norm(transformed[:, wcols], axis=1))) or 1.0
                transformed[:, wcols] *= before / after
                key = f"{variant}:g={gamma:g}"
                results.setdefault(key, {})[held] = pds(transformed, truth)
                print(
                    f"    gamma={gamma:g} pds={results[key][held]:.4f}",
                    flush=True,
                )
            del whitener, basis
            gc.collect()

    out = {}
    for key, fold_values in results.items():
        out[key] = {**fold_values, "macro": float(np.mean(list(fold_values.values())))}
    ordered = sorted(out.items(), key=lambda kv: -kv[1]["macro"])
    for key, value in ordered:
        print(f"{key:34} macro={value['macro']:.4f} {value}", flush=True)
    suffix = args.axis
    if args.max_genes:
        suffix += f"_g{args.max_genes}_s{args.gene_seed}"
    if args.variants:
        suffix += "_" + "-".join(args.variants)
    suffix += "_" + "-".join(args.held)
    (Path(__file__).parent / f"results_{suffix}.json").write_text(
        json.dumps(out, indent=2)
    )


if __name__ == "__main__":
    main()
