#!/usr/bin/env python3
"""Leakage-aware LODO test of the PerturbMatch K562 reanalysis.

PerturbMatch and the incumbent K562 donor derive from the same Replogle
experiment, so K562 can never be a valid held-out truth for this comparison.
Only Orion HCT116 and HEK293T are held out.  Every query is ranked against the
complete target intersection, and all 300 challenge target coordinates are
removed from both scoring and whitening.
"""

from __future__ import annotations

import argparse
import csv
import gc
import hashlib
import json
import os
import resource
import sys
from dataclasses import dataclass
from pathlib import Path

os.environ.setdefault("OPENBLAS_NUM_THREADS", "2")
os.environ.setdefault("OMP_NUM_THREADS", "2")
os.environ.setdefault("MKL_NUM_THREADS", "2")

import numpy as np
import pandas as pd
from scipy.linalg import eigh

ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
FOLDS = ("orion_hct116", "orion_hek293t")
ALL_BANKS = ("k562_gwps", "perturbmatch", *FOLDS)
FRACTIONS = (0.1, 0.25, 0.5, 0.75, 0.9)
GAMMA = 0.65


def peak_rss_gib() -> float:
    value = float(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
    return value / 1024**3 if sys.platform == "darwin" else value / 1024**2


def enforce_limit(limit_gib: float, stage: str) -> None:
    if limit_gib > 0 and peak_rss_gib() >= limit_gib:
        raise MemoryError(
            f"{stage}: peak RSS {peak_rss_gib():.3f} GiB >= {limit_gib:.3f} GiB"
        )


def read_lines(path: Path) -> tuple[str, ...]:
    return tuple(line.strip() for line in path.read_text().splitlines() if line.strip())


def file_sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while block := handle.read(8 << 20):
            digest.update(block)
    return digest.hexdigest()


@dataclass(frozen=True)
class Bank:
    name: str
    path: Path
    genes: tuple[str, ...]
    perts: tuple[str, ...]
    gene_pos: dict[str, int]
    pert_pos: dict[str, int]
    delta_path: Path
    shape: tuple[int, int]
    duplicate_gene_symbols: tuple[str, ...]

    @classmethod
    def load(
        cls, name: str, path: Path, *, require_perturbmatch: bool = False
    ) -> "Bank":
        path = path.resolve()
        genes = read_lines(path / "genes.txt")
        perts = read_lines(path / "perts.txt")
        if len(perts) != len(set(perts)):
            raise ValueError(f"{name}: duplicate perturbations")
        delta_path = path / "delta.npy"
        delta = np.load(delta_path, mmap_mode="r")
        if delta.shape != (len(perts), len(genes)):
            raise ValueError(
                f"{name}: delta shape {delta.shape} != {(len(perts), len(genes))}"
            )
        b_ctrl_path = path / "b_ctrl.npy"
        b_ctrl = (
            np.load(b_ctrl_path, mmap_mode="r")
            if b_ctrl_path.exists()
            else np.zeros(len(genes), dtype=np.float32)
        )
        gene_pos: dict[str, int] = {}
        duplicates: set[str] = set()
        for index, gene in enumerate(genes):
            previous = gene_pos.get(gene)
            if previous is None:
                gene_pos[gene] = index
            else:
                duplicates.add(gene)
                if float(b_ctrl[index]) > float(b_ctrl[previous]):
                    gene_pos[gene] = index
        if require_perturbmatch:
            metadata = json.loads((path / "meta.json").read_text())
            if metadata.get("source_md5") != "0017076a5cc7b9979febf3bd19bd2157":
                raise ValueError(
                    "PerturbMatch source checksum is not the official artifact"
                )
            if metadata.get("vcc_target_coverage") != 251:
                raise ValueError("PerturbMatch VCC target coverage must be 251")
            if duplicates:
                raise ValueError("PerturbMatch compact axis must be unique")
        return cls(
            name=name,
            path=path,
            genes=genes,
            perts=perts,
            gene_pos=gene_pos,
            pert_pos={pert: index for index, pert in enumerate(perts)},
            delta_path=delta_path,
            shape=tuple(map(int, delta.shape)),
            duplicate_gene_symbols=tuple(sorted(duplicates)),
        )


@dataclass
class Prepared:
    source_raw: np.ndarray
    source_gene_norm: np.ndarray
    truth_raw: np.ndarray
    center: np.ndarray
    raw_scale: float
    gene_scale: np.ndarray
    gene_norm_scale: float
    background_rows: np.ndarray
    columns: np.ndarray


def choose_rows(rows: np.ndarray, maximum: int, seed: int) -> np.ndarray:
    if len(rows) <= maximum:
        return rows
    return np.sort(np.random.default_rng(seed).choice(rows, maximum, replace=False))


def prepare_bank(
    bank: Bank,
    targets: tuple[str, ...],
    features: tuple[str, ...],
    panel: set[str],
    *,
    chunk_rows: int,
    max_rss_gib: float,
) -> Prepared:
    columns = np.asarray([bank.gene_pos[gene] for gene in features], dtype=np.int64)
    background = np.asarray(
        [index for index, pert in enumerate(bank.perts) if pert not in panel],
        dtype=np.int64,
    )
    if not len(background):
        raise ValueError(f"{bank.name}: no non-panel background perturbations")
    delta = np.load(bank.delta_path, mmap_mode="r")
    center = np.zeros(len(features), dtype=np.float64)
    for lo in range(0, len(background), chunk_rows):
        rows = background[lo : lo + chunk_rows]
        center += np.asarray(delta[np.ix_(rows, columns)], dtype=np.float64).sum(0)
        enforce_limit(max_rss_gib, f"{bank.name} center")
    center /= len(background)

    absolute = np.zeros(len(features), dtype=np.float64)
    raw_norms = np.empty(len(background), dtype=np.float64)
    for lo in range(0, len(background), chunk_rows):
        rows = background[lo : lo + chunk_rows]
        block = np.asarray(delta[np.ix_(rows, columns)], dtype=np.float64) - center
        absolute += np.abs(block).sum(0)
        raw_norms[lo : lo + len(rows)] = np.linalg.norm(block, axis=1)
        enforce_limit(max_rss_gib, f"{bank.name} raw scale")
    raw_scale = float(np.median(raw_norms)) or 1.0
    mean_abs = absolute / len(background)
    live = mean_abs[mean_abs > 0]
    floor = float(np.percentile(live, 25)) if live.size else 1.0
    gene_scale = np.maximum(mean_abs, floor)

    normalized_norms = np.empty(len(background), dtype=np.float64)
    for lo in range(0, len(background), chunk_rows):
        rows = background[lo : lo + chunk_rows]
        block = np.asarray(delta[np.ix_(rows, columns)], dtype=np.float64) - center
        normalized_norms[lo : lo + len(rows)] = np.linalg.norm(
            block / gene_scale, axis=1
        )
        enforce_limit(max_rss_gib, f"{bank.name} gene-normalized scale")
    gene_norm_scale = float(np.median(normalized_norms)) or 1.0

    target_rows = np.asarray(
        [bank.pert_pos[target] for target in targets], dtype=np.int64
    )
    truth = np.asarray(delta[np.ix_(target_rows, columns)], dtype=np.float64)
    centered = truth - center
    return Prepared(
        source_raw=centered / raw_scale,
        source_gene_norm=(centered / gene_scale) / gene_norm_scale,
        truth_raw=truth,
        center=center,
        raw_scale=raw_scale,
        gene_scale=gene_scale,
        gene_norm_scale=gene_norm_scale,
        background_rows=background,
        columns=columns,
    )


def row_normalize(matrix: np.ndarray) -> np.ndarray:
    matrix = np.asarray(matrix, dtype=np.float64)
    return matrix / np.maximum(np.linalg.norm(matrix, axis=1, keepdims=True), 1e-12)


def pds_per_query(prediction: np.ndarray, truth: np.ndarray) -> np.ndarray:
    similarity = row_normalize(prediction) @ row_normalize(truth).T
    own = np.diag(similarity)[:, None]
    n = len(similarity)
    if similarity.shape != (n, n) or n < 2:
        raise ValueError("PDS requires a square gallery with at least two targets")
    rank = (similarity > own).sum(1).astype(np.float64)
    rank += 0.5 * ((similarity == own).sum(1) - 1)
    return 1.0 - rank / (n - 1)


def candidates(
    prepared: dict[str, Prepared],
    held: str,
    mode: str,
) -> dict[str, np.ndarray]:
    other = next(name for name in FOLDS if name != held)
    old = getattr(prepared["k562_gwps"], f"source_{mode}")
    pm = getattr(prepared["perturbmatch"], f"source_{mode}")
    companion = getattr(prepared[other], f"source_{mode}")
    output = {
        "k562_direct": old,
        "perturbmatch_direct": pm,
        "baseline_old_k562_plus_orion": 0.5 * (old + companion),
        "replacement_pm_plus_orion": 0.5 * (pm + companion),
    }
    for fraction in FRACTIONS:
        slot = (1.0 - fraction) * old + fraction * pm
        label = f"{fraction:.2f}"
        output[f"k562_ensemble_f{label}"] = slot
        output[f"full_ensemble_f{label}"] = 0.5 * (slot + companion)
    return output


def fit_whitener(
    banks: dict[str, Bank],
    prepared: dict[str, Prepared],
    held: str,
    *,
    basis_per_donor: int,
    chunk_rows: int,
    max_rss_gib: float,
) -> tuple[np.ndarray, np.ndarray, dict[str, object]]:
    source_names = ("k562_gwps", next(name for name in FOLDS if name != held))
    selected = {
        name: choose_rows(
            prepared[name].background_rows,
            basis_per_donor,
            int.from_bytes(
                hashlib.sha256(f"058|{held}|{name}".encode()).digest()[:4], "little"
            ),
        )
        for name in source_names
    }
    n_features = prepared[source_names[0]].source_raw.shape[1]
    total = sum(len(rows) for rows in selected.values())
    basis = np.empty((total, n_features), dtype=np.float32)
    cursor = 0
    for name in source_names:
        bank = banks[name]
        item = prepared[name]
        delta = np.load(bank.delta_path, mmap_mode="r")
        rows = selected[name]
        for lo in range(0, len(rows), chunk_rows):
            rr = rows[lo : lo + chunk_rows]
            block = np.asarray(delta[np.ix_(rr, item.columns)], dtype=np.float32)
            block -= item.center.astype(np.float32)
            block /= np.float32(item.raw_scale)
            basis[cursor + lo : cursor + lo + len(rr)] = block
        cursor += len(rows)
        enforce_limit(max_rss_gib, f"{held} whitening basis {name}")
    covariance = (basis.T @ basis) / float(len(basis))
    del basis
    gc.collect()
    evals, evecs = eigh(covariance, overwrite_a=True, check_finite=False, driver="evr")
    del covariance
    gc.collect()
    evals = np.maximum(evals.astype(np.float64), float(evals[-1]) * 1e-12)
    evecs = np.asarray(evecs, dtype=np.float32)
    enforce_limit(max_rss_gib, f"{held} eigendecomposition")
    return (
        evals,
        evecs,
        {
            "sources": list(source_names),
            "basis_rows": {name: len(rows) for name, rows in selected.items()},
            "features": n_features,
            "gamma": GAMMA,
            "operator_excludes_panel_perturbation_rows": True,
        },
    )


def whiten(matrix: np.ndarray, evals: np.ndarray, evecs: np.ndarray) -> np.ndarray:
    values = np.asarray(matrix, dtype=np.float64)
    weights = np.power(evals, -GAMMA / 2.0)
    return ((values @ evecs) * weights) @ evecs.T


def evaluate_fold(
    banks: dict[str, Bank],
    prepared: dict[str, Prepared],
    held: str,
    *,
    basis_per_donor: int,
    chunk_rows: int,
    max_rss_gib: float,
) -> tuple[list[dict[str, object]], dict[str, object], list[dict[str, object]]]:
    truth = prepared[held].truth_raw
    rows: list[dict[str, object]] = []
    per_target_rows: list[dict[str, object]] = []
    operator_info: dict[str, object] = {}
    for mode in ("raw", "gene_norm"):
        proposed = candidates(prepared, held, mode)
        baseline_name = "baseline_old_k562_plus_orion"
        baseline = float(np.mean(pds_per_query(proposed[baseline_name], truth)))
        direct_old = float(np.mean(pds_per_query(proposed["k562_direct"], truth)))
        for name, prediction in proposed.items():
            per_query = pds_per_query(prediction, truth)
            value = float(np.mean(per_query))
            comparator = (
                direct_old
                if "direct" in name or name.startswith("k562_ensemble")
                else baseline
            )
            rows.append(
                {
                    "held": held,
                    "mode": mode,
                    "geometry": "identity",
                    "candidate": name,
                    "pds": value,
                    "comparator_pds": comparator,
                    "delta_vs_comparator": value - comparator,
                }
            )
            if mode == "raw":
                for index, target_score in enumerate(per_query):
                    per_target_rows.append(
                        {
                            "held": held,
                            "mode": mode,
                            "geometry": "identity",
                            "candidate": name,
                            "target_index": index,
                            "pds": float(target_score),
                        }
                    )

    evals, evecs, operator_info = fit_whitener(
        banks,
        prepared,
        held,
        basis_per_donor=basis_per_donor,
        chunk_rows=chunk_rows,
        max_rss_gib=max_rss_gib,
    )
    proposed = candidates(prepared, held, "raw")
    transformed_truth = truth
    transformed = {
        name: whiten(matrix, evals, evecs) for name, matrix in proposed.items()
    }
    baseline_name = "baseline_old_k562_plus_orion"
    baseline = float(
        np.mean(pds_per_query(transformed[baseline_name], transformed_truth))
    )
    direct_old = float(
        np.mean(pds_per_query(transformed["k562_direct"], transformed_truth))
    )
    for name, prediction in transformed.items():
        per_query = pds_per_query(prediction, transformed_truth)
        value = float(np.mean(per_query))
        comparator = (
            direct_old
            if "direct" in name or name.startswith("k562_ensemble")
            else baseline
        )
        rows.append(
            {
                "held": held,
                "mode": "raw",
                "geometry": "source_whitened_g0.65",
                "candidate": name,
                "pds": value,
                "comparator_pds": comparator,
                "delta_vs_comparator": value - comparator,
            }
        )
        for index, target_score in enumerate(per_query):
            per_target_rows.append(
                {
                    "held": held,
                    "mode": "raw",
                    "geometry": "source_whitened_g0.65",
                    "candidate": name,
                    "target_index": index,
                    "pds": float(target_score),
                }
            )
    del evals, evecs, transformed
    gc.collect()
    enforce_limit(max_rss_gib, f"{held} completed")
    return rows, operator_info, per_target_rows


def summarize(rows: list[dict[str, object]]) -> list[dict[str, object]]:
    grouped: dict[tuple[str, str, str], list[dict[str, object]]] = {}
    for row in rows:
        key = (str(row["mode"]), str(row["geometry"]), str(row["candidate"]))
        grouped.setdefault(key, []).append(row)
    output = []
    for (mode, geometry, candidate), values in grouped.items():
        if len(values) != len(FOLDS):
            raise ValueError(f"{mode}/{geometry}/{candidate}: incomplete folds")
        scores = [float(row["pds"]) for row in values]
        deltas = [float(row["delta_vs_comparator"]) for row in values]
        output.append(
            {
                "mode": mode,
                "geometry": geometry,
                "candidate": candidate,
                "macro_pds": float(np.mean(scores)),
                "macro_delta_vs_comparator": float(np.mean(deltas)),
                "min_delta_vs_comparator": float(np.min(deltas)),
                "positive_folds": int(sum(delta > 0 for delta in deltas)),
                "fold_pds": {str(row["held"]): float(row["pds"]) for row in values},
                "fold_delta": {
                    str(row["held"]): float(row["delta_vs_comparator"])
                    for row in values
                },
            }
        )
    return sorted(
        output,
        key=lambda row: (
            str(row["mode"]),
            str(row["geometry"]),
            -float(row["macro_pds"]),
        ),
    )


def write_csv(path: Path, rows: list[dict[str, object]]) -> None:
    if not rows:
        return
    fields = sorted({key for row in rows for key in row})
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def run(args: argparse.Namespace) -> dict[str, object]:
    panel_values = tuple(
        pd.read_csv(args.panel_csv)["target_gene"].dropna().astype(str).str.strip()
    )
    official = tuple(
        pd.read_csv(args.genes_csv)["gene_name"].dropna().astype(str).str.strip()
    )
    if len(panel_values) != 300 or len(set(panel_values)) != 300:
        raise ValueError("challenge panel must contain 300 unique targets")
    panel = set(panel_values)
    paths = {
        "k562_gwps": Path(args.k562),
        "perturbmatch": Path(args.perturbmatch),
        "orion_hct116": Path(args.hct116),
        "orion_hek293t": Path(args.hek293t),
    }
    banks = {
        name: Bank.load(name, path, require_perturbmatch=name == "perturbmatch")
        for name, path in paths.items()
    }
    targets = tuple(
        target
        for target in panel_values
        if all(target in bank.pert_pos for bank in banks.values())
    )
    features = tuple(
        gene
        for gene in official
        if gene not in panel and all(gene in bank.gene_pos for bank in banks.values())
    )
    if len(targets) < args.min_shared_targets:
        raise ValueError(f"only {len(targets)} targets shared across all banks")
    if len(features) < args.min_shared_genes:
        raise ValueError(
            f"only {len(features)} non-panel genes shared across all banks"
        )
    if panel & set(features):
        raise AssertionError("panel target coordinate entered feature axis")

    prepared = {
        name: prepare_bank(
            bank,
            targets,
            features,
            panel,
            chunk_rows=args.chunk_rows,
            max_rss_gib=args.max_rss_gib,
        )
        for name, bank in banks.items()
    }
    all_rows: list[dict[str, object]] = []
    all_per_target: list[dict[str, object]] = []
    operators: dict[str, object] = {}
    for held in FOLDS:
        rows, info, per_target = evaluate_fold(
            banks,
            prepared,
            held,
            basis_per_donor=args.basis_per_donor,
            chunk_rows=args.chunk_rows,
            max_rss_gib=args.max_rss_gib,
        )
        all_rows.extend(rows)
        all_per_target.extend(per_target)
        operators[held] = info

    summary = summarize(all_rows)
    old = prepared["k562_gwps"].source_raw
    pm = prepared["perturbmatch"].source_raw
    agreement = np.sum(row_normalize(old) * row_normalize(pm), axis=1)
    result: dict[str, object] = {
        "experiment": "058-perturbmatch",
        "status": "completed",
        "protocol": {
            "valid_holdouts": list(FOLDS),
            "invalid_holdout_omitted": (
                "K562: PerturbMatch and incumbent K562 use the same Replogle assay"
            ),
            "retrieval_gallery": (
                "complete VCC-panel target intersection across all four banks"
            ),
            "target_coordinate_policy": (
                "all 300 VCC panel genes excluded from score and operator axes"
            ),
            "truth": "uncentered raw held-Orion effect vectors",
            "source_center": (
                "gene mean over source perturbations excluding all 300 panel targets"
            ),
            "whitener": (
                "candidate-independent source-only covariance from incumbent "
                "K562 plus the non-held Orion donor"
            ),
            "gamma": GAMMA,
            "fractions": list(FRACTIONS),
            "submission": False,
            "cell_eval2": False,
        },
        "axis": {
            "targets": len(targets),
            "target_order": list(targets),
            "target_order_sha256": hashlib.sha256(
                ("\n".join(targets) + "\n").encode()
            ).hexdigest(),
            "features": len(features),
            "feature_axis_sha256": hashlib.sha256(
                ("\n".join(features) + "\n").encode()
            ).hexdigest(),
        },
        "banks": {
            name: {
                "path": str(bank.path),
                "shape": list(bank.shape),
                "duplicate_gene_symbols": list(bank.duplicate_gene_symbols),
                "delta_sha256": file_sha256(bank.delta_path),
                "background_rows": len(prepared[name].background_rows),
                "raw_scale": prepared[name].raw_scale,
                "gene_norm_scale": prepared[name].gene_norm_scale,
            }
            for name, bank in banks.items()
        },
        "k562_reanalysis_agreement": {
            "mean_cosine": float(np.mean(agreement)),
            "median_cosine": float(np.median(agreement)),
            "p10_cosine": float(np.percentile(agreement, 10)),
            "p90_cosine": float(np.percentile(agreement, 90)),
        },
        "operators": operators,
        "fold_rows": all_rows,
        "summary": summary,
        "peak_rss_gib": peak_rss_gib(),
    }
    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(result, indent=2) + "\n")
    write_csv(output.with_suffix(".csv"), all_rows)
    write_csv(output.parent / "per_target_pds.csv", all_per_target)
    print(json.dumps(result["axis"], indent=2), flush=True)
    for row in summary:
        if (
            row["candidate"]
            in {
                "k562_direct",
                "perturbmatch_direct",
                "baseline_old_k562_plus_orion",
                "replacement_pm_plus_orion",
                "full_ensemble_f0.10",
                "full_ensemble_f0.25",
                "full_ensemble_f0.50",
            }
            and row["mode"] == "raw"
        ):
            print(json.dumps(row), flush=True)
    print(f"peak_rss_gib={peak_rss_gib():.3f}", flush=True)
    return result


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--k562", default=str(ROOT / "data/donors/k562_gwps"))
    parser.add_argument("--hct116", default=str(ROOT / "data/donors/orion_hct116"))
    parser.add_argument("--hek293t", default=str(ROOT / "data/donors/orion_hek293t"))
    parser.add_argument(
        "--perturbmatch", default=str(ROOT / "data/external/perturbmatch")
    )
    parser.add_argument(
        "--panel-csv", default=str(ROOT / "data/challenge/pert_counts.csv")
    )
    parser.add_argument(
        "--genes-csv", default=str(ROOT / "data/challenge/gene_names.csv")
    )
    parser.add_argument("--basis-per-donor", type=int, default=6000)
    parser.add_argument("--chunk-rows", type=int, default=128)
    parser.add_argument("--max-rss-gib", type=float, default=12.0)
    parser.add_argument("--min-shared-targets", type=int, default=220)
    parser.add_argument("--min-shared-genes", type=int, default=2400)
    parser.add_argument("--output", default=str(HERE / "results.json"))
    return parser


if __name__ == "__main__":
    run(build_parser().parse_args())
