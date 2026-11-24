#!/usr/bin/env python3
"""057 — fixed KOLF fraction under source-only union whitening.

The experiment keeps experiment 055's chosen rule frozen: raw KOLF effects,
mixed at fraction 0.1.  HCT116 and HEK293T are the only predictive sources and
the only rows used to fit the whitening operator.  K562 is opened as truth only
after the operator and all predictions have been frozen.

Large arrays are memory-mapped and every generated artifact stays below the
experiment's output directory.  Real-data execution is intended for AWS.
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
from typing import Iterable, Sequence

os.environ.setdefault("OPENBLAS_NUM_THREADS", "4")
os.environ.setdefault("OMP_NUM_THREADS", "4")
os.environ.setdefault("MKL_NUM_THREADS", "4")

import numpy as np
import pandas as pd
from scipy.linalg import eigh


ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
SOURCE_NAMES = ("orion_hct116", "orion_hek293t")
HOLDOUT_NAME = "k562_gwps"
AUXILIARY_NAME = "kolf2"
KOLF_FRACTION = 0.1
EFFECT_SCALE = 6.817
DEFAULT_GAMMAS = (0.0, 0.55, 0.65, 0.75, 0.85)
PRIMARY_GAMMA = 0.65


def rss_gib() -> float:
    rss = float(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
    return rss / 1024**3 if sys.platform == "darwin" else rss / 1024**2


def enforce_rss(limit_gib: float, stage: str) -> None:
    if limit_gib > 0 and rss_gib() >= limit_gib:
        raise MemoryError(
            f"{stage}: peak RSS {rss_gib():.3f} GiB reached {limit_gib:.3f} GiB"
        )


def read_lines(path: Path, *, unique: bool = True) -> tuple[str, ...]:
    values = tuple(
        line.strip() for line in path.read_text().splitlines() if line.strip()
    )
    if unique and len(values) != len(set(values)):
        raise ValueError(f"{path}: duplicate labels")
    return values


def read_csv_column(path: Path, column: str) -> tuple[str, ...]:
    frame = pd.read_csv(path)
    if column not in frame:
        raise ValueError(f"{path}: missing column {column!r}")
    values = tuple(
        value
        for value in frame[column].dropna().astype(str).str.strip().tolist()
        if value
    )
    if len(values) != len(set(values)):
        raise ValueError(f"{path}: duplicate values in {column!r}")
    return values


def sha256_lines(values: Iterable[str]) -> str:
    return hashlib.sha256(("\n".join(values) + "\n").encode()).hexdigest()


def sha256_file(path: Path, block_size: int = 8 << 20) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while block := handle.read(block_size):
            digest.update(block)
    return digest.hexdigest()


def close_memmap(array: np.ndarray) -> None:
    mmap = getattr(array, "_mmap", None)
    if mmap is not None:
        mmap.close()


@dataclass(frozen=True)
class Bank:
    name: str
    path: Path
    genes: tuple[str, ...]
    perts: tuple[str, ...]
    pert_pos: dict[str, int]
    delta_path: Path
    b_ctrl_path: Path
    shape: tuple[int, int]

    @classmethod
    def load(cls, name: str, path: Path, *, require_kolf: bool = False) -> "Bank":
        path = path.expanduser().resolve()
        genes = read_lines(path / "genes.txt", unique=False)
        perts = read_lines(path / "perts.txt")
        delta_path = path / "delta.npy"
        b_ctrl_path = path / "b_ctrl.npy"
        delta = np.load(delta_path, mmap_mode="r")
        b_ctrl = np.load(b_ctrl_path, mmap_mode="r")
        if delta.shape != (len(perts), len(genes)):
            raise ValueError(
                f"{name}: delta shape {delta.shape} != {(len(perts), len(genes))}"
            )
        if b_ctrl.shape != (len(genes),):
            raise ValueError(f"{name}: b_ctrl shape {b_ctrl.shape} != {(len(genes),)}")
        if not np.issubdtype(delta.dtype, np.floating):
            raise TypeError(f"{name}: delta dtype must be floating point")
        close_memmap(delta)
        close_memmap(b_ctrl)
        if require_kolf:
            metadata = json.loads((path / "meta.json").read_text())
            if "int8 source value * 0.04" not in str(
                metadata.get("effect_definition", "")
            ):
                raise ValueError(f"{name}: unexpected KOLF effect definition")
            if int(metadata.get("vcc_target_coverage", -1)) != len(perts):
                raise ValueError(f"{name}: KOLF target coverage metadata mismatch")
        return cls(
            name=name,
            path=path,
            genes=genes,
            perts=perts,
            pert_pos={pert: index for index, pert in enumerate(perts)},
            delta_path=delta_path,
            b_ctrl_path=b_ctrl_path,
            shape=(len(perts), len(genes)),
        )

    def fingerprint(self) -> dict[str, object]:
        delta = np.load(self.delta_path, mmap_mode="r")
        fingerprint = {
            "path": str(self.path),
            "shape": list(map(int, delta.shape)),
            "dtype": str(delta.dtype),
            "delta_bytes": self.delta_path.stat().st_size,
            "genes_sha256": sha256_lines(self.genes),
            "perts_sha256": sha256_lines(self.perts),
        }
        close_memmap(delta)
        return fingerprint


@dataclass(frozen=True)
class GeneMap:
    destination: np.ndarray
    source: np.ndarray
    duplicate_symbols: tuple[str, ...]


def map_to_axis(
    bank: Bank,
    axis: Sequence[str],
    *,
    duplicate_policy: str,
) -> GeneMap:
    """Map one source column to each axis symbol.

    ``max_b_ctrl`` reproduces experiment 055. ``last`` reproduces the
    MultiDonorModel assignment used by production whitening.
    """
    if duplicate_policy not in {"max_b_ctrl", "last"}:
        raise ValueError(f"unknown duplicate policy {duplicate_policy!r}")
    axis_pos = {gene: index for index, gene in enumerate(axis)}
    selected: dict[int, int] = {}
    duplicates: set[str] = set()
    control = np.load(bank.b_ctrl_path, mmap_mode="r")
    try:
        for source, gene in enumerate(bank.genes):
            destination = axis_pos.get(gene)
            if destination is None:
                continue
            previous = selected.get(destination)
            if previous is None:
                selected[destination] = source
            else:
                duplicates.add(gene)
                if duplicate_policy == "last" or float(control[source]) > float(
                    control[previous]
                ):
                    selected[destination] = source
    finally:
        close_memmap(control)
    destination = np.asarray(sorted(selected), dtype=np.int64)
    source = np.asarray([selected[index] for index in destination], dtype=np.int64)
    return GeneMap(
        destination=destination,
        source=source,
        duplicate_symbols=tuple(sorted(duplicates)),
    )


def mapped_source_columns_with_duplicates(
    bank: Bank, official_set: set[str]
) -> np.ndarray:
    return np.fromiter(
        (index for index, gene in enumerate(bank.genes) if gene in official_set),
        dtype=np.int64,
    )


def stream_mean(
    bank: Bank,
    columns: np.ndarray,
    *,
    chunk_rows: int,
    max_rss_gib: float,
) -> np.ndarray:
    delta = np.load(bank.delta_path, mmap_mode="r")
    total = np.zeros(len(columns), dtype=np.float64)
    try:
        for lo in range(0, bank.shape[0], chunk_rows):
            hi = min(bank.shape[0], lo + chunk_rows)
            total += np.asarray(delta[lo:hi, columns], dtype=np.float64).sum(axis=0)
            enforce_rss(max_rss_gib, f"{bank.name} mean")
    finally:
        close_memmap(delta)
    return total / bank.shape[0]


def stream_median_centered_norm(
    bank: Bank,
    columns: np.ndarray,
    mean: np.ndarray,
    *,
    chunk_rows: int,
    max_rss_gib: float,
) -> float:
    delta = np.load(bank.delta_path, mmap_mode="r")
    norms = np.empty(bank.shape[0], dtype=np.float64)
    try:
        for lo in range(0, bank.shape[0], chunk_rows):
            hi = min(bank.shape[0], lo + chunk_rows)
            block = np.asarray(delta[lo:hi, columns], dtype=np.float64)
            block -= mean
            norms[lo:hi] = np.linalg.norm(block, axis=1)
            enforce_rss(max_rss_gib, f"{bank.name} norm")
    finally:
        close_memmap(delta)
    return float(np.median(norms)) or 1.0


@dataclass
class ProductionStats:
    bank: Bank
    mapping: GeneMap
    mean: np.ndarray
    donor_scale: float
    all_mapped_columns: int


def production_stats(
    bank: Bank,
    official: Sequence[str],
    *,
    chunk_rows: int,
    max_rss_gib: float,
) -> ProductionStats:
    """Compute the exact source-wide scale used by MultiDonorModel."""
    official_set = set(official)
    all_columns = mapped_source_columns_with_duplicates(bank, official_set)
    all_mean = stream_mean(
        bank,
        all_columns,
        chunk_rows=chunk_rows,
        max_rss_gib=max_rss_gib,
    )
    donor_scale = stream_median_centered_norm(
        bank,
        all_columns,
        all_mean,
        chunk_rows=chunk_rows,
        max_rss_gib=max_rss_gib,
    )
    mapping = map_to_axis(bank, official, duplicate_policy="last")
    source_to_mean = {
        int(source): float(mean)
        for source, mean in zip(all_columns.tolist(), all_mean.tolist(), strict=True)
    }
    unique_mean = np.asarray(
        [source_to_mean[int(source)] for source in mapping.source], dtype=np.float32
    )
    return ProductionStats(
        bank=bank,
        mapping=mapping,
        mean=unique_mean,
        donor_scale=donor_scale,
        all_mapped_columns=len(all_columns),
    )


def choose_basis_rows(bank: Bank, maximum: int, rng: np.random.Generator) -> np.ndarray:
    rows = np.arange(bank.shape[0], dtype=np.int64)
    if len(rows) > maximum:
        rows = np.sort(rng.choice(rows, maximum, replace=False).astype(np.int64))
    return rows


def build_source_only_operator(
    source_stats: Sequence[ProductionStats],
    official: Sequence[str],
    panel_set: set[str],
    output_dir: Path,
    *,
    basis_per_donor: int,
    operator_seed: int,
    chunk_rows: int,
    max_rss_gib: float,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, dict[str, object]]:
    """Fit the production-style union covariance using source donors only."""
    if tuple(stats.bank.name for stats in source_stats) != SOURCE_NAMES:
        raise ValueError("operator sources or source order changed")
    source_union = set().union(*(set(stats.bank.genes) for stats in source_stats))
    wcols = np.asarray(
        [
            index
            for index, gene in enumerate(official)
            if gene in source_union and gene not in panel_set
        ],
        dtype=np.int64,
    )
    if any(official[index] in panel_set for index in wcols):
        raise AssertionError("VCC target coordinate entered whitening axis")
    local = {int(official_pos): index for index, official_pos in enumerate(wcols)}

    rng = np.random.default_rng(operator_seed)
    selections = {
        stats.bank.name: choose_basis_rows(stats.bank, basis_per_donor, rng)
        for stats in source_stats
    }
    total_rows = sum(len(rows) for rows in selections.values())
    scratch = output_dir / "scratch"
    scratch.mkdir(parents=True, exist_ok=True)
    basis_path = scratch / "source_basis.npy"
    covariance_path = scratch / "union_covariance.npy"
    basis = np.lib.format.open_memmap(
        basis_path,
        mode="w+",
        dtype=np.float32,
        shape=(total_rows, len(wcols)),
    )
    basis[:] = 0.0
    reference = float(np.median([stats.donor_scale for stats in source_stats]))
    selection_audit: dict[str, object] = {}
    cursor = 0
    for stats in source_stats:
        rows = selections[stats.bank.name]
        selected_pairs = [
            (local[int(destination)], int(source), float(center))
            for destination, source, center in zip(
                stats.mapping.destination,
                stats.mapping.source,
                stats.mean,
                strict=True,
            )
            if int(destination) in local
        ]
        basis_dest = np.asarray([item[0] for item in selected_pairs], dtype=np.int64)
        source_cols = np.asarray([item[1] for item in selected_pairs], dtype=np.int64)
        center = np.asarray([item[2] for item in selected_pairs], dtype=np.float32)
        factor = reference / stats.donor_scale
        delta = np.load(stats.bank.delta_path, mmap_mode="r")
        try:
            for lo in range(0, len(rows), chunk_rows):
                hi = min(len(rows), lo + chunk_rows)
                block = np.asarray(
                    delta[np.ix_(rows[lo:hi], source_cols)], dtype=np.float32
                )
                block -= center
                view = basis[cursor + lo : cursor + hi]
                view[:, basis_dest] = block * factor
                enforce_rss(max_rss_gib, f"{stats.bank.name} basis")
        finally:
            close_memmap(delta)
        selection_audit[stats.bank.name] = {
            "rows": len(rows),
            "row_indices_sha256": hashlib.sha256(rows.tobytes()).hexdigest(),
            "mapped_operator_genes": len(source_cols),
            "donor_scale": stats.donor_scale,
            "basis_relative_scale": factor,
        }
        cursor += len(rows)
    basis.flush()

    covariance = np.lib.format.open_memmap(
        covariance_path,
        mode="w+",
        dtype=np.float32,
        shape=(len(wcols), len(wcols)),
        fortran_order=True,
    )
    print(
        f"covariance: basis={basis.shape}, operator={covariance.shape}, "
        f"peak_rss={rss_gib():.2f} GiB",
        flush=True,
    )
    with np.errstate(all="ignore"):
        np.matmul(basis.T, basis, out=covariance)
    covariance /= float(len(basis))
    covariance.flush()
    close_memmap(basis)
    del basis
    basis_path.unlink(missing_ok=True)
    gc.collect()
    enforce_rss(max_rss_gib, "covariance")

    print("eigendecomposition: float32 scipy.linalg.eigh(driver='evr')", flush=True)
    evals, evecs = eigh(
        covariance,
        overwrite_a=True,
        check_finite=False,
        driver="evr",
    )
    close_memmap(covariance)
    del covariance
    covariance_path.unlink(missing_ok=True)
    gc.collect()
    evals = np.maximum(evals.astype(np.float64), float(evals[-1]) * 1e-12)
    evecs = np.asarray(evecs, dtype=np.float32)
    if not np.isfinite(evals).all() or not np.isfinite(evecs[::257, ::257]).all():
        raise FloatingPointError("non-finite whitening operator")
    enforce_rss(max_rss_gib, "eigendecomposition")

    operator_path = output_dir / "whitener_union.npz"
    np.savez(operator_path, evals=evals, evecs=evecs)
    operator_genes = tuple(official[index] for index in wcols)
    audit = {
        "axis": "union",
        "source_donors_only": list(SOURCE_NAMES),
        "forbidden_fit_donors": [HOLDOUT_NAME, AUXILIARY_NAME],
        "basis_per_donor_cap": basis_per_donor,
        "basis_rows": total_rows,
        "operator_seed": operator_seed,
        "operator_genes": len(wcols),
        "operator_genes_sha256": sha256_lines(operator_genes),
        "operator_wcols_sha256": hashlib.sha256(wcols.tobytes()).hexdigest(),
        "excluded_target_coordinates": len(panel_set),
        "excluded_targets_sha256": sha256_lines(sorted(panel_set)),
        "panel_overlap_with_operator": len(panel_set & set(operator_genes)),
        "source_selections": selection_audit,
        "operator_npz": str(operator_path),
        "operator_npz_sha256": sha256_file(operator_path),
    }
    (output_dir / "whitener_union.json").write_text(json.dumps(audit, indent=2) + "\n")
    return wcols, evals, evecs, audit


def centered_panel_on_axis(
    bank: Bank,
    targets: Sequence[str],
    axis: Sequence[str],
    *,
    duplicate_policy: str,
    chunk_rows: int,
    max_rss_gib: float,
) -> tuple[np.ndarray, GeneMap]:
    mapping = map_to_axis(bank, axis, duplicate_policy=duplicate_policy)
    mean = stream_mean(
        bank,
        mapping.source,
        chunk_rows=chunk_rows,
        max_rss_gib=max_rss_gib,
    )
    rows = np.asarray([bank.pert_pos[target] for target in targets], dtype=np.int64)
    output = np.zeros((len(targets), len(axis)), dtype=np.float32)
    delta = np.load(bank.delta_path, mmap_mode="r")
    try:
        for lo in range(0, len(rows), chunk_rows):
            hi = min(len(rows), lo + chunk_rows)
            block = np.asarray(
                delta[np.ix_(rows[lo:hi], mapping.source)], dtype=np.float32
            )
            block -= mean.astype(np.float32)
            output[lo:hi, mapping.destination] = block
            enforce_rss(max_rss_gib, f"{bank.name} panel extraction")
    finally:
        close_memmap(delta)
    return output, mapping


def build_fixed_predictions(
    banks: dict[str, Bank],
    targets: Sequence[str],
    operator_genes: Sequence[str],
    shared_score_genes: Sequence[str],
    *,
    chunk_rows: int,
    max_rss_gib: float,
) -> tuple[np.ndarray, np.ndarray, dict[str, object]]:
    """Embed the exact exp055 raw-0.1 rule into the union operator axis."""
    matrices: dict[str, np.ndarray] = {}
    mappings: dict[str, GeneMap] = {}
    for name in (*SOURCE_NAMES, AUXILIARY_NAME):
        matrices[name], mappings[name] = centered_panel_on_axis(
            banks[name],
            targets,
            operator_genes,
            duplicate_policy="max_b_ctrl",
            chunk_rows=chunk_rows,
            max_rss_gib=max_rss_gib,
        )

    operator_pos = {gene: index for index, gene in enumerate(operator_genes)}
    shared_local = np.asarray(
        [operator_pos[gene] for gene in shared_score_genes], dtype=np.int64
    )
    panel_norms = {
        name: float(np.median(np.linalg.norm(matrix[:, shared_local], axis=1))) or 1.0
        for name, matrix in matrices.items()
    }
    reference = float(np.median([panel_norms[name] for name in SOURCE_NAMES]))
    scaled_sources = [
        matrices[name] * (reference / panel_norms[name]) for name in SOURCE_NAMES
    ]
    baseline = sum(scaled_sources) / len(scaled_sources)
    kolf = matrices[AUXILIARY_NAME] * (reference / panel_norms[AUXILIARY_NAME])
    available = np.zeros(len(operator_genes), dtype=bool)
    available[mappings[AUXILIARY_NAME].destination] = True
    candidate = baseline.copy()
    candidate[:, available] = (1.0 - KOLF_FRACTION) * baseline[
        :, available
    ] + KOLF_FRACTION * kolf[:, available]
    return (
        baseline,
        candidate,
        {
            "rule": "experiment 055 raw KOLF fraction 0.1",
            "source_panel_norms_on_exp055_shared_axis": {
                name: panel_norms[name] for name in SOURCE_NAMES
            },
            "kolf_panel_norm_on_exp055_shared_axis": panel_norms[AUXILIARY_NAME],
            "reference_panel_norm": reference,
            "kolf_fraction": KOLF_FRACTION,
            "kolf_available_operator_genes": int(available.sum()),
            "kolf_missing_policy": "preserve baseline",
            "duplicate_policy": {
                name: "max_b_ctrl" for name in (*SOURCE_NAMES, AUXILIARY_NAME)
            },
            "duplicate_symbols": {
                name: list(mappings[name].duplicate_symbols)
                for name in (*SOURCE_NAMES, AUXILIARY_NAME)
            },
        },
    )


def pds_per_query(prediction: np.ndarray, truth: np.ndarray) -> np.ndarray:
    if prediction.shape != truth.shape:
        raise ValueError(f"PDS shape mismatch {prediction.shape} vs {truth.shape}")
    prediction = np.asarray(prediction, dtype=np.float64).copy()
    truth = np.asarray(truth, dtype=np.float64).copy()
    prediction /= np.maximum(np.linalg.norm(prediction, axis=1, keepdims=True), 1e-12)
    truth /= np.maximum(np.linalg.norm(truth, axis=1, keepdims=True), 1e-12)
    with np.errstate(all="ignore"):
        similarities = prediction @ truth.T
    if not np.isfinite(similarities).all():
        raise FloatingPointError("non-finite cosine similarities")
    own = np.diag(similarities)[:, None]
    n = len(prediction)
    if n < 2:
        raise ValueError("PDS requires at least two targets")
    ranks = (similarities > own).sum(axis=1).astype(np.float64)
    ranks += 0.5 * ((similarities == own).sum(axis=1) - 1)
    return 1.0 - ranks / (n - 1)


def query_fold_summary(
    values: np.ndarray, targets: Sequence[str], folds: int
) -> list[dict[str, object]]:
    output = []
    for fold in range(folds):
        rows = np.arange(fold, len(targets), folds)
        output.append(
            {
                "fold": fold,
                "n_queries": len(rows),
                "pds": float(np.mean(values[rows])),
                "targets_sha256": sha256_lines(targets[index] for index in rows),
                "gallery_size": len(targets),
            }
        )
    return output


def score_gamma(
    baseline: np.ndarray,
    candidate: np.ndarray,
    truth_full: np.ndarray,
    truth_shared: np.ndarray,
    operator_score_positions: np.ndarray,
    evals: np.ndarray,
    evecs: np.ndarray,
    gamma: float,
    targets: Sequence[str],
    folds: int,
) -> tuple[dict[str, object], list[dict[str, object]]]:
    """Transform only observed K562 coordinates; PDS is scale invariant."""
    if gamma == 0:
        transformed = np.concatenate(
            [
                baseline[:, operator_score_positions],
                candidate[:, operator_score_positions],
            ],
            axis=0,
        )
    else:
        combined = np.concatenate([baseline, candidate], axis=0)
        with np.errstate(all="ignore"):
            latent = combined @ evecs
            weights = np.power(evals, -gamma / 2.0).astype(np.float32)
            score_vectors = np.asarray(
                evecs[operator_score_positions, :], dtype=np.float32
            )
            transformed = (latent * weights) @ score_vectors.T
        if not np.isfinite(transformed).all():
            raise FloatingPointError(f"gamma {gamma}: non-finite transformed values")

    n = len(targets)
    baseline_shared = transformed[:n]
    candidate_shared = transformed[n:]
    # The 136 K562 genes outside the HCT/HEK union are exactly zero in both
    # predictions.  They remain in truth_full's norm, matching target-observed
    # scoring rather than silently narrowing the target gene axis.
    baseline_full = np.zeros_like(truth_full)
    candidate_full = np.zeros_like(truth_full)
    baseline_full[:, : baseline_shared.shape[1]] = baseline_shared
    candidate_full[:, : candidate_shared.shape[1]] = candidate_shared

    base_per = pds_per_query(baseline_full, truth_full)
    cand_per = pds_per_query(candidate_full, truth_full)
    base_shared_per = pds_per_query(baseline_shared, truth_shared)
    cand_shared_per = pds_per_query(candidate_shared, truth_shared)
    baseline_folds = query_fold_summary(base_per, targets, folds)
    candidate_folds = query_fold_summary(cand_per, targets, folds)
    rows = []
    for index, target in enumerate(targets):
        rows.append(
            {
                "gamma": gamma,
                "target_gene": target,
                "query_fold": index % folds,
                "baseline_pds": float(base_per[index]),
                "kolf_pds": float(cand_per[index]),
                "delta": float(cand_per[index] - base_per[index]),
                "baseline_shared_axis_pds": float(base_shared_per[index]),
                "kolf_shared_axis_pds": float(cand_shared_per[index]),
                "shared_axis_delta": float(
                    cand_shared_per[index] - base_shared_per[index]
                ),
            }
        )
    return (
        {
            "gamma": gamma,
            "is_primary": gamma == PRIMARY_GAMMA,
            "target_observed_axis": {
                "baseline_pds": float(np.mean(base_per)),
                "kolf_pds": float(np.mean(cand_per)),
                "delta": float(np.mean(cand_per - base_per)),
                "baseline_query_folds": baseline_folds,
                "kolf_query_folds": candidate_folds,
                "fold_deltas": [
                    {
                        "fold": left["fold"],
                        "n_queries": left["n_queries"],
                        "delta": float(right["pds"] - left["pds"]),
                    }
                    for left, right in zip(baseline_folds, candidate_folds, strict=True)
                ],
            },
            "exp055_shared_axis": {
                "baseline_pds": float(np.mean(base_shared_per)),
                "kolf_pds": float(np.mean(cand_shared_per)),
                "delta": float(np.mean(cand_shared_per - base_shared_per)),
            },
        },
        rows,
    )


def write_csv(path: Path, rows: Sequence[dict[str, object]]) -> None:
    if not rows:
        return
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=sorted({key for row in rows for key in row})
        )
        writer.writeheader()
        writer.writerows(rows)


def run(args: argparse.Namespace) -> dict[str, object]:
    output_dir = Path(args.output_dir).expanduser().resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    if HERE.resolve() not in (output_dir, *output_dir.parents):
        raise ValueError("output directory must stay inside experiment 057")
    panel = read_csv_column(Path(args.panel_csv), "target_gene")
    official = read_csv_column(Path(args.genes_csv), "gene_name")
    panel_set = set(panel)
    banks = {
        "orion_hct116": Bank.load("orion_hct116", Path(args.hct116)),
        "orion_hek293t": Bank.load("orion_hek293t", Path(args.hek293t)),
        "k562_gwps": Bank.load("k562_gwps", Path(args.k562)),
        "kolf2": Bank.load("kolf2", Path(args.kolf), require_kolf=True),
    }
    targets = tuple(
        target
        for target in panel
        if all(target in bank.pert_pos for bank in banks.values())
    )
    if len(targets) < args.min_shared_targets:
        raise ValueError(f"only {len(targets)} targets shared by all four banks")

    # The operator is fitted and predictions are frozen before K562 effect
    # values are opened. K562 labels are used only to define the shared gallery.
    source_stats = [
        production_stats(
            banks[name],
            official,
            chunk_rows=args.chunk_rows,
            max_rss_gib=args.max_rss_gib,
        )
        for name in SOURCE_NAMES
    ]
    wcols, evals, evecs, operator_audit = build_source_only_operator(
        source_stats,
        official,
        panel_set,
        output_dir,
        basis_per_donor=args.basis_per_donor,
        operator_seed=args.operator_seed,
        chunk_rows=args.chunk_rows,
        max_rss_gib=args.max_rss_gib,
    )
    operator_genes = tuple(official[index] for index in wcols)
    source_union = set(operator_genes)
    k562_genes = set(banks[HOLDOUT_NAME].genes)
    kolf_genes = set(banks[AUXILIARY_NAME].genes)
    target_score_genes = tuple(
        gene for gene in official if gene in k562_genes and gene not in panel_set
    )
    shared_score_genes = tuple(
        gene
        for gene in target_score_genes
        if gene in source_union and gene in kolf_genes
    )
    off_operator_score_genes = tuple(
        gene for gene in target_score_genes if gene not in source_union
    )
    if any(gene in kolf_genes for gene in off_operator_score_genes):
        raise ValueError(
            "KOLF predicts K562-observed genes outside the source operator; "
            "the scale-invariant partial transform shortcut is invalid"
        )
    if len(shared_score_genes) < args.min_shared_genes:
        raise ValueError(f"only {len(shared_score_genes)} shared score genes")
    expected_score_order = shared_score_genes + off_operator_score_genes
    if set(expected_score_order) != set(target_score_genes):
        raise AssertionError("score-axis partition is incomplete")

    baseline, candidate, blend_audit = build_fixed_predictions(
        banks,
        targets,
        operator_genes,
        shared_score_genes,
        chunk_rows=args.chunk_rows,
        max_rss_gib=args.max_rss_gib,
    )
    # Predictions are now fixed. Only now load the K562 effect values.
    truth_full, truth_map = centered_panel_on_axis(
        banks[HOLDOUT_NAME],
        targets,
        expected_score_order,
        duplicate_policy="max_b_ctrl",
        chunk_rows=args.chunk_rows,
        max_rss_gib=args.max_rss_gib,
    )
    if len(truth_map.destination) != len(expected_score_order):
        raise AssertionError("K562 truth did not cover its declared score axis")
    truth_shared = truth_full[:, : len(shared_score_genes)]
    operator_pos = {gene: index for index, gene in enumerate(operator_genes)}
    operator_score_positions = np.asarray(
        [operator_pos[gene] for gene in shared_score_genes], dtype=np.int64
    )

    gammas = tuple(sorted(set(map(float, args.gammas))))
    if PRIMARY_GAMMA not in gammas or 0.0 not in gammas:
        raise ValueError("gamma sweep must include 0.0 and the fixed primary 0.65")
    gamma_results = []
    per_target_rows: list[dict[str, object]] = []
    for gamma in gammas:
        summary, rows = score_gamma(
            baseline,
            candidate,
            truth_full,
            truth_shared,
            operator_score_positions,
            evals,
            evecs,
            gamma,
            targets,
            args.query_folds,
        )
        gamma_results.append(summary)
        per_target_rows.extend(rows)
        primary = summary["target_observed_axis"]
        print(
            f"gamma={gamma:.2f} baseline={primary['baseline_pds']:.6f} "
            f"kolf={primary['kolf_pds']:.6f} delta={primary['delta']:+.6f}",
            flush=True,
        )
        enforce_rss(args.max_rss_gib, f"gamma {gamma}")

    expected_055 = json.loads(Path(args.exp055_results).read_text())["frozen_k562"]
    gamma_zero = next(item for item in gamma_results if item["gamma"] == 0.0)
    reproduction = {
        "expected_baseline_pds": float(expected_055["baseline_pds"]),
        "expected_kolf_pds": float(expected_055["pds"]),
        "observed_baseline_pds": gamma_zero["exp055_shared_axis"]["baseline_pds"],
        "observed_kolf_pds": gamma_zero["exp055_shared_axis"]["kolf_pds"],
    }
    reproduction["max_abs_error"] = max(
        abs(
            reproduction["observed_baseline_pds"]
            - reproduction["expected_baseline_pds"]
        ),
        abs(reproduction["observed_kolf_pds"] - reproduction["expected_kolf_pds"]),
    )
    reproduction["passes"] = reproduction["max_abs_error"] <= args.reproduction_tol
    if not reproduction["passes"]:
        raise AssertionError(f"experiment 055 reproduction failed: {reproduction}")

    primary = next(item for item in gamma_results if item["gamma"] == PRIMARY_GAMMA)
    result = {
        "experiment": "057-kolf-union-whiten",
        "status": "completed",
        "decision": {
            "primary_gamma": PRIMARY_GAMMA,
            "primary_delta": primary["target_observed_axis"]["delta"],
            "passes_positive_gate": primary["target_observed_axis"]["delta"] > 0,
            "rule": (
                "retain KOLF only if the a-priori gamma=0.65 target-observed-axis "
                "delta is positive; diagnostic gammas are not selectors"
            ),
        },
        "protocol": {
            "metric": "effect-level rank PDS",
            "held_out_target_role": HOLDOUT_NAME,
            "predictive_sources": list(SOURCE_NAMES),
            "auxiliary": AUXILIARY_NAME,
            "auxiliary_mode": "raw",
            "auxiliary_fraction": KOLF_FRACTION,
            "targets": len(targets),
            "full_retrieval_gallery": len(targets),
            "query_folds": args.query_folds,
            "query_fold_semantics": (
                "disjoint round-robin query subsets; every query ranks against "
                "the complete target gallery"
            ),
            "gammas": list(gammas),
            "primary_gamma_fixed_before_holdout": PRIMARY_GAMMA,
            "gamma_selection_on_k562": False,
            "operator_fit_rows": "source donors only",
            "count_emission": False,
            "cell_eval2": False,
            "submission": False,
            "renormalization": (
                "production preserve is a positive scalar per candidate and is "
                "omitted because effect-level cosine PDS is exactly scale invariant"
            ),
        },
        "axes": {
            "official_genes": len(official),
            "challenge_targets": len(panel),
            "operator_genes": len(operator_genes),
            "target_observed_score_genes": len(target_score_genes),
            "exp055_shared_score_genes": len(shared_score_genes),
            "target_score_genes_outside_operator": len(off_operator_score_genes),
            "target_score_genes_outside_operator_and_kolf": len(
                set(off_operator_score_genes) - kolf_genes
            ),
            "target_order_sha256": sha256_lines(targets),
            "target_observed_score_axis_sha256": sha256_lines(expected_score_order),
            "exp055_shared_score_axis_sha256": sha256_lines(shared_score_genes),
        },
        "leakage_audit": {
            "operator_sources_exact": list(SOURCE_NAMES),
            "k562_rows_in_operator": 0,
            "kolf_rows_in_operator": 0,
            "panel_coordinate_overlap_with_operator": operator_audit[
                "panel_overlap_with_operator"
            ],
            "predictions_frozen_before_k562_delta_open": True,
            "k562_used_for_candidate_selection": False,
            "fixed_rule_source": "experiments/055-kolf-direct/results.json",
            "exp055_reproduction": reproduction,
        },
        "blend_audit": blend_audit,
        "operator_audit": operator_audit,
        "input_fingerprints": {
            name: bank.fingerprint() for name, bank in banks.items()
        },
        "gamma_results": gamma_results,
        "peak_rss_gib": rss_gib(),
    }
    (output_dir / "results.json").write_text(json.dumps(result, indent=2) + "\n")
    write_csv(output_dir / "per_target_pds.csv", per_target_rows)
    write_csv(
        output_dir / "summary.csv",
        [
            {
                "gamma": item["gamma"],
                "is_primary": item["is_primary"],
                "baseline_pds": item["target_observed_axis"]["baseline_pds"],
                "kolf_pds": item["target_observed_axis"]["kolf_pds"],
                "delta": item["target_observed_axis"]["delta"],
                "shared_baseline_pds": item["exp055_shared_axis"]["baseline_pds"],
                "shared_kolf_pds": item["exp055_shared_axis"]["kolf_pds"],
                "shared_delta": item["exp055_shared_axis"]["delta"],
            }
            for item in gamma_results
        ],
    )
    print(
        f"wrote {output_dir / 'results.json'}; peak RSS={rss_gib():.2f} GiB",
        flush=True,
    )
    return result


def parser() -> argparse.ArgumentParser:
    ap = argparse.ArgumentParser(description=__doc__)
    inputs = HERE / "inputs"
    ap.add_argument("--hct116", default=str(inputs / "orion_hct116"))
    ap.add_argument("--hek293t", default=str(inputs / "orion_hek293t"))
    ap.add_argument("--k562", default=str(inputs / "k562_gwps"))
    ap.add_argument("--kolf", default=str(inputs / "kolf2"))
    ap.add_argument("--genes-csv", default=str(inputs / "gene_names.csv"))
    ap.add_argument("--panel-csv", default=str(inputs / "pert_counts.csv"))
    ap.add_argument("--exp055-results", default=str(inputs / "exp055-results.json"))
    ap.add_argument("--output-dir", default=str(HERE / "out"))
    ap.add_argument("--basis-per-donor", type=int, default=6000)
    ap.add_argument("--operator-seed", type=int, default=0)
    ap.add_argument("--gammas", nargs="+", type=float, default=list(DEFAULT_GAMMAS))
    ap.add_argument("--query-folds", type=int, default=5)
    ap.add_argument("--chunk-rows", type=int, default=64)
    ap.add_argument("--max-rss-gib", type=float, default=28.0)
    ap.add_argument("--min-shared-targets", type=int, default=240)
    ap.add_argument("--min-shared-genes", type=int, default=7000)
    ap.add_argument("--reproduction-tol", type=float, default=1e-9)
    return ap


def main() -> None:
    args = parser().parse_args()
    if args.basis_per_donor < 1 or args.chunk_rows < 1 or args.query_folds < 1:
        raise SystemExit("basis cap, chunk size, and query folds must be positive")
    if not 0 < args.max_rss_gib <= 29:
        raise SystemExit("--max-rss-gib must be in (0, 29]")
    run(args)


if __name__ == "__main__":
    main()
