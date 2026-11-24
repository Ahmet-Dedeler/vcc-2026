#!/usr/bin/env python3
"""Leak-aware, bounded-memory LODO screen of pooled Feng effects.

This is a screening experiment, not a submission or full cell-eval2 run. It
uses the three established somatic donor banks as leave-one-donor-out truth
and asks whether a small Feng blend improves raw-space PDS. The HCT116 and
HEK293T folds select one candidate; K562 is evaluated once after selection.

Large .npy donor matrices are read with ``os.pread`` in tiny contiguous row
blocks, avoiding multi-gigabyte memory mappings on the local Mac.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import resource
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Any

os.environ.setdefault("OPENBLAS_NUM_THREADS", "2")
os.environ.setdefault("VECLIB_MAXIMUM_THREADS", "2")
os.environ.setdefault("OMP_NUM_THREADS", "2")
os.environ.setdefault("MKL_NUM_THREADS", "2")

import numpy as np


ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
MAIN_PATHS = {
    "k562_gwps": ROOT / "data/donors/k562_gwps",
    "orion_hct116": ROOT / "data/donors/orion_hct116",
    "orion_hek293t": ROOT / "data/donors/orion_hek293t",
}
DEVELOPMENT_HELD = ("orion_hct116", "orion_hek293t")
FROZEN_HELD = "k562_gwps"
CONFIDENCE_POLICIES = (
    "all",
    "depth_ge_25",
    "depth_ge_50",
    "self_sig_0_1",
)


def peak_rss_gib() -> float:
    rss = float(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
    return rss / 1024**3 if sys.platform == "darwin" else rss / 1024**2


def enforce_memory_limit(limit_gib: float, stage: str) -> None:
    if limit_gib > 0 and peak_rss_gib() >= limit_gib:
        raise MemoryError(
            f"{stage}: peak RSS {peak_rss_gib():.3f} GiB reached {limit_gib:.3f} GiB"
        )


def read_lines(path: Path) -> tuple[str, ...]:
    return tuple(line.strip() for line in path.read_text().splitlines() if line.strip())


def read_csv_axis(path: Path, column: str) -> tuple[str, ...]:
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle)
        if reader.fieldnames is None or column not in reader.fieldnames:
            raise ValueError(f"{path}: missing {column!r}")
        values = tuple(row[column].strip() for row in reader if row[column].strip())
    if len(values) != len(set(values)):
        raise ValueError(f"{path}: duplicate {column!r} values")
    return values


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while block := handle.read(1 << 20):
            digest.update(block)
    return digest.hexdigest()


@dataclass(frozen=True)
class NpyLayout:
    path: Path
    shape: tuple[int, int]
    dtype: np.dtype[Any]
    data_offset: int
    row_bytes: int

    @classmethod
    def inspect(cls, path: Path) -> "NpyLayout":
        with path.open("rb") as handle:
            version = np.lib.format.read_magic(handle)
            if version == (1, 0):
                shape, fortran, dtype = np.lib.format.read_array_header_1_0(handle)
            elif version in ((2, 0), (3, 0)):
                shape, fortran, dtype = np.lib.format.read_array_header_2_0(handle)
            else:
                raise ValueError(f"{path}: unsupported NPY version {version}")
            offset = handle.tell()
        if fortran:
            raise ValueError(f"{path}: Fortran-order arrays are unsupported")
        if len(shape) != 2:
            raise ValueError(f"{path}: expected a matrix, got {shape}")
        dtype = np.dtype(dtype)
        expected_size = offset + int(np.prod(shape)) * dtype.itemsize
        if path.stat().st_size != expected_size:
            raise ValueError(
                f"{path}: byte size {path.stat().st_size} != {expected_size}"
            )
        return cls(
            path=path,
            shape=(int(shape[0]), int(shape[1])),
            dtype=dtype,
            data_offset=offset,
            row_bytes=int(shape[1]) * dtype.itemsize,
        )

    def _read_exact(self, descriptor: int, offset: int, size: int) -> bytes:
        chunks: list[bytes] = []
        remaining = size
        while remaining:
            payload = os.pread(descriptor, remaining, offset + size - remaining)
            if not payload:
                raise EOFError(f"{self.path}: short read at byte {offset}")
            chunks.append(payload)
            remaining -= len(payload)
        return b"".join(chunks)

    def read_block(self, descriptor: int, start: int, stop: int) -> np.ndarray:
        if not 0 <= start < stop <= self.shape[0]:
            raise IndexError((start, stop, self.shape))
        payload = self._read_exact(
            descriptor,
            self.data_offset + start * self.row_bytes,
            (stop - start) * self.row_bytes,
        )
        return np.frombuffer(payload, dtype=self.dtype).reshape(
            stop - start, self.shape[1]
        )

    def read_row(self, descriptor: int, row: int) -> np.ndarray:
        return self.read_block(descriptor, row, row + 1)[0]


@dataclass(frozen=True)
class Bank:
    name: str
    path: Path
    genes: tuple[str, ...]
    perts: tuple[str, ...]
    gene_pos: dict[str, int]
    pert_pos: dict[str, int]
    layout: NpyLayout
    duplicate_gene_symbols: int

    @classmethod
    def load(cls, name: str, path: Path) -> "Bank":
        genes = read_lines(path / "genes.txt")
        perts = read_lines(path / "perts.txt")
        layout = NpyLayout.inspect(path / "delta.npy")
        if layout.shape != (len(perts), len(genes)):
            raise ValueError(
                f"{name}: delta shape {layout.shape} != {(len(perts), len(genes))}"
            )
        b_ctrl = np.load(path / "b_ctrl.npy")
        if b_ctrl.shape != (len(genes),):
            raise ValueError(f"{name}: malformed b_ctrl.npy")
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
        if len(perts) != len(set(perts)):
            raise ValueError(f"{name}: duplicate perturbation labels")
        return cls(
            name=name,
            path=path,
            genes=genes,
            perts=perts,
            gene_pos=gene_pos,
            pert_pos={pert: index for index, pert in enumerate(perts)},
            layout=layout,
            duplicate_gene_symbols=len(duplicates),
        )

    def centered_selected(
        self,
        targets: tuple[str, ...],
        features: tuple[str, ...],
        *,
        chunk_rows: int,
        max_rss_gib: float,
    ) -> np.ndarray:
        columns = np.asarray([self.gene_pos[gene] for gene in features], dtype=np.int64)
        mean = np.zeros(len(features), dtype=np.float64)
        descriptor = os.open(self.layout.path, os.O_RDONLY)
        try:
            for start in range(0, self.layout.shape[0], chunk_rows):
                stop = min(start + chunk_rows, self.layout.shape[0])
                block = self.layout.read_block(descriptor, start, stop)
                mean += np.asarray(block[:, columns], dtype=np.float64).sum(axis=0)
                enforce_memory_limit(max_rss_gib, f"{self.name} center at row {stop}")
            mean /= self.layout.shape[0]

            selected = np.empty((len(targets), len(features)), dtype=np.float64)
            for output_row, target in enumerate(targets):
                row = self.pert_pos[target]
                selected[output_row] = (
                    np.asarray(
                        self.layout.read_row(descriptor, row)[columns],
                        dtype=np.float64,
                    )
                    - mean
                )
        finally:
            os.close(descriptor)
        enforce_memory_limit(max_rss_gib, f"{self.name} selected effects")
        return selected


@dataclass(frozen=True)
class Feng:
    path: Path
    genes: tuple[str, ...]
    targets: tuple[str, ...]
    target_pos: dict[str, int]
    centered_lfc: np.ndarray
    pval_adj: np.ndarray
    observed: np.ndarray
    metadata: dict[str, Any]
    cell_depth: dict[str, int]

    @classmethod
    def load(cls, path: Path, coverage_path: Path) -> "Feng":
        genes = read_lines(path / "genes.txt")
        targets = read_lines(path / "targets.txt")
        centered_lfc = np.load(path / "centered_lfc.npy", mmap_mode="r")
        pval_adj = np.load(path / "pval_adj.npy", mmap_mode="r")
        observed = np.load(path / "observed.npy", mmap_mode="r")
        metadata = json.loads((path / "meta.json").read_text())
        coverage = json.loads(coverage_path.read_text())
        cell_depth = coverage["genomewide_merged"]["vcc_effect_cells"]
        expected = (len(targets), len(genes))
        for name, array in (
            ("centered_lfc", centered_lfc),
            ("pval_adj", pval_adj),
            ("observed", observed),
        ):
            if array.shape != expected:
                raise ValueError(f"Feng {name} shape {array.shape} != {expected}")
        if metadata["selection"]["covered_targets"] != len(targets):
            raise ValueError("Feng metadata target count mismatch")
        return cls(
            path=path,
            genes=genes,
            targets=targets,
            target_pos={target: index for index, target in enumerate(targets)},
            centered_lfc=centered_lfc,
            pval_adj=pval_adj,
            observed=observed,
            metadata=metadata,
            cell_depth={str(k): int(v) for k, v in cell_depth.items()},
        )


def row_normalize(matrix: np.ndarray) -> np.ndarray:
    return matrix / (np.linalg.norm(matrix, axis=1, keepdims=True) + 1e-12)


def pds_rows(
    prediction: np.ndarray,
    truth: np.ndarray,
    query_mask: np.ndarray | None = None,
) -> float:
    scores = row_normalize(prediction) @ row_normalize(truth).T
    own = np.diag(scores)[:, None]
    n = scores.shape[1]
    ranks = (scores > own).sum(axis=1) + 0.5 * ((scores == own).sum(axis=1) - 1)
    values = 1.0 - ranks / (n - 1)
    if query_mask is not None:
        values = values[query_mask]
    if not values.size:
        return float("nan")
    return float(np.mean(values))


def scale_matrix(matrix: np.ndarray, reference_norm: float) -> np.ndarray:
    norm = float(np.median(np.linalg.norm(matrix, axis=1))) or 1.0
    return matrix * (reference_norm / norm)


def make_baseline(parts: list[np.ndarray]) -> tuple[np.ndarray, float]:
    norms = [float(np.median(np.linalg.norm(part, axis=1))) for part in parts]
    reference = float(np.median(norms)) or 1.0
    return sum(scale_matrix(part, reference) for part in parts) / len(parts), reference


def prepare_feng(
    feng: Feng,
    targets: tuple[str, ...],
    features: tuple[str, ...],
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    official_pos = {gene: index for index, gene in enumerate(feng.genes)}
    feature_columns = np.asarray(
        [official_pos[gene] for gene in features], dtype=np.int64
    )
    values = np.zeros((len(targets), len(features)), dtype=np.float64)
    available = np.zeros(values.shape, dtype=np.bool_)
    covered = np.zeros(len(targets), dtype=np.bool_)
    self_significant = np.zeros(len(targets), dtype=np.bool_)
    for output_row, target in enumerate(targets):
        source_row = feng.target_pos.get(target)
        if source_row is None:
            continue
        row_observed = np.asarray(
            feng.observed[source_row, feature_columns], dtype=np.bool_
        )
        values[output_row, row_observed] = np.asarray(
            feng.centered_lfc[source_row, feature_columns[row_observed]],
            dtype=np.float64,
        )
        available[output_row, row_observed] = True
        covered[output_row] = True
        target_column = official_pos.get(target)
        if (
            target_column is not None
            and bool(feng.observed[source_row, target_column])
            and float(feng.pval_adj[source_row, target_column]) < 0.1
        ):
            self_significant[output_row] = True
    return values, available, covered, self_significant


def policy_mask(
    policy: str,
    targets: tuple[str, ...],
    covered: np.ndarray,
    self_significant: np.ndarray,
    cell_depth: dict[str, int],
) -> np.ndarray:
    if policy == "all":
        return covered.copy()
    if policy == "self_sig_0_1":
        return covered & self_significant
    if policy.startswith("depth_ge_"):
        threshold = int(policy.removeprefix("depth_ge_"))
        return covered & np.asarray(
            [cell_depth.get(target, 0) >= threshold for target in targets],
            dtype=np.bool_,
        )
    raise ValueError(f"unknown confidence policy {policy!r}")


def blend_feng(
    baseline: np.ndarray,
    feng_values: np.ndarray,
    available: np.ndarray,
    active_targets: np.ndarray,
    fraction: float,
) -> tuple[np.ndarray, dict[str, float]]:
    active_available = available & active_targets[:, None]
    if not active_available.any():
        raise ValueError("candidate has no active Feng coordinates")
    active_rows = np.flatnonzero(active_targets)
    feature_mask = available[active_rows].any(axis=0)
    baseline_reference = float(
        np.median(np.linalg.norm(baseline[active_rows][:, feature_mask], axis=1))
    )
    feng_reference = float(
        np.median(np.linalg.norm(feng_values[active_rows][:, feature_mask], axis=1))
    )
    if feng_reference <= 0:
        raise ValueError("Feng median effect norm is zero")
    scale = baseline_reference / feng_reference
    candidate = baseline.copy()
    candidate[active_available] = (1.0 - fraction) * baseline[
        active_available
    ] + fraction * scale * feng_values[active_available]
    return candidate, {
        "baseline_overlap_median_norm": baseline_reference,
        "feng_overlap_median_norm": feng_reference,
        "feng_scale": scale,
    }


def evaluate_fold(
    held: str,
    effects: dict[str, np.ndarray],
    feng_values: np.ndarray,
    available: np.ndarray,
    active_targets: np.ndarray,
    fraction: float,
) -> dict[str, Any]:
    sources = [
        effects[name]
        for name in ("k562_gwps", "orion_hct116", "orion_hek293t")
        if name != held
    ]
    baseline, _ = make_baseline(sources)
    truth = effects[held]
    candidate, scaling = blend_feng(
        baseline, feng_values, available, active_targets, fraction
    )
    baseline_all = pds_rows(baseline, truth)
    candidate_all = pds_rows(candidate, truth)
    baseline_active = pds_rows(baseline, truth, active_targets)
    candidate_active = pds_rows(candidate, truth, active_targets)
    inactive = ~active_targets
    return {
        "held": held,
        "fraction": fraction,
        "active_targets": int(active_targets.sum()),
        "inactive_targets": int(inactive.sum()),
        "baseline_pds_all": baseline_all,
        "candidate_pds_all": candidate_all,
        "delta_all": candidate_all - baseline_all,
        "baseline_pds_active": baseline_active,
        "candidate_pds_active": candidate_active,
        "delta_active": candidate_active - baseline_active,
        "inactive_exactly_unchanged": bool(
            np.array_equal(candidate[inactive], baseline[inactive])
        ),
        "scaling": scaling,
    }


def choose_candidate(rows: list[dict[str, Any]]) -> dict[str, Any]:
    grouped: dict[tuple[str, float], list[dict[str, Any]]] = {}
    for row in rows:
        grouped.setdefault((str(row["policy"]), float(row["fraction"])), []).append(row)
    candidates = []
    for (policy, fraction), group in grouped.items():
        if {row["held"] for row in group} != set(DEVELOPMENT_HELD):
            raise ValueError("candidate does not have both development folds")
        deltas = [float(row["delta_all"]) for row in group]
        candidates.append(
            {
                "policy": policy,
                "fraction": fraction,
                "development_mean_delta": float(np.mean(deltas)),
                "development_min_delta": float(np.min(deltas)),
                "active_targets": min(int(row["active_targets"]) for row in group),
            }
        )
    return max(
        candidates,
        key=lambda row: (
            row["development_mean_delta"],
            row["development_min_delta"],
            row["active_targets"],
            -row["fraction"],
        ),
    )


def run(args: argparse.Namespace) -> dict[str, Any]:
    panel = read_csv_axis(args.panel_csv, "target_gene")
    official = read_csv_axis(args.genes_csv, "gene_name")
    panel_set = set(panel)
    banks = {name: Bank.load(name, path) for name, path in MAIN_PATHS.items()}
    feng = Feng.load(args.feng, args.coverage)
    if feng.genes != official:
        raise ValueError("Feng compact artifact is not on the canonical VCC gene axis")

    targets = tuple(
        target
        for target in panel
        if all(target in bank.pert_pos for bank in banks.values())
    )
    features = tuple(
        gene
        for gene in official
        if gene not in panel_set
        and all(gene in bank.gene_pos for bank in banks.values())
    )
    if len(targets) < args.min_gallery_targets:
        raise ValueError(f"only {len(targets)} shared main-donor targets")
    if len(features) < args.min_features:
        raise ValueError(f"only {len(features)} shared non-panel genes")

    effects: dict[str, np.ndarray] = {}
    for name, bank in banks.items():
        effects[name] = bank.centered_selected(
            targets,
            features,
            chunk_rows=args.chunk_rows,
            max_rss_gib=args.max_rss_gib,
        )
        print(
            f"prepared {name}: shape={effects[name].shape}, "
            f"peak_rss={peak_rss_gib():.3f} GiB",
            flush=True,
        )

    feng_values, available, covered, self_significant = prepare_feng(
        feng, targets, features
    )
    development_rows: list[dict[str, Any]] = []
    for policy in CONFIDENCE_POLICIES:
        active = policy_mask(
            policy, targets, covered, self_significant, feng.cell_depth
        )
        if active.sum() < args.min_active_targets:
            continue
        for held in DEVELOPMENT_HELD:
            for fraction in args.fractions:
                row = evaluate_fold(
                    held,
                    effects,
                    feng_values,
                    available,
                    active,
                    fraction,
                )
                row["policy"] = policy
                development_rows.append(row)

    selected = choose_candidate(development_rows)
    frozen_active = policy_mask(
        str(selected["policy"]),
        targets,
        covered,
        self_significant,
        feng.cell_depth,
    )
    frozen = evaluate_fold(
        FROZEN_HELD,
        effects,
        feng_values,
        available,
        frozen_active,
        float(selected["fraction"]),
    )
    frozen["policy"] = selected["policy"]
    three_fold_delta = float(
        (2.0 * float(selected["development_mean_delta"]) + float(frozen["delta_all"]))
        / 3.0
    )
    selected = {
        **selected,
        "selection_rule": (
            "maximize mean full-gallery PDS delta on held HCT116 and HEK293T; "
            "tie-break by worst fold, active-target count, then smaller fraction"
        ),
        "frozen_k562_delta": frozen["delta_all"],
        "three_fold_macro_delta": three_fold_delta,
        "passes_screen": (
            float(selected["development_min_delta"]) >= -0.002
            and float(frozen["delta_all"]) >= 0.0
            and three_fold_delta > 0.0
        ),
    }

    result = {
        "protocol": {
            "metric": "raw-space effect-level PDS",
            "gallery_targets": len(targets),
            "non_panel_features": len(features),
            "feng_covered_gallery_targets": int(covered.sum()),
            "feng_measured_features": int(available.any(axis=0).sum()),
            "source_centering": "per-gene mean across every perturbation in each donor",
            "feng_centering": "per-gene mean across all 6,673 Feng targets",
            "missing_target_policy": "preserve baseline exactly",
            "missing_gene_policy": "preserve baseline exactly",
            "fractions": args.fractions,
            "confidence_policies": list(CONFIDENCE_POLICIES),
            "development_folds": list(DEVELOPMENT_HELD),
            "frozen_holdout": FROZEN_HELD,
            "no_whitening": True,
            "no_count_emission": True,
            "no_cell_eval2": True,
            "no_submission": True,
            "large_array_io": "bounded contiguous os.pread blocks",
            "duplicate_gene_symbols_resolved_by_max_b_ctrl": {
                name: bank.duplicate_gene_symbols for name, bank in banks.items()
            },
        },
        "inputs": {
            "feng_meta_sha256": sha256_file(args.feng / "meta.json"),
            "coverage_sha256": sha256_file(args.coverage),
        },
        "targets": list(targets),
        "feature_axis_sha256": hashlib.sha256(
            ("\n".join(features) + "\n").encode()
        ).hexdigest(),
        "development_rows": development_rows,
        "selected_on_hct_and_hek_only": selected,
        "frozen_k562": frozen,
        "peak_rss_gib": peak_rss_gib(),
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2) + "\n")
    with args.output.with_suffix(".csv").open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "policy",
                "held",
                "fraction",
                "active_targets",
                "inactive_targets",
                "baseline_pds_all",
                "candidate_pds_all",
                "delta_all",
                "baseline_pds_active",
                "candidate_pds_active",
                "delta_active",
                "inactive_exactly_unchanged",
            ],
            extrasaction="ignore",
        )
        writer.writeheader()
        writer.writerows(development_rows)
    enforce_memory_limit(args.max_rss_gib, "completed LODO")
    print(json.dumps(selected, indent=2))
    return result


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--feng",
        type=Path,
        default=ROOT / "data/external/feng-ipsc/genomewide_vcc_panel",
    )
    parser.add_argument(
        "--coverage",
        type=Path,
        default=ROOT / "data/external/feng-ipsc/metadata_coverage.json",
    )
    parser.add_argument(
        "--panel-csv",
        type=Path,
        default=ROOT / "data/challenge/pert_counts.csv",
    )
    parser.add_argument(
        "--genes-csv",
        type=Path,
        default=ROOT / "data/challenge/gene_names.csv",
    )
    parser.add_argument(
        "--fractions",
        nargs="+",
        type=float,
        default=[0.02, 0.05, 0.10, 0.20, 0.35],
    )
    parser.add_argument("--chunk-rows", type=int, default=4)
    parser.add_argument("--max-rss-gib", type=float, default=0.90)
    parser.add_argument("--min-gallery-targets", type=int, default=200)
    parser.add_argument("--min-features", type=int, default=7_000)
    parser.add_argument("--min-active-targets", type=int, default=10)
    parser.add_argument("--output", type=Path, default=HERE / "lodo_results.json")
    return parser


def main() -> None:
    args = build_parser().parse_args()
    if args.chunk_rows < 1:
        raise SystemExit("--chunk-rows must be positive")
    if not 0 < args.max_rss_gib <= 1.0:
        raise SystemExit("--max-rss-gib must be in (0, 1]")
    if any(not 0 < fraction < 1 for fraction in args.fractions):
        raise SystemExit("all fractions must be in (0, 1)")
    run(args)


if __name__ == "__main__":
    main()
