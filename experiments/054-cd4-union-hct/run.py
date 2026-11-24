"""054 — fixed CD4 mean3 blend under the clean HCT union whitener.

This is a deliberately narrow, effect-level PDS experiment. It compares the
current K562-gwps + Orion-HEK293T baseline with the already selected 50% CD4
three-condition mean, blended *before* the exact clean-HCT union operator from
experiment 040.

The script never fits a whitener, emits counts, invokes cell-eval2, or submits
anything. Large donor arrays and the extracted whitener eigenvectors are
memory-mapped; source statistics and whitening are computed in bounded blocks.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import resource
import shutil
import sys
import tempfile
import zipfile
from collections.abc import Iterable, Sequence
from dataclasses import dataclass
from pathlib import Path

# Bound BLAS before importing NumPy. The external process-tree guard remains
# required for a real-data invocation.
os.environ.setdefault("OPENBLAS_NUM_THREADS", "2")
os.environ.setdefault("VECLIB_MAXIMUM_THREADS", "2")
os.environ.setdefault("OMP_NUM_THREADS", "2")
os.environ.setdefault("MKL_NUM_THREADS", "2")

import anndata as ad
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent

SOURCE_NAMES = ("k562_gwps", "orion_hek293t")
CD4_NAMES = ("cd4_rest", "cd4_stim8hr", "cd4_stim48hr")
CD4_MODES = ("raw", "gene_norm")
MISSING_GENE_POLICIES = ("preserve_baseline", "zero_vote")
RENORMALIZATION_POLICIES = ("production_self", "baseline_fixed")

CONTROL = "non-targeting"
PERT_COL = "target_gene"
TARGET_SUM = 5e4
EFFECT_SCALE = 6.817
BLEND_FRACTION = 0.5
WHITEN_GAMMA = 0.65


def peak_rss_gib() -> float:
    """Return peak process RSS in GiB on macOS or Linux."""
    rss = float(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
    return rss / 1024**3 if sys.platform == "darwin" else rss / 1024**2


def enforce_memory_limit(limit_gib: float, stage: str) -> None:
    """Fail before a process with an unexpectedly large live set can continue."""
    if limit_gib > 0 and peak_rss_gib() >= limit_gib:
        raise MemoryError(
            f"{stage}: peak RSS {peak_rss_gib():.3f} GiB reached the "
            f"{limit_gib:.3f} GiB internal limit"
        )


def read_lines(path: Path, *, require_unique: bool = True) -> tuple[str, ...]:
    values = tuple(
        line.strip() for line in path.read_text().splitlines() if line.strip()
    )
    if require_unique and len(values) != len(set(values)):
        raise ValueError(f"{path}: duplicate labels")
    return values


def sha256_lines(values: Iterable[str]) -> str:
    return hashlib.sha256("\n".join(values).encode()).hexdigest()


def sha256_file(path: Path, block_bytes: int = 8 << 20) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while block := handle.read(block_bytes):
            digest.update(block)
    return digest.hexdigest()


def read_csv_column(path: Path, column: str) -> tuple[str, ...]:
    frame = pd.read_csv(path)
    if column not in frame:
        raise ValueError(f"{path}: required column {column!r} is missing")
    values = tuple(
        value
        for value in frame[column].dropna().astype(str).str.strip().tolist()
        if value
    )
    if len(values) != len(set(values)):
        raise ValueError(f"{path}: {column!r} contains duplicate labels")
    return values


@dataclass(frozen=True)
class Bank:
    """A memory-mapped perturbation-effect bank."""

    name: str
    path: Path
    genes: tuple[str, ...]
    perts: tuple[str, ...]
    pert_pos: dict[str, int]
    delta_path: Path
    b_ctrl_path: Path
    shape: tuple[int, int]

    @classmethod
    def load(cls, name: str, path: Path, *, require_cd4_scale: bool = False) -> Bank:
        path = path.expanduser().resolve()
        genes = read_lines(path / "genes.txt", require_unique=False)
        perts = read_lines(path / "perts.txt")
        delta_path = path / "delta.npy"
        b_ctrl_path = path / "b_ctrl.npy"
        delta = np.load(delta_path, mmap_mode="r")
        b_ctrl = np.load(b_ctrl_path, mmap_mode="r")
        if not isinstance(delta, np.memmap):
            raise TypeError(f"{name}: delta.npy is not memory mapped")
        if delta.shape != (len(perts), len(genes)):
            raise ValueError(
                f"{name}: delta shape {delta.shape} != "
                f"labels {(len(perts), len(genes))}"
            )
        if b_ctrl.shape != (len(genes),):
            raise ValueError(f"{name}: b_ctrl shape {b_ctrl.shape} != {(len(genes),)}")
        if not np.issubdtype(delta.dtype, np.floating):
            raise TypeError(f"{name}: delta dtype {delta.dtype} is not floating")
        if require_cd4_scale:
            metadata = json.loads((path / "meta.json").read_text())
            definition = str(metadata.get("effect_definition", ""))
            if "ln(2)" not in definition or "log2 fold change" not in definition:
                raise ValueError(
                    f"{name}: unrecognized CD4 effect scale {definition!r}"
                )
        return cls(
            name=name,
            path=path,
            genes=genes,
            perts=perts,
            pert_pos={pert: index for index, pert in enumerate(perts)},
            delta_path=delta_path,
            b_ctrl_path=b_ctrl_path,
            shape=tuple(map(int, delta.shape)),
        )


@dataclass(frozen=True)
class AxisMap:
    """Donor columns mapped to the official HCT axis.

    ``all_*`` preserves duplicate donor symbols because experiment 031 counted
    duplicate mapped columns in its donor norm. ``unique_*`` uses the last
    source occurrence for assignment, matching NumPy's repeated-index
    assignment in the existing baseline implementation.
    """

    all_source_columns: np.ndarray
    all_target_positions: np.ndarray
    unique_source_columns: np.ndarray
    unique_target_positions: np.ndarray
    duplicate_symbols: tuple[str, ...]


def map_bank_to_axis(bank: Bank, official_genes: Sequence[str]) -> AxisMap:
    target_pos = {gene: index for index, gene in enumerate(official_genes)}
    all_source: list[int] = []
    all_target: list[int] = []
    last_source: dict[int, int] = {}
    seen: set[int] = set()
    duplicates: set[str] = set()
    for source_column, gene in enumerate(bank.genes):
        destination = target_pos.get(gene)
        if destination is None:
            continue
        all_source.append(source_column)
        all_target.append(destination)
        if destination in seen:
            duplicates.add(gene)
        seen.add(destination)
        last_source[destination] = source_column
    if not all_source:
        raise ValueError(f"{bank.name}: no genes map to the official HCT axis")
    unique_target = np.asarray(sorted(last_source), dtype=np.int64)
    unique_source = np.asarray(
        [last_source[position] for position in unique_target], dtype=np.int64
    )
    return AxisMap(
        all_source_columns=np.asarray(all_source, dtype=np.int64),
        all_target_positions=np.asarray(all_target, dtype=np.int64),
        unique_source_columns=unique_source,
        unique_target_positions=unique_target,
        duplicate_symbols=tuple(sorted(duplicates)),
    )


def _stream_column_mean(
    bank: Bank,
    columns: np.ndarray,
    *,
    chunk_rows: int,
    max_rss_gib: float,
) -> np.ndarray:
    delta = np.load(bank.delta_path, mmap_mode="r")
    total = np.zeros(len(columns), dtype=np.float64)
    for lo in range(0, bank.shape[0], chunk_rows):
        hi = min(lo + chunk_rows, bank.shape[0])
        block = np.asarray(delta[lo:hi, columns], dtype=np.float64)
        total += block.sum(axis=0, dtype=np.float64)
        enforce_memory_limit(max_rss_gib, f"{bank.name} source centering")
    return total / bank.shape[0]


def _stream_centered_statistics(
    bank: Bank,
    columns: np.ndarray,
    mean: np.ndarray,
    *,
    need_gene_scale: bool,
    chunk_rows: int,
    max_rss_gib: float,
) -> tuple[float, np.ndarray | None, float | None]:
    """Return raw median row norm and optional mean-absolute gene normalization."""
    delta = np.load(bank.delta_path, mmap_mode="r")
    raw_norms = np.empty(bank.shape[0], dtype=np.float64)
    absolute = np.zeros(len(columns), dtype=np.float64) if need_gene_scale else None
    for lo in range(0, bank.shape[0], chunk_rows):
        hi = min(lo + chunk_rows, bank.shape[0])
        centered = np.asarray(delta[lo:hi, columns], dtype=np.float64) - mean
        raw_norms[lo:hi] = np.linalg.norm(centered, axis=1)
        if absolute is not None:
            absolute += np.abs(centered).sum(axis=0, dtype=np.float64)
        enforce_memory_limit(max_rss_gib, f"{bank.name} source scale")
    raw_median = float(np.median(raw_norms)) or 1.0
    if absolute is None:
        return raw_median, None, None

    mean_abs = absolute / bank.shape[0]
    live = mean_abs[mean_abs > 0]
    floor = float(np.percentile(live, 25.0)) if live.size else 1.0
    gene_scale = np.maximum(mean_abs, floor)
    normalized_norms = np.empty(bank.shape[0], dtype=np.float64)
    for lo in range(0, bank.shape[0], chunk_rows):
        hi = min(lo + chunk_rows, bank.shape[0])
        centered = np.asarray(delta[lo:hi, columns], dtype=np.float64) - mean
        normalized_norms[lo:hi] = np.linalg.norm(centered / gene_scale, axis=1)
        enforce_memory_limit(max_rss_gib, f"{bank.name} normalized source scale")
    normalized_median = float(np.median(normalized_norms)) or 1.0
    return raw_median, gene_scale, normalized_median


@dataclass
class PreparedSource:
    name: str
    effects: np.ndarray
    target_available: np.ndarray
    mean_on_axis: np.ndarray
    raw_median_norm: float
    mapped_gene_count: int
    duplicate_symbols: tuple[str, ...]


def prepare_main_source(
    bank: Bank,
    targets: Sequence[str],
    official_genes: Sequence[str],
    *,
    chunk_rows: int,
    max_rss_gib: float,
) -> PreparedSource:
    """Reproduce experiment 031's source centering and donor norm exactly."""
    mapping = map_bank_to_axis(bank, official_genes)
    mean_all = _stream_column_mean(
        bank,
        mapping.all_source_columns,
        chunk_rows=chunk_rows,
        max_rss_gib=max_rss_gib,
    )
    raw_median, _, _ = _stream_centered_statistics(
        bank,
        mapping.all_source_columns,
        mean_all,
        need_gene_scale=False,
        chunk_rows=chunk_rows,
        max_rss_gib=max_rss_gib,
    )

    # Pick the last occurrence of each duplicate symbol, exactly as the old
    # repeated-index assignment did, while retaining duplicates in raw_median.
    all_column_pos = {
        int(source): index
        for index, source in enumerate(mapping.all_source_columns.tolist())
    }
    unique_mean = np.asarray(
        [
            mean_all[all_column_pos[int(source)]]
            for source in mapping.unique_source_columns
        ],
        dtype=np.float64,
    )
    mean_on_axis = np.zeros(len(official_genes), dtype=np.float64)
    mean_on_axis[mapping.unique_target_positions] = unique_mean

    effects = np.zeros((len(targets), len(official_genes)), dtype=np.float64)
    available = np.zeros(len(targets), dtype=bool)
    delta = np.load(bank.delta_path, mmap_mode="r")
    for target_index, target in enumerate(targets):
        row = bank.pert_pos.get(target)
        if row is None:
            continue
        values = np.asarray(delta[row, mapping.unique_source_columns], dtype=np.float64)
        effects[target_index, mapping.unique_target_positions] = (
            values - unique_mean
        ) * (EFFECT_SCALE / raw_median)
        available[target_index] = True
    enforce_memory_limit(max_rss_gib, f"{bank.name} target extraction")
    return PreparedSource(
        name=bank.name,
        effects=effects,
        target_available=available,
        mean_on_axis=mean_on_axis,
        raw_median_norm=raw_median,
        mapped_gene_count=len(mapping.unique_target_positions),
        duplicate_symbols=mapping.duplicate_symbols,
    )


def build_current_baseline(
    sources: Sequence[PreparedSource],
) -> tuple[np.ndarray, np.ndarray]:
    """Build the exact equal-source experiment-031 baseline shape."""
    if not sources:
        raise ValueError("at least one source donor is required")
    shape = sources[0].effects.shape
    if any(source.effects.shape != shape for source in sources):
        raise ValueError("source effect matrices have different shapes")

    baseline = np.zeros(shape, dtype=np.float64)
    votes = np.zeros(shape[0], dtype=np.int64)
    for source in sources:
        rows = source.target_available
        baseline[rows] += source.effects[rows]
        votes += rows.astype(np.int64)
    covered = votes > 0
    baseline[covered] /= votes[covered, None]

    # Exact MultiDonorModel fallback for a target absent from every source.
    missing = ~covered
    if missing.any():
        mean_response = sum(source.mean_on_axis for source in sources) / len(sources)
        mean_norm = float(np.linalg.norm(mean_response)) or 1.0
        baseline[missing] = EFFECT_SCALE * mean_response / mean_norm
    return baseline, votes


@dataclass
class CD4Prepared:
    mode: str
    mean_effect: np.ndarray
    gene_votes: np.ndarray
    target_condition_count: np.ndarray
    condition_details: list[dict[str, object]]


def prepare_cd4_mean(
    banks: Sequence[Bank],
    targets: Sequence[str],
    official_genes: Sequence[str],
    panel_targets: set[str],
    mode: str,
    *,
    chunk_rows: int,
    max_rss_gib: float,
) -> CD4Prepared:
    """Center each condition on all its perturbations, then form mean3.

    Partially covered targets average the available conditions. Targets absent
    from all three conditions remain explicitly unavailable and are never
    allowed to change the baseline.
    """
    if mode not in CD4_MODES:
        raise ValueError(f"unknown CD4 mode {mode!r}")
    official_pos = {gene: index for index, gene in enumerate(official_genes)}
    value_sum = np.zeros((len(targets), len(official_genes)), dtype=np.float64)
    gene_votes = np.zeros(value_sum.shape, dtype=np.uint8)
    condition_count = np.zeros(len(targets), dtype=np.uint8)
    details: list[dict[str, object]] = []

    for bank in banks:
        source_columns: list[int] = []
        target_positions: list[int] = []
        for source_column, gene in enumerate(bank.genes):
            destination = official_pos.get(gene)
            if destination is None or gene in panel_targets:
                continue
            source_columns.append(source_column)
            target_positions.append(destination)
        if not source_columns:
            raise ValueError(f"{bank.name}: no non-panel genes map to HCT")
        source_columns_array = np.asarray(source_columns, dtype=np.int64)
        target_positions_array = np.asarray(target_positions, dtype=np.int64)
        if len(np.unique(target_positions_array)) != len(target_positions_array):
            raise ValueError(f"{bank.name}: compact CD4 gene symbols are not unique")

        mean = _stream_column_mean(
            bank,
            source_columns_array,
            chunk_rows=chunk_rows,
            max_rss_gib=max_rss_gib,
        )
        raw_median, gene_scale, normalized_median = _stream_centered_statistics(
            bank,
            source_columns_array,
            mean,
            need_gene_scale=mode == "gene_norm",
            chunk_rows=chunk_rows,
            max_rss_gib=max_rss_gib,
        )
        if mode == "gene_norm":
            assert gene_scale is not None and normalized_median is not None
            source_multiplier = (raw_median / normalized_median) * (
                EFFECT_SCALE / raw_median
            )
        else:
            source_multiplier = EFFECT_SCALE / raw_median

        delta = np.load(bank.delta_path, mmap_mode="r")
        covered = 0
        for target_index, target in enumerate(targets):
            row = bank.pert_pos.get(target)
            if row is None:
                continue
            centered = (
                np.asarray(delta[row, source_columns_array], dtype=np.float64) - mean
            )
            if mode == "gene_norm":
                centered = centered / gene_scale
            value_sum[target_index, target_positions_array] += (
                centered * source_multiplier
            )
            gene_votes[target_index, target_positions_array] += 1
            condition_count[target_index] += 1
            covered += 1
        details.append(
            {
                "condition": bank.name,
                "source_background_perturbations": bank.shape[0],
                "mapped_non_panel_genes": len(source_columns_array),
                "covered_hct_targets": covered,
                "raw_median_centered_norm": raw_median,
                "gene_normalized_median_before_rescale": normalized_median,
                "effect_scale_after_rescale": EFFECT_SCALE,
            }
        )
        enforce_memory_limit(max_rss_gib, f"{bank.name} CD4 target extraction")

    available = gene_votes > 0
    mean_effect = np.zeros_like(value_sum)
    mean_effect[available] = value_sum[available] / gene_votes[available]
    return CD4Prepared(
        mode=mode,
        mean_effect=mean_effect,
        gene_votes=gene_votes,
        target_condition_count=condition_count,
        condition_details=details,
    )


def blend_cd4(
    baseline: np.ndarray,
    cd4: CD4Prepared,
    policy: str,
) -> np.ndarray:
    """Apply the fixed 50% blend with an explicit missing-gene policy."""
    if policy not in MISSING_GENE_POLICIES:
        raise ValueError(f"unknown missing-gene policy {policy!r}")
    candidate = baseline.copy()
    target_covered = cd4.target_condition_count > 0
    if policy == "zero_vote":
        candidate[target_covered] = (1.0 - BLEND_FRACTION) * baseline[
            target_covered
        ] + BLEND_FRACTION * cd4.mean_effect[target_covered]
    else:
        available = cd4.gene_votes > 0
        candidate[available] = (1.0 - BLEND_FRACTION) * baseline[
            available
        ] + BLEND_FRACTION * cd4.mean_effect[available]
    # This is an invariant, not a numerical expectation.
    candidate[~target_covered] = baseline[~target_covered]
    return candidate


@dataclass(frozen=True)
class WhitenerMetadata:
    path: Path
    metadata_path: Path
    metadata: dict[str, object]
    wcols: np.ndarray
    panel_order_sha256: str
    excluded_sorted_sha256: str


def audit_whitener(
    whitener_path: Path,
    official_genes: Sequence[str],
    panel: Sequence[str],
    source_banks: Sequence[Bank],
) -> WhitenerMetadata:
    """Derive experiment 040's wcols and require exact adjacent metadata."""
    whitener_path = whitener_path.expanduser().resolve()
    metadata_path = whitener_path.with_suffix(".json")
    if not whitener_path.is_file():
        raise FileNotFoundError(
            f"clean HCT whitener is absent: {whitener_path}; never fall back "
            "to experiments/040-union-gamma/out/whitener_union.npz"
        )
    if not metadata_path.is_file():
        raise FileNotFoundError(
            f"clean HCT whitener metadata is absent: {metadata_path}"
        )
    if tuple(bank.name for bank in source_banks) != SOURCE_NAMES:
        raise ValueError(
            f"source order must be exactly {SOURCE_NAMES}, got "
            f"{tuple(bank.name for bank in source_banks)}"
        )

    donor_union = set().union(*(set(bank.genes) for bank in source_banks))
    panel_set = set(panel)
    wcols = np.asarray(
        sorted(
            index
            for index, gene in enumerate(official_genes)
            if gene in donor_union and gene not in panel_set
        ),
        dtype=np.int64,
    )
    expected = {
        "axis": "union",
        "donors": list(SOURCE_NAMES),
        "n_cols": len(wcols),
        "wcols_sha256": hashlib.sha256(wcols.tobytes()).hexdigest(),
        "excluded_targets": sorted(panel_set),
    }
    actual = json.loads(metadata_path.read_text())
    if actual != expected:
        differing = sorted(
            key
            for key in set(actual) | set(expected)
            if actual.get(key) != expected.get(key)
        )
        raise ValueError(
            f"{metadata_path}: clean-HCT whitener metadata mismatch in {differing}"
        )
    return WhitenerMetadata(
        path=whitener_path,
        metadata_path=metadata_path,
        metadata=actual,
        wcols=wcols,
        panel_order_sha256=sha256_lines(panel),
        excluded_sorted_sha256=sha256_lines(expected["excluded_targets"]),
    )


@dataclass
class ExtractedWhitener:
    evals_path: Path
    evecs_path: Path
    size: int

    @classmethod
    def from_npz(
        cls,
        source: Path,
        scratch: Path,
        expected_size: int,
        *,
        max_rss_gib: float,
    ) -> ExtractedWhitener:
        """Stream the uncompressed .npy members out of npz for mmap access."""
        scratch.mkdir(parents=True, exist_ok=True)
        with zipfile.ZipFile(source) as archive:
            members = {item.filename: item for item in archive.infolist()}
            required = {"evals.npy", "evecs.npy"}
            if set(members) != required:
                raise ValueError(
                    f"{source}: expected exactly {sorted(required)}, "
                    f"found {sorted(members)}"
                )
            for member in sorted(required):
                info = members[member]
                if info.is_dir() or Path(member).name != member:
                    raise ValueError(f"{source}: unsafe npz member {member!r}")
                destination = scratch / member
                with archive.open(info) as reader, destination.open("wb") as writer:
                    shutil.copyfileobj(reader, writer, length=8 << 20)
                enforce_memory_limit(max_rss_gib, f"extracting {member}")

        evals = np.load(scratch / "evals.npy", mmap_mode="r")
        evecs = np.load(scratch / "evecs.npy", mmap_mode="r")
        if not isinstance(evals, np.memmap) or not isinstance(evecs, np.memmap):
            raise TypeError("extracted whitener arrays are not memory mapped")
        if evals.shape != (expected_size,) or evecs.shape != (
            expected_size,
            expected_size,
        ):
            raise ValueError(
                f"{source}: operator shapes {evals.shape}, {evecs.shape} "
                f"!= {(expected_size,), (expected_size, expected_size)}"
            )
        if evals.dtype != np.float64 or evecs.dtype != np.float32:
            raise TypeError(
                f"{source}: expected float64 evals and float32 evecs, got "
                f"{evals.dtype} and {evecs.dtype}"
            )
        if not np.isfinite(evals).all() or np.any(evals <= 0):
            raise ValueError(f"{source}: eigenvalues must be finite and positive")
        return cls(
            evals_path=scratch / "evals.npy",
            evecs_path=scratch / "evecs.npy",
            size=expected_size,
        )

    def linear_apply(
        self,
        matrix: np.ndarray,
        *,
        feature_block: int,
        max_rss_gib: float,
    ) -> np.ndarray:
        """Apply V diag(lambda^-gamma/2) V^T without a full float64 V copy."""
        matrix = np.asarray(matrix, dtype=np.float64)
        if matrix.ndim != 2 or matrix.shape[1] != self.size:
            raise ValueError(
                f"whitener input shape {matrix.shape} does not end in {self.size}"
            )
        evals = np.load(self.evals_path, mmap_mode="r")
        evecs = np.load(self.evecs_path, mmap_mode="r")
        weights = np.asarray(evals, dtype=np.float64) ** (-WHITEN_GAMMA / 2.0)
        latent = np.empty_like(matrix)
        for lo in range(0, self.size, feature_block):
            hi = min(lo + feature_block, self.size)
            eigenvectors = np.asarray(evecs[:, lo:hi], dtype=np.float64)
            latent[:, lo:hi] = matrix @ eigenvectors
            enforce_memory_limit(max_rss_gib, "union whitener forward projection")
        latent *= weights

        output = np.empty_like(matrix)
        for lo in range(0, self.size, feature_block):
            hi = min(lo + feature_block, self.size)
            eigenvectors = np.asarray(evecs[lo:hi, :], dtype=np.float64)
            output[:, lo:hi] = latent @ eigenvectors.T
            enforce_memory_limit(max_rss_gib, "union whitener reverse projection")
        return output


def apply_whitener(
    matrix: np.ndarray,
    operator: ExtractedWhitener,
    wcols: np.ndarray,
    *,
    renormalization_policy: str,
    baseline_scale: float | None,
    feature_block: int,
    max_rss_gib: float,
) -> tuple[np.ndarray, dict[str, float | str]]:
    """Apply the fixed operator with production or diagnostic renormalization."""
    if renormalization_policy not in RENORMALIZATION_POLICIES:
        raise ValueError(f"unknown renormalization policy {renormalization_policy!r}")
    output = np.asarray(matrix, dtype=np.float64).copy()
    before = float(np.median(np.linalg.norm(output[:, wcols], axis=1))) or 1.0
    transformed = operator.linear_apply(
        output[:, wcols],
        feature_block=feature_block,
        max_rss_gib=max_rss_gib,
    )
    after = float(np.median(np.linalg.norm(transformed, axis=1))) or 1.0
    candidate_self_scale = before / after
    if renormalization_policy == "production_self":
        scale = candidate_self_scale
    else:
        if baseline_scale is None:
            raise ValueError("baseline_fixed requires the baseline whitening scale")
        scale = baseline_scale
    output[:, wcols] = transformed * scale
    return output, {
        "policy": renormalization_policy,
        "before_operator_median_norm": before,
        "after_operator_median_norm": after,
        "candidate_self_scale": candidate_self_scale,
        "applied_scale": scale,
    }


def _sum_rows(
    matrix: object,
    rows: np.ndarray,
    n_genes: int,
    *,
    chunk_rows: int,
    max_rss_gib: float,
    stage: str,
) -> np.ndarray:
    total = np.zeros(n_genes, dtype=np.float64)
    for lo in range(0, len(rows), chunk_rows):
        hi = min(lo + chunk_rows, len(rows))
        block = matrix[rows[lo:hi]]
        total += np.asarray(block.sum(axis=0)).ravel()
        enforce_memory_limit(max_rss_gib, stage)
    return total


@dataclass
class HCTTruth:
    targets: tuple[str, ...]
    truth: np.ndarray
    control_profile: np.ndarray
    cells_per_target: dict[str, int]
    control_cells: int


def load_hct_truth(
    h5ad_path: Path,
    official_genes: Sequence[str],
    expected_targets: Sequence[str],
    *,
    chunk_rows: int,
    max_rss_gib: float,
) -> HCTTruth:
    """Compute real HCT log-pseudobulk effects without loading the count matrix."""
    h5ad_path = h5ad_path.expanduser().resolve()
    dataset = ad.read_h5ad(h5ad_path, backed="r")
    try:
        observed_genes = tuple(map(str, dataset.var_names))
        if observed_genes != tuple(official_genes):
            raise ValueError(f"{h5ad_path}: gene axis differs from gene_names.csv")
        if PERT_COL not in dataset.obs:
            raise ValueError(f"{h5ad_path}: obs[{PERT_COL!r}] is missing")
        labels = dataset.obs[PERT_COL].astype(str).to_numpy()
        actual_targets = tuple(sorted(set(labels) - {CONTROL}))
        if actual_targets != tuple(sorted(expected_targets)):
            raise ValueError(
                f"{h5ad_path}: target set differs from its metadata sidecar"
            )
        targets = actual_targets
        control_rows = np.flatnonzero(labels == CONTROL)
        if not len(control_rows):
            raise ValueError(f"{h5ad_path}: no {CONTROL!r} cells")
        control_sum = _sum_rows(
            dataset.X,
            control_rows,
            len(official_genes),
            chunk_rows=chunk_rows,
            max_rss_gib=max_rss_gib,
            stage="HCT control pseudobulk",
        )
        if control_sum.sum() <= 0:
            raise ValueError(f"{h5ad_path}: zero HCT control library")
        control_profile = np.log1p(TARGET_SUM * control_sum / control_sum.sum())

        truth = np.empty((len(targets), len(official_genes)), dtype=np.float64)
        cells_per_target: dict[str, int] = {}
        for index, target in enumerate(targets):
            rows = np.flatnonzero(labels == target)
            cells_per_target[target] = len(rows)
            pert_sum = _sum_rows(
                dataset.X,
                rows,
                len(official_genes),
                chunk_rows=chunk_rows,
                max_rss_gib=max_rss_gib,
                stage=f"HCT truth {target}",
            )
            if pert_sum.sum() <= 0:
                raise ValueError(f"{h5ad_path}: zero library for {target}")
            truth[index] = (
                np.log1p(TARGET_SUM * pert_sum / pert_sum.sum()) - control_profile
            )
        return HCTTruth(
            targets=targets,
            truth=truth,
            control_profile=control_profile,
            cells_per_target=cells_per_target,
            control_cells=len(control_rows),
        )
    finally:
        dataset.file.close()


def pds_per_target(
    prediction: np.ndarray,
    truth: np.ndarray,
    score_columns: np.ndarray,
) -> np.ndarray:
    """Return official rank-style PDS per query against the full HCT gallery."""
    predicted = np.asarray(prediction[:, score_columns], dtype=np.float64)
    observed = np.asarray(truth[:, score_columns], dtype=np.float64)
    predicted /= np.maximum(np.linalg.norm(predicted, axis=1, keepdims=True), 1e-12)
    observed /= np.maximum(np.linalg.norm(observed, axis=1, keepdims=True), 1e-12)
    similarities = predicted @ observed.T
    own = np.diag(similarities)[:, None]
    n = similarities.shape[0]
    if similarities.shape != (n, n) or n < 2:
        raise ValueError("PDS requires equal prediction/truth panels with n >= 2")
    ranks = (similarities > own).sum(axis=1)
    ranks = ranks + 0.5 * ((similarities == own).sum(axis=1) - 1)
    return 1.0 - ranks / (n - 1)


def summarize_pds(
    per_target: np.ndarray,
    target_condition_count: np.ndarray,
) -> dict[str, dict[str, float | int | None]]:
    subsets = {
        "all_hct_targets": np.ones(len(per_target), dtype=bool),
        "cd4_covered_targets": target_condition_count > 0,
        "shared_all_three_condition_targets": target_condition_count == len(CD4_NAMES),
    }
    output: dict[str, dict[str, float | int | None]] = {}
    for name, mask in subsets.items():
        output[name] = {
            "n_targets": int(mask.sum()),
            "pds": float(np.mean(per_target[mask])) if mask.any() else None,
        }
    return output


def _write_rows(path: Path, rows: Sequence[dict[str, object]]) -> None:
    if not rows:
        return
    fieldnames = sorted({key for row in rows for key in row})
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def run(args: argparse.Namespace) -> dict[str, object]:
    panel = read_csv_column(Path(args.panel_csv), "target_gene")
    official = read_csv_column(Path(args.genes_csv), "gene_name")
    panel_set = set(panel)
    score_columns = np.asarray(
        [index for index, gene in enumerate(official) if gene not in panel_set],
        dtype=np.int64,
    )

    hct_meta_path = Path(args.hct_meta)
    hct_meta = json.loads(hct_meta_path.read_text())
    expected_hct_targets = tuple(map(str, hct_meta.get("perts", ())))
    if len(expected_hct_targets) < args.min_hct_targets:
        raise ValueError(
            f"{hct_meta_path}: only {len(expected_hct_targets)} HCT targets; "
            f"require {args.min_hct_targets}"
        )
    if not set(expected_hct_targets) <= panel_set:
        raise ValueError("HCT sidecar contains targets outside the challenge panel")

    truth = load_hct_truth(
        Path(args.hct_real),
        official,
        expected_hct_targets,
        chunk_rows=args.hct_chunk_rows,
        max_rss_gib=args.max_rss_gib,
    )
    targets = truth.targets
    source_banks = (
        Bank.load("k562_gwps", Path(args.k562)),
        Bank.load("orion_hek293t", Path(args.hek293t)),
    )
    prepared_sources = [
        prepare_main_source(
            bank,
            targets,
            official,
            chunk_rows=args.source_chunk_rows,
            max_rss_gib=args.max_rss_gib,
        )
        for bank in source_banks
    ]
    baseline_pre, source_votes = build_current_baseline(prepared_sources)

    audited = audit_whitener(
        Path(args.whitener),
        official,
        panel,
        source_banks,
    )
    if len(audited.wcols) < args.min_whitener_cols:
        raise ValueError(
            f"only {len(audited.wcols)} clean-HCT union columns; "
            f"require {args.min_whitener_cols}"
        )

    cd4_root = Path(args.cd4_root)
    cd4_banks = tuple(
        Bank.load(name, cd4_root / name, require_cd4_scale=True) for name in CD4_NAMES
    )
    cd4_by_mode = {
        mode: prepare_cd4_mean(
            cd4_banks,
            targets,
            official,
            panel_set,
            mode,
            chunk_rows=args.source_chunk_rows,
            max_rss_gib=args.max_rss_gib,
        )
        for mode in CD4_MODES
    }
    coverage_reference = cd4_by_mode["raw"].target_condition_count
    if any(
        not np.array_equal(prepared.target_condition_count, coverage_reference)
        for prepared in cd4_by_mode.values()
    ):
        raise AssertionError("CD4 target coverage changed between modes")
    if int((coverage_reference > 0).sum()) < args.min_cd4_covered_targets:
        raise ValueError(
            f"only {int((coverage_reference > 0).sum())} HCT targets have CD4; "
            f"require {args.min_cd4_covered_targets}"
        )

    output_path = Path(args.output)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    scratch_parent = Path(args.scratch_parent)
    scratch_parent.mkdir(parents=True, exist_ok=True)
    rows: list[dict[str, object]] = []
    per_target_rows: list[dict[str, object]] = []

    with tempfile.TemporaryDirectory(
        prefix="vcc054-whitener-", dir=scratch_parent
    ) as scratch_name:
        operator = ExtractedWhitener.from_npz(
            audited.path,
            Path(scratch_name),
            len(audited.wcols),
            max_rss_gib=args.max_rss_gib,
        )
        baseline_post, baseline_whitening = apply_whitener(
            baseline_pre,
            operator,
            audited.wcols,
            renormalization_policy="production_self",
            baseline_scale=None,
            feature_block=args.whitening_feature_block,
            max_rss_gib=args.max_rss_gib,
        )
        baseline_scale = float(baseline_whitening["applied_scale"])
        baseline_per_target = pds_per_target(baseline_post, truth.truth, score_columns)
        baseline_summary = summarize_pds(baseline_per_target, coverage_reference)
        rows.append(
            {
                "candidate": "baseline",
                "cd4_mode": "",
                "missing_gene_policy": "",
                "whitening_renormalization_policy": "production_self",
                "whitening_candidate_self_scale": baseline_whitening[
                    "candidate_self_scale"
                ],
                "whitening_applied_scale": baseline_whitening["applied_scale"],
                **{
                    f"{subset}_n": values["n_targets"]
                    for subset, values in baseline_summary.items()
                },
                **{
                    f"{subset}_pds": values["pds"]
                    for subset, values in baseline_summary.items()
                },
            }
        )

        for index, target in enumerate(targets):
            per_target_rows.append(
                {
                    "target_gene": target,
                    "candidate": "baseline",
                    "cd4_condition_count": int(coverage_reference[index]),
                    "whitening_renormalization_policy": "production_self",
                    "pds": float(baseline_per_target[index]),
                }
            )

        candidate_results: list[dict[str, object]] = []
        absent = coverage_reference == 0
        for mode in CD4_MODES:
            prepared = cd4_by_mode[mode]
            for policy in MISSING_GENE_POLICIES:
                candidate_pre = blend_cd4(baseline_pre, prepared, policy)
                pre_absent_diff = (
                    float(np.max(np.abs(candidate_pre[absent] - baseline_pre[absent])))
                    if absent.any()
                    else 0.0
                )
                if pre_absent_diff != 0.0:
                    raise AssertionError(
                        f"{mode}:{policy}: absent targets changed before whitening"
                    )
                for renormalization_policy in RENORMALIZATION_POLICIES:
                    candidate_name = (
                        f"cd4_mean3_f0.5:{mode}:{policy}:{renormalization_policy}"
                    )
                    candidate_post, whitening = apply_whitener(
                        candidate_pre,
                        operator,
                        audited.wcols,
                        renormalization_policy=renormalization_policy,
                        baseline_scale=baseline_scale,
                        feature_block=args.whitening_feature_block,
                        max_rss_gib=args.max_rss_gib,
                    )
                    post_absent_diff = (
                        float(
                            np.max(
                                np.abs(candidate_post[absent] - baseline_post[absent])
                            )
                        )
                        if absent.any()
                        else 0.0
                    )
                    if (
                        renormalization_policy == "baseline_fixed"
                        and post_absent_diff != 0.0
                    ):
                        raise AssertionError(
                            f"{candidate_name}: baseline_fixed changed absent targets"
                        )
                    per_target = pds_per_target(
                        candidate_post, truth.truth, score_columns
                    )
                    summary = summarize_pds(per_target, coverage_reference)
                    result = {
                        "candidate": candidate_name,
                        "cd4_mode": mode,
                        "missing_gene_policy": policy,
                        "whitening_renormalization_policy": (renormalization_policy),
                        "fraction": BLEND_FRACTION,
                        "is_primary": (
                            mode == "gene_norm"
                            and policy == "preserve_baseline"
                            and renormalization_policy == "production_self"
                        ),
                        "whitening": whitening,
                        "pds": summary,
                        "delta_vs_baseline": {
                            subset: (
                                None
                                if values["pds"] is None
                                else float(values["pds"])
                                - float(baseline_summary[subset]["pds"])
                            )
                            for subset, values in summary.items()
                        },
                        "pre_whiten_absent_target_max_abs_diff": pre_absent_diff,
                        "post_whiten_absent_target_max_abs_diff": post_absent_diff,
                        "post_whiten_absent_target_exact_equal": (
                            post_absent_diff == 0.0
                        ),
                    }
                    candidate_results.append(result)
                    rows.append(
                        {
                            "candidate": candidate_name,
                            "cd4_mode": mode,
                            "missing_gene_policy": policy,
                            "whitening_renormalization_policy": (
                                renormalization_policy
                            ),
                            "whitening_candidate_self_scale": whitening[
                                "candidate_self_scale"
                            ],
                            "whitening_applied_scale": whitening["applied_scale"],
                            **{
                                f"{subset}_n": values["n_targets"]
                                for subset, values in summary.items()
                            },
                            **{
                                f"{subset}_pds": values["pds"]
                                for subset, values in summary.items()
                            },
                            **{
                                f"{subset}_delta": result["delta_vs_baseline"][subset]
                                for subset in summary
                            },
                        },
                    )
                    for index, target in enumerate(targets):
                        per_target_rows.append(
                            {
                                "target_gene": target,
                                "candidate": candidate_name,
                                "cd4_condition_count": int(coverage_reference[index]),
                                "whitening_renormalization_policy": (
                                    renormalization_policy
                                ),
                                "pds": float(per_target[index]),
                            }
                        )
                    enforce_memory_limit(
                        args.max_rss_gib, f"candidate {candidate_name}"
                    )

    target_any = coverage_reference > 0
    target_all = coverage_reference == len(CD4_NAMES)
    primary_cd4 = cd4_by_mode["gene_norm"]
    result = {
        "experiment": "054-cd4-union-hct",
        "status": "completed",
        "protocol": {
            "metric": "effect-level PDS against real HCT log-pseudobulk truth",
            "count_emission": False,
            "cell_eval2": False,
            "submission": False,
            "source_donors": list(SOURCE_NAMES),
            "held_out_context": "Orion HCT116",
            "cd4_conditions": list(CD4_NAMES),
            "cd4_fraction": BLEND_FRACTION,
            "primary_candidate": (
                "cd4_mean3_f0.5:gene_norm:preserve_baseline:production_self"
            ),
            "selection_on_hct": False,
            "whitening_gamma": WHITEN_GAMMA,
            "whitening_renormalization_policies": {
                "production_self": (
                    "production preserve behavior: each candidate uses its own "
                    "before_operator / after_operator median-norm scalar"
                ),
                "baseline_fixed": (
                    "diagnostic isolation: every candidate reuses the baseline "
                    "before_operator / after_operator scalar"
                ),
            },
            "pds_subset_semantics": (
                "query-row subset averaged against the full HCT target gallery"
            ),
        },
        "inputs": {
            "hct_real": str(Path(args.hct_real).resolve()),
            "hct_meta": str(hct_meta_path.resolve()),
            "k562": str(Path(args.k562).resolve()),
            "hek293t": str(Path(args.hek293t).resolve()),
            "cd4_root": str(cd4_root.resolve()),
            "panel_csv": str(Path(args.panel_csv).resolve()),
            "genes_csv": str(Path(args.genes_csv).resolve()),
            "whitener": str(audited.path),
            "whitener_metadata": str(audited.metadata_path),
        },
        "axis_audit": {
            "official_genes": len(official),
            "challenge_panel_targets": len(panel),
            "score_genes_after_panel_exclusion": len(score_columns),
            "wcols": len(audited.wcols),
            "wcols_sha256": audited.metadata["wcols_sha256"],
            "panel_order_sha256": audited.panel_order_sha256,
            "excluded_sorted_sha256": audited.excluded_sorted_sha256,
            "whitener_metadata": audited.metadata,
            "whitener_npz_sha256": sha256_file(audited.path),
        },
        "hct_truth": {
            "targets": len(targets),
            "target_order_sha256": sha256_lines(targets),
            "control_cells": truth.control_cells,
            "cells_per_target_min": min(truth.cells_per_target.values()),
            "cells_per_target_median": float(
                np.median(list(truth.cells_per_target.values()))
            ),
            "cells_per_target_max": max(truth.cells_per_target.values()),
        },
        "baseline": {
            "source_vote_counts": {
                str(votes): int(np.sum(source_votes == votes))
                for votes in sorted(set(source_votes.tolist()))
            },
            "source_details": [
                {
                    "name": source.name,
                    "raw_median_centered_norm": source.raw_median_norm,
                    "mapped_genes": source.mapped_gene_count,
                    "duplicate_symbols": list(source.duplicate_symbols),
                }
                for source in prepared_sources
            ],
            "whitening": baseline_whitening,
            "pds": baseline_summary,
        },
        "cd4_coverage": {
            "any_condition_targets": int(target_any.sum()),
            "all_three_condition_targets": int(target_all.sum()),
            "absent_targets": int((~target_any).sum()),
            "absent_target_symbols": [
                target
                for target, covered in zip(targets, target_any, strict=True)
                if not covered
            ],
            "primary_mode_target_gene_entries_by_vote_count": {
                "zero": int(np.sum(primary_cd4.gene_votes == 0)),
                "one": int(np.sum(primary_cd4.gene_votes == 1)),
                "two": int(np.sum(primary_cd4.gene_votes == 2)),
                "three": int(np.sum(primary_cd4.gene_votes == 3)),
            },
            "condition_details_by_mode": {
                mode: prepared.condition_details
                for mode, prepared in cd4_by_mode.items()
            },
        },
        "missing_gene_control": {
            "preserve_baseline": (
                "blend only target-gene coordinates with at least one CD4 vote; "
                "leave the baseline unchanged where CD4 has no measured gene"
            ),
            "zero_vote": (
                "map missing CD4 genes to zero, then apply the 50% convex blend "
                "globally for every CD4-covered target"
            ),
        },
        "candidates": candidate_results,
        "peak_rss_gib": peak_rss_gib(),
    }
    output_path.write_text(json.dumps(result, indent=2))
    _write_rows(output_path.with_name("summary.csv"), rows)
    _write_rows(output_path.with_name("per_target_pds.csv"), per_target_rows)
    return result


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--hct-real",
        default=ROOT / "data/harness/orion-hct116.h5ad",
        type=Path,
    )
    parser.add_argument(
        "--hct-meta",
        default=ROOT / "data/harness/orion-hct116.h5ad.meta.json",
        type=Path,
    )
    parser.add_argument("--k562", default=ROOT / "data/donors/k562_gwps", type=Path)
    parser.add_argument(
        "--hek293t", default=ROOT / "data/donors/orion_hek293t", type=Path
    )
    parser.add_argument(
        "--cd4-root",
        default=ROOT / "data/external/cd4t/cloud-compact/cd4_compact",
        type=Path,
    )
    parser.add_argument(
        "--panel-csv",
        default=ROOT / "data/challenge/pert_counts.csv",
        type=Path,
    )
    parser.add_argument(
        "--genes-csv",
        default=ROOT / "data/challenge/gene_names.csv",
        type=Path,
    )
    parser.add_argument(
        "--whitener",
        default=(ROOT / "experiments/040-union-gamma/out_hct_clean/whitener_union.npz"),
        type=Path,
    )
    parser.add_argument("--output", default=HERE / "results.json", type=Path)
    parser.add_argument("--scratch-parent", default=HERE / "scratch", type=Path)
    parser.add_argument("--source-chunk-rows", type=int, default=64)
    parser.add_argument("--hct-chunk-rows", type=int, default=128)
    parser.add_argument("--whitening-feature-block", type=int, default=256)
    parser.add_argument("--max-rss-gib", type=float, default=7.5)
    parser.add_argument("--min-hct-targets", type=int, default=100)
    parser.add_argument("--min-cd4-covered-targets", type=int, default=100)
    parser.add_argument("--min-whitener-cols", type=int, default=17_000)
    parser.add_argument(
        "--guarded",
        action="store_true",
        help="acknowledge execution under scripts/run_with_memory_guard.py",
    )
    return parser


def validate_args(args: argparse.Namespace) -> None:
    if not args.guarded:
        raise SystemExit(
            "Refusing an unguarded real-data run. Use the exact cloud command "
            "in experiments/054-cd4-union-hct/README.md."
        )
    for name in (
        "source_chunk_rows",
        "hct_chunk_rows",
        "whitening_feature_block",
    ):
        if getattr(args, name) <= 0:
            raise SystemExit(f"--{name.replace('_', '-')} must be positive")
    if args.max_rss_gib <= 0:
        raise SystemExit("--max-rss-gib must be positive")


def main() -> int:
    args = build_parser().parse_args()
    validate_args(args)
    result = run(args)
    primary = next(row for row in result["candidates"] if row["is_primary"])
    baseline = result["baseline"]["pds"]["all_hct_targets"]["pds"]
    candidate = primary["pds"]["all_hct_targets"]["pds"]
    print(
        f"baseline={baseline:.6f} primary={candidate:.6f} "
        f"delta={candidate - baseline:+.6f} "
        f"peak_rss={result['peak_rss_gib']:.3f} GiB"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
