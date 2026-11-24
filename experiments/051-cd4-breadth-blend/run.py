"""051 — leak-free CD4 breadth-adaptive blending under donor-held-out PDS.

The evaluator tunes CD4 blend rules on held-out HCT116 and HEK293T, freezes
the selected configurations, and evaluates K562 exactly once as the
lineage-analogue holdout. Adaptive weights are functions only of metadata in
the compact CD4 donor directories; K562 comparison tables and challenge
outcomes are never read.

This is an effect-level PDS experiment. It does not run cell-eval2, emit
single-cell counts, edit production code, or submit anything.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import resource
import sys
import tempfile
from collections.abc import Iterable, Sequence
from dataclasses import asdict, dataclass
from pathlib import Path

# Bound BLAS concurrency before NumPy is imported. The outer process-tree guard
# is still mandatory for real-data execution.
os.environ.setdefault("OPENBLAS_NUM_THREADS", "2")
os.environ.setdefault("VECLIB_MAXIMUM_THREADS", "2")
os.environ.setdefault("OMP_NUM_THREADS", "2")
os.environ.setdefault("MKL_NUM_THREADS", "2")

import numpy as np
import pandas as pd
from scipy.special import expit

ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
MAIN_PATHS = {
    "k562_gwps": ROOT / "data/donors/k562_gwps",
    "orion_hct116": ROOT / "data/donors/orion_hct116",
    "orion_hek293t": ROOT / "data/donors/orion_hek293t",
}
MAIN_NAMES = tuple(MAIN_PATHS)
DEVELOPMENT_HELD = ("orion_hct116", "orion_hek293t")
HOLDOUT = "k562_gwps"
CD4_NAMES = ("cd4_rest", "cd4_stim8hr", "cd4_stim48hr")
CD4_SOURCE_NAMES = ("cd4_stim8hr", "cd4_mean3")
MODES = ("raw", "gene_norm")
WHITENER_SCHEMA = "vcc2026-exp051-fold-union-whitener-v1"


def peak_rss_gib() -> float:
    """Return peak process RSS in GiB on macOS or Linux."""
    rss = float(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
    return rss / 1024**3 if sys.platform == "darwin" else rss / 1024**2


def enforce_memory_limit(limit_gib: float, stage: str) -> None:
    if limit_gib > 0 and peak_rss_gib() >= limit_gib:
        raise MemoryError(
            f"{stage}: peak RSS {peak_rss_gib():.2f} GiB reached the "
            f"{limit_gib:.2f} GiB internal stop"
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


def stable_seed(parts: Iterable[object], base: int) -> int:
    payload = "|".join(map(str, parts)) + f"|{base}"
    return int.from_bytes(hashlib.sha256(payload.encode()).digest()[:8], "little")


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

    @classmethod
    def load(cls, name: str, path: Path, *, cd4: bool = False) -> Bank:
        path = path.expanduser().resolve()
        genes = read_lines(path / "genes.txt", require_unique=False)
        perts = read_lines(path / "perts.txt")
        delta_path = path / "delta.npy"
        delta = np.load(delta_path, mmap_mode="r")
        if not isinstance(delta, np.memmap):
            raise TypeError(f"{name}: delta.npy did not open as a memory map")
        if delta.shape != (len(perts), len(genes)):
            raise ValueError(
                f"{name}: delta shape {delta.shape} != "
                f"labels {(len(perts), len(genes))}"
            )
        if not np.issubdtype(delta.dtype, np.floating):
            raise TypeError(f"{name}: delta dtype {delta.dtype} is not floating")
        b_ctrl = np.load(path / "b_ctrl.npy", mmap_mode="r")
        if b_ctrl.shape != (len(genes),):
            raise ValueError(
                f"{name}: b_ctrl shape {b_ctrl.shape} != gene axis {(len(genes),)}"
            )
        gene_pos: dict[str, int] = {}
        for index, gene in enumerate(genes):
            previous = gene_pos.get(gene)
            if previous is None or float(b_ctrl[index]) > float(b_ctrl[previous]):
                gene_pos[gene] = index
        if cd4:
            meta = json.loads((path / "meta.json").read_text())
            definition = str(meta.get("effect_definition", ""))
            if "ln(2)" not in definition or "log2 fold change" not in definition:
                raise ValueError(
                    f"{name}: unrecognized CD4 effect scale: {definition!r}"
                )
            if not (path / "obs.csv").exists():
                raise FileNotFoundError(path / "obs.csv")
        return cls(
            name=name,
            path=path,
            genes=genes,
            perts=perts,
            gene_pos=gene_pos,
            pert_pos={pert: i for i, pert in enumerate(perts)},
            delta_path=delta_path,
            shape=tuple(map(int, delta.shape)),
        )


@dataclass(frozen=True)
class PreparedBank:
    name: str
    raw_path: Path
    gene_norm_path: Path
    available: np.ndarray
    shape: tuple[int, int]

    def load(self, mode: str) -> np.memmap:
        if mode not in MODES:
            raise ValueError(f"unknown mode {mode!r}")
        path = self.raw_path if mode == "raw" else self.gene_norm_path
        matrix = np.load(path, mmap_mode="r")
        if not isinstance(matrix, np.memmap):
            raise TypeError(f"{path}: prepared effect matrix is not memory mapped")
        if matrix.shape != self.shape:
            raise ValueError(
                f"{path}: shape changed from {self.shape} to {matrix.shape}"
            )
        return matrix


@dataclass(frozen=True)
class ConfidenceData:
    source: str
    breadth_log: np.ndarray
    quality: np.ndarray
    provenance: dict[str, object]


@dataclass(frozen=True)
class Candidate:
    family: str
    mode: str
    cd4_source: str
    geometry: str
    fraction: float
    slope: float = 0.0
    qc_strength: float = 0.0

    @property
    def key(self) -> str:
        fields = [
            self.family,
            self.mode,
            self.cd4_source,
            self.geometry,
            f"f={self.fraction:.8g}",
        ]
        if self.family == "adaptive":
            fields.extend([f"slope={self.slope:.8g}", f"qc={self.qc_strength:.8g}"])
        return ":".join(fields)

    @classmethod
    def from_dict(cls, value: dict[str, object]) -> Candidate:
        return cls(
            family=str(value["family"]),
            mode=str(value["mode"]),
            cd4_source=str(value["cd4_source"]),
            geometry=str(value["geometry"]),
            fraction=float(value["fraction"]),
            slope=float(value.get("slope", 0.0)),
            qc_strength=float(value.get("qc_strength", 0.0)),
        )


@dataclass
class FoldComponents:
    held: str
    mode: str
    geometry: str
    truth: np.ndarray
    baseline: np.ndarray
    cd4: dict[str, np.ndarray]
    score_genes: tuple[str, ...]


def _stream_center_scale(
    bank: Bank,
    columns: np.ndarray,
    *,
    chunk_rows: int,
    max_rss_gib: float,
) -> tuple[np.ndarray, np.ndarray]:
    """Compute source-only center and mean-absolute scale in bounded blocks."""
    delta = np.load(bank.delta_path, mmap_mode="r")
    total = np.zeros(len(columns), dtype=np.float64)
    for lo in range(0, bank.shape[0], chunk_rows):
        hi = min(lo + chunk_rows, bank.shape[0])
        block = np.asarray(delta[lo:hi, columns], dtype=np.float32)
        total += block.sum(axis=0, dtype=np.float64)
        enforce_memory_limit(max_rss_gib, f"{bank.name} center")
    mean = total / bank.shape[0]

    absolute = np.zeros(len(columns), dtype=np.float64)
    for lo in range(0, bank.shape[0], chunk_rows):
        hi = min(lo + chunk_rows, bank.shape[0])
        block = np.asarray(delta[lo:hi, columns], dtype=np.float32)
        absolute += np.abs(block - mean).sum(axis=0, dtype=np.float64)
        enforce_memory_limit(max_rss_gib, f"{bank.name} scale")
    mean_abs = absolute / bank.shape[0]
    live = mean_abs[mean_abs > 0]
    floor = float(np.percentile(live, 25.0)) if live.size else 1.0
    return mean.astype(np.float32), np.maximum(mean_abs, floor).astype(np.float32)


def prepare_bank(
    bank: Bank,
    perturbations: Sequence[str],
    master_features: Sequence[str],
    scratch: Path,
    *,
    chunk_rows: int,
    max_rss_gib: float,
) -> PreparedBank:
    """Stream one donor into panel-row, official-union memory maps."""
    destination_columns: list[int] = []
    source_columns: list[int] = []
    for destination, gene in enumerate(master_features):
        source = bank.gene_pos.get(gene)
        if source is not None:
            destination_columns.append(destination)
            source_columns.append(source)
    dst = np.asarray(destination_columns, dtype=np.int64)
    src = np.asarray(source_columns, dtype=np.int64)
    if not len(src):
        raise ValueError(f"{bank.name}: no genes overlap the master feature axis")

    rows = np.asarray([bank.pert_pos[pert] for pert in perturbations], dtype=np.int64)
    mean, scale = _stream_center_scale(
        bank,
        src,
        chunk_rows=chunk_rows,
        max_rss_gib=max_rss_gib,
    )

    bank_dir = scratch / bank.name
    bank_dir.mkdir(parents=True)
    shape = (len(perturbations), len(master_features))
    raw = np.lib.format.open_memmap(
        bank_dir / "raw.npy", mode="w+", dtype=np.float32, shape=shape
    )
    gene_norm = np.lib.format.open_memmap(
        bank_dir / "gene_norm.npy", mode="w+", dtype=np.float32, shape=shape
    )
    raw[:] = 0
    gene_norm[:] = 0
    delta = np.load(bank.delta_path, mmap_mode="r")
    for lo in range(0, len(rows), chunk_rows):
        hi = min(lo + chunk_rows, len(rows))
        block = np.asarray(delta[np.ix_(rows[lo:hi], src)], dtype=np.float32)
        centered = block - mean
        raw[lo:hi, dst] = centered
        gene_norm[lo:hi, dst] = centered / scale
        enforce_memory_limit(max_rss_gib, f"{bank.name} panel extraction")

    before = float(np.median(np.linalg.norm(raw, axis=1))) or 1.0
    after = float(np.median(np.linalg.norm(gene_norm, axis=1))) or 1.0
    gene_norm *= np.float32(before / after)
    raw.flush()
    gene_norm.flush()
    del raw, gene_norm, delta

    available = np.zeros(len(master_features), dtype=bool)
    available[dst] = True
    return PreparedBank(
        name=bank.name,
        raw_path=bank_dir / "raw.npy",
        gene_norm_path=bank_dir / "gene_norm.npy",
        available=available,
        shape=shape,
    )


def _as_bool(value: object) -> float:
    if value is None or (isinstance(value, float) and np.isnan(value)):
        return float("nan")
    if isinstance(value, (bool, np.bool_)):
        return float(bool(value))
    text = str(value).strip().casefold()
    if text in {"1", "true", "t", "yes", "y"}:
        return 1.0
    if text in {"0", "false", "f", "no", "n"}:
        return 0.0
    return float("nan")


def _rank01(values: np.ndarray) -> np.ndarray:
    series = pd.Series(np.asarray(values, dtype=np.float64))
    finite = np.isfinite(series.to_numpy())
    if finite.sum() < 2:
        return np.full(len(series), 0.5, dtype=np.float64)
    ranked = series.rank(method="average", pct=True).to_numpy(dtype=np.float64)
    fill = float(np.nanmedian(ranked[finite]))
    return np.where(np.isfinite(ranked), ranked, fill)


def _condition_confidence(
    bank: Bank,
    perturbations: Sequence[str],
) -> tuple[np.ndarray, np.ndarray, dict[str, object]]:
    """Read only compact-CD4 breadth/QC fields for one condition."""
    frame = pd.read_csv(bank.path / "obs.csv")
    if "target_gene" not in frame:
        raise ValueError(f"{bank.name}/obs.csv has no target_gene column")
    if frame["target_gene"].astype(str).duplicated().any():
        raise ValueError(f"{bank.name}/obs.csv has duplicate target_gene rows")
    frame["target_gene"] = frame["target_gene"].astype(str)
    frame = frame.set_index("target_gene").reindex(list(perturbations))
    if frame.index.has_duplicates or frame.isna().all(axis=1).any():
        missing = frame.index[frame.isna().all(axis=1)].tolist()
        raise ValueError(f"{bank.name}: metadata missing targets {missing[:8]}")

    breadth_column = ""
    breadth = None
    for column in ("n_total_de_genes", "n_downstream"):
        if column in frame:
            values = pd.to_numeric(frame[column], errors="coerce").to_numpy(float)
            if np.isfinite(values).sum() >= max(3, len(values) // 2):
                breadth_column = column
                breadth = values
                break
    if breadth is None and {"n_up_genes", "n_down_genes"} <= set(frame.columns):
        up = pd.to_numeric(frame["n_up_genes"], errors="coerce").to_numpy(float)
        down = pd.to_numeric(frame["n_down_genes"], errors="coerce").to_numpy(float)
        breadth = up + down
        breadth_column = "n_up_genes+n_down_genes"
    if breadth is None:
        raise ValueError(f"{bank.name}: no usable CD4 effect-breadth field")
    breadth = np.maximum(np.nan_to_num(breadth, nan=0.0), 0.0)

    quality_parts: list[np.ndarray] = []
    quality_columns: list[str] = []
    if "ontarget_significant" in frame:
        quality_parts.append(
            np.asarray([_as_bool(value) for value in frame["ontarget_significant"]])
        )
        quality_columns.append("ontarget_significant")
    for column in (
        "low_target_gex",
        "distal_offtarget_flag",
        "neighboring_gene_KD",
        "offtarget_flag",
    ):
        if column in frame:
            bad = np.asarray([_as_bool(value) for value in frame[column]])
            quality_parts.append(1.0 - bad)
            quality_columns.append(f"not_{column}")
    for column in ("n_guides", "n_cells_target"):
        if column in frame:
            values = pd.to_numeric(frame[column], errors="coerce").to_numpy(float)
            quality_parts.append(_rank01(np.log1p(np.maximum(values, 0.0))))
            quality_columns.append(f"rank_log1p_{column}")

    if quality_parts:
        stacked = np.vstack(quality_parts)
        with np.errstate(invalid="ignore"):
            quality = np.nanmean(stacked, axis=0)
        quality = np.nan_to_num(quality, nan=0.5)
    else:
        quality = np.full(len(perturbations), 0.5, dtype=np.float64)

    return (
        np.log1p(breadth),
        np.clip(quality, 0.0, 1.0),
        {
            "condition": bank.name,
            "breadth_column": breadth_column,
            "quality_columns": quality_columns,
        },
    )


def load_confidence_data(
    cd4_banks: dict[str, Bank],
    perturbations: Sequence[str],
) -> dict[str, ConfidenceData]:
    per_condition = {
        name: _condition_confidence(bank, perturbations)
        for name, bank in cd4_banks.items()
    }
    stim_breadth, stim_quality, stim_provenance = per_condition["cd4_stim8hr"]
    all_breadth = np.vstack([per_condition[name][0] for name in CD4_NAMES])
    all_quality = np.vstack([per_condition[name][1] for name in CD4_NAMES])
    return {
        "cd4_stim8hr": ConfidenceData(
            source="cd4_stim8hr",
            breadth_log=stim_breadth,
            quality=stim_quality,
            provenance=stim_provenance,
        ),
        "cd4_mean3": ConfidenceData(
            source="cd4_mean3",
            breadth_log=np.mean(all_breadth, axis=0),
            quality=np.mean(all_quality, axis=0),
            provenance={
                "condition": "mean of Rest, Stim8hr, Stim48hr",
                "breadth_columns": {
                    name: per_condition[name][2]["breadth_column"] for name in CD4_NAMES
                },
                "quality_columns": {
                    name: per_condition[name][2]["quality_columns"]
                    for name in CD4_NAMES
                },
            },
        ),
    }


def robust_z(values: np.ndarray) -> np.ndarray:
    values = np.asarray(values, dtype=np.float64)
    median = float(np.median(values))
    q25, q75 = np.quantile(values, [0.25, 0.75])
    scale = float(q75 - q25)
    if scale <= 1e-12:
        scale = float(np.std(values))
    if scale <= 1e-12:
        return np.zeros_like(values)
    return (values - median) / scale


def _bounded_mean_scale(
    shape: np.ndarray,
    target_mean: float,
    lower: float,
    upper: float,
) -> np.ndarray:
    """Scale positive scores to an exact bounded mean by monotone bisection."""
    if not lower <= target_mean <= upper:
        raise ValueError(f"target mean {target_mean} must lie in [{lower}, {upper}]")
    shape = np.maximum(np.asarray(shape, dtype=np.float64), 1e-12)
    lo, hi = 0.0, max(1.0, upper / float(np.min(shape)))
    for _ in range(80):
        middle = 0.5 * (lo + hi)
        mean = float(np.clip(middle * shape, lower, upper).mean())
        if mean < target_mean:
            lo = middle
        else:
            hi = middle
    return np.clip(0.5 * (lo + hi) * shape, lower, upper)


def adaptive_weights(
    confidence: ConfidenceData,
    *,
    fraction: float,
    slope: float,
    qc_strength: float,
    lower: float,
    upper: float,
    breadth_override: np.ndarray | None = None,
) -> np.ndarray:
    """Smooth monotone CD4 weights with the same mean as a fixed comparator."""
    breadth = (
        confidence.breadth_log
        if breadth_override is None
        else np.asarray(breadth_override, dtype=np.float64)
    )
    if breadth.shape != confidence.breadth_log.shape:
        raise ValueError("breadth override shape mismatch")
    score = robust_z(breadth) + qc_strength * robust_z(confidence.quality)
    shape = expit(slope * score)
    return _bounded_mean_scale(shape, fraction, lower, upper)


def row_normalize(matrix: np.ndarray) -> np.ndarray:
    matrix = np.asarray(matrix, dtype=np.float32)
    norms = np.linalg.norm(matrix, axis=1, keepdims=True)
    return matrix / np.maximum(norms, np.float32(1e-12))


def pds(prediction: np.ndarray, truth: np.ndarray) -> float:
    """Exact VCC PDS rank statistic on a square perturbation panel."""
    scores = row_normalize(prediction) @ row_normalize(truth).T
    own = np.diag(scores)[:, None]
    n = scores.shape[1]
    if scores.shape != (n, n) or n < 2:
        raise ValueError("PDS needs equally sized prediction/truth panels with n >= 2")
    rank = (scores > own).sum(axis=1) + 0.5 * ((scores == own).sum(axis=1) - 1)
    return float(np.mean(1.0 - rank / (n - 1)))


def _scaled_to_reference(matrix: np.ndarray, reference_norm: float) -> np.ndarray:
    matrix = np.asarray(matrix, dtype=np.float32)
    norm = float(np.median(np.linalg.norm(matrix, axis=1))) or 1.0
    return matrix * np.float32(reference_norm / norm)


class AuditedWhitener:
    """Memory-mapped fold-specific production-union whitening operator."""

    def __init__(
        self,
        held: str,
        positions: np.ndarray,
        evals_path: Path,
        evecs_path: Path,
        gamma: float,
        max_rss_gib: float,
        block_rows: int,
    ):
        self.held = held
        self.positions = positions
        self.evals_path = evals_path
        self.evecs_path = evecs_path
        self.gamma = gamma
        self.max_rss_gib = max_rss_gib
        self.block_rows = block_rows

    @classmethod
    def load(
        cls,
        root: Path,
        held: str,
        score_genes: Sequence[str],
        expected_operator_genes: Sequence[str],
        expected_sources: Sequence[str],
        panel_digest: str,
        *,
        gamma: float,
        max_rss_gib: float,
        block_rows: int,
    ) -> AuditedWhitener:
        fold_dir = root / held
        meta = json.loads((fold_dir / "meta.json").read_text())
        operator_genes = read_lines(fold_dir / "genes.txt")
        expected_operator_genes = tuple(expected_operator_genes)
        checks = {
            "schema": meta.get("schema") == WHITENER_SCHEMA,
            "held_out_donor": meta.get("held_out_donor") == held,
            "fit_uses_held_truth": meta.get("fit_uses_held_truth") is False,
            "fit_uses_challenge_outcomes": (
                meta.get("fit_uses_challenge_outcomes") is False
            ),
            "axis_kind": meta.get("axis_kind") == "production_union",
            "source_donors": sorted(meta.get("source_donors", []))
            == sorted(expected_sources),
            "excluded_targets_sha256": (
                meta.get("excluded_targets_sha256") == panel_digest
            ),
            "operator_genes": operator_genes == expected_operator_genes,
            "operator_genes_sha256": (
                meta.get("operator_genes_sha256")
                == sha256_lines(expected_operator_genes)
            ),
        }
        failed = [name for name, passed in checks.items() if not passed]
        if failed:
            raise ValueError(f"{fold_dir}: whitener provenance audit failed: {failed}")

        evals_path = fold_dir / "evals.npy"
        evecs_path = fold_dir / "evecs.npy"
        evals = np.load(evals_path, mmap_mode="r")
        evecs = np.load(evecs_path, mmap_mode="r")
        n = len(operator_genes)
        if evals.shape != (n,) or evecs.shape != (n, n):
            raise ValueError(
                f"{fold_dir}: operator shapes {evals.shape}, {evecs.shape} "
                f"do not match {n} genes"
            )
        if not isinstance(evecs, np.memmap):
            raise TypeError(f"{evecs_path}: evecs must be a .npy memory map")
        score_pos = {gene: i for i, gene in enumerate(score_genes)}
        positions = np.asarray([score_pos[gene] for gene in operator_genes])
        return cls(
            held=held,
            positions=positions,
            evals_path=evals_path,
            evecs_path=evecs_path,
            gamma=gamma,
            max_rss_gib=max_rss_gib,
            block_rows=block_rows,
        )

    def apply(self, matrix: np.ndarray) -> np.ndarray:
        output = np.asarray(matrix, dtype=np.float32).copy()
        if self.gamma == 0 or not len(self.positions):
            return output
        evals = np.asarray(np.load(self.evals_path, mmap_mode="r"), dtype=np.float64)
        evecs = np.load(self.evecs_path, mmap_mode="r")
        floor = max(float(np.max(evals)) * 1e-12, np.finfo(np.float64).tiny)
        weights = np.maximum(evals, floor) ** (-self.gamma / 2.0)
        weights = weights.astype(np.float32)
        for lo in range(0, len(output), self.block_rows):
            hi = min(lo + self.block_rows, len(output))
            block = output[lo:hi, self.positions]
            latent = np.asarray(block @ evecs, dtype=np.float32)
            output[lo:hi, self.positions] = (latent * weights) @ evecs.T
            enforce_memory_limit(
                self.max_rss_gib, f"{self.held} union-whitener application"
            )
        return output


def _fold_components(
    held: str,
    mode: str,
    master_features: Sequence[str],
    banks: dict[str, Bank],
    prepared: dict[str, PreparedBank],
    whitener_root: Path | None,
    panel_digest: str,
    *,
    whitening_gamma: float,
    max_rss_gib: float,
    whitening_block_rows: int,
) -> dict[str, FoldComponents]:
    held_mask = prepared[held].available
    score_positions = np.flatnonzero(held_mask)
    score_genes = tuple(master_features[i] for i in score_positions)
    truth = np.asarray(prepared[held].load("raw")[:, score_positions], dtype=np.float32)
    source_names = [name for name in MAIN_NAMES if name != held]

    source_parts: dict[str, np.ndarray] = {}
    source_norms: dict[str, float] = {}
    for name in source_names:
        matrix = np.asarray(
            prepared[name].load(mode)[:, score_positions], dtype=np.float32
        )
        source_parts[name] = matrix
        source_norms[name] = float(np.median(np.linalg.norm(matrix, axis=1))) or 1.0
    reference_norm = float(np.median(list(source_norms.values()))) or 1.0
    baseline = sum(
        _scaled_to_reference(source_parts[name], reference_norm)
        for name in source_names
    ) / len(source_names)

    cd4_conditions = {
        name: _scaled_to_reference(
            np.asarray(prepared[name].load(mode)[:, score_positions], dtype=np.float32),
            reference_norm,
        )
        for name in CD4_NAMES
    }
    cd4 = {
        "cd4_stim8hr": cd4_conditions["cd4_stim8hr"],
        "cd4_mean3": sum(cd4_conditions.values()) / len(cd4_conditions),
    }
    geometries = {
        "none": FoldComponents(
            held=held,
            mode=mode,
            geometry="none",
            truth=truth,
            baseline=baseline,
            cd4=cd4,
            score_genes=score_genes,
        )
    }

    if whitener_root is not None:
        expected_operator_genes = tuple(
            gene
            for position, gene in zip(score_positions, score_genes, strict=True)
            if any(prepared[name].available[position] for name in source_names)
        )
        whitener = AuditedWhitener.load(
            whitener_root,
            held,
            score_genes,
            expected_operator_genes,
            source_names,
            panel_digest,
            gamma=whitening_gamma,
            max_rss_gib=max_rss_gib,
            block_rows=whitening_block_rows,
        )
        geometry_name = f"audited_union:g={whitening_gamma:g}"
        geometries[geometry_name] = FoldComponents(
            held=held,
            mode=mode,
            geometry=geometry_name,
            truth=truth,
            baseline=whitener.apply(baseline),
            cd4={name: whitener.apply(matrix) for name, matrix in cd4.items()},
            score_genes=score_genes,
        )
    enforce_memory_limit(max_rss_gib, f"{held} {mode} fold preparation")
    return geometries


def _candidate_prediction(
    candidate: Candidate,
    fold: FoldComponents,
    confidence: ConfidenceData,
    *,
    weight_lower: float,
    weight_upper: float,
    breadth_override: np.ndarray | None = None,
) -> tuple[np.ndarray, np.ndarray]:
    if candidate.family == "fixed":
        weights = np.full(len(fold.baseline), candidate.fraction, dtype=np.float64)
    elif candidate.family == "adaptive":
        weights = adaptive_weights(
            confidence,
            fraction=candidate.fraction,
            slope=candidate.slope,
            qc_strength=candidate.qc_strength,
            lower=weight_lower,
            upper=weight_upper,
            breadth_override=breadth_override,
        )
    else:
        raise ValueError(f"unknown candidate family {candidate.family!r}")
    prediction = (1.0 - weights[:, None]) * fold.baseline + weights[:, None] * fold.cd4[
        candidate.cd4_source
    ]
    return prediction.astype(np.float32), weights


def _candidate_grid(
    args: argparse.Namespace, geometries: Sequence[str]
) -> list[Candidate]:
    candidates: list[Candidate] = []
    for mode in args.modes:
        for source in args.cd4_sources:
            for geometry in geometries:
                for fraction in args.fractions:
                    candidates.append(
                        Candidate(
                            family="fixed",
                            mode=mode,
                            cd4_source=source,
                            geometry=geometry,
                            fraction=fraction,
                        )
                    )
                    for slope in args.slopes:
                        for qc_strength in args.qc_strengths:
                            candidates.append(
                                Candidate(
                                    family="adaptive",
                                    mode=mode,
                                    cd4_source=source,
                                    geometry=geometry,
                                    fraction=fraction,
                                    slope=slope,
                                    qc_strength=qc_strength,
                                )
                            )
    return candidates


def evaluate_development(
    args: argparse.Namespace,
    master_features: Sequence[str],
    banks: dict[str, Bank],
    prepared: dict[str, PreparedBank],
    confidence: dict[str, ConfidenceData],
    panel_digest: str,
) -> tuple[list[dict[str, object]], dict[tuple[str, str, str], FoldComponents]]:
    """Evaluate the full grid on HCT116 and HEK293T only."""
    fold_cache: dict[tuple[str, str, str], FoldComponents] = {}
    geometry_names: list[str] | None = None
    for held in DEVELOPMENT_HELD:
        for mode in args.modes:
            geometries = _fold_components(
                held,
                mode,
                master_features,
                banks,
                prepared,
                args.whitener_dir,
                panel_digest,
                whitening_gamma=args.whitening_gamma,
                max_rss_gib=args.max_rss_gib,
                whitening_block_rows=args.whitening_block_rows,
            )
            if geometry_names is None:
                geometry_names = list(geometries)
            elif list(geometries) != geometry_names:
                raise RuntimeError(
                    "geometry availability differs across development folds"
                )
            for geometry, fold in geometries.items():
                fold_cache[(held, mode, geometry)] = fold

    assert geometry_names is not None
    candidates = _candidate_grid(args, geometry_names)
    rows: list[dict[str, object]] = []
    for candidate in candidates:
        for held in DEVELOPMENT_HELD:
            fold = fold_cache[(held, candidate.mode, candidate.geometry)]
            baseline_score = pds(fold.baseline, fold.truth)
            prediction, weights = _candidate_prediction(
                candidate,
                fold,
                confidence[candidate.cd4_source],
                weight_lower=args.weight_lower,
                weight_upper=args.weight_upper,
            )
            score = pds(prediction, fold.truth)
            rows.append(
                {
                    "held": held,
                    "candidate": candidate.key,
                    **asdict(candidate),
                    "pds": score,
                    "baseline_pds": baseline_score,
                    "delta_vs_baseline": score - baseline_score,
                    "weight_min": float(weights.min()),
                    "weight_mean": float(weights.mean()),
                    "weight_max": float(weights.max()),
                }
            )
        enforce_memory_limit(args.max_rss_gib, f"candidate {candidate.key}")
    return rows, fold_cache


def _complexity(candidate: Candidate) -> tuple[object, ...]:
    """Deterministic preference among statistically tied development results."""
    return (
        0 if candidate.geometry == "none" else 1,
        0 if candidate.mode == "raw" else 1,
        0 if candidate.cd4_source == "cd4_stim8hr" else 1,
        abs(candidate.fraction - 1 / 3),
        candidate.slope,
        candidate.qc_strength,
        candidate.key,
    )


def summarize_and_select(
    rows: Sequence[dict[str, object]],
    family: str,
    *,
    tie_tolerance: float,
) -> tuple[Candidate, list[dict[str, object]]]:
    grouped: dict[str, list[dict[str, object]]] = {}
    for row in rows:
        if row["family"] == family:
            grouped.setdefault(str(row["candidate"]), []).append(dict(row))
    summaries: list[dict[str, object]] = []
    candidates: dict[str, Candidate] = {}
    for key, values in grouped.items():
        held = {str(value["held"]) for value in values}
        if held != set(DEVELOPMENT_HELD):
            raise RuntimeError(f"{key}: development folds are {held}")
        candidate = Candidate.from_dict(values[0])
        candidates[key] = candidate
        deltas = np.asarray(
            [float(value["delta_vs_baseline"]) for value in values],
            dtype=np.float64,
        )
        scores = np.asarray(
            [float(value["pds"]) for value in values],
            dtype=np.float64,
        )
        summaries.append(
            {
                "candidate": key,
                **asdict(candidate),
                "development_mean_pds": float(scores.mean()),
                "development_min_pds": float(scores.min()),
                "development_mean_delta": float(deltas.mean()),
                "development_min_delta": float(deltas.min()),
                "development_max_delta": float(deltas.max()),
                "development_fold_deltas": {
                    str(value["held"]): float(value["delta_vs_baseline"])
                    for value in values
                },
            }
        )
    if not summaries:
        raise RuntimeError(f"no development candidates for family {family}")
    best_mean = max(float(item["development_mean_pds"]) for item in summaries)
    near = [
        item
        for item in summaries
        if float(item["development_mean_pds"]) >= best_mean - tie_tolerance
    ]
    selected_summary = max(
        near,
        key=lambda item: float(item["development_min_pds"]),
    )
    best_min = float(selected_summary["development_min_pds"])
    equally_robust = [
        item
        for item in near
        if float(item["development_min_pds"]) >= best_min - tie_tolerance
    ]
    selected_summary = min(
        equally_robust,
        key=lambda item: _complexity(candidates[str(item["candidate"])]),
    )
    selected = candidates[str(selected_summary["candidate"])]
    ordered = sorted(
        summaries,
        key=lambda item: (
            -float(item["development_mean_pds"]),
            -float(item["development_min_pds"]),
            _complexity(candidates[str(item["candidate"])]),
        ),
    )
    return selected, ordered


def select_overall(
    selected_fixed: Candidate,
    fixed_summary: Sequence[dict[str, object]],
    selected_adaptive: Candidate,
    adaptive_summary: Sequence[dict[str, object]],
    *,
    tie_tolerance: float,
) -> Candidate:
    summary_by_key = {
        str(row["candidate"]): row for row in [*fixed_summary, *adaptive_summary]
    }
    fixed_mean = float(summary_by_key[selected_fixed.key]["development_mean_pds"])
    adaptive_mean = float(summary_by_key[selected_adaptive.key]["development_mean_pds"])
    if adaptive_mean > fixed_mean + tie_tolerance:
        return selected_adaptive
    return selected_fixed


def _shuffle_scores(
    candidate: Candidate,
    fold: FoldComponents,
    confidence: ConfidenceData,
    *,
    replicates: int,
    seed: int,
    weight_lower: float,
    weight_upper: float,
) -> list[float]:
    scores: list[float] = []
    for replicate in range(replicates):
        rng = np.random.default_rng(
            stable_seed(
                ("breadth-shuffle", candidate.key, replicate),
                base=seed,
            )
        )
        shuffled_breadth = confidence.breadth_log[
            rng.permutation(len(confidence.breadth_log))
        ]
        prediction, _ = _candidate_prediction(
            candidate,
            fold,
            confidence,
            weight_lower=weight_lower,
            weight_upper=weight_upper,
            breadth_override=shuffled_breadth,
        )
        scores.append(pds(prediction, fold.truth))
    return scores


def evaluate_selected_controls(
    selected_adaptive: Candidate,
    fold_cache: dict[tuple[str, str, str], FoldComponents],
    confidence: dict[str, ConfidenceData],
    args: argparse.Namespace,
) -> list[dict[str, object]]:
    """Evaluate predeclared controls on development folds after selection."""
    rows: list[dict[str, object]] = []
    matched = Candidate(
        family="fixed",
        mode=selected_adaptive.mode,
        cd4_source=selected_adaptive.cd4_source,
        geometry=selected_adaptive.geometry,
        fraction=selected_adaptive.fraction,
    )
    for held in DEVELOPMENT_HELD:
        fold = fold_cache[(held, selected_adaptive.mode, selected_adaptive.geometry)]
        adaptive_prediction, adaptive_weight = _candidate_prediction(
            selected_adaptive,
            fold,
            confidence[selected_adaptive.cd4_source],
            weight_lower=args.weight_lower,
            weight_upper=args.weight_upper,
        )
        constant_prediction, _ = _candidate_prediction(
            matched,
            fold,
            confidence[selected_adaptive.cd4_source],
            weight_lower=args.weight_lower,
            weight_upper=args.weight_upper,
        )
        baseline_score = pds(fold.baseline, fold.truth)
        adaptive_score = pds(adaptive_prediction, fold.truth)
        constant_score = pds(constant_prediction, fold.truth)
        shuffle_scores = _shuffle_scores(
            selected_adaptive,
            fold,
            confidence[selected_adaptive.cd4_source],
            replicates=args.shuffle_replicates,
            seed=args.seed,
            weight_lower=args.weight_lower,
            weight_upper=args.weight_upper,
        )
        rows.append(
            {
                "held": held,
                "adaptive_pds": adaptive_score,
                "matched_constant_pds": constant_score,
                "baseline_pds": baseline_score,
                "adaptive_vs_constant": adaptive_score - constant_score,
                "adaptive_vs_baseline": adaptive_score - baseline_score,
                "shuffle_mean_pds": float(np.mean(shuffle_scores)),
                "shuffle_p95_pds": float(np.quantile(shuffle_scores, 0.95)),
                "adaptive_vs_shuffle_mean": adaptive_score
                - float(np.mean(shuffle_scores)),
                "shuffle_scores": shuffle_scores,
                "weight_min": float(adaptive_weight.min()),
                "weight_mean": float(adaptive_weight.mean()),
                "weight_max": float(adaptive_weight.max()),
            }
        )
    return rows


def evaluate_k562_once(
    selected_fixed: Candidate,
    selected_adaptive: Candidate,
    selected_overall: Candidate,
    args: argparse.Namespace,
    master_features: Sequence[str],
    banks: dict[str, Bank],
    prepared: dict[str, PreparedBank],
    confidence: dict[str, ConfidenceData],
    panel_digest: str,
) -> dict[str, object]:
    """Run one frozen K562 holdout stage; no K562 candidate sweep exists."""
    needed = {selected_fixed, selected_adaptive}
    fold_cache: dict[tuple[str, str], FoldComponents] = {}
    for candidate in needed:
        cache_key = (candidate.mode, candidate.geometry)
        if cache_key in fold_cache:
            continue
        geometries = _fold_components(
            HOLDOUT,
            candidate.mode,
            master_features,
            banks,
            prepared,
            args.whitener_dir,
            panel_digest,
            whitening_gamma=args.whitening_gamma,
            max_rss_gib=args.max_rss_gib,
            whitening_block_rows=args.whitening_block_rows,
        )
        fold_cache[cache_key] = geometries[candidate.geometry]

    rows: list[dict[str, object]] = []
    for role, candidate in (
        ("selected_fixed", selected_fixed),
        ("selected_adaptive", selected_adaptive),
    ):
        fold = fold_cache[(candidate.mode, candidate.geometry)]
        prediction, weights = _candidate_prediction(
            candidate,
            fold,
            confidence[candidate.cd4_source],
            weight_lower=args.weight_lower,
            weight_upper=args.weight_upper,
        )
        baseline_score = pds(fold.baseline, fold.truth)
        score = pds(prediction, fold.truth)
        rows.append(
            {
                "role": role,
                "candidate": candidate.key,
                "pds": score,
                "baseline_pds": baseline_score,
                "delta_vs_baseline": score - baseline_score,
                "weight_min": float(weights.min()),
                "weight_mean": float(weights.mean()),
                "weight_max": float(weights.max()),
            }
        )

    adaptive_fold = fold_cache[(selected_adaptive.mode, selected_adaptive.geometry)]
    matched = Candidate(
        family="fixed",
        mode=selected_adaptive.mode,
        cd4_source=selected_adaptive.cd4_source,
        geometry=selected_adaptive.geometry,
        fraction=selected_adaptive.fraction,
    )
    matched_prediction, _ = _candidate_prediction(
        matched,
        adaptive_fold,
        confidence[selected_adaptive.cd4_source],
        weight_lower=args.weight_lower,
        weight_upper=args.weight_upper,
    )
    matched_score = pds(matched_prediction, adaptive_fold.truth)
    adaptive_score = next(
        float(row["pds"]) for row in rows if row["role"] == "selected_adaptive"
    )
    fixed_score = next(
        float(row["pds"]) for row in rows if row["role"] == "selected_fixed"
    )
    shuffle_scores = _shuffle_scores(
        selected_adaptive,
        adaptive_fold,
        confidence[selected_adaptive.cd4_source],
        replicates=args.shuffle_replicates,
        seed=args.seed,
        weight_lower=args.weight_lower,
        weight_upper=args.weight_upper,
    )
    empirical_p = (1 + sum(score >= adaptive_score for score in shuffle_scores)) / (
        len(shuffle_scores) + 1
    )
    return {
        "held": HOLDOUT,
        "evaluation_count": 1,
        "selection_was_frozen_before_truth_scoring": True,
        "selected_overall_role": (
            "selected_adaptive"
            if selected_overall == selected_adaptive
            else "selected_fixed"
        ),
        "rows": rows,
        "adaptive_matched_constant": {
            "candidate": matched.key,
            "pds": matched_score,
            "adaptive_vs_matched_constant": adaptive_score - matched_score,
            "adaptive_vs_selected_fixed": adaptive_score - fixed_score,
        },
        "adaptive_shuffled_breadth_control": {
            "replicates": len(shuffle_scores),
            "scores": shuffle_scores,
            "mean_pds": float(np.mean(shuffle_scores)),
            "p95_pds": float(np.quantile(shuffle_scores, 0.95)),
            "adaptive_vs_shuffle_mean": adaptive_score - float(np.mean(shuffle_scores)),
            "one_sided_empirical_p": float(empirical_p),
        },
    }


def _write_csv(path: Path, rows: Sequence[dict[str, object]]) -> None:
    if not rows:
        return
    scalar_rows = []
    for row in rows:
        scalar_rows.append(
            {
                key: json.dumps(value, sort_keys=True)
                if isinstance(value, (dict, list))
                else value
                for key, value in row.items()
            }
        )
    fieldnames = sorted({key for row in scalar_rows for key in row})
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(scalar_rows)


def run(args: argparse.Namespace) -> dict[str, object]:
    panel = (
        pd.read_csv(args.panel_csv)["target_gene"]
        .dropna()
        .astype(str)
        .str.strip()
        .tolist()
    )
    official = (
        pd.read_csv(args.genes_csv)["gene_name"]
        .dropna()
        .astype(str)
        .str.strip()
        .tolist()
    )
    if len(panel) != len(set(panel)):
        raise ValueError("panel has duplicate target symbols")
    if len(official) != len(set(official)):
        raise ValueError("official gene axis has duplicate symbols")
    panel_set = set(panel)
    panel_digest = sha256_lines(panel)

    main = {
        "k562_gwps": Bank.load("k562_gwps", Path(args.k562)),
        "orion_hct116": Bank.load("orion_hct116", Path(args.hct116)),
        "orion_hek293t": Bank.load("orion_hek293t", Path(args.hek293t)),
    }
    cd4_root = Path(args.cd4_root)
    cd4 = {name: Bank.load(name, cd4_root / name, cd4=True) for name in CD4_NAMES}
    banks = {**main, **cd4}
    perturbations = [
        target
        for target in panel
        if all(target in bank.pert_pos for bank in banks.values())
    ]
    if len(perturbations) < args.min_shared_targets:
        raise ValueError(
            f"only {len(perturbations)} panel targets are shared by all banks; "
            f"require {args.min_shared_targets}"
        )

    # Exclude all 300 panel target coordinates before source centering,
    # gene normalization, blending, or optional dense whitening.
    master_features = [
        gene
        for gene in official
        if gene not in panel_set
        and any(gene in bank.gene_pos for bank in banks.values())
    ]
    if panel_set & set(master_features):
        raise AssertionError("panel target coordinate leaked into master feature axis")
    if len(master_features) < args.min_master_genes:
        raise ValueError(
            f"only {len(master_features)} non-panel union genes; "
            f"require {args.min_master_genes}"
        )

    confidence = load_confidence_data(cd4, perturbations)
    output = Path(args.output).expanduser().resolve()
    output.parent.mkdir(parents=True, exist_ok=True)
    scratch_parent = Path(args.scratch_parent).expanduser().resolve()
    scratch_parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(
        prefix="exp051-", dir=scratch_parent
    ) as scratch_name:
        scratch = Path(scratch_name)
        prepared = {}
        for name, bank in banks.items():
            prepared[name] = prepare_bank(
                bank,
                perturbations,
                master_features,
                scratch,
                chunk_rows=args.chunk_rows,
                max_rss_gib=args.max_rss_gib,
            )
            print(
                f"prepared {name:16} panel={len(perturbations):3d} "
                f"union_genes={len(master_features):5d} "
                f"peak={peak_rss_gib():.2f} GiB",
                flush=True,
            )

        development_rows, fold_cache = evaluate_development(
            args,
            master_features,
            banks,
            prepared,
            confidence,
            panel_digest,
        )
        selected_fixed, fixed_summary = summarize_and_select(
            development_rows,
            "fixed",
            tie_tolerance=args.tie_tolerance,
        )
        selected_adaptive, adaptive_summary = summarize_and_select(
            development_rows,
            "adaptive",
            tie_tolerance=args.tie_tolerance,
        )
        selected_overall = select_overall(
            selected_fixed,
            fixed_summary,
            selected_adaptive,
            adaptive_summary,
            tie_tolerance=args.tie_tolerance,
        )

        development_controls = evaluate_selected_controls(
            selected_adaptive,
            fold_cache,
            confidence,
            args,
        )

        # This is deliberately after every hyperparameter has been frozen.
        k562_holdout = evaluate_k562_once(
            selected_fixed,
            selected_adaptive,
            selected_overall,
            args,
            master_features,
            banks,
            prepared,
            confidence,
            panel_digest,
        )

        result: dict[str, object] = {
            "status": "complete",
            "protocol": {
                "development_folds": list(DEVELOPMENT_HELD),
                "lineage_analogue_holdout": HOLDOUT,
                "k562_holdout_evaluations": 1,
                "selection_uses_k562_pds": False,
                "adaptive_inputs": "compact CD4 obs.csv breadth and QC fields only",
                "forbidden_inputs_not_read": [
                    "K562_comparison.suppl_table.csv",
                    "challenge outcomes",
                    "leaderboard outcomes",
                ],
                "perturbations": len(perturbations),
                "perturbation_policy": "VCC-panel intersection across all six banks",
                "master_union_genes": len(master_features),
                "gene_policy": (
                    "official-axis union across banks minus all 300 VCC target genes "
                    "before every transform"
                ),
                "truth_space": "held donor centered raw effects",
                "prediction_modes": list(args.modes),
                "cd4_sources": list(args.cd4_sources),
                "fractions": list(args.fractions),
                "slopes": list(args.slopes),
                "qc_strengths": list(args.qc_strengths),
                "weight_bounds": [args.weight_lower, args.weight_upper],
                "tie_tolerance": args.tie_tolerance,
                "shuffle_replicates": args.shuffle_replicates,
                "whitening": (
                    {
                        "enabled": True,
                        "kind": "audited fold-specific production union",
                        "gamma": args.whitening_gamma,
                        "root": str(args.whitener_dir),
                    }
                    if args.whitener_dir is not None
                    else {"enabled": False}
                ),
                "no_cell_eval2": True,
                "no_submission": True,
            },
            "digests": {
                "panel_sha256": panel_digest,
                "master_feature_axis_sha256": sha256_lines(master_features),
            },
            "confidence_provenance": {
                name: value.provenance for name, value in confidence.items()
            },
            "selected_fixed": asdict(selected_fixed),
            "selected_adaptive": asdict(selected_adaptive),
            "selected_overall": asdict(selected_overall),
            "fixed_development_summary": fixed_summary,
            "adaptive_development_summary": adaptive_summary,
            "development_rows": development_rows,
            "selected_adaptive_development_controls": development_controls,
            "k562_holdout": k562_holdout,
            "peak_rss_gib": peak_rss_gib(),
        }

    output.write_text(json.dumps(result, indent=2) + "\n")
    _write_csv(output.with_name("development_rows.csv"), development_rows)
    _write_csv(
        output.with_name("k562_holdout_rows.csv"),
        list(k562_holdout["rows"]),
    )
    print(f"wrote {output}", flush=True)
    print(
        json.dumps(
            {
                "selected_fixed": asdict(selected_fixed),
                "selected_adaptive": asdict(selected_adaptive),
                "selected_overall": asdict(selected_overall),
                "k562_holdout": k562_holdout,
                "peak_rss_gib": peak_rss_gib(),
            },
            indent=2,
        ),
        flush=True,
    )
    return result


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cd4-root", required=True, type=Path)
    parser.add_argument("--k562", default=MAIN_PATHS["k562_gwps"], type=Path)
    parser.add_argument("--hct116", default=MAIN_PATHS["orion_hct116"], type=Path)
    parser.add_argument("--hek293t", default=MAIN_PATHS["orion_hek293t"], type=Path)
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
        "--output",
        default=HERE / "results.json",
        type=Path,
    )
    parser.add_argument(
        "--scratch-parent",
        default=HERE,
        type=Path,
        help="temporary memory-map parent; files are removed after the run",
    )
    parser.add_argument(
        "--fractions",
        nargs="+",
        type=float,
        default=[0.1, 0.2, 1 / 3, 0.5, 0.67],
    )
    parser.add_argument(
        "--slopes",
        nargs="+",
        type=float,
        default=[0.5, 1.0, 2.0],
    )
    parser.add_argument(
        "--qc-strengths",
        nargs="+",
        type=float,
        default=[0.0, 0.5, 1.0],
    )
    parser.add_argument(
        "--modes",
        nargs="+",
        choices=MODES,
        default=list(MODES),
    )
    parser.add_argument(
        "--cd4-sources",
        nargs="+",
        choices=CD4_SOURCE_NAMES,
        default=list(CD4_SOURCE_NAMES),
    )
    parser.add_argument("--weight-lower", type=float, default=0.02)
    parser.add_argument("--weight-upper", type=float, default=0.80)
    parser.add_argument("--tie-tolerance", type=float, default=0.001)
    parser.add_argument("--shuffle-replicates", type=int, default=31)
    parser.add_argument("--seed", type=int, default=20260822)
    parser.add_argument("--chunk-rows", type=int, default=128)
    parser.add_argument("--whitening-block-rows", type=int, default=32)
    parser.add_argument(
        "--whitener-dir",
        type=Path,
        default=None,
        help="optional root of audited fold-specific production-union operators",
    )
    parser.add_argument("--whitening-gamma", type=float, default=0.65)
    parser.add_argument(
        "--max-rss-gib",
        type=float,
        default=3.75,
        help="internal soft stop; outer process-tree guard must be 4 GiB",
    )
    parser.add_argument("--min-shared-targets", type=int, default=240)
    parser.add_argument("--min-master-genes", type=int, default=9_000)
    parser.add_argument(
        "--guarded",
        action="store_true",
        help="confirm this CLI is running inside scripts/run_with_memory_guard.py",
    )
    return parser


def validate_args(args: argparse.Namespace) -> None:
    if not args.guarded:
        raise SystemExit(
            "Refusing an unguarded real-data run. Invoke through "
            "scripts/run_with_memory_guard.py with --max-rss-gib 4 and pass "
            "--guarded to this script."
        )
    if args.max_rss_gib <= 0 or args.max_rss_gib >= 4:
        raise SystemExit("--max-rss-gib must be > 0 and < 4")
    if args.chunk_rows < 1 or args.whitening_block_rows < 1:
        raise SystemExit("chunk sizes must be positive")
    if args.shuffle_replicates < 1:
        raise SystemExit("--shuffle-replicates must be positive")
    if not 0 <= args.weight_lower < args.weight_upper <= 1:
        raise SystemExit("weight bounds must satisfy 0 <= lower < upper <= 1")
    if any(
        not args.weight_lower <= fraction <= args.weight_upper
        for fraction in args.fractions
    ):
        raise SystemExit("every fraction must lie within the weight bounds")
    if any(slope <= 0 for slope in args.slopes):
        raise SystemExit("all slopes must be positive")
    if any(value < 0 for value in args.qc_strengths):
        raise SystemExit("QC strengths must be non-negative")
    if args.whitener_dir is not None:
        args.whitener_dir = args.whitener_dir.expanduser().resolve()


def main() -> int:
    args = build_parser().parse_args()
    validate_args(args)
    run(args)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
