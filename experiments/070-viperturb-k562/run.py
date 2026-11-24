#!/usr/bin/env python3
"""Evaluate a frozen VIPerturb addition to experiment-064 arm D.

The full held HCT116/HEK293T gallery is retained.  Only targets available in
incumbent K562, PerturbMatch, and the official filtered VIPerturb donor are
modified.  The candidate K562 slot is:

    0.5 * (0.5 * incumbent + 0.5 * PerturbMatch) + 0.5 * VIPerturb

on genes shared by all three K562 views.  Every other target/gene remains the
frozen full-gallery arm-D baseline.
"""

from __future__ import annotations

import argparse
import csv
import gc
import gzip
import hashlib
import importlib.util
import json
import os
import resource
import sys
from collections.abc import Iterable, Sequence
from dataclasses import dataclass
from pathlib import Path

os.environ.setdefault("OPENBLAS_NUM_THREADS", "8")
os.environ.setdefault("OMP_NUM_THREADS", "8")
os.environ.setdefault("MKL_NUM_THREADS", "8")

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
EXP064_PATH = ROOT / "experiments/064-pm-cd4-composition/run.py"
SPEC = importlib.util.spec_from_file_location("exp070_exp064", EXP064_PATH)
if SPEC is None or SPEC.loader is None:
    raise ImportError(f"cannot import frozen experiment-064 code: {EXP064_PATH}")
EXP064 = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = EXP064
SPEC.loader.exec_module(EXP064)
BASE = EXP064.BASE

FOLDS = ("orion_hct116", "orion_hek293t")
CD4_NAMES = ("cd4_rest", "cd4_stim8hr", "cd4_stim48hr")
PM_FRACTION = 0.50
VIP_FRACTION = 0.50
CD4_MODE = "gene_norm"
CD4_FRACTION = 0.05
WHITEN_GAMMA = 0.65
BOOTSTRAP_REPS = 20_000
BOOTSTRAP_SEED = 7_002_026
NULL_REPS = 1_000
MIN_NULL_REPS = 100
TARGET_NULL_SEED = 70_100_001
GENE_NULL_SEED = 70_200_001
GENE_NULL_CHUNK_REPS = 8
MACRO_GATE = 0.003
MIN_MODIFIED_TARGETS = 70
MIN_NON_LOW_TARGETS = 20
LOW_CELL_THRESHOLD = 30
EXPECTED_VIP_COVERAGE = 124
VIP_SOURCE_MD5 = "ac52c41233fc7a5a5c3810e0a075007f"
VIP_MANIFEST_MD5 = "c095c8f2d16136753c22ef5d9b30a22b"
EXPECTED_WHITENER_SHA256 = dict(EXP064.EXPECTED_WHITENER_SHA256)
JSONLITE_FOUR_DECIMAL_ATOL = 5.0001e-5


def peak_rss_gib() -> float:
    value = float(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
    return value / 1024**3 if sys.platform == "darwin" else value / 1024**2


def enforce_limit(limit_gib: float, stage: str) -> None:
    if limit_gib > 0 and peak_rss_gib() >= limit_gib:
        raise MemoryError(
            f"{stage}: peak RSS {peak_rss_gib():.3f} GiB >= {limit_gib:.3f} GiB"
        )


def sha256_file(path: Path, block_bytes: int = 8 << 20) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while block := handle.read(block_bytes):
            digest.update(block)
    return digest.hexdigest()


def md5_file(path: Path, block_bytes: int = 8 << 20) -> str:
    digest = hashlib.md5()
    with path.open("rb") as handle:
        while block := handle.read(block_bytes):
            digest.update(block)
    return digest.hexdigest()


def sha256_lines(values: Iterable[str]) -> str:
    return hashlib.sha256("\n".join(values).encode()).hexdigest()


def read_lines(path: Path) -> tuple[str, ...]:
    values = tuple(
        line.strip() for line in path.read_text().splitlines() if line.strip()
    )
    if not values or len(values) != len(set(values)):
        raise ValueError(f"{path}: empty or duplicate labels")
    return values


def parse_bool(value: object) -> bool:
    normalized = str(value).strip().lower()
    if normalized in {"true", "t", "1"}:
        return True
    if normalized in {"false", "f", "0"}:
        return False
    raise ValueError(f"invalid boolean value {value!r}")


def read_effect_table(
    path: Path,
    expected_genes: Sequence[str],
    expected_perts: Sequence[str],
) -> np.ndarray:
    with gzip.open(path, "rt", newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        header = next(reader)
        if not header or header[0] != "target_gene":
            raise ValueError(f"{path}: first column must be target_gene")
        genes = tuple(header[1:])
        if genes != tuple(expected_genes):
            raise ValueError(f"{path}: effect-table gene axis mismatch")
        output = np.empty((len(expected_perts), len(genes)), dtype=np.float64)
        seen: list[str] = []
        for index, row in enumerate(reader):
            if index >= len(expected_perts):
                raise ValueError(f"{path}: too many perturbation rows")
            if len(row) != len(header):
                raise ValueError(f"{path}: row {index + 2} has wrong width")
            seen.append(row[0])
            output[index] = np.asarray(row[1:], dtype=np.float64)
    if tuple(seen) != tuple(expected_perts):
        raise ValueError(f"{path}: perturbation row order mismatch")
    if not np.isfinite(output).all():
        raise ValueError(f"{path}: non-finite effect values")
    return output


def read_gene_value_table(
    path: Path,
    expected_genes: Sequence[str],
    value_column: str,
) -> np.ndarray:
    with gzip.open(path, "rt", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames != ["gene", value_column]:
            raise ValueError(f"{path}: expected gene/{value_column} columns")
        genes: list[str] = []
        values: list[float] = []
        for row in reader:
            genes.append(row["gene"])
            values.append(float(row[value_column]))
    if tuple(genes) != tuple(expected_genes):
        raise ValueError(f"{path}: gene axis mismatch")
    output = np.asarray(values, dtype=np.float64)
    if not np.isfinite(output).all():
        raise ValueError(f"{path}: non-finite {value_column} values")
    return output


def panel_median_centered_norm(
    effects: np.ndarray,
    panel_center: np.ndarray,
) -> float:
    return float(np.median(np.linalg.norm(effects - panel_center, axis=1)))


def legacy_r_recycled_panel_norm(
    effects: np.ndarray,
    panel_center: np.ndarray,
) -> float:
    """Reproduce R's accidental column-major matrix/vector recycling."""
    flat = np.asarray(effects, dtype=np.float64).ravel(order="F")
    recycled = np.resize(np.asarray(panel_center, dtype=np.float64), flat.size)
    centered = (flat - recycled).reshape(effects.shape, order="F")
    return float(np.median(np.linalg.norm(centered, axis=1)))


def resolve_panel_median_norm(
    effects: np.ndarray,
    panel_center: np.ndarray,
    stored: float,
) -> tuple[float, str, float | None]:
    canonical = panel_median_centered_norm(effects, panel_center)
    if np.isclose(stored, canonical, rtol=1e-12, atol=1e-12):
        return canonical, "canonical_full_precision", None
    if np.isclose(
        stored,
        canonical,
        rtol=0.0,
        atol=JSONLITE_FOUR_DECIMAL_ATOL,
    ):
        return canonical, "canonical_jsonlite_four_decimal", None
    legacy = legacy_r_recycled_panel_norm(effects, panel_center)
    if np.isclose(
        stored,
        legacy,
        rtol=0.0,
        atol=JSONLITE_FOUR_DECIMAL_ATOL,
    ):
        return canonical, "legacy_r_recycling_jsonlite_four_decimal", legacy
    raise ValueError(
        "VIPerturb panel-centered scale artifact is inconsistent: "
        f"stored={stored:.17g} canonical={canonical:.17g} "
        f"legacy_r_recycling={legacy:.17g}"
    )


@dataclass(frozen=True)
class VIPAudit:
    metadata: dict[str, object]
    target_rows: dict[str, dict[str, object]]
    compact_hashes: dict[str, str]
    source_variance: np.ndarray


@dataclass(frozen=True)
class FoldCache:
    held: str
    baseline_logits: np.ndarray
    vip_projection: np.ndarray
    row_factors: np.ndarray
    baseline_per_target: np.ndarray


def _align_scaled_effects(
    raw_effects: np.ndarray,
    center: np.ndarray,
    median_norm: float,
    perts: Sequence[str],
    targets: Sequence[str],
    official_count: int,
    mapped_source: np.ndarray,
    mapped_official: np.ndarray,
    panel_positions: np.ndarray,
    *,
    name: str,
) -> BASE.PreparedSource:
    if not np.isfinite(median_norm) or median_norm <= 0:
        raise ValueError(f"{name}: median centered norm must be positive")
    scaled = (raw_effects - center) * (BASE.EFFECT_SCALE / median_norm)
    target_position = {target: index for index, target in enumerate(targets)}
    output = np.zeros((len(targets), official_count), dtype=np.float64)
    available = np.zeros(len(targets), dtype=bool)
    for source_row, target in enumerate(perts):
        destination_row = target_position.get(target)
        if destination_row is None:
            continue
        output[destination_row, mapped_official] = scaled[source_row, mapped_source]
        available[destination_row] = True
    output[:, panel_positions] = 0.0
    return BASE.PreparedSource(
        name=name,
        effects=output,
        target_available=available,
        mean_on_axis=np.zeros(official_count, dtype=np.float64),
        raw_median_norm=median_norm,
        mapped_gene_count=len(mapped_official),
        duplicate_symbols=(),
    )


def load_viperturb_compact(
    path: Path,
    targets: Sequence[str],
    official: Sequence[str],
    panel_positions: np.ndarray,
    *,
    expected_coverage: int,
    expected_source_md5: str,
    expected_manifest_md5: str,
) -> tuple[BASE.PreparedSource, BASE.PreparedSource, VIPAudit]:
    path = path.expanduser().resolve()
    metadata = json.loads((path / "meta.json").read_text())
    checks = {
        "experiment": "070-viperturb-k562",
        "zenodo_record": "18460279",
        "source_file": "genome_wide_filtered.rds",
        "source_md5": expected_source_md5,
        "manifest_file": "genome_wide_manifest.txt",
        "manifest_md5": expected_manifest_md5,
        "manifest_filtered_vcc_coverage": expected_coverage,
        "normalization_scope": "all_filtered_source_perturbations",
        "low_cell_threshold": LOW_CELL_THRESHOLD,
    }
    for key, expected in checks.items():
        if metadata.get(key) != expected:
            raise ValueError(
                f"VIPerturb metadata {key}={metadata.get(key)!r} != {expected!r}"
            )

    genes = read_lines(path / "genes.txt")
    perts = read_lines(path / "perts.txt")
    source_perts = read_lines(path / "source_perts.txt")
    if len(perts) != expected_coverage:
        raise ValueError(
            f"VIPerturb compact coverage {len(perts)} != {expected_coverage}"
        )
    if tuple(metadata.get("filtered_targets", ())) != perts:
        raise ValueError("VIPerturb metadata target order disagrees with perts.txt")
    if int(metadata.get("n_perturbations", -1)) != len(perts):
        raise ValueError("VIPerturb metadata perturbation count mismatch")
    if int(metadata.get("n_genes", -1)) != len(genes):
        raise ValueError("VIPerturb metadata gene count mismatch")
    if int(metadata.get("n_source_perturbations", -1)) != len(source_perts):
        raise ValueError("VIPerturb source-wide perturbation count mismatch")
    if len(source_perts) <= len(perts):
        raise ValueError("VIPerturb normalization catalogue is not source-wide")
    if not set(perts).issubset(source_perts):
        raise ValueError("VIPerturb panel rows are absent from source catalogue")
    source_median_norm = float(metadata.get("source_median_centered_norm", 0.0))
    panel_median_norm = float(metadata.get("panel_median_centered_norm", 0.0))
    if source_median_norm <= 0 or panel_median_norm <= 0:
        raise ValueError("VIPerturb normalization scales must be positive")

    expected_artifacts = {
        str(name): str(value)
        for name, value in dict(metadata.get("artifacts_md5", {})).items()
    }
    required_artifacts = {
        "genes.txt",
        "perts.txt",
        "source_perts.txt",
        "effects.tsv.gz",
        "control.tsv.gz",
        "source_center.tsv.gz",
        "source_variance.tsv.gz",
        "panel_center.tsv.gz",
        "target_audit.csv",
        "metadata_column_audit.csv",
        "object_audit.txt",
        "r_session_info.txt",
    }
    if set(expected_artifacts) != required_artifacts:
        raise ValueError("VIPerturb compact artifact hash set is incomplete")
    compact_hashes: dict[str, str] = {}
    for name in sorted(required_artifacts):
        actual = md5_file(path / name)
        if actual != expected_artifacts[name]:
            raise ValueError(f"VIPerturb compact hash mismatch for {name}")
        compact_hashes[name] = actual
    if (
        metadata.get("source_perturbation_order_md5")
        != compact_hashes["source_perts.txt"]
    ):
        raise ValueError("VIPerturb source perturbation order hash mismatch")

    effects = read_effect_table(path / "effects.tsv.gz", genes, perts)
    source_center = read_gene_value_table(
        path / "source_center.tsv.gz", genes, "source_center"
    )
    source_variance = read_gene_value_table(
        path / "source_variance.tsv.gz", genes, "source_variance"
    )
    panel_center = read_gene_value_table(
        path / "panel_center.tsv.gz", genes, "panel_center"
    )
    if np.any(source_variance < 0):
        raise ValueError("VIPerturb source variance contains negative values")
    if not np.allclose(
        panel_center,
        effects.mean(axis=0),
        rtol=1e-12,
        atol=1e-12,
    ):
        raise ValueError(
            "VIPerturb panel-center artifact disagrees with raw panel rows"
        )
    (
        computed_panel_norm,
        panel_norm_validation,
        legacy_panel_norm,
    ) = resolve_panel_median_norm(
        effects,
        panel_center,
        panel_median_norm,
    )
    metadata = dict(metadata)
    metadata["panel_median_centered_norm_stored"] = panel_median_norm
    metadata["panel_median_centered_norm"] = computed_panel_norm
    metadata["panel_median_centered_norm_validation"] = panel_norm_validation
    if legacy_panel_norm is not None:
        metadata["panel_median_centered_norm_legacy_r_recycling"] = legacy_panel_norm
    panel_median_norm = computed_panel_norm
    gene_pos = {gene: index for index, gene in enumerate(official)}
    mapped_source: list[int] = []
    mapped_official: list[int] = []
    for source, gene in enumerate(genes):
        destination = gene_pos.get(gene)
        if destination is not None:
            mapped_source.append(source)
            mapped_official.append(destination)
    if len(mapped_source) != len(genes):
        raise ValueError("VIPerturb compact axis contains a non-official gene")

    mapped_source_array = np.asarray(mapped_source, dtype=np.int64)
    mapped_official_array = np.asarray(mapped_official, dtype=np.int64)
    source_wide = _align_scaled_effects(
        effects,
        source_center,
        source_median_norm,
        perts,
        targets,
        len(official),
        mapped_source_array,
        mapped_official_array,
        panel_positions,
        name="viperturb_k562_source_wide",
    )
    panel_centered = _align_scaled_effects(
        effects,
        panel_center,
        panel_median_norm,
        perts,
        targets,
        len(official),
        mapped_source_array,
        mapped_official_array,
        panel_positions,
        name="viperturb_k562_panel_centered_ablation",
    )

    target_rows: dict[str, dict[str, object]] = {}
    with (path / "target_audit.csv").open(newline="") as handle:
        reader = csv.DictReader(handle)
        required = {
            "target_gene",
            "manifest_cell_count",
            "observed_cell_count",
            "count_difference",
            "low_cell_manifest",
            "low_cell_observed",
        }
        if reader.fieldnames is None or set(reader.fieldnames) != required:
            raise ValueError("VIPerturb target_audit.csv schema mismatch")
        for row in reader:
            target = row["target_gene"].strip()
            target_rows[target] = {
                "manifest_cell_count": int(row["manifest_cell_count"]),
                "observed_cell_count": int(row["observed_cell_count"]),
                "count_difference": int(row["count_difference"]),
                "low_cell_manifest": parse_bool(row["low_cell_manifest"]),
                "low_cell_observed": parse_bool(row["low_cell_observed"]),
            }
    if tuple(target_rows) != perts:
        raise ValueError("VIPerturb target audit order disagrees with perts.txt")

    source_variance_official = np.zeros(len(official), dtype=np.float64)
    source_variance_official[mapped_official_array] = source_variance[
        mapped_source_array
    ]
    return (
        source_wide,
        panel_centered,
        VIPAudit(
            metadata,
            target_rows,
            compact_hashes,
            source_variance_official,
        ),
    )


def build_pm_slot_with_fallback(
    incumbent: BASE.PreparedSource,
    perturbmatch: BASE.PreparedSource,
    shared_columns: np.ndarray,
) -> tuple[BASE.PreparedSource, dict[str, object]]:
    if incumbent.effects.shape != perturbmatch.effects.shape:
        raise ValueError("incumbent and PerturbMatch matrices differ in shape")
    columns = np.asarray(shared_columns, dtype=np.int64)
    rows = incumbent.target_available & perturbmatch.target_available
    effects = np.asarray(incumbent.effects, dtype=np.float64).copy()
    if rows.any() and len(columns):
        effects[np.ix_(rows, columns)] = (1.0 - PM_FRACTION) * incumbent.effects[
            np.ix_(rows, columns)
        ] + PM_FRACTION * perturbmatch.effects[np.ix_(rows, columns)]
    return (
        BASE.PreparedSource(
            name="k562_gwps_pm50_fallback",
            effects=effects,
            target_available=np.asarray(incumbent.target_available, dtype=bool).copy(),
            mean_on_axis=np.asarray(incumbent.mean_on_axis, dtype=np.float64).copy(),
            raw_median_norm=incumbent.raw_median_norm,
            mapped_gene_count=incumbent.mapped_gene_count,
            duplicate_symbols=incumbent.duplicate_symbols,
        ),
        {
            "targets_blended": int(rows.sum()),
            "targets_incumbent_fallback": int((~rows).sum()),
            "shared_non_panel_genes": int(len(columns)),
        },
    )


def build_vip_slot(
    incumbent: BASE.PreparedSource,
    perturbmatch: BASE.PreparedSource,
    pm_slot: BASE.PreparedSource,
    viperturb: BASE.PreparedSource,
    shared_columns: np.ndarray,
    targets: Sequence[str],
) -> tuple[BASE.PreparedSource, np.ndarray, dict[str, object]]:
    shapes = {
        incumbent.effects.shape,
        perturbmatch.effects.shape,
        pm_slot.effects.shape,
        viperturb.effects.shape,
    }
    if len(shapes) != 1:
        raise ValueError("K562 view matrices differ in shape")
    modified = (
        incumbent.target_available
        & perturbmatch.target_available
        & viperturb.target_available
    )
    columns = np.asarray(shared_columns, dtype=np.int64)
    effects = np.asarray(pm_slot.effects, dtype=np.float64).copy()
    if modified.any() and len(columns):
        effects[np.ix_(modified, columns)] = (1.0 - VIP_FRACTION) * pm_slot.effects[
            np.ix_(modified, columns)
        ] + VIP_FRACTION * viperturb.effects[np.ix_(modified, columns)]
    modified_targets = [
        target for target, keep in zip(targets, modified, strict=True) if keep
    ]
    return (
        BASE.PreparedSource(
            name="k562_pm50_vip50_fallback",
            effects=effects,
            target_available=np.asarray(pm_slot.target_available, dtype=bool).copy(),
            mean_on_axis=np.asarray(pm_slot.mean_on_axis, dtype=np.float64).copy(),
            raw_median_norm=pm_slot.raw_median_norm,
            mapped_gene_count=pm_slot.mapped_gene_count,
            duplicate_symbols=pm_slot.duplicate_symbols,
        ),
        modified,
        {
            "formula": "0.5*(0.5*incumbent+0.5*PerturbMatch)+0.5*VIPerturb",
            "modified_target_count": int(modified.sum()),
            "modified_targets": modified_targets,
            "shared_non_panel_genes": int(len(columns)),
            "baseline_preserved_target_count": int((~modified).sum()),
        },
    )


def freeze_complete_cd4_mean3(cd4: BASE.CD4Prepared) -> BASE.CD4Prepared:
    complete = np.asarray(cd4.target_condition_count) == len(CD4_NAMES)
    effects = np.asarray(cd4.mean_effect, dtype=np.float64).copy()
    votes = np.asarray(cd4.gene_votes, dtype=np.uint8).copy()
    effects[~complete] = 0.0
    votes[~complete] = 0
    return BASE.CD4Prepared(
        mode=cd4.mode,
        mean_effect=effects,
        gene_votes=votes,
        target_condition_count=np.asarray(
            cd4.target_condition_count, dtype=np.uint8
        ).copy(),
        condition_details=list(cd4.condition_details),
    )


def apply_operator_with_frozen_scale(
    matrix: np.ndarray,
    operator: BASE.Operator,
    wcols: np.ndarray,
    scale: float,
) -> np.ndarray:
    output = np.asarray(matrix, dtype=np.float64).copy()
    weights = operator.evals ** (-WHITEN_GAMMA / 2.0)
    transformed = ((output[:, wcols] @ operator.evecs) * weights) @ operator.evecs.T
    output[:, wcols] = transformed * scale
    return output


def add_cd4_with_frozen_scale(
    baseline_post: np.ndarray,
    cd4: BASE.CD4Prepared,
    scale: float,
) -> np.ndarray:
    output = np.asarray(baseline_post, dtype=np.float64).copy()
    available = cd4.gene_votes > 0
    output[available] += CD4_FRACTION * cd4.mean_effect[available]
    output *= scale
    return output


def paired_target_cluster_bootstrap(
    per_fold_delta: dict[str, np.ndarray],
    *,
    reps: int,
    seed: int,
) -> dict[str, object]:
    arrays = [np.asarray(per_fold_delta[held], dtype=np.float64) for held in FOLDS]
    if len({len(values) for values in arrays}) != 1:
        raise ValueError("held folds do not share the frozen full gallery")
    cluster_delta = np.mean(np.stack(arrays, axis=1), axis=1)
    rng = np.random.default_rng(seed)
    draws = np.empty(reps, dtype=np.float64)
    for lo in range(0, reps, 1_000):
        hi = min(lo + 1_000, reps)
        indices = rng.integers(
            0, len(cluster_delta), size=(hi - lo, len(cluster_delta))
        )
        draws[lo:hi] = cluster_delta[indices].mean(axis=1)
    return {
        "method": (
            "paired target-cluster bootstrap over the common full gallery; "
            "one target resample preserves baseline/candidate pairing and both folds"
        ),
        "reps": int(reps),
        "seed": int(seed),
        "estimate": float(cluster_delta.mean()),
        "ci95_percentile": [
            float(np.percentile(draws, 2.5)),
            float(np.percentile(draws, 97.5)),
        ],
        "positive_probability": float(np.mean(draws > 0)),
        "targets": int(len(cluster_delta)),
    }


def pds_per_target_from_logits(logits: np.ndarray) -> np.ndarray:
    values = np.asarray(logits, dtype=np.float64)
    if values.ndim != 2 or values.shape[0] != values.shape[1]:
        raise ValueError("PDS logits must be a square target-by-target matrix")
    own = np.diag(values)[:, None]
    ranks = (values > own).sum(axis=1)
    ranks = ranks + 0.5 * ((values == own).sum(axis=1) - 1)
    return 1.0 - ranks / (values.shape[0] - 1)


def balanced_rank_strata(values: np.ndarray, bins: int) -> np.ndarray:
    """Assign deterministic, near-equal quantile bins even when values tie."""
    array = np.asarray(values, dtype=np.float64)
    if array.ndim != 1 or len(array) < bins:
        raise ValueError(f"need at least {bins} finite values for {bins} strata")
    if not np.isfinite(array).all():
        raise ValueError("stratification values must be finite")
    order = np.lexsort((np.arange(len(array)), array))
    labels = np.empty(len(array), dtype=np.int64)
    labels[order] = np.minimum(
        (np.arange(len(array), dtype=np.int64) * bins) // len(array),
        bins - 1,
    )
    return labels


def stratified_permutation(labels: np.ndarray, rng: np.random.Generator) -> np.ndarray:
    strata = np.asarray(labels, dtype=np.int64)
    permutation = np.arange(len(strata), dtype=np.int64)
    for label in np.unique(strata):
        members = np.flatnonzero(strata == label)
        permutation[members] = rng.permutation(members)
    return permutation


def candidate_pds_from_vip_rows(
    cache: FoldCache,
    vip_rows: np.ndarray,
    pm_rows: np.ndarray,
) -> np.ndarray:
    if vip_rows.shape != pm_rows.shape:
        raise ValueError("VIPerturb and PM null matrices differ in shape")
    logits = cache.baseline_logits.copy()
    contribution = (vip_rows - pm_rows) @ cache.vip_projection
    modified_rows = np.flatnonzero(cache.row_factors != 0)
    if len(modified_rows) != len(vip_rows):
        raise ValueError("fold cache modified-row mapping is inconsistent")
    logits[modified_rows] += contribution * cache.row_factors[modified_rows, None]
    return pds_per_target_from_logits(logits)


def summarize_null(
    values: np.ndarray,
    observed: float,
    *,
    name: str,
    seed: int,
    method: str,
) -> dict[str, object]:
    array = np.asarray(values, dtype=np.float64)
    if array.ndim != 1 or not len(array) or not np.isfinite(array).all():
        raise ValueError(f"{name}: null values must be a finite non-empty vector")
    return {
        "name": name,
        "method": method,
        "reps": int(len(array)),
        "seed": int(seed),
        "mean": float(array.mean()),
        "std": float(array.std(ddof=1)) if len(array) > 1 else 0.0,
        "p95": float(np.percentile(array, 95)),
        "maximum": float(array.max()),
        "empirical_upper_tail_p": float(
            (1 + np.sum(array >= observed)) / (len(array) + 1)
        ),
        "observed_gain": float(observed),
    }


def run_target_label_null(
    caches: dict[str, FoldCache],
    vip_rows: np.ndarray,
    pm_rows: np.ndarray,
    cell_count_quartiles: np.ndarray,
    *,
    reps: int,
    seed: int,
    max_rss_gib: float,
) -> np.ndarray:
    rng = np.random.default_rng(seed)
    values = np.empty(reps, dtype=np.float64)
    baseline_means = {
        held: float(cache.baseline_per_target.mean()) for held, cache in caches.items()
    }
    for rep in range(reps):
        permutation = stratified_permutation(cell_count_quartiles, rng)
        fold_gains = [
            float(
                candidate_pds_from_vip_rows(
                    cache, vip_rows[permutation], pm_rows
                ).mean()
                - baseline_means[held]
            )
            for held, cache in caches.items()
        ]
        values[rep] = float(np.mean(fold_gains))
        if rep % 100 == 0:
            enforce_limit(max_rss_gib, f"target-label null replicate {rep}")
    return values


def run_gene_coordinate_null(
    caches: dict[str, FoldCache],
    vip_rows: np.ndarray,
    pm_rows: np.ndarray,
    variance_deciles: np.ndarray,
    *,
    reps: int,
    seed: int,
    max_rss_gib: float,
) -> np.ndarray:
    rng = np.random.default_rng(seed)
    values = np.empty(reps, dtype=np.float64)
    baseline_means = {
        held: float(cache.baseline_per_target.mean()) for held, cache in caches.items()
    }
    modified_rows = {
        held: np.flatnonzero(cache.row_factors != 0) for held, cache in caches.items()
    }
    for lo in range(0, reps, GENE_NULL_CHUNK_REPS):
        hi = min(lo + GENE_NULL_CHUNK_REPS, reps)
        permutations = [
            stratified_permutation(variance_deciles, rng) for _ in range(lo, hi)
        ]
        stacked_delta = np.concatenate(
            [vip_rows[:, permutation] - pm_rows for permutation in permutations],
            axis=0,
        )
        fold_gains = np.empty((len(caches), hi - lo), dtype=np.float64)
        for fold_index, (held, cache) in enumerate(caches.items()):
            contributions = (stacked_delta @ cache.vip_projection).reshape(
                hi - lo,
                len(vip_rows),
                cache.baseline_logits.shape[1],
            )
            rows = modified_rows[held]
            for offset in range(hi - lo):
                logits = cache.baseline_logits.copy()
                logits[rows] += contributions[offset] * cache.row_factors[rows, None]
                fold_gains[fold_index, offset] = (
                    pds_per_target_from_logits(logits).mean() - baseline_means[held]
                )
        values[lo:hi] = fold_gains.mean(axis=0)
        enforce_limit(max_rss_gib, f"gene-coordinate null replicates {lo}:{hi}")
    return values


def low_cell_sensitivity(
    per_fold_delta: dict[str, np.ndarray],
    modified: np.ndarray,
    low_cell: np.ndarray,
    *,
    minimum_non_low_targets: int,
    low_cell_threshold: int = LOW_CELL_THRESHOLD,
) -> dict[str, object]:
    modified = np.asarray(modified, dtype=bool)
    low_cell = np.asarray(low_cell, dtype=bool)
    if modified.shape != low_cell.shape:
        raise ValueError("modified and low-cell masks differ in shape")
    cluster_delta = np.mean(
        np.stack([np.asarray(per_fold_delta[held]) for held in FOLDS], axis=1),
        axis=1,
    )
    low = modified & low_cell
    non_low = modified & ~low_cell
    statistically_defined = int(non_low.sum()) >= minimum_non_low_targets
    non_low_mean = float(cluster_delta[non_low].mean()) if non_low.any() else None
    full_gallery_reset = float(np.where(non_low, cluster_delta, 0.0).mean())
    return {
        "low_cell_threshold": int(low_cell_threshold),
        "modified_low_cell_targets": int(low.sum()),
        "modified_low_cell_mean_gain": (
            float(cluster_delta[low].mean()) if low.any() else None
        ),
        "modified_non_low_cell_targets": int(non_low.sum()),
        "minimum_non_low_targets": int(minimum_non_low_targets),
        "statistically_defined": statistically_defined,
        "modified_non_low_cell_mean_gain": non_low_mean,
        "full_gallery_gain_with_low_cell_modified_targets_reset": full_gallery_reset,
        "gain_persists": bool(
            statistically_defined and non_low_mean is not None and non_low_mean > 0.0
        ),
    }


def stopping_gate(
    modified_targets: int,
    fold_deltas: dict[str, float],
    macro_delta: float,
    bootstrap: dict[str, object],
    target_null: dict[str, object],
    gene_null: dict[str, object],
    low_cell: dict[str, object],
    *,
    minimum_modified_targets: int = MIN_MODIFIED_TARGETS,
) -> dict[str, object]:
    lower = float(bootstrap["ci95_percentile"][0])
    checks = {
        "at_least_70_modified_targets": modified_targets >= minimum_modified_targets,
        "both_folds_nonnegative": all(
            float(fold_deltas[held]) >= 0.0 for held in FOLDS
        ),
        "macro_gain_at_least_0_003": float(macro_delta) >= MACRO_GATE,
        "bootstrap_lower_bound_positive": lower > 0.0,
        "observed_gain_above_target_label_null_p95": (
            float(macro_delta) > float(target_null["p95"])
        ),
        "observed_gain_above_gene_coordinate_null_p95": (
            float(macro_delta) > float(gene_null["p95"])
        ),
        "gain_persists_without_low_cell_modified_targets": bool(
            low_cell["gain_persists"]
        ),
    }
    passed = all(checks.values())
    return {
        "thresholds": {
            "minimum_modified_targets": minimum_modified_targets,
            "minimum_fold_delta": 0.0,
            "minimum_macro_delta": MACRO_GATE,
            "minimum_bootstrap_lower_bound": 0.0,
            "target_label_null_p95": float(target_null["p95"]),
            "gene_coordinate_null_p95": float(gene_null["p95"]),
            "minimum_non_low_modified_targets": int(
                low_cell["minimum_non_low_targets"]
            ),
        },
        "checks": checks,
        "passes": passed,
        "action": (
            "eligible for one count-level confirmation"
            if passed
            else "stop; do not run count-level confirmation"
        ),
    }


def evaluate_fold(
    held: str,
    banks: dict[str, BASE.Bank],
    prepared: dict[str, BASE.PreparedSource],
    pm_slot: BASE.PreparedSource,
    viperturb: BASE.PreparedSource,
    vip_slot: BASE.PreparedSource,
    modified: np.ndarray,
    triplet_columns: np.ndarray,
    cd4: BASE.CD4Prepared,
    official: Sequence[str],
    panel: Sequence[str],
    targets: Sequence[str],
    vip_audit: VIPAudit,
    whitener_path: Path,
    *,
    expected_whitener_sha256: str,
    max_rss_gib: float,
    min_operator_genes: int,
    min_retrieval_genes: int,
) -> tuple[dict[str, object], list[dict[str, object]], np.ndarray, FoldCache]:
    companion_name = next(name for name in FOLDS if name != held)
    source_banks = (banks["k562_gwps"], banks[companion_name])
    operator, wcols = EXP064.load_frozen_operator(
        whitener_path,
        source_banks,
        official,
        set(panel),
        expected_sha256=expected_whitener_sha256,
    )
    if len(wcols) < min_operator_genes:
        raise ValueError(f"{held}: only {len(wcols)} operator genes")

    baseline_pre, baseline_votes = BASE.build_baseline(
        [pm_slot, prepared[companion_name]]
    )
    candidate_pre, candidate_votes = BASE.build_baseline(
        [vip_slot, prepared[companion_name]]
    )
    if not np.array_equal(baseline_votes, candidate_votes):
        raise AssertionError("VIPerturb changed source-vote counts")

    baseline_post, whitening = BASE.apply_operator(baseline_pre, operator, wcols)
    candidate_post = apply_operator_with_frozen_scale(
        candidate_pre,
        operator,
        wcols,
        float(whitening["applied_scale"]),
    )
    baseline_d, cd4_preserve = BASE.add_post_whiten_residual(
        baseline_post, cd4, CD4_FRACTION, wcols
    )
    candidate_d = add_cd4_with_frozen_scale(
        candidate_post,
        cd4,
        float(cd4_preserve["post_add_global_preserve_scale"]),
    )
    if (~modified).any():
        if not np.allclose(
            candidate_d[~modified],
            baseline_d[~modified],
            rtol=0.0,
            atol=1e-12,
        ):
            raise AssertionError("candidate changed a non-shared target")
        candidate_d[~modified] = baseline_d[~modified]
        if not np.array_equal(candidate_d[~modified], baseline_d[~modified]):
            raise AssertionError("non-shared arm-D rows were not preserved exactly")

    held_genes = set(banks[held].genes)
    retrieval = np.asarray(
        [index for index in wcols if official[index] in held_genes],
        dtype=np.int64,
    )
    if len(retrieval) < min_retrieval_genes:
        raise ValueError(f"{held}: only {len(retrieval)} retrieval genes")
    if set(panel) & {official[index] for index in retrieval}:
        raise AssertionError("panel target coordinate entered retrieval")

    truth = BASE.raw_truth(prepared[held])
    baseline_per_target = BASE.pds_per_target(baseline_d, truth, retrieval)
    candidate_per_target = BASE.pds_per_target(candidate_d, truth, retrieval)
    truth_retrieval = np.asarray(truth[:, retrieval], dtype=np.float64)
    truth_retrieval /= np.maximum(
        np.linalg.norm(truth_retrieval, axis=1, keepdims=True),
        1e-12,
    )
    baseline_logits = baseline_d[:, retrieval] @ truth_retrieval.T
    candidate_logits = candidate_d[:, retrieval] @ truth_retrieval.T
    np.testing.assert_array_equal(
        pds_per_target_from_logits(baseline_logits),
        baseline_per_target,
    )
    np.testing.assert_array_equal(
        pds_per_target_from_logits(candidate_logits),
        candidate_per_target,
    )

    wcol_position = {int(column): index for index, column in enumerate(wcols)}
    try:
        triplet_wpositions = np.asarray(
            [wcol_position[int(column)] for column in triplet_columns],
            dtype=np.int64,
        )
        retrieval_wpositions = np.asarray(
            [wcol_position[int(column)] for column in retrieval],
            dtype=np.int64,
        )
    except KeyError as error:
        raise AssertionError(
            "null coordinate is absent from frozen operator"
        ) from error
    truth_on_operator = np.zeros((len(targets), len(wcols)), dtype=np.float64)
    truth_on_operator[:, retrieval_wpositions] = truth_retrieval
    weights = operator.evals ** (-WHITEN_GAMMA / 2.0)
    transformed_truth = (
        (truth_on_operator @ operator.evecs) * weights
    ) @ operator.evecs.T
    vip_projection = np.asarray(
        transformed_truth[:, triplet_wpositions].T,
        dtype=np.float64,
    )
    vip_projection *= float(whitening["applied_scale"])
    vip_projection *= float(cd4_preserve["post_add_global_preserve_scale"])
    row_factors = np.zeros(len(targets), dtype=np.float64)
    if np.any(baseline_votes[modified] <= 0):
        raise AssertionError("modified target lacks a baseline source vote")
    row_factors[modified] = VIP_FRACTION / baseline_votes[modified]
    reconstructed = baseline_logits.copy()
    reconstructed[modified] += (
        (
            viperturb.effects[np.ix_(modified, triplet_columns)]
            - pm_slot.effects[np.ix_(modified, triplet_columns)]
        )
        @ vip_projection
        * row_factors[modified, None]
    )
    if not np.allclose(
        reconstructed[modified],
        candidate_logits[modified],
        rtol=1e-10,
        atol=1e-10,
    ):
        maximum_error = float(
            np.max(np.abs(reconstructed[modified] - candidate_logits[modified]))
        )
        raise AssertionError(
            f"{held}: cached null projection differs from direct pipeline "
            f"(max error {maximum_error:.3e})"
        )
    np.testing.assert_array_equal(
        pds_per_target_from_logits(reconstructed),
        candidate_per_target,
    )
    paired_delta = candidate_per_target - baseline_per_target
    baseline_pds = float(np.mean(baseline_per_target))
    candidate_pds = float(np.mean(candidate_per_target))
    delta = candidate_pds - baseline_pds

    per_target_rows: list[dict[str, object]] = []
    for target, is_modified, base, candidate, change in zip(
        targets,
        modified,
        baseline_per_target,
        candidate_per_target,
        paired_delta,
        strict=True,
    ):
        audit = vip_audit.target_rows.get(target, {})
        per_target_rows.append(
            {
                "held": held,
                "target_gene": target,
                "modified": bool(is_modified),
                "manifest_cell_count": audit.get("manifest_cell_count"),
                "observed_cell_count": audit.get("observed_cell_count"),
                "low_cell_manifest": audit.get("low_cell_manifest"),
                "baseline_arm_d_pds": float(base),
                "viperturb_arm_pds": float(candidate),
                "delta": float(change),
            }
        )
    unmodified_delta = paired_delta[~modified]
    if unmodified_delta.size and np.any(unmodified_delta != 0):
        raise AssertionError("non-shared targets changed PDS")
    modified_delta = paired_delta[modified]
    low_cell = np.asarray(
        [
            bool(
                vip_audit.target_rows.get(target, {}).get(
                    "low_cell_manifest",
                    False,
                )
                or vip_audit.target_rows.get(target, {}).get(
                    "low_cell_observed",
                    False,
                )
            )
            for target in targets
        ],
        dtype=bool,
    )
    modified_low = modified & low_cell
    modified_not_low = modified & ~low_cell
    fold = {
        "held": held,
        "companion_source": companion_name,
        "gallery_targets": len(targets),
        "modified_targets": int(modified.sum()),
        "baseline_arm_d_pds": baseline_pds,
        "viperturb_arm_pds": candidate_pds,
        "delta": delta,
        "modified_target_delta_mean": (
            float(modified_delta.mean()) if modified_delta.size else None
        ),
        "low_cell_audit": {
            "threshold_from_compact": vip_audit.metadata.get("low_cell_threshold"),
            "modified_low_cell_targets": int(modified_low.sum()),
            "modified_low_cell_delta_mean": (
                float(paired_delta[modified_low].mean()) if modified_low.any() else None
            ),
            "modified_non_low_cell_targets": int(modified_not_low.sum()),
            "modified_non_low_cell_delta_mean": (
                float(paired_delta[modified_not_low].mean())
                if modified_not_low.any()
                else None
            ),
        },
        "improved_targets": int(np.sum(paired_delta > 0)),
        "worsened_targets": int(np.sum(paired_delta < 0)),
        "tied_targets": int(np.sum(paired_delta == 0)),
        "operator_genes": len(wcols),
        "retrieval_genes": len(retrieval),
        "operator_columns_sha256": hashlib.sha256(wcols.tobytes()).hexdigest(),
        "retrieval_columns_sha256": hashlib.sha256(retrieval.tobytes()).hexdigest(),
        "whitener": {
            "path": str(whitener_path.resolve()),
            "sha256": sha256_file(whitener_path),
            "metadata_sha256": sha256_file(whitener_path.with_suffix(".json")),
            "baseline_whitening_preserve": whitening,
            "baseline_cd4_preserve": cd4_preserve,
            "candidate_reuses_baseline_scales": True,
        },
        "null_projection": {
            "method": (
                "exact PDS logits through the frozen symmetric whitening operator; "
                "prediction row norms are omitted because they cannot change within-row ranks"
            ),
            "triplet_input_genes": int(len(triplet_columns)),
            "direct_pipeline_rank_check": True,
        },
        "source_vote_counts": {
            str(value): int(np.sum(baseline_votes == value))
            for value in sorted(set(baseline_votes.tolist()))
        },
    }
    del (
        operator,
        baseline_pre,
        candidate_pre,
        baseline_post,
        candidate_post,
        baseline_d,
        candidate_d,
        truth,
        truth_retrieval,
        truth_on_operator,
        transformed_truth,
        candidate_logits,
        reconstructed,
    )
    gc.collect()
    enforce_limit(max_rss_gib, f"{held} fold completion")
    cache = FoldCache(
        held=held,
        baseline_logits=baseline_logits,
        vip_projection=vip_projection,
        row_factors=row_factors,
        baseline_per_target=baseline_per_target,
    )
    return fold, per_target_rows, paired_delta, cache


def write_csv(path: Path, rows: Sequence[dict[str, object]]) -> None:
    if not rows:
        return
    fields = sorted({key for row in rows for key in row})
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def bank_hashes(bank: BASE.Bank, *, include_b_ctrl: bool = True) -> dict[str, str]:
    names = ["delta.npy", "genes.txt", "perts.txt", "meta.json"]
    if include_b_ctrl:
        names.insert(1, "b_ctrl.npy")
    return {
        name: sha256_file(bank.path / name)
        for name in names
        if (bank.path / name).exists()
    }


def run(args: argparse.Namespace) -> dict[str, object]:
    if WHITEN_GAMMA != EXP064.WHITEN_GAMMA:
        raise AssertionError("experiment-064 whitening gamma changed")
    if args.output.exists():
        raise FileExistsError(f"refusing to overwrite {args.output}")
    panel = BASE.read_csv_column(args.panel_csv, "target_gene")
    official = BASE.read_csv_column(args.genes_csv, "gene_name")
    if len(panel) != args.expected_panel_targets:
        raise ValueError(
            f"panel has {len(panel)} targets, expected {args.expected_panel_targets}"
        )
    official_pos = {gene: index for index, gene in enumerate(official)}
    if any(target not in official_pos for target in panel):
        raise ValueError("a panel target is absent from the official gene axis")
    panel_positions = np.asarray(
        [official_pos[target] for target in panel], dtype=np.int64
    )
    panel_set = set(panel)

    banks = {
        "k562_gwps": BASE.Bank.load("k562_gwps", args.k562),
        "orion_hct116": BASE.Bank.load("orion_hct116", args.hct116),
        "orion_hek293t": BASE.Bank.load("orion_hek293t", args.hek293t),
    }
    perturbmatch = EXP064.load_perturbmatch(args.perturbmatch)
    cd4_banks = tuple(
        BASE.Bank.load(name, args.cd4_root / name, require_cd4_scale=True)
        for name in CD4_NAMES
    )
    targets = tuple(
        target
        for target in panel
        if target in banks["orion_hct116"].pert_pos
        and target in banks["orion_hek293t"].pert_pos
    )
    if len(targets) != args.expected_full_gallery_targets:
        raise ValueError(
            f"full common held gallery {len(targets)} != "
            f"{args.expected_full_gallery_targets}"
        )

    prepared = {
        name: BASE.prepare_main_source(
            bank,
            targets,
            official,
            panel_positions,
            chunk_rows=args.chunk_rows,
            max_rss_gib=args.max_rss_gib,
        )
        for name, bank in banks.items()
    }
    prepared_pm = BASE.prepare_main_source(
        perturbmatch,
        targets,
        official,
        panel_positions,
        chunk_rows=args.chunk_rows,
        max_rss_gib=args.max_rss_gib,
    )
    prepared_vip, panel_centered_vip, vip_audit = load_viperturb_compact(
        args.viperturb,
        targets,
        official,
        panel_positions,
        expected_coverage=args.expected_vip_coverage,
        expected_source_md5=args.expected_vip_source_md5,
        expected_manifest_md5=args.expected_vip_manifest_md5,
    )

    pm_columns = EXP064.shared_k562_pm_columns(
        banks["k562_gwps"], perturbmatch, official, panel_set
    )
    if len(pm_columns) < args.min_pm_shared_genes:
        raise ValueError(f"only {len(pm_columns)} shared K562/PerturbMatch genes")
    pm_slot, pm_audit = build_pm_slot_with_fallback(
        prepared["k562_gwps"], prepared_pm, pm_columns
    )
    vip_genes = set(read_lines(args.viperturb / "genes.txt"))
    triplet_columns = np.asarray(
        [index for index in pm_columns if official[int(index)] in vip_genes],
        dtype=np.int64,
    )
    if len(triplet_columns) < args.min_triplet_shared_genes:
        raise ValueError(f"only {len(triplet_columns)} shared K562/PM/VIP genes")
    vip_slot, modified, vip_slot_audit = build_vip_slot(
        prepared["k562_gwps"],
        prepared_pm,
        pm_slot,
        prepared_vip,
        triplet_columns,
        targets,
    )
    panel_vip_slot, panel_modified, panel_slot_audit = build_vip_slot(
        prepared["k562_gwps"],
        prepared_pm,
        pm_slot,
        panel_centered_vip,
        triplet_columns,
        targets,
    )
    if not np.array_equal(panel_modified, modified):
        raise AssertionError("panel-centered ablation changed target support")
    if (~modified).any() and not np.array_equal(
        panel_vip_slot.effects[~modified],
        pm_slot.effects[~modified],
    ):
        raise AssertionError("panel-centered ablation changed unsupported K562 rows")

    cd4 = freeze_complete_cd4_mean3(
        BASE.prepare_cd4_mean(
            cd4_banks,
            targets,
            official,
            panel_set,
            CD4_MODE,
            chunk_rows=args.chunk_rows,
            max_rss_gib=args.max_rss_gib,
        )
    )
    if np.any(cd4.mean_effect[:, panel_positions] != 0):
        raise AssertionError("CD4 residual contains a panel target coordinate")

    args.output.parent.mkdir(parents=True, exist_ok=True)
    folds: list[dict[str, object]] = []
    per_target_rows: list[dict[str, object]] = []
    per_fold_delta: dict[str, np.ndarray] = {}
    fold_caches: dict[str, FoldCache] = {}
    for held in FOLDS:
        expected_hash = (
            args.expected_hct_whitener_sha256
            if held == "orion_hct116"
            else args.expected_hek_whitener_sha256
        )
        fold, rows, paired, cache = evaluate_fold(
            held,
            banks,
            prepared,
            pm_slot,
            prepared_vip,
            vip_slot,
            modified,
            triplet_columns,
            cd4,
            official,
            panel,
            targets,
            vip_audit,
            args.whitener_dir / f"held_{held}.npz",
            expected_whitener_sha256=expected_hash,
            max_rss_gib=args.max_rss_gib,
            min_operator_genes=args.min_operator_genes,
            min_retrieval_genes=args.min_retrieval_genes,
        )
        folds.append(fold)
        per_target_rows.extend(rows)
        per_fold_delta[held] = paired
        fold_caches[held] = cache
        print(
            f"held={held} baseline_D={fold['baseline_arm_d_pds']:.9f} "
            f"vip={fold['viperturb_arm_pds']:.9f} "
            f"delta={fold['delta']:+.9f} "
            f"peak_rss={peak_rss_gib():.3f} GiB",
            flush=True,
        )

    baseline_macro = float(np.mean([fold["baseline_arm_d_pds"] for fold in folds]))
    candidate_macro = float(np.mean([fold["viperturb_arm_pds"] for fold in folds]))
    macro_delta = candidate_macro - baseline_macro
    bootstrap = paired_target_cluster_bootstrap(
        per_fold_delta, reps=args.bootstrap_reps, seed=args.bootstrap_seed
    )
    if not np.isclose(bootstrap["estimate"], macro_delta, atol=1e-14):
        raise AssertionError("bootstrap estimate does not match macro fold delta")

    modified_rows = np.flatnonzero(modified)
    vip_rows = prepared_vip.effects[np.ix_(modified, triplet_columns)]
    panel_vip_rows = panel_centered_vip.effects[np.ix_(modified, triplet_columns)]
    pm_rows = pm_slot.effects[np.ix_(modified, triplet_columns)]
    panel_ablation_folds: dict[str, float] = {}
    for held, cache in fold_caches.items():
        panel_score = float(
            candidate_pds_from_vip_rows(cache, panel_vip_rows, pm_rows).mean()
        )
        panel_ablation_folds[held] = panel_score
    panel_ablation_macro = float(np.mean(list(panel_ablation_folds.values())))
    panel_ablation = {
        "purpose": "diagnostic only; never eligible for shipping",
        "normalization": "center and median norm estimated on retained panel rows",
        "fold_pds": panel_ablation_folds,
        "macro_pds": panel_ablation_macro,
        "macro_delta_from_arm_d": panel_ablation_macro - baseline_macro,
        "modified_targets": int(panel_modified.sum()),
        "formula": panel_slot_audit["formula"],
    }

    cell_counts = np.asarray(
        [
            int(vip_audit.target_rows[targets[index]]["observed_cell_count"])
            for index in modified_rows
        ],
        dtype=np.float64,
    )
    cell_count_quartiles = balanced_rank_strata(cell_counts, 4)
    source_variances = vip_audit.source_variance[triplet_columns]
    variance_deciles = balanced_rank_strata(source_variances, 10)
    target_null_values = run_target_label_null(
        fold_caches,
        vip_rows,
        pm_rows,
        cell_count_quartiles,
        reps=args.null_reps,
        seed=args.target_null_seed,
        max_rss_gib=args.max_rss_gib,
    )
    gene_null_values = run_gene_coordinate_null(
        fold_caches,
        vip_rows,
        pm_rows,
        variance_deciles,
        reps=args.null_reps,
        seed=args.gene_null_seed,
        max_rss_gib=args.max_rss_gib,
    )
    target_null_path = args.output.parent / "target_label_null_macro_delta.npy"
    gene_null_path = args.output.parent / "gene_coordinate_null_macro_delta.npy"
    np.save(target_null_path, target_null_values)
    np.save(gene_null_path, gene_null_values)
    target_null = summarize_null(
        target_null_values,
        macro_delta,
        name="target_label",
        seed=args.target_null_seed,
        method=(
            "permute VIPerturb target labels only among modified targets, "
            "stratified by deterministic observed-cell-count quartile"
        ),
    )
    target_null["strata_counts"] = {
        str(label): int(np.sum(cell_count_quartiles == label)) for label in range(4)
    }
    target_null["stratification_value"] = "VIPerturb observed_cell_count"
    target_null["stratification_rule"] = (
        "stable rank split into four near-equal groups; target order breaks ties"
    )
    gene_null = summarize_null(
        gene_null_values,
        macro_delta,
        name="gene_coordinate",
        seed=args.gene_null_seed,
        method=(
            "apply one shared coordinate permutation per replicate within "
            "deterministic source-variance deciles, before the frozen operator"
        ),
    )
    gene_null["strata_counts"] = {
        str(label): int(np.sum(variance_deciles == label)) for label in range(10)
    }
    gene_null["stratification_value"] = (
        "per-gene population variance across all filtered source perturbations"
    )
    gene_null["stratification_rule"] = (
        "stable rank split into ten near-equal groups; official gene order breaks ties"
    )
    gene_null["chunk_reps"] = GENE_NULL_CHUNK_REPS
    low_cell_mask = np.asarray(
        [
            bool(
                vip_audit.target_rows.get(target, {}).get(
                    "low_cell_manifest",
                    False,
                )
                or vip_audit.target_rows.get(target, {}).get(
                    "low_cell_observed",
                    False,
                )
            )
            for target in targets
        ],
        dtype=bool,
    )
    low_cell = low_cell_sensitivity(
        per_fold_delta,
        modified,
        low_cell_mask,
        minimum_non_low_targets=args.minimum_non_low_targets,
        low_cell_threshold=int(vip_audit.metadata["low_cell_threshold"]),
    )
    fold_deltas = {str(fold["held"]): float(fold["delta"]) for fold in folds}
    gate = stopping_gate(
        int(modified.sum()),
        fold_deltas,
        macro_delta,
        bootstrap,
        target_null,
        gene_null,
        low_cell,
        minimum_modified_targets=args.minimum_modified_targets,
    )

    low_targets = {
        target
        for target, audit in vip_audit.target_rows.items()
        if bool(audit["low_cell_manifest"]) or bool(audit["low_cell_observed"])
    }
    modified_targets = set(vip_slot_audit["modified_targets"])
    result: dict[str, object] = {
        "experiment": "070-viperturb-k562",
        "status": "completed",
        "protocol": {
            "metric": "effect-level rank PDS",
            "design": "one frozen candidate versus full-gallery experiment-064 arm D",
            "valid_held_folds": list(FOLDS),
            "full_gallery": ("all challenge targets present in both held Orion donors"),
            "baseline": (
                "50:50 incumbent/PerturbMatch K562 slot with incumbent fallback, "
                "non-held Orion donor, candidate-independent gamma-0.65 whitener, "
                "then complete-three-condition CD4 gene_norm residual at 0.05"
            ),
            "candidate_k562_slot": vip_slot_audit["formula"],
            "candidate_target_policy": (
                "modify only targets present in incumbent K562, PerturbMatch, "
                "and filtered VIPerturb; preserve baseline elsewhere"
            ),
            "candidate_gene_policy": (
                "modify only shared official non-panel K562/PM/VIP coordinates"
            ),
            "viperturb_normalization": (
                "center and median centered norm over every perturbation in the "
                "filtered VIPerturb source object; retain panel rows afterward"
            ),
            "whitening_gamma": WHITEN_GAMMA,
            "operator_candidate_independent": True,
            "candidate_reuses_baseline_whitening_and_cd4_scales": True,
            "cd4_mode": CD4_MODE,
            "cd4_fraction": CD4_FRACTION,
            "count_emission": False,
            "cell_eval2": False,
            "submission": False,
            "null_reps": args.null_reps,
        },
        "axis": {
            "official_genes": len(official),
            "challenge_targets": len(panel),
            "gallery_targets": len(targets),
            "target_order": list(targets),
            "target_order_sha256": sha256_lines(targets),
            "panel_order_sha256": sha256_lines(panel),
            "k562_pm_shared_non_panel_genes": len(pm_columns),
            "k562_pm_vip_shared_non_panel_genes": len(triplet_columns),
            "triplet_columns_sha256": hashlib.sha256(
                triplet_columns.tobytes()
            ).hexdigest(),
        },
        "pm_slot": pm_audit,
        "vip_slot": vip_slot_audit,
        "viperturb": {
            "manifest_filtered_vcc_coverage": vip_audit.metadata[
                "manifest_filtered_vcc_coverage"
            ],
            "normalization_scope": vip_audit.metadata["normalization_scope"],
            "n_source_perturbations": vip_audit.metadata["n_source_perturbations"],
            "source_median_centered_norm": vip_audit.metadata[
                "source_median_centered_norm"
            ],
            "panel_median_centered_norm": vip_audit.metadata[
                "panel_median_centered_norm"
            ],
            "panel_median_centered_norm_stored": vip_audit.metadata[
                "panel_median_centered_norm_stored"
            ],
            "panel_median_centered_norm_validation": vip_audit.metadata[
                "panel_median_centered_norm_validation"
            ],
            "panel_median_centered_norm_legacy_r_recycling": (
                vip_audit.metadata.get("panel_median_centered_norm_legacy_r_recycling")
            ),
            "expression_source": vip_audit.metadata["expression_source"],
            "expression_kind": vip_audit.metadata["expression_kind"],
            "effect_definition": vip_audit.metadata["effect_definition"],
            "target_column": vip_audit.metadata["target_column"],
            "n_control_cells": vip_audit.metadata["n_control_cells"],
            "manifest_low_cell_targets": vip_audit.metadata[
                "manifest_low_cell_targets"
            ],
            "observed_low_cell_targets": vip_audit.metadata[
                "observed_low_cell_targets"
            ],
            "manifest_observed_count_mismatches": vip_audit.metadata[
                "manifest_observed_count_mismatches"
            ],
            "modified_low_cell_targets": sorted(modified_targets & low_targets),
            "compact_artifact_md5": vip_audit.compact_hashes,
        },
        "folds": folds,
        "macro": {
            "baseline_arm_d_pds": baseline_macro,
            "viperturb_arm_pds": candidate_macro,
            "delta": macro_delta,
            "bootstrap": bootstrap,
        },
        "panel_centered_ablation": panel_ablation,
        "nulls": {
            "target_label": target_null,
            "gene_coordinate": gene_null,
            "arrays": {
                "target_label": {
                    "filename": target_null_path.name,
                    "sha256": sha256_file(target_null_path),
                },
                "gene_coordinate": {
                    "filename": gene_null_path.name,
                    "sha256": sha256_file(gene_null_path),
                },
            },
        },
        "low_cell_sensitivity": low_cell,
        "stopping_gate": gate,
        "leakage_audit": {
            "held_truth_used_for_model_fitting": False,
            "held_truth_used_only_for_prespecified_evaluation": True,
            "held_truth_used_in_operator": False,
            "all_300_panel_coordinates_excluded_from_operator": True,
            "all_300_panel_coordinates_excluded_from_retrieval": True,
            "viperturb_used_in_operator": False,
            "source_wide_normalization_uses_held_truth": False,
            "null_strata_use_only_viperturb_source_metadata": True,
            "panel_centered_ablation_is_diagnostic_only": True,
            "candidate_count": 1,
            "hyperparameters_retuned": False,
            "null_definitions_and_seeds_prespecified": True,
            "k562_truth_evaluated": False,
            "unmodified_targets_preserved_exactly": True,
        },
        "input_hashes": {
            "run_py": sha256_file(Path(__file__)),
            "exp064_run_py": sha256_file(EXP064_PATH),
            "exp062_run_py": sha256_file(EXP064.BASE_PATH),
            "panel_csv": sha256_file(args.panel_csv),
            "genes_csv": sha256_file(args.genes_csv),
            "banks": {
                **{name: bank_hashes(bank) for name, bank in banks.items()},
                "perturbmatch": bank_hashes(perturbmatch, include_b_ctrl=False),
                **{bank.name: bank_hashes(bank) for bank in cd4_banks},
            },
            "viperturb_meta_json": sha256_file(args.viperturb / "meta.json"),
        },
        "infrastructure": {
            "provider": "AWS EC2",
            "instance_id": args.instance_id,
            "instance_type": args.instance_type,
            "market": args.market,
            "s3_prefix": args.s3_prefix,
        },
        "peak_rss_gib": peak_rss_gib(),
    }
    args.output.write_text(json.dumps(result, indent=2) + "\n")
    write_csv(
        args.output.parent / "fold_summary.csv",
        [
            {
                "held": fold["held"],
                "baseline_arm_d_pds": fold["baseline_arm_d_pds"],
                "viperturb_arm_pds": fold["viperturb_arm_pds"],
                "delta": fold["delta"],
                "modified_target_delta_mean": fold["modified_target_delta_mean"],
                "improved_targets": fold["improved_targets"],
                "worsened_targets": fold["worsened_targets"],
                "tied_targets": fold["tied_targets"],
            }
            for fold in folds
        ],
    )
    write_csv(args.output.parent / "per_target_pds.csv", per_target_rows)
    print(json.dumps({"macro": result["macro"], "gate": gate}, indent=2), flush=True)
    return result


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--k562", type=Path, required=True)
    parser.add_argument("--hct116", type=Path, required=True)
    parser.add_argument("--hek293t", type=Path, required=True)
    parser.add_argument("--perturbmatch", type=Path, required=True)
    parser.add_argument("--viperturb", type=Path, required=True)
    parser.add_argument("--cd4-root", type=Path, required=True)
    parser.add_argument("--panel-csv", type=Path, required=True)
    parser.add_argument("--genes-csv", type=Path, required=True)
    parser.add_argument("--whitener-dir", type=Path, required=True)
    parser.add_argument("--output", type=Path, default=HERE / "out/results.json")
    parser.add_argument("--chunk-rows", type=int, default=128)
    parser.add_argument("--max-rss-gib", type=float, default=55.0)
    parser.add_argument("--expected-panel-targets", type=int, default=300)
    parser.add_argument("--expected-full-gallery-targets", type=int, default=289)
    parser.add_argument(
        "--expected-vip-coverage", type=int, default=EXPECTED_VIP_COVERAGE
    )
    parser.add_argument(
        "--minimum-modified-targets", type=int, default=MIN_MODIFIED_TARGETS
    )
    parser.add_argument("--min-pm-shared-genes", type=int, default=6_900)
    parser.add_argument("--min-triplet-shared-genes", type=int, default=6_500)
    parser.add_argument("--min-operator-genes", type=int, default=17_000)
    parser.add_argument("--min-retrieval-genes", type=int, default=17_000)
    parser.add_argument("--bootstrap-reps", type=int, default=BOOTSTRAP_REPS)
    parser.add_argument("--bootstrap-seed", type=int, default=BOOTSTRAP_SEED)
    parser.add_argument("--null-reps", type=int, default=NULL_REPS)
    parser.add_argument("--target-null-seed", type=int, default=TARGET_NULL_SEED)
    parser.add_argument("--gene-null-seed", type=int, default=GENE_NULL_SEED)
    parser.add_argument(
        "--minimum-non-low-targets",
        type=int,
        default=MIN_NON_LOW_TARGETS,
    )
    parser.add_argument("--expected-vip-source-md5", default=VIP_SOURCE_MD5)
    parser.add_argument("--expected-vip-manifest-md5", default=VIP_MANIFEST_MD5)
    parser.add_argument(
        "--expected-hct-whitener-sha256",
        default=EXPECTED_WHITENER_SHA256["orion_hct116"],
    )
    parser.add_argument(
        "--expected-hek-whitener-sha256",
        default=EXPECTED_WHITENER_SHA256["orion_hek293t"],
    )
    parser.add_argument("--instance-id", default="unknown")
    parser.add_argument("--instance-type", default="unknown")
    parser.add_argument("--market", default="unknown")
    parser.add_argument("--s3-prefix", default="unknown")
    parser.add_argument("--guarded", action="store_true")
    return parser


def main() -> int:
    args = build_parser().parse_args()
    if not args.guarded:
        raise SystemExit("refusing unguarded real-data execution")
    if (
        args.chunk_rows < 1
        or args.max_rss_gib <= 0
        or args.bootstrap_reps < 1
        or args.null_reps < MIN_NULL_REPS
        or args.minimum_non_low_targets < 1
    ):
        raise SystemExit(
            f"chunk, memory, bootstrap, and sensitivity limits must be positive; "
            f"null reps must be at least {MIN_NULL_REPS}"
        )
    run(args)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
