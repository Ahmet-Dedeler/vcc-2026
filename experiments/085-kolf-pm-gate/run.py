#!/usr/bin/env python3
"""Frozen source-only PerturbMatch inclusion gate with one-time KOLF scoring.

The module keeps its import surface intentionally light.  Synthetic tests use
only NumPy and the pure gate/statistics helpers below.  The audited experiment
064 and production submission modules are imported lazily, after their hashes
and every non-KOLF input hash have passed.
"""

from __future__ import annotations

import argparse
import csv
import gc
import hashlib
import importlib.metadata
import importlib.util
import json
import os
import platform
import resource
import subprocess
import sys
import time
from dataclasses import asdict, dataclass
from pathlib import Path
from types import ModuleType
from typing import Callable, Iterable, Mapping, Sequence, TypeVar

os.environ.setdefault("OPENBLAS_NUM_THREADS", "4")
os.environ.setdefault("OMP_NUM_THREADS", "4")
os.environ.setdefault("MKL_NUM_THREADS", "4")
os.environ.setdefault("VECLIB_MAXIMUM_THREADS", "4")
os.environ.setdefault("NUMEXPR_NUM_THREADS", "4")

import numpy as np


ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
EXPERIMENT = "085-kolf-pm-gate"
S3_EXPERIMENT_PREFIX = "exp085-kolf-pm-gate"

FOLDS = ("orion_hct116", "orion_hek293t")
PANEL_TARGETS = 300
EXP064_TARGETS = 221
OFFICIAL_GENES = 18_533
OPERATOR_GENES = 17_942
KOLF_TARGETS = 273
KOLF_SCORE_GENES = 17_806

QUANTILE_FRACTIONS = (0.20, 0.30, 0.40, 0.50, 0.60)
TIE_TOLERANCE = 1e-12
TRAINING_MINIMUM_DELTA = 0.0015
TRAINING_PERMUTATIONS = 5_000
TRAINING_PERMUTATION_SEED = 8_502_026

CLUSTER_K = 5
CLUSTER_MINIMUM_COSINE = 0.80
MINIMUM_KOLF_CLUSTERS = 8
KOLF_MINIMUM_DELTA = 0.0015
KOLF_BOOTSTRAPS = 20_000
KOLF_BOOTSTRAP_SEED = 8_512_026
KOLF_ROUTING_PERMUTATIONS = 20_000
KOLF_ROUTING_SEED = 8_522_026

EXPECTED_PANEL_ORDER_SHA256 = (
    "cd3d00f1cd5da660609ad41d85b53e2d5557f7f5edb3e1f64559f78f7fb3cbc4"
)
EXPECTED_EXP064_TARGET_ORDER_SHA256 = (
    "59e231eeae34ac02ac336bef4828bec6bba373224051af78ddb49914c3ec383f"
)
EXPECTED_KOLF_TARGET_ORDER_SHA256 = (
    "b366592072d1e0da04aaabafbcc78882d34c97ff26eb2fca5fd9a9af0084f22f"
)
EXPECTED_KOLF_SCORE_AXIS_SHA256 = (
    "280aab34a0e86f27db8ade213e000f9066c8d5f9b89edd3ea966aba061ba9ae9"
)
EXPECTED_PRODUCTION_OPERATOR_SHA256 = (
    "5608e807230dea04f93f5e804ed26e8f7285d9487d16a2791283793ffb80fed3"
)
EXPECTED_HELD_OPERATOR_SHA256 = {
    "orion_hct116": (
        "7cecbf7f36e40d22a6bea1ff3adc86a8788b329def397044fb5b9ba675247cbc"
    ),
    "orion_hek293t": (
        "97e486e0d1a805e1f24a622f1cefbfb0917e9282c5d0491c4ddb41c8475f2af1"
    ),
}
EXPECTED_PM_SUPPORTED_TARGETS = 248
EXPECTED_CD4_COMPLETE_TARGETS = 276
EXPECTED_RUNTIME_PACKAGES = {
    "numpy": "2.5.2",
    "scipy": "1.18.0",
    "anndata": "0.13.2",
    "h5py": "3.16.0",
    "pandas": "3.0.5",
}

T = TypeVar("T")


@dataclass(frozen=True)
class GateRule:
    """One member of the frozen gate family."""

    direction: str
    quantile: float | None

    @property
    def identity(self) -> bool:
        return self.direction == "identity"

    @property
    def name(self) -> str:
        if self.identity:
            return "always_D"
        return f"{self.direction}:q={self.quantile:.2f}"


@dataclass(frozen=True)
class HeldTruth:
    values: np.ndarray
    targets: tuple[str, ...]
    score_genes: tuple[str, ...]
    panel_rows: np.ndarray
    official_columns: np.ndarray
    audit: dict[str, object]


def peak_rss_gib() -> float:
    value = float(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
    return value / 1024**3 if sys.platform == "darwin" else value / 1024**2


def memory_checkpoint(stage: str, limit_gib: float) -> None:
    peak = peak_rss_gib()
    print(f"[rss] {stage}: {peak:.3f} GiB", flush=True)
    if limit_gib > 0 and peak >= limit_gib:
        raise MemoryError(
            f"{stage}: peak RSS {peak:.3f} GiB reached {limit_gib:.3f} GiB"
        )


def sha256_file(path: Path, block_bytes: int = 8 << 20) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while block := handle.read(block_bytes):
            digest.update(block)
    return digest.hexdigest()


def sha256_array(values: np.ndarray, dtype: str = "<f8") -> str:
    array = np.ascontiguousarray(values, dtype=dtype)
    return hashlib.sha256(array.tobytes()).hexdigest()


def sha256_labels(values: Iterable[str]) -> str:
    """Hash a newline-terminated text axis, matching the pinned text files."""
    return hashlib.sha256(("\n".join(values) + "\n").encode()).hexdigest()


def sha256_joined(values: Iterable[str]) -> str:
    """Hash labels without a terminal newline, matching experiment 064."""
    return hashlib.sha256("\n".join(values).encode()).hexdigest()


def read_lines(path: Path) -> tuple[str, ...]:
    values = tuple(
        line.strip() for line in path.read_text().splitlines() if line.strip()
    )
    if not values or len(values) != len(set(values)):
        raise ValueError(f"{path}: axis is empty or contains duplicate labels")
    return values


def read_csv_column(path: Path, column: str) -> tuple[str, ...]:
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle)
        if reader.fieldnames is None or column not in reader.fieldnames:
            raise ValueError(f"{path}: required column {column!r} is missing")
        values = tuple(
            str(row[column]).strip()
            for row in reader
            if row.get(column) and str(row[column]).strip()
        )
    if not values or len(values) != len(set(values)):
        raise ValueError(f"{path}: {column!r} is empty or contains duplicates")
    return values


def read_csv_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle)
        if reader.fieldnames is None:
            raise ValueError(f"{path}: CSV has no header")
        return list(reader)


def atomic_write_text(path: Path, text: str) -> None:
    temporary = path.with_name(f".{path.name}.partial-{os.getpid()}")
    if path.exists() or temporary.exists():
        raise FileExistsError(f"refusing to overwrite {path} or {temporary}")
    try:
        temporary.write_text(text)
        os.replace(temporary, path)
    except BaseException:
        temporary.unlink(missing_ok=True)
        raise


def atomic_save_npy(path: Path, values: np.ndarray) -> None:
    temporary = path.with_name(f".{path.name}.partial-{os.getpid()}")
    if path.exists() or temporary.exists():
        raise FileExistsError(f"refusing to overwrite {path} or {temporary}")
    try:
        with temporary.open("xb") as handle:
            np.save(handle, values, allow_pickle=False)
        os.replace(temporary, path)
    except BaseException:
        temporary.unlink(missing_ok=True)
        raise


def write_csv(path: Path, rows: Sequence[Mapping[str, object]]) -> None:
    if not rows:
        raise ValueError("refusing to write an empty CSV")
    temporary = path.with_name(f".{path.name}.partial-{os.getpid()}")
    if path.exists() or temporary.exists():
        raise FileExistsError(f"refusing to overwrite {path} or {temporary}")
    try:
        fields = list(rows[0])
        with temporary.open("x", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=fields)
            writer.writeheader()
            writer.writerows(rows)
        os.replace(temporary, path)
    except BaseException:
        temporary.unlink(missing_ok=True)
        raise


def ensure_output_directory(path: Path) -> Path:
    output = path.expanduser().resolve()
    if HERE.resolve() not in (output, *output.parents):
        raise ValueError("output directory must stay inside experiment 085")
    if output.exists() and any(output.iterdir()):
        raise FileExistsError(f"refusing to reuse nonempty output directory {output}")
    output.mkdir(parents=True, exist_ok=True)
    return output


def load_module(name: str, path: Path) -> ModuleType:
    spec = importlib.util.spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise ImportError(f"cannot import frozen implementation {path}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


def verify_hash(path: Path, expected: str, label: str) -> str:
    if not path.is_file():
        raise FileNotFoundError(f"{label}: missing file {path}")
    observed = sha256_file(path)
    if observed != expected:
        raise ValueError(f"{label}: SHA-256 {observed} != frozen {expected}")
    return observed


def verify_named_files(
    root: Path,
    expected: Mapping[str, str],
    label: str,
) -> dict[str, str]:
    return {
        name: verify_hash(root / name, digest, f"{label}/{name}")
        for name, digest in expected.items()
    }


def validate_protocol(path: Path) -> dict[str, object]:
    protocol = json.loads(path.read_text())
    observed = {
        "experiment": protocol.get("experiment"),
        "frozen": protocol.get("frozen_before_kolf_evaluation"),
        "targets": protocol.get("training_truth", {}).get("targets"),
        "quantiles": tuple(protocol.get("gate", {}).get("quantile_fractions", ())),
        "production_rank_scope": protocol.get("gate", {}).get("production_rank_scope"),
        "tie_tolerance": protocol.get("gate", {}).get("tie_tolerance"),
        "training_delta": protocol.get("training_gate", {}).get(
            "minimum_macro_raw_pds_delta"
        ),
        "training_permutations": protocol.get("training_gate", {}).get(
            "full_selection_joint_target_label_permutations"
        ),
        "training_seed": protocol.get("training_gate", {}).get("permutation_seed"),
        "cluster_k": protocol.get("source_clusters", {}).get("reciprocal_k"),
        "cluster_cosine": protocol.get("source_clusters", {}).get(
            "minimum_pairwise_cosine"
        ),
        "minimum_clusters": protocol.get("source_clusters", {}).get(
            "minimum_kolf_represented_clusters"
        ),
        "kolf_targets": protocol.get("held_truth", {}).get("targets"),
        "kolf_genes": protocol.get("held_truth", {}).get("score_genes"),
        "control_used": protocol.get("held_truth", {}).get("control_used"),
        "kolf_delta": protocol.get("kolf_promotion_gate", {}).get(
            "minimum_raw_pds_delta"
        ),
        "bootstraps": protocol.get("kolf_promotion_gate", {}).get(
            "cluster_bootstrap_replicates"
        ),
        "bootstrap_seed": protocol.get("kolf_promotion_gate", {}).get(
            "cluster_bootstrap_seed"
        ),
        "routing_permutations": protocol.get("kolf_promotion_gate", {}).get(
            "matched_routing_permutations"
        ),
        "routing_seed": protocol.get("kolf_promotion_gate", {}).get(
            "matched_routing_seed"
        ),
    }
    expected = {
        "experiment": EXPERIMENT,
        "frozen": True,
        "targets": EXP064_TARGETS,
        "quantiles": QUANTILE_FRACTIONS,
        "production_rank_scope": (
            "all 248 joint K562/PerturbMatch-supported targets inside the "
            "complete 300-target VCC panel before KOLF authentication; "
            "unsupported targets use sentinel -1 and exact D"
        ),
        "tie_tolerance": TIE_TOLERANCE,
        "training_delta": TRAINING_MINIMUM_DELTA,
        "training_permutations": TRAINING_PERMUTATIONS,
        "training_seed": TRAINING_PERMUTATION_SEED,
        "cluster_k": CLUSTER_K,
        "cluster_cosine": CLUSTER_MINIMUM_COSINE,
        "minimum_clusters": MINIMUM_KOLF_CLUSTERS,
        "kolf_targets": KOLF_TARGETS,
        "kolf_genes": KOLF_SCORE_GENES,
        "control_used": False,
        "kolf_delta": KOLF_MINIMUM_DELTA,
        "bootstraps": KOLF_BOOTSTRAPS,
        "bootstrap_seed": KOLF_BOOTSTRAP_SEED,
        "routing_permutations": KOLF_ROUTING_PERMUTATIONS,
        "routing_seed": KOLF_ROUTING_SEED,
    }
    if observed != expected:
        raise ValueError(
            "preregistered protocol and code constants differ:\n"
            + json.dumps({"expected": expected, "observed": observed}, indent=2)
        )
    return {"sha256": sha256_file(path), "protocol": protocol}


def load_artifact_manifest(path: Path) -> dict[str, object]:
    if path.resolve() != (HERE / "artifact_manifest.json").resolve():
        raise ValueError("artifact manifest must be experiment 085's pinned file")
    manifest = json.loads(path.read_text())
    if manifest.get("schema_version") != 1 or manifest.get("experiment") != EXPERIMENT:
        raise ValueError("invalid experiment-085 artifact manifest")
    verify_hash(
        HERE / "preregistered_protocol.json",
        manifest["frozen_protocol_sha256"],
        "experiment-085 frozen protocol",
    )
    verify_hash(
        HERE / "README.md",
        manifest["frozen_readme_sha256"],
        "experiment-085 frozen README",
    )
    expected_axes = {
        "panel_targets": PANEL_TARGETS,
        "experiment_064_targets": EXP064_TARGETS,
        "kolf_targets": KOLF_TARGETS,
        "official_genes": OFFICIAL_GENES,
        "operator_genes": OPERATOR_GENES,
        "kolf_score_genes": KOLF_SCORE_GENES,
        "panel_order_sha256": EXPECTED_PANEL_ORDER_SHA256,
        "experiment_064_target_order_sha256": (EXPECTED_EXP064_TARGET_ORDER_SHA256),
        "kolf_target_order_sha256": EXPECTED_KOLF_TARGET_ORDER_SHA256,
        "kolf_score_axis_sha256": EXPECTED_KOLF_SCORE_AXIS_SHA256,
    }
    if manifest.get("fixed_axes") != expected_axes:
        raise ValueError("artifact manifest fixed axes differ from code")
    return manifest


def verify_pre_kolf_inputs(
    args: argparse.Namespace,
    manifest: Mapping[str, object],
) -> dict[str, object]:
    """Authenticate everything except KOLF; never touch the KOLF directory."""
    pre = manifest["pre_kolf"]
    assert isinstance(pre, dict)
    observed: dict[str, object] = {}

    observed["implementation"] = {
        name: verify_hash(ROOT / name, digest, f"implementation/{name}")
        for name, digest in pre["implementation"].items()
    }
    observed["design_provenance"] = {
        name: verify_hash(ROOT / name, digest, f"design-provenance/{name}")
        for name, digest in pre["design_provenance"].items()
    }
    exp064_root = ROOT / "experiments/064-pm-cd4-composition"
    observed["experiment_064"] = verify_named_files(
        exp064_root,
        pre["experiment_064"],
        "experiment-064",
    )
    observed["challenge"] = {
        "gene_names.csv": verify_hash(
            args.genes_csv,
            pre["challenge"]["gene_names.csv"],
            "challenge/gene_names.csv",
        ),
        "pert_counts.csv": verify_hash(
            args.panel_csv,
            pre["challenge"]["pert_counts.csv"],
            "challenge/pert_counts.csv",
        ),
    }
    donor_paths = {
        "k562_gwps": args.k562,
        "orion_hct116": args.hct116,
        "orion_hek293t": args.hek293t,
    }
    observed["donors"] = {
        name: verify_named_files(donor_paths[name], expected, f"donor/{name}")
        for name, expected in pre["donors"].items()
    }
    observed["perturbmatch"] = verify_named_files(
        args.perturbmatch,
        pre["perturbmatch"],
        "perturbmatch",
    )
    observed["cd4"] = {
        name: verify_named_files(args.cd4_root / name, expected, f"cd4/{name}")
        for name, expected in pre["cd4"].items()
    }
    observed["held_operators"] = verify_named_files(
        args.held_whitener_dir,
        pre["held_operators"],
        "held-operator",
    )
    observed["production_operator"] = {
        "whitener_union.npz": verify_hash(
            args.production_whitener,
            pre["production_operator"]["whitener_union.npz"],
            "production-operator/whitener_union.npz",
        ),
        "whitener_union.json": verify_hash(
            args.production_whitener_metadata,
            pre["production_operator"]["whitener_union.json"],
            "production-operator/whitener_union.json",
        ),
    }
    return observed


def validate_exp064_result_contract(path: Path) -> dict[str, object]:
    result = json.loads(path.read_text())
    leakage = result.get("leakage_audit", {})
    checks = {
        "experiment": result.get("experiment") == "064-pm-cd4-composition",
        "status": result.get("status") == "completed",
        "targets": result.get("axis", {}).get("gallery_targets") == EXP064_TARGETS,
        "target_hash": result.get("axis", {}).get("target_order_sha256")
        == EXPECTED_EXP064_TARGET_ORDER_SHA256,
        "panel_excluded_main": bool(
            leakage.get("all_300_panel_coordinates_zeroed_in_main_predictions")
        ),
        "panel_excluded_cd4": bool(
            leakage.get("all_300_panel_coordinates_excluded_from_cd4")
        ),
        "panel_excluded_operators": bool(
            leakage.get("all_300_panel_coordinates_excluded_from_operators")
        ),
        "panel_excluded_retrieval": bool(
            leakage.get("all_300_panel_coordinates_excluded_from_retrieval")
        ),
        "operator_candidate_independent": bool(
            leakage.get("operator_candidate_independent")
        ),
        "truth_untransformed": leakage.get("truth_transformed") is False,
        "no_retuning": leakage.get("hyperparameters_retuned") is False,
    }
    if not all(checks.values()):
        raise ValueError(
            "experiment-064 frozen result contract failed:\n"
            + json.dumps(checks, indent=2)
        )
    target_order = tuple(result["axis"]["target_order"])
    if sha256_joined(target_order) != EXPECTED_EXP064_TARGET_ORDER_SHA256:
        raise ValueError("experiment-064 target-order content drifted")
    folds = {fold["held"]: fold for fold in result["folds"]}
    if set(folds) != set(FOLDS):
        raise ValueError("experiment-064 held folds changed")
    for held in FOLDS:
        fold = folds[held]
        if fold.get("operator_genes") != OPERATOR_GENES:
            raise ValueError(f"{held}: experiment-064 operator size drifted")
        if fold.get("retrieval_genes") != KOLF_SCORE_GENES:
            raise ValueError(f"{held}: experiment-064 retrieval size drifted")
    return result


def midrank_percentile(values: np.ndarray) -> np.ndarray:
    """Return exact-tie midranks `(average_rank - 0.5) / n`."""
    array = np.asarray(values, dtype=np.float64)
    if array.ndim != 1 or len(array) < 2 or not np.isfinite(array).all():
        raise ValueError("percentile input must be finite, one-dimensional, n>=2")
    order = np.argsort(array, kind="stable")
    ranks = np.empty(len(array), dtype=np.float64)
    start = 0
    while start < len(array):
        stop = start + 1
        while stop < len(array) and array[order[stop]] == array[order[start]]:
            stop += 1
        average_rank = 0.5 * ((start + 1) + stop)
        ranks[order[start:stop]] = average_rank
        start = stop
    return (ranks - 0.5) / len(array)


def angular_displacement(left: np.ndarray, right: np.ndarray) -> np.ndarray:
    left_array = np.asarray(left, dtype=np.float64)
    right_array = np.asarray(right, dtype=np.float64)
    if left_array.shape != right_array.shape or left_array.ndim != 2:
        raise ValueError("angular displacement requires equal 2-D matrices")
    left_norm = np.linalg.norm(left_array, axis=1)
    right_norm = np.linalg.norm(right_array, axis=1)
    if np.any(left_norm <= 1e-12) or np.any(right_norm <= 1e-12):
        raise ValueError("angular displacement encountered a zero-norm query")
    cosine = np.sum(left_array * right_array, axis=1) / (left_norm * right_norm)
    return 1.0 - np.clip(cosine, -1.0, 1.0)


def gate_rules() -> tuple[GateRule, ...]:
    rules = [GateRule("identity", None)]
    for quantile in QUANTILE_FRACTIONS:
        rules.append(GateRule("low", quantile))
        rules.append(GateRule("high", quantile))
    return tuple(rules)


def route_c(feature: np.ndarray, rule: GateRule) -> np.ndarray:
    values = np.asarray(feature, dtype=np.float64)
    if rule.identity:
        return np.zeros(len(values), dtype=bool)
    assert rule.quantile is not None
    if rule.direction == "low":
        return values <= rule.quantile
    if rule.direction == "high":
        return values >= 1.0 - rule.quantile
    raise ValueError(f"unknown gate direction {rule.direction!r}")


def production_routes(
    raw_displacement: np.ndarray,
    pm_supported: np.ndarray,
    rule: GateRule,
) -> tuple[np.ndarray, np.ndarray]:
    """Rank only the support-comparable subset and force all others to D."""
    raw = np.asarray(raw_displacement, dtype=np.float64)
    supported = np.asarray(pm_supported, dtype=bool)
    if raw.ndim != 1 or raw.shape != supported.shape:
        raise ValueError("production displacement/support arrays do not align")
    if int(supported.sum()) < 2:
        raise ValueError("production routing needs at least two PM-supported targets")
    feature = np.full(len(raw), -1.0, dtype=np.float64)
    feature[supported] = midrank_percentile(raw[supported])
    assignments = np.zeros(len(raw), dtype=bool)
    assignments[supported] = route_c(feature[supported], rule)
    return feature, assignments


def _tie_key(rule: GateRule, switches: int) -> tuple[int, float, int]:
    quantile = float("inf") if rule.quantile is None else rule.quantile
    direction = {"low": 0, "high": 1, "identity": 2}[rule.direction]
    return switches, quantile, direction


def select_gate(
    features: Mapping[str, np.ndarray],
    c_minus_d: Mapping[str, np.ndarray],
    *,
    include_candidates: bool = True,
) -> dict[str, object]:
    if tuple(features) != FOLDS or tuple(c_minus_d) != FOLDS:
        raise ValueError(f"gate folds must be ordered exactly as {FOLDS}")
    lengths = {len(np.asarray(features[held])) for held in FOLDS}
    lengths |= {len(np.asarray(c_minus_d[held])) for held in FOLDS}
    if lengths != {EXP064_TARGETS}:
        raise ValueError(f"gate arrays must all have {EXP064_TARGETS} rows")

    best: dict[str, object] | None = None
    candidates: list[dict[str, object]] = []
    for rule in gate_rules():
        fold_deltas: dict[str, float] = {}
        fold_switches: dict[str, int] = {}
        for held in FOLDS:
            assignments = route_c(features[held], rule)
            fold_switches[held] = int(assignments.sum())
            fold_deltas[held] = float(
                np.mean(np.where(assignments, c_minus_d[held], 0.0))
            )
        switches = int(sum(fold_switches.values()))
        macro = float(np.mean(list(fold_deltas.values())))
        eligible = rule.identity or min(fold_deltas.values()) >= -TIE_TOLERANCE
        candidate: dict[str, object] = {
            "rule": rule,
            "name": rule.name,
            "identity": rule.identity,
            "eligible": bool(eligible),
            "macro_delta": macro,
            "fold_deltas": fold_deltas,
            "fold_switches": fold_switches,
            "switches": switches,
        }
        if include_candidates:
            candidates.append(candidate)
        if not eligible:
            continue
        if best is None:
            best = candidate
            continue
        best_macro = float(best["macro_delta"])
        if macro > best_macro + TIE_TOLERANCE:
            best = candidate
        elif abs(macro - best_macro) <= TIE_TOLERANCE:
            if _tie_key(rule, switches) < _tie_key(best["rule"], int(best["switches"])):
                best = candidate
    if best is None:
        raise AssertionError("mandatory always-D fallback was not selectable")
    return {
        "selected": best,
        "candidates": candidates if include_candidates else None,
    }


def full_selection_permutation_null(
    features: Mapping[str, np.ndarray],
    c_minus_d: Mapping[str, np.ndarray],
    *,
    repetitions: int = TRAINING_PERMUTATIONS,
    seed: int = TRAINING_PERMUTATION_SEED,
) -> dict[str, object]:
    observed = select_gate(features, c_minus_d, include_candidates=False)
    observed_delta = float(observed["selected"]["macro_delta"])
    rng = np.random.default_rng(seed)
    draws = np.empty(repetitions, dtype=np.float64)
    identity_selections = 0
    for replicate in range(repetitions):
        permutation = rng.permutation(EXP064_TARGETS)
        permuted = {held: np.asarray(c_minus_d[held])[permutation] for held in FOLDS}
        selected = select_gate(features, permuted, include_candidates=False)["selected"]
        draws[replicate] = float(selected["macro_delta"])
        identity_selections += int(bool(selected["identity"]))
    return {
        "method": (
            "same target-label permutation applied jointly to the paired "
            "HCT116/HEK293T C-minus-D outcomes; rerun complete candidate "
            "eligibility, selection, and tie-breaking"
        ),
        "repetitions": repetitions,
        "seed": seed,
        "observed_selected_macro_delta": observed_delta,
        "mean": float(draws.mean()),
        "p95": float(np.quantile(draws, 0.95)),
        "p99": float(np.quantile(draws, 0.99)),
        "one_sided_p": float((1 + np.sum(draws >= observed_delta)) / (repetitions + 1)),
        "identity_selection_fraction": float(identity_selections / repetitions),
    }


def training_promotion_gate(
    selection: Mapping[str, object],
    permutation: Mapping[str, object],
) -> dict[str, object]:
    selected = selection["selected"]
    fold_deltas = selected["fold_deltas"]
    checks = {
        "selected_nonidentity": not bool(selected["identity"]),
        "macro_delta_at_least_0_0015": (
            float(selected["macro_delta"]) >= TRAINING_MINIMUM_DELTA
        ),
        "both_donor_deltas_nonnegative": (
            min(float(value) for value in fold_deltas.values()) >= -TIE_TOLERANCE
        ),
        "full_selection_permutation_p_at_most_0_05": (
            float(permutation["one_sided_p"]) <= 0.05
        ),
    }
    return {
        **checks,
        "passed": all(checks.values()),
        "minimum_macro_raw_pds_delta": TRAINING_MINIMUM_DELTA,
    }


def row_normalize(values: np.ndarray) -> np.ndarray:
    array = np.asarray(values, dtype=np.float64)
    if array.ndim != 2 or not np.isfinite(array).all():
        raise ValueError("row normalization requires a finite 2-D matrix")
    norms = np.linalg.norm(array, axis=1, keepdims=True)
    if np.any(norms <= 1e-12):
        raise ValueError("cannot normalize a zero-norm row")
    return array / norms


def pds_per_target(prediction: np.ndarray, truth: np.ndarray) -> np.ndarray:
    predicted = np.asarray(prediction, dtype=np.float64)
    observed = np.asarray(truth, dtype=np.float64)
    if predicted.shape != observed.shape:
        raise ValueError(f"PDS shape mismatch {predicted.shape} != {observed.shape}")
    if predicted.ndim != 2 or len(predicted) < 2:
        raise ValueError("PDS requires a 2-D gallery with at least two targets")
    similarities = row_normalize(predicted) @ row_normalize(observed).T
    own = np.diag(similarities)[:, None]
    ranks = (similarities > own).sum(axis=1).astype(np.float64)
    ranks += 0.5 * ((similarities == own).sum(axis=1) - 1)
    return 1.0 - ranks / (len(predicted) - 1)


def source_effect_clusters(
    query: np.ndarray,
    *,
    reciprocal_k: int = CLUSTER_K,
    cosine_threshold: float = CLUSTER_MINIMUM_COSINE,
) -> tuple[np.ndarray, dict[str, object]]:
    """Reciprocal-kNN connected components from source queries only."""
    unit = row_normalize(query)
    similarities = unit @ unit.T
    np.fill_diagonal(similarities, -np.inf)
    k = min(max(1, reciprocal_k), len(unit) - 1)
    neighbors = np.argsort(-similarities, axis=1, kind="stable")[:, :k]
    directed = np.zeros((len(unit), len(unit)), dtype=bool)
    directed[np.arange(len(unit))[:, None], neighbors] = True
    reciprocal = directed & directed.T & (similarities >= cosine_threshold)
    parent = np.arange(len(unit), dtype=np.int64)

    def find(value: int) -> int:
        while parent[value] != value:
            parent[value] = parent[parent[value]]
            value = int(parent[value])
        return value

    def union(left: int, right: int) -> None:
        root_left, root_right = find(left), find(right)
        if root_left != root_right:
            parent[max(root_left, root_right)] = min(root_left, root_right)

    left_rows, right_rows = np.where(np.triu(reciprocal, k=1))
    for left, right in zip(left_rows.tolist(), right_rows.tolist(), strict=True):
        union(left, right)
    roots = [find(index) for index in range(len(unit))]
    mapping = {root: index for index, root in enumerate(sorted(set(roots)))}
    labels = np.asarray([mapping[root] for root in roots], dtype=np.int64)
    sizes = np.bincount(labels)
    return labels, {
        "method": "source-only arm-D reciprocal-kNN connected components",
        "reciprocal_k": k,
        "minimum_pairwise_cosine": cosine_threshold,
        "clusters": int(len(sizes)),
        "singleton_clusters": int(np.sum(sizes == 1)),
        "largest_cluster": int(np.max(sizes)),
        "cluster_sizes_sha256": sha256_array(sizes, "<i8"),
    }


def source_cluster_bootstrap(
    deltas: np.ndarray,
    labels: np.ndarray,
    *,
    repetitions: int = KOLF_BOOTSTRAPS,
    seed: int = KOLF_BOOTSTRAP_SEED,
    minimum_clusters: int = MINIMUM_KOLF_CLUSTERS,
) -> dict[str, object]:
    values = np.asarray(deltas, dtype=np.float64)
    clusters = np.asarray(labels, dtype=np.int64)
    if values.ndim != 1 or values.shape != clusters.shape:
        raise ValueError("cluster bootstrap arrays do not align")
    unique = np.unique(clusters)
    if len(unique) < minimum_clusters:
        return {
            "method": (
                "resample frozen source-query connected components with "
                "replacement, preserving every selected cluster's complete rows"
            ),
            "valid": False,
            "invalid_reason": (
                f"KOLF has {len(unique)} represented source clusters; "
                f"require {minimum_clusters}"
            ),
            "repetitions": 0,
            "seed": seed,
            "clusters": int(len(unique)),
            "point_estimate": float(values.mean()),
            "lower_95": None,
            "upper_95": None,
            "probability_positive": None,
        }
    members = [np.flatnonzero(clusters == cluster) for cluster in unique]
    cluster_sums = np.asarray([values[index].sum() for index in members])
    cluster_sizes = np.asarray([len(index) for index in members], dtype=np.int64)
    rng = np.random.default_rng(seed)
    draws = np.empty(repetitions, dtype=np.float64)
    chunk = 1_000
    for lo in range(0, repetitions, chunk):
        hi = min(lo + chunk, repetitions)
        selected = rng.integers(0, len(unique), size=(hi - lo, len(unique)))
        draws[lo:hi] = cluster_sums[selected].sum(axis=1) / cluster_sizes[selected].sum(
            axis=1
        )
    return {
        "method": (
            "resample frozen source-query connected components with "
            "replacement, preserving every selected cluster's complete rows"
        ),
        "valid": True,
        "repetitions": repetitions,
        "seed": seed,
        "clusters": int(len(unique)),
        "point_estimate": float(values.mean()),
        "lower_95": float(np.quantile(draws, 0.025)),
        "upper_95": float(np.quantile(draws, 0.975)),
        "probability_positive": float(np.mean(draws > 0.0)),
    }


def matched_routing_permutation(
    c_minus_d: np.ndarray,
    assignments_c: np.ndarray,
    pm_supported: np.ndarray,
    *,
    repetitions: int = KOLF_ROUTING_PERMUTATIONS,
    seed: int = KOLF_ROUTING_SEED,
) -> dict[str, object]:
    gains = np.asarray(c_minus_d, dtype=np.float64)
    assignments = np.asarray(assignments_c, dtype=bool)
    supported = np.asarray(pm_supported, dtype=bool)
    if (
        gains.ndim != 1
        or gains.shape != assignments.shape
        or gains.shape != supported.shape
    ):
        raise ValueError("matched-routing arrays do not align")
    if np.any(assignments & ~supported):
        raise ValueError("a PM-unsupported target was routed to C")
    supported_rows = np.flatnonzero(supported)
    c_count = int(np.sum(assignments[supported_rows]))
    if c_count <= 0 or c_count >= len(supported_rows):
        return {
            "method": (
                "uniformly reallocate the frozen C-route count among "
                "PM-supported KOLF targets; PM-unsupported targets remain exact D"
            ),
            "valid": False,
            "invalid_reason": "matched-routing null requires both C and D routes",
            "repetitions": 0,
            "seed": seed,
            "pm_supported_targets": int(len(supported_rows)),
            "c_routes": c_count,
            "d_routes": int(len(supported_rows) - c_count),
            "observed_delta": float(np.mean(np.where(assignments, gains, 0.0))),
            "mean": None,
            "p95": None,
            "p99": None,
            "one_sided_p": 1.0,
        }
    observed = float(np.mean(np.where(assignments, gains, 0.0)))
    supported_gains = gains[supported_rows]
    rng = np.random.default_rng(seed)
    draws = np.empty(repetitions, dtype=np.float64)
    chunk = 1_000
    for lo in range(0, repetitions, chunk):
        hi = min(lo + chunk, repetitions)
        random_keys = rng.random((hi - lo, len(supported_rows)))
        chosen = np.argpartition(random_keys, c_count - 1, axis=1)[:, :c_count]
        draws[lo:hi] = supported_gains[chosen].sum(axis=1) / len(gains)
    return {
        "method": (
            "uniformly reallocate the frozen C-route count among PM-supported "
            "KOLF targets; PM-unsupported targets remain exact D"
        ),
        "valid": True,
        "repetitions": repetitions,
        "seed": seed,
        "pm_supported_targets": int(len(supported_rows)),
        "c_routes": c_count,
        "d_routes": int(len(supported_rows) - c_count),
        "observed_delta": observed,
        "mean": float(draws.mean()),
        "p95": float(np.quantile(draws, 0.95)),
        "p99": float(np.quantile(draws, 0.99)),
        "one_sided_p": float((1 + np.sum(draws >= observed)) / (repetitions + 1)),
    }


def source_strata(
    deltas: np.ndarray,
    source_counts: np.ndarray,
    *,
    complete_count: int = 3,
) -> dict[str, object]:
    values = np.asarray(deltas, dtype=np.float64)
    counts = np.asarray(source_counts, dtype=np.int64)
    if values.shape != counts.shape:
        raise ValueError("source strata arrays do not align")
    masks = {
        "complete": counts == complete_count,
        "partial": (counts > 0) & (counts < complete_count),
        "direct": counts > 0,
        "missing": counts == 0,
    }
    return {
        name: {
            "targets": int(mask.sum()),
            "mean_delta": float(values[mask].mean()) if np.any(mask) else None,
        }
        for name, mask in masks.items()
    }


def execute_after_freeze(
    training_gate_passed: bool,
    freeze_builder: Callable[[], Mapping[str, object]],
    truth_loader: Callable[[Mapping[str, object]], T],
) -> tuple[Mapping[str, object] | None, T | None]:
    """Executable holdout boundary used by both the real path and tests."""
    if not training_gate_passed:
        return None, None
    frozen = freeze_builder()
    if not bool(frozen.get("sealed_before_truth")):
        raise ValueError("pre-KOLF artifact bundle is not sealed")
    return frozen, truth_loader(frozen)


def pinned_exp064_rows(
    path: Path,
    targets: Sequence[str],
) -> dict[str, dict[str, dict[str, float]]]:
    rows = read_csv_rows(path)
    if len(rows) != len(FOLDS) * len(targets):
        raise ValueError("experiment-064 per-target CSV has the wrong row count")
    output: dict[str, dict[str, dict[str, float]]] = {held: {} for held in FOLDS}
    for row in rows:
        held = row["held"]
        target = row["target_gene"]
        if held not in output or target in output[held]:
            raise ValueError("experiment-064 per-target key is invalid or duplicated")
        output[held][target] = {
            "C": float(row["pds_C"]),
            "D": float(row["pds_D"]),
        }
    for held in FOLDS:
        if tuple(output[held]) != tuple(targets):
            raise ValueError(f"{held}: experiment-064 target row order drifted")
    return output


def reconstruct_exp064_training(
    args: argparse.Namespace,
    exp064_result: Mapping[str, object],
) -> dict[str, object]:
    """Rebuild exact 064 C/D arms and return only source features/outcomes."""
    module = load_module(
        "exp085_exp064",
        ROOT / "experiments/064-pm-cd4-composition/run.py",
    )
    base = module.BASE
    panel = base.read_csv_column(args.panel_csv, "target_gene")
    official = base.read_csv_column(args.genes_csv, "gene_name")
    if len(panel) != PANEL_TARGETS or len(official) != OFFICIAL_GENES:
        raise ValueError("challenge axes have unexpected dimensions")
    official_pos = {gene: index for index, gene in enumerate(official)}
    panel_positions = np.asarray(
        [official_pos[target] for target in panel], dtype=np.int64
    )
    panel_set = set(panel)

    banks = {
        "k562_gwps": base.Bank.load("k562_gwps", args.k562),
        "orion_hct116": base.Bank.load("orion_hct116", args.hct116),
        "orion_hek293t": base.Bank.load("orion_hek293t", args.hek293t),
    }
    perturbmatch = module.load_perturbmatch(args.perturbmatch)
    cd4_banks = tuple(
        base.Bank.load(
            name,
            args.cd4_root / name,
            require_cd4_scale=True,
        )
        for name in module.CD4_NAMES
    )
    all_banks = [*banks.values(), perturbmatch, *cd4_banks]
    targets = tuple(
        target for target in panel if all(target in bank.pert_pos for bank in all_banks)
    )
    if len(targets) != EXP064_TARGETS:
        raise ValueError(f"experiment-064 gallery {len(targets)} != 221")
    if tuple(exp064_result["axis"]["target_order"]) != targets:
        raise ValueError("reconstructed experiment-064 target order differs")

    prepared = {
        name: base.prepare_main_source(
            bank,
            targets,
            official,
            panel_positions,
            chunk_rows=args.chunk_rows,
            max_rss_gib=args.max_rss_gib,
        )
        for name, bank in banks.items()
    }
    prepared_pm = base.prepare_main_source(
        perturbmatch,
        targets,
        official,
        panel_positions,
        chunk_rows=args.chunk_rows,
        max_rss_gib=args.max_rss_gib,
    )
    for source in [*prepared.values(), prepared_pm]:
        if np.any(source.effects[:, panel_positions] != 0.0):
            raise AssertionError("a panel coordinate survived source preparation")
    shared_columns = module.shared_k562_pm_columns(
        banks["k562_gwps"],
        perturbmatch,
        official,
        panel_set,
    )
    pm_slot = module.build_pm_slot(prepared["k562_gwps"], prepared_pm, shared_columns)
    cd4 = base.prepare_cd4_mean(
        cd4_banks,
        targets,
        official,
        panel_set,
        module.CD4_MODE,
        chunk_rows=args.chunk_rows,
        max_rss_gib=args.max_rss_gib,
    )
    if np.any(cd4.mean_effect[:, panel_positions] != 0.0):
        raise AssertionError("a panel coordinate entered the CD4 residual")
    pinned = pinned_exp064_rows(args.exp064_per_target, targets)

    features: dict[str, np.ndarray] = {}
    gains: dict[str, np.ndarray] = {}
    fold_audit: dict[str, object] = {}
    for held in FOLDS:
        companion = next(name for name in FOLDS if name != held)
        source_banks = (banks["k562_gwps"], banks[companion])
        operator, wcols = module.load_frozen_operator(
            args.held_whitener_dir / f"held_{held}.npz",
            source_banks,
            official,
            panel_set,
            expected_sha256=EXPECTED_HELD_OPERATOR_SHA256[held],
        )
        if len(wcols) != OPERATOR_GENES:
            raise ValueError(f"{held}: operator axis is not {OPERATOR_GENES}")
        if np.intersect1d(wcols, panel_positions).size:
            raise AssertionError(f"{held}: panel coordinate entered operator")

        baseline_a, votes_a = base.build_baseline(
            [prepared["k562_gwps"], prepared[companion]]
        )
        baseline_b, votes_b = base.build_baseline([pm_slot, prepared[companion]])
        if not np.array_equal(votes_a, votes_b):
            raise AssertionError("PerturbMatch changed experiment-064 source votes")
        arm_a, _ = base.apply_operator(baseline_a, operator, wcols)
        arm_b, _ = base.apply_operator(baseline_b, operator, wcols)
        arm_c, _ = base.add_post_whiten_residual(arm_a, cd4, module.CD4_FRACTION, wcols)
        arm_d, _ = base.add_post_whiten_residual(arm_b, cd4, module.CD4_FRACTION, wcols)
        if np.any(arm_c[:, panel_positions] != 0.0) or np.any(
            arm_d[:, panel_positions] != 0.0
        ):
            raise AssertionError(f"{held}: panel coordinate entered frozen arms")

        held_genes = set(banks[held].genes)
        retrieval = np.asarray(
            [index for index in wcols if official[index] in held_genes],
            dtype=np.int64,
        )
        if len(retrieval) != KOLF_SCORE_GENES:
            raise ValueError(f"{held}: retrieval axis is not 17,806 genes")
        truth = base.raw_truth(prepared[held])
        pds_c = base.pds_per_target(arm_c, truth, retrieval)
        pds_d = base.pds_per_target(arm_d, truth, retrieval)
        pinned_c = np.asarray([pinned[held][target]["C"] for target in targets])
        pinned_d = np.asarray([pinned[held][target]["D"] for target in targets])
        max_error_c = float(np.max(np.abs(pds_c - pinned_c)))
        max_error_d = float(np.max(np.abs(pds_d - pinned_d)))
        if max(max_error_c, max_error_d) > 1e-12:
            raise ValueError(
                f"{held}: arm C/D reproduction error "
                f"{max_error_c:.3g}/{max_error_d:.3g}"
            )

        raw_feature = angular_displacement(arm_c[:, wcols], arm_d[:, wcols])
        features[held] = midrank_percentile(raw_feature)
        gains[held] = pds_c - pds_d
        fold_audit[held] = {
            "targets": len(targets),
            "operator_genes": len(wcols),
            "retrieval_genes": len(retrieval),
            "operator_columns_sha256": hashlib.sha256(
                np.ascontiguousarray(wcols, dtype="<i8").tobytes()
            ).hexdigest(),
            "raw_feature_sha256": sha256_array(raw_feature),
            "percentile_feature_sha256": sha256_array(features[held]),
            "c_minus_d_sha256": sha256_array(gains[held]),
            "maximum_arm_c_pds_error": max_error_c,
            "maximum_arm_d_pds_error": max_error_d,
            "all_panel_coordinates_zero": True,
        }
        del (
            operator,
            baseline_a,
            baseline_b,
            arm_a,
            arm_b,
            arm_c,
            arm_d,
            truth,
            pds_c,
            pds_d,
        )
        gc.collect()
        memory_checkpoint(f"experiment-064 {held} reproduced", args.max_rss_gib)

    del banks, perturbmatch, cd4_banks, prepared, prepared_pm, pm_slot, cd4
    gc.collect()
    return {
        "targets": targets,
        "features": features,
        "c_minus_d": gains,
        "folds": fold_audit,
        "target_order_sha256": sha256_joined(targets),
    }


def build_production_freeze_inputs(
    args: argparse.Namespace,
    panel: Sequence[str],
    official: Sequence[str],
    selected_rule: GateRule,
) -> dict[str, object]:
    """Build exact three-donor production arms without opening KOLF."""
    production = load_module(
        "exp085_production",
        ROOT / "scripts/make_submission_pm_cd4.py",
    )
    production.assert_frozen_contract()
    donor_paths = {
        "k562_gwps": args.k562,
        "orion_hct116": args.hct116,
        "orion_hek293t": args.hek293t,
    }
    donors = {name: production.Donor.load(path) for name, path in donor_paths.items()}
    model = production.MultiDonorModel(
        donors,
        np.asarray(official),
        gene_norm=False,
        eb=False,
    )
    weights = production.donor_weights(model)
    pm = production.EXP061.prepare_perturbmatch(
        args.perturbmatch,
        tuple(official),
        chunk_rows=args.chunk_rows,
        max_rss_gib=args.max_rss_gib,
    )
    operator = production.load_production_operator(
        model,
        panel,
        args.production_whitener,
        args.production_whitener_metadata,
        expected_sha256=EXPECTED_PRODUCTION_OPERATOR_SHA256,
    )
    if len(operator.columns) != OPERATOR_GENES:
        raise ValueError("production operator is not the pinned 17,942-gene axis")
    official_pos = {gene: index for index, gene in enumerate(official)}
    panel_positions = np.asarray(
        [official_pos[target] for target in panel], dtype=np.int64
    )
    if np.intersect1d(operator.columns, panel_positions).size:
        raise AssertionError("a panel coordinate entered production whitening")
    cd4_banks, cd4 = production.prepare_frozen_cd4(
        args.cd4_root,
        panel,
        official,
        panel,
        chunk_rows=args.chunk_rows,
        max_rss_gib=args.max_rss_gib,
    )
    arm_c, audit_c = production.build_production_effects(
        model,
        panel,
        operator.whitener,
        operator.columns,
        panel_positions,
        pm=None,
        cd4=cd4,
        weights=weights,
    )
    arm_d, audit_d = production.build_production_effects(
        model,
        panel,
        operator.whitener,
        operator.columns,
        panel_positions,
        pm=pm,
        cd4=cd4,
        weights=weights,
    )
    if audit_c["perturbmatch"].get("enabled") is not False:
        raise AssertionError("production arm C unexpectedly enabled PerturbMatch")
    if int(audit_d["perturbmatch"]["both_supported_modified_count"]) != (
        EXPECTED_PM_SUPPORTED_TARGETS
    ):
        raise ValueError("production PerturbMatch support changed")
    if (
        int(audit_c["cd4"]["fully_supported_targets"]) != EXPECTED_CD4_COMPLETE_TARGETS
        or int(audit_d["cd4"]["fully_supported_targets"])
        != EXPECTED_CD4_COMPLETE_TARGETS
    ):
        raise ValueError("production complete-CD4 support changed")
    if (
        int(audit_c["cd4"]["panel_residual_nonzero_entries"]) != 0
        or int(audit_d["cd4"]["panel_residual_nonzero_entries"]) != 0
    ):
        raise AssertionError("panel coordinate entered production CD4 residual")

    raw_feature = angular_displacement(
        arm_c[:, operator.columns], arm_d[:, operator.columns]
    )
    pm_supported_targets = set(audit_d["perturbmatch"]["both_supported_modified"])
    pm_supported = np.asarray(
        [target in pm_supported_targets for target in panel], dtype=bool
    )
    feature, assignments_c = production_routes(raw_feature, pm_supported, selected_rule)
    gated = arm_d.copy()
    gated[assignments_c] = arm_c[assignments_c]
    if not np.array_equal(gated[~pm_supported], arm_d[~pm_supported]):
        raise AssertionError("PM-unsupported target is not exact arm D")

    source_counts = np.asarray(
        [
            sum(target in model.pert_index[name] for name in model.names)
            for target in panel
        ],
        dtype=np.uint8,
    )
    clusters, cluster_audit = source_effect_clusters(
        arm_d[:, operator.columns],
        reciprocal_k=CLUSTER_K,
        cosine_threshold=CLUSTER_MINIMUM_COSINE,
    )
    audit = {
        "donors": list(model.names),
        "weights": weights,
        "alpha": production.ALPHA,
        "beta": production.BETA,
        "gamma": production.WHITEN_GAMMA,
        "lambda": production.WHITEN_LAM,
        "operator_columns": len(operator.columns),
        "operator_columns_sha256": hashlib.sha256(
            np.ascontiguousarray(operator.columns, dtype="<i8").tobytes()
        ).hexdigest(),
        "operator_file_sha256": operator.file_sha256,
        "operator_metadata_sha256": operator.metadata_sha256,
        "arm_C": audit_c,
        "arm_D": audit_d,
        "pm_supported_targets": int(pm_supported.sum()),
        "c_routes_all_panel": int(assignments_c.sum()),
        "d_routes_all_panel": int((~assignments_c).sum()),
        "source_count_histogram": {
            str(value): int(np.sum(source_counts == value))
            for value in sorted(set(source_counts.tolist()))
        },
        "clusters": cluster_audit,
        "all_300_panel_coordinates_excluded_from_operator": True,
        "all_300_panel_coordinates_excluded_from_cd4_support": True,
        "kolf_opened": False,
    }
    del donors, model, pm, operator, cd4_banks, cd4
    gc.collect()
    memory_checkpoint("production arms and routing built", args.max_rss_gib)
    return {
        "arm_c": arm_c,
        "arm_d": arm_d,
        "gated": gated,
        "raw_feature": raw_feature,
        "feature": feature,
        "pm_supported": pm_supported,
        "assignments_c": assignments_c,
        "source_counts": source_counts,
        "clusters": clusters,
        "audit": audit,
    }


def freeze_production_artifacts(
    output: Path,
    panel: Sequence[str],
    selected_rule: GateRule,
    production: Mapping[str, object],
) -> dict[str, object]:
    freeze_dir = output / "pre_kolf"
    if freeze_dir.exists():
        raise FileExistsError(f"refusing to reuse freeze directory {freeze_dir}")
    freeze_dir.mkdir()
    arrays = {
        "arm_C.npy": np.asarray(production["arm_c"], dtype=np.float64),
        "arm_D.npy": np.asarray(production["arm_d"], dtype=np.float64),
        "gated.npy": np.asarray(production["gated"], dtype=np.float64),
        "raw_displacement.npy": np.asarray(production["raw_feature"], dtype=np.float64),
        "percentile_feature.npy": np.asarray(production["feature"], dtype=np.float64),
        "assignments_C.npy": np.asarray(production["assignments_c"], dtype=np.uint8),
        "pm_supported.npy": np.asarray(production["pm_supported"], dtype=np.uint8),
        "source_counts.npy": np.asarray(production["source_counts"], dtype=np.uint8),
        "source_clusters.npy": np.asarray(production["clusters"], dtype=np.int64),
    }
    expected_shape = (PANEL_TARGETS, OFFICIAL_GENES)
    for name in ("arm_C.npy", "arm_D.npy", "gated.npy"):
        if arrays[name].shape != expected_shape:
            raise ValueError(
                f"{name}: frozen shape {arrays[name].shape} != {expected_shape}"
            )
    for filename, values in arrays.items():
        atomic_save_npy(freeze_dir / filename, values)
    atomic_write_text(freeze_dir / "panel_targets.txt", "\n".join(panel) + "\n")
    assignment_rows = [
        {
            "target_gene": target,
            "raw_displacement": float(production["raw_feature"][index]),
            "percentile_feature": float(production["feature"][index]),
            "pm_supported": bool(production["pm_supported"][index]),
            "route": "C" if production["assignments_c"][index] else "D",
            "source_count": int(production["source_counts"][index]),
            "source_cluster": int(production["clusters"][index]),
        }
        for index, target in enumerate(panel)
    ]
    write_csv(freeze_dir / "assignments.csv", assignment_rows)

    file_hashes = {
        path.name: sha256_file(path)
        for path in sorted(freeze_dir.iterdir())
        if path.is_file()
    }
    freeze_manifest = {
        "schema_version": 1,
        "experiment": EXPERIMENT,
        "sealed_before_truth": True,
        "created_unix": time.time(),
        "targets": len(panel),
        "genes": int(np.asarray(production["arm_c"]).shape[1]),
        "selected_rule": asdict(selected_rule),
        "selected_rule_name": selected_rule.name,
        "panel_order_sha256": sha256_labels(panel),
        "files": file_hashes,
        "array_content_hashes": {
            "arm_C": sha256_array(production["arm_c"]),
            "arm_D": sha256_array(production["arm_d"]),
            "gated": sha256_array(production["gated"]),
            "assignments_C": sha256_array(production["assignments_c"], "<u1"),
            "source_clusters": sha256_array(production["clusters"], "<i8"),
        },
        "production_audit": production["audit"],
        "kolf_model_inputs_opened_before_seal": 0,
    }
    manifest_path = freeze_dir / "freeze_manifest.json"
    atomic_write_text(
        manifest_path,
        json.dumps(freeze_manifest, indent=2, sort_keys=True) + "\n",
    )
    return {
        **freeze_manifest,
        "freeze_manifest_path": str(manifest_path),
        "freeze_manifest_sha256": sha256_file(manifest_path),
    }


def parse_s3_uri(uri: str) -> tuple[str, str]:
    if not uri.startswith("s3://"):
        raise ValueError(f"not an S3 URI: {uri}")
    remainder = uri.removeprefix("s3://")
    if "/" not in remainder:
        raise ValueError(f"S3 URI has no object key: {uri}")
    bucket, key = remainder.split("/", 1)
    if not bucket or not key or "//" in key:
        raise ValueError(f"invalid S3 URI: {uri}")
    return bucket, key


def upload_durable_file(path: Path, uri: str, expected_sha256: str) -> None:
    """Atomically create one S3 object and verify its declared SHA metadata."""
    bucket, key = parse_s3_uri(uri)
    subprocess.run(
        [
            "aws",
            "s3api",
            "put-object",
            "--bucket",
            bucket,
            "--key",
            key,
            "--body",
            str(path),
            "--if-none-match",
            "*",
            "--checksum-algorithm",
            "SHA256",
            "--metadata",
            f"sha256={expected_sha256}",
            "--no-cli-pager",
        ],
        check=True,
        stdout=subprocess.DEVNULL,
    )
    observed = subprocess.run(
        [
            "aws",
            "s3api",
            "head-object",
            "--bucket",
            bucket,
            "--key",
            key,
            "--query",
            "Metadata.sha256",
            "--output",
            "text",
            "--no-cli-pager",
        ],
        check=True,
        capture_output=True,
        text=True,
    ).stdout.strip()
    if observed != expected_sha256:
        raise ValueError(
            f"durable object SHA-256 metadata {observed} != {expected_sha256}"
        )


def persist_pre_kolf_freeze(
    frozen: Mapping[str, object],
    durable_prefix: str,
    *,
    uploader: Callable[[Path, str, str], None] = upload_durable_file,
) -> dict[str, object]:
    """Persist every sealed pre-KOLF artifact before claiming truth access."""
    if not durable_prefix.endswith("/pre_kolf/"):
        raise ValueError("durable freeze URI must end with /pre_kolf/")
    _, key = parse_s3_uri(durable_prefix)
    expected_prefix = f"output/{S3_EXPERIMENT_PREFIX}/"
    if not key.startswith(expected_prefix):
        raise ValueError("durable freeze URI is outside experiment 085 output")

    manifest_path = Path(str(frozen["freeze_manifest_path"])).resolve()
    freeze_dir = manifest_path.parent
    if manifest_path.name != "freeze_manifest.json" or freeze_dir.name != "pre_kolf":
        raise ValueError("frozen manifest is not inside the pre_kolf directory")
    expected = dict(frozen["files"])
    expected[manifest_path.name] = str(frozen["freeze_manifest_sha256"])
    observed_names = {path.name for path in freeze_dir.iterdir() if path.is_file()}
    if observed_names != set(expected):
        raise ValueError(
            "pre-KOLF durable file set differs from freeze manifest: "
            f"expected={sorted(expected)}, observed={sorted(observed_names)}"
        )

    uploaded: dict[str, str] = {}
    for name in sorted(expected):
        path = freeze_dir / name
        digest = verify_hash(path, str(expected[name]), f"pre-KOLF freeze/{name}")
        uri = f"{durable_prefix}{name}"
        uploader(path, uri, digest)
        uploaded[name] = digest
    return {
        "uri": durable_prefix,
        "files": uploaded,
        "files_verified": len(uploaded),
        "freeze_manifest_sha256": str(frozen["freeze_manifest_sha256"]),
        "atomic_no_overwrite_uploads": True,
    }


def upload_durable_marker(path: Path, uri: str, expected_sha256: str) -> None:
    _, key = parse_s3_uri(uri)
    expected_key = f"output/{S3_EXPERIMENT_PREFIX}/KOLF_ACCESS_CONSUMED.json"
    if key != expected_key:
        raise ValueError(
            f"durable marker key must be the global {expected_key}, got {key}"
        )
    upload_durable_file(path, uri, expected_sha256)
    downloaded = subprocess.run(
        ["aws", "s3", "cp", uri, "-", "--only-show-errors"],
        check=True,
        capture_output=True,
    ).stdout
    observed = hashlib.sha256(downloaded).hexdigest()
    if observed != expected_sha256:
        raise ValueError(
            f"durable KOLF marker readback SHA-256 {observed} != {expected_sha256}"
        )


def mark_kolf_access(
    output: Path,
    frozen: Mapping[str, object],
    durable_freeze: Mapping[str, object],
    durable_uri: str,
    *,
    uploader: Callable[[Path, str, str], None] = upload_durable_marker,
) -> dict[str, object]:
    _, key = parse_s3_uri(durable_uri)
    expected_key = f"output/{S3_EXPERIMENT_PREFIX}/KOLF_ACCESS_CONSUMED.json"
    if key != expected_key:
        raise ValueError(
            f"durable marker key must be the global {expected_key}, got {key}"
        )
    if durable_freeze.get("freeze_manifest_sha256") != frozen.get(
        "freeze_manifest_sha256"
    ):
        raise ValueError("durable pre-KOLF freeze does not match local frozen bundle")
    marker = output / "KOLF_ACCESS_CONSUMED.json"
    record = {
        "experiment": EXPERIMENT,
        "reason": "one-time untouched KOLF truth evaluation",
        "created_unix": time.time(),
        "freeze_manifest_sha256": frozen["freeze_manifest_sha256"],
        "durable_pre_kolf_freeze": durable_freeze,
        "post_marker_refitting_allowed": False,
    }
    atomic_write_text(marker, json.dumps(record, indent=2, sort_keys=True) + "\n")
    digest = sha256_file(marker)
    uploader(marker, durable_uri, digest)
    return {
        **record,
        "path": str(marker),
        "sha256": digest,
        "durable_uri": durable_uri,
        "durable_sha256_verified_before_kolf_open": True,
    }


def verify_kolf_inputs(
    kolf_root: Path,
    manifest: Mapping[str, object],
) -> dict[str, str]:
    expected = manifest["post_freeze_kolf_truth_only"]
    assert isinstance(expected, dict)
    return verify_named_files(kolf_root, expected, "KOLF-held-truth")


def load_kolf_truth(
    args: argparse.Namespace,
    manifest: Mapping[str, object],
    official: Sequence[str],
    panel: Sequence[str],
    *,
    max_rss_gib: float,
) -> HeldTruth:
    """Authenticate and open KOLF only after the durable access marker exists."""
    marker = args.output_dir / "KOLF_ACCESS_CONSUMED.json"
    if not marker.is_file():
        raise ValueError("KOLF access marker is absent")
    hashes = verify_kolf_inputs(args.kolf, manifest)
    metadata = json.loads((args.kolf / "meta.json").read_text())
    provenance = json.loads((args.kolf / "provenance.json").read_text())
    if metadata.get("vcc_target_coverage") != KOLF_TARGETS:
        raise ValueError("KOLF target coverage changed")
    if metadata.get("vcc_gene_coverage") != 18_106:
        raise ValueError("KOLF official-gene coverage changed")
    if "delta = int8 source value * 0.04" not in str(
        metadata.get("effect_definition", "")
    ):
        raise ValueError("KOLF effect representation changed")
    control_definition = str(metadata.get("b_ctrl_definition", ""))
    if "all-ones compatibility placeholder" not in control_definition:
        raise ValueError("KOLF control placeholder contract changed")
    if provenance.get("source_repository") != metadata.get(
        "source_repository"
    ) or provenance.get("source_commit") != metadata.get("source_commit"):
        raise ValueError("KOLF provenance and metadata source identities differ")
    provenance_artifacts = provenance.get("artifacts", {})
    for filename, digest in hashes.items():
        if filename == "provenance.json":
            continue
        declared = provenance_artifacts.get(filename, {}).get("sha256")
        if declared is None:
            raise ValueError(f"KOLF provenance omits {filename}")
        if declared != digest:
            raise ValueError(f"KOLF provenance hash mismatch for {filename}")

    kolf_targets = read_lines(args.kolf / "perts.txt")
    kolf_genes = read_lines(args.kolf / "genes.txt")
    panel_set = set(panel)
    shared_targets = tuple(target for target in panel if target in set(kolf_targets))
    score_genes = tuple(
        gene for gene in official if gene in set(kolf_genes) and gene not in panel_set
    )
    if len(shared_targets) != KOLF_TARGETS or tuple(kolf_targets) != shared_targets:
        raise ValueError("KOLF target gallery is not the complete pinned panel subset")
    if len(score_genes) != KOLF_SCORE_GENES:
        raise ValueError("KOLF score axis is not 17,806 genes")
    if sha256_labels(panel) != EXPECTED_PANEL_ORDER_SHA256:
        raise ValueError("VCC panel order changed")
    if sha256_labels(shared_targets) != EXPECTED_KOLF_TARGET_ORDER_SHA256:
        raise ValueError("KOLF target order changed")
    if sha256_labels(score_genes) != EXPECTED_KOLF_SCORE_AXIS_SHA256:
        raise ValueError("KOLF score-gene order changed")
    if panel_set & set(score_genes):
        raise AssertionError("a VCC target coordinate entered KOLF scoring")

    panel_pos = {target: index for index, target in enumerate(panel)}
    official_pos = {gene: index for index, gene in enumerate(official)}
    kolf_gene_pos = {gene: index for index, gene in enumerate(kolf_genes)}
    panel_rows = np.asarray(
        [panel_pos[target] for target in shared_targets], dtype=np.int64
    )
    official_columns = np.asarray(
        [official_pos[gene] for gene in score_genes], dtype=np.int64
    )
    kolf_columns = np.asarray(
        [kolf_gene_pos[gene] for gene in score_genes], dtype=np.int64
    )

    delta = np.load(args.kolf / "delta.npy", mmap_mode="r", allow_pickle=False)
    if delta.shape != (KOLF_TARGETS, len(kolf_genes)):
        raise ValueError(f"unexpected KOLF delta shape {delta.shape}")
    raw = np.asarray(delta[:, kolf_columns], dtype=np.float64)
    values = raw - raw.mean(axis=0, keepdims=True)
    mmap = getattr(delta, "_mmap", None)
    if mmap is not None:
        mmap.close()
    del delta, raw
    gc.collect()
    if not np.isfinite(values).all():
        raise ValueError("KOLF truth contains non-finite values")
    memory_checkpoint("KOLF truth authenticated and opened", max_rss_gib)
    return HeldTruth(
        values=values,
        targets=shared_targets,
        score_genes=score_genes,
        panel_rows=panel_rows,
        official_columns=official_columns,
        audit={
            "role": "one-time untouched truth only",
            "source_repository": metadata.get("source_repository"),
            "source_commit": metadata.get("source_commit"),
            "provenance_source_commit": provenance.get("source_commit"),
            "representation": (
                "decoded raw quantized z-score, centered per gene over the "
                "complete 273-target compact gallery"
            ),
            "targets": len(shared_targets),
            "score_genes": len(score_genes),
            "target_order_sha256": sha256_labels(shared_targets),
            "score_axis_sha256": sha256_labels(score_genes),
            "hashes": hashes,
            "b_ctrl_hashed_for_provenance": True,
            "b_ctrl_loaded_or_used": False,
            "used_for_selection_or_fitting": False,
        },
    )


def evaluate_kolf_once(
    output: Path,
    held: HeldTruth,
    production: Mapping[str, object],
) -> dict[str, object]:
    rows = held.panel_rows
    columns = held.official_columns
    arm_c = np.asarray(production["arm_c"])[np.ix_(rows, columns)]
    arm_d = np.asarray(production["arm_d"])[np.ix_(rows, columns)]
    gated = np.asarray(production["gated"])[np.ix_(rows, columns)]
    if arm_c.shape != held.values.shape or arm_d.shape != held.values.shape:
        raise ValueError("production query and KOLF truth axes do not align")

    pds_c = pds_per_target(arm_c, held.values)
    pds_d = pds_per_target(arm_d, held.values)
    pds_gated = pds_per_target(gated, held.values)
    assignments = np.asarray(production["assignments_c"], dtype=bool)[rows]
    pm_supported = np.asarray(production["pm_supported"], dtype=bool)[rows]
    expected_gated = np.where(assignments, pds_c, pds_d)
    np.testing.assert_allclose(
        pds_gated,
        expected_gated,
        rtol=0.0,
        atol=1e-15,
        err_msg="gated PDS is not the frozen rowwise C/D selection",
    )
    if np.any(assignments & ~pm_supported):
        raise AssertionError("KOLF PM-unsupported target was routed to C")
    full_gated = np.asarray(production["gated"])
    full_d = np.asarray(production["arm_d"])
    full_unsupported = ~np.asarray(production["pm_supported"], dtype=bool)
    unsupported_exact = np.array_equal(
        full_gated[full_unsupported], full_d[full_unsupported]
    )
    if not unsupported_exact:
        raise AssertionError("PM-unsupported full-panel target is not exact D")

    per_target_delta = pds_gated - pds_d
    raw_delta = float(per_target_delta.mean())
    kolf_clusters = np.asarray(production["clusters"], dtype=np.int64)[rows]
    bootstrap = source_cluster_bootstrap(per_target_delta, kolf_clusters)
    routing_null = matched_routing_permutation(
        pds_c - pds_d,
        assignments,
        pm_supported,
    )
    source_counts = np.asarray(production["source_counts"], dtype=np.uint8)[rows]
    strata = source_strata(per_target_delta, source_counts)
    both_routes = bool(
        np.any(assignments[pm_supported]) and np.any(~assignments[pm_supported])
    )
    checks = {
        "complete_273_target_gallery": len(held.targets) == KOLF_TARGETS,
        "exact_17806_gene_score_axis": len(held.score_genes) == KOLF_SCORE_GENES,
        "raw_pds_delta_at_least_0_0015": raw_delta >= KOLF_MINIMUM_DELTA,
        "source_cluster_count_at_least_8": (
            int(bootstrap["clusters"]) >= MINIMUM_KOLF_CLUSTERS
        ),
        "source_cluster_bootstrap_valid": bool(bootstrap["valid"]),
        "cluster_bootstrap_lower_95_above_zero": (
            bool(bootstrap["valid"])
            and bootstrap["lower_95"] is not None
            and float(bootstrap["lower_95"]) > 0.0
        ),
        "matched_routing_null_valid": bool(routing_null["valid"]),
        "matched_routing_p_at_most_0_05": (float(routing_null["one_sided_p"]) <= 0.05),
        "both_C_and_D_used_among_pm_supported": both_routes,
        "pm_unsupported_targets_exact_D": unsupported_exact,
        "kolf_control_unused": held.audit["b_ctrl_loaded_or_used"] is False,
        "all_provenance_and_freeze_checks_passed": True,
    }
    promotion = {**checks, "passed": all(checks.values())}

    per_target_rows = []
    for index, target in enumerate(held.targets):
        count = int(source_counts[index])
        per_target_rows.append(
            {
                "target_gene": target,
                "source_cluster": int(kolf_clusters[index]),
                "source_count": count,
                "source_stratum": (
                    "complete" if count == 3 else "partial" if count > 0 else "missing"
                ),
                "direct_source": count > 0,
                "pm_supported": bool(pm_supported[index]),
                "route": "C" if assignments[index] else "D",
                "arm_C_pds": float(pds_c[index]),
                "arm_D_pds": float(pds_d[index]),
                "gated_pds": float(pds_gated[index]),
                "delta_vs_D": float(per_target_delta[index]),
            }
        )
    write_csv(output / "per_target_pds.csv", per_target_rows)
    return {
        "scores": {
            "arm_C_raw_pds": float(pds_c.mean()),
            "arm_D_raw_pds": float(pds_d.mean()),
            "gated_raw_pds": float(pds_gated.mean()),
            "gated_minus_D": raw_delta,
            "cluster_bootstrap": bootstrap,
            "matched_routing_null": routing_null,
        },
        "promotion_gate": promotion,
        "source_strata": strata,
        "routing": {
            "pm_supported_targets": int(pm_supported.sum()),
            "pm_unsupported_targets": int((~pm_supported).sum()),
            "C_routes_on_supported": int(np.sum(assignments & pm_supported)),
            "D_routes_on_supported": int(np.sum(~assignments & pm_supported)),
            "unsupported_exact_D": unsupported_exact,
        },
        "array_hashes": {
            "arm_C_eval": sha256_array(arm_c),
            "arm_D_eval": sha256_array(arm_d),
            "gated_eval": sha256_array(gated),
            "kolf_truth": sha256_array(held.values),
            "per_target_delta": sha256_array(per_target_delta),
        },
    }


def serialize_selection(selection: Mapping[str, object]) -> dict[str, object]:
    selected = dict(selection["selected"])
    selected["rule"] = asdict(selected["rule"])
    candidates = []
    for item in selection.get("candidates") or []:
        row = dict(item)
        row["rule"] = asdict(row["rule"])
        candidates.append(row)
    return {"selected": selected, "candidates": candidates}


def runtime_metadata() -> dict[str, object]:
    packages = {}
    for name in EXPECTED_RUNTIME_PACKAGES:
        try:
            packages[name] = importlib.metadata.version(name)
        except importlib.metadata.PackageNotFoundError:
            packages[name] = None
    return {
        "python": platform.python_version(),
        "platform": platform.platform(),
        "packages": packages,
    }


def verify_runtime_versions() -> dict[str, str]:
    observed = {
        name: importlib.metadata.version(name) for name in EXPECTED_RUNTIME_PACKAGES
    }
    if observed != EXPECTED_RUNTIME_PACKAGES:
        raise ValueError(
            "runtime package versions differ from the frozen environment:\n"
            + json.dumps(
                {
                    "expected": EXPECTED_RUNTIME_PACKAGES,
                    "observed": observed,
                },
                indent=2,
            )
        )
    if sys.version_info[:2] != (3, 12):
        raise ValueError(
            f"frozen runtime requires Python 3.12, got {platform.python_version()}"
        )
    return observed


def run_experiment(args: argparse.Namespace) -> dict[str, object]:
    output = ensure_output_directory(args.output_dir)
    args.output_dir = output
    frozen_runtime = verify_runtime_versions()
    protocol = validate_protocol(HERE / "preregistered_protocol.json")
    manifest = load_artifact_manifest(args.artifact_manifest)
    pre_hashes = verify_pre_kolf_inputs(args, manifest)
    exp064_result = validate_exp064_result_contract(args.exp064_results)
    memory_checkpoint("all non-KOLF inputs authenticated", args.max_rss_gib)

    training = reconstruct_exp064_training(args, exp064_result)
    selection = select_gate(training["features"], training["c_minus_d"])
    permutation = full_selection_permutation_null(
        training["features"], training["c_minus_d"]
    )
    training_gate = training_promotion_gate(selection, permutation)
    serialized_selection = serialize_selection(selection)
    preamble: dict[str, object] = {
        "schema_version": 1,
        "experiment": EXPERIMENT,
        "protocol": protocol,
        "artifact_manifest_sha256": sha256_file(args.artifact_manifest),
        "training": {
            "truth_source": "experiment-064 HCT116 and HEK293T only",
            "targets": EXP064_TARGETS,
            "target_order_sha256": training["target_order_sha256"],
            "folds": training["folds"],
            "selection": serialized_selection,
            "full_selection_permutation_null": permutation,
            "promotion_gate": training_gate,
        },
        "pre_kolf_input_hashes": pre_hashes,
        "frozen_runtime_packages": frozen_runtime,
        "kolf_opened": False,
        "count_emission": False,
        "cell_eval2": False,
        "challenge_submission": False,
    }
    if not training_gate["passed"]:
        result = {
            **preamble,
            "status": "training_gate_failed",
            "decision": "do_not_open_kolf",
            "reason": "predeclared HCT116/HEK293T gate did not pass",
            "runtime": runtime_metadata(),
            "peak_rss_gib": peak_rss_gib(),
        }
        atomic_write_text(
            output / "results.json",
            json.dumps(result, indent=2, sort_keys=True) + "\n",
        )
        print(json.dumps(result["training"]["promotion_gate"], indent=2))
        return result

    official = read_csv_column(args.genes_csv, "gene_name")
    panel = read_csv_column(args.panel_csv, "target_gene")
    if len(official) != OFFICIAL_GENES or len(panel) != PANEL_TARGETS:
        raise ValueError("challenge axes changed after authentication")
    if sha256_labels(panel) != EXPECTED_PANEL_ORDER_SHA256:
        raise ValueError("challenge panel order hash changed")
    selected_rule = selection["selected"]["rule"]
    production_holder: dict[str, object] = {}

    def freeze_builder() -> Mapping[str, object]:
        production_holder.update(
            build_production_freeze_inputs(
                args,
                panel,
                official,
                selected_rule,
            )
        )
        return freeze_production_artifacts(
            output,
            panel,
            selected_rule,
            production_holder,
        )

    def truth_loader(frozen: Mapping[str, object]) -> HeldTruth:
        durable_freeze = persist_pre_kolf_freeze(
            frozen,
            os.environ["EXP085_DURABLE_FREEZE_URI"],
        )
        access = mark_kolf_access(
            output,
            frozen,
            durable_freeze,
            os.environ["EXP085_DURABLE_MARKER_URI"],
        )
        production_holder["kolf_access_marker"] = access
        return load_kolf_truth(
            args,
            manifest,
            official,
            panel,
            max_rss_gib=args.max_rss_gib,
        )

    frozen, held = execute_after_freeze(
        bool(training_gate["passed"]),
        freeze_builder,
        truth_loader,
    )
    assert frozen is not None and held is not None
    evaluation = evaluate_kolf_once(output, held, production_holder)
    result = {
        **preamble,
        "status": "completed",
        "decision": (
            "promote_effect_level_gate"
            if evaluation["promotion_gate"]["passed"]
            else "reject_gate"
        ),
        "kolf_opened": True,
        "candidate_count_evaluated_on_kolf": 1,
        "pre_kolf_freeze": frozen,
        "kolf_access_marker": production_holder["kolf_access_marker"],
        "held_truth": held.audit,
        "evaluation": evaluation,
        "leakage_audit": {
            "exp064_C_D_reproduced_before_training": True,
            "training_labels_only_HCT116_and_HEK293T": True,
            "feature_source_only": True,
            "production_built_for_all_300_before_kolf_open": True,
            "arm_C_D_assignments_and_clusters_serialized_before_kolf_open": True,
            "kolf_rows_in_training": 0,
            "kolf_rows_in_whitening_operator": 0,
            "kolf_rows_in_source_clusters": 0,
            "kolf_control_used": False,
            "post_kolf_fit_or_selection": False,
            "complete_273_way_gallery": True,
            "all_300_panel_coordinates_excluded_from_operator_cd4_and_score": True,
        },
        "runtime": runtime_metadata(),
        "peak_rss_gib": peak_rss_gib(),
    }
    atomic_write_text(
        output / "results.json",
        json.dumps(result, indent=2, sort_keys=True) + "\n",
    )
    print(json.dumps(evaluation["scores"], indent=2), flush=True)
    print(json.dumps(evaluation["promotion_gate"], indent=2), flush=True)
    return result


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--k562", type=Path, required=True)
    parser.add_argument("--hct116", type=Path, required=True)
    parser.add_argument("--hek293t", type=Path, required=True)
    parser.add_argument("--perturbmatch", type=Path, required=True)
    parser.add_argument("--cd4-root", type=Path, required=True)
    parser.add_argument("--genes-csv", type=Path, required=True)
    parser.add_argument("--panel-csv", type=Path, required=True)
    parser.add_argument("--held-whitener-dir", type=Path, required=True)
    parser.add_argument("--production-whitener", type=Path, required=True)
    parser.add_argument("--production-whitener-metadata", type=Path, required=True)
    parser.add_argument("--exp064-results", type=Path, required=True)
    parser.add_argument("--exp064-per-target", type=Path, required=True)
    parser.add_argument("--kolf", type=Path, required=True)
    parser.add_argument(
        "--artifact-manifest",
        type=Path,
        default=HERE / "artifact_manifest.json",
    )
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--chunk-rows", type=int, default=128)
    parser.add_argument("--max-rss-gib", type=float, default=40.0)
    parser.add_argument("--guarded", action="store_true")
    return parser


def main() -> int:
    args = build_parser().parse_args()
    if not args.guarded:
        raise SystemExit("refusing real-data execution without --guarded")
    if os.environ.get("EXP085_ONE_TIME_CLEARANCE") != "YES":
        raise SystemExit(
            "refusing one-time KOLF path without EXP085_ONE_TIME_CLEARANCE=YES"
        )
    if not os.environ.get("EXP085_DURABLE_MARKER_URI"):
        raise SystemExit(
            "refusing one-time KOLF path without EXP085_DURABLE_MARKER_URI"
        )
    if not os.environ.get("EXP085_DURABLE_FREEZE_URI"):
        raise SystemExit(
            "refusing one-time KOLF path without EXP085_DURABLE_FREEZE_URI"
        )
    if args.chunk_rows < 1:
        raise SystemExit("--chunk-rows must be positive")
    if not (8.0 <= args.max_rss_gib <= 44.0):
        raise SystemExit("--max-rss-gib must be in [8, 44]")
    run_experiment(args)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
