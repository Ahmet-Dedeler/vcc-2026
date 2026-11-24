#!/usr/bin/env python3
"""Score Experiment 152 only after a future one-time scoring authorization."""

from __future__ import annotations

import hashlib
import importlib.util
import json
import os
import sys
from pathlib import Path
from types import ModuleType
from typing import Any

import numpy as np
import pandas as pd


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
PREREGISTRATION_PATH = HERE / "PREREGISTRATION.v1.json"
PREREGISTRATION_DIGEST_PATH = HERE / "PREREGISTRATION.v1.sha256"
BUILD_AUTHORIZATION_PATH = HERE / "SOURCE_BUILD_AUTHORIZATION.v1.json"
SCORING_AUTHORIZATION_PATH = HERE / "SCORING_AUTHORIZATION.v1.json"
BUILD_DIR = HERE / "work/source-build-v1"
CANDIDATE_PATH = BUILD_DIR / "candidate.h5ad"
BUILD_REPORT_PATH = BUILD_DIR / "build-report.json"
GENE_RAMP_PATH = BUILD_DIR / "gene-ramp.json"
REAL_PATH = ROOT / "data/harness/orion-hct116-thinned.h5ad"
EXPERIMENT_150_SCORING = (
    ROOT / "experiments/150-dual-moment-score-level-hct116/work/scoring-v1"
)
INCUMBENT_RESULTS = EXPERIMENT_150_SCORING / "incumbent/run/results.csv"
UNIFORM_RESULTS = EXPERIMENT_150_SCORING / "candidate/run/results.csv"
EXPERIMENT_150_RESULTS = EXPERIMENT_150_SCORING / "results.json"
EXPERIMENT_151_SCORING = (
    ROOT / "experiments/151-relaxation-consensus-selector/work/scoring-v1"
)
TARGETWISE_RESULTS = EXPERIMENT_151_SCORING / "candidate/run/results.csv"
EXPERIMENT_151_RESULTS = EXPERIMENT_151_SCORING / "results.json"
OUTDIR = HERE / "work/scoring-v1"
PERT_COL = "target_gene"
CONTROL = "non-targeting"
REQUIRED_SCORING_AUTHORIZATION_BINDINGS = frozenset(
    {
        "experiments/152-gene-consensus-ramp/PREREGISTRATION.v1.json",
        "experiments/152-gene-consensus-ramp/PREREGISTRATION.v1.sha256",
        "experiments/152-gene-consensus-ramp/SOURCE_BUILD_AUTHORIZATION.v1.json",
        "experiments/152-gene-consensus-ramp/build_candidate.py",
        "experiments/152-gene-consensus-ramp/test_experiment152.py",
        "experiments/152-gene-consensus-ramp/work/source-build-v1/build-report.json",
        "experiments/152-gene-consensus-ramp/work/source-build-v1/gene-ramp.json",
        "experiments/152-gene-consensus-ramp/work/source-build-v1/candidate.h5ad",
        "experiments/150-dual-moment-score-level-hct116/work/scoring-v1/results.json",
        "experiments/150-dual-moment-score-level-hct116/work/scoring-v1/incumbent/run/results.csv",
        "experiments/150-dual-moment-score-level-hct116/work/scoring-v1/candidate/run/results.csv",
        "experiments/151-relaxation-consensus-selector/work/scoring-v1/results.json",
        "experiments/151-relaxation-consensus-selector/work/scoring-v1/candidate/run/results.csv",
        "experiments/149-dual-moment-hek293t/score_confirmation.py",
        "scripts/run_with_memory_guard.py",
        "data/harness/orion-hct116-thinned.h5ad",
        "data/harness/orion-hct116-thinned.h5ad.meta.json",
    }
)


def load_module(name: str, path: Path) -> ModuleType:
    spec = importlib.util.spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise ImportError(path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(8 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def sha256_strings(values: tuple[str, ...] | list[str] | np.ndarray) -> str:
    payload = "\n".join(map(str, values)) + "\n"
    return hashlib.sha256(payload.encode()).hexdigest()


def atomic_write_json(path: Path, payload: dict[str, Any]) -> None:
    temporary = path.with_name(f".{path.name}.tmp-{os.getpid()}")
    temporary.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    os.replace(temporary, path)


def load_preregistration() -> dict[str, Any]:
    digest_line = PREREGISTRATION_DIGEST_PATH.read_text().strip().split()
    if len(digest_line) != 2 or digest_line[1] != PREREGISTRATION_PATH.name:
        raise ValueError("malformed immutable preregistration digest file")
    observed = sha256(PREREGISTRATION_PATH)
    if observed != digest_line[0]:
        raise ValueError("immutable preregistration hash changed")
    protocol = json.loads(PREREGISTRATION_PATH.read_text())
    if protocol.get("experiment") != "152-gene-consensus-ramp":
        raise ValueError("wrong preregistration experiment")
    return protocol


def _bound_path(relative: str) -> Path:
    path = (ROOT / relative).resolve()
    try:
        path.relative_to(ROOT.resolve())
    except ValueError as error:
        raise ValueError(f"authorization binding escapes repository: {relative}") from error
    return path


def verify_scoring_authorization(
    path: Path = SCORING_AUTHORIZATION_PATH,
) -> dict[str, Any]:
    if not path.is_file():
        raise PermissionError(
            "Experiment 152 scoring requires a future one-time "
            "SCORING_AUTHORIZATION.v1.json"
        )
    authorization = json.loads(path.read_text())
    if authorization.get("experiment") != "152-gene-consensus-ramp":
        raise ValueError("wrong scoring authorization experiment")
    if authorization.get("status") != "authorized_once_scoring":
        raise PermissionError("Experiment 152 scoring is not authorized once")
    scorer = authorization.get("scorer")
    if not isinstance(scorer, dict):
        raise ValueError("scoring authorization has no scorer binding")
    if scorer.get("sha256") != sha256(Path(__file__).resolve()):
        raise ValueError("scorer source does not match scoring authorization")
    bindings = authorization.get("bindings_sha256")
    if not isinstance(bindings, dict) or not bindings:
        raise ValueError("scoring authorization has no bindings")
    missing = sorted(REQUIRED_SCORING_AUTHORIZATION_BINDINGS - set(bindings))
    if missing:
        raise ValueError(f"scoring authorization misses bindings: {missing}")
    for relative, expected in bindings.items():
        observed = sha256(_bound_path(relative))
        if observed != expected:
            raise ValueError(f"scoring authorization binding changed: {relative}")
    return authorization


def _load_score_runtime() -> ModuleType:
    base = load_module(
        "exp152_score_base",
        ROOT / "experiments/149-dual-moment-hek293t/score_confirmation.py",
    )
    if float(base.MAX_RSS_GIB) != 24.0:
        raise ValueError("score guard is not frozen at 24 GiB")
    return base


def _metric_values(
    frame: pd.DataFrame,
    targets: tuple[str, ...],
    metric: str,
) -> np.ndarray:
    rows = frame[frame["metric"] == metric].set_index("perturbation")
    try:
        return rows.loc[list(targets), "value"].to_numpy(dtype=np.float64)
    except KeyError as error:
        raise ValueError(f"{metric}: result frame misses a primary target") from error


def projected_frame_score(
    base: ModuleType,
    frame: pd.DataFrame,
    targets: tuple[str, ...],
) -> dict[str, float]:
    arrays = base.metric_arrays(frame, np.asarray(targets))
    indices = np.arange(len(targets), dtype=np.int64)
    return base.official(base.aggregate_indices(arrays, indices))


def choose_best_frozen_comparator(
    projected_average_by_name: dict[str, float],
) -> str:
    """Choose the higher projected-average comparator with a frozen tie break."""
    required = {"uniform_lambda_1", "targetwise_selector"}
    if set(projected_average_by_name) != required:
        raise ValueError("comparator projections must contain exactly two frozen arms")
    values = projected_average_by_name
    if not all(np.isfinite(value) for value in values.values()):
        raise ValueError("comparator projected averages must be finite")
    if values["targetwise_selector"] > values["uniform_lambda_1"]:
        return "targetwise_selector"
    return "uniform_lambda_1"


def compare_frames(
    base: ModuleType,
    left: pd.DataFrame,
    right: pd.DataFrame,
    targets: tuple[str, ...],
    *,
    draws: int,
    seed: int,
) -> dict[str, Any]:
    left_arrays = base.metric_arrays(left, np.asarray(targets))
    right_arrays = base.metric_arrays(right, np.asarray(targets))
    indices = np.arange(len(targets), dtype=np.int64)
    left_raw = base.aggregate_indices(left_arrays, indices)
    right_raw = base.aggregate_indices(right_arrays, indices)
    left_scaled = base.official(left_raw)
    right_scaled = base.official(right_raw)
    left_pds = _metric_values(left, targets, "pds_cosine")
    right_pds = _metric_values(right, targets, "pds_cosine")
    return {
        "left_raw": left_raw,
        "right_raw": right_raw,
        "raw_delta": base.delta(left_raw, right_raw),
        "left_official_projection": left_scaled,
        "right_official_projection": right_scaled,
        "scaled_delta": base.delta(left_scaled, right_scaled),
        "all_per_target_pds_exact": bool(np.array_equal(left_pds, right_pds)),
        "paired_target_bootstrap": base.paired_bootstrap(
            left,
            right,
            np.asarray(targets),
            draws=draws,
            seed=seed,
        ),
    }


def evaluate_primary_gates(
    protocol: dict[str, Any],
    versus_incumbent: dict[str, Any],
    versus_best_comparator: dict[str, Any],
) -> dict[str, Any]:
    de_metrics = (
        "de_wilcoxon_lfc_nmae",
        "de_wilcoxon_direction_reach_raw",
        "de_wilcoxon_direction_fidelity_yield_raw",
        "de_wilcoxon_sig_jaccard",
    )
    improvement_metrics = (
        "de_wilcoxon_lfc_nmae",
        "de_wilcoxon_direction_reach_raw",
        "de_wilcoxon_direction_fidelity_yield_raw",
    )
    a_rules = protocol["primary_analysis"]["comparison_a_gates"]
    b_rules = protocol["primary_analysis"]["comparison_b_gates"]
    a_delta = versus_incumbent["scaled_delta"]
    b_delta = versus_best_comparator["scaled_delta"]
    improved_count = sum(a_delta[metric] > 0.0 for metric in improvement_metrics)
    minimum_de = min(a_delta[metric] for metric in de_metrics)

    comparison_a = {
        "all_per_target_pds_exact": versus_incumbent[
            "all_per_target_pds_exact"
        ],
        "projected_avg_gain": (
            a_delta["avg_score"]
            >= float(a_rules["official_anchor_projected_avg_delta_min"])
        ),
        "bootstrap_lower_bound": (
            versus_incumbent["paired_target_bootstrap"]["lower_2_5pct"]
            > float(a_rules["target_bootstrap_95pct_lower_min"])
        ),
        "minimum_improved_metrics_among_nmae_reach_fidelity": (
            improved_count
            >= int(
                a_rules["minimum_improved_metrics_among_nmae_reach_fidelity"]
            )
        ),
        "minimum_individual_de_scaled_delta": (
            minimum_de >= float(a_rules["minimum_individual_de_scaled_delta"])
        ),
        "minimum_scaled_mse_delta": (
            a_delta["expr_mse_unbiased_capped_norm"]
            >= float(a_rules["minimum_scaled_mse_delta"])
        ),
    }
    comparison_b = {
        "projected_avg_gain": (
            b_delta["avg_score"]
            >= float(b_rules["official_anchor_projected_avg_delta_min"])
        ),
        "bootstrap_lower_bound": (
            versus_best_comparator["paired_target_bootstrap"]["lower_2_5pct"]
            >= float(b_rules["target_bootstrap_95pct_lower_min"])
        ),
    }
    a_passes = all(comparison_a.values())
    b_passes = all(comparison_b.values())
    return {
        "comparison_a": {
            "checks": comparison_a,
            "passes": a_passes,
            "improved_metric_count_among_nmae_reach_fidelity": improved_count,
            "minimum_de_scaled_delta": minimum_de,
        },
        "comparison_b": {
            "checks": comparison_b,
            "passes": b_passes,
        },
        "passes": a_passes and b_passes,
    }


def main() -> int:
    protocol = load_preregistration()
    authorization = verify_scoring_authorization()
    if OUTDIR.exists():
        raise FileExistsError(f"immutable scoring output already exists: {OUTDIR}")
    gene_ramp = json.loads(GENE_RAMP_PATH.read_text())
    primary_targets = tuple(gene_ramp["primary_targets_hct_order"])
    expected_primary_hash = protocol["common_surface"][
        "primary_targets_in_hct_order_sha256_newline_terminated"
    ]
    if sha256_strings(primary_targets) != expected_primary_hash:
        raise ValueError("built primary-target hash differs from preregistration")
    if len(primary_targets) != int(protocol["common_surface"]["common_target_count"]):
        raise ValueError("built primary-target count differs")

    for path in (
        CANDIDATE_PATH,
        BUILD_REPORT_PATH,
        INCUMBENT_RESULTS,
        UNIFORM_RESULTS,
        EXPERIMENT_150_RESULTS,
        TARGETWISE_RESULTS,
        EXPERIMENT_151_RESULTS,
    ):
        if not path.is_file():
            raise FileNotFoundError(path)

    base = _load_score_runtime()
    OUTDIR.mkdir(parents=True)
    candidate_dir = OUTDIR / "candidate"
    candidate_dir.mkdir()
    memory_guard = base.run_guarded(
        [
            "uv",
            "run",
            "cell-eval2",
            "run",
            "-ap",
            str(CANDIDATE_PATH),
            "-ar",
            str(REAL_PATH),
            "-o",
            str(candidate_dir / "run"),
            "--preset",
            "vcc2026",
            "--pert-col",
            PERT_COL,
            "--control",
            CONTROL,
        ],
        candidate_dir,
    )

    incumbent_frame = pd.read_csv(INCUMBENT_RESULTS)
    uniform_frame = pd.read_csv(UNIFORM_RESULTS)
    targetwise_frame = pd.read_csv(TARGETWISE_RESULTS)
    candidate_frame = pd.read_csv(candidate_dir / "run/results.csv")
    projections = {
        "uniform_lambda_1": projected_frame_score(
            base,
            uniform_frame,
            primary_targets,
        )["avg_score"],
        "targetwise_selector": projected_frame_score(
            base,
            targetwise_frame,
            primary_targets,
        )["avg_score"],
    }
    best_name = choose_best_frozen_comparator(projections)
    best_frame = (
        targetwise_frame if best_name == "targetwise_selector" else uniform_frame
    )
    analysis = protocol["primary_analysis"]
    versus_incumbent = compare_frames(
        base,
        incumbent_frame,
        candidate_frame,
        primary_targets,
        draws=int(analysis["bootstrap_draws"]),
        seed=int(analysis["bootstrap_seed_vs_incumbent"]),
    )
    versus_best = compare_frames(
        base,
        best_frame,
        candidate_frame,
        primary_targets,
        draws=int(analysis["bootstrap_draws"]),
        seed=int(analysis["bootstrap_seed_vs_best_frozen_comparator"]),
    )
    gates = evaluate_primary_gates(protocol, versus_incumbent, versus_best)
    result = {
        "experiment": "152-gene-consensus-ramp",
        "stage": "authorized_held_hct_candidate_score",
        "scoring_authorization_sha256": sha256(SCORING_AUTHORIZATION_PATH),
        "source_build_authorization_sha256": sha256(BUILD_AUTHORIZATION_PATH),
        "preregistration_sha256": sha256(PREREGISTRATION_PATH),
        "scorer_sha256": sha256(Path(__file__).resolve()),
        "candidate_sha256": sha256(CANDIDATE_PATH),
        "candidate_results_sha256": sha256(
            candidate_dir / "run/results.csv"
        ),
        "experiment_150_results_sha256": sha256(EXPERIMENT_150_RESULTS),
        "experiment_151_results_sha256": sha256(EXPERIMENT_151_RESULTS),
        "primary_target_count": len(primary_targets),
        "primary_target_hash": sha256_strings(primary_targets),
        "best_frozen_comparator": {
            "selected": best_name,
            "projected_average_by_name": projections,
            "tie_break": "uniform_lambda_1",
        },
        "primary": {
            "versus_incumbent": versus_incumbent,
            "versus_best_frozen_comparator": versus_best,
            "gates": gates,
        },
        "diagnostics_only": {
            "requested": protocol["diagnostics_only"],
            "used_for_gating": False,
        },
        "passes": gates["passes"],
        "decision": (
            "ADVANCE_GENE_RAMP_TO_PRODUCTION_BUILD"
            if gates["passes"]
            else "CLOSE_EXACT_GENE_WISE_CONSENSUS_RAMP_FAMILY"
        ),
        "memory_guard": memory_guard,
        "only_new_candidate_scored": True,
        "parallel_scorers_used": False,
        "challenge_build_created": False,
        "challenge_submission_created": False,
        "authorization_status": authorization["status"],
    }
    atomic_write_json(OUTDIR / "results.json", result)
    print(
        json.dumps(
            {
                "decision": result["decision"],
                "passes": result["passes"],
                "best_frozen_comparator": best_name,
                "comparison_a": gates["comparison_a"],
                "comparison_b": gates["comparison_b"],
            },
            indent=2,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
