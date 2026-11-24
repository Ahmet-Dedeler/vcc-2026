#!/usr/bin/env python3
"""Score the frozen held-HCT score-level arbitration serially."""

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
AUTHORIZATION_PATH = HERE / "DEVELOPMENT_ACCESS_AUTHORIZATION.v1.json"
BUILD_DIR = HERE / "work/source-build-v1"
EFFECTS_PATH = BUILD_DIR / "effects.npz"
REAL_PATH = ROOT / "data/harness/orion-hct116-thinned.h5ad"
PREDICTIONS = {
    "incumbent": BUILD_DIR / "incumbent.h5ad",
    "candidate": BUILD_DIR / "candidate.h5ad",
}
OUTDIR = HERE / "work/scoring-v1"


def load_module(name: str, path: Path) -> ModuleType:
    spec = importlib.util.spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise ImportError(path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


BASE = load_module(
    "exp150_score_base",
    ROOT / "experiments/149-dual-moment-hek293t/score_confirmation.py",
)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(8 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def sha256_strings(values: np.ndarray) -> str:
    payload = "\n".join(map(str, values)) + "\n"
    return hashlib.sha256(payload.encode()).hexdigest()


def atomic_write_json(path: Path, payload: dict[str, Any]) -> None:
    temporary = path.with_name(f".{path.name}.tmp-{os.getpid()}")
    temporary.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    os.replace(temporary, path)


def verify_authorization() -> dict[str, Any]:
    authorization = json.loads(AUTHORIZATION_PATH.read_text())
    if authorization.get("status") != "authorized_once":
        raise ValueError("held-context scoring is not authorized")
    if authorization["scorer"]["sha256"] != sha256(Path(__file__).resolve()):
        raise ValueError("scorer source does not match authorization")
    bindings = authorization["bindings_sha256"]
    observed = {
        str((ROOT / relative).relative_to(ROOT)): sha256(ROOT / relative)
        for relative in bindings
    }
    if observed != bindings:
        raise ValueError("frozen scoring bindings do not match authorization")
    return authorization


def main() -> int:
    if OUTDIR.exists():
        raise FileExistsError(f"immutable scoring output exists: {OUTDIR}")
    authorization = verify_authorization()
    OUTDIR.mkdir(parents=True)
    memory: dict[str, dict[str, Any]] = {}
    for label, prediction in PREDICTIONS.items():
        arm_dir = OUTDIR / label
        arm_dir.mkdir()
        print(f"[score] {label}", flush=True)
        memory[label] = BASE.run_guarded(
            [
                "uv",
                "run",
                "cell-eval2",
                "run",
                "-ap",
                str(prediction),
                "-ar",
                str(REAL_PATH),
                "-o",
                str(arm_dir / "run"),
                "--preset",
                "vcc2026",
                "--pert-col",
                BASE.PERT_COL,
                "--control",
                BASE.CONTROL,
            ],
            arm_dir,
        )

    incumbent_agg = BASE.aggregate(
        OUTDIR / "incumbent/run/agg_results.csv"
    )
    candidate_agg = BASE.aggregate(
        OUTDIR / "candidate/run/agg_results.csv"
    )
    raw_delta = BASE.delta(incumbent_agg, candidate_agg)
    incumbent_scaled = BASE.official(incumbent_agg)
    candidate_scaled = BASE.official(candidate_agg)
    scaled_delta = BASE.delta(incumbent_scaled, candidate_scaled)

    incumbent_frame = pd.read_csv(
        OUTDIR / "incumbent/run/results.csv"
    )
    candidate_frame = pd.read_csv(
        OUTDIR / "candidate/run/results.csv"
    )
    with np.load(EFFECTS_PATH, allow_pickle=False) as archive:
        targets = np.asarray(archive["targets"]).astype(str)
    if len(targets) != int(authorization["target_count"]):
        raise ValueError("authorized target count differs")
    if sha256_strings(targets) != authorization["target_order_sha256"]:
        raise ValueError("authorized target order differs")
    incumbent_pds = (
        incumbent_frame[incumbent_frame["metric"] == "pds_cosine"]
        .set_index("perturbation")
        .loc[targets, "value"]
        .to_numpy(dtype=np.float64)
    )
    candidate_pds = (
        candidate_frame[candidate_frame["metric"] == "pds_cosine"]
        .set_index("perturbation")
        .loc[targets, "value"]
        .to_numpy(dtype=np.float64)
    )
    pds_exact = bool(np.array_equal(incumbent_pds, candidate_pds))
    bootstrap = BASE.paired_bootstrap(
        incumbent_frame,
        candidate_frame,
        targets,
        draws=int(authorization["bootstrap"]["draws"]),
        seed=int(authorization["bootstrap"]["seed"]),
    )

    de_metrics = (
        "de_wilcoxon_lfc_nmae",
        "de_wilcoxon_direction_reach_raw",
        "de_wilcoxon_direction_fidelity_yield_raw",
        "de_wilcoxon_sig_jaccard",
    )
    nonnegative_de = sum(scaled_delta[metric] >= 0.0 for metric in de_metrics)
    minimum_de = min(scaled_delta[metric] for metric in de_metrics)
    checks = {
        "all_per_target_pds_exact": pds_exact,
        "projected_avg_gain_at_least_0.004": (
            scaled_delta["avg_score"] >= 0.004
        ),
        "bootstrap_lower_bound_above_zero": (
            bootstrap["lower_2_5pct"] > 0.0
        ),
        "at_least_two_nonnegative_de_metrics": nonnegative_de >= 2,
        "no_de_metric_below_minus_0.015": minimum_de >= -0.015,
        "scaled_mse_does_not_decrease": (
            scaled_delta["expr_mse_unbiased_capped_norm"] >= 0.0
        ),
    }
    passed = all(checks.values())
    result = {
        "experiment": "150-dual-moment-score-level-hct116",
        "stage": "held_context_score_level_arbitration",
        "provenance_utc": "2026-08-31T18:45:00+00:00",
        "authorization_sha256": sha256(AUTHORIZATION_PATH),
        "incumbent_results_sha256": sha256(
            OUTDIR / "incumbent/run/results.csv"
        ),
        "candidate_results_sha256": sha256(
            OUTDIR / "candidate/run/results.csv"
        ),
        "incumbent_aggregate_sha256": sha256(
            OUTDIR / "incumbent/run/agg_results.csv"
        ),
        "candidate_aggregate_sha256": sha256(
            OUTDIR / "candidate/run/agg_results.csv"
        ),
        "target_count": len(targets),
        "incumbent_raw": incumbent_agg,
        "candidate_raw": candidate_agg,
        "raw_delta": raw_delta,
        "incumbent_official_projection": incumbent_scaled,
        "candidate_official_projection": candidate_scaled,
        "scaled_delta": scaled_delta,
        "nonnegative_de_metric_count": nonnegative_de,
        "minimum_de_scaled_delta": minimum_de,
        "paired_target_bootstrap": bootstrap,
        "checks": checks,
        "passes": passed,
        "decision": (
            "ADVANCE_TO_PRODUCTION_BUILD" if passed else "TERMINAL_NO_GO"
        ),
        "memory_guard": memory,
        "parallel_scorers_used": False,
        "challenge_build_created": False,
        "challenge_submission_created": False,
    }
    atomic_write_json(OUTDIR / "results.json", result)
    print(
        json.dumps(
            {
                "decision": result["decision"],
                "raw_delta": raw_delta,
                "scaled_delta": scaled_delta,
                "bootstrap": bootstrap,
                "nonnegative_de_metric_count": nonnegative_de,
                "minimum_de_scaled_delta": minimum_de,
                "checks": checks,
            },
            indent=2,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
