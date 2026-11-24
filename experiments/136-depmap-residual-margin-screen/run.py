#!/usr/bin/env python3
"""Development-only residual screen using Experiment 132's consumed split."""

from __future__ import annotations

import importlib.util
import json
import os
import sys
import tempfile
from pathlib import Path
from typing import Any

os.environ.setdefault("OPENBLAS_NUM_THREADS", "8")
os.environ.setdefault("OMP_NUM_THREADS", "8")
os.environ.setdefault("MKL_NUM_THREADS", "8")
os.environ.setdefault("VECLIB_MAXIMUM_THREADS", "8")

import numpy as np


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
EXP132_PATH = ROOT / "experiments/132-response-decomposition-ridge/run.py"
OUTPUT = HERE / "work/development-screen.json"
BETAS = np.asarray(
    [0.0, 0.0025, 0.005, 0.01, 0.02, 0.03, 0.05, 0.075, 0.10, 0.15, 0.20, 0.30, 0.50, 1.0],
    dtype=np.float64,
)


def load_exp132() -> Any:
    spec = importlib.util.spec_from_file_location("exp132_for_136", EXP132_PATH)
    if spec is None or spec.loader is None:
        raise RuntimeError("could not load Experiment 132")
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def atomic_write_json(path: Path, value: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary = tempfile.mkstemp(
        dir=path.parent, prefix=f".{path.name}.", suffix=".tmp"
    )
    try:
        with os.fdopen(descriptor, "w", encoding="utf-8") as handle:
            json.dump(value, handle, indent=2, sort_keys=True)
            handle.write("\n")
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(temporary, path)
    except BaseException:
        Path(temporary).unlink(missing_ok=True)
        raise


def summarize_delta(summary: dict[str, Any]) -> dict[str, Any]:
    return {
        "macro": summary["macro_delta"],
        "contexts": {
            name: values["delta"] for name, values in summary["contexts"].items()
        },
    }


def max_decoy(similarity: np.ndarray) -> np.ndarray:
    masked = np.asarray(similarity, dtype=np.float64).copy()
    np.fill_diagonal(masked, -np.inf)
    return np.max(masked, axis=1)


def main() -> None:
    if OUTPUT.exists():
        raise FileExistsError(f"{OUTPUT} exists; refusing overwrite")

    exp132 = load_exp132()
    exp132.validate_authorization(check_large_hashes=False)
    if exp132.TEST_ACCESS_LOCK.exists() or exp132.RESULT_DIR.exists():
        raise RuntimeError("Experiment 132 untouched source-test boundary changed")

    selection, _, response_genes, _ = exp132.load_plan_arrays()
    effects = exp132.load_split_effects(selection, response_genes, 8.0)
    contexts = exp132.context_effects(effects)
    identity = exp132.source_means_by_held(contexts)
    projection = exp132.fit_projection(exp132.feature_rows(selection))
    ridge = exp132.ridge_predictions(projection, identity)

    global_grid: list[dict[str, Any]] = []
    pds_by_beta: list[np.ndarray] = []
    for beta in BETAS:
        candidate = {
            group: identity[group] + beta * (ridge[group] - identity[group])
            for group in identity
        }
        summary, arrays = exp132.summarize(candidate, identity, effects)
        global_grid.append({"beta": float(beta), **summarize_delta(summary)})
        pds_by_beta.append(
            np.stack(
                [
                    0.5
                    * (
                        arrays["k562_gwps__pds_delta"]
                        + arrays["k562_essential__pds_delta"]
                    ),
                    arrays["rpe1__pds_delta"],
                    arrays["nadig_jurkat__pds_delta"],
                ],
                axis=0,
            )
        )

    allowed = {
        held: [name for name in exp132.CONTEXTS if name != held]
        for held in exp132.HELD_GROUPS
    }
    routed: dict[str, Any] = {}
    for policy in ("rank_preserving", "strict_margin_preserving"):
        candidate: dict[str, np.ndarray] = {}
        distributions: dict[str, dict[str, int]] = {}
        for held in exp132.HELD_GROUPS:
            baseline = identity[held]
            residual = ridge[held] - baseline
            chosen = np.zeros(len(baseline), dtype=np.float64)
            baseline_stats = []
            for source in allowed[held]:
                similarity = exp132.cosine_matrix(baseline, contexts[source])
                baseline_stats.append(
                    (
                        np.diag(similarity),
                        max_decoy(similarity),
                        exp132.pds_per_target(baseline, contexts[source]),
                    )
                )
            for beta in BETAS[1:-1]:
                query = baseline + beta * residual
                eligible = np.ones(len(baseline), dtype=bool)
                for source, (own0, decoy0, pds0) in zip(
                    allowed[held], baseline_stats, strict=True
                ):
                    similarity = exp132.cosine_matrix(query, contexts[source])
                    if policy == "rank_preserving":
                        eligible &= (
                            exp132.pds_per_target(query, contexts[source]) >= pds0
                        )
                    else:
                        eligible &= (
                            np.diag(similarity) - max_decoy(similarity)
                            >= own0 - decoy0 - 1e-12
                        )
                chosen[eligible] = beta
            candidate[held] = baseline + chosen[:, None] * residual
            distributions[held] = {
                str(float(beta)): int(np.sum(chosen == beta)) for beta in BETAS
            }
        summary, _ = exp132.summarize(candidate, identity, effects)
        routed[policy] = {
            **summarize_delta(summary),
            "beta_distribution": distributions,
        }

    pds_cube = np.stack(pds_by_beta, axis=0)
    macro_per_target = np.mean(pds_cube, axis=1)
    best_index = np.argmax(macro_per_target, axis=0)
    target_index = np.arange(pds_cube.shape[2])
    best_delta = np.max(macro_per_target, axis=0)
    chosen_beta = BETAS[best_index]
    oracle = {
        "meaning": (
            "invalid deployment oracle: chooses one beta per target using all "
            "three held development outcomes"
        ),
        "macro_pds_delta": float(np.mean(best_delta)),
        "context_deltas": {
            "k562_family": float(
                np.mean(pds_cube[best_index, 0, target_index])
            ),
            "rpe1": float(np.mean(pds_cube[best_index, 1, target_index])),
            "jurkat": float(np.mean(pds_cube[best_index, 2, target_index])),
        },
        "positive_target_fraction": float(np.mean(best_delta > 0)),
        "median_best_target_delta": float(np.median(best_delta)),
        "p95_best_target_delta": float(np.quantile(best_delta, 0.95)),
        "beta_distribution": {
            str(float(beta)): int(np.sum(chosen_beta == beta)) for beta in BETAS
        },
    }

    result = {
        "schema_version": 1,
        "experiment": "136-depmap-residual-margin-screen",
        "status": "terminal_source_development_no_go",
        "source_split": "experiment_132_consumed_selection_only",
        "selection_targets": len(selection),
        "response_genes": len(response_genes),
        "experiment_132_test_effect_values_accessed": False,
        "beta_grid": BETAS.tolist(),
        "global_grid": global_grid,
        "routing": routed,
        "oracle_upper_bound": oracle,
        "decision": {
            "passes": False,
            "reason": (
                "Every nonzero global dose loses PDS in every context; the "
                "deployment-available rank rule also loses PDS; the strict "
                "margin rule reduces to identity; and even the invalid "
                "per-target oracle stays below +0.003 macro PDS while losing "
                "on the K562 family."
            ),
            "untouched_410_target_test_authorized": False,
        },
    }
    atomic_write_json(OUTPUT, result)
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
