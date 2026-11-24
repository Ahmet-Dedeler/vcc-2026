#!/usr/bin/env python3
"""Fit/freeze experiment 083 without opening any held VCC-panel truth."""

from __future__ import annotations

import argparse
import csv
import dataclasses
import hashlib
import json
import os
import resource
import sys
from pathlib import Path
from typing import Iterable, Mapping, Sequence

os.environ.setdefault("OPENBLAS_NUM_THREADS", "8")
os.environ.setdefault("OMP_NUM_THREADS", "8")
os.environ.setdefault("MKL_NUM_THREADS", "8")
os.environ.setdefault("VECLIB_MAXIMUM_THREADS", "8")
os.environ["PYTHONDONTWRITEBYTECODE"] = "1"

import numpy as np

import model


NULL_REPETITIONS_PER_FAMILY = 100
MIN_PANEL_MULTIGUIDE = 100
TEST_BOOTSTRAP_REPETITIONS = 10_000
SPLIT_SEED = model.SPLIT_SEED


def peak_rss_gib() -> float:
    value = float(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
    return value / 1024**3 if sys.platform == "darwin" else value / 1024**2


def enforce_memory(limit_gib: float, stage: str) -> None:
    peak = peak_rss_gib()
    if limit_gib > 0 and peak >= limit_gib:
        raise MemoryError(
            f"{stage}: peak RSS {peak:.3f} GiB reached {limit_gib:.3f} GiB"
        )


def read_lines(path: Path) -> tuple[str, ...]:
    values = tuple(
        line.strip() for line in path.read_text().splitlines() if line.strip()
    )
    if len(values) != len(set(values)):
        raise ValueError(f"{path}: duplicate labels")
    return values


def read_repeated_lines(path: Path) -> tuple[str, ...]:
    return tuple(line.strip() for line in path.read_text().splitlines() if line.strip())


def read_panel(path: Path) -> tuple[str, ...]:
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle)
        if reader.fieldnames is None or "target_gene" not in reader.fieldnames:
            raise ValueError("panel CSV lacks target_gene")
        values = tuple(str(row["target_gene"]).strip() for row in reader)
    if len(values) != len(set(values)):
        raise ValueError("panel contains duplicate targets")
    return values


def sha256_file(path: Path, block_bytes: int = 8 << 20) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while block := handle.read(block_bytes):
            digest.update(block)
    return digest.hexdigest()


def sha256_lines(values: Iterable[str]) -> str:
    return hashlib.sha256("\n".join(values).encode()).hexdigest()


def load_cache(
    guide_cache: Path,
    feature_cache: Path,
) -> dict[str, object]:
    genes = read_lines(guide_cache / "genes.txt")
    guide_ids = read_lines(guide_cache / "guide_ids.txt")
    guide_targets = read_repeated_lines(guide_cache / "guide_targets.txt")
    target_names = read_lines(guide_cache / "target_names.txt")
    feature_names = read_lines(feature_cache / "feature_names.txt")
    guide_effects = np.load(
        guide_cache / "guide_delta.npy", mmap_mode="r", allow_pickle=False
    )
    pooled = np.load(
        guide_cache / "target_pooled_delta.npy",
        mmap_mode="r",
        allow_pickle=False,
    )
    guide_cells = np.load(
        guide_cache / "guide_cells.npy", mmap_mode="r", allow_pickle=False
    )
    features = np.load(
        feature_cache / "guide_features.npy", mmap_mode="r", allow_pickle=False
    )
    identifiable = np.load(
        feature_cache / "guide_identifiable.npy",
        mmap_mode="r",
        allow_pickle=False,
    )
    if guide_effects.shape != (len(guide_ids), len(genes)):
        raise ValueError("guide effect shape mismatch")
    if pooled.shape != (len(target_names), len(genes)):
        raise ValueError("pooled target shape mismatch")
    if len(guide_targets) != len(guide_ids):
        raise ValueError("guide target length mismatch")
    if features.shape != (len(guide_ids), len(feature_names)):
        raise ValueError("guide feature shape mismatch")
    if guide_cells.shape != identifiable.shape or guide_cells.shape != (
        len(guide_ids),
    ):
        raise ValueError("guide vector shape mismatch")
    return {
        "genes": genes,
        "guide_ids": guide_ids,
        "guide_targets": guide_targets,
        "target_names": target_names,
        "feature_names": feature_names,
        "guide_effects": guide_effects,
        "pooled": pooled,
        "guide_cells": guide_cells,
        "features": features,
        "identifiable": np.asarray(identifiable, dtype=bool),
    }


def clone_model(
    fitted: model.GuideLatentModel,
    candidate: model.Candidate,
) -> model.GuideLatentModel:
    return model.GuideLatentModel(
        candidate=candidate,
        feature_mean=fitted.feature_mean,
        feature_scale=fitted.feature_scale,
        feature_keep=fitted.feature_keep,
        response_mask=fitted.response_mask,
        basis=fitted.basis,
        coefficients=fitted.coefficients,
    )


def candidate_search(
    cache: Mapping[str, object],
    response_mask: np.ndarray,
    splits: Mapping[str, Sequence[str]],
    *,
    max_rss_gib: float,
    fixed_bases: Mapping[int, np.ndarray] | None = None,
) -> tuple[
    model.Candidate,
    list[dict[str, object]],
    dict[int, np.ndarray],
]:
    guide_effects = cache["guide_effects"]
    features = cache["features"]
    guide_targets = cache["guide_targets"]
    guide_cells = cache["guide_cells"]
    identifiable = cache["identifiable"]
    candidates = model.candidate_grid()
    fitted: dict[tuple[int, float], model.GuideLatentModel] = {}
    rows: list[dict[str, object]] = []
    for candidate in candidates:
        key = (candidate.rank, candidate.ridge)
        if key not in fitted:
            fitted[key] = model.fit_latent_model(
                guide_effects,
                features,
                guide_targets,
                identifiable,
                splits["fit"],
                response_mask,
                dataclasses.replace(candidate, eta=1.0, weight_power=0.0),
                seed=SPLIT_SEED,
                fixed_basis=(
                    fixed_bases.get(candidate.rank)
                    if fixed_bases is not None and candidate.rank > 0
                    else None
                ),
            )
        candidate_model = clone_model(fitted[key], candidate)
        truth, baseline, predicted, labels = model.held_guide_predictions(
            candidate_model,
            guide_effects,
            features,
            guide_targets,
            guide_cells,
            identifiable,
            splits["tune"],
            seed=SPLIT_SEED,
        )
        metrics = model.guide_prediction_metrics(truth, baseline, predicted)
        by_count: dict[str, float] = {}
        grouped = model.target_rows(guide_targets)
        if labels:
            base_cos = model.row_cosines(truth, baseline)
            pred_cos = model.row_cosines(truth, predicted)
            for count in sorted({len(grouped[target]) for target in labels}):
                selected = np.asarray(
                    [len(grouped[target]) == count for target in labels], dtype=bool
                )
                by_count[str(count)] = float(
                    np.mean(pred_cos[selected] - base_cos[selected])
                )
        rows.append(
            {
                "candidate": dataclasses.asdict(candidate),
                "tune": metrics,
                "guide_count_cosine_delta": by_count,
                "all_guide_count_strata_nonnegative": all(
                    value >= 0 for value in by_count.values()
                ),
            }
        )
        enforce_memory(max_rss_gib, f"candidate {candidate}")

    eligible = [
        row
        for row in rows
        if not model.Candidate(**row["candidate"]).is_identity
        and float(row["tune"]["cosine_delta"]) > 0
        and float(row["tune"]["nmse_delta"]) <= 0
        and bool(row["all_guide_count_strata_nonnegative"])
    ]
    if not eligible:
        return (
            candidates[0],
            rows,
            {
                rank: fitted[(rank, 0.1)].basis
                for rank in (8, 16, 32)
                if (rank, 0.1) in fitted
            },
        )
    selected_row = max(
        eligible,
        key=lambda row: (
            float(row["tune"]["cosine_delta"]),
            -float(row["tune"]["candidate_nmse"]),
            -model.Candidate(**row["candidate"]).rank,
            -model.Candidate(**row["candidate"]).eta,
        ),
    )
    return (
        model.Candidate(**selected_row["candidate"]),
        rows,
        {
            rank: fitted[(rank, 0.1)].basis
            for rank in (8, 16, 32)
            if (rank, 0.1) in fitted
        },
    )


def refit_and_test(
    cache: Mapping[str, object],
    response_mask: np.ndarray,
    splits: Mapping[str, Sequence[str]],
    selected: model.Candidate,
) -> tuple[model.GuideLatentModel, dict[str, object]]:
    fit_tune = tuple(sorted((*splits["fit"], *splits["tune"])))
    fitted = model.fit_latent_model(
        cache["guide_effects"],
        cache["features"],
        cache["guide_targets"],
        cache["identifiable"],
        fit_tune,
        response_mask,
        selected,
        seed=SPLIT_SEED,
    )
    truth, baseline, predicted, labels = model.held_guide_predictions(
        fitted,
        cache["guide_effects"],
        cache["features"],
        cache["guide_targets"],
        cache["guide_cells"],
        cache["identifiable"],
        splits["test"],
        seed=SPLIT_SEED,
    )
    metrics = model.guide_prediction_metrics(truth, baseline, predicted)
    deltas = model.row_cosines(truth, predicted) - model.row_cosines(truth, baseline)
    if len(labels) >= 2:
        pooled_pos = {
            target: index for index, target in enumerate(cache["target_names"])
        }
        source = np.stack(
            [
                np.asarray(
                    cache["pooled"][pooled_pos[target], response_mask],
                    dtype=np.float64,
                )
                for target in labels
            ]
        )
        clusters = model.source_effect_clusters(source, k=min(8, len(labels) - 1))
        bootstrap = model.clustered_bootstrap(
            deltas,
            clusters,
            repetitions=TEST_BOOTSTRAP_REPETITIONS,
            seed=model.stable_seed(SPLIT_SEED, "test-bootstrap"),
        )
    else:
        clusters = np.zeros(len(labels), dtype=np.int64)
        bootstrap = {
            "repetitions": 0,
            "clusters": int(len(labels)),
            "mean": float(np.mean(deltas)) if len(deltas) else float("nan"),
            "lower_95": float("nan"),
            "upper_95": float("nan"),
        }
    pre_gate = {
        "selected_nonidentity": not selected.is_identity,
        "test_targets": len(labels),
        "test_cosine_delta_positive": float(metrics["cosine_delta"]) > 0,
        "test_nmse_nonworse": float(metrics["nmse_delta"]) <= 0,
        "test_cluster_bootstrap_lower_above_zero": (float(bootstrap["lower_95"]) > 0),
    }
    return fitted, {
        "metrics": metrics,
        "target_order": list(labels),
        "target_order_sha256": sha256_lines(labels),
        "cluster_ids": clusters.tolist(),
        "cluster_bootstrap": bootstrap,
        "pre_gate": pre_gate,
        "passes": all(pre_gate.values()),
    }


def save_model(path: Path, fitted: model.GuideLatentModel) -> None:
    np.savez(
        path,
        feature_mean=fitted.feature_mean,
        feature_scale=fitted.feature_scale,
        feature_keep=fitted.feature_keep,
        response_mask=fitted.response_mask,
        basis=fitted.basis,
        coefficients=fitted.coefficients,
    )


def panel_adjustments(
    fitted: model.GuideLatentModel,
    cache: Mapping[str, object],
    panel: Sequence[str],
    allowed_targets: set[str],
) -> tuple[np.ndarray, np.ndarray, list[dict[str, object]], np.ndarray]:
    grouped = model.target_rows(cache["guide_targets"])
    pooled_pos = {target: index for index, target in enumerate(cache["target_names"])}
    adjustment = np.zeros((len(panel), len(cache["genes"])), dtype=np.float32)
    shrink = np.zeros_like(adjustment)
    audit: list[dict[str, object]] = []
    eligible = np.zeros(len(panel), dtype=bool)
    for panel_index, target in enumerate(panel):
        if target not in allowed_targets:
            audit.append(
                {
                    "target_gene": target,
                    "modified": False,
                    "reason": "absent_from_exact_arm_d_k562_perturbmatch_slot",
                    "guide_count": 0,
                    "identifiable_guides": 0,
                    "correction_norm": 0.0,
                }
            )
            continue
        rows = grouped.get(target)
        pooled_index = pooled_pos.get(target)
        if rows is None or pooled_index is None:
            audit.append(
                {
                    "target_gene": target,
                    "modified": False,
                    "reason": "absent_from_k562",
                    "guide_count": 0,
                    "identifiable_guides": 0,
                    "correction_norm": 0.0,
                }
            )
            continue
        original = np.asarray(cache["pooled"][pooled_index], dtype=np.float64)
        cleaned, detail = model.infer_cleaned_target(
            fitted,
            np.asarray(cache["guide_effects"][rows]),
            np.asarray(cache["features"][rows]),
            np.asarray(cache["guide_cells"][rows]),
            np.asarray(cache["identifiable"][rows]),
            original,
        )
        adjustment[panel_index] = (cleaned - original).astype(np.float32)
        shrink[panel_index] = (
            model.shrink_control(original, cleaned) - original
        ).astype(np.float32)
        eligible[panel_index] = bool(detail["modified"])
        audit.append({"target_gene": target, **detail})
    return adjustment, shrink, audit, eligible


def fit_frozen_nulls(
    output: Path,
    cache: Mapping[str, object],
    splits: Mapping[str, Sequence[str]],
    all_nonpanel: Sequence[str],
    response_mask: np.ndarray,
    panel: Sequence[str],
    allowed_targets: set[str],
    eligible_panel: np.ndarray,
    fit_bases: Mapping[int, np.ndarray],
    *,
    repetitions: int,
    max_rss_gib: float,
) -> dict[str, object]:
    eligible_indices = np.flatnonzero(eligible_panel)
    target_norms = {
        target: float(np.linalg.norm(cache["pooled"][index, response_mask]))
        for index, target in enumerate(cache["target_names"])
    }
    outputs: dict[str, object] = {}
    selection_audit: dict[str, list[dict[str, object]]] = {}
    for family, bundle in (
        ("guide_assignment", True),
        ("seed_features", False),
    ):
        path = output / f"null_{family}_adjustments.npy"
        array = np.lib.format.open_memmap(
            path,
            mode="w+",
            dtype=np.float32,
            shape=(
                repetitions,
                len(eligible_indices),
                len(cache["genes"]),
            ),
        )
        family_selection: list[dict[str, object]] = []
        for replicate in range(repetitions):
            permuted = model.stratified_feature_permutation(
                cache["features"],
                cache["guide_targets"],
                target_norms,
                seed=model.stable_seed(SPLIT_SEED, family, replicate),
                bundle=bundle,
            )
            null_cache = dict(cache)
            null_cache["features"] = permuted
            null_selected, _, _ = candidate_search(
                null_cache,
                response_mask,
                splits,
                max_rss_gib=max_rss_gib,
                fixed_bases=fit_bases,
            )
            _, null_test = refit_and_test(
                null_cache, response_mask, splits, null_selected
            )
            if null_test["passes"]:
                null_model = model.fit_latent_model(
                    cache["guide_effects"],
                    permuted,
                    cache["guide_targets"],
                    cache["identifiable"],
                    all_nonpanel,
                    response_mask,
                    null_selected,
                    seed=SPLIT_SEED,
                    fixed_basis=fit_bases.get(null_selected.rank),
                )
                adjustments, _, _, _ = panel_adjustments(
                    null_model, null_cache, panel, allowed_targets
                )
                array[replicate] = adjustments[eligible_indices]
            else:
                array[replicate] = 0.0
            family_selection.append(
                {
                    "replicate": replicate,
                    "selected": dataclasses.asdict(null_selected),
                    "nonpanel_test_passes": bool(null_test["passes"]),
                    "nonpanel_test_cosine_delta": float(
                        null_test["metrics"]["cosine_delta"]
                    ),
                }
            )
            array.flush()
            enforce_memory(max_rss_gib, f"{family} null {replicate}")
        del array
        selection_audit[family] = family_selection
        outputs[family] = {
            "repetitions": repetitions,
            "path": path.name,
            "sha256": sha256_file(path),
            "complete_grid_selection_repeated": True,
            "nonpanel_test_gate_repeated": True,
            "basis_refit_per_null": False,
            "basis_policy": (
                "frozen real-response basis; gives nulls the same favorable "
                "response subspace and tests only assignment/feature linkage"
            ),
        }
    (output / "null-selection.json").write_text(
        json.dumps(selection_audit, indent=2) + "\n"
    )
    return outputs


def run(args: argparse.Namespace) -> dict[str, object]:
    source_audit = json.loads(args.source_audit.read_text())
    preflight = source_audit["promotion_preflight"]
    if not bool(preflight["passes"]):
        raise RuntimeError(
            "source coverage preflight failed; refusing fit/heavy compute: "
            f"{preflight['observed_panel_arm_d_evaluable_identifiable_multi_guide_targets']} "
            f"< {preflight['minimum_panel_multi_guide_targets']}"
        )
    if args.output.exists() and any(args.output.iterdir()):
        raise FileExistsError(f"refusing to overwrite nonempty {args.output}")
    args.output.mkdir(parents=True)
    cache = load_cache(args.guide_cache, args.feature_cache)
    panel = read_panel(args.panel_csv)
    panel_set = set(panel)
    pm_targets = set(read_lines(args.perturbmatch_perts))
    eligible = model.identifiable_targets(
        cache["guide_targets"], cache["identifiable"], minimum_guides=2
    )
    nonpanel = tuple(target for target in eligible if target not in panel_set)
    panel_eligible = tuple(
        target for target in eligible if target in panel_set and target in pm_targets
    )
    response_mask = model.excluded_coordinate_mask(cache["genes"], (*panel, *nonpanel))
    if len(panel_eligible) < MIN_PANEL_MULTIGUIDE:
        raise RuntimeError(
            f"only {len(panel_eligible)} identifiable multi-guide panel targets; "
            f"require {MIN_PANEL_MULTIGUIDE}"
        )
    splits = model.deterministic_target_split(nonpanel, seed=SPLIT_SEED)
    model.assert_split_contract(splits, panel)
    selected, search_rows, fit_bases = candidate_search(
        cache, response_mask, splits, max_rss_gib=args.max_rss_gib
    )
    _, test_result = refit_and_test(cache, response_mask, splits, selected)
    if not test_result["passes"]:
        result = {
            "experiment": "083-guide-latent-decontamination",
            "status": "stopped_before_panel",
            "decision": "REJECT",
            "reason": "nonpanel held-guide pre-gate failed",
            "selected": dataclasses.asdict(selected),
            "splits": {key: list(value) for key, value in splits.items()},
            "candidate_search": search_rows,
            "nonpanel_test": test_result,
            "held_panel_truth_opened": False,
            "peak_rss_gib": peak_rss_gib(),
        }
        (args.output / "fit-results.json").write_text(
            json.dumps(result, indent=2) + "\n"
        )
        return result

    # The complete non-panel refit occurs before any held panel truth is opened.
    final_model = model.fit_latent_model(
        cache["guide_effects"],
        cache["features"],
        cache["guide_targets"],
        cache["identifiable"],
        nonpanel,
        response_mask,
        selected,
        seed=SPLIT_SEED,
    )
    save_model(args.output / "frozen_model.npz", final_model)
    adjustment, shrink, target_audit, eligible_mask = panel_adjustments(
        final_model, cache, panel, pm_targets
    )
    if int(eligible_mask.sum()) < MIN_PANEL_MULTIGUIDE:
        raise AssertionError("panel coverage changed after preflight")
    np.save(args.output / "panel_adjustment.npy", adjustment, allow_pickle=False)
    np.save(
        args.output / "panel_shrink_control_adjustment.npy",
        shrink,
        allow_pickle=False,
    )
    np.save(
        args.output / "eligible_panel_mask.npy",
        eligible_mask,
        allow_pickle=False,
    )
    (args.output / "panel-target-audit.json").write_text(
        json.dumps(target_audit, indent=2) + "\n"
    )
    nulls = fit_frozen_nulls(
        args.output,
        cache,
        splits,
        nonpanel,
        response_mask,
        panel,
        pm_targets,
        eligible_mask,
        fit_bases,
        repetitions=args.null_repetitions,
        max_rss_gib=args.max_rss_gib,
    )
    frozen_files = (
        "frozen_model.npz",
        "panel_adjustment.npy",
        "panel_shrink_control_adjustment.npy",
        "eligible_panel_mask.npy",
        "panel-target-audit.json",
        "null_guide_assignment_adjustments.npy",
        "null_seed_features_adjustments.npy",
        "null-selection.json",
    )
    result = {
        "experiment": "083-guide-latent-decontamination",
        "status": "frozen_before_panel",
        "decision": "PANEL_EVALUATION_AUTHORIZED",
        "selected": dataclasses.asdict(selected),
        "splits": {
            key: {
                "count": len(value),
                "targets": list(value),
                "sha256": sha256_lines(value),
            }
            for key, value in splits.items()
        },
        "candidate_search": search_rows,
        "nonpanel_test": test_result,
        "panel_coverage": {
            "required": MIN_PANEL_MULTIGUIDE,
            "identifiable_multi_guide": int(eligible_mask.sum()),
            "passes": int(eligible_mask.sum()) >= MIN_PANEL_MULTIGUIDE,
        },
        "nulls": nulls,
        "leakage_audit": {
            "all_300_panel_coordinates_removed_before_fit": True,
            "all_nonpanel_development_target_coordinates_removed_before_fit": True,
            "final_retrieval_axis_excludes_exactly_panel_targets": True,
            "panel_targets_absent_from_fit_tune_test": True,
            "fit_tune_test_disjoint": True,
            "panel_adjustments_frozen_before_panel_truth": True,
            "held_panel_truth_opened": False,
            "post_panel_refit_possible": False,
            "experiment_080_artifacts_used": False,
        },
        "frozen_artifacts": {
            filename: sha256_file(args.output / filename) for filename in frozen_files
        },
        "input_hashes": {
            "source_audit": sha256_file(args.source_audit),
            "panel_csv": sha256_file(args.panel_csv),
            "guide_cache_metadata": sha256_file(args.guide_cache / "metadata.json"),
            "feature_provenance": sha256_file(
                args.feature_cache / "seed-feature-provenance.json"
            ),
        },
        "peak_rss_gib": peak_rss_gib(),
    }
    (args.output / "fit-results.json").write_text(json.dumps(result, indent=2) + "\n")
    (args.output / "panel-evaluation-authorized.json").write_text(
        json.dumps(
            {
                "authorized": True,
                "fit_results_sha256": sha256_file(args.output / "fit-results.json"),
                "frozen_artifacts": result["frozen_artifacts"],
                "no_panel_truth_opened": True,
            },
            indent=2,
        )
        + "\n"
    )
    return result


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--guide-cache", type=Path, required=True)
    parser.add_argument("--feature-cache", type=Path, required=True)
    parser.add_argument("--panel-csv", type=Path, required=True)
    parser.add_argument("--perturbmatch-perts", type=Path, required=True)
    parser.add_argument("--source-audit", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument(
        "--null-repetitions",
        type=int,
        default=NULL_REPETITIONS_PER_FAMILY,
    )
    parser.add_argument("--max-rss-gib", type=float, default=48.0)
    parser.add_argument("--guarded", action="store_true")
    return parser


def main() -> int:
    args = build_parser().parse_args()
    if not args.guarded:
        raise SystemExit("refusing real-data fit without --guarded")
    if args.null_repetitions != NULL_REPETITIONS_PER_FAMILY:
        raise SystemExit(
            f"null repetitions are frozen at {NULL_REPETITIONS_PER_FAMILY}"
        )
    result = run(args)
    print(
        json.dumps(
            {
                "status": result["status"],
                "decision": result["decision"],
                "peak_rss_gib": result["peak_rss_gib"],
            },
            indent=2,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
