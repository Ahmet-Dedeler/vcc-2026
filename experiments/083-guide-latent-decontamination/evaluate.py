#!/usr/bin/env python3
"""Evaluate frozen guide decontamination against held donors and arm D.

This phase cannot fit or refit a model. It accepts only precomputed,
hash-authorized K562 adjustments produced before held panel truth was opened.
Experiment-080 artifacts are intentionally unsupported.
"""

from __future__ import annotations

import argparse
import csv
import gc
import hashlib
import importlib.util
import json
import os
import resource
import sys
from pathlib import Path
from typing import Any, Mapping, Sequence

os.environ.setdefault("OPENBLAS_NUM_THREADS", "8")
os.environ.setdefault("OMP_NUM_THREADS", "8")
os.environ.setdefault("MKL_NUM_THREADS", "8")
os.environ.setdefault("VECLIB_MAXIMUM_THREADS", "8")
os.environ["PYTHONDONTWRITEBYTECODE"] = "1"

import numpy as np

import model


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
EXP068_PATH = ROOT / "experiments/068-seed-offtarget-ablation/run.py"
EXP064_RESULTS = ROOT / "experiments/064-pm-cd4-composition/cloud-output/results.json"
EXP064_PER_TARGET = (
    ROOT / "experiments/064-pm-cd4-composition/cloud-output/per_target_pds.csv"
)
FLAGGED_TARGETS = frozenset({"PQBP1", "SHPRH", "FZD2", "GLYR1"})
MIN_PANEL_MULTIGUIDE = 100
MIN_MACRO_PDS_DELTA = 0.003
BOOTSTRAP_REPETITIONS = 20_000


def load_module(name: str, path: Path):
    spec = importlib.util.spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise ImportError(f"cannot import {path}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


EXP068 = load_module("exp083_exp068", EXP068_PATH)


def peak_rss_gib() -> float:
    value = float(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
    return value / 1024**3 if sys.platform == "darwin" else value / 1024**2


def enforce_memory(limit_gib: float, stage: str) -> None:
    peak = peak_rss_gib()
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


def read_lines(path: Path, *, unique: bool = True) -> tuple[str, ...]:
    values = tuple(
        line.strip() for line in path.read_text().splitlines() if line.strip()
    )
    if unique and len(values) != len(set(values)):
        raise ValueError(f"{path}: duplicate labels")
    return values


def write_csv(path: Path, rows: Sequence[Mapping[str, object]]) -> None:
    if not rows:
        return
    fields: list[str] = []
    for row in rows:
        for field in row:
            if field not in fields:
                fields.append(field)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def verify_frozen_authorization(frozen_dir: Path) -> dict[str, Any]:
    authorization_path = frozen_dir / "panel-evaluation-authorized.json"
    fit_path = frozen_dir / "fit-results.json"
    authorization = json.loads(authorization_path.read_text())
    if authorization.get("authorized") is not True:
        raise ValueError("panel evaluation is not authorized")
    if authorization.get("no_panel_truth_opened") is not True:
        raise ValueError("fit phase did not attest closed panel truth")
    if authorization.get("fit_results_sha256") != sha256_file(fit_path):
        raise ValueError("fit-results hash changed after authorization")
    fit_result = json.loads(fit_path.read_text())
    if fit_result.get("status") != "frozen_before_panel":
        raise ValueError("fit result is not frozen_before_panel")
    if fit_result["leakage_audit"].get("post_panel_refit_possible") is not False:
        raise ValueError("post-panel refit guard is absent")
    if fit_result["leakage_audit"].get("experiment_080_artifacts_used") is not False:
        raise ValueError("experiment-080 artifacts are forbidden")
    for filename, expected in authorization["frozen_artifacts"].items():
        path = frozen_dir / filename
        if sha256_file(path) != expected:
            raise ValueError(f"frozen artifact changed: {filename}")
    return fit_result


def map_adjustment_to_official(
    raw_adjustment: np.ndarray,
    k562_bank: Any,
    official: Sequence[str],
    panel_positions: np.ndarray,
    raw_median_norm: float,
    pm_shared: np.ndarray,
    pm_available: np.ndarray,
) -> np.ndarray:
    mapping = EXP068.map_bank(k562_bank, official)
    raw = np.asarray(raw_adjustment, dtype=np.float64)
    if raw.shape[1] != len(k562_bank.genes):
        raise ValueError("K562 adjustment gene axis mismatch")
    output = np.zeros((len(raw), len(official)), dtype=np.float64)
    output[:, mapping.unique_official_positions] = raw[:, mapping.unique_source_columns]
    output *= EXP068.EFFECT_SCALE / raw_median_norm
    output[:, panel_positions] = 0.0
    shared_mask = np.zeros(len(official), dtype=bool)
    shared_mask[pm_shared] = True
    output[:, shared_mask] *= 1.0 - EXP068.PM_FRACTION
    output[~pm_available] = 0.0
    return output


def build_predictions(
    pm_slot: np.ndarray,
    pm_available: np.ndarray,
    k562_mean: np.ndarray,
    companion: Any,
    cd4: Any,
    operator: Any,
    slot_adjustment: np.ndarray,
    *,
    evec_block_rows: int,
) -> np.ndarray:
    indices = tuple(range(len(pm_slot)))
    baseline = EXP068.build_baseline_rows(
        indices,
        pm_slot + slot_adjustment,
        pm_available,
        k562_mean,
        companion,
        removed=frozenset(),
    )
    transformed = EXP068.transform_rows(
        baseline[:, operator.wcols],
        operator,
        block_rows=evec_block_rows,
    )
    transformed *= operator.preserve_scale
    transformed += EXP068.CD4_FRACTION * cd4.mean_effect[:, operator.wcols]
    return transformed


def score_case(
    predictions: np.ndarray,
    held: Any,
    panel: Sequence[str],
    official: Sequence[str],
    operator: Any,
    held_bank: Any,
) -> tuple[tuple[str, ...], np.ndarray]:
    panel_pos = {target: index for index, target in enumerate(panel)}
    gallery = tuple(target for target in panel if held.available[panel_pos[target]])
    indices = np.asarray([panel_pos[target] for target in gallery], dtype=np.int64)
    held_gene_set = set(held_bank.genes)
    retrieval = np.asarray(
        [
            position
            for position in operator.wcols
            if official[position] in held_gene_set
        ],
        dtype=np.int64,
    )
    lookup = {int(position): index for index, position in enumerate(operator.wcols)}
    local = np.asarray(
        [lookup[int(position)] for position in retrieval], dtype=np.int64
    )
    query = predictions[np.ix_(indices, local)]
    truth = held.effects[np.ix_(indices, retrieval)]
    return gallery, model.pds_per_target(query, truth)


def target_deltas_across_cases(
    case_rows: Sequence[Mapping[str, object]],
) -> dict[str, float]:
    grouped: dict[str, list[float]] = {}
    for row in case_rows:
        grouped.setdefault(str(row["target_gene"]), []).append(float(row["pds_delta"]))
    return {target: float(np.mean(values)) for target, values in grouped.items()}


def fold_source_cluster_bootstrap(
    case_rows: Sequence[Mapping[str, object]],
    baseline_predictions: Mapping[str, np.ndarray],
    panel: Sequence[str],
    prepared: Mapping[str, Any],
    *,
    repetitions: int,
    seed: int,
) -> dict[str, object]:
    """Resample source-only clusters independently inside each outer fold."""
    panel_pos = {target: index for index, target in enumerate(panel)}
    fold_specs = (
        ("orion_hct116", "hct_geometry"),
        ("orion_hek293t", "hek_geometry"),
    )
    fold_data: dict[str, dict[str, object]] = {}
    for held, geometry in fold_specs:
        rows = [
            row
            for row in case_rows
            if row["held_truth"] == held and row["geometry"] == geometry
        ]
        targets = tuple(str(row["target_gene"]) for row in rows)
        expected = tuple(
            target for target in panel if prepared[held].available[panel_pos[target]]
        )
        if targets != expected:
            raise ValueError(f"{held}: rows do not match the complete held gallery")
        indices = np.asarray([panel_pos[target] for target in targets], dtype=np.int64)
        clusters = model.source_effect_clusters(
            baseline_predictions[geometry][indices], k=8
        )
        deltas = np.asarray([float(row["pds_delta"]) for row in rows])
        members = [
            np.flatnonzero(clusters == cluster) for cluster in np.unique(clusters)
        ]
        fold_data[held] = {
            "geometry": geometry,
            "target_order": list(targets),
            "cluster_ids": clusters.tolist(),
            "clusters": len(members),
            "deltas": deltas,
            "members": members,
            "leave_one_cluster_out": {
                str(cluster): float(np.mean(deltas[clusters != cluster]))
                for cluster in np.unique(clusters)
                if np.any(clusters != cluster)
            },
        }
    rng = np.random.default_rng(seed)
    draws = np.empty(repetitions, dtype=np.float64)
    for replicate in range(repetitions):
        fold_means = []
        for held, _ in fold_specs:
            members = fold_data[held]["members"]
            deltas = fold_data[held]["deltas"]
            selected = rng.integers(0, len(members), size=len(members))
            sampled_rows = np.concatenate([members[index] for index in selected])
            fold_means.append(float(np.mean(deltas[sampled_rows])))
        draws[replicate] = float(np.mean(fold_means))
    return {
        "method": (
            "fold-specific reciprocal-8NN source-query clusters; cluster "
            "samples are independent across the two held-donor folds"
        ),
        "uses_held_truth_for_clustering": False,
        "repetitions": repetitions,
        "seed": seed,
        "mean": float(
            np.mean([np.mean(fold_data[held]["deltas"]) for held, _ in fold_specs])
        ),
        "lower_95": float(np.percentile(draws, 2.5)),
        "upper_95": float(np.percentile(draws, 97.5)),
        "folds": {
            held: {
                key: value
                for key, value in details.items()
                if key not in {"members", "deltas"}
            }
            for held, details in fold_data.items()
        },
    }


def evaluate_adjustment(
    raw_adjustment: np.ndarray,
    *,
    panel: Sequence[str],
    official: Sequence[str],
    panel_positions: np.ndarray,
    banks: Mapping[str, Any],
    prepared: Mapping[str, Any],
    cd4: Any,
    pm_slot: np.ndarray,
    pm_available: np.ndarray,
    pm_shared: np.ndarray,
    operators: Mapping[str, Any],
    baseline_predictions: Mapping[str, np.ndarray],
    evec_block_rows: int,
    max_rss_gib: float,
    label: str,
) -> tuple[dict[str, float], list[dict[str, object]], dict[str, np.ndarray]]:
    mapped = map_adjustment_to_official(
        raw_adjustment,
        banks["k562_gwps"],
        official,
        panel_positions,
        prepared["k562_gwps"].raw_median_norm,
        pm_shared,
        pm_available,
    )
    rows: list[dict[str, object]] = []
    predictions: dict[str, np.ndarray] = {}
    cases = (
        ("orion_hct116", "hct_geometry"),
        ("orion_hek293t", "hek_geometry"),
        ("kolf", "hct_geometry"),
        ("kolf", "hek_geometry"),
    )
    per_case_mean: dict[str, float] = {}
    for truth_name, config in cases:
        companion_name = "orion_hek293t" if config == "hct_geometry" else "orion_hct116"
        # This is a correction to the K562/PerturbMatch source slot.
        # build_baseline_rows performs the equal-source averaging exactly once.
        slot_adjustment = mapped
        candidate = build_predictions(
            pm_slot,
            pm_available,
            prepared["k562_gwps"].mean_on_axis,
            prepared[companion_name],
            cd4,
            operators[config],
            slot_adjustment,
            evec_block_rows=evec_block_rows,
        )
        if config not in predictions:
            predictions[config] = candidate
        gallery, baseline_pds = score_case(
            baseline_predictions[config],
            prepared[truth_name],
            panel,
            official,
            operators[config],
            banks[truth_name],
        )
        candidate_gallery, candidate_pds = score_case(
            candidate,
            prepared[truth_name],
            panel,
            official,
            operators[config],
            banks[truth_name],
        )
        if gallery != candidate_gallery:
            raise AssertionError("candidate changed held gallery")
        case_id = f"{truth_name}__{config}"
        per_case_mean[case_id] = float(np.mean(candidate_pds - baseline_pds))
        for target, before, after in zip(
            gallery, baseline_pds, candidate_pds, strict=True
        ):
            rows.append(
                {
                    "variant": label,
                    "case_id": case_id,
                    "held_truth": truth_name,
                    "geometry": config,
                    "target_gene": target,
                    "baseline_pds": float(before),
                    "candidate_pds": float(after),
                    "pds_delta": float(after - before),
                }
            )
        del candidate
        gc.collect()
        enforce_memory(max_rss_gib, f"{label} {case_id}")
    donor_means = {
        "orion_hct116": per_case_mean["orion_hct116__hct_geometry"],
        "orion_hek293t": per_case_mean["orion_hek293t__hek_geometry"],
        "kolf": float(
            np.mean(
                [
                    per_case_mean["kolf__hct_geometry"],
                    per_case_mean["kolf__hek_geometry"],
                ]
            )
        ),
    }
    donor_means["macro"] = float(np.mean(list(donor_means.values())))
    return donor_means, rows, predictions


def reconstruct_arm_d(
    panel: Sequence[str],
    official: Sequence[str],
    banks: Mapping[str, Any],
    prepared: Mapping[str, Any],
    cd4: Any,
    pm_slot: np.ndarray,
    pm_available: np.ndarray,
    operators: Mapping[str, Any],
    *,
    evec_block_rows: int,
) -> dict[str, np.ndarray]:
    output = {}
    for config, companion_name in (
        ("hct_geometry", "orion_hek293t"),
        ("hek_geometry", "orion_hct116"),
    ):
        output[config] = build_predictions(
            pm_slot,
            pm_available,
            prepared["k562_gwps"].mean_on_axis,
            prepared[companion_name],
            cd4,
            operators[config],
            np.zeros((len(panel), len(official)), dtype=np.float64),
            evec_block_rows=evec_block_rows,
        )
    return output


def arm_d_fidelity(
    baseline_predictions: Mapping[str, np.ndarray],
    panel: Sequence[str],
    official: Sequence[str],
    prepared: Mapping[str, Any],
    banks: Mapping[str, Any],
    operators: Mapping[str, Any],
) -> dict[str, object]:
    experiment = json.loads(EXP064_RESULTS.read_text())
    reference_targets = tuple(experiment["axis"]["target_order"])
    panel_pos = {target: index for index, target in enumerate(panel)}
    indices = np.asarray(
        [panel_pos[target] for target in reference_targets], dtype=np.int64
    )
    expected_rows = list(csv.DictReader(EXP064_PER_TARGET.open(newline="")))
    expected = {
        (str(row["held"]), str(row["target_gene"])): float(row["pds_D"])
        for row in expected_rows
    }
    errors = []
    rows = []
    for truth_name, config in (
        ("orion_hct116", "hct_geometry"),
        ("orion_hek293t", "hek_geometry"),
    ):
        held_genes = set(banks[truth_name].genes)
        retrieval = np.asarray(
            [
                position
                for position in operators[config].wcols
                if official[position] in held_genes
            ],
            dtype=np.int64,
        )
        lookup = {
            int(position): index
            for index, position in enumerate(operators[config].wcols)
        }
        local = np.asarray([lookup[int(value)] for value in retrieval], dtype=np.int64)
        query = baseline_predictions[config][np.ix_(indices, local)]
        truth = prepared[truth_name].effects[np.ix_(indices, retrieval)]
        observed = model.pds_per_target(query, truth)
        for target, value in zip(reference_targets, observed, strict=True):
            wanted = expected[(truth_name, target)]
            error = abs(float(value) - wanted)
            errors.append(error)
            rows.append(
                {
                    "held": truth_name,
                    "target_gene": target,
                    "expected_pds_D": wanted,
                    "observed_pds_D": float(value),
                    "absolute_error": error,
                }
            )
    maximum = max(errors)
    return {
        "rows": rows,
        "maximum_absolute_error": maximum,
        "exact": maximum == 0.0,
    }


def run(args: argparse.Namespace) -> dict[str, object]:
    fit_result = verify_frozen_authorization(args.frozen_dir)
    source_audit = json.loads(args.source_audit.read_text())
    if not source_audit["promotion_preflight"]["passes"]:
        raise RuntimeError("source coverage preflight is not passing")
    if (
        int(fit_result["panel_coverage"]["identifiable_multi_guide"])
        < MIN_PANEL_MULTIGUIDE
    ):
        raise RuntimeError("frozen fit lacks required panel coverage")
    if args.output.exists() and any(args.output.iterdir()):
        raise FileExistsError(f"refusing to overwrite nonempty {args.output}")
    args.output.mkdir(parents=True)

    EXP068.MAX_RSS_GIB = args.max_rss_gib
    EXP068.STOP_RSS_GIB = args.max_rss_gib
    panel = EXP068.read_csv_column(args.panel_csv, "target_gene")
    official = EXP068.read_csv_column(args.genes_csv, "gene_name")
    panel_pos = {target: index for index, target in enumerate(panel)}
    official_pos = {gene: index for index, gene in enumerate(official)}
    panel_positions = np.asarray(
        [official_pos[target] for target in panel], dtype=np.int64
    )
    panel_set = set(panel)
    guide_genes = read_lines(args.guide_cache / "genes.txt")
    raw_adjustment = np.load(
        args.frozen_dir / "panel_adjustment.npy",
        mmap_mode="r",
        allow_pickle=False,
    )
    shrink_adjustment = np.load(
        args.frozen_dir / "panel_shrink_control_adjustment.npy",
        mmap_mode="r",
        allow_pickle=False,
    )
    eligible = np.asarray(
        np.load(
            args.frozen_dir / "eligible_panel_mask.npy",
            mmap_mode="r",
            allow_pickle=False,
        ),
        dtype=bool,
    )
    if raw_adjustment.shape != (len(panel), len(guide_genes)):
        raise ValueError("frozen adjustment shape mismatch")
    if np.any(np.asarray(raw_adjustment[~eligible]) != 0):
        raise ValueError("unsupported target was modified")
    local_panel_positions = [
        guide_genes.index(target) for target in panel if target in set(guide_genes)
    ]
    if local_panel_positions and np.any(
        np.asarray(raw_adjustment[:, local_panel_positions]) != 0
    ):
        raise ValueError("a VCC target coordinate was modified")

    banks = {
        "k562_gwps": EXP068.Bank.load("k562_gwps", args.k562),
        "orion_hct116": EXP068.Bank.load("orion_hct116", args.hct116),
        "orion_hek293t": EXP068.Bank.load("orion_hek293t", args.hek293t),
        "perturbmatch": EXP068.Bank.load("perturbmatch", args.perturbmatch),
        "kolf": EXP068.Bank.load("kolf", args.kolf),
    }
    cd4_banks = [
        EXP068.Bank.load(name, args.cd4_root / name) for name in EXP068.CD4_NAMES
    ]
    prepared = {
        name: EXP068.prepare_main(
            bank,
            panel,
            official,
            panel_positions,
            chunk_rows=args.chunk_rows,
        )
        for name, bank in banks.items()
    }
    cd4 = EXP068.prepare_cd4(
        cd4_banks,
        panel,
        official,
        panel_set,
        chunk_rows=args.chunk_rows,
    )
    pm_slot, pm_available, pm_shared = EXP068.prepare_pm_slot(
        prepared["k562_gwps"],
        prepared["perturbmatch"],
        banks["k562_gwps"],
        banks["perturbmatch"],
        official,
        panel_set,
    )
    exp064 = json.loads(EXP064_RESULTS.read_text())
    preserve = {
        "hct_geometry": float(
            next(fold for fold in exp064["folds"] if fold["held"] == "orion_hct116")[
                "production_preserve"
            ]["B"]["applied_scale"]
        ),
        "hek_geometry": float(
            next(fold for fold in exp064["folds"] if fold["held"] == "orion_hek293t")[
                "production_preserve"
            ]["B"]["applied_scale"]
        ),
    }
    operators = {
        "hct_geometry": EXP068.load_operator(
            "hct_geometry",
            args.operator_dir / "held_orion_hct116.npz",
            [banks["k562_gwps"], banks["orion_hek293t"]],
            official,
            panel_set,
            preserve["hct_geometry"],
        ),
        "hek_geometry": EXP068.load_operator(
            "hek_geometry",
            args.operator_dir / "held_orion_hek293t.npz",
            [banks["k562_gwps"], banks["orion_hct116"]],
            official,
            panel_set,
            preserve["hek_geometry"],
        ),
    }
    if not np.array_equal(
        operators["hct_geometry"].wcols, operators["hek_geometry"].wcols
    ):
        raise ValueError("frozen arm-D operators have different axes")
    baseline_predictions = reconstruct_arm_d(
        panel,
        official,
        banks,
        prepared,
        cd4,
        pm_slot,
        pm_available,
        operators,
        evec_block_rows=args.evec_block_rows,
    )
    fidelity = arm_d_fidelity(
        baseline_predictions, panel, official, prepared, banks, operators
    )
    if not fidelity["exact"]:
        raise ValueError(
            "arm-D baseline reconstruction is not bit-exact in PDS: "
            f"{fidelity['maximum_absolute_error']}"
        )

    candidate_summary, candidate_rows, _ = evaluate_adjustment(
        raw_adjustment,
        panel=panel,
        official=official,
        panel_positions=panel_positions,
        banks=banks,
        prepared=prepared,
        cd4=cd4,
        pm_slot=pm_slot,
        pm_available=pm_available,
        pm_shared=pm_shared,
        operators=operators,
        baseline_predictions=baseline_predictions,
        evec_block_rows=args.evec_block_rows,
        max_rss_gib=args.max_rss_gib,
        label="candidate",
    )
    shrink_summary, shrink_rows, _ = evaluate_adjustment(
        shrink_adjustment,
        panel=panel,
        official=official,
        panel_positions=panel_positions,
        banks=banks,
        prepared=prepared,
        cd4=cd4,
        pm_slot=pm_slot,
        pm_available=pm_available,
        pm_shared=pm_shared,
        operators=operators,
        baseline_predictions=baseline_predictions,
        evec_block_rows=args.evec_block_rows,
        max_rss_gib=args.max_rss_gib,
        label="shrink_control",
    )

    case_rows = candidate_rows
    per_target = target_deltas_across_cases(case_rows)
    bootstrap = fold_source_cluster_bootstrap(
        case_rows,
        baseline_predictions,
        panel,
        prepared,
        repetitions=BOOTSTRAP_REPETITIONS,
        seed=model.stable_seed(model.SPLIT_SEED, "panel-bootstrap"),
    )
    all_leave_one_cluster_out = [
        float(value)
        for fold in bootstrap["folds"].values()
        for value in fold["leave_one_cluster_out"].values()
    ]

    null_summaries: dict[str, object] = {}
    eligible_indices = np.flatnonzero(eligible)
    for family in ("guide_assignment", "seed_features"):
        null_path = args.frozen_dir / f"null_{family}_adjustments.npy"
        null_values = np.load(null_path, mmap_mode="r", allow_pickle=False)
        if null_values.shape[1:] != (len(eligible_indices), len(guide_genes)):
            raise ValueError(f"{family}: null adjustment shape mismatch")
        draws = np.empty(len(null_values), dtype=np.float64)
        for replicate in range(len(null_values)):
            full = np.zeros_like(raw_adjustment, dtype=np.float32)
            full[eligible_indices] = null_values[replicate]
            summary, _, _ = evaluate_adjustment(
                full,
                panel=panel,
                official=official,
                panel_positions=panel_positions,
                banks=banks,
                prepared=prepared,
                cd4=cd4,
                pm_slot=pm_slot,
                pm_available=pm_available,
                pm_shared=pm_shared,
                operators=operators,
                baseline_predictions=baseline_predictions,
                evec_block_rows=args.evec_block_rows,
                max_rss_gib=args.max_rss_gib,
                label=f"{family}_{replicate}",
            )
            draws[replicate] = summary["macro"]
        null_summaries[family] = {
            "repetitions": len(draws),
            "mean": float(np.mean(draws)),
            "p95": float(np.percentile(draws, 95)),
            "one_sided_p_ge_observed": float(
                (1 + np.sum(draws >= candidate_summary["macro"])) / (len(draws) + 1)
            ),
        }

    flagged_values = [
        delta for target, delta in per_target.items() if target in FLAGGED_TARGETS
    ]
    unflagged_multiguide_values = [
        per_target[target]
        for target in panel
        if eligible[panel_pos[target]]
        and target not in FLAGGED_TARGETS
        and target in per_target
    ]
    decision_checks = {
        "coverage_at_least_100": int(eligible.sum()) >= MIN_PANEL_MULTIGUIDE,
        "nonnegative_each_held_donor": all(
            candidate_summary[name] >= 0
            for name in ("orion_hct116", "orion_hek293t", "kolf")
        ),
        "macro_delta_at_least_0_003": (
            candidate_summary["macro"] >= MIN_MACRO_PDS_DELTA
        ),
        "cluster_bootstrap_lower_above_zero": bootstrap["lower_95"] > 0,
        "above_guide_assignment_null_p95": (
            candidate_summary["macro"] > null_summaries["guide_assignment"]["p95"]
        ),
        "above_seed_feature_null_p95": (
            candidate_summary["macro"] > null_summaries["seed_features"]["p95"]
        ),
        "target_coverage_unchanged": True,
        "all_unsupported_targets_exactly_unchanged": bool(
            np.all(np.asarray(raw_adjustment[~eligible]) == 0)
        ),
        "flagged_seed_evidence_stratum_nonnegative": (
            bool(flagged_values) and float(np.mean(flagged_values)) >= 0
        ),
        "unflagged_multiguide_stratum_nonnegative": (
            bool(unflagged_multiguide_values)
            and float(np.mean(unflagged_multiguide_values)) >= 0
        ),
        "not_explained_by_global_shrinkage": (
            candidate_summary["macro"] > shrink_summary["macro"]
        ),
        "all_leave_one_cluster_out_nonnegative": all(
            value >= 0 for value in all_leave_one_cluster_out
        ),
        "exact_arm_d_baseline_fidelity": fidelity["exact"],
    }
    decision = "PROMOTE" if all(decision_checks.values()) else "REJECT"
    result = {
        "experiment": "083-guide-latent-decontamination",
        "status": "completed",
        "decision": decision,
        "candidate": candidate_summary,
        "shrink_control": shrink_summary,
        "cluster_bootstrap": bootstrap,
        "source_derived_clusters": bootstrap["folds"],
        "strata": {
            "flagged_seed_evidence": {
                "targets": len(flagged_values),
                "mean_delta": (
                    float(np.mean(flagged_values)) if flagged_values else None
                ),
            },
            "unflagged_identifiable_multi_guide": {
                "targets": len(unflagged_multiguide_values),
                "mean_delta": (
                    float(np.mean(unflagged_multiguide_values))
                    if unflagged_multiguide_values
                    else None
                ),
            },
        },
        "nulls": null_summaries,
        "decision_checks": decision_checks,
        "arm_d_fidelity": {
            key: value for key, value in fidelity.items() if key != "rows"
        },
        "leakage_audit": {
            "experiment_080_artifacts_used": False,
            "held_truth_opened_only_after_frozen_hash_authorization": True,
            "model_fit_or_refit_in_evaluation": False,
            "all_vcc_target_coordinates_excluded": True,
            "complete_held_truth_galleries_used": True,
            "whole_source_rows_deleted": False,
            "unsupported_targets_preserved_exactly": True,
            "operators_frozen_from_experiment_062": True,
            "comparator_exact_experiment_064_arm_d": True,
        },
        "input_hashes": {
            "fit_results": sha256_file(args.frozen_dir / "fit-results.json"),
            "authorization": sha256_file(
                args.frozen_dir / "panel-evaluation-authorized.json"
            ),
            "source_audit": sha256_file(args.source_audit),
            "exp068_run": sha256_file(EXP068_PATH),
            "exp064_results": sha256_file(EXP064_RESULTS),
            "exp064_per_target": sha256_file(EXP064_PER_TARGET),
            "operators": {
                name: sha256_file(operator.path) for name, operator in operators.items()
            },
        },
        "peak_rss_gib": peak_rss_gib(),
    }
    write_csv(args.output / "per_target_pds.csv", candidate_rows)
    write_csv(args.output / "shrink_control_pds.csv", shrink_rows)
    write_csv(args.output / "arm_d_reconstruction.csv", fidelity["rows"])
    (args.output / "results.json").write_text(json.dumps(result, indent=2) + "\n")
    return result


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--frozen-dir", type=Path, required=True)
    parser.add_argument("--source-audit", type=Path, required=True)
    parser.add_argument("--guide-cache", type=Path, required=True)
    parser.add_argument("--genes-csv", type=Path, required=True)
    parser.add_argument("--panel-csv", type=Path, required=True)
    parser.add_argument("--k562", type=Path, required=True)
    parser.add_argument("--hct116", type=Path, required=True)
    parser.add_argument("--hek293t", type=Path, required=True)
    parser.add_argument("--perturbmatch", type=Path, required=True)
    parser.add_argument("--kolf", type=Path, required=True)
    parser.add_argument("--cd4-root", type=Path, required=True)
    parser.add_argument("--operator-dir", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--chunk-rows", type=int, default=4)
    parser.add_argument("--evec-block-rows", type=int, default=128)
    parser.add_argument("--max-rss-gib", type=float, default=48.0)
    parser.add_argument("--guarded", action="store_true")
    return parser


def main() -> int:
    args = build_parser().parse_args()
    if not args.guarded:
        raise SystemExit("refusing held-panel evaluation without --guarded")
    result = run(args)
    print(
        json.dumps(
            {
                "decision": result["decision"],
                "candidate": result["candidate"],
                "decision_checks": result["decision_checks"],
                "peak_rss_gib": result["peak_rss_gib"],
            },
            indent=2,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
