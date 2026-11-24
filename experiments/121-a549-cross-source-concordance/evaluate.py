#!/usr/bin/env python3
"""Run the one-shot preregistered Experiment 121 source-only evaluation."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import os
import shutil
import tempfile
from pathlib import Path
from typing import Any, Iterable, Sequence

import numpy as np

import freeze_plans


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
PROTOCOL_PATH = HERE / "preregistered_protocol.json"
AUTHORIZATION_PATH = HERE / "AUTHORIZATION_MANIFEST.json"
COMMON_TARGETS_PATH = HERE / "common_targets.txt"
FREEZE_DIR = HERE / "work" / "preregistered-freeze-v1"
DEFAULT_OUTPUT = HERE / "work" / "source-only-evaluation-v1"


class EvaluationError(RuntimeError):
    """A frozen input, metric, safety, or output invariant failed."""


def require(condition: bool, message: str) -> None:
    if not condition:
        raise EvaluationError(message)


def sha256_file(path: Path, block_bytes: int = 8 << 20) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while block := handle.read(block_bytes):
            digest.update(block)
    return digest.hexdigest()


def sha256_lines(values: Iterable[str]) -> str:
    payload = "".join(f"{value}\n" for value in values).encode()
    return hashlib.sha256(payload).hexdigest()


def json_bytes(payload: dict[str, Any]) -> bytes:
    return (json.dumps(payload, indent=2, sort_keys=True) + "\n").encode()


def read_json(path: Path) -> dict[str, Any]:
    value = json.loads(path.read_text(encoding="utf-8"))
    require(isinstance(value, dict), f"JSON root is not an object: {path}")
    return value


def read_text_axis(path: Path) -> tuple[str, ...]:
    require(path.read_bytes().endswith(b"\n"), f"axis lacks final newline: {path}")
    values = tuple(path.read_text(encoding="utf-8").splitlines())
    require(all(value and value == value.strip() for value in values), "invalid axis")
    require(len(values) == len(set(values)), f"duplicate labels in {path}")
    return values


def read_csv_axis(path: Path, column: str) -> tuple[str, ...]:
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        require(reader.fieldnames == [column], f"header changed: {path}")
        values = tuple(str(row[column]) for row in reader)
    require(all(value and value == value.strip() for value in values), "invalid axis")
    require(len(values) == len(set(values)), f"duplicate labels in {path}")
    return values


def safe_output(path: Path) -> Path:
    resolved = path.expanduser().resolve()
    require(
        resolved != HERE and resolved.is_relative_to(HERE),
        f"output must be below {HERE}",
    )
    return resolved


def authenticate_source_files(authorization: dict[str, Any]) -> None:
    for source in authorization["source_bindings"].values():
        records = (source,) if "path" in source else source.values()
        for record in records:
            path = ROOT / record["path"]
            require(path.is_file() and not path.is_symlink(), f"missing input: {path}")
            require(sha256_file(path) == record["sha256"], f"hash changed: {path}")


def target_indices(axis: Sequence[str], targets: Sequence[str]) -> np.ndarray:
    positions = {label: index for index, label in enumerate(axis)}
    require(all(target in positions for target in targets), "target missing from axis")
    return np.asarray([positions[target] for target in targets], dtype=np.int64)


def build_response_mask(
    genes: Sequence[str],
    altair_measured: np.ndarray,
    myllia_measured: np.ndarray,
    panel: Sequence[str],
    targets: Sequence[str],
) -> np.ndarray:
    require(
        altair_measured.shape == myllia_measured.shape,
        "selected measured-mask shapes differ",
    )
    require(
        altair_measured.ndim == 2
        and altair_measured.shape[1] == len(genes),
        "selected measured-mask shape changed",
    )
    mask = np.all(altair_measured, axis=0) & np.all(myllia_measured, axis=0)
    excluded = set(panel) | set(targets)
    mask &= np.asarray([gene not in excluded for gene in genes], dtype=bool)
    return mask


def cosine_matrix(
    predictor: np.ndarray,
    reference: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    predictor = np.asarray(predictor, dtype=np.float64)
    reference = np.asarray(reference, dtype=np.float64)
    require(predictor.ndim == 2, "predictor must be two-dimensional")
    require(reference.ndim == 2, "reference must be two-dimensional")
    require(predictor.shape[1] == reference.shape[1], "coordinate shapes differ")
    predictor_norm = np.linalg.norm(predictor, axis=1)
    reference_norm = np.linalg.norm(reference, axis=1)
    denominator = predictor_norm[:, None] * reference_norm[None, :]
    similarity = np.zeros((predictor.shape[0], reference.shape[0]), dtype=np.float64)
    np.divide(
        predictor @ reference.T,
        denominator,
        out=similarity,
        where=denominator > 0,
    )
    return similarity, predictor_norm == 0, reference_norm == 0


def pds_from_assignment(
    similarity: np.ndarray,
    assignment: np.ndarray,
    predictor_zero: np.ndarray,
    reference_zero: np.ndarray,
) -> np.ndarray:
    similarity = np.asarray(similarity, dtype=np.float64)
    assignment = np.asarray(assignment, dtype=np.int64)
    count = similarity.shape[0]
    require(similarity.shape == (count, count), "PDS gallery must be square")
    require(assignment.shape == (count,), "assignment shape changed")
    require(
        np.array_equal(np.sort(assignment), np.arange(count)),
        "assignment must be a permutation",
    )
    scores = np.empty(count, dtype=np.float64)
    for row in range(count):
        own = int(assignment[row])
        if predictor_zero[row] or reference_zero[own]:
            scores[row] = 0.5
            continue
        own_similarity = similarity[row, own]
        greater = int(np.count_nonzero(similarity[row] > own_similarity))
        equal = int(np.count_nonzero(similarity[row] == own_similarity))
        midrank = 1.0 + greater + 0.5 * (equal - 1)
        scores[row] = 1.0 - (midrank - 1.0) / (count - 1)
    return scores


def sign_metrics(
    predictor: np.ndarray,
    reference: np.ndarray,
    fraction: float = 0.2,
) -> tuple[np.ndarray, np.ndarray, int]:
    predictor = np.asarray(predictor, dtype=np.float64)
    reference = np.asarray(reference, dtype=np.float64)
    require(predictor.shape == reference.shape, "sign metric shapes differ")
    require(predictor.ndim == 2, "sign metric inputs must be matrices")
    selected_count = max(1, math.ceil(fraction * predictor.shape[1]))
    canonical_index = np.arange(predictor.shape[1])
    precision = np.empty(predictor.shape[0], dtype=np.float64)
    credit = np.empty(predictor.shape[0], dtype=np.float64)
    for row in range(predictor.shape[0]):
        order = np.lexsort((canonical_index, -np.abs(predictor[row])))
        selected = order[:selected_count]
        predictor_sign = np.sign(predictor[row, selected])
        reference_sign = np.sign(reference[row, selected])
        correct = (
            (predictor_sign != 0)
            & (reference_sign != 0)
            & (predictor_sign == reference_sign)
        )
        precision[row] = float(np.mean(correct))
        weights = np.abs(predictor[row, selected])
        denominator = float(np.sum(weights))
        if denominator == 0:
            credit[row] = 0.0
        else:
            credit[row] = float(
                np.sum(weights * predictor_sign * reference_sign) / denominator
            )
    return precision, credit, selected_count


def inverse_permutation(permutation: np.ndarray) -> np.ndarray:
    inverse = np.empty_like(permutation)
    inverse[permutation] = np.arange(len(permutation), dtype=permutation.dtype)
    return inverse


def null_macro_pds(
    altair_to_myllia_similarity: np.ndarray,
    myllia_to_altair_similarity: np.ndarray,
    altair_zero: np.ndarray,
    myllia_zero: np.ndarray,
    derangements: np.ndarray,
) -> np.ndarray:
    values = np.empty(len(derangements), dtype=np.float64)
    for index, permutation in enumerate(derangements):
        inverse = inverse_permutation(permutation)
        forward = pds_from_assignment(
            altair_to_myllia_similarity,
            permutation,
            altair_zero,
            myllia_zero,
        )
        reverse = pds_from_assignment(
            myllia_to_altair_similarity,
            inverse,
            myllia_zero,
            altair_zero,
        )
        values[index] = 0.5 * (float(np.mean(forward)) + float(np.mean(reverse)))
    return values


def bootstrap_macro_pds(
    forward_pds: np.ndarray,
    reverse_pds: np.ndarray,
    bootstrap_indices: np.ndarray,
) -> np.ndarray:
    forward = np.mean(forward_pds[bootstrap_indices], axis=1)
    reverse = np.mean(reverse_pds[bootstrap_indices], axis=1)
    return 0.5 * (forward + reverse)


def write_per_target(
    path: Path,
    targets: Sequence[str],
    own_cosine: np.ndarray,
    forward_pds: np.ndarray,
    reverse_pds: np.ndarray,
    forward_precision: np.ndarray,
    reverse_precision: np.ndarray,
    forward_credit: np.ndarray,
    reverse_credit: np.ndarray,
) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(
            (
                "target",
                "own_target_cosine",
                "altair_to_myllia_pds",
                "myllia_to_altair_pds",
                "altair_to_myllia_top20_sign_precision",
                "myllia_to_altair_top20_sign_precision",
                "altair_to_myllia_directional_sign_credit",
                "myllia_to_altair_directional_sign_credit",
            )
        )
        for values in zip(
            targets,
            own_cosine,
            forward_pds,
            reverse_pds,
            forward_precision,
            reverse_precision,
            forward_credit,
            reverse_credit,
            strict=True,
        ):
            writer.writerow(
                (values[0],) + tuple(f"{float(value):.17g}" for value in values[1:])
            )


def prepare_inputs() -> dict[str, Any]:
    freeze_plans.verify(FREEZE_DIR)
    protocol = read_json(PROTOCOL_PATH)
    authorization = read_json(AUTHORIZATION_PATH)
    require(
        protocol["status"]
        == "frozen_before_either_source_expression_array_was_opened",
        "protocol status changed",
    )
    authenticate_source_files(authorization)

    targets = read_text_axis(COMMON_TARGETS_PATH)
    expected_targets = tuple(protocol["target_axis"]["labels"])
    require(targets == expected_targets, "target axis changed")
    altair_perts = read_text_axis(
        ROOT / authorization["source_bindings"]["altair_experiment_119"][
            "perturbations"
        ]["path"]
    )
    myllia_perts = read_text_axis(
        ROOT / authorization["source_bindings"]["myllia_experiment_120"][
            "perturbations"
        ]["path"]
    )
    require(
        tuple(sorted(set(altair_perts) & set(myllia_perts))) == targets,
        "exact target intersection changed",
    )
    altair_genes = read_text_axis(
        ROOT
        / authorization["source_bindings"]["altair_experiment_119"]["genes"]["path"]
    )
    myllia_genes = read_text_axis(
        ROOT
        / authorization["source_bindings"]["myllia_experiment_120"]["genes"]["path"]
    )
    require(altair_genes == myllia_genes, "canonical source gene axes differ")
    canonical = read_csv_axis(
        ROOT / authorization["source_bindings"]["canonical_gene_axis"]["path"],
        "gene_name",
    )
    require(altair_genes == canonical, "source and canonical gene axes differ")
    panel = read_csv_axis(
        ROOT / authorization["source_bindings"]["vcc_panel_axis"]["path"],
        "target_gene",
    )
    require(len(panel) == 300, "VCC panel count changed")
    return {
        "protocol": protocol,
        "authorization": authorization,
        "targets": targets,
        "altair_perts": altair_perts,
        "myllia_perts": myllia_perts,
        "genes": canonical,
        "panel": panel,
    }


def run(output: Path = DEFAULT_OUTPUT) -> Path:
    output = safe_output(output)
    require(not output.exists(), f"refusing to overwrite evaluation: {output}")
    inputs = prepare_inputs()
    authorization = inputs["authorization"]
    targets = inputs["targets"]
    genes = inputs["genes"]
    panel = inputs["panel"]

    altair_record = authorization["source_bindings"]["altair_experiment_119"]
    myllia_record = authorization["source_bindings"]["myllia_experiment_120"]
    altair_values = np.load(
        ROOT / altair_record["expression_array"]["path"],
        mmap_mode="r",
        allow_pickle=False,
    )
    altair_measured = np.load(
        ROOT / altair_record["measured_mask"]["path"],
        mmap_mode="r",
        allow_pickle=False,
    )
    myllia_values = np.load(
        ROOT / myllia_record["expression_array"]["path"],
        mmap_mode="r",
        allow_pickle=False,
    )
    myllia_measured = np.load(
        ROOT / myllia_record["measured_mask"]["path"],
        mmap_mode="r",
        allow_pickle=False,
    )

    altair_rows = target_indices(inputs["altair_perts"], targets)
    myllia_rows = target_indices(inputs["myllia_perts"], targets)
    altair_selected_mask = np.asarray(altair_measured[altair_rows], dtype=bool)
    myllia_selected_mask = np.asarray(myllia_measured[myllia_rows], dtype=bool)
    response_mask = build_response_mask(
        genes,
        altair_selected_mask,
        myllia_selected_mask,
        panel,
        targets,
    )
    require(np.count_nonzero(response_mask) > 0, "response mask is empty")
    altair = np.asarray(
        altair_values[np.ix_(altair_rows, np.flatnonzero(response_mask))],
        dtype=np.float64,
    )
    myllia = np.asarray(
        myllia_values[np.ix_(myllia_rows, np.flatnonzero(response_mask))],
        dtype=np.float64,
    )
    require(np.all(np.isfinite(altair)), "Altair selected values are nonfinite")
    require(np.all(np.isfinite(myllia)), "Myllia selected values are nonfinite")

    forward_similarity, altair_zero, myllia_zero = cosine_matrix(altair, myllia)
    reverse_similarity, myllia_zero_again, altair_zero_again = cosine_matrix(
        myllia,
        altair,
    )
    require(np.array_equal(myllia_zero, myllia_zero_again), "zero-row drift")
    require(np.array_equal(altair_zero, altair_zero_again), "zero-row drift")
    require(
        np.allclose(forward_similarity.T, reverse_similarity, atol=1e-15, rtol=0),
        "bidirectional cosine matrix lost symmetry",
    )

    identity = np.arange(len(targets), dtype=np.int64)
    forward_pds = pds_from_assignment(
        forward_similarity,
        identity,
        altair_zero,
        myllia_zero,
    )
    reverse_pds = pds_from_assignment(
        reverse_similarity,
        identity,
        myllia_zero,
        altair_zero,
    )
    own_cosine = np.diag(forward_similarity)
    forward_precision, forward_credit, selected_count = sign_metrics(
        altair,
        myllia,
    )
    reverse_precision, reverse_credit, reverse_selected_count = sign_metrics(
        myllia,
        altair,
    )
    require(selected_count == reverse_selected_count, "sign selection count drift")

    derangements = np.load(
        FREEZE_DIR / "target-derangements.npy",
        allow_pickle=False,
    )
    bootstrap_indices = np.load(
        FREEZE_DIR / "whole-target-bootstrap.npy",
        allow_pickle=False,
    )
    null_macro = null_macro_pds(
        forward_similarity,
        reverse_similarity,
        altair_zero,
        myllia_zero,
        derangements,
    )
    bootstrap_macro = bootstrap_macro_pds(
        forward_pds,
        reverse_pds,
        bootstrap_indices,
    )

    forward_mean = float(np.mean(forward_pds))
    reverse_mean = float(np.mean(reverse_pds))
    macro = 0.5 * (forward_mean + reverse_mean)
    null_p95 = float(np.quantile(null_macro, 0.95, method="linear"))
    bootstrap_lower = float(np.quantile(bootstrap_macro, 0.025, method="linear"))
    bootstrap_upper = float(np.quantile(bootstrap_macro, 0.975, method="linear"))
    forward_credit_mean = float(np.mean(forward_credit))
    reverse_credit_mean = float(np.mean(reverse_credit))
    gates = {
        "altair_to_myllia_pds_delta_positive": forward_mean > 0.5,
        "myllia_to_altair_pds_delta_positive": reverse_mean > 0.5,
        "macro_pds_at_least_0p60": macro >= 0.60,
        "macro_pds_above_null_p95": macro > null_p95,
        "bootstrap_lower_95_above_0p50": bootstrap_lower > 0.50,
        "altair_to_myllia_sign_credit_positive": forward_credit_mean > 0.0,
        "myllia_to_altair_sign_credit_positive": reverse_credit_mean > 0.0,
    }
    go = all(gates.values())
    response_genes = tuple(
        gene for gene, selected in zip(genes, response_mask, strict=True) if selected
    )
    result = {
        "schema_version": 1,
        "experiment": inputs["protocol"]["experiment"],
        "status": "complete_source_only_evaluation",
        "target_axis": {
            "count": len(targets),
            "labels": list(targets),
            "sha256": sha256_lines(targets),
        },
        "response_axis": {
            "count": len(response_genes),
            "sha256": sha256_lines(response_genes),
            "global_complete_case": True,
            "excluded_vcc_panel_count": len(panel),
            "excluded_common_target_count": len(targets),
            "excluded_union_count": len(set(panel) | set(targets)),
        },
        "zero_rows": {
            "altair": int(np.count_nonzero(altair_zero)),
            "myllia": int(np.count_nonzero(myllia_zero)),
        },
        "metrics": {
            "altair_to_myllia": {
                "pds": forward_mean,
                "pds_delta_from_chance": forward_mean - 0.5,
                "own_target_cosine": float(np.mean(own_cosine)),
                "top20_predictor_magnitude_sign_precision": float(
                    np.mean(forward_precision)
                ),
                "directional_sign_credit": forward_credit_mean,
            },
            "myllia_to_altair": {
                "pds": reverse_mean,
                "pds_delta_from_chance": reverse_mean - 0.5,
                "own_target_cosine": float(np.mean(own_cosine)),
                "top20_predictor_magnitude_sign_precision": float(
                    np.mean(reverse_precision)
                ),
                "directional_sign_credit": reverse_credit_mean,
            },
            "macro_pds": macro,
            "top20_selected_coordinates_per_target": selected_count,
        },
        "derangement_null": {
            "replicates": len(null_macro),
            "macro_pds_mean": float(np.mean(null_macro)),
            "macro_pds_p95": null_p95,
            "macro_pds_maximum": float(np.max(null_macro)),
            "one_sided_p_value": float(
                (1 + np.count_nonzero(null_macro >= macro)) / (len(null_macro) + 1)
            ),
        },
        "whole_target_bootstrap": {
            "replicates": len(bootstrap_macro),
            "macro_pds_lower_95": bootstrap_lower,
            "macro_pds_upper_95": bootstrap_upper,
        },
        "go_gate": {
            "criteria": gates,
            "go": go,
            "interpretation": (
                "justifies_a_later_held_context_donor_experiment"
                if go
                else "does_not_justify_a_held_context_donor_experiment"
            ),
        },
        "safety": {
            "held_truth_opened": False,
            "existing_outcome_directories_opened": False,
            "cell_eval2": False,
            "counts_generated": False,
            "submission_paths_opened": False,
        },
    }

    output.parent.mkdir(parents=True, exist_ok=True)
    temporary = Path(
        tempfile.mkdtemp(prefix=f".{output.name}.", dir=output.parent)
    )
    try:
        (temporary / "result.json").write_bytes(json_bytes(result))
        write_per_target(
            temporary / "per-target.tsv",
            targets,
            own_cosine,
            forward_pds,
            reverse_pds,
            forward_precision,
            reverse_precision,
            forward_credit,
            reverse_credit,
        )
        (temporary / "response-genes.txt").write_text(
            "".join(f"{gene}\n" for gene in response_genes),
            encoding="utf-8",
        )
        np.save(temporary / "null-macro-pds.npy", null_macro)
        np.save(temporary / "bootstrap-macro-pds.npy", bootstrap_macro)
        primary_files = (
            "bootstrap-macro-pds.npy",
            "null-macro-pds.npy",
            "per-target.tsv",
            "response-genes.txt",
            "result.json",
        )
        artifact_manifest = {
            "schema_version": 1,
            "experiment": result["experiment"],
            "artifact": "immutable_source_only_cross_source_evaluation",
            "files": {
                name: {
                    "sha256": sha256_file(temporary / name),
                    "size_bytes": (temporary / name).stat().st_size,
                }
                for name in primary_files
            },
            "bindings": {
                "authorization_manifest_sha256": sha256_file(AUTHORIZATION_PATH),
                "evaluator_sha256": sha256_file(Path(__file__).resolve()),
                "freeze_manifest_sha256": sha256_file(
                    FREEZE_DIR / "freeze-manifest.json"
                ),
                "preregistered_protocol_sha256": sha256_file(PROTOCOL_PATH),
            },
            "immutable_by_default": True,
        }
        (temporary / "artifact-manifest.json").write_bytes(
            json_bytes(artifact_manifest)
        )
        receipt = {
            "schema_version": 1,
            "experiment": result["experiment"],
            "status": "verified_source_only_evaluation",
            "artifact_manifest_sha256": sha256_file(
                temporary / "artifact-manifest.json"
            ),
            "result_sha256": sha256_file(temporary / "result.json"),
            "source_artifact_and_receipt_bindings": authorization["source_bindings"],
            "authorization": authorization["authorization"],
        }
        (temporary / "evaluation-receipt.json").write_bytes(json_bytes(receipt))
        (temporary / "evaluation-receipt.sha256").write_text(
            f"{sha256_file(temporary / 'evaluation-receipt.json')}\n",
            encoding="utf-8",
        )
        os.replace(temporary, output)
    except BaseException:
        shutil.rmtree(temporary, ignore_errors=True)
        raise
    return output


def verify_output(output: Path = DEFAULT_OUTPUT) -> dict[str, Any]:
    output = safe_output(output)
    require(output.is_dir() and not output.is_symlink(), "evaluation output missing")
    manifest = read_json(output / "artifact-manifest.json")
    receipt = read_json(output / "evaluation-receipt.json")
    require(
        sha256_file(output / "artifact-manifest.json")
        == receipt["artifact_manifest_sha256"],
        "artifact manifest hash changed",
    )
    require(
        sha256_file(output / "result.json") == receipt["result_sha256"],
        "result hash changed",
    )
    for name, record in manifest["files"].items():
        path = output / name
        require(path.stat().st_size == record["size_bytes"], f"{name} size changed")
        require(sha256_file(path) == record["sha256"], f"{name} hash changed")
    require(
        (output / "evaluation-receipt.sha256").read_text(encoding="utf-8").strip()
        == sha256_file(output / "evaluation-receipt.json"),
        "receipt hash sidecar changed",
    )
    return read_json(output / "result.json")


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("command", choices=("run", "verify"))
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    return parser.parse_args(argv)


def main(argv: Sequence[str] | None = None) -> None:
    args = parse_args(argv)
    if args.command == "run":
        output = run(args.output)
        print(output)
    else:
        result = verify_output(args.output)
        print(json.dumps(result["go_gate"], indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
