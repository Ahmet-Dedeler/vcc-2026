#!/usr/bin/env python3
"""Audit exact production integration of Experiment 122's panel residual."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import os
from pathlib import Path
import shutil
from typing import Any

import h5py
import numpy as np


ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
PROTOCOL_PATH = HERE / "preregistered_protocol.json"
FREEZE_PATH = HERE / "IMPLEMENTATION_FREEZE.json"
TS = 50_000.0


class AuditError(RuntimeError):
    pass


def require(condition: bool, message: str) -> None:
    if not condition:
        raise AuditError(message)


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while chunk := handle.read(1 << 20):
            digest.update(chunk)
    return digest.hexdigest()


def read_axis(path: Path, column: str) -> tuple[str, ...]:
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        require(reader.fieldnames == [column], f"{path}: header changed")
        values = tuple(str(row[column]) for row in reader)
    require(len(values) == len(set(values)), f"{path}: duplicate labels")
    return values


def verify_inputs(protocol: dict[str, Any]) -> None:
    items = [
        protocol["inputs"]["panel_residual"],
        protocol["inputs"]["staged_a15_effects"],
        protocol["inputs"]["gene_names"],
        protocol["inputs"]["panel_targets"],
        *protocol["inputs"]["released_controls"].values(),
    ]
    for item in items:
        path = ROOT / item["path"]
        require(path.is_file(), f"missing input: {path}")
        require(sha256_file(path) == item["sha256"], f"input hash changed: {path}")


def verify_freeze() -> None:
    require(FREEZE_PATH.is_file(), "implementation freeze is absent")
    freeze = json.loads(FREEZE_PATH.read_text(encoding="utf-8"))
    for relative, expected in freeze["files"].items():
        path = ROOT / relative
        require(sha256_file(path) == expected, f"implementation hash changed: {path}")


def decode_strings(values: np.ndarray) -> tuple[str, ...]:
    return tuple(
        value.decode("utf-8") if isinstance(value, (bytes, np.bytes_)) else str(value)
        for value in values
    )


def read_h5ad_index(node: h5py.Dataset | h5py.Group) -> tuple[str, ...]:
    if isinstance(node, h5py.Dataset):
        return decode_strings(np.asarray(node))
    require(
        str(node.attrs.get("encoding-type", "")) == "nullable-string-array",
        "unsupported H5AD index encoding",
    )
    require(set(node.keys()) == {"mask", "values"}, "nullable string fields changed")
    mask = np.asarray(node["mask"], dtype=bool)
    require(not bool(np.any(mask)), "H5AD gene index contains missing values")
    return decode_strings(np.asarray(node["values"]))


def control_profile(path: Path, expected_genes: tuple[str, ...]) -> np.ndarray:
    """Read a CSR h5ad in bounded chunks and return log1p(TS * pseudobulk)."""
    with h5py.File(path, "r") as handle:
        require("X" in handle, f"{path}: X missing")
        matrix = handle["X"]
        require(
            str(matrix.attrs.get("encoding-type", "")) == "csr_matrix",
            f"{path}: X is not CSR",
        )
        shape = tuple(map(int, matrix.attrs["shape"]))
        require(shape[1] == len(expected_genes), f"{path}: gene dimension changed")
        var_names = read_h5ad_index(handle["var"]["_index"])
        require(var_names == expected_genes, f"{path}: gene order changed")
        data = matrix["data"]
        indices = matrix["indices"]
        totals = np.zeros(shape[1], dtype=np.float64)
        chunk = 8_000_000
        for start in range(0, len(data), chunk):
            stop = min(start + chunk, len(data))
            values = np.asarray(data[start:stop], dtype=np.float64)
            columns = np.asarray(indices[start:stop], dtype=np.int64)
            np.add.at(totals, columns, values)
    require(bool(np.all(np.isfinite(totals))), f"{path}: nonfinite totals")
    require(float(totals.sum()) > 0.0, f"{path}: zero total")
    return np.log1p(TS * totals / totals.sum())


def profile_from_effect(
    effect: np.ndarray,
    b_ctrl: np.ndarray,
    *,
    min_cpm: float,
) -> np.ndarray:
    expression = np.expm1(b_ctrl)
    cpm = 1e6 * expression / max(float(expression.sum()), 1e-12)
    gate = (cpm > min_cpm) & (b_ctrl > 0)
    gated = np.where(gate, effect, 0.0)
    profile = np.clip(np.expm1(b_ctrl + gated), 0.0, None) / TS
    total = float(profile.sum())
    require(total > 0.0, "degenerate profile")
    return profile / total


def repair_to_sum(
    values: np.ndarray,
    target_sum: float,
    mutable: np.ndarray | None = None,
) -> None:
    if mutable is None:
        indices = np.arange(values.size)
    else:
        mutable = np.asarray(mutable, dtype=bool)
        require(mutable.shape == values.shape, "sum-repair mask shape changed")
        indices = np.flatnonzero(mutable)
    require(indices.size > 0, "no mutable coordinates for sum repair")
    order = indices[np.argsort(-values[indices], kind="stable")]
    for index_value in order:
        index = int(index_value)
        for _ in range(6):
            residual = math.fsum(map(float, values)) - target_sum
            if residual == 0.0:
                require(bool(np.all(values >= 0.0)), "sum repair made a negative value")
                return
            before = float(values[index])
            after = before - residual
            if after < 0.0 or after == before:
                break
            values[index] = after
    raise AuditError("could not enforce exact nonnegative compensated sum")


def panel_simplex_transport(
    baseline_profile: np.ndarray,
    b_ctrl: np.ndarray,
    residual_row: np.ndarray,
    panel_columns: np.ndarray,
    target_index: int,
    eta: float,
    *,
    min_cpm: float,
) -> tuple[np.ndarray, dict[str, Any]]:
    baseline = np.asarray(baseline_profile, dtype=np.float64)
    candidate = baseline.copy()
    expression = np.expm1(b_ctrl)
    cpm = 1e6 * expression / max(float(expression.sum()), 1e-12)
    panel_gate = (cpm[panel_columns] > min_cpm) & (b_ctrl[panel_columns] > 0)
    mutable = panel_gate & np.isfinite(residual_row)
    mutable[target_index] = False
    mutable &= baseline[panel_columns] > 0.0
    mutable_columns = panel_columns[mutable]
    if mutable_columns.size:
        original = baseline[mutable_columns]
        target_sum = math.fsum(map(float, original))
        logits = eta * np.asarray(residual_row[mutable], dtype=np.float64)
        logits -= float(np.max(logits))
        transported = original * np.exp(logits)
        transport_sum = math.fsum(map(float, transported))
        require(transport_sum > 0.0, "transport has zero mass")
        transported *= target_sum / transport_sum
        repair_to_sum(transported, target_sum)
        candidate[mutable_columns] = transported

    panel_mask = np.zeros(baseline.size, dtype=bool)
    panel_mask[panel_columns] = True
    nonpanel = ~panel_mask
    base_panel_sum = math.fsum(map(float, baseline[panel_mask]))
    if mutable_columns.size:
        panel_candidate = candidate[panel_columns].copy()
        repair_to_sum(panel_candidate, base_panel_sum, mutable)
        candidate[panel_columns] = panel_candidate
    candidate_panel_sum = math.fsum(map(float, candidate[panel_mask]))
    base_total = math.fsum((base_panel_sum, math.fsum(map(float, baseline[nonpanel]))))
    candidate_total = math.fsum(
        (candidate_panel_sum, math.fsum(map(float, candidate[nonpanel])))
    )
    base_effect = np.log1p(TS * baseline) - b_ctrl
    candidate_effect = np.log1p(TS * candidate) - b_ctrl
    invariants = {
        "finite": bool(np.all(np.isfinite(candidate))),
        "nonnegative": bool(np.all(candidate >= 0.0)),
        "exact_nonpanel_relative_abundance_identity": bool(
            np.array_equal(candidate[nonpanel], baseline[nonpanel])
        ),
        "exact_target_self_relative_abundance_identity": bool(
            candidate[panel_columns[target_index]]
            == baseline[panel_columns[target_index]]
        ),
        "exact_panel_mass_conservation": bool(candidate_panel_sum == base_panel_sum),
        "exact_total_mass_conservation": bool(candidate_total == base_total),
        "exact_nonpanel_log_effect_identity": bool(
            np.array_equal(candidate_effect[nonpanel], base_effect[nonpanel])
        ),
    }
    return candidate, {
        "mutable_coordinates": int(mutable_columns.size),
        "invariants": invariants,
        "residual_sign_retention": (
            float(
                np.mean(
                    np.sign(candidate[mutable_columns] - baseline[mutable_columns])
                    == np.sign(residual_row[mutable])
                )
            )
            if mutable_columns.size
            else 1.0
        ),
        "maximum_absolute_panel_change": float(
            np.max(np.abs(candidate[panel_columns] - baseline[panel_columns]))
        ),
    }


def emitted_nonpanel_totals_identical(
    baseline: np.ndarray,
    candidate: np.ndarray,
    nonpanel: np.ndarray,
    libraries: tuple[int, ...],
    cells: int,
) -> bool:
    for library in libraries:
        base_totals = np.floor(cells * library * baseline[nonpanel]).astype(np.int64)
        candidate_totals = np.floor(cells * library * candidate[nonpanel]).astype(
            np.int64
        )
        if not np.array_equal(base_totals, candidate_totals):
            return False
    return True


def cosine_rows(left: np.ndarray, right: np.ndarray) -> np.ndarray:
    numerator = np.sum(left * right, axis=1)
    denominator = np.linalg.norm(left, axis=1) * np.linalg.norm(right, axis=1)
    return np.divide(
        numerator,
        denominator,
        out=np.ones_like(numerator),
        where=denominator > 0,
    )


def summarize(values: list[float]) -> dict[str, float]:
    array = np.asarray(values, dtype=np.float64)
    return {
        "minimum": float(np.min(array)),
        "median": float(np.median(array)),
        "maximum": float(np.max(array)),
    }


def run(output: Path) -> dict[str, Any]:
    protocol = json.loads(PROTOCOL_PATH.read_text(encoding="utf-8"))
    verify_inputs(protocol)
    verify_freeze()
    genes = read_axis(ROOT / protocol["inputs"]["gene_names"]["path"], "gene_name")
    panel = read_axis(ROOT / protocol["inputs"]["panel_targets"]["path"], "target_gene")
    require(len(genes) == 18_533, "official gene count changed")
    require(len(panel) == 300, "panel target count changed")
    gene_position = {gene: index for index, gene in enumerate(genes)}
    panel_columns = np.asarray([gene_position[gene] for gene in panel], dtype=np.int64)
    panel_mask = np.zeros(len(genes), dtype=bool)
    panel_mask[panel_columns] = True
    nonpanel = ~panel_mask

    effects = np.load(
        ROOT / protocol["inputs"]["staged_a15_effects"]["path"],
        allow_pickle=False,
    ).astype(np.float64)
    residual = np.load(
        ROOT / protocol["inputs"]["panel_residual"]["path"],
        allow_pickle=False,
    ).astype(np.float64)
    require(effects.shape == (300, 18_533), "a15 effect shape changed")
    require(residual.shape == (300, 300), "panel residual shape changed")
    require(bool(np.all(np.diag(residual) == 0.0)), "panel residual self changed")

    row_sums = residual.sum(axis=1)
    context_results: dict[str, Any] = {}
    all_transport_pass = True
    for context, item in protocol["inputs"]["released_controls"].items():
        b_ctrl = control_profile(ROOT / item["path"], genes)
        baseline_profiles = np.stack(
            [
                profile_from_effect(row, b_ctrl, min_cpm=10.0)
                for row in effects
            ]
        )
        baseline_effects = np.log1p(TS * baseline_profiles) - b_ctrl[None, :]
        dose_results: dict[str, Any] = {}
        for eta in map(float, protocol["fixed_doses"]):
            expanded = np.zeros_like(effects)
            expanded[:, panel_columns] = eta * residual
            naive_profiles = np.stack(
                [
                    profile_from_effect(row, b_ctrl, min_cpm=10.0)
                    for row in effects + expanded
                ]
            )
            naive_effects = np.log1p(TS * naive_profiles) - b_ctrl[None, :]
            naive_delta = naive_profiles[:, nonpanel] - baseline_profiles[:, nonpanel]
            naive_cosine = cosine_rows(
                baseline_effects[:, nonpanel],
                naive_effects[:, nonpanel],
            )

            transport_profiles = np.empty_like(baseline_profiles)
            transport_rows: list[dict[str, Any]] = []
            for target_index in range(len(panel)):
                transported, audit = panel_simplex_transport(
                    baseline_profiles[target_index],
                    b_ctrl,
                    residual[target_index],
                    panel_columns,
                    target_index,
                    eta,
                    min_cpm=10.0,
                )
                audit["exact_nonpanel_emitted_pseudobulk_total_identity_at_all_representative_libraries"] = emitted_nonpanel_totals_identical(
                    baseline_profiles[target_index],
                    transported,
                    nonpanel,
                    tuple(protocol["production_semantics"]["representative_median_libraries"]),
                    int(protocol["production_semantics"]["cells_per_perturbation"]),
                )
                audit["all_invariants_pass"] = bool(
                    all(audit["invariants"].values())
                    and audit[
                        "exact_nonpanel_emitted_pseudobulk_total_identity_at_all_representative_libraries"
                    ]
                )
                all_transport_pass &= audit["all_invariants_pass"]
                transport_profiles[target_index] = transported
                transport_rows.append(audit)
            transport_effects = np.log1p(TS * transport_profiles) - b_ctrl[None, :]
            transport_cosine = cosine_rows(
                baseline_effects[:, nonpanel],
                transport_effects[:, nonpanel],
            )
            dose_results[f"{eta:.2f}"] = {
                "naive_log_effect": {
                    "rows_with_exact_nonpanel_profile_identity": int(
                        np.sum(np.all(naive_delta == 0.0, axis=1))
                    ),
                    "maximum_absolute_nonpanel_profile_change": float(
                        np.max(np.abs(naive_delta))
                    ),
                    "nonpanel_effect_cosine": summarize(naive_cosine.tolist()),
                },
                "panel_simplex": {
                    "all_rows_pass": bool(
                        all(row["all_invariants_pass"] for row in transport_rows)
                    ),
                    "minimum_mutable_coordinates": int(
                        min(row["mutable_coordinates"] for row in transport_rows)
                    ),
                    "residual_sign_retention": summarize(
                        [row["residual_sign_retention"] for row in transport_rows]
                    ),
                    "maximum_absolute_panel_change": float(
                        max(
                            row["maximum_absolute_panel_change"]
                            for row in transport_rows
                        )
                    ),
                    "nonpanel_effect_cosine": summarize(transport_cosine.tolist()),
                },
            }
        context_results[context] = dose_results

    result = {
        "schema_version": 1,
        "experiment": "124-panel-simplex-transport",
        "status": "integration_audit_complete",
        "panel_residual_row_sums": {
            "exact_zero_rows": int(np.sum(row_sums == 0.0)),
            "nonzero_rows": int(np.sum(row_sums != 0.0)),
            "absolute": summarize(np.abs(row_sums).tolist()),
        },
        "contexts": context_results,
        "decision": {
            "naive_log_effect_exact_pds_null": bool(
                all(
                    arm["naive_log_effect"][
                        "rows_with_exact_nonpanel_profile_identity"
                    ]
                    == 300
                    for context in context_results.values()
                    for arm in context.values()
                )
            ),
            "panel_simplex_all_invariants_pass": bool(all_transport_pass),
            "production_promotion_authorized": False,
            "independent_confirmation_required": True,
        },
        "forbidden_actions_respected": {
            "challenge_truth_accessed": False,
            "held_perturbation_outcome_accessed": False,
            "cell_eval2_run": False,
            "count_or_submission_built": False,
            "submission_sent": False,
        },
    }
    require(
        not result["decision"]["naive_log_effect_exact_pds_null"],
        "naive diagnostic unexpectedly remained exactly PDS-null",
    )
    require(
        result["decision"]["panel_simplex_all_invariants_pass"],
        "panel-simplex invariant failed",
    )

    if output.exists():
        raise FileExistsError(f"refusing to overwrite {output}")
    stage = output.parent / f".{output.name}.tmp-{os.getpid()}"
    if stage.exists():
        shutil.rmtree(stage)
    stage.mkdir(parents=True)
    result_path = stage / "results.json"
    result_path.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    manifest = {
        "schema_version": 1,
        "experiment": "124-panel-simplex-transport",
        "files": {
            "results.json": {
                "bytes": result_path.stat().st_size,
                "sha256": sha256_file(result_path),
            }
        },
        "protocol_sha256": sha256_file(PROTOCOL_PATH),
        "implementation_freeze_sha256": sha256_file(FREEZE_PATH),
    }
    (stage / "artifact-manifest.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n"
    )
    os.replace(stage, output)
    for path in output.iterdir():
        path.chmod(0o444)
    output.chmod(0o555)
    return result


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--audit", action="store_true")
    parser.add_argument(
        "--output",
        type=Path,
        default=HERE / "work/integration-audit-v1",
    )
    args = parser.parse_args()
    require(args.audit, "run requires --audit")
    result = run(args.output.resolve())
    print(json.dumps(result["decision"], indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
