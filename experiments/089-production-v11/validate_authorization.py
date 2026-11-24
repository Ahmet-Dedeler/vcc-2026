#!/usr/bin/env python3
"""Validate the exact terminal experiment-065 evidence for experiment 071.

Passing this validator authorizes exploratory packaging only. It never
authorizes a Virtual Cell Challenge submission.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import os
import re
from pathlib import Path, PurePosixPath
from typing import Any, Mapping

EXPERIMENT = "065-pm-cd4-hct-count"
S3_EXPERIMENT = "exp065-pm-cd4-hct-count"
EXP061_AVG6 = 0.12173082252852911
RUN_ID_RE = re.compile(r"^[0-9]{8}T[0-9]{6}Z$")
SHA256_RE = re.compile(r"^[0-9a-f]{64}$")
SHA256_LINE_RE = re.compile(r"^([0-9a-f]{64}) [ *](.+)$")

METRICS = (
    "pds_cosine",
    "expr_mse_unbiased_capped_norm",
    "de_wilcoxon_direction_fidelity_yield_raw",
    "de_wilcoxon_direction_reach_raw",
    "de_wilcoxon_sig_jaccard",
    "de_wilcoxon_lfc_nmae",
)
LOWER_IS_BETTER = {
    "expr_mse_unbiased_capped_norm",
    "de_wilcoxon_lfc_nmae",
}
ANCHORS = {
    "pds_cosine": (0.500, 0.9525),
    "expr_mse_unbiased_capped_norm": (0.989, 0.0365),
    "de_wilcoxon_lfc_nmae": (1.0013, 0.400),
    "de_wilcoxon_direction_fidelity_yield_raw": (0.5135, 0.8135),
    "de_wilcoxon_direction_reach_raw": (0.072, 0.968),
    "de_wilcoxon_sig_jaccard": (0.029, 0.399),
}

REQUIRED_TERMINAL_ARTIFACTS = (
    "exit-code.txt",
    "artifact-sha256.txt",
    "results.json",
    "build-report.json",
    "baseline-agg_results.csv",
    "candidate-run/agg_results.csv",
    "candidate.h5ad.sha256",
    "score-memory-guard.json",
    "runtime.json",
)
REQUIRED_INTERNAL_MANIFEST_ENTRIES = (
    "results.json",
    "build-report.json",
    "baseline-agg_results.csv",
    "candidate-run/agg_results.csv",
    "candidate.h5ad",
    "candidate.h5ad.sha256",
    "score-memory-guard.json",
    "runtime.json",
)

EXPECTED_PROTOCOL = {
    "held_out_truth": "Orion HCT116",
    "candidate_independent_operator": True,
    "operator_fit_sources": ["k562_gwps", "orion_hek293t"],
    "all_300_vcc_target_coordinates_excluded_from_operator_and_cd4": True,
    "gamma": 0.65,
    "alpha": 1.08,
    "cpm_gate": 5.0,
    "emission": "quasi exact-total",
    "spread": 0.25,
    "seed": 0,
    "cd4_conditions": ["cd4_rest", "cd4_stim8hr", "cd4_stim48hr"],
    "cd4_mode": "gene_norm",
    "cd4_fraction": 0.05,
    "submission": False,
}
EXPECTED_BUILD_INPUT_HASHES = {
    "real": "d20a1e9da59debca8db7d83b955d12c3e1384b530ae6154191d95039b1655863",
    "k562_delta": ("0c9422280fed7afb81ca474e8b9024af788c9ffbde32e11f96cb52422e296aa9"),
    "hek_delta": ("d95202d984d9322dc3830f23b6ef55001ba8562dea1d72c43fb0c9d5064c7d1b"),
    "perturbmatch_delta": (
        "efa20944f6f6485a7cc7f484c7ab70bd2ba491a6e07b231f306e104b2e1e4cf2"
    ),
    "whitener": ("a46264ebb6cb4b877032a6e4c68e76069eb9112c6f1bb7f15e1133fae04a8a77"),
    "whitener_metadata": (
        "cf5c1de31d03875b15f6fd11d3f628aadfccf04e3346493c831c3a6f441b4793"
    ),
    "panel_csv": ("f57edd7b912ebd718efc7ee9d0f334772513e7cc418d133ce525470e373b3276"),
    "genes_csv": ("25bfa66715e186bebabce7ac788bbcea47e2bf59ca70be1f8f3a06f2f0e47201"),
    "exp061_build_report": (
        "be855dd8174ba055616260e5c31fe9e2db0ea0b933f0e336fc25bc4cef367c47"
    ),
    "exp061_source": (
        "77e6fbf7edf1aeea6d39be4fe801591d0b0c2867f95042b697f0fc4042a2fae4"
    ),
    "exp062_source": (
        "16417a9b985d9aad8993bf2f24ce550ea1e6871d9401bb4672d11ec6d0046176"
    ),
}
EXPECTED_CD4_HASHES = {
    "cd4_rest": {
        "delta.npy": (
            "b2bafd6e24c2aa34b03ec17836db34fcd9512f78920aed207413a9d727f9e765"
        ),
        "b_ctrl.npy": (
            "5c47ef6b15a00a867b7b24034273bf3aba118d7cdf87ce58ec58427751313d81"
        ),
        "genes.txt": (
            "cbada1fd0eaeeb0385fa9c905c08ea90f1d472f5f10e832c47836ed9d7fe45fd"
        ),
        "perts.txt": (
            "ddcd8ff07e2f38d93fc035ce2aadebe4bcaeb097ef95e1f26e8e4860b641374b"
        ),
        "meta.json": (
            "cece9564be470b133ac6399012af42bcfa1934776c1820838c6aa8ab6a15b9cc"
        ),
    },
    "cd4_stim8hr": {
        "delta.npy": (
            "eea7996d6694450fb1167b9ebbb6b709a3477e6d0347b49a4475e13ff7800c99"
        ),
        "b_ctrl.npy": (
            "75e35371ec33677462acdbc56938c2cc266d583dc2e75ed50b6b825c4a20161e"
        ),
        "genes.txt": (
            "cbada1fd0eaeeb0385fa9c905c08ea90f1d472f5f10e832c47836ed9d7fe45fd"
        ),
        "perts.txt": (
            "639feb3fba1d6cf3b5bb5e2ef9dfdc0ff3044adb4c7477b4588c82ff096d828e"
        ),
        "meta.json": (
            "9ca1e9b8ccbf785d88b2a31bc47818479327c205145faa78d9b8eaa215bff4cc"
        ),
    },
    "cd4_stim48hr": {
        "delta.npy": (
            "568119ac25c8d16116d8953d3ee61762fe3c45d9ee8c22a6e7b2669c39837a9a"
        ),
        "b_ctrl.npy": (
            "a35f6c639937332cc483c82b9717dec9ff8a3fc70e6030e1b9f6df5fcc6658fd"
        ),
        "genes.txt": (
            "cbada1fd0eaeeb0385fa9c905c08ea90f1d472f5f10e832c47836ed9d7fe45fd"
        ),
        "perts.txt": (
            "4bf6468927fc11f0709bb0e63d14852f5a49f4f577455f0fa92e88719777d618"
        ),
        "meta.json": (
            "62390f52782a592b2a144349bebebffc2675309a614bff0eef99c0b315668ef8"
        ),
    },
}
EXPECTED_EFFECT_HASHES = {
    "baseline_pre_sha256": (
        "05a26dbe3a1e0f1d2e4f420dc3431c8f1d76291f47a8c940f9389b999fe26cac"
    ),
    "candidate_pre_sha256": (
        "5bdad760fa24bfb55250200afd5c553b24872cef3c295fc518739bd7be058793"
    ),
    "baseline_final_sha256": (
        "6a5b722548c3808d2d5adc5f91aa8da7a5e0a48470ac70fa4be85ee479740b8f"
    ),
    "candidate_final_sha256": (
        "a805ff771f55a2097a0bff93d0614fe9d6141c7af329501cad1c4072bcf5b83b"
    ),
    "joint_final_sha256": (
        "709fc341b8ca5ea2d75c7f21b6bfbcdcd76874a8eb49339d3aaba8a26a877b1c"
    ),
}


def file_sha256(path: Path, block_bytes: int = 8 << 20) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while block := handle.read(block_bytes):
            digest.update(block)
    return digest.hexdigest()


def load_json(path: Path) -> Any:
    with path.open() as handle:
        value = json.load(handle)
    reject_nonfinite(value, str(path))
    return value


def reject_nonfinite(value: Any, location: str) -> None:
    if isinstance(value, bool) or value is None or isinstance(value, str):
        return
    if isinstance(value, (int, float)):
        if not math.isfinite(float(value)):
            raise ValueError(f"{location}: non-finite numeric value")
        return
    if isinstance(value, list):
        for index, item in enumerate(value):
            reject_nonfinite(item, f"{location}[{index}]")
        return
    if isinstance(value, dict):
        for key, item in value.items():
            reject_nonfinite(item, f"{location}.{key}")


def require_sha256(value: Any, location: str) -> str:
    if not isinstance(value, str) or SHA256_RE.fullmatch(value) is None:
        raise ValueError(f"{location}: invalid SHA-256")
    return value


def require_regular_file(root: Path, relative: str) -> Path:
    pure = PurePosixPath(relative)
    if pure.is_absolute() or ".." in pure.parts or "." in pure.parts:
        raise ValueError(f"unsafe artifact path: {relative!r}")
    path = root.joinpath(*pure.parts)
    if path.is_symlink() or not path.is_file():
        raise FileNotFoundError(f"missing regular artifact: {relative}")
    if not path.resolve().is_relative_to(root.resolve()):
        raise ValueError(f"artifact escapes root: {relative!r}")
    return path


def normalize_exp065_output_path(raw: str) -> str:
    if "\x00" in raw or "\\" in raw:
        raise ValueError(f"unsafe manifest path: {raw!r}")
    path = PurePosixPath(raw)
    if path.is_absolute():
        base = PurePosixPath("/work/exp065/output")
        try:
            path = path.relative_to(base)
        except ValueError as error:
            raise ValueError(f"unexpected absolute exp065 path: {raw!r}") from error
    elif path.parts[:1] == (".",):
        path = PurePosixPath(*path.parts[1:])
    if (
        path.is_absolute()
        or not path.parts
        or any(part in ("", ".", "..") for part in path.parts)
    ):
        raise ValueError(f"unsafe exp065 output path: {raw!r}")
    return path.as_posix()


def parse_sha256_manifest(path: Path) -> dict[str, str]:
    entries: dict[str, str] = {}
    for line_number, line in enumerate(path.read_text().splitlines(), start=1):
        match = SHA256_LINE_RE.fullmatch(line)
        if match is None:
            raise ValueError(f"{path}:{line_number}: invalid sha256sum line")
        digest, raw_name = match.groups()
        name = normalize_exp065_output_path(raw_name)
        if name in entries:
            raise ValueError(f"{path}: duplicate manifest entry {name}")
        entries[name] = digest
    for required in REQUIRED_INTERNAL_MANIFEST_ENTRIES:
        if required not in entries:
            raise ValueError(f"{path}: missing manifest entry {required}")
    return entries


def parse_candidate_checksum(path: Path) -> str:
    lines = path.read_text().splitlines()
    if len(lines) != 1:
        raise ValueError(f"{path}: expected exactly one checksum line")
    match = SHA256_LINE_RE.fullmatch(lines[0])
    if match is None:
        raise ValueError(f"{path}: invalid candidate checksum")
    digest, raw_name = match.groups()
    if normalize_exp065_output_path(raw_name) != "candidate.h5ad":
        raise ValueError(f"{path}: checksum does not name candidate.h5ad")
    return digest


def validate_terminal_manifest(
    root: Path,
    manifest_path: Path,
    *,
    run_id: str,
    expected_prefix: str,
) -> tuple[dict[str, Any], dict[str, dict[str, Any]]]:
    manifest = load_json(manifest_path)
    if manifest.get("schema") != 1:
        raise ValueError("terminal S3 manifest schema must be 1")
    if manifest.get("experiment") != EXPERIMENT:
        raise ValueError("terminal S3 manifest experiment mismatch")
    if manifest.get("run_id") != run_id:
        raise ValueError("terminal S3 manifest run_id mismatch")
    s3 = manifest.get("s3")
    if not isinstance(s3, dict) or s3.get("prefix") != expected_prefix:
        raise ValueError("terminal S3 manifest prefix mismatch")
    if not isinstance(s3.get("bucket"), str) or not s3["bucket"]:
        raise ValueError("terminal S3 manifest bucket is missing")

    raw_entries = manifest.get("artifacts")
    if not isinstance(raw_entries, list):
        raise ValueError("terminal S3 manifest artifacts must be a list")
    entries: dict[str, dict[str, Any]] = {}
    keys: set[str] = set()
    for entry in raw_entries:
        if not isinstance(entry, dict):
            raise ValueError("terminal S3 manifest artifact must be an object")
        relative = entry.get("path")
        if not isinstance(relative, str):
            raise ValueError("terminal S3 manifest artifact path is missing")
        path = require_regular_file(root, relative)
        if relative in entries:
            raise ValueError(f"duplicate terminal artifact path: {relative}")
        expected_key = f"{expected_prefix}/{relative}"
        if entry.get("key") != expected_key:
            raise ValueError(f"{relative}: terminal S3 key mismatch")
        if expected_key in keys:
            raise ValueError(f"duplicate terminal S3 key: {expected_key}")
        keys.add(expected_key)
        digest = require_sha256(entry.get("sha256"), f"{relative}.sha256")
        if file_sha256(path) != digest:
            raise ValueError(f"{relative}: staged artifact hash mismatch")
        if entry.get("size") != path.stat().st_size:
            raise ValueError(f"{relative}: staged artifact size mismatch")
        entries[relative] = entry
    if set(entries) != set(REQUIRED_TERMINAL_ARTIFACTS):
        raise ValueError(
            "terminal artifact set mismatch: "
            f"missing={sorted(set(REQUIRED_TERMINAL_ARTIFACTS) - set(entries))}, "
            f"unexpected={sorted(set(entries) - set(REQUIRED_TERMINAL_ARTIFACTS))}"
        )
    return manifest, entries


def read_mean(path: Path) -> dict[str, float]:
    with path.open(newline="") as handle:
        rows = list(csv.DictReader(handle))
    mean_rows = [row for row in rows if row.get("statistic") == "mean"]
    if len(mean_rows) != 1:
        raise ValueError(f"{path}: expected exactly one mean row")
    values: dict[str, float] = {}
    for metric in METRICS:
        try:
            value = float(mean_rows[0][metric])
        except (KeyError, TypeError, ValueError) as error:
            raise ValueError(f"{path}: invalid {metric}") from error
        if not math.isfinite(value):
            raise ValueError(f"{path}: non-finite {metric}")
        values[metric] = value
    return values


def scaled(metric: str, value: float) -> float:
    baseline, replicate = ANCHORS[metric]
    result = (value - baseline) / (replicate - baseline)
    if metric == "expr_mse_unbiased_capped_norm":
        return max(0.0, min(1.0, result))
    return result


def require_equal(actual: Any, expected: Any, location: str) -> None:
    if actual != expected:
        raise ValueError(f"{location}: {actual!r} != {expected!r}")


def require_close(actual: Any, expected: float, location: str) -> None:
    if isinstance(actual, bool):
        raise ValueError(f"{location}: boolean is not a metric")
    try:
        value = float(actual)
    except (TypeError, ValueError) as error:
        raise ValueError(f"{location}: not numeric") from error
    if not math.isfinite(value):
        raise ValueError(f"{location}: non-finite")
    if not math.isclose(value, expected, rel_tol=1e-12, abs_tol=1e-12):
        raise ValueError(f"{location}: {value} != recomputed {expected}")


def validate_build_report(
    report: Mapping[str, Any],
    internal_manifest: Mapping[str, str],
    candidate_digest: str,
) -> None:
    require_equal(report.get("experiment"), EXPERIMENT, "build-report experiment")
    require_equal(report.get("status"), "candidate_built", "build-report status")
    protocol = report.get("protocol")
    if not isinstance(protocol, dict):
        raise ValueError("build-report protocol is missing")
    for key, expected in EXPECTED_PROTOCOL.items():
        require_equal(protocol.get(key), expected, f"build-report protocol.{key}")

    require_equal(
        report.get("inputs_sha256"),
        EXPECTED_BUILD_INPUT_HASHES,
        "build-report inputs_sha256",
    )
    require_equal(
        report.get("cd4_inputs_sha256"),
        EXPECTED_CD4_HASHES,
        "build-report cd4_inputs_sha256",
    )
    effect_audit = report.get("effect_audit")
    if not isinstance(effect_audit, dict):
        raise ValueError("build-report effect_audit is missing")
    for key, expected in EXPECTED_EFFECT_HASHES.items():
        require_equal(effect_audit.get(key), expected, f"effect_audit.{key}")
    require_equal(
        effect_audit.get("panel_cd4_residual_nonzero_entries"),
        0,
        "effect_audit.panel_cd4_residual_nonzero_entries",
    )
    require_equal(
        effect_audit.get("missing_pm_rows_equal_to_baseline_before_cd4"),
        True,
        "effect_audit.missing_pm_rows_equal_to_baseline_before_cd4",
    )

    axis = report.get("axis")
    if not isinstance(axis, dict):
        raise ValueError("build-report axis is missing")
    require_equal(axis.get("hct_targets"), 132, "build-report axis.hct_targets")
    require_equal(axis.get("official_genes"), 18533, "build-report axis.official_genes")
    require_equal(axis.get("operator_genes"), 17942, "build-report axis.operator_genes")

    cd4_audit = report.get("cd4_audit")
    if not isinstance(cd4_audit, dict):
        raise ValueError("build-report cd4_audit is missing")
    require_equal(cd4_audit.get("mode"), "gene_norm", "build-report cd4 mode")
    require_equal(
        cd4_audit.get("fully_supported_targets"),
        125,
        "build-report fully supported CD4 targets",
    )

    candidate = report.get("candidate_h5ad")
    if not isinstance(candidate, dict):
        raise ValueError("build-report candidate_h5ad is missing")
    require_equal(
        normalize_exp065_output_path(str(candidate.get("path"))),
        "candidate.h5ad",
        "build-report candidate path",
    )
    require_equal(candidate.get("sha256"), candidate_digest, "candidate build hash")
    require_equal(
        internal_manifest.get("candidate.h5ad"),
        candidate_digest,
        "candidate artifact-manifest hash",
    )
    if not isinstance(candidate.get("bytes"), int) or candidate["bytes"] <= 0:
        raise ValueError("build-report candidate byte count is invalid")


def validate_results(
    result: Mapping[str, Any],
    baseline_path: Path,
    candidate_path: Path,
) -> dict[str, float]:
    require_equal(result.get("experiment"), EXPERIMENT, "results experiment")
    require_equal(result.get("status"), "scored", "results status")
    require_equal(
        result.get("comparison"),
        "joint frozen candidate minus exact clean gamma-0.65 HCT baseline",
        "results comparison",
    )
    metrics = result.get("metrics")
    if not isinstance(metrics, dict) or set(metrics) != set(METRICS):
        raise ValueError("results metrics set mismatch")

    baseline = read_mean(baseline_path)
    candidate = read_mean(candidate_path)
    baseline_scaled_values: list[float] = []
    candidate_scaled_values: list[float] = []
    for metric in METRICS:
        block = metrics[metric]
        if not isinstance(block, dict):
            raise ValueError(f"results metric block is invalid: {metric}")
        baseline_scaled = scaled(metric, baseline[metric])
        candidate_scaled = scaled(metric, candidate[metric])
        raw_delta = candidate[metric] - baseline[metric]
        expected = {
            "preferred": "lower" if metric in LOWER_IS_BETTER else "higher",
            "baseline_raw": baseline[metric],
            "candidate_raw": candidate[metric],
            "raw_candidate_minus_baseline": raw_delta,
            "raw_favorable_delta": (
                -raw_delta if metric in LOWER_IS_BETTER else raw_delta
            ),
            "baseline_scaled": baseline_scaled,
            "candidate_scaled": candidate_scaled,
            "scaled_candidate_minus_baseline": candidate_scaled - baseline_scaled,
        }
        require_equal(
            block.get("preferred"), expected["preferred"], f"{metric}.preferred"
        )
        for key, value in expected.items():
            if key != "preferred":
                require_close(block.get(key), float(value), f"{metric}.{key}")
        baseline_scaled_values.append(baseline_scaled)
        candidate_scaled_values.append(candidate_scaled)

    baseline_avg6 = sum(baseline_scaled_values) / len(METRICS)
    candidate_avg6 = sum(candidate_scaled_values) / len(METRICS)
    delta = candidate_avg6 - baseline_avg6
    avg6 = result.get("avg6")
    if not isinstance(avg6, dict):
        raise ValueError("results avg6 is missing")
    require_close(avg6.get("baseline"), baseline_avg6, "avg6.baseline")
    require_close(avg6.get("candidate"), candidate_avg6, "avg6.candidate")
    require_close(avg6.get("candidate_minus_baseline"), delta, "avg6.delta")

    score_files = result.get("score_files")
    if not isinstance(score_files, dict):
        raise ValueError("results score_files is missing")
    expected_files = {
        "baseline": ("baseline-agg_results.csv", baseline_path),
        "candidate": ("candidate-run/agg_results.csv", candidate_path),
    }
    for role, (relative, path) in expected_files.items():
        item = score_files.get(role)
        if not isinstance(item, dict):
            raise ValueError(f"results score_files.{role} is missing")
        require_equal(
            normalize_exp065_output_path(str(item.get("path"))),
            relative,
            f"results score_files.{role}.path",
        )
        require_equal(
            item.get("sha256"),
            file_sha256(path),
            f"results score_files.{role}.sha256",
        )

    pds_delta = candidate_scaled_values[0] - baseline_scaled_values[0]
    if delta <= 0:
        raise ValueError(f"exploratory packaging gate failed: avg6 delta {delta} <= 0")
    if pds_delta < 0:
        raise ValueError(
            f"exploratory packaging gate failed: PDS scaled delta {pds_delta} < 0"
        )
    if candidate_avg6 <= EXP061_AVG6:
        raise ValueError(
            "exploratory packaging gate failed: "
            f"candidate avg6 {candidate_avg6} <= exp061 {EXP061_AVG6}"
        )
    return {
        "baseline_avg6": baseline_avg6,
        "candidate_avg6": candidate_avg6,
        "candidate_minus_baseline": delta,
        "pds_scaled_candidate_minus_baseline": pds_delta,
        "exp061_candidate_avg6": EXP061_AVG6,
    }


def validate_runtime(path: Path) -> None:
    runtime = load_json(path)
    require_equal(runtime.get("cell_eval2"), "0.16.0", "runtime cell_eval2")
    if not isinstance(runtime.get("python"), str) or not runtime["python"]:
        raise ValueError("runtime Python version is missing")


def validate_score_guard(path: Path) -> None:
    guard = load_json(path)
    require_equal(guard.get("return_code"), 0, "score guard return_code")
    require_equal(
        guard.get("killed_for_memory"), False, "score guard killed_for_memory"
    )
    maximum = guard.get("max_rss_gib")
    peak = guard.get("peak_process_tree_rss_gib")
    require_close(maximum, 55.0, "score guard max_rss_gib")
    if not isinstance(peak, (int, float)) or isinstance(peak, bool):
        raise ValueError("score guard peak RSS is invalid")
    if float(peak) < 0 or float(peak) > float(maximum):
        raise ValueError("score guard peak RSS lies outside the enforced limit")
    command = guard.get("command")
    if not isinstance(command, list):
        raise ValueError("score guard command is missing")
    joined = " ".join(map(str, command))
    for token in ("cell-eval2", "--preset", "vcc2026"):
        if token not in joined:
            raise ValueError(f"score guard command is missing {token}")


def validate_authorization(
    artifact_root: Path,
    terminal_manifest_path: Path,
    *,
    exp065_run_id: str,
    expected_prefix: str,
) -> dict[str, Any]:
    if RUN_ID_RE.fullmatch(exp065_run_id) is None:
        raise ValueError(f"invalid EXP065_RUN_ID: {exp065_run_id!r}")
    exact_prefix = f"output/{S3_EXPERIMENT}/{exp065_run_id}"
    if expected_prefix != exact_prefix:
        raise ValueError(f"expected prefix {expected_prefix!r} != {exact_prefix!r}")

    root = artifact_root.expanduser().resolve()
    terminal_manifest_path = terminal_manifest_path.expanduser().resolve()
    terminal_manifest, terminal_entries = validate_terminal_manifest(
        root,
        terminal_manifest_path,
        run_id=exp065_run_id,
        expected_prefix=expected_prefix,
    )
    exit_code = require_regular_file(root, "exit-code.txt").read_text().strip()
    if exit_code != "0":
        raise ValueError(f"experiment 065 exit code is {exit_code!r}, not '0'")

    internal_manifest_path = require_regular_file(root, "artifact-sha256.txt")
    internal_manifest = parse_sha256_manifest(internal_manifest_path)
    for relative in REQUIRED_INTERNAL_MANIFEST_ENTRIES:
        if relative == "candidate.h5ad":
            continue
        local = require_regular_file(root, relative)
        require_equal(
            internal_manifest[relative],
            file_sha256(local),
            f"artifact-sha256 {relative}",
        )

    candidate_digest = parse_candidate_checksum(
        require_regular_file(root, "candidate.h5ad.sha256")
    )
    require_equal(
        candidate_digest,
        internal_manifest["candidate.h5ad"],
        "candidate checksum vs artifact manifest",
    )

    build_report = load_json(require_regular_file(root, "build-report.json"))
    validate_build_report(build_report, internal_manifest, candidate_digest)
    validate_runtime(require_regular_file(root, "runtime.json"))
    validate_score_guard(require_regular_file(root, "score-memory-guard.json"))

    results_path = require_regular_file(root, "results.json")
    baseline_path = require_regular_file(root, "baseline-agg_results.csv")
    candidate_path = require_regular_file(root, "candidate-run/agg_results.csv")
    result = load_json(results_path)
    gate = validate_results(result, baseline_path, candidate_path)

    for relative, entry in terminal_entries.items():
        require_equal(
            entry["sha256"],
            file_sha256(require_regular_file(root, relative)),
            f"terminal manifest {relative}",
        )

    return {
        "schema": 1,
        "experiment": "071-production-pm-cd4",
        "status": "exploratory_packaging_gate_passed",
        "authorization_scope": "exploratory_packaging_only",
        "challenge_submission_authorized": False,
        "exp065": {
            "experiment": EXPERIMENT,
            "run_id": exp065_run_id,
            "s3_bucket": terminal_manifest["s3"]["bucket"],
            "s3_prefix": expected_prefix,
            "terminal_manifest_sha256": file_sha256(terminal_manifest_path),
            "artifact_manifest_sha256": file_sha256(internal_manifest_path),
            "candidate_h5ad_sha256": candidate_digest,
        },
        "packaging_gate": gate,
        "semantic_limitations": {
            "exact_v9_equivalence": False,
            "normalization": (
                "experiment 065 builds at alpha 1.08 and applies a frozen "
                "baseline-derived operator/global rescaling; exact v9 production "
                "uses alpha 1.0796 inside MultiDonorModel.effect and "
                "candidate-specific on-axis preservation"
            ),
            "random_streams": (
                "experiment 065 uses one RNG for control-library selection and "
                "count emission; exact v9 production uses separate deterministic "
                "RNG streams"
            ),
            "interpretation": (
                "a positive experiment-065 result permits construction and "
                "inspection of an exploratory package only; it does not establish "
                "an exact-v9 challenge-context gain and does not authorize submission"
            ),
        },
    }


def atomic_write_json(path: Path, value: Mapping[str, Any]) -> None:
    path = path.expanduser().resolve()
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.name}.partial-{os.getpid()}")
    try:
        temporary.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n")
        os.replace(temporary, path)
    except BaseException:
        temporary.unlink(missing_ok=True)
        raise


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--artifact-root", type=Path, required=True)
    parser.add_argument("--terminal-manifest", type=Path, required=True)
    parser.add_argument("--exp065-run-id", required=True)
    parser.add_argument("--expected-prefix", required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    authorization = validate_authorization(
        args.artifact_root,
        args.terminal_manifest,
        exp065_run_id=args.exp065_run_id,
        expected_prefix=args.expected_prefix,
    )
    atomic_write_json(args.output, authorization)
    print(json.dumps(authorization, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
