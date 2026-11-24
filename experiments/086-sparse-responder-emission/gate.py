#!/usr/bin/env python3
"""Fail-closed authentication of a terminal experiment-077 result."""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import re
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Mapping


EXPECTED_ALPHAS = {0.35, 0.40}
ANCHORS = {
    "pds_cosine": (0.500, 0.9525),
    "expr_mse_unbiased_capped_norm": (0.989, 0.0365),
    "de_wilcoxon_lfc_nmae": (1.0013, 0.400),
    "de_wilcoxon_direction_fidelity_yield_raw": (0.5135, 0.8135),
    "de_wilcoxon_direction_reach_raw": (0.072, 0.968),
    "de_wilcoxon_sig_jaccard": (0.029, 0.399),
}
FLOAT_ABS_TOLERANCE = 1e-12
RUN_ID_RE = re.compile(r"[0-9]{8}T[0-9]{6}Z(?:-[a-z0-9][a-z0-9-]{0,31})?")
SHA256_RE = re.compile(r"[0-9a-f]{64}")
SSM_COMMAND_ID_RE = re.compile(
    r"[0-9a-f]{8}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{12}"
)
REQUIRED_ATTESTED_HASHES = {
    "source_archive",
    "runner",
    "grid",
    "summarizer",
    "exit_code",
    "summary",
    "results",
    "artifact_manifest",
}


def sha256_file(path: Path, block_bytes: int = 8 << 20) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while block := handle.read(block_bytes):
            digest.update(block)
    return digest.hexdigest()


def _regular_file(path: Path, *, label: str) -> Path:
    if path.is_symlink() or not path.is_file():
        raise ValueError(f"{label} must be a regular non-symlink file: {path}")
    return path


def _load_json_file(path: Path, *, label: str) -> object:
    _regular_file(path, label=label)
    try:
        return json.loads(path.read_text())
    except json.JSONDecodeError as error:
        raise ValueError(f"{label} is not valid JSON") from error


def _sha256(value: object, *, label: str) -> str:
    text = str(value)
    if SHA256_RE.fullmatch(text) is None:
        raise ValueError(f"{label} is not a lowercase SHA-256")
    return text


def parse_artifact_checksums(path: Path) -> dict[str, str]:
    records: dict[str, str] = {}
    for raw in _regular_file(path, label="artifact manifest").read_text().splitlines():
        if not raw.strip():
            continue
        pieces = raw.split(maxsplit=1)
        if len(pieces) != 2 or SHA256_RE.fullmatch(pieces[0].lower()) is None:
            raise ValueError(f"malformed checksum line: {raw!r}")
        name = pieces[1].strip()
        if name.startswith("*"):
            name = name[1:]
        if name.startswith("./"):
            name = name[2:]
        relative = Path(name)
        if relative.is_absolute() or ".." in relative.parts:
            raise ValueError(f"unsafe checksum path {name!r}")
        if name in records:
            raise ValueError(f"duplicate checksum path {name}")
        records[name] = pieces[0].lower()
    return records


def alpha_from_tag(tag: str) -> float:
    match = re.search(r"(?:^|_)a(0\.35|0\.4)(?:_|$)", str(tag))
    if match is None:
        raise ValueError(f"cannot identify supported alpha from tag {tag!r}")
    return float(match.group(1))


def scaled(raw: Mapping[str, float]) -> dict[str, float]:
    result: dict[str, float] = {}
    for metric, (baseline, replicate) in ANCHORS.items():
        value = (float(raw[metric]) - baseline) / (replicate - baseline)
        if metric == "expr_mse_unbiased_capped_norm":
            value = max(0.0, min(1.0, value))
        result[metric] = float(value)
    result["avg6"] = float(sum(result.values()) / 6)
    return result


def _finite_metrics(record: Mapping[str, object], *, where: str) -> dict[str, float]:
    output: dict[str, float] = {}
    for metric in ANCHORS:
        try:
            value = float(record[metric])
        except (KeyError, TypeError, ValueError) as error:
            raise ValueError(f"{where} lacks numeric {metric}") from error
        if not math.isfinite(value):
            raise ValueError(f"{where} {metric} is non-finite")
        output[metric] = value
    return output


def _assert_close(observed: float, expected: float, *, where: str) -> None:
    if not math.isclose(
        float(observed),
        float(expected),
        rel_tol=0.0,
        abs_tol=FLOAT_ABS_TOLERANCE,
    ):
        raise ValueError(f"{where}: observed {observed!r}, expected {expected!r}")


def _validate_scaled(
    observed: Mapping[str, object],
    expected: Mapping[str, float],
    *,
    where: str,
) -> None:
    for metric, value in expected.items():
        try:
            got = float(observed[metric])
        except (KeyError, TypeError, ValueError) as error:
            raise ValueError(f"{where} lacks numeric {metric}") from error
        _assert_close(got, value, where=f"{where}.{metric}")


def _single_ssm_command(document: object, *, label: str) -> Mapping[str, object]:
    if not isinstance(document, dict):
        raise ValueError(f"{label} must be a JSON object")
    commands = document.get("Commands")
    if not isinstance(commands, list) or len(commands) != 1:
        raise ValueError(f"{label} must contain exactly one SSM command")
    command = commands[0]
    if not isinstance(command, dict):
        raise ValueError(f"{label} SSM command is malformed")
    return command


def _expected_launch_commands(launch: Mapping[str, object]) -> list[str]:
    bucket = str(launch["bucket"])
    source = str(launch["source_prefix"])
    output = str(launch["output_prefix"])
    run_id = str(launch["run_id"])
    return [
        "set -Eeuo pipefail",
        f"aws s3 cp s3://{bucket}/{source}/run.sh "
        "/root/exp077-run.sh --only-show-errors",
        "chmod +x /root/exp077-run.sh",
        f"systemd-run --unit=vcc-exp077-{run_id} "
        f"--setenv=BUCKET={bucket} --setenv=RUN_ID={run_id} "
        f"--setenv=SOURCE_PREFIX={source} --setenv=OUTPUT_PREFIX={output} "
        "/bin/bash /root/exp077-run.sh",
    ]


def _expected_terminal_commands(
    launch: Mapping[str, object],
    required: Mapping[str, object],
) -> list[str]:
    bucket = str(launch["bucket"])
    attestor_key = str(required["attestor_s3_key"])
    attestor_sha = str(required["attestor_sha256"])
    run_id = str(launch["run_id"])
    return [
        "set -Eeuo pipefail",
        f"aws s3 cp s3://{bucket}/{attestor_key} "
        "/root/exp086-attest-exp077.py --only-show-errors",
        f"printf '%s  %s\\n' '{attestor_sha}' "
        "/root/exp086-attest-exp077.py | sha256sum --check - >/dev/null",
        "python3 /root/exp086-attest-exp077.py "
        f"--run-id {run_id} "
        f"--instance-id {launch['instance_id']} "
        f"--region {launch['region']} "
        f"--source-prefix {launch['source_prefix']} "
        f"--output-prefix {launch['output_prefix']} "
        f"--work-root /work/exp077-{run_id}",
    ]


def _validate_ssm_command(
    document: object,
    *,
    command_id: str,
    instance_id: str,
    expected_commands: list[str],
    label: str,
) -> None:
    command = _single_ssm_command(document, label=label)
    if command.get("CommandId") != command_id:
        raise ValueError(f"{label} command id does not match trusted provenance")
    if command.get("DocumentName") != "AWS-RunShellScript":
        raise ValueError(f"{label} used an unexpected SSM document")
    if command.get("Status") != "Success" or command.get("StatusDetails") != "Success":
        raise ValueError(f"{label} is not successful")
    if command.get("InstanceIds") != [instance_id]:
        raise ValueError(f"{label} instance does not match trusted provenance")
    parameters = command.get("Parameters")
    if (
        not isinstance(parameters, dict)
        or parameters.get("commands") != expected_commands
    ):
        raise ValueError(f"{label} command body does not match the frozen command")


def _validate_provenance(
    root: Path,
    *,
    run_id: str,
    trusted_provenance_path: Path,
    launch_metadata_snapshot_path: Path,
    launch_command_evidence: Path | None,
    terminal_command_evidence: Path | None,
    terminal_invocation_evidence: Path | None,
) -> tuple[dict[str, object], dict[str, str], str]:
    trusted = _load_json_file(
        trusted_provenance_path, label="trusted exp077 provenance"
    )
    if not isinstance(trusted, dict) or trusted.get("schema_version") != 1:
        raise ValueError("trusted exp077 provenance has an unsupported schema")
    launch = trusted.get("launch")
    required = trusted.get("required_terminal_attestation")
    if not isinstance(launch, dict) or not isinstance(required, dict):
        raise ValueError("trusted exp077 provenance lacks launch/attestation policy")
    if launch.get("run_id") != run_id:
        raise ValueError("run id does not match trusted exp077 launch metadata")
    if RUN_ID_RE.fullmatch(run_id) is None:
        raise ValueError("run id is not a canonical timestamp identifier")
    for field in (
        "launch_metadata_sha256",
        "source_archive_sha256",
        "runner_sha256",
        "grid_sha256",
        "summarizer_sha256",
    ):
        _sha256(launch.get(field), label=f"trusted launch {field}")
    snapshot = _load_json_file(
        launch_metadata_snapshot_path, label="exp077 launch metadata snapshot"
    )
    if not isinstance(snapshot, dict):
        raise ValueError("exp077 launch metadata snapshot must be a JSON object")
    if sha256_file(launch_metadata_snapshot_path) != launch["launch_metadata_sha256"]:
        raise ValueError("exp077 launch metadata snapshot hash changed")
    expected_snapshot = {
        "run_id": run_id,
        "region": launch["region"],
        "instance_id": launch["instance_id"],
        "source_s3_prefix": (f"s3://{launch['bucket']}/{launch['source_prefix']}/"),
        "output_s3_prefix": (f"s3://{launch['bucket']}/{launch['output_prefix']}/"),
        "ssm_command_id": launch["ssm_command_id"],
        "cell_eval2_jobs": 2,
        "submission": False,
    }
    if snapshot != expected_snapshot:
        raise ValueError(
            "exp077 launch metadata snapshot disagrees with trusted launch"
        )
    attestor_sha = _sha256(
        required.get("attestor_sha256"), label="trusted attestor_sha256"
    )
    expected_attestor_key = (
        f"input/exp086-provenance-attestor/{attestor_sha}/attest_exp077.py"
    )
    if required.get("attestor_s3_key") != expected_attestor_key:
        raise ValueError("trusted attestor S3 key is not content-addressed")
    if set(required.get("required_hashes", [])) != REQUIRED_ATTESTED_HASHES:
        raise ValueError("trusted terminal attestation hash set changed")

    terminal = trusted.get("terminal_attestation")
    if not isinstance(terminal, dict):
        raise ValueError(
            "immutable exp077 terminal attestation is absent; after exp077 "
            "completes, run the pinned read-only attestor through AWS SSM and "
            "record its command id in exp077-trusted-provenance.json"
        )
    terminal_command_id = str(terminal.get("ssm_command_id", ""))
    if SSM_COMMAND_ID_RE.fullmatch(terminal_command_id) is None:
        raise ValueError("trusted terminal SSM command id is invalid")
    launch_command_id = str(launch.get("ssm_command_id", ""))
    if SSM_COMMAND_ID_RE.fullmatch(launch_command_id) is None:
        raise ValueError("trusted launch SSM command id is invalid")
    if terminal_command_id == launch_command_id:
        raise ValueError("terminal SSM command must be distinct from launch command")
    evidence_paths = {
        "launch SSM evidence": launch_command_evidence,
        "terminal SSM command evidence": terminal_command_evidence,
        "terminal SSM invocation evidence": terminal_invocation_evidence,
    }
    missing_evidence = [name for name, path in evidence_paths.items() if path is None]
    if missing_evidence:
        raise ValueError(f"missing immutable provenance evidence: {missing_evidence}")

    launch_document = _load_json_file(
        launch_command_evidence, label="launch SSM evidence"
    )
    terminal_document = _load_json_file(
        terminal_command_evidence, label="terminal SSM command evidence"
    )
    invocation = _load_json_file(
        terminal_invocation_evidence, label="terminal SSM invocation evidence"
    )
    _validate_ssm_command(
        launch_document,
        command_id=launch_command_id,
        instance_id=str(launch["instance_id"]),
        expected_commands=_expected_launch_commands(launch),
        label="launch SSM evidence",
    )
    _validate_ssm_command(
        terminal_document,
        command_id=terminal_command_id,
        instance_id=str(launch["instance_id"]),
        expected_commands=_expected_terminal_commands(launch, required),
        label="terminal SSM evidence",
    )
    if not isinstance(invocation, dict):
        raise ValueError("terminal SSM invocation evidence must be a JSON object")
    if (
        invocation.get("CommandId") != terminal_command_id
        or invocation.get("InstanceId") != launch["instance_id"]
        or invocation.get("Status") != "Success"
        or invocation.get("StatusDetails") != "Success"
        or invocation.get("ResponseCode") != 0
        or invocation.get("StandardErrorContent") not in ("", None)
    ):
        raise ValueError(
            "terminal SSM invocation is not the trusted successful command"
        )
    stdout = invocation.get("StandardOutputContent")
    if not isinstance(stdout, str) or not stdout.strip():
        raise ValueError("terminal SSM invocation has no attestation stdout")
    try:
        attestation = json.loads(stdout)
    except json.JSONDecodeError as error:
        raise ValueError(
            "terminal SSM stdout is not one canonical JSON record"
        ) from error
    canonical = json.dumps(attestation, sort_keys=True, separators=(",", ":"))
    if stdout.strip() != canonical:
        raise ValueError("terminal SSM attestation is not canonical JSON")
    if not isinstance(attestation, dict):
        raise ValueError("terminal SSM attestation must be a JSON object")
    expected_identity = {
        "schema_version": 1,
        "kind": "exp077-terminal-attestation",
        "run_id": run_id,
        "instance_id": launch["instance_id"],
        "region": launch["region"],
        "source_prefix": launch["source_prefix"],
        "output_prefix": launch["output_prefix"],
        "systemd_unit": f"vcc-exp077-{run_id}.service",
    }
    for key, expected in expected_identity.items():
        if attestation.get(key) != expected:
            raise ValueError(
                f"terminal attestation {key} does not match launch evidence"
            )
    state = attestation.get("systemd_state")
    if not isinstance(state, dict) or (
        state.get("ActiveState") != "inactive"
        or state.get("Result") != "success"
        or state.get("ExecMainStatus") != "0"
    ):
        raise ValueError(
            "terminal attestation does not prove successful unit completion"
        )
    hashes = attestation.get("sha256")
    if not isinstance(hashes, dict) or set(hashes) != REQUIRED_ATTESTED_HASHES:
        raise ValueError("terminal attestation has an incomplete hash set")
    authenticated_hashes = {
        key: _sha256(value, label=f"terminal attestation {key}")
        for key, value in hashes.items()
    }
    expected_source_hashes = {
        "source_archive": launch["source_archive_sha256"],
        "runner": launch["runner_sha256"],
        "grid": launch["grid_sha256"],
        "summarizer": launch["summarizer_sha256"],
    }
    for key, expected in expected_source_hashes.items():
        if authenticated_hashes[key] != expected:
            raise ValueError(
                f"terminal attestation {key} changed from the frozen launch"
            )
    output_files = {
        "exit_code": root / "exit-code.txt",
        "summary": root / "summary.json",
        "results": root / "results.json",
        "artifact_manifest": root / "artifact-sha256.txt",
    }
    for key, path in output_files.items():
        _regular_file(path, label=f"exp077 {key}")
        observed = sha256_file(path)
        if observed != authenticated_hashes[key]:
            raise ValueError(
                f"{key} SHA-256 {observed} != immutable SSM attestation "
                f"{authenticated_hashes[key]}"
            )
    return launch, authenticated_hashes, terminal_command_id


@dataclass(frozen=True)
class Exp077Gate:
    run_id: str
    selected_tag: str
    selected_alpha: float
    summary_sha256: str
    results_sha256: str
    artifact_manifest_sha256: str
    source_archive_sha256: str
    runner_sha256: str
    grid_sha256: str
    summarizer_sha256: str
    launch_ssm_command_id: str
    terminal_attestation_ssm_command_id: str
    candidate_alphas: tuple[float, ...]
    passes_promotion_gate: bool

    def to_dict(self) -> dict[str, object]:
        return asdict(self)


def validate_terminal_exp077(
    root: Path,
    *,
    run_id: str,
    trusted_provenance_path: Path,
    launch_metadata_snapshot_path: Path,
    launch_command_evidence: Path | None = None,
    terminal_command_evidence: Path | None = None,
    terminal_invocation_evidence: Path | None = None,
) -> Exp077Gate:
    launch, attested_hashes, terminal_command_id = _validate_provenance(
        root,
        run_id=run_id,
        trusted_provenance_path=trusted_provenance_path,
        launch_metadata_snapshot_path=launch_metadata_snapshot_path,
        launch_command_evidence=launch_command_evidence,
        terminal_command_evidence=terminal_command_evidence,
        terminal_invocation_evidence=terminal_invocation_evidence,
    )
    if (root / "exit-code.txt").read_text().strip() != "0":
        raise ValueError("experiment 077 terminal exit code is not zero")

    checksums = parse_artifact_checksums(root / "artifact-sha256.txt")
    for name in ("summary.json", "results.json"):
        expected = checksums.get(name)
        if expected is None:
            raise ValueError(f"artifact manifest does not list {name}")
        observed = sha256_file(root / name)
        if observed != expected:
            raise ValueError(f"{name} SHA-256 {observed} != manifest {expected}")
        attested_key = "summary" if name == "summary.json" else "results"
        if observed != attested_hashes[attested_key]:
            raise ValueError(f"{name} is not bound by immutable SSM evidence")

    summary = json.loads((root / "summary.json").read_text())
    results = json.loads((root / "results.json").read_text())
    if not isinstance(summary, dict):
        raise ValueError("experiment 077 summary must be a JSON object")
    if not isinstance(results, list) or len(results) != 2:
        raise ValueError("experiment 077 results must contain exactly two runs")
    candidates = summary.get("candidates")
    if not isinstance(candidates, list) or len(candidates) != 2:
        raise ValueError("experiment 077 summary must contain exactly two candidates")
    if any(not isinstance(candidate, dict) for candidate in candidates):
        raise ValueError("experiment 077 candidates must be JSON objects")
    if any(not isinstance(result, dict) for result in results):
        raise ValueError("experiment 077 results must be JSON objects")
    by_tag = {str(candidate.get("tag")): candidate for candidate in candidates}
    if len(by_tag) != 2:
        raise ValueError("experiment 077 candidate tags are not unique")
    result_by_tag = {str(result.get("tag")): result for result in results}
    if len(result_by_tag) != 2 or set(result_by_tag) != set(by_tag):
        raise ValueError("experiment 077 results and summary tags disagree")
    alphas = tuple(sorted(alpha_from_tag(tag) for tag in by_tag))
    if set(alphas) != EXPECTED_ALPHAS:
        raise ValueError(f"experiment 077 alpha set {alphas} is not {EXPECTED_ALPHAS}")

    baseline = summary.get("baseline")
    if not isinstance(baseline, dict):
        raise ValueError("experiment 077 summary lacks baseline")
    baseline_raw_record = baseline.get("raw")
    baseline_scaled_record = baseline.get("scaled")
    if not isinstance(baseline_raw_record, dict) or not isinstance(
        baseline_scaled_record, dict
    ):
        raise ValueError("experiment 077 baseline lacks raw/scaled metrics")
    baseline_raw = _finite_metrics(baseline_raw_record, where="baseline.raw")
    baseline_scaled = scaled(baseline_raw)
    _validate_scaled(baseline_scaled_record, baseline_scaled, where="baseline.scaled")

    for tag, candidate in by_tag.items():
        result = result_by_tag[tag]
        if "error" in result:
            raise ValueError(f"experiment 077 result {tag} is failed or malformed")
        raw_record = candidate.get("raw")
        result_raw_record = result.get("raw")
        if not isinstance(raw_record, dict) or not isinstance(result_raw_record, dict):
            raise ValueError(f"experiment 077 candidate {tag} lacks raw metrics")
        raw = _finite_metrics(raw_record, where=f"candidate[{tag}].raw")
        result_raw = _finite_metrics(result_raw_record, where=f"results[{tag}].raw")
        for metric in ANCHORS:
            _assert_close(
                raw[metric],
                result_raw[metric],
                where=f"candidate/results[{tag}].{metric}",
            )
        expected_scaled = scaled(raw)
        candidate_scaled = candidate.get("scaled")
        result_scaled = result.get("official")
        if not isinstance(candidate_scaled, dict) or not isinstance(
            result_scaled, dict
        ):
            raise ValueError(f"experiment 077 candidate {tag} lacks scaled metrics")
        _validate_scaled(
            candidate_scaled,
            expected_scaled,
            where=f"candidate[{tag}].scaled",
        )
        _validate_scaled(
            result_scaled,
            expected_scaled,
            where=f"results[{tag}].official",
        )
        expected_delta_avg6 = expected_scaled["avg6"] - baseline_scaled["avg6"]
        expected_delta_pds = raw["pds_cosine"] - baseline_raw["pds_cosine"]
        expected_positive_mse = (
            raw["expr_mse_unbiased_capped_norm"]
            < ANCHORS["expr_mse_unbiased_capped_norm"][0]
        )
        expected_individual = bool(
            expected_positive_mse
            and expected_delta_pds >= -0.003
            and expected_delta_avg6 >= 0.003
        )
        _assert_close(
            float(candidate.get("delta_avg6")),
            expected_delta_avg6,
            where=f"candidate[{tag}].delta_avg6",
        )
        _assert_close(
            float(candidate.get("delta_pds_raw")),
            expected_delta_pds,
            where=f"candidate[{tag}].delta_pds_raw",
        )
        if candidate.get("positive_mse") is not expected_positive_mse:
            raise ValueError(f"candidate {tag} positive_mse is inconsistent")
        if candidate.get("passes_individual_gate") is not expected_individual:
            raise ValueError(f"candidate {tag} individual gate is inconsistent")

    selected_tag = str(summary.get("selected"))
    if selected_tag not in by_tag:
        raise ValueError("experiment 077 selected tag is absent")
    expected_selected = max(
        candidates,
        key=lambda candidate: float(candidate["scaled"]["avg6"]),
    )["tag"]
    if selected_tag != expected_selected:
        raise ValueError("experiment 077 selected tag is not the best candidate")
    selected = by_tag[selected_tag]
    if selected.get("passes_individual_gate") is not True:
        raise ValueError("experiment 077 selected candidate failed its individual gate")
    adjacent_positive = all(float(item["delta_avg6"]) > 0 for item in candidates)
    if summary.get("adjacent_alphas_both_positive") is not adjacent_positive:
        raise ValueError("experiment 077 adjacent-alpha gate is inconsistent")
    expected_promotion = bool(selected["passes_individual_gate"] and adjacent_positive)
    if summary.get("passes_promotion_gate") is not expected_promotion:
        raise ValueError("experiment 077 promotion gate is inconsistent")
    if not expected_promotion:
        raise ValueError("experiment 077 did not pass its promotion gate")
    if summary.get("submission_authorized") is not False:
        raise ValueError("unexpected experiment 077 authorization semantics")

    return Exp077Gate(
        run_id=run_id,
        selected_tag=selected_tag,
        selected_alpha=alpha_from_tag(selected_tag),
        summary_sha256=attested_hashes["summary"],
        results_sha256=attested_hashes["results"],
        artifact_manifest_sha256=attested_hashes["artifact_manifest"],
        source_archive_sha256=attested_hashes["source_archive"],
        runner_sha256=attested_hashes["runner"],
        grid_sha256=attested_hashes["grid"],
        summarizer_sha256=attested_hashes["summarizer"],
        launch_ssm_command_id=str(launch["ssm_command_id"]),
        terminal_attestation_ssm_command_id=terminal_command_id,
        candidate_alphas=alphas,
        passes_promotion_gate=True,
    )


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--exp077-output", type=Path, required=True)
    parser.add_argument("--run-id", required=True)
    parser.add_argument("--trusted-provenance", type=Path, required=True)
    parser.add_argument("--launch-metadata-snapshot", type=Path, required=True)
    parser.add_argument("--launch-command-evidence", type=Path)
    parser.add_argument("--terminal-command-evidence", type=Path)
    parser.add_argument("--terminal-invocation-evidence", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    gate = validate_terminal_exp077(
        args.exp077_output,
        run_id=args.run_id,
        trusted_provenance_path=args.trusted_provenance,
        launch_metadata_snapshot_path=args.launch_metadata_snapshot,
        launch_command_evidence=args.launch_command_evidence,
        terminal_command_evidence=args.terminal_command_evidence,
        terminal_invocation_evidence=args.terminal_invocation_evidence,
    )
    if args.output.is_symlink() or args.output.exists():
        raise FileExistsError(f"refusing existing gate output {args.output}")
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(gate.to_dict(), indent=2) + "\n")
    print(json.dumps(gate.to_dict(), indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
