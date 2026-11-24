#!/usr/bin/env python3
"""Bounded tests for the sparse-responder emitter and exp077 gate."""

from __future__ import annotations

import hashlib
import importlib.util
import json
import os
import subprocess
import sys
import tempfile
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import scipy.sparse as sp
from cell_eval2.moments import jackknife_correction as official_jackknife_correction


HERE = Path(__file__).resolve().parent


def load_module(name: str, filename: str):
    spec = importlib.util.spec_from_file_location(name, HERE / filename)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


EMITTER = load_module("exp086_emitter_test", "emitter.py")
GATE = load_module("exp086_gate_test", "gate.py")
RUN = load_module("exp086_run_test", "run.py")


def test_exact_columns_nonnegative_rows_and_determinism() -> None:
    rng = np.random.default_rng(860)
    template = rng.poisson(3.0, size=(40, 71)).astype(np.int64)
    template[rng.random(template.shape) < 0.45] = 0
    totals = np.asarray(template.sum(axis=0), dtype=np.int64)
    moved = np.minimum(totals[:10] // 4, totals[:10])
    totals[:10] -= moved
    totals[10:20] += moved
    first, first_diag = EMITTER.emit_sparse_responders(
        sp.csr_matrix(template),
        column_totals=totals,
        responder_fraction=0.40,
        key="synthetic/T0",
        seed=861,
        control_reference_totals=template.sum(axis=0),
        scored_gene_mask=np.arange(template.shape[1]) < 50,
    )
    second, second_diag = EMITTER.emit_sparse_responders(
        sp.csr_matrix(template),
        column_totals=totals,
        responder_fraction=0.40,
        key="synthetic/T0",
        seed=861,
        control_reference_totals=template.sum(axis=0),
        scored_gene_mask=np.arange(template.shape[1]) < 50,
    )
    np.testing.assert_array_equal(first.toarray(), second.toarray())
    assert first_diag == second_diag
    np.testing.assert_array_equal(np.asarray(first.sum(axis=0)).ravel(), totals)
    assert first.data.min(initial=0) >= 0
    rows = np.asarray(first.sum(axis=1)).ravel()
    assert rows.min() > 0
    assert rows.max() <= EMITTER.MAX_COUNTS_PER_CELL
    assert first_diag.max_abs_column_error == 0
    np.testing.assert_array_equal(
        np.asarray(first.sum(axis=1)).ravel(),
        template.sum(axis=1),
    )
    assert (
        first_diag.input_column_totals_sha256 == first_diag.output_column_totals_sha256
    )


def test_negative_effect_infeasibility_is_repaired_and_reported() -> None:
    template = np.full((10, 4), 10, dtype=np.int64)
    totals = np.asarray([20, 180, 100, 100], dtype=np.int64)
    emitted, diag = EMITTER.emit_sparse_responders(
        template,
        column_totals=totals,
        responder_fraction=0.30,
        key="negative",
        seed=862,
        control_reference_totals=template.sum(axis=0),
        scored_gene_mask=np.ones(4, dtype=bool),
    )
    np.testing.assert_array_equal(np.asarray(emitted.sum(axis=0)).ravel(), totals)
    assert emitted.data.min(initial=0) >= 0
    assert diag.nonresponder_clipped_counts > 0
    assert diag.nonresponder_absolute_change_counts > 0
    assert diag.scored_nonresponder_clipped_genes >= 1
    assert diag.representative_subset_screening_fraction_scored == 0.8
    assert not diag.selected_subset_exact_control_feasible_scored


def test_selected_subset_feasibility_is_not_overridden_by_aggregate_bound() -> None:
    key = "nonrepresentative-subset"
    seed = 863
    responders = EMITTER.deterministic_responder_indices(4, 0.50, key=key, seed=seed)
    nonresponders = np.setdiff1d(np.arange(4), responders)
    template = np.zeros((4, 2), dtype=np.int64)
    template[responders, 0] = 10
    template[responders, 1] = 1
    template[nonresponders, 1] = 10
    baseline = template.sum(axis=0)
    totals = baseline + np.asarray([-15, 15], dtype=np.int64)
    emitted, diag = EMITTER.emit_sparse_responders(
        template,
        column_totals=totals,
        responder_fraction=0.50,
        key=key,
        seed=seed,
    )
    np.testing.assert_array_equal(
        emitted.toarray()[nonresponders], template[nonresponders]
    )
    assert diag.representative_subset_screening_fraction_all == 0.75
    assert diag.actual_responder_fraction == 0.50
    assert diag.nonresponder_absolute_change_counts == 0
    assert diag.scored_nonresponder_absolute_change_counts == 0
    assert diag.selected_subset_exact_control_feasible_all
    assert diag.selected_subset_exact_control_feasible_scored


def test_nonresponders_are_exact_control_rows_when_no_clipping_is_needed() -> None:
    rng = np.random.default_rng(863)
    template = rng.integers(0, 5, size=(20, 12), dtype=np.int64)
    totals = np.asarray(template.sum(axis=0), dtype=np.int64)
    moved = min(int(totals[0] // 4), 4)
    totals[0] -= moved
    totals[1] += moved
    emitted, diag = EMITTER.emit_sparse_responders(
        template,
        column_totals=totals,
        responder_fraction=0.50,
        key="positive-only",
        seed=864,
    )
    responders = EMITTER.deterministic_responder_indices(
        len(template), 0.50, key="positive-only", seed=864
    )
    nonresponders = np.setdiff1d(np.arange(len(template)), responders)
    np.testing.assert_array_equal(
        emitted.toarray()[nonresponders], template[nonresponders]
    )
    assert diag.nonresponder_clipped_counts == 0
    assert diag.nonresponder_absolute_change_counts == 0
    assert diag.unchanged_nonresponder_row_fraction == 1.0


def test_exact_margins_fix_pseudobulk_but_not_jackknife_or_cpm_mean() -> None:
    rng = np.random.default_rng(865)
    template = rng.poisson(4.0, size=(50, 40)).astype(np.int64)
    template[rng.random(template.shape) < 0.40] = 0
    totals = np.asarray(template.sum(axis=0), dtype=np.int64)
    moved = np.minimum(totals[:8] // 4, totals[:8])
    totals[:8] -= moved
    totals[8:16] += moved
    sparse, _ = EMITTER.emit_sparse_responders(
        template,
        column_totals=totals,
        responder_fraction=0.50,
        key="metric-separation",
        seed=866,
    )
    # A different exact-margin layout: even integer allocation.
    even = np.zeros_like(template)
    for gene, total in enumerate(totals):
        base, remainder = divmod(int(total), len(template))
        even[:, gene] = base
        even[:remainder, gene] += 1
    even = sp.csr_matrix(even)
    np.testing.assert_array_equal(
        EMITTER.scorer_pseudobulk(sparse), EMITTER.scorer_pseudobulk(even)
    )
    assert not np.array_equal(
        EMITTER.arithmetic_mean_cpm(sparse),
        EMITTER.arithmetic_mean_cpm(even),
    )
    assert not np.isclose(
        EMITTER.jackknife_correction(sparse, chunk_genes=11),
        EMITTER.jackknife_correction(even, chunk_genes=11),
    )


def test_jackknife_implementation_matches_cell_eval2() -> None:
    rng = np.random.default_rng(867)
    counts = rng.poisson(2.0, size=(37, 83)).astype(np.int32)
    counts[rng.random(counts.shape) < 0.55] = 0
    matrix = sp.csr_matrix(counts)
    observed = EMITTER.jackknife_correction(matrix, chunk_genes=13)
    expected = official_jackknife_correction(
        matrix,
        np.zeros(matrix.shape[0], dtype=np.intp),
        1,
        EMITTER.BULK_TARGET_SUM,
        chunk=13,
    )[0]
    np.testing.assert_allclose(observed, expected, rtol=1e-13, atol=1e-13)


def test_row_cap_and_zero_row_fail_closed() -> None:
    template = np.ones((10, 3), dtype=np.int64)
    try:
        EMITTER.emit_sparse_responders(
            template,
            column_totals=np.asarray([10, 10, 10]),
            responder_fraction=0.10,
            key="cap",
            max_counts_per_cell=2,
        )
    except ValueError as error:
        assert "exceeds cap" in str(error)
    else:
        raise AssertionError("row-cap violation was accepted")

    try:
        EMITTER.emit_sparse_responders(
            template,
            column_totals=np.asarray([10, 10, 10]),
            responder_fraction=0.50,
            key="zero-row",
            minimum_counts_per_cell=4,
        )
    except ValueError as error:
        assert "below minimum" in str(error)
    else:
        raise AssertionError("zero-library rows were accepted")


def _write_terminal_fixture(
    root: Path,
    *,
    passes: bool,
    exit_code: int = 0,
) -> dict[str, object]:
    tags = [
        "g0.65_q1.0_a0.35_spread0.25_cpm5.0_seed0",
        "g0.65_q1.0_a0.4_spread0.25_cpm5.0_seed0",
    ]
    baseline_raw = {
        "pds_cosine": 0.790,
        "expr_mse_unbiased_capped_norm": 0.995,
        "de_wilcoxon_lfc_nmae": 0.970,
        "de_wilcoxon_direction_fidelity_yield_raw": 0.510,
        "de_wilcoxon_direction_reach_raw": 0.090,
        "de_wilcoxon_sig_jaccard": 0.035,
    }
    metrics = {
        "pds_cosine": 0.80,
        "expr_mse_unbiased_capped_norm": 0.98,
        "de_wilcoxon_lfc_nmae": 0.95,
        "de_wilcoxon_direction_fidelity_yield_raw": 0.52,
        "de_wilcoxon_direction_reach_raw": 0.10,
        "de_wilcoxon_sig_jaccard": 0.04,
    }
    baseline_scaled = GATE.scaled(baseline_raw)
    candidates = []
    results = []
    for tag in tags:
        score = GATE.scaled(metrics)
        delta_avg6 = score["avg6"] - baseline_scaled["avg6"]
        delta_pds = metrics["pds_cosine"] - baseline_raw["pds_cosine"]
        positive_mse = (
            metrics["expr_mse_unbiased_capped_norm"]
            < GATE.ANCHORS["expr_mse_unbiased_capped_norm"][0]
        )
        individual = bool(positive_mse and delta_pds >= -0.003 and delta_avg6 >= 0.003)
        candidates.append(
            {
                "tag": tag,
                "raw": metrics,
                "scaled": score,
                "delta_avg6": delta_avg6,
                "delta_pds_raw": delta_pds,
                "positive_mse": positive_mse,
                "passes_individual_gate": individual if passes else False,
            }
        )
        results.append({"tag": tag, "raw": metrics, "official": score})
    summary = {
        "baseline": {"raw": baseline_raw, "scaled": baseline_scaled},
        "candidates": candidates,
        "adjacent_alphas_both_positive": passes,
        "selected": tags[0],
        "passes_promotion_gate": passes,
        "submission_authorized": False,
    }
    (root / "summary.json").write_text(json.dumps(summary))
    (root / "results.json").write_text(json.dumps(results))
    (root / "exit-code.txt").write_text(f"{exit_code}\n")
    lines = []
    for name in ("summary.json", "results.json"):
        digest = hashlib.sha256((root / name).read_bytes()).hexdigest()
        lines.append(f"{digest}  ./{name}")
    (root / "artifact-sha256.txt").write_text("\n".join(lines) + "\n")
    run_id = "20260823T060938Z"
    launch_command_id = "11111111-1111-4111-8111-111111111111"
    terminal_command_id = "22222222-2222-4222-8222-222222222222"
    launch = {
        "run_id": run_id,
        "region": "us-east-1",
        "instance_id": "i-0123456789abcdef0",
        "bucket": "test-vcc-bucket",
        "source_prefix": f"input/exp077-low-alpha-union/{run_id}",
        "output_prefix": f"output/exp077-low-alpha-union/{run_id}",
        "ssm_command_id": launch_command_id,
        "launch_metadata_sha256": "",
        "source_archive_sha256": "2" * 64,
        "runner_sha256": "3" * 64,
        "grid_sha256": "4" * 64,
        "summarizer_sha256": "5" * 64,
    }
    launch_snapshot = {
        "run_id": run_id,
        "region": launch["region"],
        "instance_id": launch["instance_id"],
        "source_s3_prefix": (f"s3://{launch['bucket']}/{launch['source_prefix']}/"),
        "output_s3_prefix": (f"s3://{launch['bucket']}/{launch['output_prefix']}/"),
        "ssm_command_id": launch_command_id,
        "cell_eval2_jobs": 2,
        "submission": False,
    }
    snapshot_path = root / "launch-metadata.snapshot.json"
    snapshot_path.write_text(json.dumps(launch_snapshot, indent=2) + "\n")
    launch["launch_metadata_sha256"] = GATE.sha256_file(snapshot_path)
    required = {
        "authority": "AWS Systems Manager Run Command",
        "attestor_sha256": "6" * 64,
        "attestor_s3_key": (
            f"input/exp086-provenance-attestor/{'6' * 64}/attest_exp077.py"
        ),
        "required_hashes": sorted(GATE.REQUIRED_ATTESTED_HASHES),
    }
    trusted = {
        "schema_version": 1,
        "launch": launch,
        "terminal_attestation": {"ssm_command_id": terminal_command_id},
        "required_terminal_attestation": required,
    }
    trusted_path = root / "trusted.json"
    trusted_path.write_text(json.dumps(trusted))

    attested_hashes = {
        "source_archive": launch["source_archive_sha256"],
        "runner": launch["runner_sha256"],
        "grid": launch["grid_sha256"],
        "summarizer": launch["summarizer_sha256"],
        "exit_code": GATE.sha256_file(root / "exit-code.txt"),
        "summary": GATE.sha256_file(root / "summary.json"),
        "results": GATE.sha256_file(root / "results.json"),
        "artifact_manifest": GATE.sha256_file(root / "artifact-sha256.txt"),
    }
    attestation = {
        "schema_version": 1,
        "kind": "exp077-terminal-attestation",
        "run_id": run_id,
        "instance_id": launch["instance_id"],
        "region": launch["region"],
        "source_prefix": launch["source_prefix"],
        "output_prefix": launch["output_prefix"],
        "systemd_unit": f"vcc-exp077-{run_id}.service",
        "systemd_state": {
            "ActiveState": "inactive",
            "SubState": "dead",
            "Result": "success",
            "ExecMainStatus": "0",
        },
        "sha256": attested_hashes,
    }
    launch_document = {
        "Commands": [
            {
                "CommandId": launch_command_id,
                "DocumentName": "AWS-RunShellScript",
                "Status": "Success",
                "StatusDetails": "Success",
                "InstanceIds": [launch["instance_id"]],
                "Parameters": {"commands": GATE._expected_launch_commands(launch)},
            }
        ]
    }
    terminal_document = {
        "Commands": [
            {
                "CommandId": terminal_command_id,
                "DocumentName": "AWS-RunShellScript",
                "Status": "Success",
                "StatusDetails": "Success",
                "InstanceIds": [launch["instance_id"]],
                "Parameters": {
                    "commands": GATE._expected_terminal_commands(launch, required)
                },
            }
        ]
    }
    invocation = {
        "CommandId": terminal_command_id,
        "InstanceId": launch["instance_id"],
        "Status": "Success",
        "StatusDetails": "Success",
        "ResponseCode": 0,
        "StandardOutputContent": json.dumps(
            attestation, sort_keys=True, separators=(",", ":")
        )
        + "\n",
    }
    launch_path = root / "launch-command.json"
    terminal_path = root / "terminal-command.json"
    invocation_path = root / "terminal-invocation.json"
    launch_path.write_text(json.dumps(launch_document))
    terminal_path.write_text(json.dumps(terminal_document))
    invocation_path.write_text(json.dumps(invocation))
    return {
        "run_id": run_id,
        "trusted_provenance_path": trusted_path,
        "launch_metadata_snapshot_path": snapshot_path,
        "launch_command_evidence": launch_path,
        "terminal_command_evidence": terminal_path,
        "terminal_invocation_evidence": invocation_path,
    }


def test_terminal_gate_accepts_only_authenticated_passing_result() -> None:
    with tempfile.TemporaryDirectory(dir=HERE) as directory:
        root = Path(directory)
        evidence = _write_terminal_fixture(root, passes=True)
        gate = GATE.validate_terminal_exp077(root, **evidence)
        assert gate.selected_alpha == 0.35
        assert gate.candidate_alphas == (0.35, 0.4)

        data = json.loads((root / "summary.json").read_text())
        data["passes_promotion_gate"] = False
        (root / "summary.json").write_text(json.dumps(data))
        try:
            GATE.validate_terminal_exp077(root, **evidence)
        except ValueError as error:
            assert "immutable SSM attestation" in str(error)
        else:
            raise AssertionError("tampered summary was accepted")


def test_terminal_gate_rejects_coherent_rehash_of_colocated_manifest() -> None:
    with tempfile.TemporaryDirectory(dir=HERE) as directory:
        root = Path(directory)
        evidence = _write_terminal_fixture(root, passes=True)
        data = json.loads((root / "summary.json").read_text())
        data["forged_but_semantically_ignored"] = True
        (root / "summary.json").write_text(json.dumps(data))
        lines = []
        for name in ("summary.json", "results.json"):
            digest = hashlib.sha256((root / name).read_bytes()).hexdigest()
            lines.append(f"{digest}  ./{name}")
        (root / "artifact-sha256.txt").write_text("\n".join(lines) + "\n")
        try:
            GATE.validate_terminal_exp077(root, **evidence)
        except ValueError as error:
            assert "immutable SSM attestation" in str(error)
        else:
            raise AssertionError("coherently rehashed S3 output was accepted")


def test_terminal_gate_rejects_absent_or_failed_terminal_marker() -> None:
    with tempfile.TemporaryDirectory(dir=HERE) as directory:
        root = Path(directory)
        evidence = _write_terminal_fixture(root, passes=True, exit_code=1)
        try:
            GATE.validate_terminal_exp077(root, **evidence)
        except ValueError as error:
            assert "exit code" in str(error)
        else:
            raise AssertionError("failed exp077 terminal state was accepted")
        (root / "exit-code.txt").unlink()
        try:
            GATE.validate_terminal_exp077(root, **evidence)
        except ValueError as error:
            assert "regular non-symlink file" in str(error)
        else:
            raise AssertionError("nonterminal exp077 state was accepted")


def test_terminal_gate_rejects_tampered_ssm_command_body() -> None:
    with tempfile.TemporaryDirectory(dir=HERE) as directory:
        root = Path(directory)
        evidence = _write_terminal_fixture(root, passes=True)
        command_path = Path(evidence["terminal_command_evidence"])
        document = json.loads(command_path.read_text())
        document["Commands"][0]["Parameters"]["commands"][-1] += " --forged"
        command_path.write_text(json.dumps(document))
        try:
            GATE.validate_terminal_exp077(root, **evidence)
        except ValueError as error:
            assert "frozen command" in str(error)
        else:
            raise AssertionError("tampered SSM command body was accepted")


def test_production_provenance_fails_closed_until_terminal_attestation() -> None:
    with tempfile.TemporaryDirectory(dir=HERE) as directory:
        try:
            GATE.validate_terminal_exp077(
                Path(directory),
                run_id="20260823T060938Z",
                trusted_provenance_path=HERE / "exp077-trusted-provenance.json",
                launch_metadata_snapshot_path=(
                    HERE / "exp077-launch-metadata.snapshot.json"
                ),
            )
        except ValueError as error:
            assert "immutable exp077 terminal attestation is absent" in str(error)
        else:
            raise AssertionError("missing immutable terminal attestation was accepted")


def test_run_protocol_rejects_fraction_seed_and_scorer_changes() -> None:
    valid = SimpleNamespace(
        fractions=[0.35, 0.50, 0.65],
        maximum_scorers=2,
        seed=8_602_026,
    )
    RUN.validate_frozen_protocol(valid)
    for field, value in (
        ("fractions", [0.35, 0.65]),
        ("maximum_scorers", 1),
        ("maximum_scorers", 2.1),
        ("seed", 0),
        ("seed", 8_602_026.5),
    ):
        changed = SimpleNamespace(**vars(valid))
        setattr(changed, field, value)
        try:
            RUN.validate_frozen_protocol(changed)
        except ValueError:
            pass
        else:
            raise AssertionError(f"changed frozen {field} was accepted")


def test_run_rejects_existing_empty_or_symlink_output() -> None:
    base = {
        "fractions": [0.35, 0.50, 0.65],
        "maximum_scorers": 2,
        "seed": 8_602_026,
    }
    with tempfile.TemporaryDirectory() as directory:
        root = Path(directory).resolve()
        existing = root / "existing"
        existing.mkdir()
        args = SimpleNamespace(**base, output=existing)
        try:
            RUN.run(args)
        except FileExistsError:
            pass
        else:
            raise AssertionError("existing empty output directory was accepted")

        target = root / "target"
        target.mkdir()
        symlink = root / "output-link"
        symlink.symlink_to(target, target_is_directory=True)
        args = SimpleNamespace(**base, output=symlink)
        try:
            RUN.run(args)
        except FileExistsError:
            pass
        else:
            raise AssertionError("symlink output directory was accepted")


def _runner_environment(*, run_id: str, work_root: Path) -> dict[str, str]:
    environment = os.environ.copy()
    environment.update(
        {
            "BUCKET": "<your-s3-bucket>",
            "RUN_ID": run_id,
            "EXP077_RUN_ID": "20260823T060938Z",
            "EXPECTED_SOURCE_SHA256": "1" * 64,
            "EXPECTED_RUNNER_SHA256": hashlib.sha256(
                (HERE / "aws/run.sh").read_bytes()
            ).hexdigest(),
            "EXPECTED_RESERVATION_SHA256": "2" * 64,
            "WORK_ROOT": str(work_root),
        }
    )
    return environment


def test_shell_identifier_and_reused_work_guards_are_early() -> None:
    invalid = subprocess.run(
        ["/bin/bash", str(HERE / "aws/run.sh")],
        env=_runner_environment(
            run_id="../../tmp/victim",
            work_root=Path("/tmp/unused-exp086"),
        ),
        check=False,
        capture_output=True,
        text=True,
    )
    assert invalid.returncode == 67
    assert "unsafe run id" in invalid.stderr

    with tempfile.TemporaryDirectory() as directory:
        work_root = Path(directory).resolve() / "exp086"
        work_root.mkdir()
        run_id = "20990101T010203Z"
        reused = work_root / run_id
        reused.mkdir()
        sentinel = reused / "must-survive.txt"
        sentinel.write_text("preserve me")
        process = subprocess.run(
            ["/bin/bash", str(HERE / "aws/run.sh")],
            env=_runner_environment(run_id=run_id, work_root=work_root),
            check=False,
            capture_output=True,
            text=True,
        )
        assert process.returncode == 79
        assert sentinel.read_text() == "preserve me"

    with tempfile.TemporaryDirectory() as directory:
        parent = Path(directory).resolve()
        work_root = parent / "exp086"
        work_root.mkdir()
        run_id = "20990101T010204Z"
        target = parent / "elsewhere"
        target.mkdir()
        sentinel = target / "must-survive.txt"
        sentinel.write_text("preserve me too")
        (work_root / run_id).symlink_to(target, target_is_directory=True)
        process = subprocess.run(
            ["/bin/bash", str(HERE / "aws/run.sh")],
            env=_runner_environment(run_id=run_id, work_root=work_root),
            check=False,
            capture_output=True,
            text=True,
        )
        assert process.returncode == 79
        assert sentinel.read_text() == "preserve me too"


def test_launcher_rejects_unsafe_identifiers_before_aws() -> None:
    base = os.environ.copy()
    base.update(
        {
            "I_HAVE_EXPLICIT_CLEARANCE": "YES",
            "INSTANCE_ID": "i-0123456789abcdef0",
            "AWS_REGION": "us-east-1",
            "BUCKET": "<your-s3-bucket>",
            "EXP077_RUN_ID": "20260823T060938Z",
            "RUN_ID": "20990101T010205Z",
        }
    )
    cases = (
        ("RUN_ID", "../../tmp/victim"),
        ("EXP077_RUN_ID", "bad;touch-pwned"),
        ("INSTANCE_ID", "i-good;whoami"),
        ("AWS_REGION", "us-east-1;whoami"),
        ("BUCKET", "bucket/escape"),
    )
    for field, value in cases:
        environment = dict(base)
        environment[field] = value
        process = subprocess.run(
            ["/bin/bash", str(HERE / "aws/launch.sh")],
            env=environment,
            check=False,
            capture_output=True,
            text=True,
        )
        assert process.returncode == 67, (field, process.stderr)
        assert "unsafe" in process.stderr


def test_launcher_rejects_existing_local_metadata_and_missing_attestation() -> None:
    run_id = "20990101T010206Z"
    environment = os.environ.copy()
    environment.update(
        {
            "I_HAVE_EXPLICIT_CLEARANCE": "YES",
            "INSTANCE_ID": "i-0123456789abcdef0",
            "AWS_REGION": "us-east-1",
            "BUCKET": "<your-s3-bucket>",
            "EXP077_RUN_ID": "20260823T060938Z",
            "RUN_ID": run_id,
        }
    )
    metadata = HERE / f"launch-metadata-{run_id}.json"
    try:
        metadata.write_text("do not overwrite")
        collision = subprocess.run(
            ["/bin/bash", str(HERE / "aws/launch.sh")],
            env=environment,
            check=False,
            capture_output=True,
            text=True,
        )
        assert collision.returncode == 79
        assert metadata.read_text() == "do not overwrite"
    finally:
        metadata.unlink(missing_ok=True)

    blocked = subprocess.run(
        ["/bin/bash", str(HERE / "aws/launch.sh")],
        env=environment,
        check=False,
        capture_output=True,
        text=True,
    )
    assert blocked.returncode != 0
    assert "immutable exp077 terminal attestation is absent" in blocked.stderr
    assert not metadata.exists()


def test_aws_runner_is_disjoint_gated_and_non_submitting() -> None:
    launch = (HERE / "aws/launch.sh").read_text()
    runner = (HERE / "aws/run.sh").read_text()
    assert "I_HAVE_EXPLICIT_CLEARANCE" in launch
    assert "run-instances" not in launch
    assert "start-instances" not in launch
    assert "vcc submit" not in launch
    assert "vcc submit" not in runner
    assert "rm -rf" not in launch
    assert "rm -rf" not in runner
    assert "--if-none-match '*'" in launch
    assert 'ARCHIVE="${PREFLIGHT}/source.${RUN_ID}.tar.gz"' in launch
    assert "launch-metadata-${RUN_ID}.json" in launch
    assert "experiments/080-" not in launch
    assert "experiments/080-" not in runner
    gate_position = runner.index("gate.py")
    donor_position = runner.index("input/clean-hct")
    scorer_position = runner.index("experiments/086-sparse-responder-emission/run.py")
    assert gate_position < donor_position
    assert gate_position < scorer_position


def test_local_memory_guard_terminates_an_over_limit_child() -> None:
    with tempfile.TemporaryDirectory(dir=HERE) as directory:
        report = Path(directory) / "memory.json"
        process = subprocess.run(
            [
                sys.executable,
                str(HERE / "memory_guard.py"),
                "--max-rss-gib",
                "0.05",
                "--poll-seconds",
                "0.02",
                "--report",
                str(report),
                "--",
                sys.executable,
                "-c",
                (
                    "import time; "
                    "payload = bytearray(96 * 1024 * 1024); "
                    "time.sleep(5); "
                    "print(len(payload))"
                ),
            ],
            check=False,
            capture_output=True,
            text=True,
            timeout=15,
        )
        assert process.returncode == 137
        record = json.loads(report.read_text())
        assert record["killed_for_memory"] is True
        assert record["peak_process_tree_rss_gib"] > 0.05


def main() -> int:
    tests = [
        test_exact_columns_nonnegative_rows_and_determinism,
        test_negative_effect_infeasibility_is_repaired_and_reported,
        test_selected_subset_feasibility_is_not_overridden_by_aggregate_bound,
        test_nonresponders_are_exact_control_rows_when_no_clipping_is_needed,
        test_exact_margins_fix_pseudobulk_but_not_jackknife_or_cpm_mean,
        test_jackknife_implementation_matches_cell_eval2,
        test_row_cap_and_zero_row_fail_closed,
        test_terminal_gate_accepts_only_authenticated_passing_result,
        test_terminal_gate_rejects_coherent_rehash_of_colocated_manifest,
        test_terminal_gate_rejects_absent_or_failed_terminal_marker,
        test_terminal_gate_rejects_tampered_ssm_command_body,
        test_production_provenance_fails_closed_until_terminal_attestation,
        test_run_protocol_rejects_fraction_seed_and_scorer_changes,
        test_run_rejects_existing_empty_or_symlink_output,
        test_shell_identifier_and_reused_work_guards_are_early,
        test_launcher_rejects_unsafe_identifiers_before_aws,
        test_launcher_rejects_existing_local_metadata_and_missing_attestation,
        test_aws_runner_is_disjoint_gated_and_non_submitting,
        test_local_memory_guard_terminates_an_over_limit_child,
    ]
    for test in tests:
        test()
        print(f"PASS {test.__name__}", flush=True)
    print(
        f"{len(tests)} tests passed; peak RSS={EMITTER.process_peak_rss_gib():.3f} GiB",
        flush=True,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
