#!/usr/bin/env python3
"""Synthetic tests for experiment-110 artifact and manifest contracts."""

from __future__ import annotations

import hashlib
import sys
import tempfile
from pathlib import Path

import numpy as np


HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

import arbitration as ARB  # noqa: E402
import artifact_contract as CONTRACT  # noqa: E402


def fake_selection_arrays() -> dict[str, np.ndarray]:
    selected_index = 0
    model = ARB.RidgeLogisticModel(
        ARB.FIXED_GRID[selected_index].ridge_lambda,
        ARB.Standardizer(
            np.zeros(len(ARB.FEATURE_NAMES)),
            np.ones(len(ARB.FEATURE_NAMES)),
        ),
        np.zeros(len(ARB.FEATURE_NAMES)),
        0.0,
    )
    rotation_scores = np.zeros((3, len(ARB.FIXED_GRID)))
    rotation_scores[:, selected_index] = 1.0
    grid_scores = ARB.GridScores(
        ARB.DEVELOPMENT_CONTEXTS,
        ARB.FIXED_GRID,
        rotation_scores,
        np.full(3, selected_index, dtype=np.int64),
        np.ones(3, dtype=bool),
    )
    selection = ARB.LocoSelection(
        grid_scores,
        selected_index,
        ARB.FIXED_GRID[selected_index],
        model,
    )
    null = ARB.FullSelectionNull(
        selected_config_indices=np.zeros(100, dtype=np.int64),
        selected_configs=np.repeat(
            ARB.FIXED_GRID[selected_index].as_array()[None, :], 100, axis=0
        ),
        selected_macro_scores=np.ones(100),
        rotation_scores=np.repeat(rotation_scores[None, :, :], 100, axis=0),
        rotation_winner_indices=np.zeros((100, 3), dtype=np.int64),
        rotation_unique_winner=np.ones((100, 3), dtype=bool),
        common_rotation_winner=np.ones(100, dtype=bool),
        label_sha256=np.full((100, 3), "a" * 64, dtype="U64"),
        scaler_mean=np.zeros((100, len(ARB.FEATURE_NAMES))),
        scaler_scale=np.ones((100, len(ARB.FEATURE_NAMES))),
        coefficient=np.zeros((100, len(ARB.FEATURE_NAMES))),
        intercept=np.zeros(100),
    )
    return CONTRACT.selection_arrays(selection, null)


def selection_arrays_aligned_to_held() -> dict[str, np.ndarray]:
    arrays = fake_selection_arrays()
    logit_first = np.log(0.9 / 0.1)
    logit_second = np.log(0.2 / 0.8)
    coefficient = np.zeros(len(ARB.FEATURE_NAMES))
    coefficient[0] = (logit_second - logit_first) / 0.9
    arrays["coefficient"] = coefficient
    arrays["intercept"] = np.asarray([logit_first])
    arrays["null_coefficient"] = np.repeat(coefficient[None, :], 100, axis=0)
    arrays["null_intercept"] = np.full(100, logit_first)
    return arrays


def fake_held_arrays() -> dict[str, np.ndarray]:
    genes = np.asarray(["G0", "G1", "UNSUPPORTED", "FORBIDDEN"])
    targets = np.asarray(["P0", "P1"])
    incumbent = np.asarray(
        [
            [1.0, -0.5, 0.25, 0.0],
            [-0.8, 0.6, -0.2, 0.0],
        ],
        dtype=np.float64,
    )
    state_mean = np.asarray(
        [
            [-0.7, 0.4, -0.1, 0.0],
            [0.5, -0.3, 0.1, 0.0],
        ],
        dtype=np.float64,
    )
    offsets = np.asarray([-0.03, -0.01, 0.01, 0.03])[:, None, None]
    state_banks = state_mean[None, :, :] + offsets
    state_banks[:, :, 3] = 0.0
    # The symmetric offsets have exact zero mean in float64.
    state_mean = state_banks.mean(axis=0)
    supported = np.asarray([True, True, False, False])
    forbidden = np.asarray([False, False, False, True])
    geometry = ~forbidden
    rows = np.asarray([0, 1], dtype=np.int64)
    cols = np.asarray([0, 1], dtype=np.int64)
    features = np.arange(18, dtype=np.float64).reshape(2, 9) / 10
    probabilities = np.asarray([0.9, 0.2])
    config = ARB.FIXED_GRID[0]
    replacement = np.zeros_like(incumbent, dtype=bool)
    selected = probabilities >= config.probability_threshold
    replacement[rows[selected], cols[selected]] = True
    uncapped = ARB.sparse_state_replacement(
        incumbent, state_mean, replacement, config.state_blend_fraction
    )
    candidate, row_scales = ARB.apply_row_cosine_distance_cap(
        incumbent,
        uncapped,
        geometry,
        config.maximum_cosine_distance,
    )
    mutable = np.broadcast_to(supported & ~forbidden, incumbent.shape).copy()
    candidate = ARB.enforce_coordinate_contract(
        candidate, incumbent, mutable, forbidden
    )
    null_configs = np.repeat(config.as_array()[None, :], 100, axis=0)
    null_pair_selected = np.repeat(selected[None, :], 100, axis=0)
    null_scales = np.repeat(row_scales[None, :], 100, axis=0)
    return {
        "feature_names": np.asarray(ARB.FEATURE_NAMES, dtype="U"),
        "genes": genes,
        "targets": targets,
        "incumbent": incumbent,
        "state_banks": state_banks,
        "state_mean": state_mean,
        "state_supported": supported,
        "forbidden_gene_mask": forbidden,
        "geometry_mask": geometry,
        "pair_rows": rows,
        "pair_cols": cols,
        "features": features,
        "selected_probability": probabilities,
        "replacement_mask": replacement,
        "row_cap_scale": row_scales,
        "candidate": candidate,
        "selected_config": config.as_array(),
        "null_selected_configs": null_configs,
        "null_pair_selected": null_pair_selected,
        "null_row_cap_scale": null_scales,
    }


def fake_development_context_arrays(
    context_name: str = "rpe1",
) -> dict[str, np.ndarray]:
    held = fake_held_arrays()
    incumbent = held["incumbent"]
    state_banks = held["state_banks"]
    state_mean = held["state_mean"]
    genes = held["genes"]
    targets = held["targets"]
    reference = incumbent.copy()
    reference[0, 0] = state_mean[0, 0]
    reference[1, 1] = state_mean[1, 1]
    significant = np.ones_like(reference, dtype=bool)
    own_target = np.zeros_like(reference, dtype=bool)
    supported = held["state_supported"]
    forbidden = held["forbidden_gene_mask"]
    cross_donor = np.asarray([[0.8, -0.4, 0.2, 0.0], [-0.7, 0.5, -0.1, 0.0]])
    control_log_cpm = np.asarray([6.0, 5.0, 4.0, 3.0])
    rebuilt = ARB.make_development_context(
        name=context_name,
        incumbent_lfc=incumbent,
        state_bank_lfc=state_banks,
        a14_cross_donor_agreement=cross_donor,
        control_log_cpm=control_log_cpm,
        reference_lfc=reference,
        reference_significant=significant,
        state_supported_gene_mask=supported,
        forbidden_gene_mask=forbidden,
        own_target_coordinate_mask=own_target,
        target_blocks_by_row=targets,
    )
    return {
        "context": np.asarray([context_name]),
        "feature_names": np.asarray(ARB.FEATURE_NAMES, dtype="U"),
        "genes": genes,
        "targets": targets,
        "target_blocks_by_row": targets.copy(),
        "incumbent": incumbent,
        "state_banks": state_banks,
        "state_mean": state_mean,
        "a14_cross_donor_agreement": cross_donor,
        "control_log_cpm": control_log_cpm,
        "reference_lfc": reference,
        "reference_significant": significant,
        "own_target_coordinate_mask": own_target,
        "state_supported": supported,
        "forbidden_gene_mask": forbidden,
        "geometry_mask": ~forbidden,
        "direction_mask": rebuilt.direction_mask,
        "mutable_mask": rebuilt.mutable_mask,
        "pair_rows": rebuilt.pair_rows,
        "pair_cols": rebuilt.pair_cols,
        "pair_target_blocks": rebuilt.target_blocks,
        "features": rebuilt.features,
        "labels": rebuilt.labels,
    }


def base_manifest(stage: str, *, context_name: str = "rpe1") -> dict[str, object]:
    if stage == CONTRACT.DEVELOPMENT_SELECTION_STAGE:
        inputs = [
            {
                "path": f"development/{context}/manifest.json",
                "sha256": "1" * 64,
                "role": f"development-context-{context}",
            }
            for context in ARB.DEVELOPMENT_CONTEXTS
        ]
    elif stage == CONTRACT.DEVELOPMENT_CONTEXT_STAGE:
        inputs = [
            {
                "path": "checkpoints/rpe1-out/manifest.json",
                "sha256": "1" * 64,
                "role": "state-checkpoint-manifest",
            }
        ]
    else:
        inputs = [
            {
                "path": "selection/manifest.json",
                "sha256": "1" * 64,
                "role": "selection-freeze",
            }
        ]
    common: dict[str, object] = {
        "schema_version": 1,
        "experiment": CONTRACT.EXPERIMENT,
        "stage": stage,
        "artifact": CONTRACT.ARTIFACT_NAME,
        "artifact_sha256": "0" * 64,
        "protocol": "preregistered_protocol.json",
        "protocol_sha256": hashlib.sha256(
            (HERE / "preregistered_protocol.json").read_bytes()
        ).hexdigest(),
        "implementation_sha256": {
            "arbitration.py": hashlib.sha256(
                (HERE / "arbitration.py").read_bytes()
            ).hexdigest(),
            "artifact_contract.py": hashlib.sha256(
                (HERE / "artifact_contract.py").read_bytes()
            ).hexdigest(),
        },
        "input_artifacts": inputs,
        "held_truth_read": False,
    }
    if stage == CONTRACT.DEVELOPMENT_CONTEXT_STAGE:
        checkpoint_path = f"checkpoints/{context_name}-out/manifest.json"
        common["input_artifacts"] = [
            {
                "path": checkpoint_path,
                "sha256": "1" * 64,
                "role": "state-checkpoint-manifest",
            }
        ]
        common.update(
            {
                "context": context_name,
                "contains_source_reference_truth": True,
                "source_truth_scope": "public-development-only",
                "state_training_exclusion": context_name,
                "in_context_trained_state_forbidden": True,
                "control_banks": 4,
                "gene_axis_sha256": "2" * 64,
                "target_axis_sha256": "3" * 64,
                "state_checkpoint_manifest": checkpoint_path,
                "state_checkpoint_manifest_sha256": "1" * 64,
            }
        )
    elif stage == CONTRACT.DEVELOPMENT_SELECTION_STAGE:
        common.update(
            {
                "development_contexts": list(ARB.DEVELOPMENT_CONTEXTS),
                "feature_names": list(ARB.FEATURE_NAMES),
                "grid_sha256": ARB.sha256_array(ARB.fixed_grid_array()),
                "null_repetitions": 100,
                "checkpoint_policy": {
                    "name": "last.ckpt",
                    "global_step": 20_000,
                    "best_ckpt_forbidden": True,
                    "strict_context_holdout_required": True,
                },
            }
        )
    else:
        common.update(
            {
                "surface": "synthetic",
                "contains_reference_truth": False,
                "sealed_before_held_truth": True,
                "control_banks": 4,
                "gene_axis_sha256": "2" * 64,
                "target_axis_sha256": "3" * 64,
                "selection_manifest": "selection/manifest.json",
                "selection_manifest_sha256": "4" * 64,
            }
        )
    return common


def materialize_manifest_inputs(manifest: dict[str, object], root: Path) -> None:
    inputs = manifest["input_artifacts"]
    assert isinstance(inputs, list)
    for item in inputs:
        assert isinstance(item, dict)
        path = root / str(item["path"])
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_bytes(b"bound test input")
        item["sha256"] = CONTRACT.sha256_file(path)
        if item["role"] == "selection-freeze":
            manifest["selection_manifest_sha256"] = item["sha256"]
        if item["role"] == "state-checkpoint-manifest":
            manifest["state_checkpoint_manifest_sha256"] = item["sha256"]


def test_selection_schema_rejects_feature_or_grid_drift() -> None:
    arrays = fake_selection_arrays()
    CONTRACT.validate_selection_arrays(arrays)
    changed = dict(arrays)
    changed["feature_names"] = arrays["feature_names"].copy()
    changed["feature_names"][0] = "gene_identity"
    try:
        CONTRACT.validate_selection_arrays(changed)
    except CONTRACT.ArtifactContractError:
        pass
    else:
        raise AssertionError("forbidden feature-name drift was accepted")

    changed = dict(arrays)
    changed["grid"] = arrays["grid"].copy()
    changed["grid"][0, 0] = 99
    try:
        CONTRACT.validate_selection_arrays(changed)
    except CONTRACT.ArtifactContractError:
        pass
    else:
        raise AssertionError("candidate-grid drift was accepted")


def test_development_context_schema_recomputes_features_and_labels() -> None:
    arrays = fake_development_context_arrays()
    CONTRACT.validate_development_context_arrays(arrays)
    changed = dict(arrays)
    changed["features"] = arrays["features"].copy()
    changed["features"][0, 0] += 1e-8
    try:
        CONTRACT.validate_development_context_arrays(changed)
    except CONTRACT.ArtifactContractError:
        pass
    else:
        raise AssertionError("source feature drift was accepted")


def test_atomic_development_context_artifact_round_trip() -> None:
    arrays = fake_development_context_arrays()
    manifest = base_manifest(CONTRACT.DEVELOPMENT_CONTEXT_STAGE)
    with tempfile.TemporaryDirectory() as temporary:
        root = Path(temporary)
        materialize_manifest_inputs(manifest, root)
        output = root / "development/rpe1"
        CONTRACT.write_artifact(
            output,
            stage=CONTRACT.DEVELOPMENT_CONTEXT_STAGE,
            arrays=arrays,
            manifest=manifest,
            provenance_root=root,
        )
        loaded, loaded_manifest = CONTRACT.load_artifact(
            output,
            expected_stage=CONTRACT.DEVELOPMENT_CONTEXT_STAGE,
            provenance_root=root,
        )
        assert loaded_manifest["state_training_exclusion"] == "rpe1"
        np.testing.assert_array_equal(loaded["features"], arrays["features"])


def test_selection_schema_requires_one_identical_winner() -> None:
    arrays = fake_selection_arrays()
    arrays["rotation_winner_indices"] = np.asarray([0, 1, 0])
    try:
        CONTRACT.validate_selection_arrays(arrays)
    except CONTRACT.ArtifactContractError:
        pass
    else:
        raise AssertionError("non-identical rotation winners were accepted")


def test_held_schema_rejects_truth_extra_and_fallback_drift() -> None:
    arrays = fake_held_arrays()
    CONTRACT.validate_held_prediction_arrays(arrays)

    extra = dict(arrays)
    extra["reference_truth"] = np.zeros_like(arrays["incumbent"])
    try:
        CONTRACT.validate_held_prediction_arrays(extra)
    except CONTRACT.ArtifactContractError:
        pass
    else:
        raise AssertionError("held truth array was accepted")

    changed = dict(arrays)
    changed["candidate"] = arrays["candidate"].copy()
    changed["candidate"][0, 2] += 1e-6
    try:
        CONTRACT.validate_held_prediction_arrays(changed)
    except CONTRACT.ArtifactContractError:
        pass
    else:
        raise AssertionError("unsupported fallback drift was accepted")


def test_held_schema_rejects_forbidden_replacement_and_bad_bank_mean() -> None:
    arrays = fake_held_arrays()
    changed = dict(arrays)
    changed["replacement_mask"] = arrays["replacement_mask"].copy()
    changed["replacement_mask"][0, 3] = True
    try:
        CONTRACT.validate_held_prediction_arrays(changed)
    except CONTRACT.ArtifactContractError:
        pass
    else:
        raise AssertionError("forbidden replacement was accepted")

    changed = dict(arrays)
    changed["state_mean"] = arrays["state_mean"].copy()
    changed["state_mean"][0, 0] += 1e-5
    try:
        CONTRACT.validate_held_prediction_arrays(changed)
    except CONTRACT.ArtifactContractError:
        pass
    else:
        raise AssertionError("STATE bank/mean mismatch was accepted")


def test_held_masks_are_bound_to_frozen_observed_and_null_models() -> None:
    held = fake_held_arrays()
    selection = selection_arrays_aligned_to_held()
    CONTRACT.validate_held_against_selection(held, selection)
    changed = dict(held)
    changed["null_pair_selected"] = held["null_pair_selected"].copy()
    changed["null_pair_selected"][0, 0] = False
    try:
        CONTRACT.validate_held_against_selection(changed, selection)
    except CONTRACT.ArtifactContractError:
        pass
    else:
        raise AssertionError("held null mask drift was accepted")


def test_manifest_fails_closed_on_truth_or_checkpoint_policy() -> None:
    manifest = base_manifest(CONTRACT.HELD_PREDICTION_STAGE)
    CONTRACT.validate_manifest(
        manifest,
        stage=CONTRACT.HELD_PREDICTION_STAGE,
        artifact_sha256="0" * 64,
    )
    manifest["held_truth_read"] = True
    try:
        CONTRACT.validate_manifest(
            manifest,
            stage=CONTRACT.HELD_PREDICTION_STAGE,
            artifact_sha256="0" * 64,
        )
    except CONTRACT.ArtifactContractError:
        pass
    else:
        raise AssertionError("held-truth-read manifest was accepted")

    selection_manifest = base_manifest(CONTRACT.DEVELOPMENT_SELECTION_STAGE)
    selection_manifest["checkpoint_policy"] = {
        "name": "best.ckpt",
        "global_step": 20_000,
        "best_ckpt_forbidden": False,
        "strict_context_holdout_required": True,
    }
    try:
        CONTRACT.validate_manifest(
            selection_manifest,
            stage=CONTRACT.DEVELOPMENT_SELECTION_STAGE,
            artifact_sha256="0" * 64,
        )
    except CONTRACT.ArtifactContractError:
        pass
    else:
        raise AssertionError("best.ckpt policy was accepted")


def test_atomic_selection_artifact_round_trip() -> None:
    arrays = fake_selection_arrays()
    manifest = base_manifest(CONTRACT.DEVELOPMENT_SELECTION_STAGE)
    with tempfile.TemporaryDirectory() as temporary:
        root = Path(temporary)
        materialize_manifest_inputs(manifest, root)
        output = root / "selection-output"
        path = CONTRACT.write_artifact(
            output,
            stage=CONTRACT.DEVELOPMENT_SELECTION_STAGE,
            arrays=arrays,
            manifest=manifest,
            provenance_root=root,
        )
        assert path == (output / "manifest.json").resolve()
        loaded, loaded_manifest = CONTRACT.load_artifact(
            output,
            expected_stage=CONTRACT.DEVELOPMENT_SELECTION_STAGE,
            provenance_root=root,
        )
        assert loaded_manifest["artifact_sha256"] == CONTRACT.sha256_file(
            output / CONTRACT.ARTIFACT_NAME
        )
        assert set(loaded) == set(arrays)
        try:
            CONTRACT.write_artifact(
                output,
                stage=CONTRACT.DEVELOPMENT_SELECTION_STAGE,
                arrays=arrays,
                manifest=manifest,
                provenance_root=root,
            )
        except FileExistsError:
            pass
        else:
            raise AssertionError("artifact writer overwrote an existing freeze")


def test_atomic_held_artifact_round_trip_without_truth() -> None:
    arrays = fake_held_arrays()
    manifest = base_manifest(CONTRACT.HELD_PREDICTION_STAGE)
    with tempfile.TemporaryDirectory() as temporary:
        root = Path(temporary)
        materialize_manifest_inputs(manifest, root)
        output = root / "held"
        CONTRACT.write_artifact(
            output,
            stage=CONTRACT.HELD_PREDICTION_STAGE,
            arrays=arrays,
            manifest=manifest,
            provenance_root=root,
        )
        loaded, loaded_manifest = CONTRACT.load_artifact(
            output,
            expected_stage=CONTRACT.HELD_PREDICTION_STAGE,
            provenance_root=root,
        )
        assert loaded_manifest["contains_reference_truth"] is False
        assert loaded_manifest["sealed_before_held_truth"] is True
        assert not any(
            token in name
            for name in loaded
            for token in ("truth", "reference", "label", "significant")
        )


def test_selection_artifact_requires_three_valid_context_artifacts() -> None:
    with tempfile.TemporaryDirectory() as temporary:
        root = Path(temporary)
        context_dirs: dict[str, Path] = {}
        for context in ARB.DEVELOPMENT_CONTEXTS:
            manifest = base_manifest(
                CONTRACT.DEVELOPMENT_CONTEXT_STAGE,
                context_name=context,
            )
            materialize_manifest_inputs(manifest, root)
            directory = root / f"development/{context}"
            CONTRACT.write_artifact(
                directory,
                stage=CONTRACT.DEVELOPMENT_CONTEXT_STAGE,
                arrays=fake_development_context_arrays(context),
                manifest=manifest,
                provenance_root=root,
            )
            context_dirs[context] = directory

        selection_manifest = base_manifest(CONTRACT.DEVELOPMENT_SELECTION_STAGE)
        for item in selection_manifest["input_artifacts"]:
            path = root / str(item["path"])
            item["sha256"] = CONTRACT.sha256_file(path)
        selection_dir = root / "selection"
        CONTRACT.write_artifact(
            selection_dir,
            stage=CONTRACT.DEVELOPMENT_SELECTION_STAGE,
            arrays=fake_selection_arrays(),
            manifest=selection_manifest,
            provenance_root=root,
        )
        _, _, loaded_contexts = CONTRACT.load_selection_with_contexts(
            selection_dir,
            context_dirs,
            provenance_root=root,
        )
        assert set(loaded_contexts) == set(ARB.DEVELOPMENT_CONTEXTS)


def test_linked_held_and_selection_artifacts_cross_validate() -> None:
    with tempfile.TemporaryDirectory() as temporary:
        root = Path(temporary)
        selection_manifest = base_manifest(CONTRACT.DEVELOPMENT_SELECTION_STAGE)
        materialize_manifest_inputs(selection_manifest, root)
        selection_dir = root / "selection"
        CONTRACT.write_artifact(
            selection_dir,
            stage=CONTRACT.DEVELOPMENT_SELECTION_STAGE,
            arrays=selection_arrays_aligned_to_held(),
            manifest=selection_manifest,
            provenance_root=root,
        )

        held_manifest = base_manifest(CONTRACT.HELD_PREDICTION_STAGE)
        held_manifest["selection_manifest"] = "selection/manifest.json"
        held_manifest["input_artifacts"] = [
            {
                "path": "selection/manifest.json",
                "sha256": CONTRACT.sha256_file(selection_dir / "manifest.json"),
                "role": "selection-freeze",
            }
        ]
        held_manifest["selection_manifest_sha256"] = held_manifest["input_artifacts"][
            0
        ]["sha256"]
        held_dir = root / "held"
        CONTRACT.write_artifact(
            held_dir,
            stage=CONTRACT.HELD_PREDICTION_STAGE,
            arrays=fake_held_arrays(),
            manifest=held_manifest,
            provenance_root=root,
        )
        held, _, selection, _ = CONTRACT.load_held_with_selection(
            held_dir,
            selection_dir,
            provenance_root=root,
        )
        CONTRACT.validate_held_against_selection(held, selection)


def main() -> None:
    tests = [
        value
        for name, value in sorted(globals().items())
        if name.startswith("test_") and callable(value)
    ]
    for test in tests:
        test()
        print(f"PASS {test.__name__}")
    print(f"{len(tests)} tests passed")


if __name__ == "__main__":
    main()
