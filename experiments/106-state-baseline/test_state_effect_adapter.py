#!/usr/bin/env python3
"""Synthetic and local-axis tests for the STATE effect adapter."""

from __future__ import annotations

import importlib.util
import inspect
import sys
import tempfile
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import scipy.sparse as sp


HERE = Path(__file__).resolve().parent
SPEC = importlib.util.spec_from_file_location(
    "state_effect_adapter", HERE / "state_effect_adapter.py"
)
if SPEC is None or SPEC.loader is None:
    raise ImportError(HERE / "state_effect_adapter.py")
ADAPTER = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = ADAPTER
SPEC.loader.exec_module(ADAPTER)


def test_checkpoint_contract_rejects_wrong_shape() -> None:
    genes = tuple(f"G{i}" for i in range(ADAPTER.EXPECTED_INPUT_DIM))
    good = SimpleNamespace(
        input_dim=ADAPTER.EXPECTED_INPUT_DIM,
        output_dim=ADAPTER.EXPECTED_INPUT_DIM,
        pert_dim=ADAPTER.EXPECTED_PERT_DIM,
        cell_sentence_len=ADAPTER.SET_LEN,
        hidden_dim=ADAPTER.EXPECTED_HIDDEN_DIM,
        output_space="all",
        predict_residual=True,
        batch_encoder=None,
        gene_decoder=None,
        hparams={"gene_names": genes},
    )
    ADAPTER.validate_model_contract(good, genes)
    bad = SimpleNamespace(**good.__dict__)
    bad.output_dim -= 1
    try:
        ADAPTER.validate_model_contract(bad, genes)
    except ValueError:
        pass
    else:
        raise AssertionError("wrong model shape was accepted")


def test_frozen_runtime_artifacts_are_rehashed() -> None:
    with tempfile.TemporaryDirectory(dir=HERE) as temp:
        run_dir = Path(temp)
        names = (
            "config.yaml",
            "hparams.yaml",
            "var_dims.pkl",
            "pert_onehot_map.pt",
            "batch_onehot_map.torch",
            "cell_type_onehot_map.torch",
            "uv.lock",
            "starter.local.toml",
        )
        hashes = {}
        for name in names:
            path = run_dir / name
            path.write_bytes(name.encode())
            hashes[name] = ADAPTER.sha256_file(path)
        checkpoint_path = run_dir / "checkpoint.ckpt"
        checkpoint_path.write_bytes(b"checkpoint")
        source_tree = run_dir / "source/state"
        source_tree.mkdir(parents=True)
        (source_tree / "__init__.py").write_text("VALUE = 1\n")
        packages = ["package==1"]
        packages_path = run_dir / "runtime-packages.txt"
        packages_path.write_text("\n".join(packages) + "\n")
        runtime = {
            "python_executable": "/test/python",
            "python_version": "3.test",
            "torch_version": "test",
            "transformers_version": "test",
            "packages_sha256": ADAPTER.sha256_file(packages_path),
        }
        spec = ADAPTER.CheckpointSpec(
            path=checkpoint_path,
            sha256=ADAPTER.sha256_file(checkpoint_path),
            global_step=ADAPTER.EXPECTED_STEP,
            profile=ADAPTER.EXPECTED_PROFILE,
            source_commit=ADAPTER.EXPECTED_SOURCE_COMMIT,
            source_patch_sha256="patch",
            source_tree="source/state",
            source_tree_sha256=ADAPTER.tree_sha256(source_tree),
            config_sha256=hashes["config.yaml"],
            hparams_sha256=hashes["hparams.yaml"],
            var_dims_sha256=hashes["var_dims.pkl"],
            pert_map_sha256=hashes["pert_onehot_map.pt"],
            batch_map_sha256=hashes["batch_onehot_map.torch"],
            cell_type_map_sha256=hashes["cell_type_onehot_map.torch"],
            uv_lock_sha256=hashes["uv.lock"],
            support_archive_sha256="archive",
            state_gene_axis_sha256="axis",
            starter_toml_sha256=hashes["starter.local.toml"],
            runtime_packages_sha256=ADAPTER.sha256_file(packages_path),
            runtime=runtime,
            training_recipe={"max_steps": ADAPTER.EXPECTED_STEP, "train_seed": 42},
            support_files={"competition_train.h5": {"sha256": "frozen"}},
        )
        original_axis_hash = ADAPTER.sha256_strings
        original_runtime = ADAPTER.runtime_environment
        try:
            ADAPTER.sha256_strings = lambda values: spec.state_gene_axis_sha256
            ADAPTER.runtime_environment = lambda: runtime
            ADAPTER.validate_runtime_provenance(spec, run_dir, ("A",))
            (run_dir / "config.yaml").write_bytes(b"changed")
            try:
                ADAPTER.validate_runtime_provenance(spec, run_dir, ("A",))
            except ValueError:
                pass
            else:
                raise AssertionError("changed frozen runtime artifact was accepted")
        finally:
            ADAPTER.sha256_strings = original_axis_hash
            ADAPTER.runtime_environment = original_runtime


def test_incumbent_manifest_binds_surface_controls_and_axes() -> None:
    import json

    with tempfile.TemporaryDirectory(dir=HERE) as temp:
        root = Path(temp)
        incumbent = root / "frozen.npz"
        incumbent.write_bytes(b"incumbent")
        controls = root / "controls.h5ad"
        controls.write_bytes(b"controls")
        controls_manifest = root / "controls.h5ad.manifest.json"
        controls_manifest.write_text(
            json.dumps(
                {
                    "stage": "held-control-only-input",
                    "surface": "hct116",
                    "contains_perturbation_rows": False,
                    "control_rows": ADAPTER.EXPECTED_CONTROL_ROWS,
                    "expression_read_policy": (
                        "only predicate-filtered control rows are materialized "
                        "by experiment code"
                    ),
                    "accepted_control_label_counts": {"Non-Targeting": 165_777},
                }
            )
        )
        operator = root / "operator.npz"
        operator_metadata = operator.with_suffix(".json")
        operator.write_bytes(b"operator")
        operator_metadata.write_bytes(b"metadata")
        original_operator_paths = ADAPTER.HELD_OPERATOR_PATHS
        ADAPTER.HELD_OPERATOR_PATHS = {
            **original_operator_paths,
            "hct116": operator,
        }
        manifest = root / "manifest.json"
        manifest.write_text(
            json.dumps(
                {
                    "experiment": "106-state-baseline",
                    "stage": "held-incumbent-freeze",
                    "surface": "hct116",
                    "sealed_before_perturbation_truth": True,
                    "endpoint": "a14-shaped",
                    "alpha": 0.982,
                    "cd4_fraction": 0.05,
                    "basis_prune": "essential",
                    "frozen_sha256": ADAPTER.sha256_file(incumbent),
                    "controls_h5ad_sha256": ADAPTER.sha256_file(controls),
                    "controls_manifest_sha256": ADAPTER.sha256_file(controls_manifest),
                    "held_operator": str(operator.relative_to(ADAPTER.ROOT)),
                    "held_operator_sha256": ADAPTER.sha256_file(operator),
                    "held_operator_metadata_sha256": ADAPTER.sha256_file(
                        operator_metadata
                    ),
                    "target_order_sha256": "targets",
                    "gene_axis_sha256": "genes",
                }
            )
        )
        try:
            spec = ADAPTER.IncumbentSpec.load(manifest, incumbent, controls, "hct116")
            assert spec.surface == "hct116"
            controls.write_bytes(b"changed")
            try:
                ADAPTER.IncumbentSpec.load(manifest, incumbent, controls, "hct116")
            except ValueError:
                pass
            else:
                raise AssertionError("changed controls were accepted")
        finally:
            ADAPTER.HELD_OPERATOR_PATHS = original_operator_paths


def test_h1_axis_is_exact_identity() -> None:
    genes = ADAPTER.read_axis(ADAPTER.STATE_GENE_PATH)
    alignment = ADAPTER.GeneAlignment.build(genes, genes)
    np.testing.assert_array_equal(alignment.source_cols, np.arange(len(genes)))
    np.testing.assert_array_equal(alignment.state_cols, np.arange(len(genes)))
    assert alignment.source_only == ()
    assert alignment.state_only == ()


def test_2026_alignment_is_18077_456_3() -> None:
    source = ADAPTER.read_axis(
        ADAPTER.ROOT / "data/challenge/gene_names.csv", "gene_name"
    )
    state = ADAPTER.read_axis(ADAPTER.STATE_GENE_PATH)
    alignment = ADAPTER.GeneAlignment.build(source, state)
    assert len(alignment.source_cols) == 18_077
    assert len(alignment.source_only) == 456
    assert set(alignment.state_only) == {"HSPA14-1", "TBCE-1", "TMSB15B-1"}


def test_duplicate_aliases_are_not_merged() -> None:
    alignment = ADAPTER.GeneAlignment.build(
        ("HSPA14", "TBCE", "X"), ("HSPA14-1", "TBCE-1", "X")
    )
    assert alignment.source_cols.tolist() == [2]
    assert alignment.state_cols.tolist() == [2]
    assert alignment.source_only == ("HSPA14", "TBCE")
    assert alignment.state_only == ("HSPA14-1", "TBCE-1")


def test_missing_or_zero_target_embedding_fails_closed() -> None:
    good = np.ones(ADAPTER.EXPECTED_PERT_DIM, dtype=np.float32)
    try:
        ADAPTER.load_target_embeddings(("MISSING",), {"GOOD": good})
    except KeyError:
        pass
    else:
        raise AssertionError("missing target embedding was accepted")
    try:
        ADAPTER.load_target_embeddings(("ZERO",), {"ZERO": np.zeros_like(good)})
    except ValueError:
        pass
    else:
        raise AssertionError("zero target embedding was accepted")


def test_target_order_matches_allowlist() -> None:
    a = np.full(ADAPTER.EXPECTED_PERT_DIM, 1, dtype=np.float32)
    b = np.full(ADAPTER.EXPECTED_PERT_DIM, 2, dtype=np.float32)
    observed = ADAPTER.load_target_embeddings(("B", "A"), {"A": a, "B": b})
    assert observed[:, 0].tolist() == [2.0, 1.0]


def test_control_sampling_is_seed_stable_and_reused_across_targets() -> None:
    rows = np.arange(1_000)
    first = ADAPTER.sample_control_sets(rows, "HCT116")
    second = ADAPTER.sample_control_sets(rows, "HCT116")
    other = ADAPTER.sample_control_sets(rows, "HEK293T")
    np.testing.assert_array_equal(first, second)
    assert not np.array_equal(first, other)
    assert first.shape == (ADAPTER.N_CONTROL_SETS, ADAPTER.SET_LEN)


def test_each_transformer_sequence_is_target_homogeneous() -> None:
    embedding = np.arange(ADAPTER.EXPECTED_PERT_DIM, dtype=np.float32)
    batch = ADAPTER.build_perturbation_batch(embedding, 4)
    reshaped = batch.reshape(4, ADAPTER.SET_LEN, -1)
    for block in reshaped:
        np.testing.assert_array_equal(block, np.broadcast_to(embedding, block.shape))


def test_streamed_accumulation_matches_materialized_accumulation() -> None:
    rng = np.random.default_rng(4)
    values = rng.uniform(0, 5, size=(512, 17)).astype(np.float32)
    streamed = ADAPTER.streamed_mean_counts(
        [values[:111], values[111:300], values[300:]]
    )
    materialized = np.expm1(values.astype(np.float64)).mean(axis=0)
    np.testing.assert_allclose(streamed, materialized, rtol=1e-12, atol=1e-12)


def test_streamed_accumulation_clips_like_upstream_inference() -> None:
    values = np.asarray([[-2.0, 3.0, 17.0], [1.0, 4.0, 14.0]])
    observed = ADAPTER.streamed_mean_counts([values])
    expected = np.expm1(np.clip(values, 0.0, ADAPTER.MAX_LOG_EXPRESSION)).mean(axis=0)
    np.testing.assert_allclose(observed, expected, rtol=1e-12, atol=1e-12)


def test_identity_model_produces_exact_zero_effect() -> None:
    counts = np.asarray([1.0, 2.0, 0.0, 7.0])
    effect = ADAPTER.scorer_effect(counts, counts)
    np.testing.assert_array_equal(effect, np.zeros_like(effect))


def test_challenge_only_genes_preserve_control_expression() -> None:
    alignment = ADAPTER.GeneAlignment.build(("A", "B", "C"), ("A", "C", "D"))
    state_mean = np.asarray([10.0, 30.0, 5.0])
    control = np.asarray([2.0, 7.0, 8.0])
    projected = alignment.project_mean_to_source(state_mean, control)
    assert projected[1] == control[1]
    effect = ADAPTER.scorer_effect(projected, control)
    assert effect[1] == 0.0


def test_log_pseudobulk_is_not_mean_log_expression() -> None:
    counts = np.asarray([[1.0, 9.0], [9.0, 1.0]])
    correct = ADAPTER.bulk_profile(counts.mean(axis=0))
    wrong = np.log1p(counts).mean(axis=0)
    assert not np.allclose(correct, wrong)


def test_contexts_are_never_pooled() -> None:
    rows = np.arange(1_000)
    one = ADAPTER.sample_control_sets(rows, "HCT116")
    two = ADAPTER.sample_control_sets(rows, "HEK293T")
    assert not np.array_equal(one, two)


def test_panel_coordinates_are_removed_from_full_pds_gallery() -> None:
    panel = ADAPTER.read_axis(ADAPTER.PANEL_PATH, "target_gene")
    genes = (panel[0], "SAFE", panel[1])
    effect = np.ones((2, 3))
    masked = ADAPTER.apply_coordinate_mask(effect, genes, set(panel))
    np.testing.assert_array_equal(masked[:, [0, 2]], 0)
    np.testing.assert_array_equal(masked[:, 1], 1)


def test_sparse_control_alignment_uses_log1p_raw_counts() -> None:
    alignment = ADAPTER.GeneAlignment.build(("B", "A"), ("A", "C", "B"))
    raw = sp.csr_matrix(np.asarray([[3.0, 8.0], [0.0, 2.0]]))
    observed = alignment.align_and_log_controls(raw)
    expected = np.asarray(
        [[np.log1p(8.0), 0.0, np.log1p(3.0)], [np.log1p(2.0), 0.0, 0.0]],
        dtype=np.float32,
    )
    np.testing.assert_allclose(observed, expected)


def test_inference_api_has_no_perturbed_truth_argument() -> None:
    parameters = inspect.signature(ADAPTER.infer_effects).parameters
    forbidden = {"truth", "reference", "perturbed", "real_effect"}
    assert forbidden.isdisjoint(parameters)


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
