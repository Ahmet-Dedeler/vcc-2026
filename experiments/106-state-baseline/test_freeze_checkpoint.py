#!/usr/bin/env python3
"""Pure contract tests for fixed-step STATE checkpoint freezing."""

from __future__ import annotations

import importlib.util
import sys
import tempfile
from pathlib import Path
from types import SimpleNamespace


HERE = Path(__file__).resolve().parent
SPEC = importlib.util.spec_from_file_location(
    "freeze_checkpoint", HERE / "freeze_checkpoint.py"
)
if SPEC is None or SPEC.loader is None:
    raise ImportError(HERE / "freeze_checkpoint.py")
FREEZE = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = FREEZE
SPEC.loader.exec_module(FREEZE)


def valid_metadata(genes: tuple[str, ...]) -> dict[str, object]:
    return {
        "global_step": 20_000,
        "input_dim": 18_080,
        "output_dim": 18_080,
        "pert_dim": 5_120,
        "hidden_dim": 672,
        "cell_set_len": 128,
        "output_space": "all",
        "predict_residual": True,
        "gene_decoder_bool": False,
        "batch_encoder": False,
        "gene_names": genes,
    }


def test_checkpoint_contract_accepts_exact_fixed_step() -> None:
    genes = ("A", "B")
    FREEZE.validate_checkpoint_metadata(valid_metadata(genes), genes)


def test_checkpoint_contract_rejects_best_step() -> None:
    genes = ("A", "B")
    metadata = valid_metadata(genes)
    metadata["global_step"] = 18_000
    try:
        FREEZE.validate_checkpoint_metadata(metadata, genes)
    except ValueError:
        pass
    else:
        raise AssertionError("non-20k checkpoint was accepted")


def test_checkpoint_contract_rejects_axis_drift() -> None:
    genes = ("A", "B")
    try:
        FREEZE.validate_checkpoint_metadata(valid_metadata(genes), ("B", "A"))
    except ValueError:
        pass
    else:
        raise AssertionError("gene-axis drift was accepted")


def test_training_detection_accepts_split_output_dir_and_name() -> None:
    original_run = FREEZE.subprocess.run
    try:
        run_dir = Path("/tmp/runs/official-state_sm")
        resolved = run_dir.resolve()
        command = (
            "27064 state tx train "
            f"output_dir={resolved.parent} name={resolved.name} model=state_sm"
        )
        FREEZE.subprocess.run = lambda *args, **kwargs: SimpleNamespace(
            stdout=command + "\n"
        )
        observed = FREEZE.state_training_processes(run_dir)
        assert observed == [command]
    finally:
        FREEZE.subprocess.run = original_run


def test_tree_hash_changes_with_source_content() -> None:
    import tempfile

    with tempfile.TemporaryDirectory() as temp:
        root = Path(temp)
        (root / "state").mkdir()
        source = root / "state/model.py"
        source.write_text("x = 1\n")
        first = FREEZE.tree_sha256(root)
        source.write_text("x = 2\n")
        assert FREEZE.tree_sha256(root) != first


def test_tree_hash_ignores_bytecode_caches_excluded_from_snapshot() -> None:
    with tempfile.TemporaryDirectory() as temp:
        root = Path(temp)
        (root / "state").mkdir()
        (root / "state/model.py").write_text("x = 1\n")
        first = FREEZE.tree_sha256(root)
        cache = root / "state/__pycache__"
        cache.mkdir()
        (cache / "model.cpython-314.pyc").write_bytes(b"cache")
        assert FREEZE.tree_sha256(root) == first


def test_training_recipe_binds_seed_split_and_paths() -> None:
    with tempfile.TemporaryDirectory() as temp:
        root = Path(temp)
        starter = root / "starter.local.toml"
        support = FREEZE.SUPPORT_ROOT
        starter.write_text(
            "\n".join(
                (
                    "[datasets]",
                    f'replogle_h1 = "{support}/'
                    '{competition_train,k562_gwps,rpe1,jurkat,k562,hepg2}.h5"',
                    "[training]",
                    'replogle_h1 = "train"',
                    "[zeroshot]",
                    '"replogle_h1.hepg2" = "test"',
                    "[fewshot]",
                    "",
                )
            )
        )
        config = root / "config.yaml"
        config.write_text(
            f"""
data:
  name: PerturbationDataModule
  kwargs:
    toml_config_path: {starter.resolve()}
    num_workers: 8
    pin_memory: false
    batch_col: batch_var
    pert_col: target_gene
    cell_type_key: cell_type
    control_pert: non-targeting
    perturbation_features_file: {support / "ESM2_pert_features.pt"}
training:
  batch_size: 16
  lr: 0.0001
  max_steps: 20000
  train_seed: 42
  val_freq: 2000
  ckpt_every_n_steps: 20000
name: official-state_sm
"""
        )
        observed = FREEZE.validate_training_recipe(config, starter)
        assert observed["train_seed"] == 42
        assert observed["zeroshot_context"] == "hepg2"


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
