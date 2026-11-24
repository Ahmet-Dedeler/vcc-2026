#!/usr/bin/env python3
"""Pure tests for the held-context control-only reservoir."""

from __future__ import annotations

import importlib.util
import sys
import tempfile
from pathlib import Path

import numpy as np


HERE = Path(__file__).resolve().parent
SPEC = importlib.util.spec_from_file_location(
    "prepare_held_controls", HERE / "prepare_held_controls.py"
)
if SPEC is None or SPEC.loader is None:
    raise ImportError(HERE / "prepare_held_controls.py")
PREPARE = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = PREPARE
SPEC.loader.exec_module(PREPARE)


def test_control_labels_are_recognized_without_panel_names() -> None:
    for label in ("non-targeting", "NTC-12", "control", "safe-harbor", "neg_1"):
        assert PREPARE.is_control(label)
    for label in ("ADNP", "MYRF", "ZNF714", "NEGR1"):
        assert not PREPARE.is_control(label)


def test_reservoir_is_deterministic() -> None:
    def run() -> list[int]:
        rng = np.random.default_rng(7)
        slots = []
        for seen in range(1, 101):
            slot = PREPARE.reservoir_slot(seen, 10, rng)
            if slot is not None:
                slots.append(slot)
        return slots

    assert run() == run()


def test_sparse_mapping_drops_nonofficial_tokens() -> None:
    mapping = np.asarray([2, -1, 0], dtype=np.int64)
    columns, values = PREPARE.map_sparse_row([0, 1, 2], [4.0, 5.0, 6.0], mapping)
    np.testing.assert_array_equal(columns, [2, 0])
    np.testing.assert_array_equal(values, [4.0, 6.0])


def test_pyarrow_column_projection_uses_supported_list_type() -> None:
    import pyarrow as pa
    import pyarrow.parquet as pq

    with tempfile.TemporaryDirectory() as temp:
        path = Path(temp) / "rows.parquet"
        pq.write_table(
            pa.table(
                {
                    "gene_token_id": [[1]],
                    "gene_expression": [[2.0]],
                    "gene_target": ["non-targeting"],
                    "pass_guide_filter": [1],
                }
            ),
            path,
        )
        columns = [
            "gene_token_id",
            "gene_expression",
            "gene_target",
            "pass_guide_filter",
        ]
        assert pq.read_table(path, columns=columns).column_names == columns


def test_metadata_preflight_never_projects_expression() -> None:
    import pyarrow as pa

    calls = []

    class FakeParquet:
        num_row_groups = 1

        def read_row_group(self, row_group, columns):
            calls.append((row_group, tuple(columns)))
            assert "gene_token_id" not in columns
            assert "gene_expression" not in columns
            return pa.table(
                {
                    "gene_target": ["Non-Targeting", "NEGR1", "ADNP"],
                    "pass_guide_filter": [1, 1, 1],
                }
            )

    original = PREPARE.pq.ParquetFile
    PREPARE.pq.ParquetFile = lambda path: FakeParquet()
    try:
        layout = PREPARE.inspect_control_layout(
            Path("synthetic.parquet"),
            allowed_labels=("Non-Targeting",),
        )
    finally:
        PREPARE.pq.ParquetFile = original
    assert calls == [(0, ("gene_target", "pass_guide_filter"))]
    assert layout.accepted_label_counts == {"Non-Targeting": 1}
    assert layout.eligible_controls == 1
    assert layout.mixed_row_groups == 1
    np.testing.assert_array_equal(layout.eligible_source_rows, [0])


def test_expression_projection_uses_exact_control_predicates() -> None:
    import pyarrow as pa

    calls = []

    def fake_read_table(path, *, columns, filters, use_threads):
        calls.append((path, tuple(columns), filters, use_threads))
        return pa.table(
            {
                "gene_token_id": [[2]],
                "gene_expression": [[4.0]],
                "gene_target": ["Non-Targeting"],
                "pass_guide_filter": [1],
            }
        )

    layout = PREPARE.FileControlLayout(
        path=Path("synthetic.parquet"),
        row_groups=1,
        mixed_row_groups=1,
        all_control_row_groups=0,
        eligible_controls=1,
        accepted_label_counts={"Non-Targeting": 1},
        eligible_source_rows=np.asarray([7]),
    )
    original = PREPARE.pq.read_table
    PREPARE.pq.read_table = fake_read_table
    try:
        tokens, expression = PREPARE.read_predicate_filtered_controls(
            layout,
            accepted_labels=("Non-Targeting",),
        )
    finally:
        PREPARE.pq.read_table = original
    assert calls == [
        (
            Path("synthetic.parquet"),
            (
                "gene_token_id",
                "gene_expression",
                "gene_target",
                "pass_guide_filter",
            ),
            [
                ("pass_guide_filter", "=", 1),
                ("gene_target", "in", ("Non-Targeting",)),
            ],
            False,
        )
    ]
    assert tokens.tolist() == [[2]]
    assert expression.tolist() == [[4.0]]


def test_mixed_parquet_materializes_only_predicate_filtered_controls() -> None:
    import pyarrow as pa
    import pyarrow.parquet as pq

    with tempfile.TemporaryDirectory() as temp:
        path = Path(temp) / "mixed.parquet"
        pq.write_table(
            pa.table(
                {
                    "gene_token_id": [[1], [2], [3], [4]],
                    "gene_expression": [[10.0], [20.0], [30.0], [40.0]],
                    "gene_target": [
                        "Non-Targeting",
                        "NEGR1",
                        "Non-Targeting",
                        "ADNP",
                    ],
                    "pass_guide_filter": [1, 1, 0, 1],
                }
            ),
            path,
        )
        layout = PREPARE.inspect_control_layout(
            path,
            allowed_labels=("Non-Targeting",),
        )
        tokens, expression = PREPARE.read_predicate_filtered_controls(
            layout,
            accepted_labels=("Non-Targeting",),
        )
    assert layout.accepted_label_counts == {"Non-Targeting": 1}
    assert layout.mixed_row_groups == 1
    np.testing.assert_array_equal(layout.eligible_source_rows, [0])
    assert tokens.tolist() == [[1]]
    assert expression.tolist() == [[10.0]]


def test_actual_orion_layout_requires_predicate_filtering() -> None:
    for line in ("HCT116", "HEK293T"):
        files = sorted(
            PREPARE.DEFAULT_DATA.glob(f"{line}_Batch*.parquet"),
            key=PREPARE.batch_number,
        )
        preflight = PREPARE.preflight_control_layout(line, files)
        expected = PREPARE.EXPECTED_LAYOUT[line]
        assert preflight.accepted_label_counts == {
            "Non-Targeting": expected["eligible_controls"]
        }
        assert preflight.row_groups == expected["row_groups"]
        assert preflight.mixed_row_groups == preflight.row_groups
        assert preflight.all_control_row_groups == 0


def test_partial_control_artifact_is_removed_before_restart() -> None:
    with tempfile.TemporaryDirectory() as temp:
        output = Path(temp) / "controls.h5ad"
        output.write_bytes(b"partial")
        observed = PREPARE.validate_existing(
            output,
            line="HCT116",
            n_control=PREPARE.DEFAULT_N_CONTROL,
            seed=PREPARE.DEFAULT_SEED,
        )
        assert observed is None
        assert not output.exists()
        assert not output.with_suffix(".h5ad.manifest.json").exists()
        assert callable(PREPARE.shutil.rmtree)


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
