#!/usr/bin/env python3
"""Source-only tests for Experiment 107 strict artifact construction."""

from __future__ import annotations

import importlib.util
import csv
import hashlib
import sys
import unittest
from pathlib import Path
from unittest import mock

import numpy as np


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]


def load_module(name: str, path: Path):
    spec = importlib.util.spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise ImportError(path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


BUILDER = load_module(
    "exp107_test_strict_artifact_builder", HERE / "strict_artifact_builder.py"
)


class StrictArtifactBuilderTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.plans, cls.contract, cls.source = BUILDER.build_surface_plans()

    def test_metadata_plan_never_opens_or_hashes_nadig_delta(self) -> None:
        forbidden = (ROOT / "data/donors/nadig_jurkat/delta.npy").resolve()
        real_open = Path.open
        real_load = np.load

        def guarded_open(path: Path, *args: object, **kwargs: object):
            if Path(path).resolve() == forbidden:
                raise AssertionError("Nadig delta was opened")
            return real_open(path, *args, **kwargs)

        def guarded_load(path: object, *args: object, **kwargs: object):
            if isinstance(path, (str, Path)) and Path(path).resolve() == forbidden:
                raise AssertionError("Nadig delta was loaded")
            return real_load(path, *args, **kwargs)

        with (
            mock.patch.object(Path, "open", guarded_open),
            mock.patch.object(BUILDER.np, "load", side_effect=guarded_load),
        ):
            report = BUILDER.metadata_report()
        self.assertFalse(report["nadig_outcome_opened_or_hashed"])
        self.assertEqual(len(report["surfaces"]), 8)

    def test_all_frozen_surface_geometry_matches(self) -> None:
        expected = {
            "nadig.gallery_1": (300, 230, 17661, 29295, 298),
            "nadig.gallery_2": (300, 226, 17675, 29282, 295),
            "nadig.gallery_3": (300, 224, 17668, 29293, 298),
            "nadig.gallery_4": (300, 229, 17670, 29285, 294),
            "nadig.gallery_5": (300, 226, 17668, 29298, 293),
            "nadig.gallery_6": (300, 228, 17665, 29276, 287),
            "h1": (85, 343, 17862, 29244, 83),
            "orion_hek293t": (132, 349, 17942, 18370, 131),
        }
        observed = {
            key: (
                len(plan.targets),
                len(plan.score_genes),
                len(plan.operator_columns),
                sum(len(rows) for rows in plan.basis_rows.values()),
                int(plan.source_transferable.sum()),
            )
            for key, plan in self.plans.items()
        }
        self.assertEqual(observed, expected)

    def test_nonexportable_mrpl36_is_never_scored(self) -> None:
        self.assertIn("MRPL36", set(self.source.response_symbols))
        self.assertTrue(
            all("MRPL36" not in plan.score_genes for plan in self.plans.values())
        )

    def test_gallery_source_coverage_anchors(self) -> None:
        combined = {
            target
            for number in range(1, 7)
            for target in self.plans[f"nadig.gallery_{number}"].targets
        }
        pm = set(
            (ROOT / "data/external/perturbmatch/perts.txt").read_text().splitlines()
        )
        k562 = set((ROOT / "data/donors/k562_gwps/perts.txt").read_text().splitlines())
        cd4 = [
            set(
                (
                    ROOT
                    / f"data/external/cd4t/cloud-compact/cd4_compact/{condition}/perts.txt"
                )
                .read_text()
                .splitlines()
            )
            for condition in ("cd4_rest", "cd4_stim8hr", "cd4_stim48hr")
        ]
        self.assertEqual(len(combined & pm), 1734)
        self.assertEqual(len(combined & pm & k562), 1724)
        self.assertEqual(len(combined & set.intersection(*cd4)), 0)

    def test_vcc_only_operator_plan_reproduces_production_pruned_geometry(self) -> None:
        def read_csv_column(path: Path, column: str) -> tuple[str, ...]:
            with path.open(newline="", encoding="utf-8") as handle:
                return tuple(row[column] for row in csv.DictReader(handle))

        official = read_csv_column(ROOT / "data/challenge/gene_names.csv", "gene_name")
        panel = read_csv_column(ROOT / "data/challenge/pert_counts.csv", "target_gene")
        donor_genes = {
            name: (ROOT / f"data/donors/{name}/genes.txt").read_text().splitlines()
            for name in BUILDER.STRICT.DONOR_ORDER
        }
        donor_perts = {
            name: (ROOT / f"data/donors/{name}/perts.txt").read_text().splitlines()
            for name in BUILDER.STRICT.DONOR_ORDER
        }
        essential = set(
            (ROOT / "data/donors/k562_essential/perts.txt").read_text().splitlines()
        )
        columns = BUILDER.STRICT.whitener_columns(official, donor_genes, panel, ())
        rows = BUILDER.STRICT.select_whitener_basis_rows(donor_perts, essential, set())
        self.assertEqual(len(columns), 17942)
        self.assertEqual(
            hashlib.sha256(columns.tobytes()).hexdigest(),
            "b2cae0d381c3cd7065a3de33f1a6da5514577f2b3187b6dc84786ee789e27e4d",
        )
        self.assertEqual(
            {name: len(value) for name, value in rows.items()},
            {
                "k562_gwps": 7784,
                "orion_hct116": 10807,
                "orion_hek293t": 10803,
            },
        )

    def test_held_donor_sets_and_no_hct116_promotion_role_are_frozen(self) -> None:
        self.assertEqual(
            self.plans["h1"].donor_order,
            ("k562_gwps", "orion_hct116", "orion_hek293t"),
        )
        self.assertEqual(
            self.plans["orion_hek293t"].donor_order,
            ("k562_gwps", "orion_hct116"),
        )
        self.assertFalse(
            any(
                role.endswith(".orion_hct116") and not role.startswith("donor.")
                for role in BUILDER.STRICT.REQUIRED_ARTIFACT_ROLES
            )
        )
        self.assertIsNone(
            self.contract["held_context_gate"]["diagnostic_directional_surface"]
        )

    def test_scored_coordinate_missing_from_operator_is_retained_as_fallback(self) -> None:
        official = BUILDER._read_csv_column(  # noqa: SLF001
            ROOT / "data/challenge/gene_names.csv", "gene_name"
        )
        btbd11 = official.index("BTBD11")
        for surface in ("h1", "orion_hek293t"):
            plan = self.plans[surface]
            supported, positions, missing = BUILDER.partition_score_axis_by_operator(
                plan.score_indices_in_official,
                plan.operator_columns,
            )
            self.assertEqual(missing.tolist(), [btbd11])
            self.assertEqual(int(supported.sum()), len(plan.score_genes) - 1)
            self.assertEqual(len(positions), int(supported.sum()))
            self.assertEqual(
                plan.score_genes[
                    int(np.flatnonzero(~supported)[0])
                ],
                "BTBD11",
            )

    def test_score_operator_partition_rejects_duplicate_axes(self) -> None:
        with self.assertRaisesRegex(ValueError, "unique vectors"):
            BUILDER.partition_score_axis_by_operator(
                np.asarray([1, 1], dtype=np.int64),
                np.asarray([0, 1], dtype=np.int64),
            )
        with self.assertRaisesRegex(ValueError, "unique vectors"):
            BUILDER.partition_score_axis_by_operator(
                np.asarray([0, 1], dtype=np.int64),
                np.asarray([1, 1], dtype=np.int64),
            )

    def test_top_level_role_paths_are_pairwise_unique(self) -> None:
        paths = [
            path
            for role in BUILDER.STRICT.REQUIRED_ARTIFACT_ROLES
            for path in BUILDER._role_relpaths(role)  # noqa: SLF001
        ]
        self.assertEqual(len(paths), 2 * len(BUILDER.STRICT.REQUIRED_ARTIFACT_ROLES))
        self.assertEqual(len(paths), len(set(paths)))

    def test_frozen_transitive_input_sets_never_include_nadig_delta(self) -> None:
        forbidden = (ROOT / "data/donors/nadig_jurkat/delta.npy").resolve()
        for role in BUILDER.STRICT.REQUIRED_ARTIFACT_ROLES:
            inputs = {
                path.resolve()
                for path in BUILDER.STRICT.expected_transitive_paths(
                    role, self.contract
                )
            }
            self.assertNotIn(forbidden, inputs, role)
            if ".nadig.gallery_" in role:
                self.assertIn(
                    (ROOT / "data/donors/nadig_jurkat/perts.txt").resolve(),
                    inputs,
                )
                self.assertIn(
                    (ROOT / "data/donors/nadig_jurkat/genes.txt").resolve(),
                    inputs,
                )
        with self.assertRaisesRegex(ValueError, "refuses to open or hash"):
            BUILDER._file_bindings([forbidden])  # noqa: SLF001

    def test_tiny_eigensystem_uses_float64_values_and_float32_vectors(self) -> None:
        covariance = np.asarray([[4.0, 1.0], [1.0, 2.0]], dtype=np.float64)
        evals, evecs = BUILDER.fit_whitener_eigensystem_from_covariance(covariance)
        self.assertEqual(evals.dtype, np.dtype("<f8"))
        self.assertEqual(evecs.dtype, np.dtype("<f4"))
        self.assertTrue(np.all(evals > 0))
        np.testing.assert_allclose(
            evecs.astype(np.float64).T @ evecs.astype(np.float64),
            np.eye(2),
            rtol=1e-6,
            atol=1e-6,
        )

    def test_production_import_is_exact_a14_not_a15(self) -> None:
        production = BUILDER._load_production_module()  # noqa: SLF001
        self.assertEqual(production.ALPHA, 0.982)
        self.assertEqual(production.CD4_FRACTION, 0.05)
        self.assertIsNone(production.CD4_TARGET_F)
        self.assertEqual(production.PM_BLEND, 0.5)
        self.assertEqual(production.WHITEN_GAMMA, 0.65)
        self.assertEqual(production.WHITEN_LAM, 0.0)


if __name__ == "__main__":
    unittest.main()
