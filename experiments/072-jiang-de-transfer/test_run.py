from __future__ import annotations

import importlib.util
import sys
import tempfile
import unittest
import zipfile
from pathlib import Path

import numpy as np


HERE = Path(__file__).resolve().parent
SPEC = importlib.util.spec_from_file_location("exp072_run", HERE / "run.py")
if SPEC is None or SPEC.loader is None:
    raise ImportError("cannot load experiment 072 runner")
RUN = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = RUN
SPEC.loader.exec_module(RUN)


def table(rows: list[tuple[str, list[float | None]]]) -> str:
    header = " ".join(RUN.expected_header())
    output = [header]
    for gene, values in rows:
        encoded = ["NA" if value is None else str(value) for value in values]
        beta = ["0"] * len(RUN.CELL_LINES)
        pvalue = ["1"] * len(RUN.CELL_LINES)
        output.append(" ".join([gene, *encoded, *beta, *pvalue]))
    return "\n".join(output) + "\n"


def add_member(
    archive: zipfile.ZipFile,
    *,
    target: str,
    pathway: str,
    rows: list[tuple[str, list[float | None]]],
) -> None:
    archive.writestr(
        "DE_results_all_pathway/"
        f"Parse_{pathway}/{target}_{pathway}_pathway_DE_results.txt",
        table(rows),
    )


class JiangTransferTests(unittest.TestCase):
    def test_parse_member(self) -> None:
        info = zipfile.ZipInfo(
            "DE_results_all_pathway/Parse_IFNG/RUNX1_IFNG_pathway_DE_results.txt"
        )
        parsed = RUN.parse_member(info)
        self.assertIsNotNone(parsed)
        assert parsed is not None
        self.assertEqual(parsed.target, "RUNX1")
        self.assertEqual(parsed.pathway, "IFNG")
        self.assertIsNone(RUN.parse_member(zipfile.ZipInfo("__MACOSX/._x")))

    def test_pathway_centering_is_leave_one_target_out(self) -> None:
        genes = ("G1", "G2")
        panel = ("PANEL",)
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "fixture.zip"
            with zipfile.ZipFile(path, "w") as archive:
                add_member(
                    archive,
                    target="N1",
                    pathway="IFNG",
                    rows=[
                        ("G1", [1.0] * 6),
                        ("G2", [2.0] * 6),
                    ],
                )
                add_member(
                    archive,
                    target="N2",
                    pathway="IFNG",
                    rows=[
                        ("G1", [3.0] * 6),
                        ("G2", [4.0] * 6),
                    ],
                )
                add_member(
                    archive,
                    target="PANEL",
                    pathway="IFNG",
                    rows=[
                        ("G1", [5.0] * 6),
                        ("G2", [6.0] * 6),
                    ],
                )
            with zipfile.ZipFile(path) as archive:
                members = tuple(
                    RUN.parse_member(info)
                    for info in archive.infolist()
                    if RUN.parse_member(info) is not None
                )
                pathway_sum = np.zeros((5, 6, 2), dtype=np.float64)
                pathway_count = np.zeros((5, 6, 2), dtype=np.uint16)
                for member in members:
                    assert member is not None
                    if member.target == "PANEL":
                        continue
                    for gene, values in RUN.iter_effect_rows(archive, member):
                        column = genes.index(gene)
                        pathway_sum[1, :, column] += values
                        pathway_count[1, :, column] += 1
            audit = RUN.AuditState(
                members=tuple(member for member in members if member is not None),
                targets=("N1", "N2", "PANEL"),
                target_pathways={
                    "N1": ("IFNG",),
                    "N2": ("IFNG",),
                    "PANEL": ("IFNG",),
                },
                file_rows=(),
                panel_targets=("PANEL",),
                proxy_targets=("N1", "N2"),
                shared_panel_targets=("PANEL",),
                gene_seen=np.ones(2, dtype=bool),
                pathway_sum=pathway_sum,
                pathway_count=pathway_count,
                archive_sha256="fixture",
                archive_md5="fixture",
            )
            centered = RUN.build_centered_effects(
                archive_path=path,
                audit=audit,
                official_genes=genes,
                panel=panel,
                max_rss_gib=1.0,
            )
            np.testing.assert_allclose(
                centered.effects[centered.target_pos["N1"], :, :],
                -2.0,
            )
            np.testing.assert_allclose(
                centered.effects[centered.target_pos["N2"], :, :],
                2.0,
            )
            np.testing.assert_allclose(
                centered.effects[centered.target_pos["PANEL"], :, :],
                3.0,
            )

    def test_feature_axis_excludes_panel_and_proxy_targets(self) -> None:
        class TinyBank:
            def __init__(self) -> None:
                self.gene_pos = {"P": 0, "Q": 1, "G": 2}

        axis = RUN.build_feature_axis(
            official_genes=("P", "Q", "G"),
            panel={"P"},
            proxy_targets={"Q"},
            banks={"a": TinyBank(), "b": TinyBank(), "c": TinyBank()},
            jiang_gene_seen=np.asarray([True, True, True]),
        )
        self.assertEqual(axis, ("G",))

    def test_k562_is_omitted_only_when_k562_is_held(self) -> None:
        held_k562 = [
            RUN.CELL_LINES[index] for index in RUN.jiang_line_indices("k562_gwps")
        ]
        held_hct = [
            RUN.CELL_LINES[index] for index in RUN.jiang_line_indices("orion_hct116")
        ]
        self.assertNotIn("K562", held_k562)
        self.assertIn("K562", held_hct)

    def test_blend_preserves_missing_coordinates(self) -> None:
        baseline = np.asarray([[1.0, 2.0], [3.0, 4.0]])
        jiang = np.asarray([[2.0, 9.0], [8.0, 6.0]])
        observed = np.asarray([[True, False], [False, True]])
        candidate, scale = RUN.blend_candidate(
            baseline=baseline,
            jiang=jiang,
            observed=observed,
            fraction=0.25,
            frozen_scale=1.0,
        )
        self.assertEqual(scale, 1.0)
        self.assertEqual(candidate[0, 1], baseline[0, 1])
        self.assertEqual(candidate[1, 0], baseline[1, 0])
        self.assertEqual(candidate[0, 0], 1.25)
        self.assertEqual(candidate[1, 1], 4.5)

    def test_low_panel_coverage_fails_gate(self) -> None:
        checks = {
            "panel_targets_at_least_minimum": 9 >= 30,
            "shared_panel_targets_at_least_minimum": 9 >= 30,
            "mapped_vcc_genes_at_least_minimum": 10_000 >= 5_000,
            "proxy_targets_at_least_minimum": 200 >= 100,
        }
        self.assertFalse(all(checks.values()))


if __name__ == "__main__":
    unittest.main()
