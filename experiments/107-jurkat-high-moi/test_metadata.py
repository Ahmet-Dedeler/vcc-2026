from __future__ import annotations

import gzip
import sys
import tempfile
import unittest
from pathlib import Path

import numpy as np


HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

import metadata  # noqa: E402


class MetadataUnitTests(unittest.TestCase):
    def test_fasta_requires_unique_ids_and_nonempty_sequences(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "guides.fa.gz"
            with gzip.open(path, "wt", encoding="utf-8") as handle:
                handle.write(">A_1\nACGT\n>A_1\nTGCA\n")
            with self.assertRaisesRegex(metadata.MetadataError, "duplicate"):
                metadata.read_fasta_ids(path)

            with gzip.open(path, "wt", encoding="utf-8") as handle:
                handle.write(">A_1\n>B_1\nACGT\n")
            with self.assertRaisesRegex(metadata.MetadataError, "no sequence"):
                metadata.read_fasta_ids(path)

    def test_guide_crosswalk_marks_every_unsupported_class(self) -> None:
        targets = (
            metadata.TargetRecord("A", 1, 4),
            metadata.TargetRecord("B", 2, 5),
            metadata.TargetRecord("B", 3, 6),
        )
        guide_ids = ("A_1", "B_1", "A~B_1", "missing_1", "non-targeting_1")
        crosswalk = metadata.build_guide_crosswalk(
            guide_ids,
            targets,
            ("A", "B", "missing"),
            archive_guide_ids=guide_ids,
        )

        self.assertEqual(
            crosswalk.guide_status,
            (
                metadata.GUIDE_SUPPORTED,
                metadata.GUIDE_AMBIGUOUS,
                metadata.GUIDE_COMPOSITE,
                metadata.GUIDE_UNMAPPED,
                metadata.GUIDE_NON_TARGETING,
            ),
        )
        np.testing.assert_array_equal(
            crosswalk.guide_to_target,
            [0, -1, -1, -1, -1],
        )
        np.testing.assert_array_equal(
            crosswalk.non_targeting_guides,
            [False, False, False, False, True],
        )
        self.assertEqual(crosswalk.ambiguous_labels, ("B",))
        self.assertEqual(crosswalk.composite_labels, ("A~B",))
        self.assertEqual(crosswalk.unmapped_labels, ("missing",))
        self.assertEqual(crosswalk.panel_covered_targets, ("A",))
        self.assertEqual(crosswalk.panel_missing_targets, ("B", "missing"))
        self.assertFalse(crosswalk.guide_to_target.flags.writeable)
        self.assertFalse(crosswalk.non_targeting_guides.flags.writeable)

    def test_guide_crosswalk_hard_fails_axis_mismatch_and_bad_ids(self) -> None:
        targets = (metadata.TargetRecord("A", 1, 4),)
        with self.assertRaisesRegex(metadata.MetadataError, "feature axes differ"):
            metadata.build_guide_crosswalk(
                ("A_1",),
                targets,
                ("A",),
                archive_guide_ids=("A_2",),
            )
        with self.assertRaisesRegex(metadata.MetadataError, "positive ordinal"):
            metadata.build_guide_crosswalk(("A",), targets, ("A",))
        with self.assertRaisesRegex(metadata.MetadataError, "positive ordinal"):
            metadata.build_guide_crosswalk(("A_0",), targets, ("A",))

    def test_hash_crosswalk_uses_workbook_condition_inventory(self) -> None:
        oligos = (
            metadata.HashOligo("A" * 45, "untreated", 4),
            metadata.HashOligo("C" * 45, "activated", 5),
            metadata.HashOligo("G" * 45, "activated", 6),
        )
        crosswalk = metadata.build_hash_crosswalk(
            ("activated1_1", "activated2_1", "untreated_1"),
            oligos,
        )
        self.assertEqual(
            crosswalk.label_to_condition,
            {
                "activated1_1": "activated",
                "activated2_1": "activated",
                "untreated_1": "untreated",
            },
        )

        with self.assertRaisesRegex(metadata.MetadataError, "condition counts"):
            metadata.build_hash_crosswalk(
                ("activated1_1", "untreated_1"),
                oligos,
            )
        with self.assertRaisesRegex(metadata.MetadataError, "contiguous"):
            metadata.build_hash_crosswalk(
                ("activated2_1", "untreated_1"),
                oligos[:2],
            )
        with self.assertRaisesRegex(metadata.MetadataError, "unsupported hash"):
            metadata.build_hash_crosswalk(("mystery_1",), oligos[:1])

    def test_response_crosswalk_reports_missing_and_ambiguous_axes(self) -> None:
        responses = (
            metadata.ResponseRecord("A", 1, "group", 4),
            metadata.ResponseRecord("B", 2, "group", 5),
            metadata.ResponseRecord("C", 3, "group", 6),
        )
        crosswalk = metadata.build_response_crosswalk(
            responses,
            ("A", "B", "B"),
            ("A", "C", "C"),
            ("A",),
        )
        self.assertEqual(crosswalk.missing_source_symbols, ("C",))
        self.assertEqual(crosswalk.ambiguous_source_symbols, ("B",))
        self.assertEqual(crosswalk.missing_vcc_symbols, ("B",))
        self.assertEqual(crosswalk.ambiguous_vcc_symbols, ("C",))
        self.assertEqual(crosswalk.visible_vcc_indices, ())
        self.assertEqual(crosswalk.records[0].vcc_gene_index, 0)
        self.assertIsNone(crosswalk.records[1].vcc_gene_index)
        self.assertIsNone(crosswalk.records[2].vcc_gene_index)

    def test_csv_axis_rejects_duplicates(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "axis.csv"
            path.write_text("gene_name\nA\nA\n", encoding="utf-8")
            with self.assertRaisesRegex(metadata.MetadataError, "duplicate"):
                metadata.read_single_column_csv(path, "gene_name")


class PinnedMetadataIntegrationTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.crosswalks = metadata.load_pinned_crosswalks()

    def test_canonical_guide_axis_and_classification(self) -> None:
        guide = self.crosswalks.guide
        self.assertEqual(len(guide.guide_ids), 83_401)
        self.assertEqual(
            guide.guide_ids[:5], ("A1BG_1", "A1BG_2", "A1BG_3", "A1BG_4", "A1CF_1")
        )
        self.assertEqual(guide.guide_ids[-1], "non-targeting_250")
        self.assertEqual(
            guide.status_counts,
            {
                "ambiguous": 5,
                "composite": 201,
                "non_targeting": 250,
                "supported": 82_945,
            },
        )
        self.assertEqual(guide.ambiguous_labels, ("MEMO1",))
        self.assertEqual(len(guide.composite_labels), 135)
        self.assertEqual(guide.unmapped_labels, ())
        self.assertEqual(
            guide.duplicate_target_symbols,
            (
                metadata.DuplicateSymbol(
                    symbol="MEMO1",
                    entrez_ids=(51072, 7795),
                    workbook_rows=(8146, 12175),
                ),
            ),
        )

    def test_panel_coverage_is_exactly_299_of_300(self) -> None:
        guide = self.crosswalks.guide
        self.assertEqual(len(guide.panel_targets), 300)
        self.assertEqual(len(guide.panel_covered_targets), 299)
        self.assertEqual(guide.panel_missing_targets, ("ZNF714",))

    def test_hash_labels_map_to_one_untreated_and_six_activated(self) -> None:
        hashes = self.crosswalks.hashes
        self.assertEqual(
            tuple(record.label for record in hashes.records),
            (
                "activated1_1",
                "activated2_1",
                "activated3_1",
                "activated4_1",
                "activated5_1",
                "activated6_1",
                "untreated_1",
            ),
        )
        self.assertEqual(
            tuple(record.condition for record in hashes.records),
            (
                "activated",
                "activated",
                "activated",
                "activated",
                "activated",
                "activated",
                "untreated",
            ),
        )
        self.assertEqual(
            [record.condition for record in hashes.workbook_oligos].count("untreated"),
            1,
        )
        self.assertEqual(
            [record.condition for record in hashes.workbook_oligos].count("activated"),
            6,
        )

    def test_response_axis_and_reports_are_exact(self) -> None:
        response = self.crosswalks.response
        self.assertEqual(len(response.records), 374)
        self.assertEqual(len({record.symbol for record in response.records}), 374)
        self.assertEqual(
            len({record.entrez_id for record in response.records}),
            374,
        )
        self.assertEqual(response.duplicate_response_symbols, ())
        self.assertEqual(response.duplicate_response_entrez_ids, ())
        self.assertEqual(response.missing_source_symbols, ("ATP6",))
        self.assertEqual(
            response.ambiguous_source_symbols,
            ("HSPA14", "SOD2"),
        )
        self.assertEqual(
            response.missing_vcc_symbols,
            (
                "MAG",
                "HLA-DRB1",
                "MTRNR2L10",
                "IFTAP",
                "H1-4",
                "H1-5",
                "H2BC5",
                "IFT74",
                "H2AC14",
                "H2AC4",
                "H2AC17",
                "H2AC20",
                "H2BC4",
                "FXYD6-FXYD2",
                "MRPL36",
                "ATP6",
            ),
        )
        self.assertEqual(response.ambiguous_vcc_symbols, ())
        self.assertEqual(
            sum(record.vcc_gene_index is not None for record in response.records),
            358,
        )
        self.assertEqual(len(response.visible_vcc_indices), 349)
        self.assertEqual(len(set(response.visible_vcc_indices)), 349)

    def test_reload_is_deterministic(self) -> None:
        repeated = metadata.load_pinned_crosswalks()
        self.assertEqual(
            repeated.guide.guide_ids,
            self.crosswalks.guide.guide_ids,
        )
        np.testing.assert_array_equal(
            repeated.guide.guide_to_target,
            self.crosswalks.guide.guide_to_target,
        )
        self.assertEqual(repeated.hashes.records, self.crosswalks.hashes.records)
        self.assertEqual(repeated.response.records, self.crosswalks.response.records)


if __name__ == "__main__":
    unittest.main()
