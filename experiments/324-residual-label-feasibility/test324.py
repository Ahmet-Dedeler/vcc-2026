"""Known-answer metadata graphs and deterministic shared-noise falsifier."""

import unittest
from audit324 import design


def row(dataset, target, line, study, method="shRNA", stable=""):
    return {"Dataset": dataset, "pbgene": target, "CellLineName": line,
            "GSE": study, "method": method, "Cellosaurus": stable}


class DesignTests(unittest.TestCase):
    def test_two_studies_and_lines_are_not_a_crossed_design(self):
        rows = [row("1", "T", "A", "S1"), row("2", "T", "B", "S2")]
        summary, _, _ = design(rows, "raw", True)
        self.assertEqual(summary["same_context_accession_pairs"], 0)
        self.assertEqual(summary["crossed_rectangles"], 0)

    def test_one_accession_many_groups_is_not_replication(self):
        rows = [row(str(i), "T", "A", "S1") for i in range(3)]
        summary, _, _ = design(rows, "raw", True)
        self.assertEqual(summary["unique_target_method_context_accession_edges"], 1)
        self.assertEqual(summary["same_context_accession_pairs"], 0)

    def test_crossed_rectangle_and_exact_group_binding(self):
        rows = [row(f"{s}{l}", "T", l, s) for s in ("S1", "S2") for l in ("A", "B")]
        rows += [row("extra", "T", "A", "S1")]
        summary, pairs, rects = design(rows, "raw", True)
        self.assertEqual(len(pairs), 2)
        self.assertEqual(summary["crossed_rectangles"], 1)
        self.assertEqual(rects[0]["datasets_aa"], "S1A|extra")

    def test_method_mixing_is_not_a_same_method_rectangle(self):
        rows = [row(f"{s}{l}", "T", l, s, "shRNA" if s == "S1" else "knockout")
                for s in ("S1", "S2") for l in ("A", "B")]
        self.assertEqual(design(rows, "raw", False)[0]["crossed_rectangles"], 1)
        self.assertEqual(design(rows, "raw", True)[0]["crossed_rectangles"], 0)

    def test_stable_identity_prevents_alias_double_count(self):
        rows = [row(f"{s}{l}", "T", l, s, stable="CVCL_0001")
                for s in ("S1", "S2") for l in ("AliasA", "AliasB")]
        self.assertEqual(design(rows, "raw", True)[0]["crossed_rectangles"], 1)
        self.assertEqual(design(rows, "id_only", True)[0]["crossed_rectangles"], 0)

    def test_missing_identity_not_imputed_in_strict_subset(self):
        rows = [row(f"{s}{l}", "T", l, s)
                for s in ("S1", "S2") for l in ("A", "B")]
        self.assertEqual(design(rows, "id_only", True)[0]["crossed_rectangles"], 0)
        self.assertEqual(design(rows, "id_fallback", True)[0]["crossed_rectangles"], 1)

    def test_shared_estimate_noise_can_create_residual_repeatability(self):
        # Orthogonal noise views, zero true residual. No biological data used.
        ya, yb = [1., 0., 0., 0.], [0., 1., 0., 0.]
        ua, ub = [0., 0., 2., 0.], [0., 0., 0., 2.]
        dot = lambda a, b: sum(x * y for x, y in zip(a, b))
        subtract = lambda a, b: [x - y for x, y in zip(a, b)]
        self.assertEqual(dot(ya, yb), 0.)
        self.assertEqual(dot(subtract(ya, ua), subtract(yb, ua)), 4.)
        self.assertEqual(dot(subtract(ya, ua), subtract(yb, ub)), 0.)


if __name__ == "__main__":
    unittest.main()
