"""Offline checks of the takeover's arithmetic and retained provenance."""
import json
from pathlib import Path
import subprocess
import sys
import unittest

from audit_takeover_state import decompose

ROOT = Path(__file__).resolve().parents[1]
IDS = ("Q519I2pLWRCW7PAEx5ul", "1n5LFk6PSEHq15UHENIE",
       "2K5ppxkowMSj2ugoZqYb", "xuDelEik7zGWMXdNu0R4")


class TakeoverTests(unittest.TestCase):
    def test_published_gap_is_sum_of_six(self):
        state = json.loads((ROOT / "docs/review-2026-09-09/state.json").read_text())
        result = decompose(state["entry"], state["top_six"][0])
        self.assertAlmostEqual(result["overall_delta"], 0.0434327243975573)
        self.assertLess(result["overall_contributions"]["mse"], 0)
        self.assertEqual(state["entry"]["score_avg"], state["best_recorded_entry"]["score_avg"])

    def test_invalid_mean_fails_closed(self):
        row = {"score_avg": 0, **{f"score_{m}": 0 for m in
                                 ("pds", "mse", "nmae", "fid", "reach", "jac")}}
        with self.assertRaises(AssertionError):
            decompose(row, {**row, "score_avg": 1})

    def test_four_backfills_match_independently_extracted_hashes(self):
        ledger = json.loads((ROOT / "submissions/SUBMITTED_LEDGER.json").read_text())
        audit = json.loads((ROOT / "docs/storage-sweep-2026-09-08/audit-2026-09-09.json").read_text())
        by_id = {r.get("captured_entry_id"): r for r in audit["packages"]}
        for entry_id in IDS:
            row = by_id[entry_id]
            recorded = ledger[row["sha256"]]
            self.assertEqual(recorded["entry_id"], entry_id)
            self.assertEqual(recorded["content_sha256"], row["novelty_content_sha256"])
            self.assertEqual(recorded["h5ad_sha256"], row["source_h5ad_sha256"])
            self.assertTrue(row["payload_readable"])
            self.assertTrue(all(c["match"] for c in row["storage_manifest_checks"]))

    def test_novelty_cli_refuses_all_four_backfilled_packages(self):
        # Real CLI, real package hashes. Refusal happens before any decompression.
        for entry_id in IDS:
            cap = json.loads((ROOT / f"data/board-entries/{entry_id}.json").read_text())
            result = subprocess.run(
                [sys.executable, str(ROOT / "scripts/check_submission_novelty.py"),
                 "check", cap["package"]], cwd=ROOT, capture_output=True, text=True)
            self.assertEqual(result.returncode, 1, result.stdout + result.stderr)
            self.assertIn("ALREADY BEEN SCORED", result.stderr)
            self.assertIn(entry_id, result.stderr)


if __name__ == "__main__":
    unittest.main()
