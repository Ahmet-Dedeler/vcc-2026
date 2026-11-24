"""Regression tests for the false-completion modes seen during takeover."""

import unittest

from round5_status import classify


class ClassificationTests(unittest.TestCase):
    def state(self, **overrides):
        values = dict(
            launched=True,
            alive=False,
            events=0,
            tools=0,
            terminal=None,
            result_exists=False,
            preaccess_only=False,
        )
        values.update(overrides)
        return classify(**values)

    def test_empty_trace_and_dead_pid_is_not_done(self):
        self.assertEqual(self.state(), "LAUNCH_FAILED_NO_EVENTS")

    def test_partial_trace_is_interrupted(self):
        self.assertEqual(self.state(events=12, tools=2), "INTERRUPTED_NO_TERMINAL")

    def test_early_results_file_does_not_imply_completion(self):
        self.assertEqual(
            self.state(alive=True, events=12, tools=2, result_exists=True), "RUNNING"
        )

    def test_terminal_error_is_not_success_even_with_report(self):
        self.assertEqual(
            self.state(
                terminal={"subtype": "error", "is_error": True},
                result_exists=True,
            ),
            "TERMINAL_ERROR",
        )

    def test_unknown_error_field_fails_closed(self):
        self.assertEqual(
            self.state(terminal={"subtype": "success"}, result_exists=True),
            "TERMINAL_ERROR",
        )

    def test_success_without_artifact_needs_attention(self):
        self.assertEqual(
            self.state(terminal={"subtype": "success", "is_error": False}),
            "TERMINAL_MISSING_ARTIFACT",
        )

    def test_preaccess_success_is_not_a_scientific_result(self):
        self.assertEqual(
            self.state(
                terminal={"subtype": "success", "is_error": False},
                result_exists=True,
                preaccess_only=True,
            ),
            "PREACCESS_REPORT_READY",
        )

    def test_final_result_still_requires_audit(self):
        self.assertEqual(
            self.state(
                terminal={"subtype": "success", "is_error": False},
                result_exists=True,
            ),
            "RESULT_READY_FOR_AUDIT",
        )


if __name__ == "__main__":
    unittest.main()
