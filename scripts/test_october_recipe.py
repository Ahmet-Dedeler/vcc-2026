"""scripts/october_recipe.py: the decisions file -> def_day.sh call mapping, on temp copies (the real files are untouched).

  uv run python scripts/test_october_recipe.py
"""
import json
import os
import shutil
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts/october_recipe.py"
DD497 = ROOT / "experiments/497-october-switches/work/def_day.497.sh"
# 538 (applied 2026-10-05) adds six PENDING switches. Tests that need `command` to resolve without --allow-pending
# decide them off (10-05 CHECKLIST v2 step 5b; 650 L5).
OFF538 = [("rho_scale", "off"), ("rc_callset_cpm", "off"), ("hook_estimator", "off"), ("lodo_strength", "off"),
          ("lodo_hi", "off"), ("callset_block", "off"),
          # 690a: three more PENDING switches, decided off for the same tests
          ("de_uprn", "off"), ("mde_sign_keygate", "off"), ("de_dnmag", "off"),
          ("dump_hook", "off")]   # 690g


class RecipeTests(unittest.TestCase):
    def setUp(self):
        self._td = tempfile.TemporaryDirectory()
        self.tmp = Path(self._td.name)
        shutil.copy(ROOT / "docs/october-decisions.json", self.tmp / "d.json")

    def tearDown(self):
        self._td.cleanup()

    def run_(self, *args, def_day=None, builder=None):
        env = dict(os.environ, OCTOBER_DECISIONS=str(self.tmp / "d.json"))
        if def_day:
            env["OCTOBER_DEF_DAY"] = str(def_day)
        if builder:
            env["OCTOBER_BUILDER"] = str(builder)
        r = subprocess.run([sys.executable, str(SCRIPT), *args], env=env, capture_output=True, text=True)
        return r.returncode, r.stdout, r.stderr

    def cmd(self, *extra, **kw):
        return self.run_("command", "--bundle", "B", "--work", "W", "--name", "N", *extra, **kw)

    def decide(self, pairs):
        for sw, v in pairs:
            rc, out, err = self.run_("decide", sw, v, "--evidence", "test")
            self.assertEqual(rc, 0, err)

    def test_live_state_refuses_pending_then_resolves_with_allow(self):
        # 690a: state-aware. Plain `command` refuses (rc 2) iff a switch is PENDING; in the all-decided state (Oct 21) it
        # resolves (rc 0), and the --allow-pending checks below hold either way.
        live = json.loads((Path(__file__).resolve().parents[1] / "docs/october-decisions.json").read_text())
        pending = [k for k, v in live.get("switches", live).items() if v.get("status") == "PENDING"]
        rc, out, err = self.cmd()
        if pending:
            self.assertEqual(rc, 2)
            self.assertIn("PENDING", err)
        else:
            self.assertEqual(rc, 0, err)
        rc, out, err = self.cmd("--allow-pending")
        # de_mag was decided 2026-09-28 (dn_k250_zs). Until 502's pass-through is in scripts/def_day.sh the command is
        # refused (rc 3) and printed commented out; once it is in, the command resolves with the rider (519).
        live_dd = (Path(__file__).resolve().parent / "def_day.sh").read_text()
        if 'export VCC_DE_MAG="$DEF_DAY_DE_MAG"' in live_dd:
            self.assertEqual(rc, 0, err)
            # The live command grows with every `decide`, so pin its invariants and each decided switch's export instead
            # of one string (2026-10-02: sharpen_pt's decide broke the exact-string form).
            last = out.strip().splitlines()[-1]
            self.assertTrue(last.endswith("bash scripts/def_day.sh B W N 1.75 0.5 noprog 4"), last)
            for part in ("DEF_DAY_LODO=1", "DEF_DAY_DE_MAG=dn_k250_zs"):
                self.assertIn(part, last)
            sw = json.loads((Path(__file__).resolve().parents[1] / "docs/october-decisions.json").read_text())
            sw = sw.get("switches", sw)
            spt = sw.get("sharpen_pt", {})
            if spt.get("status") == "DECIDED" and spt.get("value") not in (None, "off"):
                self.assertIn(f"DEF_DAY_SHARPEN_PT={spt['value']}", last)
        else:
            self.assertEqual(rc, 3, err)
            self.assertIn("has no DEF_DAY_DE_MAG", err)
            self.assertTrue(out.strip().splitlines()[-1].startswith("# REFUSED DEF_DAY_LODO=1 DEF_DAY_DE_MAG=dn_k250_zs"))

    def test_all_switches_on_with_the_497_script(self):
        self.decide([("lodo_hi_pct", "96.64697522769474"), ("sharpen_pt", "stercile:2.5,2.0,1.75"),
                     ("programme_source", "h1_2025_full"), ("prior_gate", "le1"), ("de_mag", "off"), ("mde_sign", "off")] + OFF538)
        rc, out, err = self.cmd(def_day=DD497)
        self.assertEqual(rc, 0, err)
        self.assertEqual(out.strip().splitlines()[-1],
                         "DEF_DAY_LODO=1 DEF_DAY_LODO_HI_PCT=96.64697522769474 DEF_DAY_SHARPEN_PT=stercile:2.5,2.0,1.75 "
                         "DEF_DAY_PROGRAMME_SOURCE=h1_2025_full DEF_DAY_PRIOR_GATE=le1 "
                         "bash scripts/def_day.sh B W N 1.75 0.5 '' 4")
        hist = json.loads((self.tmp / "d.json").read_text())["switches"]["prior_gate"]["history"]
        self.assertEqual(hist[-1]["to"], "le1")

    def test_riders_need_502(self):
        self.decide([("lodo_hi_pct", "off"), ("sharpen_pt", "off"), ("programme_source", "off"), ("prior_gate", "off"),
                     ("de_mag", "dn_k250_zs"), ("mde_sign", "clean4:0.05")] + OFF538)
        rc, out, err = self.cmd(def_day=DD497)
        self.assertEqual(rc, 3)
        self.assertIn("has no DEF_DAY_DE_MAG", err)
        self.assertIn("has no DEF_DAY_MDE_SIGN", err)

    def test_dependency_violations(self):
        self.decide([("lodo_hi_pct", "off"), ("sharpen_pt", "off"), ("programme_source", "off"), ("prior_gate", "le1"),
                     ("de_mag", "off"), ("mde_sign", "off")] + OFF538)
        rc, out, err = self.cmd(def_day=DD497)
        self.assertEqual(rc, 3)
        self.assertIn("prior_gate is on without programme_source", err)
        self.decide([("prior_gate", "off"), ("cd4_stack", "x1.5")])
        rc, out, err = self.cmd(def_day=DD497)
        self.assertEqual(rc, 3)
        self.assertIn("cd4_stack", err)

    def test_empty_flags_need_the_497_script(self):
        self.decide([("lodo_hi_pct", "off"), ("sharpen_pt", "off"), ("programme_source", "h1_2025_full"),
                     ("prior_gate", "le1"), ("de_mag", "off"), ("mde_sign", "off")] + OFF538)
        old = self.tmp / "old_def_day.sh"
        old.write_text(DD497.read_text().replace('SLF="${6-pt:noprog}"', 'SLF="${6:-pt:noprog}"'))
        rc, out, err = self.cmd(def_day=old)
        self.assertEqual(rc, 3)
        self.assertIn("empty flags argument", err)
        rc, out, err = self.cmd(def_day=DD497)
        self.assertEqual(rc, 0, err)

    def test_decide_requires_evidence_and_known_switch(self):
        self.assertEqual(self.run_("decide", "nope", "1", "--evidence", "x")[0], 1)
        self.assertEqual(self.run_("decide", "q", "2.0", "--evidence", " ")[0], 1)


if __name__ == "__main__":
    unittest.main()
