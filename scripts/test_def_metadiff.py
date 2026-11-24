#!/usr/bin/env python3
"""Tests for scripts/def_metadiff.py, on 524dr2's finished run (658 work/dr2/run) and synthetic copies in new work dirs.
  M1 the real R-ABC final meta vs 524dr2: only de_dnmag + de_uprn differ (both shas are provenance-only)
  M2 a rebuild in a new work dir with identical arrays: no real diff, both shas provenance-only
  M3 one dump array changed -> bulk_load_sha256 is REAL;  M4 provenance text changed beyond the path -> sameline REAL
  M5 the same-line record's w changed -> sameline REAL
  M6 a dump hook: the run holds the stack dump AND dumphook/<same name>; the npz is chosen by the meta's sha
  M7 --allow / --need / --expect semantics;  M8 bad input -> rc 2;  M9 1-ULP float round-off;  M10 the near-zero noise floor (atol)
usage: uv run python scripts/test_def_metadiff.py
"""
from __future__ import annotations

import copy
import hashlib
import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
TOOL = ROOT / "scripts/def_metadiff.py"
DR2 = ROOT / "experiments/658-prereg-1007-builds/work/dr2/run"
META2 = ROOT / "submissions/524dr2-october-recipe-abc.h5ad.meta.json"
RABC = ROOT / "submissions/694rabc-october-recipe-abc.h5ad.meta.json"
DUMP, SL = "dump_h1x2_cd4x1.npz", "sameline.npz"


def sha(p: Path) -> str:
    return hashlib.sha256(p.read_bytes()).hexdigest()


def run(*args) -> subprocess.CompletedProcess:
    return subprocess.run([sys.executable, str(TOOL), *map(str, args)], capture_output=True, text=True)


@unittest.skipUnless(DR2.is_dir() and META2.is_file(), "524dr2's run dir is not on disk")
class MetaDiff(unittest.TestCase):
    def setUp(self):
        self._td = tempfile.TemporaryDirectory(prefix="t_metadiff_")
        self.td = Path(self._td.name)

    def tearDown(self):
        self._td.cleanup()

    def rebuild(self, tag, bump_array=False, bump_text=False, w=None, hook=False):
        """a copy of 524dr2's run in a NEW dir: the two npz re-saved with this dir inside their string fields"""
        r = self.td / tag / "run"
        r.mkdir(parents=True)
        new = str(r.resolve())
        shas = {}
        for name, key in ((DUMP, "env"), (SL, "provenance")):
            z = dict(np.load(DR2 / name, allow_pickle=False))
            s = str(z[key]).replace(str(DR2.resolve().relative_to(ROOT)), new)
            if bump_text and name == SL:
                s = s.replace("same-line", "programme")
            z[key] = np.array(s)
            if bump_array and name == DUMP:
                k = next(k for k in sorted(z) if z[k].dtype.kind == "f" and z[k].size)
                z[k] = z[k].copy(); z[k].flat[0] += 1.0
            np.savez(r / name, **z)
            shas[name] = sha(r / name)
        if hook:   # the hook's replacement: same name under dumphook/, one array changed; the meta records ITS sha
            z = dict(np.load(r / DUMP, allow_pickle=False))
            k = next(k for k in sorted(z) if z[k].dtype.kind == "f" and z[k].size)
            z[k] = z[k].copy(); z[k].flat[0] *= 1.5
            (r / "dumphook").mkdir()
            np.savez(r / "dumphook" / DUMP, **z)
            shas["hook"] = sha(r / "dumphook" / DUMP)
        m = json.loads(META2.read_text())
        m["config"]["bulk_load_sha256"] = shas["hook"] if hook else shas[DUMP]
        m["config"]["sameline"]["sha256"] = shas[SL]
        m["config"]["sameline"]["path"] = new + "/" + SL
        if w is not None:
            m["config"]["sameline"]["w"] = w
        mp = self.td / tag / "final.meta.json"
        mp.write_text(json.dumps(m))
        return mp, r

    @unittest.skipUnless(RABC.is_file() and (ROOT / "experiments/694-rabc/work/run" / DUMP).is_file(), "R-ABC's run not on disk")
    def test_m1_real_rabc(self):
        p = run(RABC, ROOT / "experiments/694-rabc/work/run", META2, DR2, "--expect", "de_dnmag,de_uprn")
        self.assertEqual(p.returncode, 0, p.stdout + p.stderr)
        self.assertIn("provenance-only ['bulk_load_sha256', 'sameline.sha256']", p.stdout)

    def test_m2_provenance_only(self):
        mp, r = self.rebuild("m2")
        p = run(mp, r, META2, DR2)
        self.assertEqual(p.returncode, 0, p.stdout + p.stderr)
        self.assertIn("diff [] ", p.stdout)
        self.assertIn("provenance-only ['bulk_load_sha256', 'sameline.sha256']", p.stdout)

    def test_m3_array_changed(self):
        mp, r = self.rebuild("m3", bump_array=True)
        p = run(mp, r, META2, DR2)
        self.assertEqual(p.returncode, 1)
        self.assertIn("bulk_load_sha256: REAL: array", p.stdout)
        self.assertIn("diff ['bulk_load_sha256']", p.stdout)

    def test_m4_text_beyond_path(self):
        mp, r = self.rebuild("m4", bump_text=True)
        p = run(mp, r, META2, DR2)
        self.assertEqual(p.returncode, 1)
        self.assertIn("sameline: REAL: provenance differs beyond the work dir", p.stdout)

    def test_m5_record_changed(self):
        mp, r = self.rebuild("m5", w=0.6)
        p = run(mp, r, META2, DR2)
        self.assertEqual(p.returncode, 1)
        self.assertIn("the same-line record differs beyond its sha256", p.stdout)

    def test_m6_hook_dump_chosen_by_sha(self):
        mp, r = self.rebuild("m6", hook=True)
        p = run(mp, r, META2, DR2)
        self.assertEqual(p.returncode, 1, p.stdout)
        self.assertIn("bulk_load_sha256: REAL: array", p.stdout)   # the hook's dump was compared, not the stack dump

    def test_m7_rules(self):
        mp, r = self.rebuild("m7", w=0.6)
        self.assertEqual(run(mp, r, META2, DR2, "--allow", "sameline,bulk_ctx_scale").returncode, 0)
        self.assertEqual(run(mp, r, META2, DR2, "--allow", "sameline", "--need", "sameline").returncode, 0)
        self.assertEqual(run(mp, r, META2, DR2, "--allow", "bulk_ctx_scale").returncode, 1)
        self.assertEqual(run(mp, r, META2, DR2, "--allow", "sameline", "--need", "pert_scale").returncode, 1)
        self.assertEqual(run(mp, r, META2, DR2, "--expect", "sameline").returncode, 0)
        self.assertEqual(run(mp, r, META2, DR2, "--expect", "").returncode, 1)

    def test_m9_float_roundoff(self):
        """a 1-ULP float change inside a record is not a difference (default rtol 1e-12); a real change still is; rtol 0 is exact"""
        mp, r = self.rebuild("m9")
        m = json.loads(mp.read_text()); pt = m["config"]["bulk_sharpen_pt"]
        pt["edges"] = [float(np.nextafter(e, 1.0)) for e in pt["edges"]]
        mp.write_text(json.dumps(m))
        self.assertEqual(run(mp, r, META2, DR2).returncode, 0)
        self.assertEqual(run(mp, r, META2, DR2, "--float-rtol", "0", "--float-atol", "0").returncode, 1)
        pt["edges"] = [e * 1.001 for e in pt["edges"]]
        mp.write_text(json.dumps(m))
        p = run(mp, r, META2, DR2)
        self.assertEqual(p.returncode, 1)
        self.assertIn("bulk_sharpen_pt", p.stdout)

    def test_m10_noise_floor(self):
        """698 (10-18): a round-off diagnostic near zero (tercilec key_ratio_maxdev 8.66e-15 vs 1.02e-14, 15 % relative) is not a
        difference under the default atol 1e-12; with atol 0 it is; a change above the floor (1e-9) still is"""
        mp, r = self.rebuild("m10")
        m = json.loads(mp.read_text()); pt = m["config"]["bulk_sharpen_pt"]
        pt["key_ratio_maxdev"] = 1.021405182655144e-14
        mp.write_text(json.dumps(m))
        self.assertEqual(run(mp, r, META2, DR2).returncode, 0)
        self.assertEqual(run(mp, r, META2, DR2, "--float-atol", "0").returncode, 1)
        pt["key_ratio_maxdev"] = 1e-9
        mp.write_text(json.dumps(m))
        p = run(mp, r, META2, DR2)
        self.assertEqual(p.returncode, 1)
        self.assertIn("bulk_sharpen_pt", p.stdout)

    def test_m8_bad_input(self):
        self.assertEqual(run(self.td / "nope.json", DR2, META2, DR2).returncode, 2)
        self.assertEqual(run(META2, self.td / "nodir", META2, DR2).returncode, 2)


if __name__ == "__main__":
    unittest.main(verbosity=2)
