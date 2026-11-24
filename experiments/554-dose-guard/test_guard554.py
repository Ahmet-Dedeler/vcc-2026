#!/usr/bin/env python3
"""554 -- tests for the dose-guard patch. Temp trees only: the live files the patch targets (and 538's nine) are hashed at the
start and at the end and must not change.

  A  patch mechanics on a 538-patched tree: --check clean and writes nothing; apply; a second apply refused; a partial
     tree refused; an anchor miss leaves every file untouched.
  B  dose_guard (imported from the patched def_postbuild.py): synthetic metas/dumps (healthy, one context collapsed, the
     ALLOW override, a dump that is not the meta's, the pooled fallback, a non-K-class meta) and the real dumps
     (433a healthy, 429a the known failure, a cross-paired mismatch).
  C  def_day.sh on 538's stub harness (stub uv, no data): with every knob unset the 554 day makes the same calls as the
     538 day except VCC_RC_AUDIT_DUMP on the level builds and --rc-audit on def_postbuild; def_day.final.env identical;
     the final dump is linked and the non-final ones pruned; with the two overrides set both are logged and recorded.

    uv run python experiments/554-dose-guard/test_guard554.py -v
The slow real-data identity checks (def_preflight / def_postbuild text and exit codes, base vs patched) are identity554.py.
"""
from __future__ import annotations

import hashlib
import importlib.util
import json
import os
import shutil
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
PY = ROOT / ".venv/bin/python"
PATCH = HERE / "patch_guard554.py"
APPLIER538 = ROOT / "experiments/538-cloud-bundle-extras/apply_october538.py"
sys.path.insert(0, str(ROOT / "experiments/538-cloud-bundle-extras"))
import test_extras538 as T538  # noqa: E402  (538's stub harness; its own tests are not re-run here)

REL538 = list(T538.REL)
EXTRA = [Path("scripts/def_preflight.py"), Path("scripts/frozen_inputs.json"), Path("scripts/make_submission_pm_cd4.py")]
PATCHED = [Path("scripts/def_preflight.py"), Path("scripts/def_postbuild.py"), Path("scripts/def_day.sh")]
LIVE = sorted({ROOT / p for p in REL538 + PATCHED})


def sha(p: Path) -> str:
    return hashlib.sha256(p.read_bytes()).hexdigest()


def live_hashes() -> dict:
    return {str(p): sha(p) for p in LIVE}


def make_tree(dest: Path, with554: bool) -> Path:
    """A tree holding copies of the files 538 and 554 touch (+ preflight's frozen table and the builder postbuild parses),
    538 applied, 554 applied if asked; data/ and submissions/ are read-only symlinks to the live repo."""
    for rel in REL538 + EXTRA:
        t = dest / rel
        t.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy(ROOT / rel, t)
    r = subprocess.run([str(PY), str(APPLIER538), "--root", str(dest)], capture_output=True, text=True)
    assert r.returncode == 0, r.stdout + r.stderr
    if with554:
        r = subprocess.run([str(PY), str(PATCH), "--root", str(dest)], capture_output=True, text=True)
        assert r.returncode == 0, r.stdout + r.stderr
    for d in ("data", "submissions"):
        os.symlink(ROOT / d, dest / d)
    return dest


def load_postbuild(tree: Path):
    spec = importlib.util.spec_from_file_location(f"postbuild554_{abs(hash(str(tree)))}", tree / "scripts/def_postbuild.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def meta_from_rows(rows: list[dict]) -> dict:
    rc = {"blocks": len(rows)}
    for k in ("prior_lambda", "mde_lambda", "call_set_missed"):
        v = [float(r[k]) for r in rows if k in r]
        if v:
            rc[k] = {"median": float(np.median(v)), "min": float(min(v)), "max": float(max(v))}
    return {"output": {"rc_callset": rc}}


def rows_for(ctx_vals: dict[str, tuple[float, float, float]], n: int = 4) -> list[dict]:
    return [{"context": c, "target": f"T{i}", "prior_lambda": p, "mde_lambda": m, "call_set_missed": s}
            for c, (p, m, s) in ctx_vals.items() for i in range(n)]


class Guard554(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls._live = live_hashes()
        cls.tmp = Path(tempfile.mkdtemp(prefix="guard554-"))
        cls.base = make_tree(cls.tmp / "base", with554=False)
        cls.tree = make_tree(cls.tmp / "patched", with554=True)
        cls.pb = load_postbuild(cls.tree)

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.tmp, ignore_errors=True)
        if live_hashes() != cls._live:
            raise AssertionError("a test wrote a live file")

    def _env(self, **kw):
        old = {k: os.environ.get(k) for k in kw}
        for k, v in kw.items():
            if v is None:
                os.environ.pop(k, None)
            else:
                os.environ[k] = v
        return old

    def _restore(self, old):
        for k, v in old.items():
            if v is None:
                os.environ.pop(k, None)
            else:
                os.environ[k] = v

    # ------------------------------------------------------------------ A: patch mechanics
    def test_A_check_apply_refuse_anchor(self):
        d = self.tmp / "mech"
        for rel in REL538 + EXTRA:
            (d / rel).parent.mkdir(parents=True, exist_ok=True)
            shutil.copy(ROOT / rel, d / rel)
        subprocess.run([str(PY), str(APPLIER538), "--root", str(d)], check=True, capture_output=True)
        before = {str(p): sha(d / p) for p in PATCHED}
        r = subprocess.run([str(PY), str(PATCH), "--root", str(d), "--check"], capture_output=True, text=True)
        self.assertEqual(r.returncode, 0, r.stdout + r.stderr)
        self.assertIn("patches apply cleanly (not written)", r.stdout)
        self.assertEqual(before, {str(p): sha(d / p) for p in PATCHED})
        r = subprocess.run([str(PY), str(PATCH), "--root", str(d)], capture_output=True, text=True)
        self.assertEqual(r.returncode, 0, r.stdout + r.stderr)
        after = {str(p): sha(d / p) for p in PATCHED}
        self.assertTrue(all(after[k] != before[k] for k in before))
        again = subprocess.run([str(PY), str(PATCH), "--root", str(d)], capture_output=True, text=True)
        self.assertNotEqual(again.returncode, 0)
        self.assertIn("already patched", again.stdout + again.stderr)
        self.assertEqual(after, {str(p): sha(d / p) for p in PATCHED})
        for p in PATCHED:                                           # a tree with one patched file is refused
            shutil.copy(self.base / p, d / p)
        shutil.copy(self.tree / PATCHED[0], d / PATCHED[0])
        part = subprocess.run([str(PY), str(PATCH), "--root", str(d), "--check"], capture_output=True, text=True)
        self.assertNotEqual(part.returncode, 0)
        self.assertIn("partially patched", part.stdout + part.stderr)
        for p in PATCHED:                                           # an anchor miss writes nothing
            shutil.copy(self.base / p, d / p)
        pb = d / "scripts/def_postbuild.py"
        pb.write_text(pb.read_text().replace("def main() -> None:\n", "def main():\n", 1))
        b2 = {str(p): sha(d / p) for p in PATCHED}
        miss = subprocess.run([str(PY), str(PATCH), "--root", str(d)], capture_output=True, text=True)
        self.assertNotEqual(miss.returncode, 0)
        self.assertIn("anchor found 0 times", miss.stdout + miss.stderr)
        self.assertEqual(b2, {str(p): sha(d / p) for p in PATCHED})

    def test_A_bash_and_python_syntax(self):
        self.assertEqual(subprocess.run(["bash", "-n", str(self.tree / "scripts/def_day.sh")]).returncode, 0)
        for p in ("scripts/def_preflight.py", "scripts/def_postbuild.py"):
            r = subprocess.run([str(PY), "-m", "py_compile", str(self.tree / p)], capture_output=True, text=True)
            self.assertEqual(r.returncode, 0, r.stderr)

    # ------------------------------------------------------------------ B: dose_guard
    def _dump(self, rows) -> Path:
        f = Path(tempfile.mkstemp(dir=self.tmp, suffix=".jsonl")[1])
        f.write_text("".join(json.dumps(r) + "\n" for r in rows))
        return f

    def test_B_synthetic(self):
        old = self._env(DEF_DAY_ALLOW_LOW_DOSE=None)
        try:
            ok = rows_for({"D": (1.0, 0.956, 5.0), "E": (1.0, 1.0, 17.0), "F": (1.0, 1.0, 11.0)})
            g = self.pb.dose_guard(meta_from_rows(ok), self._dump(ok), ["D", "E", "F"])
            self.assertEqual(g["rows"], [])
            self.assertAlmostEqual(g["per_context"]["D"]["margin"]["mde_lambda"], 0.056, places=9)
            bad = rows_for({"D": (1.0, 0.0, 5.0), "E": (1.0, 1.0, 17.0), "F": (1.0, 1.0, 11.0)})
            g = self.pb.dose_guard(meta_from_rows(bad), self._dump(bad), ["D", "E", "F"])
            self.assertEqual(len(g["rows"]), 1)
            name, okc, value = g["rows"][0]
            self.assertFalse(okc)
            self.assertIn("D: median mde_lambda 0 past its floor 0.9", value)
            self.assertNotIn("E:", value)
            low_prior = rows_for({"D": (0.43, 1.0, 5.0), "E": (1.0, 1.0, 17.0), "F": (1.0, 1.0, 150.0)})
            g = self.pb.dose_guard(meta_from_rows(low_prior), self._dump(low_prior), ["D", "E", "F"])
            self.assertIn("D: median prior_lambda 0.43", g["rows"][0][2])
            self.assertIn("F: median call_set_missed 150", g["rows"][0][2])
            os.environ["DEF_DAY_ALLOW_LOW_DOSE"] = "1"
            g = self.pb.dose_guard(meta_from_rows(bad), self._dump(bad), ["D", "E", "F"])
            self.assertTrue(g["rows"][0][1])
            self.assertIn("ALLOWED by DEF_DAY_ALLOW_LOW_DOSE=1", g["rows"][0][2])
            self.assertTrue(g["allowed"])
            os.environ.pop("DEF_DAY_ALLOW_LOW_DOSE")
            two = rows_for({"D": (1.0, 0.0, 5.0), "E": (1.0, 0.0, 17.0), "F": (1.0, 1.0, 11.0)})
            g = self.pb.dose_guard(meta_from_rows(ok), self._dump(two), ["D", "E", "F"])     # not this meta's dump
            self.assertIn("does not belong to this meta", g["rows"][0][2])
            # (a dump from another build with the same block count AND the same pooled medians is not detectable here:
            #  `bad` vs `ok` share both; the check guards against stale or mixed files, not a forged twin)
            g = self.pb.dose_guard(meta_from_rows(ok), self._dump(ok[:-1]), ["D", "E", "F"])  # a block short
            self.assertIn("does not belong to this meta", g["rows"][0][2])
            g = self.pb.dose_guard(meta_from_rows(ok), None, ["D", "E", "F"])                # pooled fallback, healthy
            self.assertEqual(g["rows"], [])
            self.assertIn("pooled", g["source"])
            collapse = rows_for({"D": (0.2, 1.0, 5.0), "E": (0.3, 1.0, 5.0), "F": (1.0, 1.0, 5.0)})
            g = self.pb.dose_guard(meta_from_rows(collapse), None, ["D", "E", "F"])          # pooled fallback, collapse
            self.assertIn("pooled: median prior_lambda 0.3", g["rows"][0][2])
            g = self.pb.dose_guard({"output": {}}, None, ["D", "E", "F"])                   # not a K-class build
            self.assertEqual(g["rows"], [])
            self.assertIn("none", g["source"])
        finally:
            self._restore(old)

    def test_B_real_dumps(self):
        old = self._env(DEF_DAY_ALLOW_LOW_DOSE=None)
        try:
            m433 = json.loads((ROOT / "submissions/433a-cap15-partial50-clean4q80-pertscale3-ballast-sharpen125-stackload-mde010coarse.h5ad.meta.json").read_text())
            m429 = json.loads((ROOT / "submissions/429a-cap15-partial50-clean4q80-pertscale3-ballast-sharpen125-stackload-mde010-zr5-ctxlevel-zs250.h5ad.meta.json").read_text())
            d433 = ROOT / "experiments/433-coarse-mde-arm/work/rc_audit_433a.jsonl"
            d429 = ROOT / "experiments/429-zs-production/work/rc_audit_429a.jsonl"
            g = self.pb.dose_guard(m433, d433, ["A", "B", "C"])
            self.assertEqual(g["rows"], [])
            self.assertAlmostEqual(g["per_context"]["A"]["mde_lambda"], 0.9560546875, places=9)
            g = self.pb.dose_guard(m429, d429, ["A", "B", "C"])
            v = g["rows"][0][2]
            self.assertFalse(g["rows"][0][1])
            self.assertIn("A: median mde_lambda 0 past", v)
            self.assertIn("C: median mde_lambda 0.855", v)
            self.assertNotIn("B:", v)
            g = self.pb.dose_guard(m429, d433, ["A", "B", "C"])
            self.assertIn("does not belong to this meta", g["rows"][0][2])
        finally:
            self._restore(old)

    # ------------------------------------------------------------------ C: def_day.sh on 538's stub harness
    def _day(self, kind: str, env: dict) -> dict:
        h = T538.Extras538("_assert_day_ok")
        h.star = f"{T538.rho_star():.16g}"
        dest = Path(tempfile.mkdtemp(prefix=f"day554-{kind}-", dir=self.tmp))
        repo = dest / "repo"
        for rel in REL538 + EXTRA:
            (repo / rel).parent.mkdir(parents=True, exist_ok=True)
            shutil.copy(self.tree / rel if kind == "554" else self.base / rel, repo / rel)
        h._stub_tree(repo)
        bindir = dest / "bin"                                       # 538's stub uv + one line: write the dump it is given
        bindir.mkdir()
        body = T538.UV_BODY.replace(
            '        h5.write_bytes(b"h5")\n',
            '        h5.write_bytes(b"h5")\n'
            '        if os.environ.get("VCC_RC_AUDIT_DUMP"):\n'
            '            Path(os.environ["VCC_RC_AUDIT_DUMP"]).write_text(json.dumps({"context": "D", "build": name}) + "\\n")\n', 1)
        assert body != T538.UV_BODY
        (bindir / "uv").write_text("#!" + str(PY) + "\n" + body)
        (bindir / "uv").chmod(0o755)
        if kind != "554":
            day = repo / "scripts/def_day.sh"                       # the 538 day derives REPO from its path, as 554's
            assert "REPO=\"$(CDPATH= cd" in day.read_text()
        return h._run_bash(dest, repo, env)

    @staticmethod
    def _strip(records):
        out = []
        for r in records:
            a = list(r["argv"])
            if "--rc-audit" in a:
                i = a.index("--rc-audit")
                del a[i:i + 2]
            v = {k: x for k, x in r["vcc"].items() if k != "VCC_RC_AUDIT_DUMP"}
            out.append({"argv": a, "vcc": v})
        return out

    def test_C_day_unset_and_overrides(self):
        for env in ({}, {"DEF_DAY_LODO": "1", "DEF_DAY_DE_MAG": "dn_k250_zs"}):
            b, p = self._day("538", env), self._day("554", env)
            self.assertEqual(b["rc"], 0, b["log"][-2000:] + b["err"][-1000:])
            self.assertEqual(p["rc"], 0, p["log"][-2000:] + p["err"][-1000:])
            self.assertEqual(self._strip(b["records"]), self._strip(p["records"]))
            self.assertEqual(b["final"], p["final"])
            builds = [r for r in p["records"] if any(x.endswith("def_pipeline.py") for x in r["argv"])
                      and "--stages" in r["argv"] and r["argv"][r["argv"].index("--stages") + 1] == "build"]
            for r in builds:
                name = r["argv"][r["argv"].index("--name") + 1]
                if name.endswith("-dump"):
                    self.assertNotIn("VCC_RC_AUDIT_DUMP", r["vcc"], name)
                else:
                    self.assertEqual(r["vcc"].get("VCC_RC_AUDIT_DUMP"), f"$W/{name}.rc_audit.jsonl", name)
            post = [r for r in p["records"] if any(x.endswith("def_postbuild.py") for x in r["argv"])]
            self.assertEqual(len(post), 1)
            a = post[0]["argv"]
            self.assertEqual(a[a.index("--rc-audit") + 1], "$W/day_final.rc_audit.jsonl")
            w = p["work"]
            self.assertTrue((w / "day_final.rc_audit.jsonl").is_file())
            final_from = json.loads((w / "day_final.rc_audit.jsonl").read_text())["build"]
            left = sorted(f.name for f in w.glob("day-*.rc_audit.jsonl"))
            self.assertEqual(left, [f"{final_from}.rc_audit.jsonl"])  # non-final dumps pruned with their h5ads
            self.assertNotIn("554:", p["log"])
        o = self._day("554", {"DEF_DAY_ALLOW_WIDE_LIBS": "1", "DEF_DAY_ALLOW_LOW_DOSE": "1"})
        self.assertEqual(o["rc"], 0, o["log"][-2000:])
        self.assertIn("DEF_DAY_ALLOW_WIDE_LIBS=1\n", o["final"])
        self.assertIn("DEF_DAY_ALLOW_LOW_DOSE=1\n", o["final"])
        self.assertIn("554: DEF_DAY_ALLOW_WIDE_LIBS=1", o["log"])
        self.assertIn("554: DEF_DAY_ALLOW_LOW_DOSE=1", o["log"])


if __name__ == "__main__":
    unittest.main()
