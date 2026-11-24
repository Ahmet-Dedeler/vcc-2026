#!/usr/bin/env python3
"""644 -- tests for DEF_DAY_SAMELINE_SKIP (patch_skip644.py v2: def_sameline --skip-contexts + the def_day pass-through).

Trees are 606's 10-06 chain: test598.canon_tree (HEAD + 538 -> 554 -> 590), then 641 -> 598 -> 612 -> 612dd -> 630 -> 642 ->
604 (compose606.PATCHES), with data/, submissions/ and the other experiments symlinked read-only; the skip patch goes last.
The live files are hashed at start and end.

  M  mechanics: --check clean, writes nothing; apply; second apply refused; a tree without 642 refused; an anchor miss.
  a  UNSET: on 538's stub harness the patched day makes the same calls (argv + VCC_*) and writes the same
     def_day.final.env; its def_sameline call carries no --skip-contexts; the real def_sameline.py without the flag
     writes byte-identical artifacts, with no programme source and with --programme-source orion_hct116.
  b  SKIP: the stub day passes --skip-contexts E, logs and records it; a non-manifest name stops the day. The real
     def_sameline.py on an identity that matches all three contexts (JURKAT / K562 / HCT116; unanimous, 0.75):
       no source: the kept contexts' arrays and provenance are byte-identical, and the skipped one is absent;
       --programme-source orion_hct116: the kept contexts are byte-identical, and the skipped one has exactly an
         unmatched context's form (the source's programme, w 0, no rows, no signs);
       a 642 HEK293 call on the skipped context (knob + decision set): the skip wins.
  c  builder 4daae237 (no build): its own sameline_load (ast-extracted) returns exactly the kept contexts for the
     no-source artifact (.get(skipped) is None: no KeyError), and for the source artifact a skipped record with no
     targets and w 0. Every same-line effect sits inside `if _slr is not None:`, whose else is the unmatched path.

    uv run python experiments/644-runbook-amendments/test_skip644.py -v
"""
from __future__ import annotations

import ast
import hashlib
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
for d in ("", "../606-applyday-composition", "../598-tercilec-knob", "../554-dose-guard", "../538-cloud-bundle-extras",
          "../642-hek293-contingency"):
    sys.path.insert(0, str((HERE / d).resolve()))
import patch_skip644 as S  # noqa: E402
import compose606 as C6  # noqa: E402
import test598 as T598  # noqa: E402
import test_guard554 as G  # noqa: E402
import test_extras538 as T538  # noqa: E402

PY = str(ROOT / ".venv/bin/python")
PATCH = HERE / "patch_skip644.py"
BUILDER = ROOT / "scripts/make_submission_pm_cd4.py"
LIVE = sorted({*G.LIVE, ROOT / "scripts/def_day.sh", ROOT / "scripts/def_sameline.py", BUILDER})
CALL3 = {"call": {"A": "JURKAT_HAEMATOPOIETIC_AND_LYMPHOID_TISSUE", "B": "K562_HAEMATOPOIETIC_AND_LYMPHOID_TISSUE",
                  "C": "HCT116_LARGE_INTESTINE"},
         "call_unanimous_across_methods": {"A": True, "B": True, "C": True},
         "min_top1_score": {"A": 0.75, "B": 0.75, "C": 0.75}}
RESULTS: dict = {}


def sha(p: Path) -> str:
    return hashlib.sha256(p.read_bytes()).hexdigest()


def live_hashes() -> dict:
    return {str(p): sha(p) for p in LIVE}


def clean_env(**kw) -> dict:
    env = {k: v for k, v in os.environ.items() if not k.startswith(("DEF_DAY_", "VCC_"))}
    env.update(kw)
    return env


def chain_tree(dest: Path, upto: str = "604") -> Path:
    T598.canon_tree(dest)
    (dest / "data/identity_refs").mkdir(parents=True, exist_ok=True)
    shutil.copy2(ROOT / "data/identity_refs/hek293_nonflex_ref642.npz", dest / "data/identity_refs/hek293_nonflex_ref642.npz")
    (dest / "src/vcc2026").mkdir(parents=True, exist_ok=True)
    for f in (ROOT / "src/vcc2026").glob("*.py"):
        shutil.copy2(f, dest / "src/vcc2026" / f.name)
    for name in C6.DEFAULT_ORDER[: C6.DEFAULT_ORDER.index(upto) + 1]:
        script = C6.PATCHES[name][0]
        c = subprocess.run([PY, str(ROOT / script), "--root", str(dest), "--check"], capture_output=True, text=True)
        co = (c.stdout + c.stderr).lower()
        if c.returncode != 0 and "already" in co and "partially" not in co:
            continue   # B2 ([4d0207] for the coordinator): already on HEAD after the 10-05 / 10-06 applies
        r = subprocess.run([PY, str(ROOT / script), "--root", str(dest)], capture_output=True, text=True)
        assert r.returncode == 0, (name, r.stdout[-1500:], r.stderr[-1500:])
    for top in ("data", "experiments"):
        (dest / top).mkdir(exist_ok=True)
        for e in sorted((ROOT / top).iterdir()):
            if not (dest / top / e.name).exists():
                os.symlink(e, dest / top / e.name)
    if not (dest / "submissions").exists():
        os.symlink(ROOT / "submissions", dest / "submissions")
    return dest


class Skip644(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        if S.MARK in (ROOT / "scripts/def_sameline.py").read_text():
            raise unittest.SkipTest("644 (skip) is applied on HEAD (10-05 / 10-06 apply days): this suite validates the patch against a PRE-644 (skip) tree, which cannot be built from live files any more. Run it in a worktree at the pre-apply commit with the gitignored inputs linked (606 RESULTS, the P1/P2 worktree proof); post-apply regression = experiments/644-runbook-amendments/verify644.py")
        cls.live = live_hashes()
        cls.tmp = Path(tempfile.mkdtemp(prefix="test644s-"))
        cls.base = chain_tree(cls.tmp / "base")
        cls.tree = chain_tree(cls.tmp / "patched")
        r = subprocess.run([PY, str(PATCH), "--root", str(cls.tree)], capture_output=True, text=True)
        assert r.returncode == 0, r.stdout + r.stderr
        cls.shas = {str(f): sha(cls.tree / f) for f in S.FILES}

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.tmp, ignore_errors=True)
        assert live_hashes() == cls.live, "a live file changed"
        (HERE / "work").mkdir(exist_ok=True)
        RESULTS["patched_shas"] = cls.shas
        (HERE / "work/test_skip644.json").write_text(json.dumps(RESULTS, indent=1, default=str))

    # ------------------------------------------------------------------------------------------ helpers
    def _day(self, src: Path, env: dict) -> dict:
        h = T538.Extras538("_assert_day_ok")
        h.star = f"{T538.rho_star():.16g}"
        dest = Path(tempfile.mkdtemp(prefix="day644-", dir=self.tmp))
        repo = dest / "repo"
        for rel in G.REL538 + G.EXTRA:
            (repo / rel).parent.mkdir(parents=True, exist_ok=True)
            shutil.copy(src / rel, repo / rel)
        h._stub_tree(repo)
        bindir = dest / "bin"
        bindir.mkdir()
        body = T538.UV_BODY.replace(
            '        h5.write_bytes(b"h5")\n',
            '        h5.write_bytes(b"h5")\n'
            '        if os.environ.get("VCC_RC_AUDIT_DUMP"):\n'
            '            Path(os.environ["VCC_RC_AUDIT_DUMP"]).write_text(json.dumps({"context": "D", "build": name}) + "\\n")\n', 1)
        (bindir / "uv").write_text("#!" + PY + "\n" + body)
        (bindir / "uv").chmod(0o755)
        return h._run_bash(dest, repo, env)

    @staticmethod
    def _sl_calls(records):
        return [r["argv"] for r in records if any(x.endswith("def_sameline.py") for x in r["argv"]) and "--identity" in r["argv"]]

    def _sl(self, tree: Path, ident: dict, out_name: str, extra=(), env=None, ident_name="ident.json"):
        ip = self.tmp / ident_name               # the artifact's provenance embeds this path: keep it fixed for byte-identity
        ip.write_text(json.dumps(ident))
        out = self.tmp / out_name
        p = subprocess.run([PY, str(tree / "scripts/def_sameline.py"), "--bundle", str(ROOT / "data/challenge"),
                            "--identity", str(ip), "--out", str(out), *extra], cwd=tree, env=env or clean_env(),
                           capture_output=True, text=True)
        self.assertEqual(p.returncode, 0, p.stdout[-1500:] + p.stderr[-1500:])
        z = np.load(out, allow_pickle=False)
        prov = json.loads(str(z["provenance"]))
        prov.pop("identity_source", None)
        return out, {k: z[k] for k in z.files if k != "provenance"}, prov, p.stdout

    @staticmethod
    def _eq(x, y):
        return np.array_equal(x, y, equal_nan=True) if np.issubdtype(np.asarray(x).dtype, np.floating) else np.array_equal(x, y)

    # ------------------------------------------------------------------------------------------ M
    def test_M_mechanics(self):
        t = chain_tree(self.tmp / "mech")
        b = {str(f): sha(t / f) for f in S.FILES}
        r = subprocess.run([PY, str(PATCH), "--root", str(t), "--check"], capture_output=True, text=True)
        self.assertEqual(r.returncode, 0, r.stdout + r.stderr)
        self.assertEqual(b, {str(f): sha(t / f) for f in S.FILES})
        self.assertEqual(subprocess.run([PY, str(PATCH), "--root", str(t)], capture_output=True).returncode, 0)
        r2 = subprocess.run([PY, str(PATCH), "--root", str(t)], capture_output=True, text=True)
        self.assertNotEqual(r2.returncode, 0)
        self.assertIn("refusing a second apply", r2.stdout + r2.stderr)
        self.assertEqual(subprocess.run(["bash", "-n", str(t / "scripts/def_day.sh")]).returncode, 0)
        pre = chain_tree(self.tmp / "pre642", upto="630")
        bp = {str(f): sha(pre / f) for f in S.FILES}
        r3 = subprocess.run([PY, str(PATCH), "--root", str(pre)], capture_output=True, text=True)
        self.assertNotEqual(r3.returncode, 0)
        self.assertIn("applies after 642", r3.stdout + r3.stderr)
        self.assertEqual(bp, {str(f): sha(pre / f) for f in S.FILES})
        tm = chain_tree(self.tmp / "miss")
        f = tm / S.FILES[0]
        f.write_text(f.read_text().replace('    ap.add_argument("--min-programme-rows", type=int, default=150,\n',
                                           '    ap.add_argument("--min-programme-rows", type=int,  default=150,\n'))
        bm = {str(x): sha(tm / x) for x in S.FILES}
        r4 = subprocess.run([PY, str(PATCH), "--root", str(tm)], capture_output=True, text=True)
        self.assertNotEqual(r4.returncode, 0)
        self.assertIn("anchor found 0 times", r4.stdout + r4.stderr)
        self.assertEqual(bm, {str(x): sha(tm / x) for x in S.FILES})
        RESULTS["M"] = "PASS"

    # ------------------------------------------------------------------------------------------ a
    def test_a_unset(self):
        rec = {}
        for env in ({}, {"DEF_DAY_LODO": "1", "DEF_DAY_DE_MAG": "dn_k250_zs"}):
            b, p = self._day(self.base, env), self._day(self.tree, env)
            self.assertEqual(b["rc"], 0, b["log"][-2000:] + b["err"][-1000:])
            self.assertEqual(p["rc"], 0, p["log"][-2000:] + p["err"][-1000:])
            self.assertEqual(b["records"], p["records"])
            self.assertEqual(b["final"], p["final"])
            calls = self._sl_calls(p["records"])
            self.assertEqual(len(calls), 1)
            self.assertNotIn("--skip-contexts", calls[0])
            self.assertNotIn("644:", p["log"])
            rec[json.dumps(env)] = len(p["records"])
        for tag, extra in (("nosrc", ()), ("src", ("--programme-source", "orion_hct116"))):
            ob, *_ = self._sl(self.base, CALL3, f"a_{tag}_base.npz", extra)
            op, *_ = self._sl(self.tree, CALL3, f"a_{tag}_pat.npz", extra)
            self.assertEqual(sha(ob), sha(op), tag)
            rec[f"artifact_{tag}"] = sha(op)
        RESULTS["a"] = rec

    # ------------------------------------------------------------------------------------------ b
    def test_b_skip(self):
        rec = {}
        s = self._day(self.tree, {"DEF_DAY_SAMELINE_SKIP": "E"})
        self.assertEqual(s["rc"], 0, s["log"][-2000:])
        calls = self._sl_calls(s["records"])
        self.assertEqual(calls[0][calls[0].index("--skip-contexts") + 1], "E")
        self.assertIn("644: same-line skipped for E", s["log"])
        self.assertIn("644: DEF_DAY_SAMELINE_SKIP=E", s["log"])
        self.assertIn("DEF_DAY_SAMELINE_SKIP=E\n", s["final"])
        sp = self._day(self.tree, {"DEF_DAY_SAMELINE_SKIP": "E, F"})              # 607: whitespace is ignored
        self.assertEqual(sp["rc"], 0, sp["log"][-2000:])
        spc = self._sl_calls(sp["records"])
        self.assertEqual(spc[0][spc[0].index("--skip-contexts") + 1], "E,F")
        self.assertIn("644: same-line skipped for E,F (def_sameline --skip-contexts): treated as unmatched", sp["log"])
        bad = self._day(self.tree, {"DEF_DAY_SAMELINE_SKIP": "Q"})
        self.assertNotEqual(bad["rc"], 0)
        self.assertIn("644: DEF_DAY_SAMELINE_SKIP names Q, which is not a manifest context", bad["log"])
        # real def_sameline, no programme source
        _, full, pf, _ = self._sl(self.tree, CALL3, "b_full.npz")
        out, kept, pk, so = self._sl(self.tree, CALL3, "b_skipB.npz", ("--skip-contexts", "B"))
        self.assertEqual(sorted(map(str, full["contexts"])), ["A", "B", "C"])
        self.assertEqual(sorted(map(str, kept["contexts"])), ["A", "C"])
        for c in ("A", "C"):
            for k in [k for k in full if k.endswith(f"_{c}")]:
                self.assertTrue(self._eq(full[k], kept[k]), (c, k))
            self.assertEqual(pf["contexts"][c], pk["contexts"][c])
        self.assertFalse(any(k.endswith("_B") for k in kept))
        self.assertIn("644", pk["contexts"]["B"]["mode"])
        self.assertIn("B: K562_HAEMATOPOIETIC_AND_LYMPHOID_TISSUE skipped (--skip-contexts): treated as unmatched (644)", so)
        _, kws, _, _ = self._sl(self.tree, CALL3, "b_skipBC_ws.npz", ("--skip-contexts", "B, C"))   # 607: stripped names
        self.assertEqual(sorted(map(str, kws["contexts"])), ["A"])
        for k in [k for k in full if k.endswith("_A")]:
            self.assertTrue(self._eq(full[k], kws[k]), k)
        self.__class__.art_nosrc = out
        rec["nosrc_covered_full"] = {c: pf["contexts"][c].get("covered_targets") for c in "ABC"}
        # with a programme source: the skipped context has an unmatched context's form
        src = ("--programme-source", "orion_hct116")
        _, fulls, pfs, _ = self._sl(self.tree, CALL3, "b_full_src.npz", src)
        outs, keps, pks, _ = self._sl(self.tree, CALL3, "b_skipB_src.npz", (*src, "--skip-contexts", "B"))
        for c in ("A", "C"):
            for k in [k for k in fulls if k.endswith(f"_{c}")]:
                self.assertTrue(self._eq(fulls[k], keps[k]), (c, k))
        self.assertEqual(float(keps["w_B"]), 0.0)
        self.assertEqual(keps["targets_B"].size, 0)
        self.assertEqual(keps["rows_B"].shape[0], 0)
        self.assertEqual(keps["signs_B"].shape[0], 0)
        self.assertTrue(np.array_equal(keps["programme_B"], keps["programme_A"]))   # the source's programme, as for any context
        self.assertTrue(bool(keps["programme_ok_B"]))
        self.assertIn("programme_source", pks["contexts"]["B"])
        self.assertIn("644", pks["contexts"]["B"]["mode"])
        self.__class__.art_src = outs
        # 642 on the skipped context: the skip wins even with the knob and a decision
        ident = json.loads(json.dumps(CALL3))
        per = {c: {"call": "X", "unanimous": False, "n_methods_agree": 0, "min_top1_score": 0.5, "top5": [], "row_top1_methods": 0,
                   "accepted": False} for c in "ABC"}
        per["B"] = {"call": "HEK293_REF", "unanimous": True, "n_methods_agree": 8, "min_top1_score": 0.82, "top5": [],
                    "row_top1_methods": 8, "accepted": True}
        ident["nonccle_ref"] = {"reference": "hek293", "sha256": "x", "per_context": per,
                                "call": {"A": None, "B": "HEK293", "C": None}}
        _, k642, p642, _ = self._sl(self.tree, ident, "b_642.npz", ("--skip-contexts", "B"),
                                    env=clean_env(DEF_DAY_ALLOW_NONCCLE_REF="hek293", DEF_DAY_NONCCLE_DECISION="test"),
                                    ident_name="ident642.json")
        self.assertFalse(any(k.endswith("_B") for k in k642))
        self.assertIn("644", p642["contexts"]["B"]["mode"])
        rec["skip_vs_642"] = p642["contexts"]["B"]["mode"]
        RESULTS["b"] = rec

    # ------------------------------------------------------------------------------------------ c
    def test_c_builder(self):
        if not hasattr(self.__class__, "art_src"):
            self.test_b_skip()
        src = BUILDER.read_text()
        tree = ast.parse(src)
        fn = next(n for n in tree.body if isinstance(n, ast.FunctionDef) and n.name == "sameline_load")
        res = {}
        for tag, art in (("nosrc", self.art_nosrc), ("src", self.art_src)):
            ns = {"np": np, "Path": Path, "SAMELINE": {"path": str(art), "sha256": sha(art), "w": None}}
            exec(compile(ast.Module(body=[fn], type_ignores=[]), str(BUILDER), "exec"), ns)
            z = np.load(art)
            cache = ns["sameline_load"](list(map(str, z["genes"])), list(map(str, z["perts"])))
            if tag == "nosrc":
                self.assertEqual(sorted(cache), ["A", "C"])
                self.assertIsNone(cache.get("B"))
            else:
                self.assertEqual(sorted(cache), ["A", "B", "C"])
                self.assertEqual(len(cache["B"]["targets"]), 0)
                self.assertEqual(float(cache["B"]["w"]), 0.0)
            res[tag] = sorted(cache)
        found = [n for n in ast.walk(tree) if isinstance(n, ast.If) and isinstance(n.test, ast.Compare)
                 and isinstance(n.test.left, ast.Name) and n.test.left.id == "_slr" and isinstance(n.test.ops[0], ast.IsNot)]
        self.assertEqual(len(found), 1)
        blk = found[0]
        inside = {id(n) for n in ast.walk(ast.Module(body=blk.body, type_ignores=[]))}
        effects = []
        for node in ast.walk(tree):
            if isinstance(node, ast.Call) and getattr(node.func, "id", "") == "sameline_blend":
                effects.append(("sameline_blend", node))
            if isinstance(node, (ast.Assign, ast.AugAssign)):
                for t in (node.targets if isinstance(node, ast.Assign) else [node.target]):
                    if isinstance(t, ast.Subscript) and getattr(t.value, "id", "") in ("sl_prior", "SAMELINE_AUDIT"):
                        effects.append((t.value.id, node))
        bad = [(k, n.lineno) for k, n in effects if id(n) not in inside and k in ("sameline_blend", "SAMELINE_AUDIT")]
        self.assertEqual(bad, [])
        self.assertIn("not in the artifact -- unchanged", ast.unparse(ast.Module(body=blk.orelse, type_ignores=[])))
        seg = src.splitlines()[blk.lineno - 12:blk.lineno]
        self.assertTrue(any("sl_prior: dict[str, np.ndarray] = {}" in ln for ln in seg))
        # the pt loop iterates the record's own targets: an empty target list adds no per-target sign
        pt_src = ast.get_source_segment(src, blk)
        self.assertIn('for _k, _ti in enumerate(np.asarray(_slr["targets"], dtype=np.int64)):', pt_src)
        RESULTS["c"] = dict(builder_sha256=sha(BUILDER), branch=[blk.lineno, blk.end_lineno], sameline_load=res)


if __name__ == "__main__":
    unittest.main()
