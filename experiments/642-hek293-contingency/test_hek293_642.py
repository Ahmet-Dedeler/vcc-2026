#!/usr/bin/env python3
"""642 -- tests for the HEK293 contingency patch. Temp POST-BUNDLE trees only (tree642.make_tree: copies + 538's applier +
554 + 590 + 641); the live files are hashed at the start and the end and must not change.

  T1 mechanics: --check clean and writes nothing; apply; second apply refused; a pre-590/641 tree refused; a partial tree
     refused; an anchor miss writes nothing; a modified reference file refused.
  T2 def_identity CLI on data/challenge: knob OFF -> --verify PASS and the --out JSON equals the unpatched post-bundle
     tree's (timing masked); knob ON -> every stage-1 key equal, `nonccle_ref` added, A/B/C not called HEK293.
  T3 the coordinator's pins, through the PATCHED module's nonccle_stage over 590's 423 labelled profiles + 597's three
     real-Flex profiles: (1) stage 1 (CCLE only) reproduces 590 on all 423 (0 mismatches), knob off and on alike (stage 2
     never touches stage 1); (2) 0 of the 423 non-HEK293 profiles called HEK293; (3) the held-out real-Flex HEK293 called
     HEK293, 8/8, min_top1 >= 0.65.
  T4 the Flex guard: a reference whose provenance lists a Flex source is refused (and an unpinned sha, and no sources).
  T5 def_sameline (bundle data/challenge, the real stage-1 identity + a synthetic second-stage block): no block -> npz
     and provenance identical to the unpatched tree; HEK293 on B with the knob unset / no decision / the wrong knob ->
     B unchanged; knob + decision -> B reads orion_hek293t (w 2, decision recorded), A untouched; HEK293 on A (a
     registered, unanimous CCLE call) -> the CCLE call is vetoed; --force untouched.
  T6 the def_day.sh lines (bash with a stub log): logged and recorded only when set; a decision with spaces round-trips.

    uv run python experiments/642-hek293-contingency/test_hek293_642.py -v       (T3 takes ~5 min, <= 3.6 GB)
"""
from __future__ import annotations

import importlib.util
import json
import os
import re
import shutil
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(ROOT / "experiments/554-dose-guard"))
import tree642 as TR  # noqa: E402
import patch_hek293_642 as P  # noqa: E402
import test_guard554 as G  # noqa: E402

PY = TR.PY
PATCH = HERE / "patch_hek293_642.py"
REF = ROOT / "data/identity_refs/hek293_nonflex_ref642.npz"
LIVE = sorted({*G.LIVE, ROOT / "scripts/def_identity.py", ROOT / "scripts/def_sameline.py", ROOT / "scripts/def_day.sh", REF})
RESULTS: dict = {}


def live_hashes() -> dict:
    return {str(p): G.sha(p) for p in LIVE}


def clean_env(**extra) -> dict:
    env = {k: v for k, v in os.environ.items() if not k.startswith(("DEF_DAY_", "VCC_"))}
    env.update(extra)
    return env


def apply642(tree: Path) -> subprocess.CompletedProcess:
    return subprocess.run([str(PY), str(PATCH), "--root", str(tree)], capture_output=True, text=True)


def load_module(path: Path, name: str):
    spec = importlib.util.spec_from_file_location(name, path)
    m = importlib.util.module_from_spec(spec)
    sys.modules[name] = m
    spec.loader.exec_module(m)
    return m


def strip_timing(d: dict) -> dict:
    d = json.loads(json.dumps(d))
    d.pop("timing_seconds", None)
    d.pop("nonccle_ref", None)
    return d


class Hek642(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        if P.MARK in (ROOT / "scripts/def_identity.py").read_text():
            raise unittest.SkipTest("642 is applied on HEAD (10-05 / 10-06 apply days): this suite validates the patch against a PRE-642 tree, which cannot be built from live files any more. Run it in a worktree at the pre-apply commit with the gitignored inputs linked (606 RESULTS, the P1/P2 worktree proof); post-apply regression = experiments/644-runbook-amendments/verify644.py")
        cls.live = live_hashes()
        cls.tmp = Path(tempfile.mkdtemp(prefix="test642-"))
        cls.base = TR.make_tree(cls.tmp / "base")          # post-bundle, 642 not applied
        cls.tree = TR.make_tree(cls.tmp / "patched")
        r = apply642(cls.tree)
        assert r.returncode == 0, r.stdout + r.stderr

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.tmp, ignore_errors=True)
        assert live_hashes() == cls.live, "a live file changed"
        (HERE / "work").mkdir(exist_ok=True)
        (HERE / "work/test642.json").write_text(json.dumps(RESULTS, indent=1, default=str))

    # ------------------------------------------------------------------------------------------------- T1
    def test_T1_mechanics(self):
        t = TR.make_tree(self.tmp / "mech")
        before = {f: G.sha(t / f) for f in P.FILES}
        r = subprocess.run([str(PY), str(PATCH), "--root", str(t), "--check"], capture_output=True, text=True)
        self.assertEqual(r.returncode, 0, r.stderr)
        self.assertEqual(before, {f: G.sha(t / f) for f in P.FILES})
        self.assertEqual(apply642(t).returncode, 0)
        after = {f: G.sha(t / f) for f in P.FILES}
        self.assertTrue(all(after[f] != before[f] for f in P.FILES))
        r2 = apply642(t)
        self.assertNotEqual(r2.returncode, 0)
        self.assertIn("refusing a second apply", r2.stdout + r2.stderr)
        self.assertEqual(subprocess.run(["bash", "-n", str(t / "scripts/def_day.sh")]).returncode, 0)
        # pre-590/641: 554's tree (538 + 554 only) plus def_identity copied unpatched. B2 ([4d0207] for the coordinator):
        # after the 10-05 apply the LIVE files already carry 590 / 641, so a pre-590/641 tree cannot be built from them;
        # this one sub-check (a property of the patch, verified on the pre-10-05 HEAD) is then skipped, the rest of T1 runs.
        live641 = subprocess.run([str(PY), str(ROOT / "experiments/641-depth-row/patch_depth641.py"), "--root", str(ROOT), "--check"],
                                 capture_output=True, text=True)
        if live641.returncode != 0 and "already" in (live641.stdout + live641.stderr).lower():
            print("T1: live tree already carries 641 (post-10-05): pre-590/641 refusal sub-check skipped", flush=True)
        else:
            t0 = G.make_tree(self.tmp / "pre", with554=True)
            shutil.copy(ROOT / "scripts/def_identity.py", t0 / "scripts/def_identity.py")
            b0 = {f: G.sha(t0 / f) for f in P.FILES}
            r3 = apply642(t0)
            self.assertNotEqual(r3.returncode, 0)
            self.assertIn("applies after 538 + 554 + 590 + 641", r3.stdout + r3.stderr)
            self.assertEqual(b0, {f: G.sha(t0 / f) for f in P.FILES})
        # partial
        tp = TR.make_tree(self.tmp / "partial")
        new = P.plan(tp)
        (tp / P.FILES[1]).write_text(new[P.FILES[1]])
        bp = {f: G.sha(tp / f) for f in P.FILES}
        r4 = apply642(tp)
        self.assertNotEqual(r4.returncode, 0)
        self.assertIn("partially patched", r4.stdout + r4.stderr)
        self.assertEqual(bp, {f: G.sha(tp / f) for f in P.FILES})
        # anchor miss in def_day.sh: nothing written anywhere
        tm = TR.make_tree(self.tmp / "miss")
        (tm / P.FILES[2]).write_text((tm / P.FILES[2]).read_text().replace(P.DAY_ENV_A, P.DAY_ENV_A.replace("=1\\n", "=1 \\n")))
        bm = {f: G.sha(tm / f) for f in P.FILES}
        r5 = apply642(tm)
        self.assertNotEqual(r5.returncode, 0)
        self.assertIn("anchor found 0 times", r5.stdout + r5.stderr)
        self.assertEqual(bm, {f: G.sha(tm / f) for f in P.FILES})
        # a modified reference: a tree whose data/identity_refs holds a changed copy
        tr = TR.make_tree(self.tmp / "badref")
        os.unlink(tr / "data")
        (tr / "data/identity_refs").mkdir(parents=True)
        bad = tr / "data/identity_refs/hek293_nonflex_ref642.npz"
        shutil.copy(REF, bad)
        with open(bad, "ab") as fh:
            fh.write(b"\0")
        r6 = apply642(tr)
        self.assertNotEqual(r6.returncode, 0)
        self.assertIn("not the pinned reference", r6.stdout + r6.stderr)
        RESULTS["T1"] = "PASS"

    # ------------------------------------------------------------------------------------------------- T2
    def test_T2_identity_cli(self):
        rec = {}
        outs = {}
        for name, tree, env in (("base", self.base, clean_env()), ("off", self.tree, clean_env()),
                                ("on", self.tree, clean_env(DEF_DAY_ALLOW_NONCCLE_REF="hek293"))):
            o = self.tmp / f"ident_{name}.json"
            p = subprocess.run([str(PY), str(tree / "scripts/def_identity.py"), "--bundle", str(ROOT / "data/challenge"),
                                "--verify", "--out", str(o)], cwd=tree, env=env, capture_output=True, text=True)
            self.assertEqual(p.returncode, 0, p.stdout[-1500:] + p.stderr[-1500:])
            self.assertIn("PASS", p.stdout)
            outs[name] = json.loads(o.read_text())
            rec[name] = {"rc": p.returncode, "has_nonccle": "nonccle_ref" in outs[name],
                         "642_lines": [ln for ln in p.stdout.splitlines() if ln.startswith("[642]")]}
        self.assertNotIn("nonccle_ref", outs["off"])
        self.assertEqual(strip_timing(outs["base"]), strip_timing(outs["off"]))
        self.assertEqual(strip_timing(outs["base"]), strip_timing(outs["on"]))       # stage 1 unchanged with the knob on
        n = outs["on"]["nonccle_ref"]
        self.assertEqual(n["sha256"], P.REF_SHA256)
        self.assertEqual(n["call"], {"A": None, "B": None, "C": None})
        rec["on_per_context"] = {c: {k: v for k, v in r.items() if k != "top5"} for c, r in n["per_context"].items()}
        RESULTS["T2"] = rec
        self.__class__.ident_abc = outs["base"]

    # ------------------------------------------------------------------------------------------------- T3
    def test_T3_pins(self):
        import pandas as pd
        DIp = load_module(self.tree / "scripts/def_identity.py", "def_identity_642p")
        sys.path.insert(0, str(ROOT / "experiments/590-identity-validation"))
        import match642 as MM
        Xfull, ccle_genes, samples, clinical = DIp.load_cache()
        official = pd.read_csv(ROOT / "data/challenge/gene_names.csv")["gene_name"].astype(str).tolist()
        g2c = {g: i for i, g in enumerate(ccle_genes)}
        common = [g for g in official if g in g2c]
        X = Xfull[:, np.array([g2c[g] for g in common], dtype=np.int64)]
        genes = np.asarray(common)
        del Xfull
        queries, info = MM.build_queries(genes, official)
        variance, mean = np.var(X, axis=0), np.mean(X, axis=0)
        non_hk = np.array([not re.match(r"^(RPL|RPS|MT-|HIST|H[12]B)", g) for g in genes])
        fsets = {"all": np.arange(len(genes))}
        for k in (2000, 5000, 10000):
            el = np.flatnonzero(non_hk & (mean > 0.25))
            fsets[f"hvg{k}"] = el[np.argsort(-variance[el])[:k]]
        base, qual = MM.score_all(X, samples, queries, fsets)
        ref590 = {r["query"]: r for r in json.loads((ROOT / "experiments/590-identity-validation/work/match590.json").read_text())["rows"]}
        common_q = [q for q in queries if q in ref590]
        mism = [q for q in common_q if (ref590[q]["call"] != base[q]["call"] or ref590[q]["unanimous"] != base[q]["unanimous"]
                                        or abs(ref590[q]["min_s1"] - base[q]["min_s1"]) > 1e-9)]
        self.assertEqual(len(common_q), 423)
        self.assertEqual(mism, [])                                                          # pin (1), stage 1
        st2 = DIp.nonccle_stage("hek293", X, genes, samples, fsets, qual, queries, list(queries))
        hek = {q for q in queries if MM.hek_source(q, info[q])}
        self.assertEqual(hek, {"donor:orion_hek293t", "xu:HEK293", "pilot:HEK293"})
        false = [q for q in queries if q not in hek and st2["per_context"][q]["accepted"]]
        self.assertEqual(false, [])                                                         # pin (2)
        h = st2["per_context"]["pilot:HEK293"]
        self.assertTrue(h["accepted"])                                                      # pin (3)
        self.assertEqual(h["n_methods_agree"], 8)
        self.assertEqual(DIp.NONCCLE_FLOOR, 0.75)
        self.assertGreaterEqual(h["min_top1_score"], DIp.NONCCLE_FLOOR)
        # the coordinator's 0.75: every non-HEK293 profile with the row top-1 in >= 5/8 methods is below it (blocked twice)
        strong_near = [(q, st2["per_context"][q]["min_top1_score"]) for q in queries
                       if q not in hek and st2["per_context"][q]["row_top1_methods"] >= 5]
        self.assertTrue(strong_near and all(v < DIp.NONCCLE_FLOOR for _, v in strong_near), strong_near)
        near = sorted(((q, st2["per_context"][q]["row_top1_methods"], round(st2["per_context"][q]["min_top1_score"], 4))
                       for q in queries if q not in hek and st2["per_context"][q]["row_top1_methods"] > 0), key=lambda x: -x[1])
        RESULTS["T3"] = dict(stage1_common=len(common_q), stage1_mismatches=mism, false_calls=false,
                             heldout_flex_hek293={k: h[k] for k in ("call", "n_methods_agree", "min_top1_score", "accepted")},
                             in_sample_sources={q: {k: st2["per_context"][q][k] for k in ("call", "n_methods_agree", "min_top1_score")}
                                                for q in ("donor:orion_hek293t", "xu:HEK293")},
                             near_misses_row_top1=near)

    # ------------------------------------------------------------------------------------------------- T4
    def test_T4_flex_guard(self):
        DIp = load_module(self.tree / "scripts/def_identity.py", "def_identity_642g")
        z = np.load(REF)
        prov = json.loads(str(z["provenance"]))
        res = {}
        for tag, srcs in (("flex", prov["sources"] + [{"name": "pilot", "chemistry": "10x Flex v1"}]), ("nosrc", [])):
            p = self.tmp / f"ref_{tag}.npz"
            np.savez(p, genes=z["genes"], values=z["values"], provenance=np.asarray(json.dumps({**prov, "sources": srcs})))
            DIp.NONCCLE_REFS["t_" + tag] = (str(p), DIp.sha256(p))
            with self.assertRaises(SystemExit) as cm:
                DIp.load_nonccle_ref("t_" + tag)
            self.assertIn("refusing", str(cm.exception))
            res[tag] = str(cm.exception)
        DIp.NONCCLE_REFS["t_sha"] = (str(REF), "0" * 64)
        with self.assertRaises(SystemExit) as cm:
            DIp.load_nonccle_ref("t_sha")
        self.assertIn("not the pinned", str(cm.exception))
        res["sha"] = str(cm.exception)
        g, v, pv = DIp.load_nonccle_ref("hek293")                                           # the real one loads
        self.assertEqual(len(g), 18533)
        RESULTS["T4"] = res

    # ------------------------------------------------------------------------------------------------- T5
    def test_T5_sameline(self):
        ident = getattr(self.__class__, "ident_abc", None)
        if ident is None:
            o = self.tmp / "ident_t5.json"
            subprocess.run([str(PY), str(self.base / "scripts/def_identity.py"), "--bundle", str(ROOT / "data/challenge"),
                            "--out", str(o)], cwd=self.base, env=clean_env(), capture_output=True, text=True, check=True)
            ident = json.loads(o.read_text())

        def with_stage2(ctx):
            d = json.loads(json.dumps(ident))
            per = {c: {"call": "X", "unanimous": False, "n_methods_agree": 0, "min_top1_score": 0.5, "top5": [],
                       "row_top1_methods": 0, "accepted": False} for c in "ABC"}
            per[ctx] = {"call": "HEK293_REF", "unanimous": True, "n_methods_agree": 8, "min_top1_score": 0.82, "top5": [],
                        "row_top1_methods": 8, "accepted": True}
            d["nonccle_ref"] = {"reference": "hek293", "sha256": P.REF_SHA256, "per_context": per,
                                "call": {c: ("HEK293" if c == ctx else None) for c in "ABC"}}
            return d

        def run_sl(tree, idj, tag, env, extra=()):
            ip = self.tmp / f"sl_{tag}.ident.json"
            ip.write_text(json.dumps(idj))
            o = self.tmp / f"sl_{tag}.npz"
            p = subprocess.run([str(PY), str(tree / "scripts/def_sameline.py"), "--bundle", str(ROOT / "data/challenge"),
                                "--identity", str(ip), "--out", str(o), *extra], cwd=tree, env=env, capture_output=True, text=True)
            self.assertEqual(p.returncode, 0, p.stdout[-1500:] + p.stderr[-1500:])
            z = np.load(o, allow_pickle=False)
            pv = json.loads(str(z["provenance"]))
            pv.pop("identity_source", None)          # the identity file's path differs per run by construction
            return {k: z[k] for k in z.files if k != "provenance"}, pv, p.stdout

        def same(a, b, skip=()):
            ka, kb = sorted(k for k in a if k not in skip), sorted(k for k in b if k not in skip)
            def eq(x, y):
                return np.array_equal(x, y, equal_nan=True) if np.issubdtype(np.asarray(x).dtype, np.floating) else np.array_equal(x, y)
            return ka == kb and all(eq(a[k], b[k]) for k in ka)

        rec = {}
        b0, pv0, _ = run_sl(self.base, ident, "base", clean_env())
        n0, pvn, _ = run_sl(self.tree, ident, "noblock", clean_env(DEF_DAY_ALLOW_NONCCLE_REF="hek293", DEF_DAY_NONCCLE_DECISION="x"))
        self.assertTrue(same(b0, n0, skip=("provenance",)))
        self.assertEqual(pv0, pvn)
        rec["no_block_identical"] = True
        for tag, env in (("unset", clean_env()), ("nodecision", clean_env(DEF_DAY_ALLOW_NONCCLE_REF="hek293")),
                         ("wrongknob", clean_env(DEF_DAY_ALLOW_NONCCLE_REF="rpe1", DEF_DAY_NONCCLE_DECISION="x"))):
            z, pv, out = run_sl(self.tree, with_stage2("B"), tag, env)
            self.assertNotIn("rows_B", z)
            self.assertIn("642", pv["contexts"]["B"]["mode"])
            self.assertTrue(same({k: v for k, v in b0.items() if k.endswith("_A")}, {k: v for k, v in z.items() if k.endswith("_A")}))
            rec[f"B_{tag}"] = pv["contexts"]["B"]["mode"]
        dec = "coordinator 2026-10-22: switch D to Orion HEK293T rows (642, unanimous 0.82)"
        z, pv, out = run_sl(self.tree, with_stage2("B"), "on", clean_env(DEF_DAY_ALLOW_NONCCLE_REF="hek293", DEF_DAY_NONCCLE_DECISION=dec))
        self.assertIn("rows_B", z)
        cb = pv["contexts"]["B"]
        self.assertEqual(cb["mode"], "same-line")
        self.assertEqual(cb["banks"], ["orion_hek293t"])
        self.assertEqual(cb["nonccle_decision"], dec)
        self.assertEqual(float(z["w_B"]), 2.0)
        self.assertTrue(same({k: v for k, v in b0.items() if k.endswith("_A")}, {k: v for k, v in z.items() if k.endswith("_A")}))
        rec["B_on"] = dict(covered_targets=cb["covered_targets"], w=float(z["w_B"]), call=cb["call"], banks=cb["banks"],
                           programme_bank=cb["programme_bank"])
        # HEK293 on A, whose CCLE call (JURKAT, registered, unanimous, 0.772) would switch on: vetoed without the knob
        z, pv, out = run_sl(self.tree, with_stage2("A"), "vetoA", clean_env())
        self.assertNotIn("rows_A", z)
        self.assertIn("642", pv["contexts"]["A"]["mode"])
        rec["A_veto"] = pv["contexts"]["A"]["mode"]
        # --force: untouched
        fb, fpv, _ = run_sl(self.base, ident, "fbase", clean_env(), ("--force", "A=JURKAT_HAEMATOPOIETIC_AND_LYMPHOID_TISSUE"))
        fp, fpp, _ = run_sl(self.tree, with_stage2("A"), "fpat", clean_env(DEF_DAY_ALLOW_NONCCLE_REF="hek293"),
                            ("--force", "A=JURKAT_HAEMATOPOIETIC_AND_LYMPHOID_TISSUE"))
        self.assertTrue(same(fb, fp, skip=("provenance",)))
        rec["force_identical"] = True
        RESULTS["T5"] = rec

    # ------------------------------------------------------------------------------------------------- T6
    def test_T6_day_lines(self):
        text = (self.tree / "scripts/def_day.sh").read_text()
        self.assertIn(P.DAY_LOG_B, text)
        self.assertIn(P.DAY_ENV_B, text)
        add_log, add_env = P.DAY_LOG_B[len(P.DAY_LOG_A):], P.DAY_ENV_B[len(P.DAY_ENV_A):]
        rec = {}
        dec = "coordinator: switch D (unanimous 0.82; see 642)"
        for tag, env in (("unset", {}), ("knob", {"DEF_DAY_ALLOW_NONCCLE_REF": "hek293"}),
                         ("both", {"DEF_DAY_ALLOW_NONCCLE_REF": "hek293", "DEF_DAY_NONCCLE_DECISION": dec})):
            w = self.tmp / f"day_{tag}"
            w.mkdir()
            script = f'set -euo pipefail\nW="{w}"\nlog() {{ printf "LOG %s\\n" "$*"; }}\n{add_log}{add_env}echo END\n'
            p = subprocess.run(["/bin/bash", "-c", script], env=clean_env(**env), capture_output=True, text=True)
            self.assertEqual(p.returncode, 0, p.stderr)
            envf = w / "def_day.final.env"
            logged = "642: the HEK293 contingency" in p.stdout
            body = envf.read_text() if envf.exists() else ""
            self.assertEqual(logged, tag != "unset")
            self.assertEqual("DEF_DAY_ALLOW_NONCCLE_REF=hek293" in body, tag != "unset")
            if tag == "both":
                back = subprocess.run(["/bin/bash", "-c", f'set -e; source "{envf}"; printf %s "$DEF_DAY_NONCCLE_DECISION"'],
                                      capture_output=True, text=True)
                self.assertEqual(back.stdout, dec)                                         # %q round-trips
            rec[tag] = dict(logged=logged, env=body)
        RESULTS["T6"] = rec


if __name__ == "__main__":
    unittest.main()
