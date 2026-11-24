#!/usr/bin/env python3
"""538: extras on a temp copy of the repo. Nothing here writes the live tree.

Run: uv run python experiments/538-cloud-bundle-extras/test_extras538.py
"""
from __future__ import annotations

import hashlib
import json
import os
import shutil
import subprocess
import sys
import tempfile
import textwrap
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
PY = ROOT / ".venv/bin/python"
APPLIER = ROOT / "experiments/538-cloud-bundle-extras/apply_october538.py"
PATCH538 = ROOT / "experiments/538-cloud-bundle-extras/patch_extras538.py"
PATCH530 = ROOT / "experiments/530-cloud-october-bundle/patch_bundle530.py"
PATCH528 = ROOT / "experiments/528-callprop-block/patch_defday528.py"
REF_META = ROOT / "submissions/453dr-october-path-abc-q1.75.h5ad.meta.json"
REL = [
    Path("scripts/def_day.sh"),
    Path("scripts/def_day_helpers.py"),
    Path("experiments/448-def-rehearsal-q175/relevel_iter.py"),
    Path("scripts/def_sameline.py"),
    Path("scripts/october_recipe.py"),
    Path("docs/october-decisions.json"),
    Path("experiments/524-october-recipe-dr/gate524.py"),
    Path("scripts/def_postbuild.py"),
    Path("scripts/def_pipeline.py"),
]
LIVE = [ROOT / p for p in REL]

BUILDER_STUB = """# grep tokens the day script and the resolver look for. Not the real builder.
VCC_SAMELINE
VCC_DE_MAG
VCC_DE_MDE_SIGN
programme:<tau>
file:<npz>
VCC_RC_CALLSET_CPM
VCC_RC_CALLSET_BLOCK
VCC_DE_PRIOR_GATE
def sharpen_q_by_tercile
"""
LODO_STUB = '''# stub, enough for the day's greps
ap.add_argument("--hi", type=float, default=3.0)
ap.add_argument("--hi-pct")
'''
UV_BODY = r'''
"""Stand-in for `uv`. Records argv and VCC_* then writes the artifacts def_day.sh waits on."""
import json, os, sys
from pathlib import Path

argv = sys.argv[1:]
log = Path(os.environ["STUB_LOG"])
rec = {"argv": argv, "vcc": {k: v for k, v in os.environ.items() if k.startswith("VCC_")}}
with log.open("a") as fh:
    fh.write(json.dumps(rec) + "\n")
if len(argv) < 3 or argv[0] != "run" or argv[1] != "python":
    sys.exit("stub uv: expected `run python ...`")
script, rest = argv[2], argv[3:]

def flag(name):
    return rest[rest.index(name) + 1]

if script.endswith("def_day_helpers.py") or script.endswith("relevel_iter.py"):
    os.execv(os.environ["REAL_PY"], [os.environ["REAL_PY"], script, *rest])

if script.endswith("def_pipeline.py"):
    name, stages, work = flag("--name"), flag("--stages"), flag("--work")
    w = Path(work)
    if "preflight" in stages:
        (w / "identity.json").write_text(json.dumps({"call": {"D": "X", "E": "Y", "F": "Z"}, "call_unanimous_across_methods": False}))
        (w / "preflight.json").write_text(json.dumps({"contexts": ["D", "E", "F"]}))
        (w / "census.json").write_text(json.dumps({"consensus": {"targets_with_0_donors": 2}}))
        print("--- [repin-panel-support] rc=0")
        sys.exit(0)
    if stages == "prep":
        print("dry run verified")
        sys.exit(0)
    if name.endswith("-dump"):
        dump = Path(os.environ["VCC_BULK_DUMP"])
        dump.write_bytes(b"dump")
        logs = w / "logs"
        logs.mkdir(parents=True, exist_ok=True)
        (logs / "build.log").write_text("[bulk-dump] wrote\n")
        sys.exit(0)
    if stages == "build":
        h5 = w / f"{name}.h5ad"
        cfg = {
            "bulk_sharpen_q": float(os.environ.get("VCC_BULK_SHARPEN_Q", "1.75")),
            "cd4_target_g": 0.3443655278649641,
            "rc_callset_cpm": int(os.environ["VCC_RC_CALLSET_CPM"]) if os.environ.get("VCC_RC_CALLSET_CPM", "").isdigit() else 2200,
        }
        if os.environ.get("VCC_BULK_CTX_SCALE"):
            cfg["bulk_ctx_scale"] = {k: float(v) for k, v in (p.split("=") for p in os.environ["VCC_BULK_CTX_SCALE"].split(","))}
        if os.environ.get("VCC_BULK_LOAD"):
            cfg["bulk_load"] = Path(os.environ["VCC_BULK_LOAD"].split(":")[0]).name
        if os.environ.get("VCC_DE_MAG"):
            cfg["de_mag"] = {"variant": os.environ["VCC_DE_MAG"]}
        sign = os.environ.get("VCC_DE_MDE_SIGN", "")
        if sign:
            p = sign.split(":")
            if p[0] == "clean4" and len(p) == 2:
                cfg["de_mde_sign"] = {"source": "clean4", "tau": float(p[1])}
        blk = os.environ.get("VCC_RC_CALLSET_BLOCK", "")
        if blk:
            p = blk.split(":")
            cfg["rc_callset_block"] = {
                "source": p[0],
                "sha256": p[2],
                "swap": float(p[3]) if len(p) >= 4 else 1.0,
                "steer": p[4].lower() if len(p) == 5 else "kept",
            }
        h5.write_bytes(b"h5")
        Path(str(h5) + ".meta.json").write_text(json.dumps({"config": cfg}))
        sys.exit(0)
    sys.exit(f"stub pipeline: unhandled name={name} stages={stages}")

if script.endswith("def_sameline.py"):
    import numpy as np
    out = flag("--out")
    np.savez(out, contexts=np.array([]))
    sys.exit(0)

if script.endswith("def_lodo_key.py"):
    out = Path(flag("--out"))
    out.write_bytes(b"lodo")
    hi = float(flag("--hi")) if "--hi" in rest else 3.0
    out.with_suffix(".json").write_text(json.dumps({"hi": hi, "hi_pct": None}))
    sys.exit(0)

if script.endswith("callblock528.py"):
    out = Path(flag("--out"))
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_bytes(b"callblock")
    sys.exit(0)

if script.endswith("def_pin_panel_artifacts.py"):
    sys.exit(0)

if script.endswith("ntclevel417.py"):
    h5 = rest[0]
    rho = float(os.environ["STUB_RHO"])
    per = {c: {"rho": rho, "scale": 1.0} for c in "DEF"}
    Path(h5.replace(".h5ad", ".ntclevel.json")).write_text(json.dumps({"per_context": per}))
    sys.exit(0)

if script.endswith("def_postbuild.py"):
    Path(flag("--out")).write_text(json.dumps({"checks": [{"check": "stub", "pass": True}]}))
    sys.exit(0)

if script.endswith("audit428.py"):
    print("ZERO block total (target excluded): 0 in 0 blocks")
    sys.exit(0)

sys.exit("stub uv: unhandled " + script)
'''


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def live_hashes() -> dict[str, str]:
    return {str(p): sha(p) for p in LIVE}


def copy_tree(dest: Path) -> None:
    for rel in REL:
        target = dest / rel
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy(ROOT / rel, target)


def run_applier(dest: Path, check: bool = False) -> subprocess.CompletedProcess:
    cmd = [str(PY), str(APPLIER), "--root", str(dest)]
    if check:
        cmd.append("--check")
    return subprocess.run(cmd, capture_output=True, text=True)


def rho_star() -> float:
    import math
    from scipy.stats import norm
    q = 0.31 / (0.31 + 1.78)
    return math.exp(2.0 * 0.01 * float(norm.ppf(q)))


def _sub(s: str, work: str, bundle: str) -> str:
    return s.replace(work, "$W").replace(bundle, "$B")


def _norm(records, work: str, bundle: str):
    out = []
    for rec in records:
        out.append({"argv": [_sub(a, work, bundle) for a in rec["argv"]],
                    "vcc": {k: _sub(v, work, bundle) for k, v in sorted(rec["vcc"].items())}})
    return out


def load_538():
    import importlib.util
    spec = importlib.util.spec_from_file_location("patch_extras538", PATCH538)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def constants_528() -> dict:
    src = PATCH528.read_text()
    start, end = src.index("VALIDATE ="), src.index("E = [")
    ns: dict = {}
    exec(src[start:end], ns)
    return ns


class Extras538(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls._live = live_hashes()
        cls.tmp = Path(tempfile.mkdtemp(prefix="extras538-"))
        cls.repo = cls.tmp / "repo"
        copy_tree(cls.repo)
        r = run_applier(cls.repo)
        if r.returncode != 0:
            raise RuntimeError(r.stdout + r.stderr)
        ref = cls.repo / "submissions/453dr-october-path-abc-q1.75.h5ad.meta.json"
        ref.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy(REF_META, ref)
        cls.mod = load_538()
        cls.star = f"{rho_star():.16g}"

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.tmp, ignore_errors=True)
        now = live_hashes()
        if now != cls._live:
            raise AssertionError("the test wrote the live tree")

    def test_check_clean_second_apply_refused_anchor_miss(self):
        dest = self.tmp / "check"
        copy_tree(dest)
        before = {str(p): sha(dest / p) for p in REL}
        r = run_applier(dest, check=True)
        self.assertEqual(r.returncode, 0, r.stderr + r.stdout)
        self.assertIn("patches apply cleanly (not written)", r.stdout)
        self.assertIn("528 did not apply", r.stdout)
        self.assertIn("before ", r.stdout)
        self.assertIn("after  ", r.stdout)
        self.assertEqual(before, {str(p): sha(dest / p) for p in REL})
        first = run_applier(dest)
        self.assertEqual(first.returncode, 0, first.stderr + first.stdout)
        digest = {str(p): sha(dest / p) for p in REL}
        self.assertNotEqual(digest, before)
        second = run_applier(dest)
        self.assertNotEqual(second.returncode, 0)
        self.assertIn("already patched", second.stderr + second.stdout)
        self.assertEqual(digest, {str(p): sha(dest / p) for p in REL})
        again = run_applier(dest, check=True)
        self.assertIn("already patched", again.stderr + again.stdout)
        self.assertEqual(digest, {str(p): sha(dest / p) for p in REL})
        broken = self.tmp / "broken"
        copy_tree(broken)
        sh = broken / "scripts/def_day.sh"
        sh.write_text(sh.read_text().replace(
            "REPO=~/Code/virtual-cell-challenge-2026; cd \"$REPO\" || exit 1\n",
            "REPO=/tmp/not-the-repo; cd \"$REPO\" || exit 1\n", 1))
        b2 = {str(p): sha(broken / p) for p in REL}
        miss = run_applier(broken, check=True)
        self.assertNotEqual(miss.returncode, 0, miss.stdout)
        self.assertIn("anchor count", miss.stderr + miss.stdout)
        self.assertEqual(b2, {str(p): sha(broken / p) for p in REL})

    def test_528_anchor_miss_leaves_bytes_and_snippets_match(self):
        dest = self.tmp / "only530"
        copy_tree(dest)
        r = subprocess.run([str(PY), str(PATCH530), "--root", str(dest)], capture_output=True, text=True)
        self.assertEqual(r.returncode, 0, r.stderr)
        day = dest / "scripts/def_day.sh"
        before = day.read_bytes()
        r = subprocess.run([str(PY), str(PATCH528), str(day)], capture_output=True, text=True)
        self.assertNotEqual(r.returncode, 0)
        self.assertIn("anchor found 0 times", r.stderr + r.stdout)
        self.assertEqual(day.read_bytes(), before)
        patched = (self.repo / "scripts/def_day.sh").read_text()
        ns = constants_528()
        for name in ("VALIDATE", "COMPUTE", "EXPORT", "FINALCHECK"):
            self.assertIn(ns[name], patched, name)
        self.assertIn(ns["VALIDATE"], patched)
        # re-anchored: VALIDATE sits immediately before the sharpening assignment, not after the MDE fi
        self.assertIn(ns["VALIDATE"] + 'SQ="$Q"\n', patched)
        self.assertNotIn(
            'die "the builder has no VCC_DE_MDE_SIGN (496p)"\nfi\n' + ns["VALIDATE"], patched)

    def test_bash_n_and_repo_derivation(self):
        sh = self.repo / "scripts/def_day.sh"
        bash = subprocess.run(["bash", "-n", str(sh)], capture_output=True, text=True)
        self.assertEqual(bash.returncode, 0, bash.stderr)
        text = sh.read_text()
        self.assertNotIn('REPO=~/Code/virtual-cell-challenge-2026; cd', text)
        self.assertIn('REPO="$(CDPATH= cd "$(dirname "$0")/.." && pwd)"', text)
        self.assertIn('re.fullmatch(r"[1-9][0-9]*(?:\\.[0-9]+)?", t)', text)
        wt = self.tmp / "worktree"
        (wt / "scripts").mkdir(parents=True)
        shutil.copy(sh, wt / "scripts/def_day.sh")
        traced = subprocess.run(["bash", "-x", "scripts/def_day.sh"], cwd=wt, capture_output=True, text=True)
        self.assertIn(f"REPO={wt.resolve()}", traced.stderr.replace("\n", " ") + "\n" + traced.stderr)
        # the assignment trace is `+ REPO=<path>`
        self.assertIn(f"+ REPO={wt.resolve()}", traced.stderr)
        loose = self.tmp / "loose"
        loose.mkdir()
        shutil.copy(sh, loose / "def_day.sh")
        bad = subprocess.run(["bash", str(loose / "def_day.sh")], capture_output=True, text=True)
        self.assertNotEqual(bad.returncode, 0)
        self.assertIn("no scripts/def_day.sh", bad.stderr)

    def test_unset_matches_530_plus_callset(self):
        base = self._day(kind="baseline", env={})
        full = self._day(kind="full", env={})
        self._assert_day_ok(base)
        self._assert_day_ok(full)
        self.assertEqual(base["records"], full["records"])
        self.assertEqual(base["final"], full["final"])
        blob = json.dumps(full["records"])
        self.assertNotIn("--census-max-zero", blob)
        self.assertNotIn('"--hi"', blob)
        self.assertNotIn("callblock528.py", blob)
        self.assertNotIn("VCC_RC_CALLSET_BLOCK", blob)
        self.assertNotIn("DEF_DAY_LODO_HI=", full["final"])
        self.assertNotIn("DEF_DAY_CENSUS_MAX_ZERO=", full["final"])
        self.assertNotIn("VCC_RC_CALLSET_BLOCK=", full["final"])
        october = {"DEF_DAY_LODO": "1", "DEF_DAY_DE_MAG": "dn_k250_zs"}
        b2 = self._day(kind="baseline", env=october)
        f2 = self._day(kind="full", env=october)
        self._assert_day_ok(b2)
        self._assert_day_ok(f2)
        self.assertEqual(b2["records"], f2["records"])
        self.assertEqual(b2["final"], f2["final"])
        lodo = [r for r in f2["records"] if any(a.endswith("def_lodo_key.py") for a in r["argv"])]
        self.assertEqual(len(lodo), 1)
        self.assertNotIn("--hi", lodo[0]["argv"])
        self.assertNotIn("--hi-pct", lodo[0]["argv"])
        self.assertTrue(any(a.endswith("pert_scale_lodo_def_lo0.5_hi3.npz") for a in lodo[0]["argv"]))
        # a 530 knob on, both 538 knobs off: the 4-arg riders spec must stay 530's string
        rho = {"DEF_DAY_RHO_SCALE": "1.02"}
        b3 = self._day(kind="baseline", env=rho)
        f3 = self._day(kind="full", env=rho)
        self._assert_day_ok(b3)
        self._assert_day_ok(f3)
        self.assertEqual(b3["records"], f3["records"])
        self.assertEqual(b3["final"], f3["final"])
        riders = [r for r in f3["records"] if "riders_match" in r["argv"]]
        self.assertEqual(riders[-1]["argv"][-1], "rc=-,stack=x1,rho=1.02,strength=-")
        print(f"unset identity: {len(full['records'])} uv calls match; LODO=1 {len(f2['records'])} match; "
              f"rho-only {len(f3['records'])} match", flush=True)

    def test_lodo_hi_set_and_refusals(self):
        got = self._day(kind="full", env={"DEF_DAY_LODO": "1", "DEF_DAY_LODO_HI": "6"})
        self._assert_day_ok(got)
        lodo = [r for r in got["records"] if any(a.endswith("def_lodo_key.py") for a in r["argv"])]
        argv = lodo[0]["argv"]
        self.assertEqual(argv[argv.index("--hi") + 1], "6")
        self.assertNotIn("--hi-pct", argv)
        self.assertTrue(any(a.endswith("pert_scale_lodo_def_lo0.5_hi6.npz") for a in argv))
        self.assertIn("lodo_lo0.5_hi6_def", got["log"])
        self.assertIn("DEF_DAY_LODO_HI=6", got["final"])
        riders = [r for r in got["records"] if "riders_match" in r["argv"]]
        self.assertIn("lodo_hi=6", riders[-1]["argv"][-1])
        self.assertNotIn("callset=", riders[-1]["argv"][-1])
        edge = self._day(kind="full", env={"DEF_DAY_LODO": "1", "DEF_DAY_LODO_HI": "20"})
        self._assert_day_ok(edge)
        self.assertIn("lodo_lo0.5_hi20_def", edge["log"])
        refusals = [
            ({"DEF_DAY_LODO_HI": "6"}, "needs DEF_DAY_LODO=1"),
            ({"DEF_DAY_LODO": "1", "DEF_DAY_LODO_HI": "6", "DEF_DAY_LODO_HI_PCT": "96"}, "two cap forms"),
            ({"DEF_DAY_LODO": "1", "DEF_DAY_LODO_HI": "1"}, "1 < h <= 20"),
            ({"DEF_DAY_LODO": "1", "DEF_DAY_LODO_HI": "21"}, "1 < h <= 20"),
            ({"DEF_DAY_LODO": "1", "DEF_DAY_LODO_HI": "1e1"}, "1 < h <= 20"),
            ({"DEF_DAY_LODO": "1", "DEF_DAY_LODO_HI": "06"}, "1 < h <= 20"),
            ({"DEF_DAY_LODO": "1", "DEF_DAY_LODO_HI": "0.5"}, "1 < h <= 20"),
            ({"DEF_DAY_LODO": "1", "DEF_DAY_LODO_HI": "6", "DEF_DAY_LODO_HI_PCT": "96", "DEF_DAY_LODO_STRENGTH": "0.25"},
             "two cap forms"),
            ({"DEF_DAY_LODO": "1", "DEF_DAY_LODO_HI": "6", "DEF_DAY_LODO_STRENGTH": "0.25"}, "percentile cap"),
        ]
        for env, needle in refusals:
            bad = self._day(kind="full", env=env)
            self.assertNotEqual(bad["rc"], 0, env)
            self.assertIn(needle, bad["log"] + bad["err"], env)

    def test_callset_set_and_refusals(self):
        for value, swap, steer in (("prop", "1.0", "kept"), ("prop:0.5:key", "0.5", "key")):
            got = self._day(kind="full", env={"DEF_DAY_CALLSET_BLOCK": value})
            self._assert_day_ok(got)
            calls = [r for r in got["records"] if any(a.endswith("callblock528.py") for a in r["argv"])]
            self.assertEqual(len(calls), 1, value)
            self.assertIn("--contexts", calls[0]["argv"])
            builds = [r for r in got["records"] if _name(r) and _name(r).endswith("-lvl") or (
                "--stages" in r["argv"] and r["argv"][r["argv"].index("--stages") + 1] == "build" and not _name(r).endswith("-dump"))]
            level = [r for r in got["records"] if "--stages" in r["argv"] and r["argv"][r["argv"].index("--stages") + 1] == "build"
                     and not _name(r).endswith("-dump")]
            self.assertTrue(level, value)
            self.assertTrue(level[0]["vcc"]["VCC_RC_CALLSET_BLOCK"].startswith("prop:"), value)
            self.assertIn("VCC_RC_CALLSET_BLOCK=", got["final"])
            riders = [r for r in got["records"] if "riders_match" in r["argv"]]
            self.assertIn(f"callset={value}", riders[-1]["argv"][-1])
            self.assertIn(f":{swap}:" if ":" in value else "prop:", level[0]["vcc"]["VCC_RC_CALLSET_BLOCK"])
            if steer != "kept" or ":" in value and value.count(":") == 2:
                self.assertTrue(level[0]["vcc"]["VCC_RC_CALLSET_BLOCK"].endswith(":" + steer) or steer == "kept")
        refusals = [
            ({"DEF_DAY_CALLSET_BLOCK": "banana"}, "must be prop"),
            ({"DEF_DAY_CALLSET_BLOCK": "prop:2"}, "0 < swap"),
            ({"DEF_DAY_CALLSET_BLOCK": "prop:0"}, "0 < swap"),
            ({"DEF_DAY_CALLSET_BLOCK": "prop:0.5:nope"}, "kept|key|cpm"),
        ]
        for env, needle in refusals:
            bad = self._day(kind="full", env=env)
            self.assertNotEqual(bad["rc"], 0, env)
            self.assertIn(needle, bad["log"] + bad["err"], env)
            self.assertFalse(any("callblock528.py" in a for r in bad["records"] for a in r["argv"]))

    def test_census_override_and_refusals(self):
        got = self._day(kind="full", env={"DEF_DAY_CENSUS_MAX_ZERO": "8"})
        self._assert_day_ok(got)
        setup = [r for r in got["records"] if _name(r) == "day-setup"]
        argv = setup[0]["argv"]
        self.assertEqual(argv[argv.index("--census-max-zero") + 1], "8")
        self.assertIn("census override: max zero-carrier 8 (count 2)", got["log"])
        self.assertIn("DEF_DAY_CENSUS_MAX_ZERO=8", got["final"])
        # the same work dir is now not fresh
        again = self._rerun(got, {"DEF_DAY_CENSUS_MAX_ZERO": "8"})
        self.assertNotEqual(again["rc"], 0)
        self.assertIn("fresh work dir", again["log"] + again["err"])
        for bad_value in ("5", "31", "6.0", "06", "0", "abc"):
            bad = self._day(kind="full", env={"DEF_DAY_CENSUS_MAX_ZERO": bad_value})
            self.assertNotEqual(bad["rc"], 0, bad_value)
            self.assertIn("DEF_DAY_CENSUS_MAX_ZERO must be an integer", bad["log"] + bad["err"])
            self.assertEqual(bad["records"], [], bad_value)
        edge = self._day(kind="full", env={"DEF_DAY_CENSUS_MAX_ZERO": "30"})
        self._assert_day_ok(edge)
        self.assertIn("census override: max zero-carrier 30 (count 2)", edge["log"])

    def test_pipeline_forwards_max_zero_carrier(self):
        bare = self.tmp / "pipe-bare"
        patched = self.tmp / "pipe-patched"
        for dest, src in ((bare, ROOT / "scripts/def_pipeline.py"), (patched, self.repo / "scripts/def_pipeline.py")):
            (dest / "scripts").mkdir(parents=True)
            shutil.copy(src, dest / "scripts/def_pipeline.py")
        bundle = self.tmp / "abc-manifest"
        bundle.mkdir()
        (bundle / "manifest.json").write_text(json.dumps({"contexts": ["A", "B", "C"], "panel_id": "t"}))
        work = self.tmp / "pipe-work"
        frozen = self.tmp / "frozen.json"

        def run(script, extra):
            return subprocess.run(
                [str(PY), str(script), "--bundle", str(bundle), "--work", str(work), "--name", "n",
                 "--dry-run", "--stages", "census", "--frozen-inputs", str(frozen), *extra],
                capture_output=True, text=True)

        old = run(bare / "scripts/def_pipeline.py", [])
        new = run(patched / "scripts/def_pipeline.py", [])
        self.assertEqual(old.returncode, 0, old.stderr)
        self.assertEqual(new.returncode, 0, new.stderr)
        self.assertEqual(old.stdout, new.stdout)
        self.assertNotIn("--max-zero-carrier", new.stdout)
        raised = run(patched / "scripts/def_pipeline.py", ["--census-max-zero", "8"])
        self.assertEqual(raised.returncode, 0, raised.stderr)
        self.assertIn("--max-zero-carrier 8", raised.stdout)
        self.assertNotIn("--max-zero-carrier", old.stdout)
        for bad in ("5", "31", "06"):
            refused = run(patched / "scripts/def_pipeline.py", ["--census-max-zero", bad])
            self.assertNotEqual(refused.returncode, 0, bad)
            self.assertIn("--census-max-zero must be an integer", refused.stderr + refused.stdout)

    def test_postbuild_require_config_diff(self):
        world = self.tmp / "post-world"
        meta, inc, bundle, builder = _post_world(world)
        unpatched = ROOT / "scripts/def_postbuild.py"
        patched = self.repo / "scripts/def_postbuild.py"
        script = world / "scripts/def_postbuild.py"
        shutil.copy(unpatched, script)
        base_cmd = [str(PY), str(script), "--meta", str(meta), "--incumbent", str(inc),
                    "--bundle", str(bundle), "--builder", str(builder), "--cd4-root", str(world / "cd4"),
                    "--allow-arm-drift", "cd4_deflate_basis"]
        old = subprocess.run(base_cmd, capture_output=True, text=True)
        self.assertEqual(old.returncode, 0, old.stdout + old.stderr)
        shutil.copy(patched, script)
        new = subprocess.run(base_cmd, capture_output=True, text=True)
        self.assertEqual(new.returncode, 0, new.stdout + new.stderr)
        self.assertEqual(old.stdout, new.stdout)
        self.assertEqual(old.stderr, new.stderr)

        def with_keys(meta_path, **changed):
            doc = json.loads(meta.read_text())
            doc["config"].update(changed)
            meta_path.write_text(json.dumps(doc))
            return meta_path

        ok_meta = with_keys(world / "ok.meta.json", k1="moved", k2="moved")
        inc_doc = json.loads(inc.read_text())
        inc_doc["config"].update(k1="old", k2="old")
        inc2 = world / "inc2.json"
        inc2.write_text(json.dumps(inc_doc))
        common = ["--incumbent", str(inc2), "--bundle", str(bundle), "--builder", str(builder),
                  "--cd4-root", str(world / "cd4"), "--allow-arm-drift", "cd4_deflate_basis",
                  "--allow-config-diff", "k1,k2"]
        passed = subprocess.run([str(PY), str(script), "--meta", str(ok_meta), *common,
                                 "--require-config-diff", "k1,k2"], capture_output=True, text=True)
        self.assertEqual(passed.returncode, 0, passed.stdout + passed.stderr)
        self.assertIn("every --require-config-diff key moved", passed.stdout)
        missing_meta = with_keys(world / "miss.meta.json", k1="moved", k2="old")
        missed = subprocess.run([str(PY), str(script), "--meta", str(missing_meta), *common,
                                 "--require-config-diff", "k1,k2"], capture_output=True, text=True)
        self.assertNotEqual(missed.returncode, 0)
        self.assertIn("missing=['k2']", missed.stdout)
        extra_meta = with_keys(world / "extra.meta.json", k1="moved", k2="moved", k3="surprise")
        extra = subprocess.run([str(PY), str(script), "--meta", str(extra_meta), *common,
                                "--require-config-diff", "k1,k2"], capture_output=True, text=True)
        self.assertNotEqual(extra.returncode, 0)
        self.assertIn("k3", extra.stdout)

    def test_recipe_pending_exports_nothing_new(self):
        r = self._recipe("command", "--allow-pending", "--bundle", "B", "--work", "W", "--name", "N")
        self.assertEqual(r.returncode, 0, r.stderr + r.stdout)
        lines = [ln for ln in r.stdout.splitlines() if ln and not ln.startswith("#")]
        # The live command grows with every `decide` (sharpen_pt was decided 2026-10-02), so pin invariants, not one
        # string (650 L3, 2026-10-02).
        self.assertEqual(len(lines), 1, lines)
        self.assertTrue(lines[0].endswith("bash scripts/def_day.sh B W N 1.75 0.5 noprog 4"), lines[0])
        for part in ("DEF_DAY_LODO=1", "DEF_DAY_DE_MAG=dn_k250_zs"):
            self.assertIn(part, lines[0])
        # "DEF_DAY_LODO_HI=" exactly: DEF_DAY_LODO_HI_PCT (decided 2026-10-03, 495) is the other cap form and may be exported.
        self.assertNotIn("DEF_DAY_LODO_HI=", lines[0])
        self.assertNotIn("DEF_DAY_CALLSET_BLOCK", lines[0])
        d = json.loads((self.repo / "docs/october-decisions.json").read_text())
        sw = d["switches"]
        self.assertEqual(sw["lodo_hi"]["status"], "PENDING")
        self.assertIsNone(sw["lodo_hi"]["value"])
        self.assertEqual(sw["callset_block"]["status"], "PENDING")
        self.assertIsNone(sw["callset_block"]["value"])
        self.assertEqual(sw["cd4_stack"]["status"], "DECIDED")
        self.assertEqual(sw["cd4_stack"]["value"], "x1")
        self.assertEqual(sw["lodo"]["value"], "1")
        self.assertEqual(sw["de_mag"]["value"], "dn_k250_zs")
        names = list(sw)
        self.assertLess(names.index("lodo_hi"), names.index("cd4_thresholds"))
        self.assertLess(names.index("callset_block"), names.index("cd4_thresholds"))

    def test_recipe_decided_and_refusals(self):
        ok = self._recipe("command", "--allow-pending", "--json", "--bundle", "B", "--work", "W", "--name", "N",
                          decisions=self._decide(lodo_hi="6"))
        self.assertEqual(ok.returncode, 0, ok.stderr + ok.stdout)
        env = json.loads(ok.stdout)["env"]
        self.assertEqual(env["DEF_DAY_LODO_HI"], "6")
        self.assertNotIn("DEF_DAY_CALLSET_BLOCK", env)
        cs = self._recipe("command", "--allow-pending", "--json", "--bundle", "B", "--work", "W", "--name", "N",
                          decisions=self._decide(callset="prop:0.5:key"))
        self.assertEqual(cs.returncode, 0, cs.stderr + cs.stdout)
        self.assertEqual(json.loads(cs.stdout)["env"]["DEF_DAY_CALLSET_BLOCK"], "prop:0.5:key")
        for kwargs, needle in (
            ({"lodo_hi": "1"}, "1 < h <= 20"),
            ({"lodo_hi": "21"}, "1 < h <= 20"),
            ({"lodo_hi": "1e1"}, "1 < h <= 20"),
            ({"lodo_hi": "6", "with_lodo": False}, "lodo = 1"),
            ({"lodo_hi": "6", "with_pct": True}, "two cap forms"),
            ({"callset": "prop:2"}, "callset_block"),
            ({"callset": "nope"}, "callset_block"),
            ({"callset": "prop:0.5:sideways"}, "callset_block"),
        ):
            proc = self._recipe("command", "--allow-pending", "--bundle", "B", "--work", "W", "--name", "N",
                                decisions=self._decide(**kwargs))
            self.assertEqual(proc.returncode, 3, (kwargs, proc.returncode, proc.stderr, proc.stdout))
            self.assertIn(needle, proc.stderr)

    def test_riders_match_extension(self):
        helpers = str(self.repo / "scripts/def_day_helpers.py")
        star = self.tmp / "metas"
        star.mkdir(exist_ok=True)

        def meta(name, cfg):
            p = star / name
            p.write_text(json.dumps({"config": cfg}))
            return str(p)

        def run(*args):
            return subprocess.run([str(PY), helpers, "riders_match", *args], capture_output=True, text=True)

        base = {"bulk_sharpen_q": 1.75, "cd4_target_g": 0.3443655278649641, "rc_callset_cpm": 2200,
                "bulk_load": "dump_h1x2_cd4x1.npz"}
        sha64 = "ab" * 32
        block = {"source": "prop", "sha256": sha64, "swap": 0.5, "steer": "key"}
        accept = [
            (meta("four.json", dict(base)), "-", "-", "rc=-,stack=x1,rho=-,strength=-"),
            (meta("absent.json", dict(base)), "-", "-", "rc=-,stack=x1,rho=-,strength=-,lodo_hi=6"),
            (meta("recorded.json", {**base, "lodo_hi": 6}), "-", "-", "rc=-,stack=x1,rho=-,strength=-,lodo_hi=6"),
            (meta("cs.json", {**base, "rc_callset_block": block}), "-", "-",
             "rc=-,stack=x1,rho=-,strength=-,callset=prop:0.5:key"),
            (meta("plain.json", {**base, "rc_callset_block": {"source": "prop", "sha256": sha64, "swap": 1.0, "steer": "kept"}}),
             "-", "-", "rc=-,stack=x1,rho=-,strength=-,callset=prop"),
        ]
        for args in accept:
            got = run(*args)
            self.assertEqual(got.returncode, 0, (args[-1], got.stderr))
        reject = [
            (meta("badhi.json", {**base, "lodo_hi": 4}), "-", "-", "rc=-,stack=x1,rho=-,strength=-,lodo_hi=6"),
            (meta("noblock.json", dict(base)), "-", "-", "rc=-,stack=x1,rho=-,strength=-,callset=prop"),
            (meta("badswap.json", {**base, "rc_callset_block": {**block, "swap": 1.0}}), "-", "-",
             "rc=-,stack=x1,rho=-,strength=-,callset=prop:0.5:key"),
            (meta("badsteer.json", {**base, "rc_callset_block": {**block, "steer": "kept"}}), "-", "-",
             "rc=-,stack=x1,rho=-,strength=-,callset=prop:0.5:key"),
        ]
        for args in reject:
            got = run(*args)
            self.assertNotEqual(got.returncode, 0, args[-1])
        bad = run(meta("spec.json", dict(base)), "-", "-", "rc=-,stack=x1,rho=-,strength=-,nope=1")
        self.assertEqual(bad.returncode, 2, bad.stderr)
        short = run(meta("short.json", dict(base)), "-", "-", "rc=-")
        self.assertEqual(short.returncode, 2)

    def test_gate_rows(self):
        gate = str(self.repo / "experiments/524-october-recipe-dr/gate524.py")
        sha64 = "ab" * 32

        def run(env, *, meta_edit=None, log_extra="", sidecar=None):
            work = Path(tempfile.mkdtemp(prefix="gate538-", dir=self.tmp))
            run_dir = work / "run"
            run_dir.mkdir()
            (work / "recipe524.json").write_text(json.dumps({
                "env": env, "args": ["B", "W", "N", "1.75", "auto", "pt", "4"],
                "pending_treated_off": [], "problems": [],
            }))
            (run_dir / "def_day.log").write_text(
                "CONVERGED at D=1,E=1,F=1 after 0 levelled build(s)\n"
                "final: lvl0 factors D=1,E=1,F=1 -> out\n" + log_extra)
            (run_dir / "lvl0.ntclevel.json").write_text(json.dumps({"per_context": {c: {"rho": 0.5} for c in "DEF"}}))
            (run_dir / "postbuild.json").write_text(json.dumps({"checks": [{"check": "a", "pass": True}]}))
            (run_dir / "audit428.txt").write_text("ZERO block total (target excluded): 0 in 0 blocks\n")
            if sidecar is not None:
                (run_dir / sidecar[0]).write_text(json.dumps(sidecar[1]))
            meta = json.loads(REF_META.read_text())
            if meta_edit:
                meta_edit(meta["config"])
            stem = work / "pkg"
            Path(str(stem) + ".h5ad.meta.json").write_text(json.dumps(meta))
            blob = b"vcc-synthetic"
            Path(str(stem) + ".vcc").write_bytes(blob)
            digest = hashlib.sha256(blob).hexdigest()
            chain = work / "chain.log"
            chain.write_text(f"novelty rc=0\npackage pkg.vcc sha256 {digest}\n")
            return subprocess.run([str(PY), gate, "--write", "--work", str(work), "--stem", str(stem),
                                   "--chainlog", str(chain)], capture_output=True, text=True)

        base = run({})
        self.assertEqual(base.returncode, 0, base.stdout + base.stderr)
        hi = run({"DEF_DAY_LODO_HI": "6"}, log_extra="lodo_lo0.5_hi6_def\n",
                 sidecar=("pert_scale_lodo_def_lo0.5_hi6.json", {"hi": 6}))
        self.assertEqual(hi.returncode, 0, hi.stdout + hi.stderr)
        self.assertIn("LODO fixed cap", hi.stdout)
        hi_bad = run({"DEF_DAY_LODO_HI": "6"}, log_extra="lodo_lo0.5_hi6_def\n",
                     sidecar=("pert_scale_lodo_def_lo0.5_hi6.json", {"hi": 3}))
        self.assertIn("ALL GATES FAIL", hi_bad.stdout)
        block = {"source": "prop", "sha256": sha64, "swap": 1.0, "steer": "kept"}
        cs = run({"DEF_DAY_CALLSET_BLOCK": "prop"}, meta_edit=lambda c: c.update(rc_callset_block=block))
        self.assertEqual(cs.returncode, 0, cs.stdout + cs.stderr)
        self.assertIn("call-propensity block", cs.stdout)
        cs_bad = run({"DEF_DAY_CALLSET_BLOCK": "prop:0.5:key"},
                     meta_edit=lambda c: c.update(rc_callset_block=block))
        self.assertIn("ALL GATES FAIL", cs_bad.stdout)

    # --- day harness -----------------------------------------------------------------

    def _baseline_day_text(self, repo: Path) -> str:
        """530 applied, then 528's snippets re-anchored, with the Mac REPO rewritten to this tree.
        That is the 530+528 behaviour the unset comparison is against: 528's own patch does not apply."""
        r = subprocess.run([str(PY), str(PATCH530), "--root", str(repo)], capture_output=True, text=True)
        if r.returncode != 0:
            raise AssertionError(r.stderr)
        day = repo / "scripts/def_day.sh"
        text = self.mod._apply(day.read_text(), self.mod.callset_edits(), "def_day.sh")
        old = "REPO=~/Code/virtual-cell-challenge-2026"
        if old not in text:
            raise AssertionError("baseline lost the hardcoded REPO")
        return text.replace(old, f"REPO={repo}", 1)

    def _day(self, kind: str, env: dict) -> dict:
        dest = Path(tempfile.mkdtemp(prefix="day538-", dir=self.tmp))
        repo = dest / "repo"
        copy_tree(repo)
        if kind == "baseline":
            text = self._baseline_day_text(repo)
            # helpers/relevel stay at 530 (copy_tree + patch 530 already rewrote them inside _baseline)
        elif kind == "full":
            for rel in REL:
                target = repo / rel
                target.parent.mkdir(parents=True, exist_ok=True)
                shutil.copy(self.repo / rel, target)
            text = (repo / "scripts/def_day.sh").read_text()
        else:
            raise AssertionError(kind)
        (repo / "scripts/def_day.sh").write_text(text)
        self._stub_tree(repo)
        return self._run_bash(dest, repo, env)

    def _stub_tree(self, repo: Path) -> None:
        files = {
            Path("scripts/make_submission_pm_cd4.py"): BUILDER_STUB,
            Path("scripts/def_pipeline.py"): "# stub\n",
            Path("scripts/def_lodo_key.py"): LODO_STUB,
            Path("scripts/def_pin_panel_artifacts.py"): "# stub\n",
            Path("scripts/def_postbuild.py"): "# stub\n",
            Path("experiments/522-october-bulk-key/def_strength_key522.py"): "# stub\n",
            Path("experiments/417-bulk-sharpen/ntclevel417.py"): "# stub\n",
            Path("experiments/428-zero-gene-audit/audit428.py"): "# stub\n",
            Path("experiments/528-callprop-block/callblock528.py"): "# stub\n",
            Path("src/vcc2026/rc_emit_k.py"): "block_steer=\n",
        }
        for rel, body in files.items():
            p = repo / rel
            p.parent.mkdir(parents=True, exist_ok=True)
            p.write_text(body)
        g = ["git", "-C", str(repo), "-c", "user.name=t", "-c", "user.email=t@t", "-c", "commit.gpgsign=false"]
        subprocess.run(["git", "init", "-q", str(repo)], check=True)
        subprocess.run(g + ["add", "-A"], check=True)
        subprocess.run(g + ["commit", "-q", "-m", "stub"], check=True)

    def _run_bash(self, dest: Path, repo: Path, env: dict, work: Path | None = None) -> dict:
        bindir = dest / "bin"
        if not bindir.exists():
            bindir.mkdir()
            uv = bindir / "uv"
            uv.write_text("#!" + str(PY) + "\n" + UV_BODY)
            uv.chmod(0o755)
        bundle = dest / "bundle"
        if not bundle.exists():
            bundle.mkdir()
            (bundle / "manifest.json").write_text(json.dumps({"contexts": ["D", "E", "F"]}))
        if work is None:
            work = dest / "work"
        log = dest / "uv.jsonl"
        if log.exists():
            log.unlink()
        run_env = os.environ.copy()
        for key in list(run_env):
            if key.startswith("VCC_") or key.startswith("DEF_DAY_"):
                del run_env[key]
        run_env.update(env)
        run_env["PATH"] = str(bindir) + os.pathsep + run_env.get("PATH", "")
        run_env["STUB_LOG"] = str(log)
        run_env["STUB_RHO"] = self.star
        run_env["REAL_PY"] = str(PY)
        proc = subprocess.run(
            ["bash", "scripts/def_day.sh", str(bundle), str(work), "day", "1.75", "auto", "pt:noprog", "4"],
            cwd=repo, env=run_env, capture_output=True, text=True)
        raw = []
        if log.is_file():
            raw = [json.loads(line) for line in log.read_text().splitlines() if line.strip()]
        day_log = (work / "def_day.log").read_text() if (work / "def_day.log").is_file() else proc.stdout
        final = (work / "def_day.final.env").read_text() if (work / "def_day.final.env").is_file() else ""
        return {"rc": proc.returncode, "raw": raw, "records": _norm(raw, str(work), str(bundle)),
                "final": _sub(final, str(work), str(bundle)), "log": day_log, "err": proc.stderr,
                "dest": dest, "repo": repo, "work": work, "env": env}

    def _rerun(self, previous: dict, env: dict) -> dict:
        return self._run_bash(previous["dest"], previous["repo"], env, work=previous["work"])

    def _assert_day_ok(self, got):
        if got["rc"] != 0:
            self.fail(f"rc {got['rc']}\n{got['log'][-2500:]}\n{got['err'][-1500:]}")

    def _decide(self, lodo_hi=None, callset=None, with_lodo=True, with_pct=False) -> Path:
        d = json.loads((self.repo / "docs/october-decisions.json").read_text())
        sw = d["switches"]
        if not with_lodo:
            sw["lodo"].update(status="DECIDED", value="0")
        if with_pct:
            sw["lodo_hi_pct"].update(status="DECIDED", value="96.64697522769474")
        else:
            # lodo_hi_pct was DECIDED live on 2026-10-03 (the 495 read); these fixtures exercise lodo_hi alone, so the copy
            # resets the percentile form to its pre-decide state.
            sw["lodo_hi_pct"].update(status="PENDING", value=None)
        if lodo_hi is not None:
            sw["lodo_hi"].update(status="DECIDED", value=lodo_hi)
        if callset is not None:
            sw["callset_block"].update(status="DECIDED", value=callset)
        path = self.tmp / f"dec-{hashlib.sha256(json.dumps(d).encode()).hexdigest()[:10]}.json"
        path.write_text(json.dumps(d))
        return path

    def _recipe(self, *args, decisions=None):
        env = os.environ.copy()
        env["OCTOBER_DECISIONS"] = str(decisions or (self.repo / "docs/october-decisions.json"))
        env["OCTOBER_DEF_DAY"] = str(self.repo / "scripts/def_day.sh")
        env["OCTOBER_BUILDER"] = str(ROOT / "scripts/make_submission_pm_cd4.py")
        env["OCTOBER_SAMELINE"] = str(self.repo / "scripts/def_sameline.py")
        env["OCTOBER_EMITTER"] = str(ROOT / "src/vcc2026/rc_emit_k.py")
        env["OCTOBER_CALLBLOCK"] = str(ROOT / "experiments/528-callprop-block/callblock528.py")
        return subprocess.run([str(PY), str(self.repo / "scripts/october_recipe.py"), *args],
                              capture_output=True, text=True, env=env)


def _name(rec) -> str:
    argv = rec["argv"]
    if "--name" in argv:
        return argv[argv.index("--name") + 1]
    return ""


def _post_world(world: Path):
    """A synthetic new-panel world that passes every def_postbuild check except the diffs under test.
    No challenge bundle, no donor tables from the repo: every input is written here."""
    import numpy as np
    world.mkdir(parents=True, exist_ok=True)
    genes = [f"T{i:02d}" for i in range(100)] + ["GEXTRA"]
    panel = genes[:-1]
    (world / "scripts").mkdir()
    for donor in ("k562_gwps", "orion_hct116", "orion_hek293t"):
        d = world / "data/donors" / donor
        d.mkdir(parents=True)
        (d / "genes.txt").write_text("GEXTRA\n")
        (d / "perts.txt").write_text("\n".join(panel) + "\n")
    pm = world / "data/external/perturbmatch"
    pm.mkdir(parents=True)
    (pm / "perts.txt").write_text("\n".join(panel) + "\n")
    cd4 = world / "cd4"
    for name in ("cd4_rest", "cd4_stim8hr", "cd4_stim48hr"):
        b = cd4 / name
        b.mkdir(parents=True)
        (b / "perts.txt").write_text("\n".join(panel) + "\n")
    bundle = world / "bundle"
    bundle.mkdir()
    (bundle / "manifest.json").write_text(json.dumps({
        "contexts": ["D", "E", "F"], "n_constructs": 100, "cells_per_pert": 400, "panel_id": "synth538"}))
    (bundle / "pert_counts.csv").write_text("target_gene\n" + "\n".join(panel) + "\n")
    (bundle / "gene_names.csv").write_text("gene_name\n" + "\n".join(genes) + "\n")
    op = world / "op.bin"
    op.write_bytes(b"operator-bytes-538")
    op_sha = hashlib.sha256(op.read_bytes()).hexdigest()
    meta_sha = "b" * 64
    tindex = {g: i for i, g in enumerate(genes)}
    expect_cols = np.asarray(sorted(tindex[g] for g in {"GEXTRA"} - set(panel)), dtype=np.int64)
    cols_sha = hashlib.sha256(np.ascontiguousarray(expect_cols, dtype="<i8").tobytes()).hexdigest()
    axis_sha = hashlib.sha256("\n".join(str(genes[i]) for i in expect_cols).encode()).hexdigest()
    builder = world / "scripts/builder.py"
    builder.write_text(textwrap.dedent(f"""\
        WHITENER_REGISTRY = {{
            "union_pruned_def": ("{op_sha}", "{meta_sha}"),
        }}
        LEGACY_OPERATOR_AXIS_SHA256 = {{}}
    """))
    cfg = {
        "contexts": ["D", "E", "F"], "targets": 100, "whitener_id": "union_pruned_def",
        "alpha": 0.982, "cd4_target_f": 0.34, "h1_union_g": 0.28,
        "h1_union_sources": ["gse295214_h1", "h1_2025"], "emitter": "rc_callset",
        "cd4_deflate_k": 5, "cd4_sizing": "carrying", "cd4_target_g": 0.3443655278649641,
    }
    operator = {
        "sha256": op_sha, "metadata_sha256": meta_sha, "path": str(op),
        "metadata_path": str(world / "whitener.json"), "columns": int(len(expect_cols)),
        "columns_sha256": cols_sha, "axis_sha256": axis_sha,
        "axis_sha256_checked_against": "metadata:whitener.json",
        "all_panel_coordinates_excluded": True,
        "metadata": {"n_cols": int(len(expect_cols)), "rows_kept_per_donor": {"k562_gwps": 1}, "basis": [1]},
    }
    effects = {
        "perturbmatch": {"both_supported_modified_count": 100},
        "cd4": {"fully_supported_targets": 100, "residual_frame": {"sizing_rule": "carrying"}},
        "whitening": {"after_all_median_norm": 1.0},
    }
    output = {"nnz": 10, "n_rows": 100 * 400 * 3,
              "contexts": {c: {"emitted_rows": 40000, "stored_entries": 1} for c in "DEF"}}
    frozen = {"panel_support": {"cd4_complete_targets": 100, "perturbmatch_replacement_targets": 100}}
    meta = {"config": cfg, "operator": operator, "effects": effects, "output": output, "frozen_bundle": frozen}
    meta_path = world / "meta.json"
    meta_path.write_text(json.dumps(meta))
    inc_path = world / "inc.json"
    inc_path.write_text(json.dumps({"config": dict(cfg), "effects": {"whitening": {"after_all_median_norm": 1.0}}}))
    return meta_path, inc_path, bundle, builder


if __name__ == "__main__":
    unittest.main(verbosity=2)
