"""530: the October pass-through bundle, tested on copies. Nothing here writes the live tree.

Run: uv run python experiments/530-cloud-october-bundle/test_bundle530.py
"""
from __future__ import annotations

import hashlib
import itertools
import json
import os
import shutil
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]  # was /workspace on the cloud VM
PY = ROOT / ".venv/bin/python"
PATCH = ROOT / "experiments/530-cloud-october-bundle/patch_bundle530.py"
REL = {
    "def_day": Path("scripts/def_day.sh"),
    "helpers": Path("scripts/def_day_helpers.py"),
    "relevel": Path("experiments/448-def-rehearsal-q175/relevel_iter.py"),
    "sameline": Path("scripts/def_sameline.py"),
    "recipe": Path("scripts/october_recipe.py"),
    "decisions": Path("docs/october-decisions.json"),
    "gate": Path("experiments/524-october-recipe-dr/gate524.py"),
}
LIVE = {k: ROOT / p for k, p in REL.items()}
REF_META = ROOT / "submissions/453dr-october-path-abc-q1.75.h5ad.meta.json"
X1_G = 0.3443655278649641
X15_G = "0.5165482917974462"
REFUSE526 = """# A reboot skips this script's EXIT trap, so the builder can still be the LODO pin. A fresh work
# dir would record that pin as BSHA0 and the trap would then "restore" it (526). Refuse.
_HEADB="$(git show HEAD:scripts/make_submission_pm_cd4.py | shasum -a 256 | cut -d' ' -f1)"
if [ "$BSHA0" != "$_HEADB" ]; then
  die "builder on disk ${BSHA0:0:16} is not HEAD ${_HEADB:0:16}. Restore before a fresh work dir: uv run python scripts/def_pin_panel_artifacts.py --restore --record <previous-work>/pin-artifacts.json"
fi
"""
BUILDER_STUB = """# grep tokens the day script and the resolver look for. Not the real builder.
VCC_SAMELINE
VCC_DE_MAG
VCC_DE_MDE_SIGN
programme:<tau>
file:<npz>
VCC_RC_CALLSET_CPM
VCC_DE_PRIOR_GATE
def sharpen_q_by_tercile
"""
UV_STUB = r'''#!~/Code/virtual-cell-challenge-2026/.venv/bin/python
"""Stand-in for `uv`. Records argv and VCC_* then writes the artifacts def_day.sh waits on.
Helpers and relevel_iter run for real, so the patched and unpatched arithmetic are the ones under test.
"""
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
            elif p[0] == "programme" and len(p) in (2, 3):
                cfg["de_mde_sign"] = {"source": "programme", "tau": float(p[1]), "tau_c4": float(p[2]) if len(p) == 3 else 0.05}
            elif p[0] == "file" and len(p) in (4, 5, 6):
                ms = {"source": "file", "sha256": p[2], "tau": float(p[3]), "tau_c4": float(p[4]) if len(p) >= 5 else 0.05}
                if len(p) == 6:
                    ms["le"] = int(p[5][2:])
                cfg["de_mde_sign"] = ms
            else:
                cfg["de_mde_sign"] = {"source": "unparsed", "raw": sign}
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
    out.with_suffix(".json").write_text(json.dumps({"hi": 3.0}))
    sys.exit(0)

if script.endswith("def_strength_key522.py"):
    out = Path(flag("--out"))
    out.write_bytes(b"strength")
    out.with_suffix(".json").write_text(json.dumps({"theta": float(flag("--theta"))}))
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
    return {str(p): sha(p) for p in LIVE.values()}


def copy_tree(dest: Path) -> None:
    for rel in REL.values():
        target = dest / rel
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy(ROOT / rel, target)


def run_patch(dest: Path, check: bool = False, extra: list[str] | None = None) -> subprocess.CompletedProcess:
    cmd = [str(PY), str(PATCH), "--root", str(dest)]
    if check:
        cmd.append("--check")
    if extra:
        cmd.extend(extra)
    return subprocess.run(cmd, capture_output=True, text=True)


def rho_star(scale: float = 1.0) -> float:
    import math
    from scipy.stats import norm
    q = 0.31 / (0.31 + 1.78)
    base = math.exp(2.0 * 0.01 * float(norm.ppf(q)))
    return base if scale == 1.0 else base * scale


class Bundle530(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls._live = live_hashes()
        cls.tmp = Path(tempfile.mkdtemp(prefix="bundle530-"))
        cls.repo = cls.tmp / "repo"
        copy_tree(cls.repo)
        r = run_patch(cls.repo)
        if r.returncode != 0:
            raise RuntimeError(r.stdout + r.stderr)
        ref = cls.repo / "submissions/453dr-october-path-abc-q1.75.h5ad.meta.json"
        ref.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy(REF_META, ref)

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.tmp, ignore_errors=True)

    def test_check_does_not_write(self):
        dest = self.tmp / "check"
        copy_tree(dest)
        before = {p: sha(dest / p) for p in REL.values()}
        r = run_patch(dest, check=True)
        self.assertEqual(r.returncode, 0, r.stderr)
        self.assertIn("patches apply cleanly (not written)", r.stdout)
        self.assertEqual(before, {p: sha(dest / p) for p in REL.values()})

    def test_second_apply_refuses_and_keeps_bytes(self):
        dest = self.tmp / "twice"
        copy_tree(dest)
        first = run_patch(dest)
        self.assertEqual(first.returncode, 0, first.stderr)
        digest = {p: sha(dest / p) for p in REL.values()}
        second = run_patch(dest)
        self.assertNotEqual(second.returncode, 0)
        self.assertIn("already patched", second.stderr + second.stdout)
        self.assertEqual(digest, {p: sha(dest / p) for p in REL.values()})
        again = run_patch(dest, check=True)
        self.assertIn("already patched", again.stderr + again.stdout)
        self.assertEqual(digest, {p: sha(dest / p) for p in REL.values()})

    def test_partial_bundle_refuses(self):
        dest = self.tmp / "partial"
        copy_tree(dest)
        helpers = dest / REL["helpers"]
        helpers.write_text(helpers.read_text() + "\ndef _october530_spec():\n    return None\n")
        before = {p: sha(dest / p) for p in REL.values()}
        r = run_patch(dest)
        self.assertNotEqual(r.returncode, 0)
        self.assertIn("partial bundle", r.stderr + r.stdout)
        self.assertIn("helpers", r.stderr + r.stdout)
        self.assertEqual(before, {p: sha(dest / p) for p in REL.values()})

    def test_applies_with_and_without_526_insertion(self):
        # coordinator 2026-09-29: 526's refusal is now at HEAD (600adea8), so the live copy carries it and the OTHER variant
        # is the one without it. Test whichever variant the live tree is not.
        bare = self.repo / REL["def_day"]
        self.assertIn("DEF_DAY_RHO_SCALE", bare.read_text())
        live_has_526 = "_HEADB=" in bare.read_text()
        dest = self.tmp / ("without526" if live_has_526 else "with526")
        copy_tree(dest)
        sh = dest / REL["def_day"]
        anchor = 'log "start: bundle $B work $W name $N q $Q sameline w $SLW flags $SLF; builder ${BSHA0:0:16} driver $(shasum -a 256 scripts/def_pipeline.py | cut -c1-16)"\n'
        text = sh.read_text()
        self.assertEqual(text.count(anchor), 1)
        if live_has_526:
            self.assertEqual(text.count(REFUSE526), 1)
            sh.write_text(text.replace(REFUSE526, "", 1))
        else:
            sh.write_text(text.replace(anchor, anchor + REFUSE526, 1))
        r = run_patch(dest)
        self.assertEqual(r.returncode, 0, r.stderr)
        out = sh.read_text()
        self.assertEqual(out.count(REFUSE526), 0 if live_has_526 else 1)
        self.assertIn("DEF_DAY_RHO_SCALE", out)
        self.assertIn('export VCC_DE_MDE_SIGN="$MDS"', out)
        bash = subprocess.run(["bash", "-n", str(sh)], capture_output=True, text=True)
        self.assertEqual(bash.returncode, 0, bash.stderr)

    def test_bash_n_and_strength_regex(self):
        sh = self.repo / REL["def_day"]
        bash = subprocess.run(["bash", "-n", str(sh)], capture_output=True, text=True)
        self.assertEqual(bash.returncode, 0, bash.stderr)
        text = sh.read_text()
        # The generated shell source must contain a raw regex for a decimal point, not a literal backslash.
        self.assertIn('re.fullmatch(r"-?(?:0|[1-9][0-9]*)(?:\\.[0-9]+)?", t)', text)
        self.assertNotIn("(?:\\\\.[0-9]+)?", text)
        body = _between(text, "python3 -c 'import re,sys\nt=sys.argv[1]\n", "' \"$DEF_DAY_LODO_STRENGTH\"")
        for value, ok in [("0.25", True), ("-0.25", True), ("1", True), ("-1", True), ("0.9", True),
                          ("0", False), ("0.0", False), ("1.1", False), ("1e-1", False), (".5", False), ("00.25", False)]:
            r = subprocess.run(["python3", "-c", body, value])
            self.assertEqual(r.returncode == 0, ok, value)
        hook = _between(text, "python3 -c 'import re,sys\np=sys.argv[1].split", "' \"$DEF_DAY_MDE_SIGN\"")
        for value, ok in [("hookbank:h1_2025_full:0.05", True), ("hookbank:b:0.05:0.04", True),
                          ("hookbank:b:0.05:0.04:le2", True), ("hookbank:b:0.05:0.04:le", False),
                          ("hookbank::0.05", False), ("hookbank:b:0.6", False), ("clean4:0.05", False)]:
            r = subprocess.run(["python3", "-c", hook, value])
            self.assertEqual(r.returncode == 0, ok, value)

    def test_decisions_new_switches_are_pending_and_cd4_stays_x1(self):
        d = json.loads((self.repo / REL["decisions"]).read_text())
        sw = d["switches"]
        self.assertEqual(sw["cd4_stack"]["status"], "DECIDED")
        self.assertEqual(sw["cd4_stack"]["value"], "x1")
        self.assertIn("DEF_DAY_CD4_STACK", sw["cd4_stack"]["evidence"])
        for key in ("rho_scale", "rc_callset_cpm", "hook_estimator", "lodo_strength"):
            self.assertEqual(sw[key]["status"], "PENDING", key)
            self.assertIsNone(sw[key]["value"])
        self.assertEqual(sw["mde_sign"]["status"], "PENDING")
        self.assertIn("hookbank:", sw["mde_sign"]["if_hookbank"])
        # insertion kept the day-of note and did not reorder the decided October command inputs
        self.assertEqual(sw["q"]["value"], "1.75")
        self.assertEqual(sw["lodo"]["value"], "1")
        self.assertEqual(sw["de_mag"]["value"], "dn_k250_zs")

    def test_resolver_stock_allow_pending_is_todays_command(self):
        r = self._recipe("command", "--allow-pending", "--bundle", "B", "--work", "W", "--name", "N")
        self.assertEqual(r.returncode, 0, r.stderr)
        lines = [ln for ln in r.stdout.splitlines() if ln and not ln.startswith("#")]
        self.assertEqual(lines, ["DEF_DAY_LODO=1 DEF_DAY_DE_MAG=dn_k250_zs bash scripts/def_day.sh B W N 1.75 0.5 noprog 4"])
        self.assertNotIn("DEF_DAY_CD4_STACK", r.stdout)
        self.assertNotIn("DEF_DAY_RHO_SCALE", lines[0])
        refused = self._recipe("command", "--bundle", "B", "--work", "W", "--name", "N")
        self.assertEqual(refused.returncode, 2)
        self.assertIn("# REFUSED", refused.stdout)
        self.assertIn("PENDING", refused.stderr)

    def test_resolver_every_combination_of_the_new_switches(self):
        n = 0
        for rho, cd4, rc, hook, strength in itertools.product((False, True), repeat=5):
            dest_dec = self._decide(rho=rho, cd4_x15=cd4, rc=rc, hook=hook, strength=strength)
            proc = self._recipe("command", "--allow-pending", "--json", "--bundle", "B", "--work", "W", "--name", "N",
                                decisions=dest_dec)
            self.assertEqual(proc.returncode, 0, proc.stderr + proc.stdout)
            got = json.loads(proc.stdout)
            self.assertEqual(got["problems"], [], (rho, cd4, rc, hook, strength, got["problems"]))
            self.assertEqual(list(got["env"]), _expect_keys(rho, cd4, rc, hook, strength))
            env = got["env"]
            self.assertEqual(env["DEF_DAY_LODO"], "1")
            self.assertEqual(env["DEF_DAY_DE_MAG"], "dn_k250_zs")
            if rho:
                self.assertEqual(env["DEF_DAY_RHO_SCALE"], "1.02")
            if cd4:
                self.assertEqual(env["DEF_DAY_CD4_STACK"], "x1.5")
            else:
                self.assertNotIn("DEF_DAY_CD4_STACK", env)
            if rc:
                self.assertEqual(env["DEF_DAY_RC_CALLSET_CPM"], "1500")
            if hook:
                self.assertEqual(env["DEF_DAY_HOOK_ESTIMATOR"], "meanz_ne")
                self.assertEqual(env["DEF_DAY_MDE_SIGN"], "hookbank:h1_2025_full:0.05")
            if strength:
                self.assertEqual(env["DEF_DAY_LODO_STRENGTH"], "0.25")
                self.assertEqual(env["DEF_DAY_LODO_HI_PCT"], "96.64697522769474")
            self.assertEqual(got["args"][3:], ["1.75", "0.5", "noprog", "4"])
            n += 1
        self.assertEqual(n, 32)
        print(f"resolver combinations {n} ok", flush=True)
        # boundaries that are on, but not part of the 32
        for value in ("0.9", "1.1", "1"):
            proc = self._recipe("command", "--allow-pending", "--json", "--bundle", "B", "--work", "W", "--name", "N",
                                decisions=self._decide(rho_value=value))
            self.assertEqual(proc.returncode, 0, proc.stderr)
            self.assertEqual(json.loads(proc.stdout)["env"]["DEF_DAY_RHO_SCALE"], value)
        for value in ("100", "5000", "2200"):
            proc = self._recipe("command", "--allow-pending", "--json", "--bundle", "B", "--work", "W", "--name", "N",
                                decisions=self._decide(rc_value=value))
            self.assertEqual(proc.returncode, 0, proc.stderr)
            self.assertEqual(json.loads(proc.stdout)["env"]["DEF_DAY_RC_CALLSET_CPM"], value)

    def test_resolver_rejects_bad_values_and_unpatched_cd4(self):
        cases = [
            ({"rho_value": "0.89"}, "rho_scale"),
            ({"rho_value": "1.11"}, "rho_scale"),
            ({"rc_value": "99"}, "rc_callset_cpm"),
            ({"rc_value": "5001"}, "rc_callset_cpm"),
            ({"rc_value": "1500.5"}, "rc_callset_cpm"),
            ({"cd4_value": "x2"}, "cd4_stack"),
            ({"hook": True, "hook_value": "pdown", "mde": "hookbank:b:0.05"}, "hook_estimator"),
            ({"hook": True, "mde": "clean4:0.05"}, "hookbank"),
            ({"hook": True, "mde": "programme:0.05"}, "hookbank"),
            ({"strength": True, "strength_value": "0"}, "lodo_strength"),
            ({"strength": True, "strength_value": "1e-1"}, "lodo_strength"),
            ({"strength": True, "strength_value": "1.01"}, "lodo_strength"),
            ({"strength": True, "strength_value": "-0.25", "with_lodo": False}, "lodo = 1"),
            ({"strength": True, "strength_value": "0.25", "with_pct": False}, "lodo_hi_pct"),
            ({"mde": "file:a:b:0.05"}, "mde_sign"),
        ]
        for kwargs, needle in cases:
            proc = self._recipe("command", "--allow-pending", "--bundle", "B", "--work", "W", "--name", "N",
                                decisions=self._decide(**kwargs))
            self.assertEqual(proc.returncode, 3, (kwargs, proc.returncode, proc.stderr, proc.stdout))
            self.assertIn("# REFUSED", proc.stdout)
            self.assertIn(needle, proc.stderr)
        # x1.5 is legal only once def_day.sh contains the dump override. The unpatched script still refuses.
        raw = self.tmp / "rawday"
        copy_tree(raw)
        proc = self._recipe("command", "--allow-pending", "--bundle", "B", "--work", "W", "--name", "N",
                            decisions=self._decide(cd4_x15=True), def_day=raw / REL["def_day"])
        self.assertEqual(proc.returncode, 3, proc.stderr)
        self.assertIn("cd4_stack", proc.stderr)
        # a patched day accepts it
        proc = self._recipe("command", "--allow-pending", "--json", "--bundle", "B", "--work", "W", "--name", "N",
                            decisions=self._decide(cd4_x15=True))
        self.assertEqual(proc.returncode, 0, proc.stderr)
        self.assertEqual(json.loads(proc.stdout)["env"]["DEF_DAY_CD4_STACK"], "x1.5")

    def test_riders_match_accepts_and_rejects(self):
        helpers = str(self.repo / REL["helpers"])
        star = self.tmp / "metas"
        star.mkdir(exist_ok=True)

        def meta(name, cfg):
            p = star / name
            p.write_text(json.dumps({"config": cfg}))
            return str(p)

        def run(*args):
            return subprocess.run([str(PY), helpers, "riders_match", *args], capture_output=True, text=True)

        base = {"bulk_sharpen_q": 1.75, "cd4_target_g": X1_G, "rc_callset_cpm": 2200, "bulk_load": "dump_h1x2_cd4x1.npz"}
        sha64 = "a" * 64
        accept = [
            (meta("none.json", dict(base)), "-", "-"),
            (meta("mag.json", {**base, "de_mag": {"variant": "dn_k250_zs"}}), "dn_k250_zs", "-"),
            (meta("c4.json", {**base, "de_mde_sign": {"source": "clean4", "tau": 0.05}}), "-", "clean4:0.05"),
            (meta("prog.json", {**base, "de_mde_sign": {"source": "programme", "tau": 0.05, "tau_c4": 0.05}}), "-", "programme:0.05"),
            (meta("prog2.json", {**base, "de_mde_sign": {"source": "programme", "tau": 0.05, "tau_c4": 0.04}}), "-", "programme:0.05:0.04"),
            (meta("file.json", {**base, "de_mde_sign": {"source": "file", "sha256": sha64, "tau": 0.05, "tau_c4": 0.05}}), "-", f"file:p:{sha64}:0.05"),
            (meta("filele.json", {**base, "de_mde_sign": {"source": "file", "sha256": sha64, "tau": 0.05, "tau_c4": 0.04, "le": 2}}),
             "-", f"file:p:{sha64}:0.05:0.04:le2"),
            (meta("oct.json", dict(base)), "-", "-", "rc=-,stack=x1,rho=-,strength=-"),
            (meta("oct15.json", {**base, "bulk_load": "dump_h1x2_cd4x15_s0.25.npz", "rc_callset_cpm": 1500}),
             "-", "-", "rc=1500,stack=x1.5,rho=1.02,strength=0.25"),
            (meta("rho.json", {**base, "rho_scale": 1.02}), "-", "-", "rc=-,stack=x1,rho=1.02,strength=-"),
            # a 3-arg call does not grow a dump check: a renamed bulk_load is still today's riders_match
            (meta("renamed.json", {**base, "bulk_load": "dump_h1x2_cd4x15.npz"}), "-", "-"),
        ]
        for args in accept:
            got = run(*args)
            self.assertEqual(got.returncode, 0, (args, got.stderr))
        reject = [
            (meta("badmag.json", {**base, "de_mag": {"variant": "other"}}), "dn_k250_zs", "-"),
            (meta("signpresent.json", {**base, "de_mde_sign": {"source": "clean4", "tau": 0.05}}), "-", "-"),
            (meta("signabsent.json", dict(base)), "-", "clean4:0.05"),
            (meta("badtau.json", {**base, "de_mde_sign": {"source": "clean4", "tau": 0.40}}), "-", "clean4:0.05"),
            (meta("badc4.json", {**base, "de_mde_sign": {"source": "programme", "tau": 0.05, "tau_c4": 0.05}}), "-", "programme:0.05:0.40"),
            (meta("badsha.json", {**base, "de_mde_sign": {"source": "file", "sha256": "b" * 64, "tau": 0.05, "tau_c4": 0.05}}),
             "-", f"file:p:{sha64}:0.05"),
            (meta("badle.json", {**base, "de_mde_sign": {"source": "file", "sha256": sha64, "tau": 0.05, "tau_c4": 0.04, "le": 1}}),
             "-", f"file:p:{sha64}:0.05:0.04:le2"),
            (meta("hookraw.json", dict(base)), "-", "hookbank:b:0.05"),
            (meta("baddump.json", {**base, "bulk_load": "dump_h1x2_cd4x15.npz"}), "-", "-", "rc=-,stack=x1,rho=-,strength=-"),
            (meta("leak.json", {**base, "bulk_load": "dump_h1x2_cd4x15.npz", "cd4_target_g": float(X15_G)}),
             "-", "-", "rc=-,stack=x1.5,rho=-,strength=-"),
            (meta("badrc.json", {**base, "rc_callset_cpm": 2200}), "-", "-", "rc=1500,stack=x1,rho=-,strength=-"),
            (meta("badrho.json", {**base, "rho_scale": 1.1}), "-", "-", "rc=-,stack=x1,rho=1.02,strength=-"),
            (meta("sm.json", {**base, "bulk_load": "dump_h1x2_cd4x1_s0.25.npz"}), "-", "-", "rc=-,stack=x1,rho=-,strength=-0.25"),
        ]
        for args in reject:
            got = run(*args)
            self.assertNotEqual(got.returncode, 0, args)
        badspec = run(meta("spec.json", dict(base)), "-", "-", "rc=-")
        self.assertEqual(badspec.returncode, 2, badspec.stderr)

    def test_relevel_and_helper_scale(self):
        import importlib.util
        unpatched = self.tmp / "relevel-src"
        copy_tree(unpatched)
        spec = importlib.util.spec_from_file_location("helpers530", self.repo / REL["helpers"])
        helpers = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(helpers)
        spec_u = importlib.util.spec_from_file_location("helpers_u", unpatched / REL["helpers"])
        old = importlib.util.module_from_spec(spec_u)
        spec_u.loader.exec_module(old)
        self.assertEqual(helpers.rho_star(), old.rho_star())
        self.assertEqual(helpers.rho_star(scale=1.0), helpers.rho_star())
        self.assertEqual(helpers.rho_star(scale=1.02), helpers.rho_star() * 1.02)
        src = (self.repo / REL["helpers"]).read_text()
        self.assertIn("return base if scale == 1.0 else base * scale", src)
        self.assertEqual(helpers._dump_name530("x1", "-"), "dump_h1x2_cd4x1.npz")
        self.assertEqual(helpers._dump_name530("x1.5", "-"), "dump_h1x2_cd4x15.npz")
        self.assertEqual(helpers._dump_name530("x1.5", "0.25"), "dump_h1x2_cd4x15_s0.25.npz")
        self.assertEqual(helpers._dump_name530("x1", "-0.25"), "dump_h1x2_cd4x1_sm0.25.npz")

        points = self.tmp / "points"
        points.mkdir(exist_ok=True)
        nt = points / "a.ntclevel.json"
        nt.write_text(json.dumps({"per_context": {c: {"rho": 0.5, "scale": 1.1} for c in "DEF"}}))
        spec_arg = f"{nt}:D=1,E=1,F=1"
        new_re = str(self.repo / REL["relevel"])
        old_re = str(unpatched / REL["relevel"])
        bare = subprocess.run([str(PY), new_re, spec_arg], capture_output=True, text=True)
        one = subprocess.run([str(PY), new_re, spec_arg, "--rho-scale", "1"], capture_output=True, text=True)
        old_out = subprocess.run([str(PY), old_re, spec_arg], capture_output=True, text=True)
        self.assertEqual(bare.returncode, 0, bare.stderr)
        self.assertEqual(bare.stdout, one.stdout)
        self.assertEqual(bare.stdout, old_out.stdout)
        scaled = subprocess.run([str(PY), new_re, spec_arg, "--rho-scale", "1.02"], capture_output=True, text=True)
        self.assertEqual(scaled.returncode, 0, scaled.stderr)
        self.assertNotEqual(scaled.stdout, bare.stdout)
        self.assertIn(f"{rho_star(1.02):.4f}", scaled.stdout)
        self.assertIn(f"{rho_star():.4f}", bare.stdout)
        for bad in ("0.5", "1.2"):
            r = subprocess.run([str(PY), new_re, spec_arg, "--rho-scale", bad], capture_output=True, text=True)
            self.assertNotEqual(r.returncode, 0, bad)

        def write_nt(path: Path, rho: float):
            path.write_text(json.dumps({"per_context": {c: {"rho": rho, "scale": 1.0} for c in "DEF"}}))
            path.with_name(path.name.replace(".ntclevel.json", ".h5ad")).write_bytes(b"h5")

        at_star = points / "star.ntclevel.json"
        at_scaled = points / "scaled.ntclevel.json"
        write_nt(at_star, rho_star())
        write_nt(at_scaled, rho_star(1.02))
        self.assertLess(rho_star(1.02), 1.0)
        h = str(self.repo / REL["helpers"])

        def inband(path, scale=None):
            cmd = [str(PY), h, "inband", str(path)] + ([str(scale)] if scale is not None else [])
            return subprocess.run(cmd).returncode

        self.assertEqual(inband(at_star), 0)
        self.assertEqual(inband(at_star, "1.0"), 0)
        self.assertEqual(inband(at_star, "1.02"), 1)
        self.assertEqual(inband(at_scaled), 1)
        self.assertEqual(inband(at_scaled, "1.02"), 0)
        # one-arg inband matches the unpatched helper
        old_h = str(unpatched / REL["helpers"])
        self.assertEqual(subprocess.run([str(PY), old_h, "inband", str(at_star)]).returncode, 0)

        def best(scale, *paths):
            cmd = [str(PY), h, "best_safe"] + ([ "--scale", str(scale)] if scale is not None else []) + list(paths)
            return subprocess.run(cmd, capture_output=True, text=True).stdout.strip()

        a = f"{at_star}:D=1,E=1,F=1"
        b = f"{at_scaled}:D=1,E=1,F=1"
        self.assertTrue(best(None, a, b).endswith("star"))
        self.assertTrue(best("1.02", a, b).endswith("scaled"))
        self.assertEqual(best("1.0", a, b), best(None, a, b))

    def test_gate_rows(self):
        gate = str(self.repo / REL["gate"])
        sha64 = "ab" * 32

        def run(env, *, meta_edit=None, log_extra="", sidecar=None, drop_sameline=False):
            work = Path(tempfile.mkdtemp(prefix="gate530-", dir=self.tmp))
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
                (run_dir / "pert_scale_lodo_strength.json").write_text(json.dumps(sidecar))
            meta = json.loads(REF_META.read_text())
            if meta_edit:
                meta_edit(meta["config"])
            if drop_sameline:
                meta["config"].pop("sameline", None)
            stem = work / "pkg"
            Path(str(stem) + ".h5ad.meta.json").write_text(json.dumps(meta))
            blob = b"vcc-synthetic"
            Path(str(stem) + ".vcc").write_bytes(blob)
            digest = hashlib.sha256(blob).hexdigest()
            chain = work / "chain.log"
            chain.write_text(f"novelty rc=0\npackage pkg.vcc sha256 {digest}\n")
            proc = subprocess.run([str(PY), gate, "--write", "--work", str(work), "--stem", str(stem), "--chainlog", str(chain)],
                                  capture_output=True, text=True)
            return proc

        # minimal: no new knob. The off-check requires today's dump name, which 453dr has.
        base = run({})
        self.assertEqual(base.returncode, 0, base.stdout + base.stderr)
        self.assertIn("ALL GATES PASS", base.stdout)
        dropped = run({}, drop_sameline=True)
        self.assertEqual(dropped.returncode, 1)
        self.assertIn("ALL GATES FAIL", dropped.stdout)
        self.assertIn("same-line record", dropped.stdout)

        def sign(cfg, **ms):
            cfg["de_mde_sign"] = ms

        ok = run({"DEF_DAY_MDE_SIGN": "clean4:0.05"}, meta_edit=lambda c: sign(c, source="clean4", tau=0.05))
        self.assertEqual(ok.returncode, 0, ok.stdout)
        bad_tau = run({"DEF_DAY_MDE_SIGN": "clean4:0.05"}, meta_edit=lambda c: sign(c, source="clean4", tau=0.40))
        self.assertIn("ALL GATES FAIL", bad_tau.stdout)
        self.assertIn("uncalled-gene sign hook", bad_tau.stdout)
        # programme's tau_c4 is part of the value, not dropped by a one-colon partition
        prog_bad = run({"DEF_DAY_MDE_SIGN": "programme:0.05:0.40"},
                       meta_edit=lambda c: sign(c, source="programme", tau=0.05, tau_c4=0.05))
        self.assertIn("ALL GATES FAIL", prog_bad.stdout)
        prog_ok = run({"DEF_DAY_MDE_SIGN": "programme:0.05:0.40"},
                      meta_edit=lambda c: sign(c, source="programme", tau=0.05, tau_c4=0.40))
        self.assertEqual(prog_ok.returncode, 0, prog_ok.stdout)
        prog_default = run({"DEF_DAY_MDE_SIGN": "programme:0.05"},
                           meta_edit=lambda c: sign(c, source="programme", tau=0.05, tau_c4=0.05))
        self.assertEqual(prog_default.returncode, 0, prog_default.stdout)

        hook_env = {
            "DEF_DAY_MDE_SIGN": "hookbank:h1_2025_full:0.05:0.04:le2",
            "DEF_DAY_HOOK_ESTIMATOR": "meanz_ne",
            "DEF_DAY_RHO_SCALE": "1.02",
            "DEF_DAY_CD4_STACK": "x1.5",
            "DEF_DAY_RC_CALLSET_CPM": "1500",
            "DEF_DAY_LODO_STRENGTH": "0.25",
        }

        def all_edit(cfg):
            cfg["bulk_load"] = "dump_h1x2_cd4x15_s0.25.npz"
            cfg["bulk_load_sha256"] = "c" * 64
            cfg["rc_callset_cpm"] = 1500
            cfg["de_mde_sign"] = {"source": "file", "sha256": sha64, "tau": 0.05, "tau_c4": 0.04, "le": 2}

        all_ok = run(hook_env, meta_edit=all_edit, sidecar={"theta": 0.25},
                     log_extra="re-level rho-scale=1.02\nLODO strength rider: theta 0.25\nhook estimator meanz_ne\n")
        self.assertEqual(all_ok.returncode, 0, all_ok.stdout + all_ok.stderr)
        self.assertIn("ALL GATES PASS", all_ok.stdout)

        # 1.02 must not match a log line that only names 1.020
        rho_prefix = run({"DEF_DAY_RHO_SCALE": "1.02"}, log_extra="re-level rho-scale=1.020\n")
        self.assertIn("ALL GATES FAIL", rho_prefix.stdout)
        self.assertIn("rho scale", rho_prefix.stdout)
        rho_ok = run({"DEF_DAY_RHO_SCALE": "1.02"}, log_extra="re-level rho-scale=1.02\n")
        self.assertEqual(rho_ok.returncode, 0, rho_ok.stdout)

        leaked = run({"DEF_DAY_CD4_STACK": "x1.5"},
                     meta_edit=lambda c: c.update(bulk_load="dump_h1x2_cd4x15.npz", cd4_target_g=float(X15_G)))
        self.assertIn("ALL GATES FAIL", leaked.stdout)
        self.assertIn("stack dump name", leaked.stdout)
        silent = run({}, meta_edit=lambda c: c.update(bulk_load="dump_h1x2_cd4x15.npz"))
        self.assertIn("ALL GATES FAIL", silent.stdout)
        self.assertIn("CD4 stack x1", silent.stdout)
        wrong_rc = run({"DEF_DAY_RC_CALLSET_CPM": "1500"})
        self.assertIn("call-set CPM", wrong_rc.stdout)
        self.assertIn("ALL GATES FAIL", wrong_rc.stdout)
        hook_src = run({"DEF_DAY_MDE_SIGN": "hookbank:b:0.05", "DEF_DAY_HOOK_ESTIMATOR": "meanz_ne"},
                       meta_edit=lambda c: sign(c, source="hookbank", sha256=sha64, tau=0.05, tau_c4=0.05),
                       log_extra="estimator meanz_ne\n")
        self.assertIn("ALL GATES FAIL", hook_src.stdout)
        no_le = run({"DEF_DAY_MDE_SIGN": "hookbank:b:0.05:0.04:le2"},
                    meta_edit=lambda c: sign(c, source="file", sha256=sha64, tau=0.05, tau_c4=0.04))
        self.assertIn("ALL GATES FAIL", no_le.stdout)
        extra_le = run({"DEF_DAY_MDE_SIGN": "hookbank:b:0.05"},
                       meta_edit=lambda c: sign(c, source="file", sha256=sha64, tau=0.05, tau_c4=0.05, le=2))
        self.assertIn("ALL GATES FAIL", extra_le.stdout)
        strength_miss = run({"DEF_DAY_LODO_STRENGTH": "0.25"},
                            meta_edit=lambda c: c.update(bulk_load="dump_h1x2_cd4x1_s0.25.npz"),
                            log_extra="LODO strength rider: theta 0.25\n")
        self.assertIn("ALL GATES FAIL", strength_miss.stdout)
        neg = run({"DEF_DAY_LODO_STRENGTH": "-0.25"},
                  meta_edit=lambda c: c.update(bulk_load="dump_h1x2_cd4x1_sm0.25.npz"),
                  sidecar={"theta": -0.25}, log_extra="LODO strength rider: theta -0.25\n")
        self.assertEqual(neg.returncode, 0, neg.stdout + neg.stderr)

    def test_unset_matches_todays_exports_and_arguments(self):
        star = f"{rho_star():.16g}"
        bare_old = self._day(patched=False, env={"VCC_RC_CALLSET_CPM": "9999"}, rho=star)
        bare_new = self._day(patched=True, env={"VCC_RC_CALLSET_CPM": "9999"}, rho=star)
        self._assert_day_ok(bare_old)
        self._assert_day_ok(bare_new)
        self.assertEqual(bare_old["records"], bare_new["records"])
        self.assertEqual(bare_old["final"], bare_new["final"])
        self._assert_today(bare_old)
        october = {"DEF_DAY_LODO": "1", "DEF_DAY_DE_MAG": "dn_k250_zs", "VCC_RC_CALLSET_CPM": "9999"}
        oct_old = self._day(patched=False, env=october, rho=star)
        oct_new = self._day(patched=True, env=october, rho=star)
        self._assert_day_ok(oct_old)
        self._assert_day_ok(oct_new)
        self.assertEqual(oct_old["records"], oct_new["records"])
        self.assertEqual(oct_old["final"], oct_new["final"])
        self.assertIn("VCC_DE_MAG=dn_k250_zs", oct_new["final"])
        self.assertIn("VCC_PERT_SCALE", json.dumps(oct_new["records"]))
        # explicit x1 is the same command as leaving the knob unset
        x1 = self._day(patched=True, env={"DEF_DAY_CD4_STACK": "x1", "VCC_RC_CALLSET_CPM": "9999"}, rho=star)
        self._assert_day_ok(x1)
        self.assertEqual(bare_new["records"], x1["records"])
        self.assertEqual(bare_new["final"], x1["final"])
        print(f"unset identity: {len(bare_new['records'])} uv calls match; october {len(oct_new['records'])} match", flush=True)

    def test_knobs_on_change_only_their_exports(self):
        rho = f"{rho_star(1.02):.16g}"
        env = {
            "DEF_DAY_LODO": "1",
            "DEF_DAY_LODO_HI_PCT": "96.64697522769474",
            "DEF_DAY_RHO_SCALE": "1.02",
            "DEF_DAY_CD4_STACK": "x1.5",
            "DEF_DAY_RC_CALLSET_CPM": "1500",
            "DEF_DAY_MDE_SIGN": "hookbank:h1bank:0.05:0.04:le2",
            "DEF_DAY_HOOK_ESTIMATOR": "meanz_ne",
            "DEF_DAY_LODO_STRENGTH": "0.25",
            "VCC_RC_CALLSET_CPM": "9999",
        }
        got = self._day(patched=True, env=env, rho=rho)
        self._assert_day_ok(got)
        dump = [r for r in got["raw"] if _name(r) == "day-dump"]
        self.assertEqual(len(dump), 1)
        self.assertTrue(dump[0]["vcc"]["VCC_BULK_DUMP"].endswith("dump_h1x2_cd4x15_s0.25.npz"))
        self.assertEqual(dump[0]["vcc"]["VCC_CD4_TARGET_G"], X15_G)
        self.assertEqual(dump[0]["vcc"]["VCC_H1_UNION_G"], "0.56")
        self.assertNotIn("VCC_RC_CALLSET_CPM", dump[0]["vcc"])
        self.assertNotIn("VCC_DE_MDE_SIGN", dump[0]["vcc"])
        levels = [r for r in got["raw"] if _name(r) == "day-lvl0"]
        self.assertEqual(len(levels), 1)
        self.assertEqual(levels[0]["vcc"]["VCC_RC_CALLSET_CPM"], "1500")
        self.assertNotIn("VCC_CD4_TARGET_G", levels[0]["vcc"])
        self.assertTrue(levels[0]["vcc"]["VCC_DE_MDE_SIGN"].startswith("file:"))
        self.assertIn(":0.05:0.04:le2", levels[0]["vcc"]["VCC_DE_MDE_SIGN"])
        self.assertNotIn("9999", json.dumps(got["raw"]))
        relevel = [r for r in got["raw"] if "relevel_iter.py" in " ".join(r["argv"])]
        self.assertEqual(relevel[-1]["argv"][-2:], ["--rho-scale", "1.02"])
        strength = [r for r in got["raw"] if "def_strength_key522.py" in " ".join(r["argv"])]
        self.assertEqual(len(strength), 1)
        self.assertIn("0.25", strength[0]["argv"])
        pin = [r for r in got["raw"] if "def_pin_panel_artifacts.py" in " ".join(r["argv"]) and "--apply" in r["argv"]]
        self.assertTrue(any(a.endswith("pert_scale_lodo_strength.npz") for a in pin[0]["argv"]))
        sameline = [r for r in got["raw"] if r["argv"][2].endswith("def_sameline.py")]
        self.assertEqual(len(sameline), 2)
        self.assertIn("h1bank", sameline[1]["argv"])
        self.assertIn("meanz_ne", sameline[1]["argv"])
        riders = [r for r in got["raw"] if "riders_match" in r["argv"]]
        self.assertEqual(riders[-1]["argv"][-1], "rc=1500,stack=x1.5,rho=1.02,strength=0.25")
        self.assertIn("DEF_DAY_RHO_SCALE=1.02", got["final"])
        self.assertIn("DEF_DAY_CD4_STACK=x1.5", got["final"])
        self.assertIn("VCC_RC_CALLSET_CPM=1500", got["final"])
        self.assertIn("DEF_DAY_LODO_STRENGTH=0.25", got["final"])
        self.assertIn("DEF_DAY_HOOK_ESTIMATOR=meanz_ne", got["final"])
        log = got["log"]
        self.assertIn("re-level rho-scale=1.02", log)
        self.assertIn("LODO strength rider: theta 0.25", log)
        self.assertIn("estimator meanz_ne", log)
        self.assertIn("CD4 stack x1.5", log)

    def test_bad_scale_stops_before_any_build(self):
        got = self._day(patched=True, env={"DEF_DAY_RHO_SCALE": "1.2"}, rho=f"{rho_star():.16g}")
        self.assertNotEqual(got["rc"], 0)
        self.assertIn("DEF_DAY_RHO_SCALE", got["log"])
        self.assertEqual(got["raw"], [])

    def test_zz_live_tree_untouched(self):
        self.assertEqual(self._live, live_hashes())


    def _decide(self, rho=False, cd4_x15=False, rc=False, hook=False, strength=False,
                rho_value="1.02", rc_value="1500", cd4_value="x1.5", hook_value="meanz_ne",
                strength_value="0.25", mde="hookbank:h1_2025_full:0.05", with_lodo=True, with_pct=True):
        d = json.loads((self.repo / REL["decisions"]).read_text())
        sw = d["switches"]
        if rho or rho_value != "1.02":
            sw["rho_scale"].update(status="DECIDED", value=rho_value)
        if cd4_x15 or cd4_value != "x1.5":
            sw["cd4_stack"].update(status="DECIDED", value=cd4_value)
        if rc or rc_value != "1500":
            sw["rc_callset_cpm"].update(status="DECIDED", value=rc_value)
        if hook or hook_value != "meanz_ne" or mde != "hookbank:h1_2025_full:0.05":
            sw["hook_estimator"].update(status="DECIDED", value=hook_value)
            sw["mde_sign"].update(status="DECIDED", value=mde)
        elif mde != "hookbank:h1_2025_full:0.05" and not hook:
            sw["mde_sign"].update(status="DECIDED", value=mde)
        if strength or strength_value != "0.25":
            sw["lodo_strength"].update(status="DECIDED", value=strength_value)
            if with_pct:
                sw["lodo_hi_pct"].update(status="DECIDED", value="96.64697522769474")
            if not with_lodo:
                sw["lodo"].update(status="DECIDED", value="0")
        # a caller that only wanted a bad mde (no hook) still needs the value decided
        if mde != "hookbank:h1_2025_full:0.05" and not hook and hook_value == "meanz_ne":
            sw["mde_sign"].update(status="DECIDED", value=mde)
            sw["hook_estimator"].update(status="PENDING", value=None)
        path = self.tmp / f"dec-{sha(self._write_tmp(json.dumps(d)))[:8]}.json"
        # _write_tmp already wrote it; return that path. The hash name above double-writes. Keep it simple:
        return path if path.is_file() else self._write_json(d)

    def _write_tmp(self, text: str) -> Path:
        p = self.tmp / f"blob-{hashlib.sha256(text.encode()).hexdigest()[:10]}.json"
        if not p.exists():
            p.write_text(text)
        return p

    def _write_json(self, d: dict) -> Path:
        return self._write_tmp(json.dumps(d))

    def _recipe(self, *args, decisions=None, def_day=None, sameline=None):
        env = os.environ.copy()
        env["OCTOBER_DECISIONS"] = str(decisions or (self.repo / REL["decisions"]))
        env["OCTOBER_DEF_DAY"] = str(def_day or (self.repo / REL["def_day"]))
        env["OCTOBER_BUILDER"] = str(ROOT / "scripts/make_submission_pm_cd4.py")
        env["OCTOBER_SAMELINE"] = str(sameline or (self.repo / REL["sameline"]))
        return subprocess.run([str(PY), str(self.repo / REL["recipe"]), *args], capture_output=True, text=True, env=env)

    def _day(self, patched: bool, env: dict, rho: str) -> dict:
        dest = Path(tempfile.mkdtemp(prefix="day530-", dir=self.tmp))
        repo = dest / "repo"
        if patched:
            for rel in (REL["def_day"], REL["helpers"], REL["relevel"], REL["sameline"]):
                target = repo / rel
                target.parent.mkdir(parents=True, exist_ok=True)
                shutil.copy(self.repo / rel, target)
        else:
            copy_tree(repo)
        for rel, text in {
            Path("scripts/make_submission_pm_cd4.py"): BUILDER_STUB,
            Path("scripts/def_pipeline.py"): "# stub\n",
            Path("scripts/def_lodo_key.py"): "# --hi-pct\n",
            Path("scripts/def_pin_panel_artifacts.py"): "# stub\n",
            Path("scripts/def_postbuild.py"): "# stub\n",
            Path("experiments/522-october-bulk-key/def_strength_key522.py"): "# stub\n",
            Path("experiments/417-bulk-sharpen/ntclevel417.py"): "# stub\n",
            Path("experiments/428-zero-gene-audit/audit428.py"): "# stub\n",
        }.items():
            p = repo / rel
            p.parent.mkdir(parents=True, exist_ok=True)
            p.write_text(text)
        if not patched:
            # copy_tree already placed the real sameline; leave it. The stub intercepts the call.
            pass
        sh = repo / REL["def_day"]
        sh.write_text(sh.read_text().replace(
            "REPO=~/Code/virtual-cell-challenge-2026", f"REPO={repo}"))
        # coordinator 2026-09-29: the live def_day.sh now carries 526's refusal (`git show HEAD:<builder>` must equal the
        # builder on disk), so the stub tree is committed in a throwaway repo, with the stub builder at HEAD.
        g = ["git", "-C", str(repo), "-c", "user.name=t", "-c", "user.email=t@t", "-c", "commit.gpgsign=false"]
        subprocess.run(["git", "init", "-q", str(repo)], check=True)
        subprocess.run(g + ["add", "-A"], check=True)
        subprocess.run(g + ["commit", "-q", "-m", "stub"], check=True)
        bindir = dest / "bin"
        bindir.mkdir()
        uv = bindir / "uv"
        uv.write_text(UV_STUB)
        uv.chmod(0o755)
        bundle = dest / "bundle"
        bundle.mkdir()
        (bundle / "manifest.json").write_text(json.dumps({"contexts": ["D", "E", "F"]}))
        work = dest / "work"
        log = dest / "uv.jsonl"
        run_env = os.environ.copy()
        for key in list(run_env):
            if key.startswith("VCC_"):
                del run_env[key]
        run_env.update(env)
        run_env["PATH"] = str(bindir) + os.pathsep + run_env.get("PATH", "")
        run_env["STUB_LOG"] = str(log)
        run_env["STUB_RHO"] = rho
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
                "final": _sub(final, str(work), str(bundle)), "log": day_log, "err": proc.stderr}

    def _assert_day_ok(self, got):
        if got["rc"] != 0:
            self.fail(f"rc {got['rc']}\n{got['log'][-2000:]}\n{got['err'][-1000:]}")

    def _assert_today(self, got):
        blob = json.dumps(got["records"])
        self.assertNotIn("VCC_RC_CALLSET_CPM", blob)
        self.assertNotIn("9999", blob)
        self.assertNotIn("--rho-scale", blob)
        self.assertNotIn("def_strength_key522.py", blob)
        self.assertNotIn("VCC_CD4_TARGET_G", blob)
        self.assertNotIn("--programme-estimator", blob)
        dumps = [r for r in got["records"] if _name(r) == "day-dump"]
        self.assertEqual(set(dumps[0]["vcc"]), {"VCC_H1_UNION_G", "VCC_BULK_DUMP"})
        self.assertEqual(dumps[0]["vcc"]["VCC_H1_UNION_G"], "0.56")
        self.assertTrue(dumps[0]["vcc"]["VCC_BULK_DUMP"].endswith("dump_h1x2_cd4x1.npz"))
        riders = [r for r in got["records"] if "riders_match" in r["argv"]]
        self.assertEqual(riders[-1]["argv"][-2:], ["-", "-"])
        for needle in ("DEF_DAY_RHO_SCALE=", "DEF_DAY_CD4_STACK=", "VCC_RC_CALLSET_CPM=",
                       "DEF_DAY_LODO_STRENGTH=", "DEF_DAY_HOOK_ESTIMATOR="):
            self.assertNotIn(needle, got["final"])


def _expect_keys(rho, cd4, rc, hook, strength) -> list[str]:
    keys = ["DEF_DAY_LODO"]
    if strength:
        keys.append("DEF_DAY_LODO_HI_PCT")
    keys.append("DEF_DAY_DE_MAG")
    if hook:
        keys.append("DEF_DAY_MDE_SIGN")
    if rho:
        keys.append("DEF_DAY_RHO_SCALE")
    if rc:
        keys.append("DEF_DAY_RC_CALLSET_CPM")
    if hook:
        keys.append("DEF_DAY_HOOK_ESTIMATOR")
    if strength:
        keys.append("DEF_DAY_LODO_STRENGTH")
    if cd4:
        keys.append("DEF_DAY_CD4_STACK")
    return keys


def _between(text: str, start: str, end: str) -> str:
    i = text.find(start)
    if i < 0:
        raise AssertionError(f"missing {start!r}")
    j = text.find(end, i + len(start))
    if j < 0:
        raise AssertionError(f"missing end after {start!r}")
    return text[i + len("python3 -c '"):j] if start.startswith("python3 -c '") else text[i + len(start):j]


def _name(rec) -> str:
    argv = rec["argv"]
    if "--name" in argv:
        return argv[argv.index("--name") + 1]
    return ""


def _sub(s: str, work: str, bundle: str) -> str:
    return s.replace(work, "$W").replace(bundle, "$B")


def _norm(records, work: str, bundle: str):
    out = []
    for rec in records:
        out.append({"argv": [_sub(a, work, bundle) for a in rec["argv"]],
                    "vcc": {k: _sub(v, work, bundle) for k, v in sorted(rec["vcc"].items())}})
    return out


if __name__ == "__main__":
    unittest.main(verbosity=2)
