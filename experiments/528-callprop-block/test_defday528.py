#!/usr/bin/env python3
"""528 tests for patch_defday528.py on a def_day.sh COPY (never the production script).

usage: PYTHONDONTWRITEBYTECODE=1 uv run python experiments/528-callprop-block/test_defday528.py <def_day.sh copy>
  D1  the patch applies; bash -n passes on the patched copy; the diff is exactly the seven intended insertions
  D2  the validation block (run in isolation with stub log/die): unset passes silently; `prop` dies on the UNPATCHED
      production builder (no knob), and in a fake root holding the patched builder / emitter / helper prop, prop:0.25 and
      prop:1 pass while prop:0 / prop:1.5 / cpm / prop:x die with the knob's name
  D3  the compute block on the real challenge bundle (A/B/C as a stand-in for D/E/F): unset -> CSB empty; prop -> the
      artifact equals callblock528's A/B/C artifact byte for byte (sha 04fb24ec...) and CSB = prop:<W>/callblock_def.npz:<sha>;
      prop:0.25 -> the same with :0.25
  D4  the final-meta check accepts a meta carrying the knob and refuses a meta without it or with another sha/swap
"""
import difflib
import hashlib
import json
import os
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
SRC = Path(sys.argv[1])
W = HERE / "work/patchtest/defday"
W.mkdir(parents=True, exist_ok=True)
D0, D1 = W / "def_day_base.sh", W / "def_day_528.sh"
D0.write_text(SRC.read_text())
D1.write_text(SRC.read_text())
subprocess.run([sys.executable, str(HERE / "patch_defday528.py"), str(D1)], check=True)
ok = {}
ok["D1 bash -n"] = subprocess.run(["/bin/bash", "-n", str(D1)]).returncode == 0
a, b = D0.read_text().splitlines(), D1.read_text().splitlines()
ops = [g for g in difflib.SequenceMatcher(None, a, b).get_opcodes() if g[0] != "equal"]
ok["D1 seven insertions only"] = len(ops) == 7 and all(g[0] == "insert" for g in ops)
txt = D1.read_text()


def block(start, end):
    i = txt.index(start)
    j = txt.index(end, i)
    return txt[i:j]


VAL = block("# 528: the call-propensity CPM block (DE-only", "\nSQ=")
COMP = block("# 2c. 528: DEF_DAY_CALLSET_BLOCK", "# 2b. DEF_DAY_LODO=1")
FIN = block('if [ -n "${VCC_RC_CALLSET_BLOCK:-}" ]; then\n  python3 -', 'for P in "${POINTS[@]}"; do [ "${P%%.ntclevel.json:*}" != "$FINAL" ]')
STUB = 'log(){ echo "LOG $*"; }; die(){ echo "DIE $*"; exit 3; }; LOG=/dev/null\n'


def sh(code, env):
    e = {k: v for k, v in os.environ.items() if not k.startswith(("DEF_DAY_", "VCC_"))}
    e.update(env)
    r = subprocess.run(["/bin/bash", "-c", STUB + code], env=e, cwd=ROOT, capture_output=True, text=True)
    return r.returncode, r.stdout + r.stderr


# the validation greps the builder / emitter / helper by relative path: run it in a fake root holding the PATCHED copies
FR = W / "fakeroot"
for rel, src in (("scripts/make_submission_pm_cd4.py", HERE / "work/patchtest/bld_528.py"),
                 ("src/vcc2026/rc_emit_k.py", HERE / "k_emit528.py"),
                 ("experiments/528-callprop-block/callblock528.py", HERE / "callblock528.py")):
    (FR / rel).parent.mkdir(parents=True, exist_ok=True)
    (FR / rel).write_bytes(src.read_bytes())


def sh_in(code, env, cwd):
    e = {k: v for k, v in os.environ.items() if not k.startswith(("DEF_DAY_", "VCC_"))}
    e.update(env)
    r = subprocess.run(["/bin/bash", "-c", STUB + code], env=e, cwd=cwd, capture_output=True, text=True)
    return r.returncode, r.stdout + r.stderr


rc, out = sh(VAL, {})
ok["D2 unset passes silently (production tree)"] = rc == 0 and "DIE" not in out
rc, out = sh(VAL, {"DEF_DAY_CALLSET_BLOCK": "prop"})
ok["D2 prop on the UNPATCHED production builder dies (builder has no knob)"] = rc == 3 and "VCC_RC_CALLSET_BLOCK (528)" in out
for v in ("prop", "prop:0.25", "prop:1", "prop:1:kept", "prop:0.5:key", "prop:1:cpm"):
    rc, out = sh_in(VAL, {"DEF_DAY_CALLSET_BLOCK": v}, FR)
    ok[f"D2 {v} passes (patched tree)"] = rc == 0 and "DIE" not in out
for v in ("prop:0", "prop:1.5", "cpm", "prop:x", "prop:1:both", "prop:1:kept:x"):
    rc, out = sh_in(VAL, {"DEF_DAY_CALLSET_BLOCK": v}, FR)
    ok[f"D2 {v} dies"] = rc == 3 and "DEF_DAY_CALLSET_BLOCK" in out
art = HERE / "work/callblock528_abc.npz"
want_sha = hashlib.sha256(art.read_bytes()).hexdigest()
TW = W / "daywork"
TW.mkdir(exist_ok=True)
env = {"B": str(ROOT / "data/challenge"), "W": str(TW), "CTXS": "A,B,C", "PYTHONDONTWRITEBYTECODE": "1"}
rc, out = sh(COMP + 'echo "CSB=$CSB"\n', env)
ok["D3 unset -> empty"] = rc == 0 and "CSB=\n" in out + "\n"
rc, out = sh(COMP + 'echo "CSB=$CSB"\n', {**env, "DEF_DAY_CALLSET_BLOCK": "prop"})
got_sha = hashlib.sha256((TW / "callblock_def.npz").read_bytes()).hexdigest() if (TW / "callblock_def.npz").exists() else ""
ok["D3 prop: the artifact == callblock528's A/B/C artifact"] = got_sha == want_sha
ok["D3 prop: CSB"] = f"CSB=prop:{TW}/callblock_def.npz:{want_sha}\n" in out
rc, out = sh(COMP + 'echo "CSB=$CSB"\n', {**env, "DEF_DAY_CALLSET_BLOCK": "prop:0.25"})
ok["D3 prop:0.25: CSB"] = f"CSB=prop:{TW}/callblock_def.npz:{want_sha}:0.25\n" in out
rc, out = sh(COMP + 'echo "CSB=$CSB"\n', {**env, "DEF_DAY_CALLSET_BLOCK": "prop:1:key"})
ok["D3 prop:1:key: CSB"] = f"CSB=prop:{TW}/callblock_def.npz:{want_sha}:1:key\n" in out
meta = TW / "final.h5ad.meta.json"
for label, cfg, knob, expect in (
        ("accepts the matching meta", {"rc_callset_block": {"source": "prop", "sha256": want_sha, "swap": 1.0, "steer": "kept"}},
         f"prop:{TW}/callblock_def.npz:{want_sha}", 0),
        ("accepts swap 0.25", {"rc_callset_block": {"source": "prop", "sha256": want_sha, "swap": 0.25, "steer": "kept"}},
         f"prop:{TW}/callblock_def.npz:{want_sha}:0.25", 0),
        ("accepts steer key", {"rc_callset_block": {"source": "prop", "sha256": want_sha, "swap": 1.0, "steer": "key"}},
         f"prop:{TW}/callblock_def.npz:{want_sha}:1:key", 0),
        ("refuses a meta without it", {}, f"prop:{TW}/callblock_def.npz:{want_sha}", 3),
        ("refuses another sha", {"rc_callset_block": {"source": "prop", "sha256": "0" * 64, "swap": 1.0, "steer": "kept"}},
         f"prop:{TW}/callblock_def.npz:{want_sha}", 3),
        ("refuses another swap", {"rc_callset_block": {"source": "prop", "sha256": want_sha, "swap": 1.0, "steer": "kept"}},
         f"prop:{TW}/callblock_def.npz:{want_sha}:0.5", 3),
        ("refuses another steer", {"rc_callset_block": {"source": "prop", "sha256": want_sha, "swap": 1.0, "steer": "key"}},
         f"prop:{TW}/callblock_def.npz:{want_sha}", 3)):
    meta.write_text(json.dumps({"config": cfg}))
    rc, out = sh(FIN, {"FINAL": str(TW / "final"), "VCC_RC_CALLSET_BLOCK": knob})
    ok[f"D4 {label}"] = rc == expect
rc, out = sh(FIN, {"FINAL": str(TW / "final")})
ok["D4 unset: no check"] = rc == 0
for k, v in ok.items():
    print(("PASS " if v else "FAIL ") + k)
print(f"{sum(map(bool, ok.values()))}/{len(ok)} PASS (def_day copy sha {hashlib.sha256(SRC.read_bytes()).hexdigest()[:16]})")
sys.exit(0 if all(ok.values()) else 1)
