#!/usr/bin/env python3
"""528 unit tests for patch_builder528.py on builder COPIES (never the production file).

usage: PYTHONDONTWRITEBYTECODE=1 uv run python experiments/528-callprop-block/test_patch528.py <builder-copy.py>
It copies <builder-copy.py> to work/patchtest/bld_base.py and bld_528.py, patches the latter, then checks:
  P1  unset and `cpm` parse to None and leave every other call-set constant equal to the base's
  P2  prop:<npz>:<sha256> -> swap 1.0, block_keep 0; :0.25 -> block_keep 1650; with 523's VCC_RC_CALLSET_CPM=1600, 0.25 -> 1200
  P2b steer kept (default) / key / cpm parse (case-insensitive)
  P3  refusals: bad source, short sha, missing path, swap 0 / 1.5 / 'x', a bad steer, 6 fields
  P4  the patched builder differs from the base in exactly the six intended places (difflib hunks), compiles, and pins the
      patched emitter's sha (k_emit528.py) in DE_SIGN_EMITTER_SHA256
  P5  the patched emitter (k_emit528.py) and the production emitter differ only by patch_emitter528's replacements
"""
import difflib
import hashlib
import importlib.util
import os
import py_compile
import shutil
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
SRC = Path(sys.argv[1])
W = HERE / "work/patchtest"
W.mkdir(parents=True, exist_ok=True)
B0, B1 = W / "bld_base.py", W / "bld_528.py"
_txt = SRC.read_text()
_root = "ROOT = Path(__file__).resolve().parents[1]\n"
assert _txt.count(_root) == 1
_txt = _txt.replace(_root, f'ROOT = Path("{ROOT}")  # 528 test copy only (the builder resolves its inputs from ROOT)\n')
B0.write_text(_txt)
B1.write_text(_txt)
subprocess.run([sys.executable, str(HERE / "patch_builder528.py"), str(B1)], check=True)
ART = HERE / "work/callblock528_abc.npz"
SHA = hashlib.sha256(ART.read_bytes()).hexdigest()
EM = hashlib.sha256((HERE / "k_emit528.py").read_bytes()).hexdigest()


def load(path, env):
    for k in [k for k in os.environ if k.startswith("VCC_")]:
        del os.environ[k]
    os.environ.update(env)
    name = f"b528_{path.stem}_{abs(hash(tuple(sorted(env.items()))))}"
    spec = importlib.util.spec_from_file_location(name, path)
    m = importlib.util.module_from_spec(spec)
    sys.modules[name] = m
    spec.loader.exec_module(m)
    return m


def refused(path, env):
    code = ("import importlib.util,sys; s=importlib.util.spec_from_file_location('x', sys.argv[1]); "
            "m=importlib.util.module_from_spec(s); sys.modules['x']=m; s.loader.exec_module(m)")
    e = {k: v for k, v in os.environ.items() if not k.startswith("VCC_")}
    e.update(env)
    r = subprocess.run([sys.executable, "-c", code, str(path)], env=e, capture_output=True, text=True)
    return r.returncode != 0, (r.stderr.strip().splitlines() or [""])[-1][:160]


ok = {}
m0, m1 = load(B0, {}), load(B1, {})
ok["P1 unset -> None"] = m1.RC_CALLSET_BLOCK is None and m1.RC_CALLSET_BLOCK_AUDIT == {}
ok["P1 unset: call-set constants equal"] = all(getattr(m0, k) == getattr(m1, k) for k in
                                               ("RC_CALLSET_HEAD", "RC_CALLSET_CPM", "RC_CALLSET_Z_CAP", "RC_SCORER_MIN_CPM"))
ok["P1 cpm -> None"] = load(B1, {"VCC_RC_CALLSET_BLOCK": "cpm"}).RC_CALLSET_BLOCK is None
m = load(B1, {"VCC_RC_CALLSET_BLOCK": f"prop:{ART}:{SHA}"})
ok["P2 prop full"] = m.RC_CALLSET_BLOCK == {"source": "prop", "path": str(ART), "sha256": SHA, "swap": 1.0, "block_keep": 0,
                                           "steer": "kept"}
for st in ("kept", "key", "cpm", "KEY"):
    m = load(B1, {"VCC_RC_CALLSET_BLOCK": f"prop:{ART}:{SHA}:1:{st}"})
    ok[f"P2 steer {st}"] = m.RC_CALLSET_BLOCK["steer"] == st.lower() and m.RC_CALLSET_BLOCK["block_keep"] == 0
m = load(B1, {"VCC_RC_CALLSET_BLOCK": f"prop:{ART}:{SHA}:0.25"})
ok["P2 prop 0.25 -> keep 1650"] = m.RC_CALLSET_BLOCK["block_keep"] == 1650 and m.RC_CALLSET_BLOCK["swap"] == 0.25
m = load(B1, {"VCC_RC_CALLSET_BLOCK": f"prop:{ART}:{SHA}:0.25", "VCC_RC_CALLSET_CPM": "1600"})
ok["P2 with VCC_RC_CALLSET_CPM=1600, 0.25 -> keep 1200"] = m.RC_CALLSET_BLOCK["block_keep"] == 1200
for label, val in (("bad source", f"cpmx:{ART}:{SHA}"), ("short sha", f"prop:{ART}:{SHA[:40]}"), ("missing path", f"prop::{SHA}"),
                   ("swap 0", f"prop:{ART}:{SHA}:0"), ("swap 1.5", f"prop:{ART}:{SHA}:1.5"), ("swap x", f"prop:{ART}:{SHA}:x"),
                   ("bad steer", f"prop:{ART}:{SHA}:0.5:both"), ("6 fields", f"prop:{ART}:{SHA}:0.5:kept:x")):
    r, msg = refused(B1, {"VCC_RC_CALLSET_BLOCK": val})
    ok[f"P3 refuses {label}"] = r and "VCC_RC_CALLSET_BLOCK" in msg
a, b = B0.read_text().splitlines(), B1.read_text().splitlines()
hunks = [g for g in difflib.SequenceMatcher(None, a, b).get_opcodes() if g[0] != "equal"]
ok["P4 six hunks"] = len(hunks) == 6
py_compile.compile(str(B1), cfile=str(W / "bld_528.pyc"), doraise=True)
ok["P4 compiles"] = True
ok["P4 emitter pin = k_emit528"] = f'DE_SIGN_EMITTER_SHA256 = "{EM}"' in B1.read_text()
pe = (ROOT / "src/vcc2026/rc_emit_k.py").read_text()
sys.path.insert(0, str(HERE))
from patch_emitter528 import patch_text  # noqa: E402
ok["P5 k_emit528 == patch_emitter528(production)"] = patch_text(pe) == (HERE / "k_emit528.py").read_text()
for k, v in ok.items():
    print(("PASS " if v else "FAIL ") + k)
print(f"{sum(map(bool, ok.values()))}/{len(ok)} PASS  (builder copy sha {hashlib.sha256(SRC.read_bytes()).hexdigest()[:16]}; "
      f"patched {hashlib.sha256(B1.read_bytes()).hexdigest()[:16]}; emitter pin {EM[:16]})")
sys.exit(0 if all(ok.values()) else 1)
