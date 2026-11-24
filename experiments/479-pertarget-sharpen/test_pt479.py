"""479 unit tests on a builder file (default: the production builder): parse of VCC_BULK_SHARPEN_Q's per-target forms, and
sharpen_q_by_tercile on 446bj2's LODO stack dump (third session's known answer: key terciles 100/100/100 with the 97 floor
targets in LOW; strength terciles 100/100/100).
usage: uv run python experiments/479-pertarget-sharpen/test_pt479.py [builder.py]"""
import importlib.util
import os
import subprocess
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
BUILDER = Path(sys.argv[1]) if len(sys.argv) > 1 else ROOT / "scripts/make_submission_pm_cd4.py"
DUMP = ROOT / "experiments/417-bulk-sharpen/work/dump_h1x2_cd4x15_lodo05.npz"
BASE = {"VCC_PERT_SCALE": "wshare_lo0.3_hi3", "VCC_PERT_SCALE_BALLAST": "1"}


def load(env):
    for k in [k for k in os.environ if k.startswith("VCC_")]:
        del os.environ[k]
    os.environ.update(env)
    name = f"b479_{abs(hash(tuple(sorted(env.items()))))}"
    spec = importlib.util.spec_from_file_location(name, BUILDER)
    m = importlib.util.module_from_spec(spec); sys.modules[name] = m; spec.loader.exec_module(m)
    return m


def refused(env):
    code = ("import importlib.util,sys; s=importlib.util.spec_from_file_location('x', sys.argv[1]); "
            "m=importlib.util.module_from_spec(s); sys.modules['x']=m; s.loader.exec_module(m)")
    e = {k: v for k, v in os.environ.items() if not k.startswith("VCC_")}; e.update(env)
    r = subprocess.run([sys.executable, "-c", code, str(BUILDER)], env=e, capture_output=True, text=True)
    return r.returncode != 0, (r.stderr.strip().splitlines() or [""])[-1][:150]


ok = True
m0 = load({})
assert m0.BULK_SHARPEN_PT is None and m0.BULK_SHARPEN_Q is None
m1 = load({**BASE, "VCC_BULK_SHARPEN_Q": "1.75"})
assert m1.BULK_SHARPEN_PT is None and m1.BULK_SHARPEN_Q == 1.75
m2 = load({**BASE, "VCC_BULK_SHARPEN_Q": "tercile:2.5,2.0,1.75"})
assert m2.BULK_SHARPEN_PT == {"mode": "tercile", "q": [2.5, 2.0, 1.75]} and m2.BULK_SHARPEN_Q is None
m3 = load({**BASE, "VCC_BULK_SHARPEN_Q": "stercile:2.5,2,1.75"})
assert m3.BULK_SHARPEN_PT == {"mode": "stercile", "q": [2.5, 2.0, 1.75]}
print("PASS parse: unset, scalar 1.75, tercile, stercile")
for env, why in (({**BASE, "VCC_BULK_SHARPEN_Q": "tercile:2.5,2.0"}, "two values"),
                 ({**BASE, "VCC_BULK_SHARPEN_Q": "tercile:1.0,2,2"}, "q <= 1"),
                 ({**BASE, "VCC_BULK_SHARPEN_Q": "stercile:4,2,2"}, "q > 3.5"),
                 ({**BASE, "VCC_BULK_SHARPEN_Q": "tercile:a,b,c"}, "not numbers"),
                 ({**BASE, "VCC_BULK_SHARPEN_Q": "2.5"}, "scalar > 2 (unchanged rule)"),
                 ({"VCC_PERT_SCALE": "wshare_lo0.3_hi3", "VCC_BULK_SHARPEN_Q": "tercile:2.5,2,1.75"}, "no ballast")):
    r, msg = refused(env); ok &= r
    print(("PASS" if r else "FAIL") + f" refused: {why} ({msg})")
z = np.load(DUMP, allow_pickle=False)
genes, perts = list(map(str, z["genes"])), list(map(str, z["perts"]))
for mode, want in (("tercile", ([100, 100, 100], 97)), ("stercile", ([100, 100, 100], 0))):
    q, a = m2.sharpen_q_by_tercile(z["bulk"], z["carriers"], genes, perts, {"mode": mode, "q": [2.5, 2.0, 1.75]})
    good = a["counts"] == want[0] and a["floor_to_low"] == want[1] and a["tie_to_mid"] == 0 and a["nonfinite"] == 0
    good &= set(np.unique(q)) == {2.5, 2.0, 1.75} and int((q == 2.5).sum()) == 100
    ok &= good
    print(("PASS" if good else "FAIL") + f" {mode} on 446bj2's LODO dump: {a}")
q_all, a_all = m2.sharpen_q_by_tercile(z["bulk"], z["carriers"], genes, perts, {"mode": "tercile", "q": [1.75, 1.75, 1.75]})
good = bool(np.all(q_all == 1.75)); ok &= good
print(("PASS" if good else "FAIL") + " all-1.75 gives q == 1.75 for every target")
print("ALL PASS" if ok else "SOME FAILED"); sys.exit(0 if ok else 1)
