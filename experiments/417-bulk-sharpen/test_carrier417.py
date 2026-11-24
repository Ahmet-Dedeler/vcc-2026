"""417 window 2c: VCC_CARRIER_GAMMA parses and refuses (the second pass itself is verified by the KA build: gamma 0.65 must
reproduce the shipped carriers, and so 417a's X, bit for bit). TEST_BUILDER=<path> tests a dry copy."""
import os
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
BUILDER = Path(os.environ.get("TEST_BUILDER", ROOT / "scripts/make_submission_pm_cd4.py"))
BASE = dict(VCC_PERT_SCALE="wshare_lo0.3_hi3", VCC_PERT_SCALE_BALLAST="1")
CODE = (f"import importlib.util as u, sys; s=u.spec_from_file_location('b', {str(BUILDER)!r}); m=u.module_from_spec(s); "
        "sys.modules['b']=m; s.loader.exec_module(m); print('CG', m.CARRIER_GAMMA)")


def run(env):
    e = {k: v for k, v in os.environ.items() if not k.startswith("VCC_")}
    e.update(env)
    return subprocess.run([sys.executable, "-c", CODE], env=e, capture_output=True, text=True)


r = run(dict(BASE)); assert r.returncode == 0 and "CG None" in r.stdout, r.stderr[-300:]
r = run(dict(BASE, VCC_CARRIER_GAMMA="0")); assert r.returncode == 0 and "CG 0.0" in r.stdout, r.stderr[-300:]
for bad, want in ((dict(VCC_PERT_SCALE="wshare_lo0.3_hi3", VCC_CARRIER_GAMMA="0"), "needs VCC_PERT_SCALE_BALLAST"),
                  (dict(BASE, VCC_CARRIER_GAMMA="4"), "must be in [0, 3]")):
    r = run(bad)
    assert r.returncode != 0 and want in (r.stderr + r.stdout), (bad, r.returncode, r.stderr[-300:])
print("417 carrier-gamma tests PASS (unset None; 0 parses; 2 refusals)")
