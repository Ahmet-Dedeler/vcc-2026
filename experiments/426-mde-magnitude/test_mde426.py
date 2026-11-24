"""426: VCC_DE_MDE parses and refuses (the emission itself is verified by the H1 known answers: run426.py --emitter prod
must reproduce 416's base bit-for-bit with the knob off and k_emit426's MDE run bit-for-bit with it on).
TEST_BUILDER=<path> tests a dry copy."""
import os
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
BUILDER = Path(os.environ.get("TEST_BUILDER", ROOT / "scripts/make_submission_pm_cd4.py"))
BASE = dict(VCC_DE_SIGN_PRIOR="clean4q80", VCC_PERT_SCALE="wshare_lo0.3_hi3", VCC_PERT_SCALE_BALLAST="1")
CODE = (f"import importlib.util as u, sys; s=u.spec_from_file_location('b', {str(BUILDER)!r}); m=u.module_from_spec(s); "
        "sys.modules['b']=m; s.loader.exec_module(m); print('MDE', m.DE_MDE)")


def run(env):
    e = {k: v for k, v in os.environ.items() if not k.startswith("VCC_")}
    e.update(env)
    return subprocess.run([sys.executable, "-c", CODE], env=e, capture_output=True, text=True, cwd=str(ROOT))


r = run(dict(BASE)); assert r.returncode == 0 and "MDE None" in r.stdout, r.stderr[-400:]
r = run(dict(BASE, VCC_DE_MDE="0.1,0.1,0.3")); assert r.returncode == 0 and "MDE (0.1, 0.1, 0.3)" in r.stdout, r.stderr[-400:]
for bad, want in ((dict(VCC_DE_MDE="0.1,0.1,0.3"), "needs VCC_DE_SIGN_PRIOR"),
                  (dict(BASE, VCC_DE_MDE="0.1,0.1"), "must be 'c,cap,min_count'"),
                  (dict(BASE, VCC_DE_MDE="0.9,0.1,0.3"), "must be 'c,cap,min_count'"),
                  (dict(BASE, VCC_DE_MDE="a,b,c"), "floats")):
    r = run(bad)
    assert r.returncode != 0 and want in (r.stderr + r.stdout), (bad, r.returncode, (r.stderr + r.stdout)[-400:])
print("426 VCC_DE_MDE tests PASS (unset None; 0.1,0.1,0.3 parses; 4 refusals)")
