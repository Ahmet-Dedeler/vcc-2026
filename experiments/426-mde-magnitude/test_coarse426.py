"""426h: VCC_DE_MDE's optional 4th field `coarse` parses and refuses. The emission itself is verified by the H1 known
answers (run426.py with RUN426_EMITTER=path on the patched emitter reproduces k_emit426 bit-for-bit: knob off vs 416's
base, plain MDE, coarse MDE). TEST_BUILDER=<path> tests a dry copy."""
import os
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
BUILDER = Path(os.environ.get("TEST_BUILDER", ROOT / "scripts/make_submission_pm_cd4.py"))
BASE = dict(VCC_DE_SIGN_PRIOR="clean4q80", VCC_PERT_SCALE="wshare_lo0.3_hi3", VCC_PERT_SCALE_BALLAST="1")
CODE = (f"import importlib.util as u, sys; s=u.spec_from_file_location('b', {str(BUILDER)!r}); m=u.module_from_spec(s); "
        "sys.modules['b']=m; s.loader.exec_module(m); print('MDE', m.DE_MDE, 'COARSE', m.DE_MDE_COARSE)")


def run(env):
    e = {k: v for k, v in os.environ.items() if not k.startswith("VCC_")}
    e.update(env)
    return subprocess.run([sys.executable, "-c", CODE], env=e, capture_output=True, text=True, cwd=str(ROOT))


r = run(dict(BASE)); assert r.returncode == 0 and "MDE None COARSE False" in r.stdout, r.stderr[-400:]
r = run(dict(BASE, VCC_DE_MDE="0.1,0.1,0.3"))
assert r.returncode == 0 and "MDE (0.1, 0.1, 0.3) COARSE False" in r.stdout, r.stderr[-400:]
r = run(dict(BASE, VCC_DE_MDE="0.1,0.1,0.3,coarse"))
assert r.returncode == 0 and "MDE (0.1, 0.1, 0.3) COARSE True" in r.stdout, r.stderr[-400:]
for bad, want in ((dict(BASE, VCC_DE_MDE="0.1,0.1,0.3,fine"), "or 'c,cap,min_count,coarse'"),
                  (dict(BASE, VCC_DE_MDE="0.1,0.1,0.3,coarse,x"), "or 'c,cap,min_count,coarse'"),
                  (dict(BASE, VCC_DE_MDE="0.1,0.1,0.3,0.5"), "or 'c,cap,min_count,coarse'"),
                  (dict(BASE, VCC_DE_MDE="coarse"), "or 'c,cap,min_count,coarse'"),
                  (dict(BASE, VCC_DE_MDE="0.9,0.1,0.3,coarse"), "0 < c <= 0.5"),
                  (dict(VCC_DE_MDE="0.1,0.1,0.3,coarse"), "needs VCC_DE_SIGN_PRIOR")):
    r = run(bad)
    assert r.returncode != 0 and want in (r.stderr + r.stdout), (bad, r.returncode, (r.stderr + r.stdout)[-400:])
print("426h VCC_DE_MDE coarse tests PASS (unset; 3-field; 4-field coarse; 6 refusals)")
