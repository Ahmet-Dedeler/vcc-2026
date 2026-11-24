"""422: VCC_CARRIER_SIGN_ONLY parses and refuses; the sign-only arithmetic keeps |final| and takes the pass's signs
(shipped sign where the pass is 0). TEST_BUILDER=<path> tests a dry copy."""
import os
import subprocess
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
BUILDER = Path(os.environ.get("TEST_BUILDER", ROOT / "scripts/make_submission_pm_cd4.py"))
BASE = dict(VCC_PERT_SCALE="wshare_lo0.3_hi3", VCC_PERT_SCALE_BALLAST="1")
CODE = (f"import importlib.util as u, sys; s=u.spec_from_file_location('b', {str(BUILDER)!r}); m=u.module_from_spec(s); "
        "sys.modules['b']=m; s.loader.exec_module(m); print('SO', m.CARRIER_SIGN_ONLY)")


def run(env):
    e = {k: v for k, v in os.environ.items() if not k.startswith("VCC_")}
    e.update(env)
    return subprocess.run([sys.executable, "-c", CODE], env=e, capture_output=True, text=True)


r = run(dict(BASE)); assert r.returncode == 0 and "SO None" in r.stdout, r.stderr[-300:]
r = run(dict(BASE, VCC_CARRIER_GAMMA="0", VCC_CARRIER_SIGN_ONLY="1")); assert r.returncode == 0 and "SO 1" in r.stdout, r.stderr[-300:]
for bad, want in ((dict(BASE, VCC_CARRIER_SIGN_ONLY="1"), "needs VCC_CARRIER_GAMMA"),
                  (dict(BASE, VCC_CARRIER_GAMMA="0", VCC_CARRIER_SIGN_ONLY="yes"), "must be 1 when set")):
    r = run(bad)
    assert r.returncode != 0 and want in (r.stderr + r.stdout), (bad, r.returncode, r.stderr[-300:])
# the arithmetic, as patched
rng = np.random.default_rng(422)
final = rng.normal(0, 1, (5, 40)); final[0, 3] = 0.0
fc = rng.normal(0, 1, (5, 40)); fc[1, 7] = 0.0
cf = np.abs(final) * np.where(fc != 0, np.sign(fc), np.sign(final))
assert np.allclose(np.abs(cf), np.abs(final)) and np.allclose(np.linalg.norm(cf, axis=1), np.linalg.norm(final, axis=1))
assert np.sign(cf[1, 7]) == np.sign(final[1, 7]) and cf[0, 3] == 0.0
src = BUILDER.read_text()
assert "carrier_final = np.abs(final) * np.where(fc != 0, np.sign(fc), np.sign(final))" in src
print("422 carrier-sign-only tests PASS (unset None; parse; 2 refusals; arithmetic: |final| kept, pass signs, zero fallback)")
