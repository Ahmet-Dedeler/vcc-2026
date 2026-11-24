"""417 window 2b: VCC_BULK_DUMP / VCC_BULK_LOAD parse and refuse correctly (the I/O itself is verified in the first real
use: a LOAD of a same-env DUMP must log max|loaded-own| = 0). TEST_BUILDER=<path> tests a dry copy."""
import os
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
BUILDER = Path(os.environ.get("TEST_BUILDER", ROOT / "scripts/make_submission_pm_cd4.py"))
SHA = "0" * 64
BASE = dict(VCC_PERT_SCALE="wshare_lo0.3_hi3", VCC_PERT_SCALE_BALLAST="1")
CODE = (f"import importlib.util as u, sys; s=u.spec_from_file_location('b', {str(BUILDER)!r}); m=u.module_from_spec(s); "
        "sys.modules['b']=m; s.loader.exec_module(m); print('DUMP', m.BULK_DUMP); print('LOAD', m.BULK_LOAD)")


def run(env):
    e = {k: v for k, v in os.environ.items() if not k.startswith("VCC_")}
    e.update(env)
    return subprocess.run([sys.executable, "-c", CODE], env=e, capture_output=True, text=True)


r = run(dict(BASE))
assert r.returncode == 0 and "DUMP None" in r.stdout and "LOAD None" in r.stdout, r.stderr[-400:]
r = run(dict(BASE, VCC_BULK_LOAD=f"/x/y.npz:{SHA}"))
assert r.returncode == 0 and f"LOAD (PosixPath('/x/y.npz'), '{SHA}')" in r.stdout, (r.stdout, r.stderr[-400:])
r = run(dict(BASE, VCC_BULK_DUMP="/x/d.npz"))
assert r.returncode == 0 and "DUMP /x/d.npz" in r.stdout, r.stderr[-400:]
for bad, want in ((dict(VCC_PERT_SCALE="wshare_lo0.3_hi3", VCC_BULK_LOAD=f"/x/y.npz:{SHA}"), "needs VCC_PERT_SCALE_BALLAST"),
                  (dict(BASE, VCC_BULK_LOAD="/x/y.npz"), "must be '<path.npz>:<sha256>'"),
                  (dict(BASE, VCC_BULK_LOAD="/x/y.npz:abc"), "full sha256"),
                  (dict(BASE, VCC_BULK_LOAD=f"/x/y.npz:{SHA}", VCC_BULK_DUMP="/x/d.npz"), "exclusive")):
    r = run(bad)
    assert r.returncode != 0 and want in (r.stderr + r.stdout), (bad, r.returncode, r.stderr[-300:])
print("417 bulk-io tests PASS (unset None; LOAD/DUMP parse; 4 refusals)")
