"""428: zero_repair invariants on synthetic blocks + VCC_DE_ZERO_REPAIR parse/refusals. TEST_BUILDER / TEST_MODULE test
dry copies."""
import importlib.util
import os
import subprocess
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
MOD = Path(os.environ.get("TEST_MODULE", ROOT / "src/vcc2026/zero_repair.py"))
spec = importlib.util.spec_from_file_location("zr", MOD); zr = importlib.util.module_from_spec(spec); spec.loader.exec_module(zr)
rng = np.random.default_rng(0)
G, n, K = 60, 400, 3
gated = np.zeros(G, bool); gated[:40] = True
ctrl = rng.uniform(6, 500, 40)
blk = rng.poisson(3.0, (n, G)).astype(np.int64)
blk[:, 5] = 0; blk[:, 7] = 0; blk[:, 50] = 0          # two gated zeros (5, 7) and one sub-gate zero (50)
blk[:, 9] = 0                                         # the target's own column: never repaired
out, d = zr.zero_repair(blk, n - K, gated, ctrl, 9, np.random.default_rng(1), floor_log2=-5.0)
assert d["zero_repaired"] == 2, d
assert out[:, 5].sum() > 0 and out[:, 7].sum() > 0 and out[:, 9].sum() == 0
assert np.array_equal(out.sum(1), blk.sum(1)), "libraries moved"
assert np.array_equal(out[n - K:], blk[n - K:]), "ballast rows moved"
others = np.setdiff1d(np.flatnonzero(gated), [5, 7])
assert np.array_equal(out[:, others], blk[:, others]), "another tested gene moved"
assert out.min() >= 0
ok_blk = blk.copy(); ok_blk[:, 5] = 1; ok_blk[:, 7] = 1
same, d2 = zr.zero_repair(ok_blk, n - K, gated, ctrl, 9, np.random.default_rng(1))
assert same is ok_blk and d2["zero_repaired"] == 0, "no-op must return the input object"
print("zero_repair invariants PASS")

BUILDER = Path(os.environ.get("TEST_BUILDER", ROOT / "scripts/make_submission_pm_cd4.py"))
BASE = dict(VCC_DE_SIGN_PRIOR="clean4q80", VCC_PERT_SCALE="wshare_lo0.3_hi3", VCC_PERT_SCALE_BALLAST="1")
CODE = (f"import importlib.util as u, sys; s=u.spec_from_file_location('b', {str(BUILDER)!r}); m=u.module_from_spec(s); "
        "sys.modules['b']=m; s.loader.exec_module(m); print('ZR', m.DE_ZERO_REPAIR)")


def run(env):
    e = {k: v for k, v in os.environ.items() if not k.startswith("VCC_")}
    e.update(env)
    return subprocess.run([sys.executable, "-c", CODE], env=e, capture_output=True, text=True, cwd=str(ROOT))


r = run(dict(BASE)); assert r.returncode == 0 and "ZR None" in r.stdout, r.stderr[-400:]
r = run(dict(BASE, VCC_DE_ZERO_REPAIR="-5")); assert r.returncode == 0 and "ZR -5.0" in r.stdout, r.stderr[-400:]
for bad, want in ((dict(VCC_DE_SIGN_PRIOR="clean4q80", VCC_PERT_SCALE="wshare_lo0.3_hi3", VCC_DE_ZERO_REPAIR="-5"), "needs VCC_PERT_SCALE_BALLAST"),
                  (dict(BASE, VCC_DE_ZERO_REPAIR="-0.5"), "must be in [-10, -1]"),
                  (dict(BASE, VCC_DE_ZERO_REPAIR="x"), "must be a float")):
    r = run(bad)
    assert r.returncode != 0 and want in (r.stderr + r.stdout), (bad, (r.stderr + r.stdout)[-400:])
print("428 VCC_DE_ZERO_REPAIR tests PASS (unset None; -5 parses; 3 refusals)")
