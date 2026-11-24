"""417: bulk sharpening restores each row's gated non-panel L2 norm, leaves panel/ungated coordinates alone, equals 169's
sharpen(mode='row') on the panel-zeroed row 414 priced, renders through profile_from_effect, and the knob refuses to run
without ballast. TEST_BUILDER=<path> tests a dry-run copy instead of scripts/make_submission_pm_cd4.py."""
import importlib.util
import os
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
BUILDER = Path(os.environ.get("TEST_BUILDER", ROOT / "scripts/make_submission_pm_cd4.py"))
sys.path.insert(0, str(ROOT / "scripts"))
sys.path.insert(0, str(ROOT / "experiments/169-sharpening-six-member"))
from sharp import sharpen  # noqa: E402

env = dict(VCC_PERT_SCALE="wshare_lo0.3_hi3", VCC_PERT_SCALE_BALLAST="1", VCC_BULK_SHARPEN_Q="1.25")
os.environ.update(env)
spec_ = importlib.util.spec_from_file_location("builder417", BUILDER)
M = importlib.util.module_from_spec(spec_)
sys.modules["builder417"] = M
spec_.loader.exec_module(M)
assert M.BULK_SHARPEN_Q == 1.25, M.BULK_SHARPEN_Q

perts = [str(t) for t in pd.read_csv(ROOT / "data/challenge/pert_counts.csv")["target_gene"]]
genes = [str(g) for g in pd.read_csv(ROOT / "data/challenge/gene_names.csv")["gene_name"]]
panel = np.array([g in set(perts) for g in genes])
assert panel.sum() == 300, panel.sum()
rng = np.random.default_rng(417)
b_ctrl = np.log1p(rng.gamma(0.3, 2.0, len(genes)))            # a control log-profile with a realistic gated fraction
expr = np.expm1(b_ctrl)
gate = (1e6 * expr / expr.sum() > 10.0) & (b_ctrl > 0)
worst = 0.0
for i in range(50):
    h = rng.standard_t(3, len(genes)) * 0.05
    t = np.flatnonzero(panel)[i]
    h[t] = -3.0                                                 # the row's own knockdown, the largest |h| in the row
    g = np.where(gate, h, 0.0)
    out = M.bulk_sharpen(g, 1.25, gate, panel)
    on = gate & ~panel
    assert abs(np.linalg.norm(out[on]) / np.linalg.norm(g[on]) - 1) < 1e-12
    assert np.array_equal(out[panel], g[panel]) and np.array_equal(out[~gate], g[~gate])
    ref = sharpen(np.where(panel, 0.0, g), 1.25, gate, mode="row")    # 414's priced object (panel zeroed)
    worst = max(worst, float(np.max(np.abs(out[on] - ref[on]))))
    assert np.array_equal(np.sign(out[on]), np.sign(g[on]))           # monotone: signs and |h| order kept
    assert np.array_equal(np.argsort(np.abs(out[on]), kind="stable"), np.argsort(np.abs(g[on]), kind="stable"))
assert worst < 1e-12, worst
# profile_from_effect: sharpen_q=None is the unpatched render; sharpen_q renders the sharpened gated row
h = rng.standard_t(3, len(genes)) * 0.05
p0 = M.profile_from_effect(h, b_ctrl, min_cpm=10.0)
p1 = M.profile_from_effect(h, b_ctrl, min_cpm=10.0, sharpen_q=1.25, protect=panel)
g = np.where(gate, h, 0.0)
man = np.clip(np.expm1(b_ctrl + M.bulk_sharpen(g, 1.25, gate, panel)), 0, None) / M.TS
assert np.max(np.abs(p1 - man / man.sum())) < 1e-15
man0 = np.clip(np.expm1(b_ctrl + g), 0, None) / M.TS
assert np.array_equal(p0, man0 / man0.sum())
try:
    M.profile_from_effect(h, b_ctrl, min_cpm=10.0, sharpen_q=1.25)
    raise AssertionError("sharpen_q without protect must raise")
except ValueError:
    pass
# the knob refuses without ballast, and refuses q outside (1, 2]
for bad, want in ((dict(VCC_BULK_SHARPEN_Q="1.25", VCC_PERT_SCALE="wshare_lo0.3_hi3"), "needs VCC_PERT_SCALE_BALLAST"),
                  (dict(env, VCC_BULK_SHARPEN_Q="0.8"), "must be in (1, 2]")):
    e = {k: v for k, v in os.environ.items() if not k.startswith("VCC_")}
    e.update(bad)
    r = subprocess.run([sys.executable, "-c", f"import importlib.util as u; s=u.spec_from_file_location('b', {str(BUILDER)!r}); "
                        "import sys; m=u.module_from_spec(s); sys.modules['b']=m; s.loader.exec_module(m)"], env=e, capture_output=True, text=True)
    assert r.returncode != 0 and want in (r.stderr + r.stdout), (bad, r.returncode, r.stderr[-300:])
print(f"417 bulk-sharpen tests PASS (50 rows: norm, panel/ungated untouched, == 169 sharpen(row) to {worst:.1e}, "
      "order/sign kept; render path; refusals)")
