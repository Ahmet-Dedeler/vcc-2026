"""410b ballast helper tests (run after patch_ballast410.py --write): exact totals, per-cell cap, row count, K bounds, unset no-op."""
import importlib
import os
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts"))
os.environ.pop("VCC_PERT_SCALE_BALLAST", None)
os.environ.pop("VCC_PERT_SCALE", None)
import make_submission_pm_cd4 as M  # noqa: E402

assert M.PERT_SCALE_BALLAST is None and M.PERT_SCALE is None
# refuses ballast without a scale
os.environ["VCC_PERT_SCALE_BALLAST"] = "1"
try:
    importlib.reload(M)
    raise AssertionError("ballast without VCC_PERT_SCALE must be refused")
except SystemExit:
    pass
os.environ["VCC_PERT_SCALE"] = "wshare_lo0.3_hi3"
M = importlib.reload(M)
assert M.PERT_SCALE_BALLAST == "1"

rng = np.random.default_rng(3)
G = 2000
base = rng.gamma(0.5, 1.0, G)
base /= base.sum()
eff = rng.normal(0, 0.4, G)
prof_de = base * np.exp(eff); prof_de /= prof_de.sum()
for w in (0.19, 0.62, 1.0, 1.87):
    prof_bk = base * np.exp(w * eff); prof_bk /= prof_bk.sum()
    libs = rng.integers(12000, 30000, 400).astype(float)
    k, ratio = M.ballast_count(prof_de, prof_bk, libs, 400)
    assert 1 <= k <= M.BALLAST_KMAX
    carriers = np.vstack([rng.multinomial(int(l), prof_de) for l in libs[: 400 - k]])
    block, diag = M.ballast_rows(carriers, prof_de, prof_bk, k)
    assert block.shape == (400, G), block.shape
    assert np.array_equal(block[: 400 - k], carriers)
    assert block.sum(1).max() <= 1_000_000 and block.min() >= 0
    S = int(carriers.sum()) + k * M.BALLAST_LB
    tc = carriers.sum(0); tg = tc / tc.sum() * np.where(prof_de > 0, prof_bk / prof_de, 1.0)
    raw = tg / tg.sum() * S; tb = np.floor(raw).astype(np.int64); rem = int(S - tb.sum())
    tb[np.argsort(-(raw - tb), kind="stable")[:rem]] += 1
    tot = block.sum(0)
    ok = tb >= carriers.sum(0)
    if diag["ballast_trim"] == 0:
        assert np.array_equal(tot[ok], tb[ok]), "exact totals where feasible"
    assert np.array_equal(tot[~ok], carriers.sum(0)[~ok]), "clipped genes sit at the carriers' totals"
    jk_c = M.block_jackknife(carriers); jk_b = M.block_jackknife(block)
    assert k < 2 or jk_b >= 0.999 * jk_c or diag["ballast_jk_rho"] >= 0.95, (jk_b, jk_c, diag)
    print(f"   jk carriers {jk_c:.4f} -> block {jk_b:.4f} (ratio {diag['ballast_jk_ratio']:.3f}, rho {diag['ballast_jk_rho']:.3f})")
    comp = tot / tot.sum()
    rel = np.abs(np.log(comp[prof_bk > 1e-4] / prof_bk[prof_bk > 1e-4]))
    print(f"w {w}: K {k} ratio {ratio:.3f} clipped {diag['ballast_clip_genes']} max|log comp/bulk| {rel.max():.4f} "
          f"ballast libs {diag['ballast_lib_min']}..{diag['ballast_lib_max']}")
print("410b ballast tests PASS")
