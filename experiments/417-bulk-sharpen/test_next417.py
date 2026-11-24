"""417 window 2: the 415 LODO key artifacts load through apply_pert_scale (weights == the npz's, non-panel energy kept);
VCC_BULK_CTX_SCALE parses, refuses without ballast / out of range / malformed. TEST_BUILDER=<path> tests a dry copy."""
import importlib.util
import os
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
BUILDER = Path(os.environ.get("TEST_BUILDER", ROOT / "scripts/make_submission_pm_cd4.py"))


def load(env):
    for k in [k for k in os.environ if k.startswith("VCC_")]:
        os.environ.pop(k)
    os.environ.update(env)
    s = importlib.util.spec_from_file_location("builder417n", BUILDER)
    m = importlib.util.module_from_spec(s)
    sys.modules["builder417n"] = m
    s.loader.exec_module(m)
    return m


perts = [str(t) for t in pd.read_csv(ROOT / "data/challenge/pert_counts.csv")["target_gene"]]
genes = [str(g) for g in pd.read_csv(ROOT / "data/challenge/gene_names.csv")["gene_name"]]
nonp = np.array([g not in set(perts) for g in genes])
E = np.random.default_rng(5).normal(0, 0.1, (300, len(genes)))
for key, npz in (("lodo415_lo0.3_hi3", "415-lodo-key/work/pert_scale_415pp_lo0.3_hi3.npz"), ("lodo415_lo0.3_inf", "415-lodo-key/work/pert_scale_415pp_lo0.3_hi1e+09.npz"),
                 ("wshare_lo0.5_hi3", "410-pertscale-revisit/work/pert_scale_410_lo0.5_hi3.npz"), ("wshare_lo0.5_inf", "410-pertscale-revisit/work/pert_scale_410_lo0.5_hi1e+09.npz")):
    M = load(dict(VCC_PERT_SCALE=key, VCC_PERT_SCALE_BALLAST="1"))
    out, audit = M.apply_pert_scale(E, perts, genes)
    z = np.load(ROOT / "experiments" / npz)
    w = dict(zip(map(str, z["targets"]), z["w"]))
    ratio = np.linalg.norm(out, axis=1) / np.linalg.norm(E, axis=1)
    want = np.array([w[t] for t in perts])
    assert np.allclose(ratio / ratio.mean(), want / want.mean(), rtol=1e-10), key
    assert abs(np.linalg.norm(out[:, nonp]) / np.linalg.norm(E[:, nonp]) - 1) < 1e-12, key
M = load(dict(VCC_PERT_SCALE="wshare_lo0.3_hi3", VCC_PERT_SCALE_BALLAST="1", VCC_BULK_CTX_SCALE="A=0.978,B=1.0136,C=1.0386"))
assert M.BULK_CTX_SCALE == {"A": 0.978, "B": 1.0136, "C": 1.0386}, M.BULK_CTX_SCALE
M = load(dict(VCC_PERT_SCALE="wshare_lo0.3_hi3", VCC_PERT_SCALE_BALLAST="1"))
assert M.BULK_CTX_SCALE is None and M.BULK_SHARPEN_Q is None
for bad, want in ((dict(VCC_PERT_SCALE="wshare_lo0.3_hi3", VCC_BULK_CTX_SCALE="A=1,B=1,C=1"), "needs VCC_PERT_SCALE_BALLAST"),
                  (dict(VCC_PERT_SCALE="wshare_lo0.3_hi3", VCC_PERT_SCALE_BALLAST="1", VCC_BULK_CTX_SCALE="A=1.5,B=1,C=1"), "[0.7, 1.4]"),
                  (dict(VCC_PERT_SCALE="wshare_lo0.3_hi3", VCC_PERT_SCALE_BALLAST="1", VCC_BULK_CTX_SCALE="A:1,B:1"), "must look like")):
    e = {k: v for k, v in os.environ.items() if not k.startswith("VCC_")}
    e.update(bad)
    r = subprocess.run([sys.executable, "-c", f"import importlib.util as u, sys; s=u.spec_from_file_location('b', {str(BUILDER)!r}); "
                        "m=u.module_from_spec(s); sys.modules['b']=m; s.loader.exec_module(m)"], env=e, capture_output=True, text=True)
    assert r.returncode != 0 and want in (r.stderr + r.stdout), (bad, r.returncode, r.stderr[-300:])
for k in [k for k in os.environ if k.startswith("VCC_")]:
    os.environ.pop(k)
print("417 window-2 tests PASS (415 LODO + 410 floor-0.5 key artifacts load and scale as their npz; ctx-scale parse + 3 refusals)")
