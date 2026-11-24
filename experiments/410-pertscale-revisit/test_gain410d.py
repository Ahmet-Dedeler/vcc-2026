"""410d knob tests: refused without ballast, unset no-op, set multiplies the bulk per gene energy-preservingly."""
import os, sys, importlib, hashlib
from pathlib import Path
import numpy as np, pandas as pd
ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts"))
for k in ("VCC_BULK_GENE_GAIN", "VCC_PERT_SCALE_BALLAST", "VCC_PERT_SCALE"):
    os.environ.pop(k, None)
import make_submission_pm_cd4 as M
assert M.BULK_GENE_GAIN is None
os.environ["VCC_BULK_GENE_GAIN"] = "g228reg_t025"
try:
    importlib.reload(M); raise AssertionError("gene gain without ballast must be refused")
except SystemExit:
    pass
os.environ["VCC_PERT_SCALE"] = "wshare_lo0.3_inf"; os.environ["VCC_PERT_SCALE_BALLAST"] = "1"
M = importlib.reload(M)
assert M.BULK_GENE_GAIN == "g228reg_t025"
path, pin = M.BULK_GENE_GAIN_ARTIFACTS[M.BULK_GENE_GAIN]
assert hashlib.sha256(path.read_bytes()).hexdigest() == pin
perts = [str(t) for t in pd.read_csv(ROOT / "data/challenge/pert_counts.csv")["target_gene"]]
genes = [str(g) for g in pd.read_csv(ROOT / "data/challenge/gene_names.csv")["gene_name"]]
E = np.random.default_rng(1).normal(0, 0.1, (300, len(genes)))
out, audit = M.apply_bulk_gene_gain(E, perts, genes)
nonp = np.array([g not in set(perts) for g in genes])
assert abs(np.linalg.norm(out[:, nonp]) - np.linalg.norm(E[:, nonp])) < 1e-9 * np.linalg.norm(E[:, nonp])
m = np.load(path)["m"]
ratio = out[:, 5] / E[:, 5]
assert np.allclose(ratio, m[5] * audit["energy_constant"])
for k in ("VCC_BULK_GENE_GAIN", "VCC_PERT_SCALE_BALLAST", "VCC_PERT_SCALE"):
    os.environ.pop(k, None)
M = importlib.reload(M); assert M.BULK_GENE_GAIN is None
print("410d gene-gain knob tests PASS", audit["m_quantiles"])
