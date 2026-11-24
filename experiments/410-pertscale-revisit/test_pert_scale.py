"""410 knob tests: unset is a no-op (no key, no module state), set scales rows energy-preservingly and refuses bad pins."""
import os, sys, importlib, hashlib
from pathlib import Path
import numpy as np
ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts"))
os.environ.pop("VCC_PERT_SCALE", None)
import make_submission_pm_cd4 as M
assert M.PERT_SCALE is None
os.environ["VCC_PERT_SCALE"] = "wshare_lo0.3_hi3"
M2 = importlib.reload(M)
assert M2.PERT_SCALE == "wshare_lo0.3_hi3"
path, pin = M2.PERT_SCALE_ARTIFACTS[M2.PERT_SCALE]
assert hashlib.sha256(path.read_bytes()).hexdigest() == pin
import pandas as pd
perts = [str(t) for t in pd.read_csv(ROOT / "data/challenge/pert_counts.csv")["target_gene"]]
genes = [str(g) for g in pd.read_csv(ROOT / "data/challenge/gene_names.csv")["gene_name"]]
rng = np.random.default_rng(0)
E = rng.normal(0, 0.1, (300, len(genes)))
out, audit = M2.apply_pert_scale(E, perts, genes)
panel = set(perts); nonp = np.array([g not in panel for g in genes])
assert abs(np.linalg.norm(out[:, nonp]) - np.linalg.norm(E[:, nonp])) < 1e-9 * np.linalg.norm(E[:, nonp])
w = np.load(path)["w"]; ratio = np.linalg.norm(out, axis=1) / np.linalg.norm(E, axis=1)
assert np.allclose(ratio / ratio.mean(), w / w.mean())
# direction of every row unchanged
cos = (out * E).sum(1) / np.linalg.norm(out, axis=1) / np.linalg.norm(E, axis=1)
assert np.allclose(cos, 1.0)
try:
    M2.apply_pert_scale(E, perts[::-1], genes); raise AssertionError("order check missing")
except SystemExit:
    pass
os.environ.pop("VCC_PERT_SCALE")
M3 = importlib.reload(M2); assert M3.PERT_SCALE is None
print("410 knob tests PASS", audit["final_scale_quantiles"])
