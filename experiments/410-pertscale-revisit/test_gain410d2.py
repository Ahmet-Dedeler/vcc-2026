"""410d2: the render-gated normalisation preserves the summed gated energy exactly; 410d's value is unchanged."""
import os, sys, importlib
from pathlib import Path
import numpy as np, pandas as pd
ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts"))
os.environ.update(VCC_PERT_SCALE="wshare_lo0.3_inf", VCC_PERT_SCALE_BALLAST="1", VCC_BULK_GENE_GAIN="g228reg_t025rg")
import make_submission_pm_cd4 as M
M = importlib.reload(M)
perts = [str(t) for t in pd.read_csv(ROOT / "data/challenge/pert_counts.csv")["target_gene"]]
genes = [str(g) for g in pd.read_csv(ROOT / "data/challenge/gene_names.csv")["gene_name"]]
gates = [M.control_render_gate(ROOT / "data/challenge", c, genes, 10.0) for c in "ABC"]
print("gated genes per context", [int(g.sum()) for g in gates])
E = np.random.default_rng(2).normal(0, 0.1, (300, len(genes)))
out, audit = M.apply_bulk_gene_gain(E, perts, genes, gates=gates)
nonp = np.array([g not in set(perts) for g in genes])
e_before = sum((E[:, g & nonp] ** 2).sum() for g in gates); e_after = sum((out[:, g & nonp] ** 2).sum() for g in gates)
assert abs(e_after / e_before - 1) < 1e-12, e_after / e_before
out0, audit0 = M.apply_bulk_gene_gain(E, perts, genes, gates=None)
assert abs(np.linalg.norm(out0[:, nonp]) / np.linalg.norm(E[:, nonp]) - 1) < 1e-12
print("410d2 gated-normalisation tests PASS; constants gated", round(audit["energy_constant"], 6), "effect-space", round(audit0["energy_constant"], 6), audit["normalisation"])
for k in ("VCC_BULK_GENE_GAIN", "VCC_PERT_SCALE_BALLAST", "VCC_PERT_SCALE"):
    os.environ.pop(k, None)
