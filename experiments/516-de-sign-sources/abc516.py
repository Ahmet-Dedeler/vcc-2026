"""516 (rider-class convenience, NOT a lever: the source did not pass the mandate's +0.02 H1-and-fold bar).
The Stouffer-6 sign source on the A/B/C panel, context-agnostic, on the challenge gene axis (gene_names.csv order):
Z = sum_i clip(x_i / (mad_g,i tau_p,i), +-3) / sqrt(n_i) over K562-GWPS, Orion HCT116, Orion HEK293T, CD4 Rest / Stim8hr /
Stim48hr (same code as vote516 / sources516; noise models from 2,000 random non-panel perturbations per bank).
Confidence: the H1 out-of-fold calibration of |Z| (work/stack516 eq rows) is written beside it as P(sign right | |Z| bin).
Known answer: the same function on the H1 axis must reproduce vote516's H1 Z exactly (checked before writing).
    PYTHONDONTWRITEBYTECODE=1 uv run python abc516.py
"""
import hashlib
import json

import numpy as np
import pandas as pd

import sources516 as SRC
from common516 import ROOT, WORK

NAMES = ["k562_gwps", "orion_hct116", "orion_hek293t", "cd4_Rest", "cd4_Stim8hr", "cd4_Stim48hr"]


def stz(targets, genes, exclude):
    parts = {}
    for bn in NAMES[:3]:
        r = SRC.bank(bn, targets, genes, exclude)
        parts[bn] = (r["x"].astype(float), r["mad"], r["tau"])
    c = SRC.cd4(targets, genes, exclude)
    for cn in ("Rest", "Stim8hr", "Stim48hr"):
        parts[f"cd4_{cn}"] = (c[cn]["log_fc"].astype(float), c[cn]["mad"], c[cn]["tau"])
    num = 0; den = 0
    for bn in NAMES:
        x, mad, tau = parts[bn]
        s = np.clip(x / (mad[None, :] * tau[:, None]), -3, 3)
        num = num + np.nan_to_num(s); den = den + np.isfinite(s)
    Z = np.where(den > 0, num / np.sqrt(np.maximum(den, 1)), np.nan)
    return Z, den


# known answer: H1 axis vs vote516's construction from the cached sources
import importlib.util as u
s = u.spec_from_file_location("v", "vote516.py"); V = u.module_from_spec(s); s.loader.exec_module(V)
Zv = V.stouffer(["K", "H", "E", "R", "S8", "S48"], False, 3.0)
Zh, _ = stz(V.d["targets"], V.d["genes"], set(V.d["targets"]))
dz = np.nanmax(np.abs(np.nan_to_num(Zh) - np.nan_to_num(Zv)))
print(f"KNOWN ANSWER {'PASS' if dz < 1e-5 else 'FAIL'}: H1 Z from abc516's path vs vote516's, max |dZ| {dz:.2e}")
assert dz < 1e-5
targets = pd.read_csv(ROOT / "data/challenge/pert_counts.csv")["target_gene"].astype(str).to_numpy()
genes = pd.read_csv(ROOT / "data/challenge/gene_names.csv")["gene_name"].astype(str).to_numpy()
Z, n = stz(targets, genes, set(targets))
Z = np.where(genes[None, :] == targets[:, None], np.nan, Z)          # the target's own gene is never scored
calib = {"|Z| [0,0.5)": 0.515, "|Z| [0.5,1)": 0.536, "|Z| [1,1.5)": 0.562, "|Z| [1.5,2)": 0.594, "|Z| [2,3)": 0.626,
         "|Z| >=3": 0.744}   # decomp516 all85 Stouffer-6 accuracies on H1 pool genes by |Z| bin (in-sample, H1 only)
out = WORK / "abc-stz6-signs.npz"
np.savez_compressed(out, targets=targets, genes=genes, Z=Z.astype(np.float32), n_sources=n.astype(np.int8))
sha = hashlib.sha256(open(out, "rb").read()).hexdigest()
meta = dict(file=str(out.relative_to(ROOT)), sha256=sha, construction="Stouffer-6 uncentred cap 3 (vote516/sources516)",
            status="rider-class only: H1 +0.016 pooled uncalled (no_ess +0.012), cancer folds ~0 vs the production analog",
            coverage_targets=int(np.isfinite(Z).any(1).sum()), coverage_cells=float(np.isfinite(Z).mean()),
            h1_accuracy_by_absZ=calib, recommended_use="sign for genes production does not call, |Z| >= 0.5 only")
json.dump(meta, open(WORK / "abc-stz6-signs.json", "w"), indent=1)
print(json.dumps(meta, indent=1))
