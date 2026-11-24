"""462e — does 458d's denoised programme recover the H1 recipient programme better than 458's bank? h1-test gate signs
(441's instrument, 58 targets) from P_down over the OTHER H1 perturbations (h1-test's 85 targets excluded), for:
  bank   data/donors/h1_2025 (harness: 400 cells/target, 18,399 control rows, thinned)
  full   data/donors/h1_2025_full (458d: all public cells, 38,176 controls)
against the LOTO oracle (0.6245) and CLEAN4 (0.5733). Same-line and same-lab, so an UPPER reference for what either
programme can carry; the cross-line question is the board's (458d). usage: uv run python .../h1prog462.py"""
import glob, json, sys
from pathlib import Path
import numpy as np
ROOT = Path(__file__).resolve().parents[2]; HERE = Path(__file__).resolve().parent
z = np.load(glob.glob(str(ROOT / "experiments/323-*/work/reference.npz"))[0], allow_pickle=True)
genes, keep, targets = z["genes"].astype(str), z["keep"], list(z["targets"].astype(str))
q, lfc = z["q"], z["lfc"]; kg = genes[keep]
own = np.array([[g == t for g in kg] for t in targets])
gate = (q < 0.05) & np.isfinite(lfc) & (lfc != 0) & ~own
use = np.flatnonzero(gate.sum(1) >= 10)
res = {}
for name in ("h1_2025", "h1_2025_full"):
    d = ROOT / "data/donors" / name
    delta = np.load(d / "delta.npy", mmap_mode="r"); perts = [l.strip() for l in open(d / "perts.txt")]
    bg = [l.strip() for l in open(d / "genes.txt")]
    rows = [i for i, p in enumerate(perts) if p not in set(targets)]
    M = np.asarray(delta[rows], dtype=np.float64)
    nz = (M != 0).sum(0); dn = (M < 0).sum(0)
    with np.errstate(invalid="ignore", divide="ignore"):
        pdn = np.where(nz > 0, dn / np.maximum(nz, 1), np.nan)
    pos = {g: i for i, g in enumerate(bg)}; col = np.array([pos.get(g, -1) for g in kg])
    pk = np.where(col >= 0, pdn[np.maximum(col, 0)], np.nan)
    sign = np.where(np.isfinite(pk), np.where(pk > 0.5, -1.0, 1.0), np.nan)
    hit = n = 0; whit = wsum = 0.0; cov = []
    for t in use:
        g = gate[t] & np.isfinite(sign)
        h = sign[g] == np.sign(lfc[t][g]); w = np.abs(lfc[t][g])
        hit += h.sum(); n += g.sum(); whit += (h * w).sum(); wsum += w.sum(); cov.append(g.sum() / gate[t].sum())
    res[name] = dict(rows=len(rows), gate_acc=round(hit / n, 4), gate_acc_weighted=round(whit / wsum, 4), coverage=round(float(np.mean(cov)), 3))
res["reference"] = {"oracle_loto": [0.6245, 0.6324], "clean4": [0.5733, 0.5713]}
print(json.dumps(res, indent=1)); (HERE / "work/h1prog462.json").write_text(json.dumps(res, indent=1))
