"""450l: where does ARC's own programme sit? P_down from Arc's 2025 H1 screen (data/donors/h1_2025: 255 perturbations, Flex,
400 cells each; A/B/C panel rows excluded) against the other labs' programmes (labline450's banks) and Pan's H1 (same line,
Whitehead). If programmes are lab signatures, Arc-H1 should agree with nothing cross-lab except through the H1 line itself
(Pan). usage: uv run python experiments/450-same-line-programme/arcprog450.py"""
import csv
import sys
from pathlib import Path

import h5py
import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts"))
import def_sameline as S  # noqa: E402

S.BANK_DIRS.setdefault("orion_hek293t", ROOT / "data/donors/orion_hek293t")
S.BANK_DIRS.setdefault("rpe1", ROOT / "data/donors/rpe1")
S.BANK_DIRS.setdefault("h1_2025", ROOT / "data/donors/h1_2025")
official = [r["gene_name"] for r in csv.DictReader(open(ROOT / "data/challenge/gene_names.csv"))]
panel = {r["target_gene"] for r in csv.DictReader(open(ROOT / "data/challenge/pert_counts.csv"))}
P = {}
for b in ("h1_2025", "k562_gwps", "k562_essential", "nadig_jurkat", "nadig_hepg2", "rpe1", "orion_hct116", "orion_hek293t", "viperturb_k562"):
    delta, perts, genes, _ = S.load_bank(b)
    _, pdown, nrows, _ = S.programme(delta, perts, genes, official, panel)
    P[b] = pdown
    print(f"{b:15s} rows {nrows}", flush=True)
# Pan H1 (Whitehead, 10x 5'), from the pseudobulk database, A/B/C panel rows excluded
db = ROOT / "data/external/gse295214/pseudo_bulk_perturbation_database.h5ad"
def cat(g):
    return g["categories"][:].astype(str)[g["codes"][:]] if isinstance(g, h5py.Group) else g[:].astype(str)
with h5py.File(db, "r") as f:
    ds, pg, vg = cat(f["obs/dataset"]), cat(f["obs/perturbed_gene_name"]), cat(f["var/gene_name"])
    sel = np.sort(np.where((ds == "Pan2025_hESC") & ~np.isin(pg, list(panel)))[0])
    X = np.where(f["layers/measure_mask"][sel, :] > 0, f["X"][sel, :], np.nan)
with np.errstate(invalid="ignore"):
    pp = np.nanmean(np.where(np.isnan(X), np.nan, (X < 0).astype(float)), axis=0)
vi = {g: i for i, g in enumerate(vg)}
col = np.array([vi.get(g, -1) for g in official])
P["pan_h1"] = np.where(col >= 0, pp[np.maximum(col, 0)], np.nan)
a = P["h1_2025"] - 0.5
print("\nArc H1 programme vs:")
for k, v in P.items():
    if k == "h1_2025":
        continue
    b = v - 0.5
    m = np.isfinite(a) & np.isfinite(b)
    s = m & (np.abs(a) >= 0.05) & (np.abs(b) >= 0.05)
    print(f"  {k:15s} n {m.sum():6d}  r {np.corrcoef(a[m], b[m])[0, 1]:+.3f}  sign {np.mean(np.sign(a[s]) == np.sign(b[s])):.3f} (n {s.sum()})")
