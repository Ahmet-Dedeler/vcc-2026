"""Donor inventory: coverage of the 300 challenge targets and of the 18,533 gene axis."""
import json, itertools
from pathlib import Path
import numpy as np

ROOT = Path("~/Code/virtual-cell-challenge-2026")
DON = ROOT / "data/donors"
NAMES = ["k562_gwps", "orion_hct116", "orion_hek293t", "h1_2025",
         "k562_essential", "nadig_jurkat", "rpe1"]

targets = [l.strip() for l in open(ROOT / "data/challenge/pert_counts.csv")][1:]
targets = [t for t in targets if t]
axis = [l.strip() for l in open(ROOT / "data/challenge/gene_names.csv")]
axis = [g for g in axis if g and g != "gene"]
print(f"targets={len(targets)} axis={len(axis)}")
TS_ = set(targets); AX = set(axis)

rows = {}
for n in NAMES:
    d = DON / n
    perts = [l.strip() for l in open(d / "perts.txt") if l.strip()]
    genes = [l.strip() for l in open(d / "genes.txt") if l.strip()]
    meta = json.loads((d / "meta.json").read_text())
    cov = sorted(TS_ & set(perts))
    rows[n] = dict(n_perts=len(perts), n_genes=len(genes),
                   axis_overlap=len(AX & set(genes)),
                   cov=cov, n_cov=len(cov),
                   ctrl_cells=meta.get("n_control_cells"),
                   med_cells=meta.get("cells_per_pert", {}).get("median"),
                   src=meta.get("source"))
    print(f"{n:16s} perts={len(perts):6d} genes={len(genes):6d} "
          f"axis_ov={len(AX & set(genes)):6d} targets={len(cov):3d}/300 "
          f"ctrl={meta.get('n_control_cells')} medcells={meta.get('cells_per_pert',{}).get('median')}")

np.save(ROOT / "experiments/099-donors/work/cov.npy",
        np.array([[t in set(rows[n]['cov']) for t in targets] for n in NAMES]))
json.dump({n: rows[n]['cov'] for n in NAMES},
          open(ROOT / "experiments/099-donors/work/cov.json", "w"))

print("\n-- pairwise target intersection --")
print("            " + " ".join(f"{n[:9]:>10s}" for n in NAMES))
for a in NAMES:
    print(f"{a:12s}" + " ".join(f"{len(set(rows[a]['cov'])&set(rows[b]['cov'])):10d}" for b in NAMES))

print("\n-- how many donors carry each target --")
cnt = np.zeros(len(targets), int)
for i, t in enumerate(targets):
    cnt[i] = sum(t in set(rows[n]['cov']) for n in NAMES)
for k in range(8):
    print(f"  carried by {k} donors: {(cnt==k).sum()}")
prod = ["k562_gwps", "orion_hct116", "orion_hek293t"]
cntp = np.array([sum(t in set(rows[n]['cov']) for n in prod) for t in targets])
print("\n-- production 3 --")
for k in range(4):
    print(f"  carried by {k}: {(cntp==k).sum()}")
