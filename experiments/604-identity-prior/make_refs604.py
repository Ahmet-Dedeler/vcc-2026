#!/usr/bin/env python3
"""604: pin the Flex reference profiles the identity prior uses (622: Flex only; 3' references hurt).
  cpm_A / cpm_B / cpm_C : per-cell CPM mean of the challenge contexts' controls (data/challenge, float64, every gene)
  cpm_h1                : 480's gate-h1_full_ka.npz `ctrl_cpm` (Arc Flex 2025 H1 controls; its > 5 CPM genes; NaN elsewhere)
  calls                 : the identity call per profile (A/B/C: def_identity's consensus on the challenge bundle; h1: H1_ESC,
                          not a CCLE line), used by the same-line-reference guard.
Writes work/flex_refs604.npz (+ .json with sha256 and the input shas)."""
import csv, hashlib, json, sys
from pathlib import Path
import numpy as np
HERE = Path(__file__).resolve().parent; ROOT = HERE.parents[1]
sys.path.insert(0, str(HERE))
import def_identity_prior as D
DATA = ROOT / "data/challenge"
H1 = ROOT / "experiments/480-nmae-instrument/work/gate-h1_full_ka.npz"
CALLS = {"A": "JURKAT_HAEMATOPOIETIC_AND_LYMPHOID_TISSUE", "B": "HELA_CERVIX", "C": "CAL33_UPPER_AERODIGESTIVE_TRACT", "h1": "H1_ESC"}
with open(DATA / "gene_names.csv") as fh:
    genes = [r["gene_name"] for r in csv.DictReader(fh)]
out = {"genes": np.asarray(genes), "names": np.asarray(["A", "B", "C", "h1"]), "calls": np.asarray([CALLS[k] for k in ("A", "B", "C", "h1")])}
inputs = {}
for c in "ABC":
    p = DATA / f"context_{c}.h5ad"
    out[f"cpm_{c}"] = D.ctrl_cpm(p, genes)
    inputs[f"context_{c}.h5ad"] = D.file_sha256(p)
z = np.load(H1, allow_pickle=True)
pos = dict(zip(z["genes"].astype(str), z["ctrl_cpm"].astype(np.float64)))
out["cpm_h1"] = np.array([pos.get(g, np.nan) for g in genes])
inputs["gate-h1_full_ka.npz"] = D.file_sha256(H1)
dst = HERE / "work/flex_refs604.npz"
np.savez_compressed(dst, **out)
sha = D.file_sha256(dst)
meta = {"sha256": sha, "inputs_sha256": inputs, "calls": CALLS, "h1_genes": int(np.isfinite(out["cpm_h1"]).sum())}
Path(str(dst) + ".json").write_text(json.dumps(meta, indent=1))
print(json.dumps(meta, indent=1))
