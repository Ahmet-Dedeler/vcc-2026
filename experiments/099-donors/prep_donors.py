"""Materialise `Donor`-loadable directories for two candidates that ship without a `b_ctrl`.

Neither `b_ctrl` is used by anything this experiment runs: `eb=False` and the weighting is
`uniform`, so `MultiDonorModel` touches `b_ctrl` only in `basal_weights()`, which is never
called. A zero vector is written so `Donor.load` succeeds, and the fact that it is a
PLACEHOLDER is recorded here rather than inferred later.
"""
import numpy as np, os
from pathlib import Path
ROOT = Path("~/Code/virtual-cell-challenge-2026")
W = ROOT / "experiments/099-donors/work"

# --- PerturbMatch: Fortran-order (pert, gene) lfc, no control profile in the release --------
src = ROOT / "data/external/perturbmatch"
dst = W / "perturbmatch"; dst.mkdir(parents=True, exist_ok=True)
d = np.load(src / "delta.npy", mmap_mode="r")
print("perturbmatch delta", d.shape, d.dtype, "F" if d.flags.f_contiguous else "C")
for f in ("delta.npy", "genes.txt", "perts.txt"):
    t = dst / f
    if not t.exists():
        t.symlink_to(src / f)
np.save(dst / "b_ctrl.npy", np.zeros(d.shape[1], dtype=np.float32))

# --- Feng iPSC: pooled genome-wide screen lfc on the VCC axis, 6,151 genes measured ---------
src = ROOT / "data/external/feng-ipsc/genomewide_vcc_panel"
dst = W / "feng_ipsc"; dst.mkdir(parents=True, exist_ok=True)
lfc = np.load(src / "lfc.npy")
print("feng lfc", lfc.shape, lfc.dtype, "measured genes",
      int(np.load(src / "measured_gene_mask.npy").sum()))
np.save(dst / "delta.npy", lfc.astype(np.float32))
(dst / "genes.txt").write_text((src / "genes.txt").read_text())
(dst / "perts.txt").write_text((src / "targets.txt").read_text())
np.save(dst / "b_ctrl.npy", np.zeros(lfc.shape[1], dtype=np.float32))
print("done")
