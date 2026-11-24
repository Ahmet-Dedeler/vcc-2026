"""642: build the pinned non-CCLE HEK293 reference row for def_identity's second stage.

Sources: every HEK293 / 293T non-targeting control profile on disk that is NOT Flex (642's rule: a Flex-built reference
pulls Satija's Flex HAP1 to a unanimous 8/8 call at 0.818; same-platform references attract same-platform queries):
  xu      XuCao 2023 HEK293, 10x 3', 2,758 control cells (590's profile, official-axis counts)
  orion   X-Atlas/Orion HEK293T b_ctrl (log1p CP50k), Orion FiCS
Each source -> log2(1 + CPM) over its own present genes (def_identity's query arithmetic: counts_query / bctrl_query),
then the mean in log2 space where both are present, else the one present; NaN where neither. Stored on the official
gene axis (data/challenge/gene_names.csv order); def_identity maps it onto its common axis at match time.

    uv run python experiments/642-hek293-contingency/build_ref642.py   -> data/identity_refs/hek293_nonflex_ref642.npz
"""
from __future__ import annotations

import hashlib
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
OUT = ROOT / "data/identity_refs/hek293_nonflex_ref642.npz"
XU = ROOT / "experiments/590-identity-validation/work/profiles/xu.npz"
ORION = ROOT / "data/donors/orion_hek293t"


def sha(p: Path) -> str:
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(1 << 24), b""):
            h.update(b)
    return h.hexdigest()


def main() -> None:
    official = pd.read_csv(ROOT / "data/challenge/gene_names.csv")["gene_name"].astype(str).tolist()
    z = np.load(XU)
    meta = json.loads(str(z["meta"]))
    assert len(meta) == 1 and meta[0]["query"] == "xu:HEK293" and "Flex" not in meta[0]["chemistry"], meta
    c = np.asarray(z["counts"][0], dtype=np.float64)
    m = np.isfinite(c)
    xu = np.full(len(official), np.nan)
    xu[m] = np.log2(1.0 + 1e6 * c[m] / c[m].sum())                      # = match590.counts_query on the official axis
    og = np.asarray((ORION / "genes.txt").read_text().splitlines())
    b = np.expm1(np.load(ORION / "b_ctrl.npy").astype(np.float64))
    ol = np.log2(1.0 + 1e6 * b / b.sum())                              # = def_identity.log2_cpm_from_bctrl
    pos = {g: i for i, g in enumerate(og)}
    orion = np.array([ol[pos[g]] if g in pos else np.nan for g in official])
    both = np.isfinite(xu) & np.isfinite(orion)
    val = np.where(both, (xu + orion) / 2, np.where(np.isfinite(xu), xu, orion))
    prov = {
        "name": "hek293_nonflex_ref642", "line": "HEK293 (HEK293 + HEK293T controls)", "built_by": "experiments/642-hek293-contingency/build_ref642.py",
        "units": "log2(1 + CPM) per source over its own present genes, then the log2-space mean (official gene axis)",
        "sources": [
            {"name": "xu", "dataset": "XuCao2023 (scPerturb XuCao2023.h5ad, cell_line HEK293, control cells)", "chemistry": "10x 3'",
             "n_cells": meta[0]["n_cells"], "file": str(XU.relative_to(ROOT)), "sha256": sha(XU)},
            {"name": "orion", "dataset": "X-Atlas/Orion HEK293T non-targeting b_ctrl (data/donors/orion_hek293t)",
             "chemistry": "Orion FiCS (10x GEM-X 5')", "file": str((ORION / "b_ctrl.npy").relative_to(ROOT)),
             "sha256": sha(ORION / "b_ctrl.npy"), "genes_sha256": sha(ORION / "genes.txt")},
        ],
        "rule": "642: never build a non-CCLE reference from Flex data (same-platform references attract same-platform queries)",
        "n_genes": len(official), "n_present": int(np.isfinite(val).sum()), "n_both": int(both.sum()),
    }
    OUT.parent.mkdir(parents=True, exist_ok=True)
    np.savez(OUT, genes=np.asarray(official), values=val.astype(np.float32), provenance=np.asarray(json.dumps(prov)))
    print(json.dumps({**prov, "out": str(OUT.relative_to(ROOT)), "out_sha256": sha(OUT)}, indent=1))


if __name__ == "__main__":
    main()
