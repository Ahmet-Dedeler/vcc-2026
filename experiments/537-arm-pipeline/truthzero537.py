#!/usr/bin/env python3
"""537 — is there a "truth-zero" nmae pot? A gate gene whose reference perturbed mean is exactly 0 carries lfc = log2(1e-9 /
ctrl CPM) ~ -32 and would dominate its target's nmae (a cis-silenced, low-expressed neighbour under CRISPRi, say). Predicting
the zero would win that target outright; the 433h zero repair shows the reverse error's size. Counted on h1-test's reference
(323's reference.npz: the scorer's own p, q, lfc on the 5-CPM universe). Writes work/truthzero537.json. < 0.5 GB.
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]


def main():
    z = np.load(ROOT / "experiments/323-current-family-validation/work/reference.npz", allow_pickle=False)
    lfc, q = z["lfc"], z["q"]
    genes, tg = z["genes"][z["keep"]].astype(str), z["targets"].astype(str)
    sig = (q < 0.05) & np.isfinite(lfc) & ~np.equal.outer(tg, genes)
    out = {"gate_genes_total": int(sig.sum()), "gate_genes_per_target_median": float(np.median(sig.sum(1)))}
    for thr in (-15, -10, -5):
        m = sig & (lfc < thr)
        out[f"lt{-thr}"] = int(m.sum())
        out[f"targets_lt{-thr}"] = int((m.sum(1) > 0).sum())
    out["gt5"] = int((sig & (lfc > 5)).sum())
    out["examples_lt5"] = [(tg[i], genes[j], round(float(lfc[i, j]), 2)) for i, j in zip(*np.nonzero(sig & (lfc < -5)))]
    print(json.dumps(out, indent=1))
    (HERE / "work/truthzero537.json").write_text(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()
