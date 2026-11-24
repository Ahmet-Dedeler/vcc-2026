"""Shared-target counts between every candidate line and the donor tables we already hold.

Panel overlap decides whether a line can serve at INFERENCE. It is not the only thing a line
is worth: experiment 161's residual model learns from line PAIRS, and a pair only needs shared
*non-panel* targets. The brief makes this point about Jiang 2025 (9/300 panel, 171+ shared
non-panel targets). This script computes that second number for every candidate.

Writes `work/shared.json`.
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd

from census import PANEL, REPO, inspect, read_col  # noqa: F401  (PANEL/REPO reused)
import h5py

OUT = Path(__file__).resolve().parent / "work"

DONORS = REPO / "data/donors"
# Candidate lines: every on-disk source that is NOT already a donor table, has a usable
# single-gene target vocabulary, and is a cell line we do not already hold.
CANDIDATES = {
    "nadig_hepg2": ("data/external/nadig2025/hepg2_raw_singlecell.h5ad", "liver / hepatoblastoma"),
    "frangieh_melanoma": ("data/external/scperturb/FrangiehIzar2021_RNA.h5ad", "melanoma"),
    "xucao_hek293": ("data/external/scperturb/XuCao2023.h5ad", "kidney (HEK293)"),
    "tian_ipsc_neuron": ("data/external/scperturb/TianKampmann2021_CRISPRi.h5ad", "iPSC-derived neuron"),
    "sunshine_calu3": ("data/external/scperturb/SunshineHein2023.h5ad", "lung (Calu-3)"),
    "papalexi_thp1": ("data/external/scperturb/PapalexiSatija2021_eccite_RNA.h5ad", "monocyte (THP-1)"),
    "wessels_thp1": ("data/external/scperturb/WesselsSatija2023.h5ad", "monocyte (THP-1)"),
    "norman_k562": ("data/external/scperturb/NormanWeissman2019_filtered.h5ad", "CML (K562)"),
    "datlinger_jurkat": ("data/external/scperturb/DatlingerBock2017.h5ad", "T-ALL (Jurkat)"),
    "dixit_k562": ("data/external/scperturb/DixitRegev2016.h5ad", "CML (K562)"),
}

PERT_COLS = ("perturbation", "gene", "target_gene", "gene_target")


def labels(path: Path) -> set[str]:
    with h5py.File(path, "r") as h:
        col = next(c for c in PERT_COLS if c in h["obs"])
        p = read_col(h, "obs", col)
    out: set[str] = set()
    for s in set(p.tolist()):
        s = str(s)
        if s.lower() in {"control", "non-targeting", "nan", "*", "none", "ntc"}:
            continue
        # keep single-gene labels only; multi-guide cells are not a clean per-target effect
        if "_" in s or "+" in s or "," in s or ":" in s:
            continue
        out.add(s)
    return out


def main() -> None:
    held = {}
    for d in sorted(DONORS.iterdir()):
        if (d / "perts.txt").exists():
            held[d.name] = set((d / "perts.txt").read_text().split("\n"))
    print("held donor tables:", {k: len(v) for k, v in held.items()}, "\n")

    rows = []
    for name, (rel, tissue) in CANDIDATES.items():
        p = REPO / rel
        if not p.exists():
            print(f"MISSING {rel}")
            continue
        L = labels(p)
        rec = {
            "name": name, "tissue": tissue, "file": rel,
            "size_gb": round(p.stat().st_size / 1e9, 2),
            "n_single_gene_targets": len(L),
            "panel_overlap": len(PANEL & L),
            "shared_with": {k: len(L & v) for k, v in held.items()},
        }
        rec["shared_nonpanel_max"] = max(
            len((L & v) - PANEL) for v in held.values()) if held else 0
        rows.append(rec)
        print(f"{name:<20} {tissue:<26} targets={len(L):>5} panel={rec['panel_overlap']:>3} "
              f"| shared: " + " ".join(f"{k.split('_')[0][:6]}={v}"
                                       for k, v in rec["shared_with"].items()))

    (OUT / "shared.json").write_text(json.dumps(
        {"held": {k: len(v) for k, v in held.items()}, "candidates": rows}, indent=2))
    print(f"\nwrote {OUT/'shared.json'}")


if __name__ == "__main__":
    main()
