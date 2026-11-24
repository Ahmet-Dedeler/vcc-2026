"""DEF-RUNBOOK section 5: before relaxing 047's CD4 gates, CONFIRM THE MAPPING IS RIGHT.

Reads only the obs axis of the pinned 16.8 GB CD4 source (h5py, no X), and checks every ef-bundle
target that the extractor could not map: is it absent from the atlas, or present under another
case / a known alias?  Known-answer rows: A/B/C maps 283/287/284 (union 291); ef-bundle 186.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import h5py
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
SRC = ROOT / "experiments/199-def-runbook/work/GWCD4i.DE_stats.h5ad"


def read_col(h: h5py.File, key: str) -> np.ndarray:
    obj = h["obs"][key]
    if isinstance(obj, h5py.Group):  # categorical
        cats = obj["categories"][()].astype(str)
        codes = obj["codes"][()]
        return cats[codes]
    return obj[()].astype(str)


def main() -> None:
    with h5py.File(SRC, "r") as h:
        keys = list(h["obs"].keys())
        print("obs keys:", keys)
        obs = {k: read_col(h, k) for k in keys if k != "_index"}
        idx = h["obs"]["_index"][()].astype(str) if "_index" in h["obs"] else None
    df = pd.DataFrame(obs)
    print(df.head(3).to_string())
    print("n obs rows", len(df))
    # the target column: the extractor's own choice is in experiments/047-cd4-direct/extract_cd4.py;
    # here take every string column and report which one carries gene-symbol-like values
    cands = {}
    for c in df.columns:
        vals = set(df[c].unique())
        cands[c] = len(vals)
    print("unique per column:", cands)
    tcol = "target_contrast_gene_name"  # the extractor's symbol column (047); "index" is target x condition
    targets = set(df[tcol].unique())
    print(f"target column guess {tcol!r}: {len(targets)} unique")
    conds = [c for c in df.columns if df[c].nunique() <= 5]
    print("condition-like columns:", {c: sorted(df[c].unique()) for c in conds})
    out = {"source": str(SRC), "target_column": tcol, "n_targets": len(targets)}
    for tag, bundle in (("abc", ROOT / "data/challenge"),
                        ("ef238", ROOT / "experiments/238-essential-free-panel/work/ef-bundle")):
        panel = [str(t) for t in pd.read_csv(bundle / "pert_counts.csv")["target_gene"]]
        mapped = [t for t in panel if t in targets]
        missing = [t for t in panel if t not in targets]
        lower = {t.lower(): t for t in targets}
        case_hits = {t: lower[t.lower()] for t in missing if t.lower() in lower}
        per_cond = {}
        for c in conds:
            for cv in sorted(df[c].unique()):
                sub = set(df.loc[df[c] == cv, tcol])
                per_cond[f"{c}={cv}"] = sum(1 for t in panel if t in sub)
        out[tag] = {"n_panel": len(panel), "mapped_any_condition": len(mapped), "missing": len(missing),
                    "missing_present_under_other_case": case_hits, "per_condition": per_cond,
                    "missing_list": missing}
        print(f"[{tag}] mapped {len(mapped)}/{len(panel)}  missing {len(missing)}  case-variant hits {case_hits}")
        print(f"[{tag}] per condition {per_cond}")
    (HERE / "work/mapping244.json").write_text(json.dumps(out, indent=1))
    ok = out["abc"]["mapped_any_condition"] == 291 and out["ef238"]["mapped_any_condition"] == 186 \
        and not out["ef238"]["missing_present_under_other_case"]
    print("MAPPING", "PASS" if ok else "FAIL")
    raise SystemExit(0 if ok else 1)


if __name__ == "__main__":
    main()
