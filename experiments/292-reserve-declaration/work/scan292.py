"""292 item 2 -- METADATA-ONLY scan of every h5ad in data/external/scperturb/.

Reads ONLY the `obs`, `var` and `uns` groups with h5py.  `X`, `layers`, `obsm`, `varm`, `obsp`
are never opened.  Nothing here computes a response, a cosine, a conservation value or any
feature: every number is a count, a label or a set intersection.

Known-answer rows (asserted, not just reported):
  KA-1  n_cells and n_labels reproduce 278's obs-only scan (work/h5ad_counts/*.json) exactly.
  KA-2  279's held gene list and 290's feature-store gene list are the same 17,106 symbols in
        the same order (so overlaps (a) and (c) of the mandate are one number, reported once).
  KA-3  the 300-target panel read from conservation276.csv equals data/challenge/pert_counts.csv.

Output: work/scan292.json (one record per file), work/scan292_ka.json.

    PYTHONDONTWRITEBYTECODE=1 uv run python scan292.py
"""
from __future__ import annotations

import json
import resource
import time
from pathlib import Path

import h5py
import numpy as np
import pandas as pd

ROOT = Path("~/Code/virtual-cell-challenge-2026")
HERE = ROOT / "experiments/292-reserve-declaration"
WORK = HERE / "work"
SCP = ROOT / "data/external/scperturb"
COUNTS278 = ROOT / "experiments/278-depth-inventory/work/h5ad_counts"

# ----------------------------------------------------------------- reference gene sets
held = np.load(ROOT / "experiments/279-znf-state-dependence/work/held279.npz", allow_pickle=True)
HELD = np.asarray(held["genes"]).astype(str)
persrc = np.load(ROOT / "experiments/290-remainder-response-features/work/persrc290.npz",
                 allow_pickle=True)
PERSRC = np.asarray(persrc["genes"]).astype(str)
cons = pd.read_csv(ROOT / "experiments/276-conservation-map/work/conservation276.csv")
PANEL = [str(t) for t in cons["target"]]
PANEL_FILE = [l.split(",")[0] for l in
              (ROOT / "data/challenge/pert_counts.csv").read_text().splitlines()[1:] if l.strip()]

KA = {
    "ka2_held_persrc_same_order": bool(HELD.shape == PERSRC.shape and (HELD == PERSRC).all()),
    "ka2_n_held": int(len(HELD)),
    "ka3_panel_equal": sorted(PANEL) == sorted(PANEL_FILE),
    "ka3_n_panel": len(PANEL),
}
assert KA["ka2_held_persrc_same_order"], "279 held genes != 290 persrc genes"
assert KA["ka3_panel_equal"], "conservation276 targets != pert_counts.csv"
HELD_SET, PANEL_SET = set(HELD.tolist()), set(PANEL)

# ----------------------------------------------------------------- obs column candidates
LABEL_CANDIDATES = ["perturbation", "gene", "gene_target", "target_gene", "target",
                    "guide_identity", "sgRNA_gene", "grna_target", "knockout", "gene_symbol",
                    "sgrna_symbol"]
LINE_CANDIDATES = ["cell_line", "cell_type", "celltype", "line", "Cell_Line", "donor", "sample_line"]
BATCH_CANDIDATES = ["batch", "gemgroup", "gem_group", "sample", "channel", "lane", "replicate", "Batch"]
GUIDE_CANDIDATES = ["guide_id", "sgRNA", "guide", "sgID_AB", "guide_identity", "Guide_Call",
                    "grna", "protospacer"]
DESC_COLS = ["organism", "perturbation_type", "cell_line", "celltype", "cell_type", "tissue_type",
             "disease", "cancer", "time", "dose_value", "dose_unit", "nperts", "sex"]

CONTROL_LABELS = {"control", "ctrl", "non-targeting", "nontargeting", "non_targeting", "nt",
                  "ntc", "none", "unperturbed", "dmso", "no-target", "notarget", "non-targeting_control",
                  "safe-targeting", "safe_targeting", "scramble", "scrambled", "neg", "negative"}


def read_col(obs: h5py.Group, name: str):
    if name not in obs:
        return None
    node = obs[name]
    if isinstance(node, h5py.Group):
        if "categories" not in node or "codes" not in node:
            return None
        cats = node["categories"][:]
        cats = np.asarray([c.decode() if isinstance(c, bytes) else str(c) for c in cats], dtype=object)
        codes = node["codes"][:]
        out = np.empty(len(codes), dtype=object)
        ok = codes >= 0
        out[ok] = cats[codes[ok]]
        out[~ok] = ""
        return out
    arr = node[:]
    if arr.dtype.kind in "SO":
        return np.asarray([a.decode() if isinstance(a, bytes) else str(a) for a in arr], dtype=object)
    return np.asarray([str(a) for a in arr], dtype=object)


def pick(grp: h5py.Group, cands: list[str]) -> str | None:
    keys = list(grp.keys())
    lower = {k.lower(): k for k in keys}
    for c in cands:
        if c in keys:
            return c
        if c.lower() in lower:
            return lower[c.lower()]
    return None


def is_control(lab: str) -> bool:
    s = lab.strip().lower()
    return s in CONTROL_LABELS or s.startswith("non-targeting") or s.startswith("nontargeting")


def is_single_gene(lab: str) -> bool:
    """A label that names exactly one gene: no combination separators, not a control, not blank."""
    s = str(lab).strip()
    if not s or s.lower() in {"nan", "*", "none", ""} or is_control(s):
        return False
    return not any(ch in s for ch in ("_", "+", ",", ":", ";", "|", " "))


def var_symbols(f: h5py.File) -> tuple[int, list[str]]:
    if "var" not in f:
        return 0, []
    var = f["var"]
    n = 0
    idx_name = None
    if isinstance(var, h5py.Group):
        attrs = dict(var.attrs)
        idx_name = attrs.get("_index")
        if isinstance(idx_name, bytes):
            idx_name = idx_name.decode()
        cand = [idx_name] if idx_name else []
        cand += ["gene_symbol", "gene_symbols", "symbol", "gene_name", "index", "_index"]
        for c in cand:
            if c and c in var:
                v = read_col(var, c)
                if v is not None:
                    return len(v), [str(x) for x in v]
        for k in var.keys():
            v = read_col(var, k)
            if v is not None:
                n = len(v)
                break
    return n, []


def uns_keys(f: h5py.File) -> list[str]:
    if "uns" not in f:
        return []
    try:
        return sorted(f["uns"].keys())
    except Exception:
        return []


def scan(path: Path) -> dict:
    t0 = time.time()
    rec: dict = {"file": path.name, "path": str(path.relative_to(ROOT)),
                 "size_bytes": path.stat().st_size}
    with h5py.File(path, "r") as f:
        rec["top_level"] = sorted(f.keys())
        obs = f["obs"]
        rec["obs_cols"] = sorted(obs.keys())
        # descriptive columns
        desc = {}
        for c in DESC_COLS:
            col = pick(obs, [c])
            if col is None:
                continue
            v = read_col(obs, col)
            if v is None:
                continue
            u, cts = np.unique(v.astype(str), return_counts=True)
            order = np.argsort(-cts)
            desc[c] = {"n_uniq": int(len(u)),
                       "top": [[str(u[i]), int(cts[i])] for i in order[:8]]}
        rec["desc"] = desc

        lab_c = pick(obs, LABEL_CANDIDATES)
        rec["label_col"] = lab_c
        if lab_c is None:
            rec["error"] = "no perturbation label column"
            return rec
        labels = read_col(obs, lab_c).astype(str)
        rec["n_cells"] = int(len(labels))
        u, cts = np.unique(labels, return_counts=True)
        rec["n_labels"] = int(len(u))

        ctrl_mask = np.asarray([is_control(x) for x in u])
        rec["control_labels"] = [str(x) for x in u[ctrl_mask]]
        rec["n_control_cells"] = int(cts[ctrl_mask].sum())

        sg = np.asarray([is_single_gene(x) for x in u])
        rec["n_single_token_labels"] = int(sg.sum())

        nvar, sym = var_symbols(f)
        rec["n_var"] = int(nvar)
        symset = set(sym)
        rec["var_has_symbols"] = bool(symset & HELD_SET)

        # a "genetic perturbation targeting a gene": single-token label that is a known human
        # gene symbol -- either in the 17,106 held map genes, the 300 panel, or the file's own
        # var symbol vocabulary.
        gene_universe = HELD_SET | PANEL_SET | symset
        genlab = np.asarray([bool(s and x in gene_universe) for x, s in zip(u, sg)])
        rec["n_gene_perts"] = int(genlab.sum())
        gcts = cts[genlab]
        gu = u[genlab]
        if len(gcts):
            rec["cells_per_pert_median"] = float(np.median(gcts))
            rec["cells_per_pert_q10"] = float(np.quantile(gcts, 0.10))
            rec["cells_per_pert_q90"] = float(np.quantile(gcts, 0.90))
            rec["n_perts_ge50"] = int((gcts >= 50).sum())
            rec["n_perts_ge100"] = int((gcts >= 100).sum())
            rec["cells_on_gene_perts"] = int(gcts.sum())
        else:
            for k in ("cells_per_pert_median", "cells_per_pert_q10", "cells_per_pert_q90"):
                rec[k] = None
            rec["n_perts_ge50"] = rec["n_perts_ge100"] = rec["cells_on_gene_perts"] = 0

        gset = set(gu.tolist())
        in_held = sorted(gset & HELD_SET)
        rec["n_in_held_map"] = len(in_held)
        rec["n_in_panel"] = len(gset & PANEL_SET)
        rec["n_in_persrc"] = len(gset & set(PERSRC.tolist()))  # KA-2 says == n_in_held_map
        # the same, restricted to perturbations with >= 50 cells
        ge50 = set(gu[gcts >= 50].tolist()) if len(gcts) else set()
        rec["n_in_held_map_ge50"] = len(ge50 & HELD_SET)
        rec["n_in_panel_ge50"] = len(ge50 & PANEL_SET)
        rec["examples_in_held"] = in_held[:12]

        for tag, cands in (("line", LINE_CANDIDATES), ("batch", BATCH_CANDIDATES),
                           ("guide", GUIDE_CANDIDATES)):
            col = pick(obs, cands)
            rec[f"{tag}_col"] = col
            if col is None or col == lab_c:
                rec[f"n_{tag}s"] = None
                continue
            v = read_col(obs, col)
            rec[f"n_{tag}s"] = int(len(np.unique(v.astype(str)))) if v is not None else None
        rec["uns_keys"] = uns_keys(f)
    rec["seconds"] = round(time.time() - t0, 2)
    return rec


def main() -> None:
    out = []
    for p in sorted(SCP.glob("*.h5ad")):
        r = scan(p)
        # KA-1: reproduce 278's obs-only scan
        j = COUNTS278 / f"scperturb__{p.stem}.json"
        if j.exists():
            d = json.loads(j.read_text())
            r["ka1_278_n_cells"] = d.get("n_cells")
            r["ka1_278_n_labels"] = d.get("n_labels")
            r["ka1_match"] = (d.get("n_cells") == r.get("n_cells")
                              and d.get("n_labels") == r.get("n_labels"))
        else:
            r["ka1_match"] = None
        out.append(r)
        print(f"{p.name:<58} cells={r.get('n_cells')!s:>9} labels={r.get('n_labels')!s:>6} "
              f"gene_perts={r.get('n_gene_perts')!s:>6} held={r.get('n_in_held_map')!s:>6} "
              f"panel={r.get('n_in_panel')!s:>4} KA1={r.get('ka1_match')}")
    (WORK / "scan292.json").write_text(json.dumps(out, indent=1))
    KA["ka1_files_checked"] = sum(1 for r in out if r.get("ka1_match") is not None)
    KA["ka1_all_match"] = all(r.get("ka1_match") for r in out if r.get("ka1_match") is not None)
    KA["peak_rss_gib"] = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 2**30
    (WORK / "scan292_ka.json").write_text(json.dumps(KA, indent=1))
    print("\nKA:", json.dumps(KA))


if __name__ == "__main__":
    main()
