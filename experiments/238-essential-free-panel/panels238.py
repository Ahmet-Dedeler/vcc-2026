"""238 — the FIVE panels, and the machinery that points 162's PLOLO instrument at any of them.

`panels233.py` made "the panel" a variable across three panels.  This adds the two the mandate
asks for and keeps the other three addressable so every ratio is against the same instrument:

  abc   `data/challenge`                     today's real panel     0 essential genes
  real  199 `work/def-bundle-r`              K562 n Orion n CD4    34 essential genes
  pess  199 `work/def-bundle`                any one donor         30 essential genes
  efr   238 `work/ef-bundle-r`   NEW         real pool  MINUS the essential list
  efp   238 `work/ef-bundle`     NEW         pess pool  MINUS the essential list

The two new bundles are 199's own draws with EXACTLY ONE change -- the eligible pool has the
common-essential list removed -- and the same seeds (199 / 1990), so nothing else about the
sampling model moves.  They match Arc's one observed panel on BOTH properties that matter here:
disjoint from A/B/C, and ~0 common-essential targets.

Operators and CD4 banks: `abc`/`real`/`pess` reuse 233's (already fitted, axis-verified); the
two new panels get their own, fitted into 238's `work/`.  Every load goes through the axis
assertion, because the operator's gene axis is `union(donor coverage) - panel` and it fails
SILENTLY on the wrong panel (AGENTS.md §8.3; 233 §2 priced it at -0.021 overall on `pds`).
"""
from __future__ import annotations

import hashlib
import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
WORK = HERE / "work"
WORK.mkdir(parents=True, exist_ok=True)
R233 = ROOT / "experiments/233-october-retention"

BUNDLES = {
    "abc": ROOT / "data/challenge",
    "real": ROOT / "experiments/199-def-runbook/work/def-bundle-r",
    "pess": ROOT / "experiments/199-def-runbook/work/def-bundle",
    "efr": WORK / "ef-bundle-r",
    "efp": WORK / "ef-bundle",
}
OWN = ("efr", "efp")                      # the panels whose artefacts 238 makes
FOLD_SOURCES = {
    "k562_gwps": ("orion_hct116", "orion_hek293t"),
    "orion_hct116": ("k562_gwps", "orion_hek293t"),
    "orion_hek293t": ("k562_gwps", "orion_hct116"),
    # PLOLO's fourth fold.  `kolf2` is an iPSC line and a PANEL-FITTED library: 273 of its 273
    # targets are on today's 300 and **0 of them are on any disjoint panel** (225 §8.3;
    # re-measured here on all four).  So this fold is evaluable on `abc` ONLY, and it is
    # reported separately -- on any October-shaped panel it has zero targets and cannot run.
    "kolf2": ("k562_gwps", "orion_hct116", "orion_hek293t"),
}
CANCER_FOLDS = ("k562_gwps", "orion_hct116", "orion_hek293t")
ALL_FOLDS = tuple(FOLD_SOURCES)


def read_panel(key: str) -> tuple[str, ...]:
    path = BUNDLES[key] / "pert_counts.csv"
    lines = path.read_text().splitlines()
    assert lines[0].strip() == "target_gene", f"{path}: unexpected header {lines[0]!r}"
    panel = tuple(ln.split(",")[0].strip() for ln in lines[1:] if ln.strip())
    assert len(panel) == 300 and len(set(panel)) == 300, f"{path}: {len(panel)} targets"
    return panel


def read_genes(key: str) -> tuple[str, ...]:
    lines = (BUNDLES[key] / "gene_names.csv").read_text().splitlines()
    assert lines[0].strip() == "gene_name"
    return tuple(ln.strip() for ln in lines[1:] if ln.strip())


def axis_sha256(panel_set: set[str], source_banks) -> tuple[str, int]:
    """sha256 of the operator's gene axis = `union(donor coverage) - panel`, official order.

    Identical construction to `panels233.axis_sha256` (which is 064's `operator_columns`).
    """
    official = read_genes("abc")
    donor_union: set[str] = set()
    for bank in source_banks:
        donor_union |= set(bank.genes)
    wcols = np.asarray([i for i, g in enumerate(official)
                        if g in donor_union and g not in panel_set], dtype=np.int64)
    return hashlib.sha256(wcols.tobytes()).hexdigest(), int(len(wcols))


def op_path(panel: str, held: str) -> Path:
    """The operator this (panel, fold) must use.  238 fits its own; 233's are reused as-is.

    The `kolf2` fold's basis is all THREE production donors, which is exactly the production
    operator's basis -- so it uses 238's own re-fit of it, the one whose known-answer row
    reproduces the shipped `e2f65ecf…` (`whiten238.py --jobs ka`).
    """
    if held == "kolf2":
        return WORK / "ops" / f"{panel}__production3.npz"
    root = WORK if panel in OWN else (R233 / "work")
    return root / "ops" / f"{panel}__held_{held}.npz"


def cd4_root(panel: str) -> Path:
    root = WORK if panel in OWN else (R233 / "work")
    return root / "cd4" / panel / "cd4_compact"


def bind(panel: str) -> dict:
    """Point `plolo` (and everything importing it) at `panel`.  Returns an audit dict."""
    sys.path.insert(0, str(ROOT / "experiments/162-retrodicting-validator"))
    import plolo

    panel_t = read_panel(panel)
    plolo.PANEL = panel_t
    plolo.PANEL_SET = set(panel_t)
    plolo.WORK = WORK / "plolo" / panel
    plolo.WORK.mkdir(parents=True, exist_ok=True)
    plolo.CD4_ROOT = cd4_root(panel)
    for held in ALL_FOLDS:
        if held not in plolo.FOLDS:
            continue
        plolo.FOLDS[held]["op"] = op_path(panel, held)
        plolo.FOLDS[held]["op_kind"] = "refit233"
        plolo.FOLDS[held]["sources"] = FOLD_SOURCES[held]
    return dict(panel=panel, bundle=str(BUNDLES[panel]).replace(str(ROOT) + "/", ""),
                n_targets=len(panel_t), cd4_root=str(cd4_root(panel)))


def load_frozen_operator_238(path, source_banks, official, panel_set, *, expected_sha256):
    """`panels233.load_frozen_operator_233`, verbatim: verify the AXIS, not just its length."""
    sys.path.insert(0, str(ROOT / "experiments/162-retrodicting-validator"))
    import plolo
    from plolo import E062

    path = Path(path)
    meta = json.loads(path.with_suffix(".json").read_text())
    sha, n = axis_sha256(set(panel_set), source_banks)
    actual = E062.sha256_file(path)
    if expected_sha256 and actual != expected_sha256:
        raise ValueError(f"{path}: operator sha256 {actual[:12]} != {expected_sha256[:12]}")
    if int(meta["n_cols"]) != n:
        raise ValueError(f"{path}: n_cols {meta['n_cols']} != {n}")
    if meta.get("axis_sha256") != sha:
        raise ValueError(
            f"{path}: OPERATOR AXIS MISMATCH -- fitted for panel {meta.get('panel')!r} "
            f"(axis {str(meta.get('axis_sha256'))[:12]}), asked for axis {sha[:12]}")
    if meta.get("donors") != [b.name for b in source_banks]:
        raise ValueError(f"{path}: donors {meta.get('donors')} != "
                         f"{[b.name for b in source_banks]}")
    with np.load(path) as z:
        evals = np.asarray(z["evals"], dtype=np.float64)
        evecs = np.asarray(z["evecs"], dtype=np.float32)
    official_pos = {g: i for i, g in enumerate(official)}
    donor_union: set[str] = set()
    for b in source_banks:
        donor_union |= set(b.genes)
    wcols = np.asarray(sorted(official_pos[g] for g in donor_union
                              if g in official_pos and g not in panel_set), dtype=np.int64)
    if evecs.shape != (len(wcols), len(wcols)):
        raise ValueError(f"{path}: eigenvector shape {evecs.shape} invalid")
    return plolo.E064.BASE.Operator(evals, evecs, meta, path, path.with_suffix(".json")), wcols
