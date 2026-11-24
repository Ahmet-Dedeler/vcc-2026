"""233 — the three panels, and the machinery that points 162's PLOLO instrument at any of them.

THE QUESTION.  On 2026-10-22 the 300-target panel is replaced by one DISJOINT from today's.
Everything we have measured is measured on today's panel.  This module makes "the panel" a
variable so the same instrument can be run on all three and the RATIO read off.

THE THREE PANELS (all 300 targets; the two synthetic ones are 199's, reused verbatim so no
new sampling model enters -- `docs/DEF-RUNBOOK.md` §5):

  abc   `data/challenge`                      today's real panel (the control)
  real  199 `work/def-bundle-r`               realistic draw: K562 n Orion n CD4
  pess  199 `work/def-bundle`                 pessimistic draw: any one production donor

Both synthetic bundles carry the SAME control h5ads (A/B/C relabelled D/E/F, byte-identical)
and the SAME `gene_names.csv`.  Only `pert_counts.csv` differs.  So the contrast isolates
exactly one variable -- WHICH 300 GENES ARE THE PANEL -- and holds the cell lines, the gene
axis and the control cells fixed.  What it therefore does NOT measure is the change of
CONTEXT (A/B/C -> D/E/F cell lines); that half is unmeasurable before the day.  PLOLO's own
folds hold out a whole donor LINE as truth, so cross-line transfer is still exercised.

THE AXIS ASSERTION (AGENTS.md §8.3, the Oct-22 trap).  The whitening operator's gene axis is
`union(donor coverage) - panel`.  Production's operator metadata records `n_cols` but NOT a
column hash, so pointing the A/B/C operator at a new panel passes every check the builder
runs while 98.8% of coordinates are paired with the wrong gene.  `axis_sha256` below is the
assertion that fires; `trap233.py` demonstrates it firing, and demonstrates what passes
without it.
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

BUNDLES = {
    "abc": ROOT / "data/challenge",
    "real": ROOT / "experiments/199-def-runbook/work/def-bundle-r",
    "pess": ROOT / "experiments/199-def-runbook/work/def-bundle",
}
FOLD_SOURCES = {
    "k562_gwps": ("orion_hct116", "orion_hek293t"),
    "orion_hct116": ("k562_gwps", "orion_hek293t"),
    "orion_hek293t": ("k562_gwps", "orion_hct116"),
}
CANCER_FOLDS = tuple(FOLD_SOURCES)


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
    """sha256 of the operator's gene axis = `union(donor coverage) - panel`, on the official axis.

    Identical construction to `experiments/064-pm-cd4-composition/run.py:operator_columns`.
    """
    official = read_genes("abc")
    donor_union: set[str] = set()
    for bank in source_banks:
        donor_union |= set(bank.genes)
    wcols = np.asarray([i for i, g in enumerate(official)
                        if g in donor_union and g not in panel_set], dtype=np.int64)
    return hashlib.sha256(wcols.tobytes()).hexdigest(), int(len(wcols))


def op_path(panel: str, held: str) -> Path:
    """The operator this (panel, fold) must use.  Every one is re-fitted for its own panel."""
    return WORK / "ops" / f"{panel}__held_{held}.npz"


def cd4_root(panel: str) -> Path:
    return WORK / "cd4" / panel / "cd4_compact"


def bind(panel: str) -> dict:
    """Point `plolo` (and everything importing it) at `panel`.  Returns an audit dict.

    Patches the module GLOBALS `plolo` reads at call time -- `PANEL`, `PANEL_SET`, `WORK`
    (so each panel gets its own fold cache), `CD4_ROOT`, and each fold's operator path.
    """
    sys.path.insert(0, str(ROOT / "experiments/162-retrodicting-validator"))
    import plolo

    panel_t = read_panel(panel)
    plolo.PANEL = panel_t
    plolo.PANEL_SET = set(panel_t)
    plolo.WORK = WORK / "plolo" / panel
    plolo.WORK.mkdir(parents=True, exist_ok=True)
    plolo.CD4_ROOT = cd4_root(panel)
    for held in CANCER_FOLDS:
        plolo.FOLDS[held]["op"] = op_path(panel, held)
        plolo.FOLDS[held]["op_kind"] = "refit233"
        plolo.FOLDS[held]["sources"] = FOLD_SOURCES[held]
    plolo.FOLDS.pop("kolf2", None)          # iPSC; the mandate's three folds are the cancer ones
    return dict(panel=panel, bundle=str(BUNDLES[panel].relative_to(ROOT)),
                n_targets=len(panel_t), cd4_root=str(plolo.CD4_ROOT.relative_to(ROOT)))


def load_frozen_operator_233(path, source_banks, official, panel_set, *, expected_sha256):
    """Drop-in for `E064.load_frozen_operator` that verifies the AXIS, not just its length.

    Production's own operator metadata has no `wcols_sha256`; 062's held operators do.  Ours
    always does, and this refuses to load an operator whose axis was fitted for another panel
    even when `n_cols` coincides -- which it does, exactly, between A/B/C and both synthetic
    panels (17,942 columns all three times).  That coincidence is the trap.
    """
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
