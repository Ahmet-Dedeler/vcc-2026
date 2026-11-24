"""250 — point 233's instrument (fold233 / score233, verbatim) at the five panels, writing
NOTHING outside this experiment.

238's `bind238.py` pattern: rebind the attributes `fold233` / `score233` read at CALL time on
their shared `panels233` module object, BEFORE `score233` is imported (its `SHARED` path is
computed from `panels233.WORK` at import).

  abc   `data/challenge`                          233's folds / ops / CD4 banks
  real  199 `work/def-bundle-r`                   233's
  pess  199 `work/def-bundle`                     233's
  efr   238 `work/ef-bundle-r` (essential-free)   238's folds / ops / CD4 banks
  efp   238 `work/ef-bundle`   (essential-free)   238's

Every artefact is READ from 233's / 238's directories; every cache write (`plolo.WORK`)
lands in `experiments/250-deflation-october/work/`.
"""
from __future__ import annotations

import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
R233 = ROOT / "experiments/233-october-retention"
R238 = ROOT / "experiments/238-essential-free-panel"
WORK = HERE / "work"
WORK.mkdir(parents=True, exist_ok=True)

for p in (str(HERE), str(R233),
          str(ROOT / "experiments/162-retrodicting-validator"),
          str(ROOT / "experiments/188-magnitude-calibration"),
          str(ROOT / "experiments/203-reopen-125"),
          str(ROOT / "experiments/218-independent-donors"),
          str(ROOT / "experiments/228-theta-gain-plolo")):
    if p not in sys.path:
        sys.path.insert(0, p)

if "score233" in sys.modules or "fold233" in sys.modules:
    raise SystemExit("bind250 must be imported BEFORE score233/fold233")

import numpy as np  # noqa: E402
import panels233 as P233  # noqa: E402

BUNDLES = {
    "abc": ROOT / "data/challenge",
    "real": ROOT / "experiments/199-def-runbook/work/def-bundle-r",
    "pess": ROOT / "experiments/199-def-runbook/work/def-bundle",
    "efr": R238 / "work/ef-bundle-r",
    "efp": R238 / "work/ef-bundle",
}
OWN238 = ("efr", "efp")
DISJOINT = ("real", "pess", "efp", "efr")
FOLD_SOURCES = {
    "k562_gwps": ("orion_hct116", "orion_hek293t"),
    "orion_hct116": ("k562_gwps", "orion_hek293t"),
    "orion_hek293t": ("k562_gwps", "orion_hct116"),
    "kolf2": ("k562_gwps", "orion_hct116", "orion_hek293t"),   # abc only (238 §4)
}
CANCER_FOLDS = ("k562_gwps", "orion_hct116", "orion_hek293t")


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


def _root(panel: str) -> Path:
    return (R238 / "work") if panel in OWN238 else (R233 / "work")


def op_path(panel: str, held: str) -> Path:
    if held == "kolf2":
        assert panel == "abc", "kolf2 has 0 targets on every disjoint panel (238 §4)"
        return R238 / "work/ops" / f"{panel}__production3.npz"
    return _root(panel) / "ops" / f"{panel}__held_{held}.npz"


def cd4_root(panel: str) -> Path:
    return _root(panel) / "cd4" / panel / "cd4_compact"


def fold_cache(panel: str, held: str) -> Path:
    if held == "kolf2":                       # 238 built PLOLO's fourth fold on abc (233 dropped it)
        return R238 / "work" / "folds" / f"{panel}__{held}.npz"
    return _root(panel) / "folds" / f"{panel}__{held}.npz"


def bind(panel: str) -> dict:
    """`panels233.bind`, but every cache path is under 250's `work/`."""
    import plolo
    panel_t = read_panel(panel)
    plolo.PANEL = panel_t
    plolo.PANEL_SET = set(panel_t)
    plolo.WORK = WORK / "plolo" / panel
    plolo.WORK.mkdir(parents=True, exist_ok=True)
    plolo.CD4_ROOT = cd4_root(panel)
    for held in FOLD_SOURCES:
        if held == "kolf2" and panel != "abc":
            plolo.FOLDS.pop("kolf2", None)
            continue
        plolo.FOLDS.setdefault(held, {})
        plolo.FOLDS[held]["op"] = op_path(panel, held)
        plolo.FOLDS[held]["op_kind"] = "refit233"
        plolo.FOLDS[held]["sources"] = FOLD_SOURCES[held]
    return dict(panel=panel, bundle=str(BUNDLES[panel].relative_to(ROOT)),
                n_targets=len(panel_t), cd4_root=str(cd4_root(panel).relative_to(ROOT)))


P233.BUNDLES = BUNDLES
P233.WORK = WORK
P233.read_panel = read_panel
P233.read_genes = read_genes
P233.op_path = op_path
P233.cd4_root = cd4_root
P233.bind = bind
P233.FOLD_SOURCES = FOLD_SOURCES
# `panels233.axis_sha256` / `load_frozen_operator_233` read `read_genes` from their own module
# globals, which are now rebound; CANCER_FOLDS is unchanged.

# panel disjointness, asserted at import (233 falsifier 1)
_abc = set(read_panel("abc"))
for _p in DISJOINT:
    assert not (_abc & set(read_panel(_p))), f"panel {_p} intersects abc"
