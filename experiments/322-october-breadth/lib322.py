"""322 — shared setup.  233's instrument, imported not rewritten.

WHAT IS NEW HERE.  Nothing about the scoring path: `panels233`, `fold233`, `score233` and
`oct233`'s `run()` shape are 233's, and every operator / CD4 bank / fold this experiment
reads was written by 233 and is still on disk (`experiments/233-october-retention/work/`,
12 GB).  The only new thing is (a) a wider SOURCE REGISTRY (225's `_reg` trick: register
extra donor-shaped tables into `sources218.SOURCES` at runtime, read-only) and (b) two
new source TRANSFORMS -- coverage-restriction and coverage-matched random re-assignment --
which are the mandate's two controls.

THE TWO SENSES OF "BROAD" THAT THIS EXPERIMENT SEPARATES.  A source has

  * PANEL coverage   -- how many of the 300 CURRENT targets it carries; and
  * SOURCE breadth   -- how many perturbations it contains AT ALL.

They are not the same axis and the mandate's candidate is the case where they disagree
most: `kolf2` is 273/300 on panel coverage and 273 perturbations TOTAL.  Which of the two
predicts October's coverage is measured in `census322.py`.

NOTHING IS WRITTEN OUTSIDE `experiments/322-october-breadth/`.  233's work directory is
opened read-only; `plolo.WORK` is re-pointed at 233's SHARED cache exactly as `score233`
does (that cache is panel-independent: control blocks and challenge b_ctrl).
"""
from __future__ import annotations

import os
import sys
from pathlib import Path

os.environ.setdefault("OPENBLAS_NUM_THREADS", "6")
os.environ.setdefault("OMP_NUM_THREADS", "6")

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
WORK = HERE / "work"
WORK.mkdir(parents=True, exist_ok=True)

E233 = ROOT / "experiments/233-october-retention"
assert E233.is_dir(), E233

# 233's own modules come first on the path so `import panels233 / score233 / fold233 / oct233`
# resolve to 233's files, unchanged.
for p in (str(E233),
          str(ROOT / "experiments/162-retrodicting-validator"),
          str(ROOT / "experiments/218-independent-donors"),
          str(ROOT / "experiments/223-silent-zero-audit"),
          str(ROOT / "experiments/228-nmae-anatomy"),
          str(ROOT / "experiments/188-real-de"),
          str(ROOT / "experiments/197-mse-inversion")):
    if p not in sys.path:
        sys.path.insert(0, p)

_T223 = "experiments/223-silent-zero-audit/work/tables"
_T225 = "experiments/225-source-sweep/work/tables"

# --------------------------------------------------------------- the extra registrations
# 225's `harness225._reg`, reproduced here so 322 does not have to import 225's harness
# (which re-points `plolo.WORK` at 225's directory).  Paths and roles copied verbatim from
# `experiments/225-source-sweep/harness225.py`; the tables themselves are read-only.
EXTRA_SOURCES: dict[str, dict] = {
    "dld1_lowmoi": dict(
        path=f"{_T223}/dld1_lowmoi", lineage="epithelial (colorectal, DLD-1)",
        lab="116's acquisition", semantics="author log2FC, split-mean",
        role="223's live source"),
    "jiang_5line_noK562": dict(
        path=f"{_T225}/jiang_5line_noK562", lineage="mixed (5 lines, K562 removed)",
        lab="Jiang 2024", semantics="log-pseudobulk delta", role="225 candidate"),
    "igvf_cardiomyocyte": dict(
        path=f"{_T225}/igvf_cardiomyocyte", lineage="cardiomyocyte (WTC11-derived)",
        lab="IGVF / Hon", semantics="published log2FC, significant hits only",
        role="225 candidate"),
    "cmap_xpr": dict(
        path=f"{_T223}/cmap_xpr", lineage="mixed cancer lines", lab="LINCS 2020",
        semantics="MODZ consensus z", role="223 death"),
    "cmap_sh_cgs": dict(
        path=f"{_T223}/cmap_sh_cgs", lineage="mixed cancer lines", lab="LINCS 2020",
        semantics="MODZ consensus z", role="223 death"),
    "altair_a549": dict(
        path=f"{_T223}/altair_a549", lineage="epithelial (lung adeno)", lab="Altair",
        semantics="official OPS log2FC", role="223 death"),
}


def register_extra() -> list[str]:
    """Add the extra tables to 218's registry.  Idempotent; returns the names added."""
    import sources218
    sources218.SOURCES = dict(sources218.SOURCES)
    added = []
    for name, meta in EXTRA_SOURCES.items():
        if name in sources218.SOURCES:
            continue
        p = ROOT / meta["path"]
        if not (p / "perts.txt").exists():
            continue
        sources218.SOURCES[name] = dict(meta)
        added.append(name)
    return added


def bank_root(name: str) -> Path:
    import sources218
    return (ROOT / sources218.SOURCES[name]["path"]).resolve()


def source_perts(name: str) -> list[str]:
    return [ln.strip() for ln in (bank_root(name) / "perts.txt").read_text().splitlines()
            if ln.strip()]


# --------------------------------------------------------------- the two control transforms
def restrict_residual(res, keep_rows: np.ndarray):
    """A source residual with every row OUTSIDE `keep_rows` zeroed and de-voted.

    This is the *coverage held fixed* control: the source's own biology on a chosen subset
    of the targets it carries, so "more targets covered" and "different biology" separate.
    Same object type `oct233.deranged` returns, so it drops straight into `run()`.
    """
    from plolo import E062
    n = np.asarray(res.mean_effect).shape[0]
    mask = np.zeros(n, dtype=bool)
    mask[np.asarray(keep_rows, dtype=np.int64)] = True
    eff = np.asarray(res.mean_effect, dtype=np.float64).copy()
    votes = np.asarray(res.gene_votes).copy()
    cnt = np.asarray(res.target_condition_count).copy()
    eff[~mask] = 0.0
    votes[~mask] = 0
    cnt[~mask] = 0
    return E062.CD4Prepared(mode=res.mode + "|restricted", mean_effect=eff,
                            gene_votes=votes, target_condition_count=cnt,
                            condition_details=[])


def carrying_rows(res) -> np.ndarray:
    return np.flatnonzero(np.asarray(res.target_condition_count) > 0)


def coverage_matched_random(res, n_keep: int, rng) -> object:
    """Keep a RANDOM `n_keep` of the source's carrying rows.

    The mandate's *coverage-matched random source assignment*: it holds the source's
    identity and its per-row biology fixed and varies only HOW MANY targets it reaches, so
    a contribution that tracks coverage and one that tracks identity come apart.
    """
    rows = carrying_rows(res)
    if n_keep >= rows.size:
        return restrict_residual(res, rows)
    keep = rng.choice(rows, size=int(n_keep), replace=False)
    return restrict_residual(res, keep)


def memory_ok(min_free_pct: float = 40.0, min_free_inactive_gib: float = 20.0) -> tuple[bool, str]:
    """233's gate (`whiten233.memory_ok`), same thresholds, restated here.

    `vm.swapusage free` is the unused tail of a DYNAMICALLY sized swap file, not a pressure
    reading, so the gate is `memory_pressure` free % AND vm_stat free+inactive.
    """
    import subprocess
    mp = subprocess.run(["memory_pressure"], capture_output=True, text=True).stdout
    pct = None
    for ln in mp.splitlines():
        if "free percentage" in ln:
            pct = float(ln.rsplit(":", 1)[1].strip().rstrip("%"))
    vs = subprocess.run(["vm_stat"], capture_output=True, text=True).stdout
    page = 16384
    vals = {}
    for ln in vs.splitlines():
        if ":" in ln and "Pages" in ln:
            k, v = ln.split(":", 1)
            try:
                vals[k.strip()] = int(v.strip().rstrip("."))
            except ValueError:
                pass
    fi = (vals.get("Pages free", 0) + vals.get("Pages inactive", 0)) * page / 2**30
    ok = (pct is not None and pct >= min_free_pct) and fi >= min_free_inactive_gib
    return ok, f"memory_pressure free {pct}% ; free+inactive {fi:.1f} GiB"
