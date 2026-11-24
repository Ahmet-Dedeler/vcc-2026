"""233 Phase 0/2 — coverage of every source under the three panels, and the CARRYING-ROW
counts each PLOLO fold actually gets.

Two jobs.

1. KNOWN-ANSWER.  Re-derive the per-source panel coverage independently and check it against
   `docs/DEF-RUNBOOK.md` §5.-2's published table (225's numbers).  Any disagreement means my
   panel files or my source registry are not the ones those numbers were measured on, and
   everything downstream would be a ratio against the wrong control.

2. THE ATTRIBUTION SUBSTRATE.  §5's coverage table is stated per SOURCE against the whole
   300.  What a PLOLO fold sees is different and is what actually carries the score: for the
   fold that holds out donor X, each evaluated target has 0, 1 or 2 SOURCE carriers among the
   remaining two, and the CD4/H1 residuals carry their own row sets.  Falsifier 4 of the
   mandate: report CARRYING rows per source per panel, never the all-rows median.

    uv run python census233.py
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(ROOT / "experiments/218-independent-donors"))
sys.path.insert(0, str(HERE))

import panels233 as P  # noqa: E402
from sources218 import SOURCES  # noqa: E402

# `docs/DEF-RUNBOOK.md` §5.-2, the pessimistic column, and §5's realistic column.
KA_PESS = {"orion_hek293t": 293, "orion_hct116": 291, "k562_gwps": 150,
           "gse295214_h1": 34, "h1_2025": 3, "nadig_hepg2": 33, "nadig_jurkat": 33,
           "rpe1": 32, "k562_essential": 30, "kolf2": 0, "feng_ipsc": 0}
KA_TODAY = {"orion_hek293t": 294, "orion_hct116": 290, "k562_gwps": 269,
            "gse295214_h1": 64, "h1_2025": 20, "nadig_hepg2": 0, "nadig_jurkat": 0,
            "rpe1": 0, "k562_essential": 0, "kolf2": 273, "feng_ipsc": 182}
KA_REAL = {"k562_gwps": 300, "orion_hct116": 300, "orion_hek293t": 300}

ORDER = ("k562_gwps", "orion_hct116", "orion_hek293t", "kolf2", "feng_ipsc",
         "gse295214_h1", "h1_2025", "nadig_hepg2", "nadig_jurkat", "rpe1",
         "k562_essential", "telohaec", "myllia_a549", "frangieh_melanoma", "perturbmatch")


def targets_of(name: str) -> set[str] | None:
    meta = SOURCES.get(name)
    if meta is None:
        return None
    p = ROOT / meta["path"] / "perts.txt"
    if not p.exists():
        return None
    return {x.strip() for x in p.read_text().split("\n") if x.strip()}


def main() -> None:
    panels = {k: set(P.read_panel(k)) for k in ("abc", "real", "pess")}
    assert not (panels["abc"] & panels["real"]), "realistic draw overlaps today's panel"
    assert not (panels["abc"] & panels["pess"]), "pessimistic draw overlaps today's panel"
    print("FALSIFIER 1 — panel disjointness, asserted not assumed:")
    for k in ("real", "pess"):
        print(f"  |abc n {k}| = {len(panels['abc'] & panels[k])}   "
              f"|{k}| = {len(panels[k])}   |abc u {k}| = {len(panels['abc'] | panels[k])}")
    print(f"  |real n pess| = {len(panels['real'] & panels['pess'])}\n")

    cov, fails = {}, []
    print(f"{'source':22} {'total':>7} {'abc':>6} {'real':>6} {'pess':>6}   KA")
    for name in ORDER:
        t = targets_of(name)
        if t is None:
            print(f"{name:22} {'-- not on disk --':>28}")
            continue
        c = {k: len(t & panels[k]) for k in panels}
        cov[name] = dict(total=len(t), **c)
        marks = []
        for tbl, key in ((KA_TODAY, "abc"), (KA_REAL, "real"), (KA_PESS, "pess")):
            if name in tbl:
                ok = tbl[name] == c[key]
                marks.append(f"{key}={'PASS' if ok else f'FAIL({tbl[name]})'}")
                if not ok:
                    fails.append((name, key, tbl[name], c[key]))
        print(f"{name:22} {len(t):7d} {c['abc']:6d} {c['real']:6d} {c['pess']:6d}   "
              + " ".join(marks))

    # fold-level carriers: for each fold, how many EVALUATED targets each source carries
    print("\nPer-fold carrier structure (evaluated targets = panel n held-donor coverage):")
    fold_rows = {}
    for panel_key, panel in panels.items():
        for held, sources in P.FOLD_SOURCES.items():
            held_t = targets_of(held)
            ev = panel & held_t
            per = {s: len(ev & targets_of(s)) for s in sources}
            both = len(ev & targets_of(sources[0]) & targets_of(sources[1]))
            one = sum(1 for g in ev
                      if (g in targets_of(sources[0])) ^ (g in targets_of(sources[1])))
            none = len(ev) - both - one
            fold_rows[f"{panel_key}|{held}"] = dict(
                eval_targets=len(ev), per_source=per,
                carriers_2=both, carriers_1=one, carriers_0=none)
            print(f"  {panel_key:5} held {held:14} eval {len(ev):3d}  "
                  f"sources {per}  carriers 2/1/0 = {both}/{one}/{none}")

    out = dict(coverage=cov, folds=fold_rows,
               panels={k: sorted(v) for k, v in panels.items()},
               ka_failures=fails)
    (P.WORK / "census233.json").write_text(json.dumps(out, indent=2))
    print(f"\nknown-answer failures: {len(fails)}")
    if fails:
        for f in fails:
            print("  FAIL", f)
        raise SystemExit(1)
    print("-> work/census233.json")


if __name__ == "__main__":
    main()
