"""229 — the October table: what each registered source covers when the panel is REPLACED.

On 2026-10-22 the 300-target panel is replaced with 300 targets disjoint from A/B/C's, and
only the last entry after that counts.  Every source's value in this arm is a coverage
statistic, so it is projectable in advance against the two synthetic draws 199 built:

  `def-bundle-r`  the REALISTIC draw -- targets in K562 ∩ Orion ∩ CD4
  `def-bundle`    the PESSIMISTIC draw -- targets in ANY donor

Counted the way the residual builder counts: a source carries a target iff the target is a row
of the bank's own `perts.txt`.  (The builder additionally requires >= 1 mapped non-panel gene
column, which every bank here satisfies for every row it has.)  Also reported: the carrying-row
sizing consequence -- `solve_carrying_row_fraction` is coverage-blind by construction, so a
coverage drop does NOT silently zero the residual the way production's all-rows median does
(that is the whole point of 222's rule), but the CD4 residual beside it still over-doses by
~1.24x at the pessimistic draw (222 §10).

    uv run python octo229.py
"""
from __future__ import annotations

import csv
import json
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
WORK = HERE / "work"

PANELS = {
    "A/B/C (current)": ROOT / "data/challenge/pert_counts.csv",
    "realistic draw": ROOT / "experiments/199-def-runbook/work/def-bundle-r/pert_counts.csv",
    "pessimistic draw": ROOT / "experiments/199-def-runbook/work/def-bundle/pert_counts.csv",
}


def read_panel(path: Path) -> set[str]:
    with path.open() as fh:
        return {row["target_gene"].strip() for row in csv.DictReader(fh)
                if row.get("target_gene", "").strip()}


def main() -> None:
    sys.path.insert(0, str(HERE))
    import price229 as P  # noqa: E402  (imports the builder machinery, not the folds)

    MOD, builder_path = P.load_builder()
    panels = {name: read_panel(path) for name, path in PANELS.items()}
    for name, panel in panels.items():
        print(f"{name:<18} {len(panel)} targets")
    overlap = panels["A/B/C (current)"] & (panels["realistic draw"]
                                           | panels["pessimistic draw"])
    print(f"A/B/C ∩ (either draw): {len(overlap)}  (the draws are built disjoint)")

    out = {"panels": {k: len(v) for k, v in panels.items()},
           "abc_overlap_with_draws": len(overlap), "sources": []}
    for name in MOD.UNION_SOURCE_REGISTRY:
        path = MOD.UNION_SOURCE_REGISTRY[name]["path"]
        perts = {line.strip() for line in (path / "perts.txt").read_text().splitlines()
                 if line.strip()}
        row = {"source": name, "path": str(path.relative_to(ROOT)),
               "bank_rows": len(perts)}
        for pname, panel in panels.items():
            row[pname] = len(perts & panel)
        out["sources"].append(row)

    hdr = f"{'source':<22}{'rows':>8}" + "".join(f"{k:>20}" for k in panels)
    print("\n" + hdr)
    for row in out["sources"]:
        print(f"{row['source']:<22}{row['bank_rows']:>8}"
              + "".join(f"{row[k]:>20}" for k in panels))

    # the union's carrying count, which is what the sizing rule reads
    for tag, names in (("h1pair (incumbent)", MOD.UNION_SOURCES_DEFAULT),
                       ("stack4", ("dld1_lowmoi", "jiang_5line_noK562",
                                   "igvf_cardiomyocyte", "myllia_a549")),
                       ("stack6", tuple(MOD.UNION_SOURCE_REGISTRY))):
        u = {}
        for pname, panel in panels.items():
            carried = set()
            for n in names:
                p = MOD.UNION_SOURCE_REGISTRY[n]["path"]
                perts = {l.strip() for l in (p / "perts.txt").read_text().splitlines()
                         if l.strip()}
                carried |= (perts & panel)
            u[pname] = len(carried)
        out.setdefault("unions", {})[tag] = u
        print(f"{'UNION ' + tag:<22}{'':>8}" + "".join(f"{u[k]:>20}" for k in panels))

    (WORK / "octo229.json").write_text(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()
