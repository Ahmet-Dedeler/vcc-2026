"""235 addendum — redo item B's coverage on a panel that is disjoint AND ESSENTIAL-FREE.

`AGENTS.md` §8.1g landed (commit `eed4e4a`) while this experiment was running: **the real A/B/C
panel contains ZERO essential genes**, for an assay reason that should hold again in October, and
`DEF-RUNBOOK` §5's two synthetic bundles carry **34 and 30** of them because they are drawn from
donor coverage.  The rule it sets is that no synthetic-panel number may be quoted without stating
the panel's essential-gene count, and that a panel used to predict October must have **both**
properties.

235's item-B numbers were measured on those bundles.  The *structural* claim -- the carrying-row
rule sizes every source the registry can select, and the shipped all-rows rule returns exactly zero
for all of them -- does not depend on which targets the panel holds.  The *coverage* numbers do.
So this rebuilds the measurement on a third panel that has both properties:

    pool  = targets carried by any of the three production donors, ON the official gene axis
    panel = 300 drawn from  pool - (A/B/C panel) - (common-essential list),  seed 0

and reports, for every registered source, coverage and both sizing rules -- alongside the same
quantities on A/B/C (essential 0) and on 199's pessimistic bundle (essential 30).

Known-answer rows: the panel is 300 unique symbols, its intersection with A/B/C's panel is 0, its
intersection with the essential list is 0, and the A/B/C column reproduces `work/carry235.json`.

    uv run python experiments/235-october-hardening/essfree235.py
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(HERE))

import carry235 as C  # noqa: E402
import make_submission_pm_cd4 as M  # noqa: E402

ESSENTIAL = ROOT / "data/donors/k562_essential/perts.txt"
DONORS = ("k562_gwps", "orion_hct116", "orion_hek293t")
OUT = HERE / "work" / "essfree235.json"


def essential_set() -> set[str]:
    return {line.strip() for line in ESSENTIAL.read_text().splitlines() if line.strip()}


def donor_pool(genes: tuple[str, ...]) -> set[str]:
    """Targets carried by any production donor that are also on the official gene axis."""
    on_axis = set(genes)
    pool: set[str] = set()
    for name in DONORS:
        perts = (ROOT / "data/donors" / name / "perts.txt").read_text().split("\n")
        pool |= {p for p in map(str, perts) if p in on_axis}
    return pool


def main() -> None:
    panel_abc, genes = C.panel_and_genes(C.PANELS["abc"])
    panel_pess, _ = C.panel_and_genes(C.PANELS["def199_pessimistic"])
    ess = essential_set()

    pool = sorted(donor_pool(genes) - set(panel_abc) - ess)
    rng = np.random.default_rng(0)
    panel_ef = tuple(str(x) for x in rng.choice(np.array(pool), 300, replace=False))

    result: dict[str, object] = {
        "construction": (
            "300 drawn seed 0 from (targets of k562_gwps|orion_hct116|orion_hek293t on the "
            "official axis) - (A/B/C panel) - (data/donors/k562_essential/perts.txt)"
        ),
        "pool_size": len(pool),
        "essential_list_size": len(ess),
        "panels": {
            "abc": {
                "n": len(panel_abc),
                "essential_genes": len(set(panel_abc) & ess),
                "overlap_with_abc": len(set(panel_abc) & set(panel_abc)),
            },
            "def199_pessimistic": {
                "n": len(panel_pess),
                "essential_genes": len(set(panel_pess) & ess),
                "overlap_with_abc": len(set(panel_pess) & set(panel_abc)),
            },
            "essential_free_235": {
                "n": len(set(panel_ef)),
                "essential_genes": len(set(panel_ef) & ess),
                "overlap_with_abc": len(set(panel_ef) & set(panel_abc)),
            },
        },
        "panel_sha256": M.operator_axis_sha256(panel_ef),
    }
    checks = [
        ("the essential-free panel has 300 unique targets", len(set(panel_ef)) == 300),
        ("it is disjoint from A/B/C", not (set(panel_ef) & set(panel_abc))),
        ("it contains no essential gene", not (set(panel_ef) & ess)),
        ("A/B/C itself has zero essential genes (8.1g)", not (set(panel_abc) & ess)),
        ("199's pessimistic bundle has 30 (8.1g)", len(set(panel_pess) & ess) == 30),
    ]

    rows: dict[str, object] = {}
    prior = json.loads((HERE / "work" / "carry235.json").read_text())
    for name in sorted(M.UNION_SOURCE_REGISTRY):
        covered = C.coverage(name, panel_ef)
        row: dict[str, object] = {
            "covered_essential_free": covered,
            "covered_abc": prior["coverage"]["abc"][name],
            "covered_def199_pessimistic": prior["coverage"]["def199_pessimistic"][name],
        }
        if covered:
            sized = C.size_sources((name,), panel_ef, genes)
            row["solved_fraction"] = sized["solved_fraction"]
            row["all_rows_would_solve"] = sized["all_rows_would_solve"]
            row["carrying_rows"] = sized["frame"]["carrying_rows"]
            checks.append((
                f"essential-free/{name}: carrying rule sizes above zero",
                float(sized["solved_fraction"]) > 0.0,
            ))
            checks.append((
                f"essential-free/{name}: all-rows rule delivers zero",
                float(sized["all_rows_would_solve"]) == 0.0,
            ))
        rows[name] = row
        print(f"[{name:20}] essential-free {covered:3d} | abc {row['covered_abc']:3d} | "
              f"pessimistic {row['covered_def199_pessimistic']:3d}"
              + (f" | solved {row['solved_fraction']:.6f}" if covered else "  <- carries NOTHING"),
              flush=True)

    # the October four, as a union, if any of them carries anything at all
    present = tuple(n for n in C.OCTOBER_FOUR if rows[n]["covered_essential_free"])
    if present:
        union = C.size_sources(present, panel_ef, genes)
        result["october_union"] = {
            "sources": list(present),
            "solved_fraction": union["solved_fraction"],
            "all_rows_would_solve": union["all_rows_would_solve"],
            "carrying_rows": union["frame"]["carrying_rows"],
        }
    else:
        result["october_union"] = {"sources": [], "note": "none of the four carries a target"}

    result["sources"] = rows
    result["checks"] = [{"check": c, "verdict": "PASS" if ok else "FAIL"} for c, ok in checks]
    result["all_pass"] = all(ok for _, ok in checks)
    OUT.write_text(json.dumps(result, indent=1, sort_keys=True))
    for check, ok in checks:
        print(f"[{'PASS' if ok else 'FAIL'}] {check}")
    print(json.dumps(result["panels"], indent=1))
    print(f"\nESSENTIAL-FREE ADDENDUM {'PASS' if result['all_pass'] else 'FAIL'}")
    raise SystemExit(0 if result["all_pass"] else 1)


if __name__ == "__main__":
    main()
