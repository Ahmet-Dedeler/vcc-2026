"""238 — falsifier 4: the ESSENTIAL FRACTION of every source's CARRYING ROWS, per fold, on
every panel.  Pure set arithmetic; no banks are loaded and no memory is used.

A PLOLO fold holding donor `held` out scores the targets `panel n coverage(held)`.  A source
`s` carries the rows `panel n coverage(held) n coverage(s)`.  That is exactly the
`carrying_rows` count `resid218.source_residual` produces, so this table can be checked against
236 §1's measured overlap matrix -- which it is, as a known-answer row.

The question it answers: are the October four's carrying rows the panel's ESSENTIAL genes?  If
they are, then the +18.9% relative cosine 236 measured is a property of a stratum Arc's panel
does not contain, and the asset is an artefact of 199's sampling model rather than of biology.

    uv run python carry238.py
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(HERE))
import panels238 as P  # noqa: E402

SOURCE_PATHS = {
    "nadig_hepg2": "experiments/166-donor-acquisition/donors/nadig_hepg2",
    "nadig_jurkat": "data/donors/nadig_jurkat",
    "rpe1": "data/donors/rpe1",
    "k562_essential": "data/donors/k562_essential",
    "k562_gwps": "data/donors/k562_gwps",
    "orion_hct116": "data/donors/orion_hct116",
    "orion_hek293t": "data/donors/orion_hek293t",
    "gse295214_h1": "experiments/115-gse295214-h1/work/compact-donor",
    "h1_2025": "data/donors/h1_2025",
}
OCTOBER = ("nadig_hepg2", "nadig_jurkat", "rpe1", "k562_essential")
# 236 §1's published carrying-row matrix -- the known-answer row for this table
KA_236 = {
    ("pess", "orion_hct116"): dict(nadig_hepg2=28, nadig_jurkat=28, rpe1=27,
                                   k562_essential=25, UNION=30),
    ("pess", "orion_hek293t"): dict(nadig_hepg2=29, nadig_jurkat=29, rpe1=28,
                                    k562_essential=26, UNION=31),
    ("pess", "k562_gwps"): dict(nadig_hepg2=32, nadig_jurkat=32, rpe1=31, UNION=34),
    ("real", "orion_hct116"): dict(nadig_hepg2=45, nadig_jurkat=48, rpe1=47,
                                   k562_essential=34, UNION=50),
    ("real", "orion_hek293t"): dict(nadig_hepg2=45, nadig_jurkat=48, rpe1=47,
                                    k562_essential=34, UNION=50),
    ("real", "k562_gwps"): dict(nadig_hepg2=45, nadig_jurkat=48, rpe1=47, UNION=50),
}
FORBIDDEN = {"k562_gwps": {"k562_essential"}}


def perts(rel: str) -> set[str]:
    return {ln.strip() for ln in (ROOT / rel / "perts.txt").read_text().splitlines()
            if ln.strip()}


def main() -> None:
    cover = {n: perts(p) for n, p in SOURCE_PATHS.items()}
    essential = perts("data/donors/k562_essential")
    out: dict = {}
    ka_pass = ka_total = 0

    for panel in P.BUNDLES:
        if not P.BUNDLES[panel].exists():
            continue
        pset = set(P.read_panel(panel))
        print(f"\n=== panel {panel}   (|panel n essential| = {len(pset & essential)})")
        print(f"{'fold (held)':16}{'n_targets':>10}   " +
              "".join(f"{s[:13]:>15}" for s in OCTOBER) + f"{'UNION':>10}")
        for held in P.CANCER_FOLDS:
            targets = pset & cover[held]
            allowed = [s for s in OCTOBER if s not in FORBIDDEN.get(held, set())]
            carry = {s: targets & cover[s] for s in allowed}
            union = set().union(*carry.values()) if carry else set()
            row = {}
            cells = []
            for s in OCTOBER:
                if s not in carry:
                    cells.append(f"{'(leak)':>15}")
                    continue
                c, e = len(carry[s]), len(carry[s] & essential)
                row[s] = dict(carrying=c, essential=e,
                              frac=(e / c) if c else None)
                cells.append(f"{c:>7}{('/' + str(e)):>4}"
                             + (f"{100 * e / c:>4.0f}%" if c else f"{'—':>4}"))
                ka = KA_236.get((panel, held), {}).get(s)
                if ka is not None:
                    ka_total += 1
                    ka_pass += int(ka == c)
                    if ka != c:
                        print(f"  ** KA MISMATCH {panel}/{held}/{s}: {c} != 236's {ka}")
            row["UNION"] = dict(carrying=len(union), essential=len(union & essential),
                                frac=(len(union & essential) / len(union)) if union else None)
            ka = KA_236.get((panel, held), {}).get("UNION")
            if ka is not None:
                ka_total += 1
                ka_pass += int(ka == len(union))
                if ka != len(union):
                    print(f"  ** KA MISMATCH {panel}/{held}/UNION: {len(union)} != {ka}")
            out[f"{panel}|{held}"] = dict(n_targets=len(targets), sources=row)
            print(f"{held:16}{len(targets):>10}   " + "".join(cells)
                  + f"{len(union):>6}/{len(union & essential):<4}")
        print("   cells are  carrying/essential  and the essential PERCENTAGE of carrying rows")

    print(f"\nKNOWN-ANSWER (236 §1's carrying-row matrix): {ka_pass}/{ka_total} exact")
    out["ka_236_carrying_matrix"] = dict(passed=ka_pass, total=ka_total)
    (P.WORK / "carry238.json").write_text(json.dumps(out, indent=2))
    print(f"-> {P.WORK / 'carry238.json'}")


if __name__ == "__main__":
    main()
