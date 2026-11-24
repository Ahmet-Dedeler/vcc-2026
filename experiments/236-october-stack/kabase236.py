"""236 Phase 0 — the shipped stack's own levels must reproduce 233 §3.1 on ALL SIX members.

`stack236.py` already asserts 233's 157 single-source DIRECTION rows.  This asserts the other
half: that the BASE this experiment builds on -- "the shipped stack as it will exist in
October", three donors + PerturbMatch + CD4 @ f=0.34 + the H1 union @ g=0.28, with each
panel's OWN re-fitted operator and OWN re-extracted CD4 banks -- is the same object 233
scored, on all six members and both synthetic panels.  If BASE differs, every gain below is a
gain over a different incumbent and nothing composes with 233's table.

    uv run python kabase236.py
"""
from __future__ import annotations

import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
WORK = HERE / "work"
MEMBERS = ("cos_agg", "pds", "nmae", "fid", "reach")


def main() -> None:
    r233 = {(r["panel"], r["fold"], r["arm"]): r for r in json.loads(
        (ROOT / "experiments/233-october-retention/work/score233.json").read_text())}
    r236 = [r for r in json.loads((WORK / "six236.json").read_text()) if r["arm"] == "BASE"]
    bad, checked, rows = 0, 0, []
    print(f"{'panel/fold':22s}{'member':9s}{'233':>12}{'236':>12}{'|diff|':>10}")
    for a in r236:
        b = r233.get((a["panel"], a["fold"], "BASE"))
        assert b is not None, f"233 has no BASE for {a['panel']}/{a['fold']}"
        for k in MEMBERS + ("jac89",):
            bv = b["jac"]["89"] if k == "jac89" else b[k]
            av = a["jac"]["89"] if k == "jac89" else a[k]
            d = abs(bv - av)
            checked += 1
            bad += d > 1e-12
            rows.append(dict(panel=a["panel"], fold=a["fold"], member=k,
                             v233=bv, v236=av, diff=d))
            print(f"{a['panel'] + '/' + a['fold']:22s}{k:9s}{bv:12.7f}{av:12.7f}{d:10.2e}")
    verdict = "PASS" if bad == 0 else f"FAIL ({bad}/{checked})"
    print(f"\nKA-BASE: {verdict} — {checked - bad}/{checked} values bit-identical to 233 §3.1")
    (WORK / "kabase236.json").write_text(json.dumps(
        dict(checked=checked, failed=bad, verdict=verdict, rows=rows), indent=2))
    assert bad == 0, verdict


if __name__ == "__main__":
    main()
