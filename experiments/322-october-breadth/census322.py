"""322 Phase 0/1 — the coverage census, and the question it settles before any scoring.

KA-F (233's row, reproduced): 233 published a 25-row coverage table for 15 sources on the
three panels.  Every one of those rows is recomputed here from `perts.txt` and asserted.

THE NEW COLUMN is `n_perts` -- how many perturbations the source contains AT ALL.  233's
table has only panel coverage.  Putting the two side by side is the whole of Phase 1:

    panel coverage on TODAY's panel   is not   what a source will carry in October.

`kolf2` is the extreme case (the mandate's candidate): 273 perturbations total, all 273
inside Arc's 300.  It is the broadest non-production source ON THE PANEL and the narrowest
but one AS A SOURCE.

    uv run python census322.py
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import lib322 as L  # noqa: E402

import panels233 as P  # noqa: E402

# 233's published census (`census233.py:33-40` / `work/census233.json`), verbatim.
KA_F = {
    "k562_gwps": (269, 300, 150), "orion_hct116": (290, 300, 291),
    "orion_hek293t": (294, 300, 293), "kolf2": (273, 0, 0), "feng_ipsc": (182, 0, 0),
    "gse295214_h1": (64, 57, 34), "h1_2025": (20, 7, 3), "nadig_hepg2": (0, 45, 33),
    "nadig_jurkat": (0, 48, 33), "rpe1": (0, 47, 32), "k562_essential": (0, 34, 30),
    "telohaec": (65, 42, 40), "myllia_a549": (14, 5, 4), "frangieh_melanoma": (5, 0, 2),
    "perturbmatch": (251, 273, 143),
}
PANELS = ("abc", "real", "pess")


def main() -> None:
    added = L.register_extra()
    import sources218
    panels = {k: set(P.read_panel(k)) for k in PANELS}
    for a in PANELS:
        for b in PANELS:
            if a != b and (a, b) != ("real", "pess") and (b, a) != ("real", "pess"):
                inter = panels[a] & panels[b]
                assert not inter or a == b, f"{a} n {b} = {len(inter)} (expected 0)"

    rows = []
    for name in sorted(sources218.SOURCES):
        try:
            perts = L.source_perts(name)
        except FileNotFoundError:
            continue
        pset = set(perts)
        cov = {k: len(pset & panels[k]) for k in PANELS}
        rows.append(dict(source=name, n_perts=len(pset),
                         abc=cov["abc"], real=cov["real"], pess=cov["pess"],
                         role=str(sources218.SOURCES[name].get("role", "")),
                         path=str(sources218.SOURCES[name]["path"])))

    by = {r["source"]: r for r in rows}
    ka = []
    for name, (a, r, p) in KA_F.items():
        got = by.get(name)
        ok = got is not None and (got["abc"], got["real"], got["pess"]) == (a, r, p)
        ka.append(dict(source=name, expected=[a, r, p],
                       measured=None if got is None else [got["abc"], got["real"], got["pess"]],
                       pass_=bool(ok)))
    n_pass = sum(k["pass_"] for k in ka)
    print(f"KA-F (233's census): {n_pass}/{len(ka)} rows reproduce exactly")
    for k in ka:
        if not k["pass_"]:
            print("  FAIL", k)
    assert n_pass == len(ka), "233's census did not reproduce — abort (AGENTS.md rule)"

    # --- the two senses of breadth, and which one predicts October -------------------
    sel = [r for r in rows if r["n_perts"] > 0]
    x_panel = np.array([r["abc"] for r in sel], float)
    x_broad = np.array([np.log10(r["n_perts"]) for r in sel])
    for tgt in ("real", "pess"):
        y = np.array([r[tgt] for r in sel], float)
        from scipy.stats import spearmanr
        rp = spearmanr(x_panel, y).statistic
        rb = spearmanr(x_broad, y).statistic
        print(f"  Spearman(coverage on ABC, coverage on {tgt})  = {rp:+.3f}   "
              f"Spearman(log10 n_perts, coverage on {tgt}) = {rb:+.3f}   n={len(sel)}")

    print(f"\n{'source':24s} {'n_perts':>8s} {'abc':>5s} {'real':>5s} {'pess':>5s}  role")
    for r in sorted(rows, key=lambda r: -r["abc"]):
        print(f"{r['source']:24s} {r['n_perts']:8d} {r['abc']:5d} {r['real']:5d} "
              f"{r['pess']:5d}  {r['role']}")

    out = dict(ka_f=ka, ka_f_pass=n_pass, rows=rows, extra_registered=added)
    (L.WORK / "census322.json").write_text(json.dumps(out, indent=2))
    print(f"\n-> {L.WORK / 'census322.json'}")


if __name__ == "__main__":
    main()
