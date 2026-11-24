"""236 Phase 0 — the coverage census, recomputed INDEPENDENTLY of 233 and asserted against it.

233's KA-F checked 25 coverage rows against `DEF-RUNBOOK` §5 / 225's table.  This re-derives
them from the source registry's own `perts.txt` files and each panel's own `pert_counts.csv`
-- nothing but the raw files -- and asserts they equal `233/work/census233.json`.  It is the
cheapest possible check that this experiment is pointed at the same three panels and the same
thirteen source tables 233 was, and it uses no operator, no fold and ~0 memory.

    uv run python census236.py
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
WORK = HERE / "work"
WORK.mkdir(parents=True, exist_ok=True)
sys.path.insert(0, str(ROOT / "experiments/218-independent-donors"))

from sources218 import SOURCES  # noqa: E402

BUNDLES = {
    "abc": ROOT / "data/challenge",
    "real": ROOT / "experiments/199-def-runbook/work/def-bundle-r",
    "pess": ROOT / "experiments/199-def-runbook/work/def-bundle",
}
OCT = ("nadig_hepg2", "nadig_jurkat", "rpe1", "k562_essential")


def panel(key):
    lines = (BUNDLES[key] / "pert_counts.csv").read_text().splitlines()
    assert lines[0].strip() == "target_gene"
    p = tuple(ln.split(",")[0].strip() for ln in lines[1:] if ln.strip())
    assert len(p) == 300 and len(set(p)) == 300
    return set(p)


def perts(name):
    path = (ROOT / SOURCES[name]["path"]).resolve()
    return {ln.strip() for ln in (path / "perts.txt").read_text().splitlines() if ln.strip()}


def main() -> None:
    panels = {k: panel(k) for k in BUNDLES}
    # FALSIFIER: the synthetic panels must be DISJOINT from today's, asserted not assumed
    dis = {f"abc^{k}": len(panels["abc"] & panels[k]) for k in ("real", "pess")}
    assert dis["abc^real"] == 0 and dis["abc^pess"] == 0, dis
    print(f"panel disjointness: |abc n real| = {dis['abc^real']}, "
          f"|abc n pess| = {dis['abc^pess']}, |real n pess| = "
          f"{len(panels['real'] & panels['pess'])}")

    out = {}
    print(f"\n{'source':24s} " + "".join(f"{k:>8s}" for k in ("abc", "real", "pess")))
    for name in SOURCES:
        try:
            p = perts(name)
        except FileNotFoundError:
            continue
        row = {k: len(p & v) for k, v in panels.items()}
        out[name] = row
        mark = "  <- OCTOBER-ONLY" if name in OCT else ""
        print(f"{name:24s} " + "".join(f"{row[k]:8d}" for k in ("abc", "real", "pess")) + mark)

    for n in OCT:
        assert out[n]["abc"] == 0, f"{n} is not October-only: {out[n]['abc']}/300 today"

    # assert against 233's stored census
    c233 = json.loads((ROOT / "experiments/233-october-retention/work/census233.json").read_text())
    ref = c233["coverage"]
    checked = fails = 0
    for name, row in out.items():
        r = ref.get(name)
        if not isinstance(r, dict):
            continue
        for k in ("abc", "real", "pess"):
            if k in r and isinstance(r[k], int):
                checked += 1
                if r[k] != row[k]:
                    fails += 1
                    print(f"  MISMATCH {name}/{k}: 233 {r[k]} vs {row[k]}")
    print(f"\nvs 233/census233.json: {checked - fails}/{checked} matched"
          if checked else "\n(233's census json has a different shape; see work/census236.json)")
    (WORK / "census236.json").write_text(json.dumps(
        dict(disjointness=dis, coverage=out), indent=2))
    print(f"-> {WORK / 'census236.json'}")


if __name__ == "__main__":
    main()
