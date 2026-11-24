"""224 Phase 0b — choose the leave-one-LINE-out surfaces and cache them.

The slope in `k` is the decision variable, so the surface that matters most is the one with
the MOST distinct cell lines, not the one with the most targets.  This script prints the
(targets x genes) cost of every nested line set and writes the caches the ladder reads.

    uv run python surfaces224.py
"""
from __future__ import annotations

import json
import os
from pathlib import Path

import numpy as np

os.environ.setdefault("OPENBLAS_NUM_THREADS", "6")
os.environ.setdefault("OMP_NUM_THREADS", "6")

from substrate224 import (LINES, PANEL_SET, ROOT, WORK, load_basal,  # noqa: E402
                          load_effects)

# Nested candidate line sets, cheapest first.  `k562_essential` is a SECOND LIBRARY of the
# K562 line, kept so the same-line positive control (159 §5.2) is expressible; it is fused
# with `k562_gwps` into one LINE for every split.
NESTED = [
    ("S5", ["k562_gwps", "k562_essential", "orion_hct116", "orion_hek293t", "rpe1",
            "nadig_jurkat"]),
    ("S6", ["k562_gwps", "k562_essential", "orion_hct116", "orion_hek293t", "rpe1",
            "nadig_jurkat", "nadig_hepg2"]),
    ("S7", ["k562_gwps", "k562_essential", "orion_hct116", "orion_hek293t", "rpe1",
            "nadig_jurkat", "nadig_hepg2", "gse295214_h1"]),
    ("S8", ["k562_gwps", "k562_essential", "orion_hct116", "orion_hek293t", "rpe1",
            "nadig_jurkat", "nadig_hepg2", "gse295214_h1", "dld1_lowmoi"]),
    ("S9", ["k562_gwps", "k562_essential", "orion_hct116", "orion_hek293t", "rpe1",
            "nadig_jurkat", "nadig_hepg2", "gse295214_h1", "dld1_lowmoi", "telohaec"]),
]

# The surfaces actually cached.  S5 is 159's exact six-line set (the known-answer row).
CACHE = ("S5", "S7", "S8")


def sets_for(members: list[str]) -> tuple[list[str], list[str]]:
    ps, gs = [], []
    for m in members:
        p = ROOT / LINES[m]["path"]
        ps.append(set(x.strip() for x in (p / "perts.txt").read_text().splitlines()
                      if x.strip()))
        gs.append(set(x.strip() for x in (p / "genes.txt").read_text().splitlines()
                      if x.strip()))
    targets = sorted(set.intersection(*ps))
    genes = sorted(set.intersection(*gs) - set(targets) - PANEL_SET)
    return targets, genes


def main() -> None:
    rows = []
    print(f"{'set':5s} {'lines':>6s} {'tables':>7s} {'targets':>8s} {'genes':>7s} "
          f"{'panel tgt':>10s}  members")
    for name, members in NESTED:
        targets, genes = sets_for(members)
        groups = sorted({LINES[m]["group"] for m in members})
        rows.append(dict(name=name, members=members, groups=groups,
                         n_targets=len(targets), n_genes=len(genes),
                         panel_targets=len(set(targets) & PANEL_SET)))
        print(f"{name:5s} {len(groups):6d} {len(members):7d} {len(targets):8d} "
              f"{len(genes):7d} {len(set(targets) & PANEL_SET):10d}  {','.join(groups)}")

    for name in CACHE:
        members = dict(NESTED)[name]
        targets, genes = sets_for(members)
        E = np.stack([load_effects(LINES[m]["path"], targets, genes) for m in members])
        B = np.stack([load_basal(LINES[m]["path"], genes) for m in members])
        # mean-remove each table over its own targets: that mean IS the metric's zero
        C = E - E.mean(axis=1, keepdims=True)
        np.savez_compressed(
            WORK / f"surface-{name}.npz",
            C=C.astype(np.float32), basal=B.astype(np.float32),
            members=np.asarray(members), groups=np.asarray([LINES[m]["group"] for m in members]),
            targets=np.asarray(targets), genes=np.asarray(genes))
        print(f"cached surface-{name}.npz  C{C.shape}  basal{B.shape}")

    (WORK / "surfaces224.json").write_text(json.dumps(rows, indent=1))


if __name__ == "__main__":
    main()
