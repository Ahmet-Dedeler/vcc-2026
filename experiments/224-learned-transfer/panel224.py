"""224 Phase 0d — the PANEL-restricted leave-one-LINE-out surfaces.

S5/S7/S8 are built on the multi-line target intersection, which is dominated by
common-essential genes: the 2026 panel excludes those by construction (099 §1.1, and 159
§5.3 measures the challenge-panel cosines at 3-5x smaller than the broad ones).  So a rung
that clears on S7 has to be re-read on the surface production actually ships into.

These surfaces are 159 §5.3's: the four panel-covering donor lines, plus CD4-T as a fifth.
`kolf2` is a z-score matrix and 159 §4 shows it is anti-aligned and carries its own blob;
`cd4_rest` is `log2FC x ln2`.  Both are kept because production keeps them, and both are
flagged in the results.

    uv run python panel224.py
"""
from __future__ import annotations

import json

import numpy as np

from substrate224 import LINES, PANEL_SET, ROOT, WORK, load_basal, load_effects  # noqa: E402

CD4 = "data/external/cd4t/cloud-compact/cd4_compact/cd4_rest"
EXTRA = {"cd4_rest": dict(path=CD4, group="CD4T")}

SETS = {
    "P4": ["k562_gwps", "orion_hct116", "orion_hek293t", "kolf2"],
    "P5": ["k562_gwps", "orion_hct116", "orion_hek293t", "kolf2", "cd4_rest"],
    "P6": ["k562_gwps", "orion_hct116", "orion_hek293t", "kolf2", "cd4_rest", "feng_ipsc"],
    "P7": ["k562_gwps", "orion_hct116", "orion_hek293t", "kolf2", "cd4_rest", "feng_ipsc",
           "dld1_lowmoi"],
}


def meta(name):
    return EXTRA.get(name) or LINES[name]


def sets_for(members):
    ps, gs = [], []
    for m in members:
        p = ROOT / meta(m)["path"]
        ps.append(set(x.strip() for x in (p / "perts.txt").read_text().splitlines()
                      if x.strip()))
        gs.append(set(x.strip() for x in (p / "genes.txt").read_text().splitlines()
                      if x.strip()))
    targets = sorted(set.intersection(*ps) & PANEL_SET)
    genes = sorted(set.intersection(*gs) - PANEL_SET)          # drop all 300 panel coords
    return targets, genes


def main() -> None:
    rows = []
    for name, members in SETS.items():
        targets, genes = sets_for(members)
        rows.append(dict(name=name, members=members, n_targets=len(targets),
                         n_genes=len(genes)))
        print(f"{name}: {len(members)} lines, {len(targets)} PANEL targets x {len(genes)} genes"
              f"   {members}")
        if len(targets) < 40:
            print("   -> too few targets for a leave-one-line-out; not cached")
            continue
        E = np.stack([load_effects(meta(m)["path"], targets, genes) for m in members])
        B = np.stack([load_basal(meta(m)["path"], genes) for m in members])
        C = E - E.mean(axis=1, keepdims=True)
        np.savez_compressed(
            WORK / f"surface-{name}.npz", C=C.astype(np.float32),
            basal=B.astype(np.float32), members=np.asarray(members),
            groups=np.asarray([meta(m)["group"] for m in members]),
            targets=np.asarray(targets), genes=np.asarray(genes))
        print(f"   cached surface-{name}.npz  C{C.shape}")
    (WORK / "panel224.json").write_text(json.dumps(rows, indent=1))


if __name__ == "__main__":
    main()
