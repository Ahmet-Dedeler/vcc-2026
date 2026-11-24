"""224 Phase 2, falsifier 4 — the leak audit, computed rather than asserted.

States every donor <-> surface overlap this experiment could have leaked through:
  * which LINE identities are shared between tables (the LOLO split is on LINE, not table);
  * the h1-test trap (`data/donors/h1_2025` covers 85/85 of h1-test's targets, AGENTS.md
    §6.4) — is that table inside any surface used here?
  * 222's shipped union membership — which of this experiment's lines are production INPUTS,
    so that a fold holding them out is not independent of the shipped artifact;
  * the gene axis: `pds_cosine` drops all 300 panel target coordinates, so every surface must
    too, and the broad surfaces additionally drop their own intersection targets;
  * the panel-target content of each surface.

    uv run python leakaudit224.py
"""
from __future__ import annotations

import json

import numpy as np

from core224 import WORK, load_surface
from substrate224 import LINES, PANEL, PANEL_SET, ROOT

SURFACES = ("S5", "S7", "S8", "P4", "P5", "P6")
H1TEST = ROOT / "data/harness/h1-test.h5ad"
# production's shipped sources, from scripts/make_submission_pm_cd4.py / 222's union
PRODUCTION_INPUTS = {"k562_gwps", "orion_hct116", "orion_hek293t", "cd4_rest",
                     "cd4_stim8hr", "cd4_stim48hr", "perturbmatch"}
PRODUCTION_INPUTS_222 = PRODUCTION_INPUTS | {"h1_2025"}


def h1test_targets() -> list[str]:
    try:
        import anndata as ad
        a = ad.read_h5ad(H1TEST, backed="r")
        t = sorted(set(a.obs["target_gene"].astype(str)) - {"non-targeting"})
        a.file.close()
        return t
    except Exception:                                              # noqa: BLE001
        return []


def main() -> None:
    out = {}
    h1t = h1test_targets()
    print(f"h1-test targets: {len(h1t)}")
    for name, meta in LINES.items():
        p = ROOT / meta["path"]
        if not (p / "perts.txt").exists():
            continue
        have = set(x.strip() for x in (p / "perts.txt").read_text().splitlines() if x.strip())
        cov = len(have & set(h1t))
        if cov and h1t:
            print(f"   {name:18s} covers {cov:3d}/{len(h1t)} of h1-test "
                  f"({100 * cov / len(h1t):.1f}%)")
        out.setdefault("h1test_coverage", {})[name] = cov

    print("\nSURFACE AUDIT")
    for s in SURFACES:
        f = WORK / f"surface-{s}.npz"
        if not f.exists():
            continue
        C, basal, members, groups, targets, genes, dropped = load_surface(s)
        gset = set(genes)
        tset = set(targets)
        rec = dict(
            lines=members, groups=groups, dropped_same_line_tables=dropped,
            n_targets=len(targets), n_genes=len(genes),
            panel_targets=len(tset & PANEL_SET),
            panel_coords_in_gene_axis=len(gset & PANEL_SET),
            target_coords_in_gene_axis=len(gset & tset),
            duplicate_line_groups=[g for g in set(groups) if groups.count(g) > 1],
            production_inputs_present=sorted(set(members) & PRODUCTION_INPUTS_222),
            h1test_leaking_tables=sorted(
                m for m in members if out["h1test_coverage"].get(m, 0) > 0.5 * max(len(h1t), 1)))
        out.setdefault("surfaces", {})[s] = rec
        print(f"  {s}: {len(members)} LINES {members}")
        print(f"     targets {len(targets)} (panel {len(tset & PANEL_SET)}) x genes "
              f"{len(genes)};  panel coords in gene axis {len(gset & PANEL_SET)};  "
              f"own target coords in gene axis {len(gset & tset)}")
        print(f"     duplicate LINE groups: {rec['duplicate_line_groups'] or 'none'};  "
              f"production inputs present: {rec['production_inputs_present'] or 'none'}")

    print("\nWHAT THIS EXPERIMENT NEVER TOUCHED")
    for s in ("cell-eval2 / any scorer", "the whitening operator (no arm applies one)",
              "vcc / any submission package", "submissions/", "any file outside "
              "experiments/224-learned-transfer/", "the CD4 residual sizing rule "
              "(no residual is built here)"):
        print(f"   - {s}")

    (WORK / "leakaudit224.json").write_text(json.dumps(out, indent=1))
    print(f"\nwrote {WORK / 'leakaudit224.json'}")


if __name__ == "__main__":
    main()
