"""189 Phase 1a — the fitting inventory: which LINE pairs exist, how big, and what they cover.

The prior `r_g` is estimated from PAIRS of cell lines.  What is available:

  * six donor tables on disk, but only FIVE distinct cell LINES -- `k562_gwps` and
    `k562_essential` are two libraries of one line (161 `common.py`), so a leave-one-LINE-out
    fold that splits them leaks (159 measured the same-line pair at +0.2712 against
    +0.03..+0.16 for genuine cross-line pairs).
  * `h1_2025` is BANNED: 85/85 coverage of h1-test's targets (AGENTS.md trap #1).

This script prints, for every unordered pair of donor TABLES:
  shared targets, of which on the 300-gene challenge panel, and shared genes on the working
  axis; plus the h1-test gene-axis coverage of each table.

The working gene axis is h1-test's 10,779-gene DE axis intersected with the genes carried by
at least two donor tables -- the prior only has to be defined where it will be read.
"""
from __future__ import annotations

import itertools
import json
import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..",
                                "172-sign-confidence-key"))
from common189 import TABLES, WORK, load_table, panel_targets  # noqa: E402

FIT_TABLES = ["orion_hct116", "orion_hek293t", "k562_gwps", "k562_essential", "rpe1",
              "nadig_jurkat"]


def h1test_genes() -> list[str]:
    from harness172 import H1Test
    return list(H1Test().genes)


def main() -> None:
    os.makedirs(WORK, exist_ok=True)
    panel = panel_targets()
    hg = h1test_genes()
    print(f"h1-test DE gene axis: {len(hg)} genes;  challenge panel: {len(panel)} targets\n")

    tabs = {t: load_table(t) for t in FIT_TABLES}
    tabs = {k: v for k, v in tabs.items() if v is not None}

    print(f"{'table':<16}{'LINE':<9}{'targets':>9}{'genes':>8}{'panel tgt':>11}"
          f"{'genes on h1 axis':>18}")
    for t, tab in tabs.items():
        on_h1 = len(set(tab["genes"]) & set(hg))
        print(f"{t:<16}{tab['line']:<9}{len(tab['perts']):9d}{len(tab['genes']):8d}"
              f"{len(set(tab['perts']) & panel):11d}{on_h1:18d}")

    # working axis: on the h1-test DE axis and carried by >= 2 donor tables
    from collections import Counter
    cnt = Counter()
    for tab in tabs.values():
        cnt.update(set(tab["genes"]))
    work = sorted(g for g in hg if cnt[g] >= 2)
    print(f"\nworking gene axis (h1-test axis AND >=2 donor tables): {len(work)}")
    for c in range(2, len(tabs) + 1):
        print(f"    carried by >= {c} tables: {sum(1 for g in work if cnt[g] >= c)}")

    print(f"\n{'A':<16}{'B':<16}{'same LINE':<11}{'targets':>9}{'panel':>7}"
          f"{'nonpanel':>10}{'work genes':>12}")
    rows = []
    for a, b in itertools.combinations(sorted(tabs), 2):
        sh = set(tabs[a]["perts"]) & set(tabs[b]["perts"])
        gg = set(tabs[a]["genes"]) & set(tabs[b]["genes"]) & set(work)
        if not sh:
            continue
        same = tabs[a]["line"] == tabs[b]["line"]
        rows.append(dict(A=a, B=b, lineA=tabs[a]["line"], lineB=tabs[b]["line"], same=same,
                         shared=len(sh), panel=len(sh & panel), nonpanel=len(sh - panel),
                         genes=len(gg)))
    rows.sort(key=lambda r: -r["shared"])
    for r in rows:
        print(f"{r['A']:<16}{r['B']:<16}{('SAME -- excluded' if r['same'] else ''):<11}"
              f"{r['shared']:9d}{r['panel']:7d}{r['nonpanel']:10d}{r['genes']:12d}")

    json.dump(dict(work_genes=work, pairs=rows,
                   h1_axis=len(hg), n_panel=len(panel)),
              open(os.path.join(WORK, "inventory189.json"), "w"))
    np.save(os.path.join(WORK, "work_genes.npy"), np.array(work))
    print(f"\nwrote work/inventory189.json and work/work_genes.npy")


if __name__ == "__main__":
    main()
