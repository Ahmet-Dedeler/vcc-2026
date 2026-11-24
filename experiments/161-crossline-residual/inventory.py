"""161 SS1 — the pairwise shared-target inventory, and the leakage audit.

`docs/CURRENT_STATE.md` closed donor expansion because `rpe1`, `nadig_jurkat` and
`k562_essential` carry **zero** of the 300 challenge panel targets.  That closure is about
PANEL COVERAGE AT INFERENCE.  Training a residual model needs only targets shared between a
pair of lines, and it needs them to be OFF the panel (otherwise the training set and the
deployment set are the same rows).  So the inventory that matters here is the opposite one:
how many NON-panel targets does each pair share, and how many distinct cell LINES cover each
target.

Outputs
-------
  * per-table census (targets, genes, panel overlap, blob flag)
  * pairwise shared-target matrix, split panel / non-panel
  * the multi-line coverage histogram: how many targets are seen by >= k distinct LINES
  * candidate training surfaces (a set of lines, their common targets and common genes)
"""
from __future__ import annotations

import itertools
import json
import os

import numpy as np

from common import BLOBBY, TABLES, WORK, load_table, panel_targets

os.makedirs(WORK, exist_ok=True)


def main() -> None:
    panel = panel_targets()
    print(f"challenge panel: {len(panel)} targets\n")

    tabs = {}
    print(f"{'table':<16}{'line':<9}{'targets':>9}{'genes':>8}{'panel':>7}{'nonpanel':>10}"
          f"{'cells/pert':>12}  note")
    for tid in TABLES:
        t = load_table(tid)
        if t is None:
            print(f"{tid:<16}{'-':<9}{'MISSING':>9}")
            continue
        tabs[tid] = t
        ps = set(t["perts"])
        npan = len(ps & panel)
        cp = (t["meta"].get("cells_per_pert") or {}).get("median")
        note = "BLOB (159 §4: within-line cos +0.17..+0.28)" if tid in BLOBBY else ""
        print(f"{tid:<16}{t['line']:<9}{len(ps):9d}{len(t['genes']):8d}{npan:7d}"
              f"{len(ps) - npan:10d}{str(cp):>12}  {note}")

    print(f"\n{'=' * 108}\nPAIRWISE SHARED TARGETS  (non-panel is the trainable pool)")
    print(f"{'A':<16}{'B':<16}{'same line?':<11}{'shared':>8}{'panel':>7}{'nonpanel':>10}"
          f"{'genes':>8}")
    rows = []
    for a, b in itertools.combinations(sorted(tabs), 2):
        sa, sb = set(tabs[a]["perts"]), set(tabs[b]["perts"])
        sh = sa & sb
        if not sh:
            continue
        gg = set(tabs[a]["genes"]) & set(tabs[b]["genes"])
        same = tabs[a]["line"] == tabs[b]["line"]
        r = dict(A=a, B=b, lineA=tabs[a]["line"], lineB=tabs[b]["line"], same_line=same,
                 shared=len(sh), panel=len(sh & panel), nonpanel=len(sh - panel),
                 genes=len(gg))
        rows.append(r)
    rows.sort(key=lambda r: -r["nonpanel"])
    for r in rows:
        print(f"{r['A']:<16}{r['B']:<16}{('SAME' if r['same_line'] else ''):<11}"
              f"{r['shared']:8d}{r['panel']:7d}{r['nonpanel']:10d}{r['genes']:8d}")

    # ---- how many distinct LINES cover each target -----------------------------------
    from collections import Counter
    cov = Counter()
    for tid, t in tabs.items():
        for p in set(t["perts"]):
            cov[p] = cov[p]
    lines_per = {}
    for tid, t in tabs.items():
        for p in set(t["perts"]):
            lines_per.setdefault(p, set()).add(t["line"])
    hist = Counter(len(v) for v in lines_per.values())
    print(f"\n{'=' * 108}\nDISTINCT CELL LINES PER TARGET  ({len(lines_per)} targets seen at all)")
    print(f"{'>= k lines':>12}{'all targets':>14}{'non-panel':>12}{'panel':>8}")
    for k in range(1, max(hist) + 1):
        sel = [p for p, v in lines_per.items() if len(v) >= k]
        print(f"{k:12d}{len(sel):14d}{len([p for p in sel if p not in panel]):12d}"
              f"{len([p for p in sel if p in panel]):8d}")

    # ---- candidate training surfaces --------------------------------------------------
    CANDIDATES = {
        "core6": ["k562_gwps", "k562_essential", "rpe1", "nadig_jurkat",
                  "orion_hct116", "orion_hek293t"],
        "core5-noess": ["k562_gwps", "rpe1", "nadig_jurkat", "orion_hct116", "orion_hek293t"],
        "wide-h1": ["k562_gwps", "orion_hct116", "orion_hek293t", "h1_2025"],
        "wide-h1-kolf": ["k562_gwps", "orion_hct116", "orion_hek293t", "h1_2025", "kolf2"],
        "wide-cd4": ["k562_gwps", "orion_hct116", "orion_hek293t", "cd4_rest"],
        "all11": list(TABLES),
    }
    print(f"\n{'=' * 108}\nCANDIDATE TRAINING SURFACES  (strict intersection of targets AND genes)")
    print(f"{'surface':<16}{'tables':>7}{'lines':>7}{'targets':>9}{'nonpanel':>10}"
          f"{'genes':>8}  members")
    surfaces = {}
    for name, mem in CANDIDATES.items():
        mem = [m for m in mem if m in tabs]
        if len(mem) < 3:
            continue
        T = set.intersection(*[set(tabs[m]["perts"]) for m in mem])
        G = set.intersection(*[set(tabs[m]["genes"]) for m in mem])
        G = G - T                       # pds_cosine drops every panel target coordinate
        lines = sorted({tabs[m]["line"] for m in mem})
        surfaces[name] = dict(tables=mem, lines=lines, n_targets=len(T),
                              n_nonpanel=len(T - panel), n_genes=len(G))
        print(f"{name:<16}{len(mem):7d}{len(lines):7d}{len(T):9d}{len(T - panel):10d}"
              f"{len(G):8d}  {','.join(lines)}")

    # ---- leakage audit ----------------------------------------------------------------
    print(f"\n{'=' * 108}\nLEAKAGE AUDIT — table x evaluation-line overlap")
    print("A fold that evaluates recipient line L must exclude EVERY table whose line is L.")
    for line in sorted({t["line"] for t in tabs.values()}):
        mem = [tid for tid, t in tabs.items() if t["line"] == line]
        print(f"  recipient line {line:<8} -> must hold out tables: {', '.join(sorted(mem))}")
    print("  h1-test harness: `data/donors/h1_2025` covers 85/85 of its targets "
          "(AGENTS.md §6.4).  This experiment never scores on the h1-test HARNESS; it scores "
          "on the h1_2025 DONOR TABLE with that table held out of the donor pool, which is "
          "the same exclusion enforced structurally.")

    out = dict(panel=len(panel), pairs=rows, surfaces=surfaces,
               lines_per_target_hist={str(k): v for k, v in sorted(hist.items())},
               tables={tid: dict(line=t["line"], n_targets=len(t["perts"]),
                                 n_genes=len(t["genes"]),
                                 n_panel=len(set(t["perts"]) & panel),
                                 blob=tid in BLOBBY) for tid, t in tabs.items()})
    with open(os.path.join(WORK, "inventory.json"), "w") as fh:
        json.dump(out, fh, indent=1)
    print(f"\nwrote {os.path.join(WORK, 'inventory.json')}")


if __name__ == "__main__":
    main()
