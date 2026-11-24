"""184 -- aggregate a phase2 run into the tables that go in RESULTS.md."""
from __future__ import annotations

import collections
import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))

ARMS = ["M", "gain-mean", "gain-basal", "gain-basal-clip", "gain-basal/perm-line",
        "gain-ORACLE", "bilinear", "bilinear-nc", "bilinear/perm-line", "bilinear/perm-drug",
        "M/perm-drug", "KNOWN-truth", "KNOWN-zero"]
METRICS = ["sign", "sign_weak", "cos", "cos_orth"]


def main(path):
    d = json.load(open(path))
    rows = d["rows"]
    ks = sorted({r["k"] for r in rows})
    # per (held line, k, arm): average over draws first, so every line weighs the same
    byline = collections.defaultdict(dict)
    for r in rows:
        byline[(r["held"], r["k"], r["arm"])].setdefault("v", []).append(r)
    per = {}
    for (h, k, arm), v in byline.items():
        per[(h, k, arm)] = {m: float(np.mean([x[m] for x in v["v"] if m in x]))
                            for m in METRICS if any(m in x for x in v["v"])}
    lines = sorted({h for h, _, _ in per})
    print(f"{len(lines)} held lines, k in {ks}, config {d['config']}\n")

    for m in METRICS:
        print(f"===== {m}  (macro over {len(lines)} held lines)")
        hdr = "  ".join(f"k={k:<7}" for k in ks)
        print(f"{'arm':26s} {hdr}")
        for arm in ARMS:
            cells = []
            for k in ks:
                vals = [per[(h, k, arm)][m] for h in lines
                        if (h, k, arm) in per and m in per[(h, k, arm)]]
                cells.append(f"{np.mean(vals):+.4f} " if vals else "  --    ")
            print(f"{arm:26s} " + "  ".join(f"{c:<9}" for c in cells))
        # deltas vs M, with the across-line sem and the win count
        print(f"{'-- delta vs M':26s}")
        for arm in ARMS:
            if arm in ("M", "KNOWN-truth", "KNOWN-zero"):
                continue
            cells = []
            for k in ks:
                dv = [per[(h, k, arm)][m] - per[(h, k, "M")][m] for h in lines
                      if (h, k, arm) in per and m in per[(h, k, arm)] and m in per[(h, k, "M")]]
                if not dv:
                    cells.append("  --    ")
                    continue
                dv = np.asarray(dv)
                cells.append(f"{dv.mean():+.4f}({(dv > 0).sum():2d}/{len(dv)})")
            print(f"{arm:26s} " + "  ".join(f"{c:<9}" for c in cells))
        print()


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else os.path.join(HERE, "work/phase2-main.json"))
