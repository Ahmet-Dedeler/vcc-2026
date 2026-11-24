"""522 — summary of the sub-panel (panel-mix) tests: merges work/subpanel522*.json (one file per invocation; the newest
record of a surface wins) and prints, per surface and mode, the mean over draws of r vs the sub-panel's own cap 3, the
percentile cap's hi, and paired draw-level differences (pct vs cap 6; pct x S0.25 vs pct; S1 vs pct; uncapped vs pct).
Panel means are over the three cancer folds.

    uv run python experiments/522-october-bulk-key/sub522sum.py
"""
from __future__ import annotations

import sys

sys.dont_write_bytecode = True
import json  # noqa: E402
from pathlib import Path  # noqa: E402

import numpy as np  # noqa: E402

HERE = Path(__file__).resolve().parent
WORK = HERE / "work"
ARMS = ("cap6", "pctPROD", "cap8", "uncapped", "S1", "pctPRODxS0.25", "pctPRODxS0.5", "cap6xS0.5")
PAIRS = (("pctPROD", "cap6"), ("pctPRODxS0.25", "pctPROD"), ("S1", "pctPROD"), ("uncapped", "pctPROD"))


def main():
    files = sorted(WORK.glob("subpanel522*.json"), key=lambda p: p.stat().st_mtime)
    R = {}
    for f in files:
        R.update(json.loads(f.read_text()))
    lines, out = [], {}
    P = lines.append
    modes = []
    for k, v in R.items():
        for m in v["modes"]:
            if m not in modes:
                modes.append(m)
    for m in modes:
        P(f"\n## mode: {m}   (r vs the sub-panel's cap 3, mean over draws; paired differences in points, fraction of draws > 0)")
        P(f"{'surface':22s} {'n_sub':>5s} {'hi':>6s} " + " ".join(f"{a:>13s}" for a in ARMS) + "  | "
          + " | ".join(f"{a}-{b}" for a, b in PAIRS))
        panel_acc = {}
        for key in sorted(R):
            rec = R[key]["modes"].get(m)
            if rec is None:
                continue
            A = rec["arms"]
            diffs = []
            for a, b in PAIRS:
                d = np.asarray(A[a]["r"]) - np.asarray(A[b]["r"])
                diffs.append((100 * d.mean(), float(np.mean(d > 0))))
            P(f"{key:22s} {rec['info']['n_sub']:5d} {rec['hi_mean']:6.2f} " + " ".join(f"{100 * A[a]['r_mean']:+13.2f}" for a in ARMS)
              + "  | " + " | ".join(f"{d:+.2f} ({fp:.2f})" for d, fp in diffs))
            p = key.split("/")[0]
            panel_acc.setdefault(p, []).append(({a: A[a]["r_mean"] for a in ARMS}, [d for d, _ in diffs], rec["hi_mean"]))
        for p, L in panel_acc.items():
            mean = {a: float(np.mean([x[0][a] for x in L])) for a in ARMS}
            dm = [float(np.mean([x[1][i] for x in L])) for i in range(len(PAIRS))]
            P(f"{'  ' + p + ' mean (' + str(len(L)) + ' folds)':22s} {'':5s} {np.mean([x[2] for x in L]):6.2f} "
              + " ".join(f"{100 * mean[a]:+13.2f}" for a in ARMS) + "  | " + " | ".join(f"{d:+.2f}" for d in dm))
            out.setdefault(m, {})[p] = dict(n_folds=len(L), r=mean, paired_points=dict(zip([f"{a}-{b}" for a, b in PAIRS], dm)),
                                            hi=float(np.mean([x[2] for x in L])))
    txt = "\n".join(lines)
    print(txt)
    (WORK / "sub522sum.txt").write_text(txt + "\n")
    (WORK / "sub522sum.json").write_text(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()
