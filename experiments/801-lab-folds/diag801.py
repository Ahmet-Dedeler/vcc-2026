"""801 -- panel diagnostics for the WORKING panels L1-L6 only (never L7-L8), written before any score: strength-bin shares,
donor coverage, CD4 coverage, coverage classes, basal-expression quantiles; versus A/B/C and versus 534's efr / efp rules
(238's ef-bundle-r / ef-bundle panels, the tuned-on instances of those rules). The bins and the basal CPM are draw801's
(Rule), so every panel is described in the draw's own terms.

    uv run python experiments/801-lab-folds/diag801.py      -> work/diag801.json   (light, < 1 GB, seconds)
"""
from __future__ import annotations

import sys

sys.dont_write_bytecode = True
import collections  # noqa: E402
import json  # noqa: E402
from pathlib import Path  # noqa: E402

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import draw801 as D  # noqa: E402

WORKING = ("L1", "L2", "L3", "L4", "L5", "L6")
OTHER = {"abc": D.SRC / "pert_counts.csv",
         "efp": D.ROOT / "experiments/238-essential-free-panel/work/ef-bundle/pert_counts.csv",
         "efr": D.ROOT / "experiments/238-essential-free-panel/work/ef-bundle-r/pert_counts.csv"}


def describe(R, panel):
    n = len(panel)
    cpm = np.asarray([R.cpm[R.pos[t]] for t in panel if t in R.pos])
    nd = np.asarray([R.nd.get(t, np.nan) for t in panel])
    sb = collections.Counter(R.sbin(t) for t in panel)
    cv = collections.Counter(str(R.cov(t)) for t in panel)
    return dict(
        n=n, essential=len(set(panel) & R.ess),
        coverage={d: len(set(panel) & R.donor[d]) for d in D.DONORS},
        cd4_rest=int(np.isfinite(nd).sum()),
        strength_bins={b: sb.get(b, 0) for b in D.QBINS + ("NA",)},
        coverage_class={c: cv.get(c, 0) for c in D.CLASSES + ("None",)},
        n_downstream_quantiles={f"q{int(p * 100)}": float(np.nanquantile(nd, p)) for p in (0.1, 0.25, 0.5, 0.75, 0.9)},
        basal_cpm_quantiles={f"q{p}": float(np.percentile(cpm, p)) for p in (5, 10, 25, 50, 75, 90)},
        below_basal_threshold=int((cpm < R.basal_thr).sum()))


def main():
    R = D.Rule()
    draws = json.loads(D.DRAWS.read_text())["draws"]
    out = dict(note="L1-L6 only; L7-L8 are sealed (no statistic)", basal_threshold_cpm=R.basal_thr,
               strength_edges=R.edges.tolist(), panels={})
    for k in WORKING:
        out["panels"][k] = describe(R, draws[k]["targets"])
    for k, p in OTHER.items():
        out["panels"][k] = describe(R, D.read_targets(p))
    (D.WORK / "diag801.json").write_text(json.dumps(out, indent=1))
    hdr = f"{'panel':6s} {'ess':>3s} {'K562':>4s} {'HCT':>4s} {'HEK':>4s} {'CD4':>4s} " + " ".join(f"{b:>3s}" for b in D.QBINS + ("NA",)) \
        + "  " + " ".join(f"{c:>3s}" for c in D.CLASSES) + "  nd_med  cpm_q5  cpm_q50  <thr"
    print(hdr)
    for k, d in out["panels"].items():
        print(f"{k:6s} {d['essential']:3d} {d['coverage']['k562_gwps']:4d} {d['coverage']['orion_hct116']:4d} "
              f"{d['coverage']['orion_hek293t']:4d} {d['cd4_rest']:4d} " + " ".join(f"{d['strength_bins'][b]:3d}" for b in D.QBINS + ("NA",))
              + "  " + " ".join(f"{d['coverage_class'][c]:3d}" for c in D.CLASSES)
              + f"  {d['n_downstream_quantiles']['q50']:6.1f}  {d['basal_cpm_quantiles']['q5']:6.2f}  {d['basal_cpm_quantiles']['q50']:7.2f}  {d['below_basal_threshold']:4d}")


if __name__ == "__main__":
    main()
