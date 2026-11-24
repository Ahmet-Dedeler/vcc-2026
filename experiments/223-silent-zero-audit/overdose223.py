"""223 — INBOX 1 §2: the THIRD failure mode of the median rule. Above 50% coverage it does not
zero the residual, it OVER-DOSES it.

`solve_cd4_fraction` sets `fraction = target_f / (median_all_rows / base_median)`.  Below 50%
coverage `median_all_rows` is 0 and the fraction is 0 (the silent zero).  Above 50% it is a
LOW ORDER STATISTIC of the carrying rows -- the median of a distribution padded with zeros --
so the fraction comes out too LARGE by exactly

        overdose = median(carrying norms) / median(all-row norms)

and the residual is added at `overdose` times the dose the caller asked for.  This is measured
here per source per fold, together with the percentile of the all-rows median inside the
carrying-row norm distribution, which is the quantity that makes the size of the error obvious.

    uv run python overdose223.py
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
WORK = HERE / "work"
sys.path.insert(0, str(ROOT / "experiments/162-retrodicting-validator"))

import plolo  # noqa: E402

plolo.WORK = WORK / "plolo"

from plolo import FOLDS, build_fold, cd4_for  # noqa: E402

sys.path.insert(0, str(HERE))
import sweep223  # noqa: E402  (registers the extra tables)

sys.path.insert(0, str(ROOT / "experiments/218-independent-donors"))
from resid218 import source_residual  # noqa: E402

SRC = ("telohaec", "dld1_lowmoi", "altair_a549", "myllia_a549", "frangieh_melanoma",
       "feng_ipsc", "kolf2", "cmap_sh_cgs", "cmap_xpr")
FOLD = "orion_hct116"


def stats(res) -> dict:
    added = np.where(np.asarray(res.gene_votes) > 0, np.asarray(res.mean_effect), 0.0)
    n = np.linalg.norm(added, axis=1)
    live = n[n > 0]
    med_all = float(np.median(n))
    med_car = float(np.median(live)) if live.size else 0.0
    pct = (float(100.0 * (live < med_all).mean()) if live.size and med_all > 0 else 0.0)
    return dict(rows=int(n.size), carrying=int(live.size),
                coverage=float(live.size / n.size),
                median_all=med_all, median_carrying=med_car,
                overdose_factor=(med_car / med_all) if med_all > 0 else None,
                all_median_percentile_of_carrying=pct)


def main() -> None:
    fold = build_fold(FOLD)
    targets = fold["targets"]
    out = {"cd4": stats(cd4_for(targets))}
    for s in SRC:
        out[s] = stats(source_residual(s, targets, "gene_norm"))
    (WORK / "overdose223.json").write_text(json.dumps(dict(fold=FOLD, sources=out), indent=2))

    hdr = (f"{'source':18s}{'carry/rows':>13s}{'cov':>8s}{'med(all)':>11s}"
           f"{'med(carry)':>12s}{'OVERDOSE x':>12s}{'pctile of carry':>17s}  effect of the shipped rule")
    print(hdr)
    print("-" * len(hdr))
    for k, v in sorted(out.items(), key=lambda kv: -kv[1]["coverage"]):
        od = v["overdose_factor"]
        if od is None:
            note = "SILENT ZERO -- residual identically absent"
            ods, pct = "  (zero)", "      -"
        elif od < 1.05:
            note = "negligible (coverage high enough that the two medians agree)"
            ods, pct = f"{od:12.3f}", f"{v['all_median_percentile_of_carrying']:17.1f}"
        else:
            note = f"OVER-DOSED by {od:.2f}x wherever this source was median-sized"
            ods, pct = f"{od:12.3f}", f"{v['all_median_percentile_of_carrying']:17.1f}"
        print(f"{k:18s}{v['carrying']:>6d}/{v['rows']:<6d}{v['coverage']:>8.1%}"
              f"{v['median_all']:>11.4f}{v['median_carrying']:>12.4f}{ods}{pct}  {note}")
    print(f"\n-> {WORK / 'overdose223.json'}")


if __name__ == "__main__":
    main()
