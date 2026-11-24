"""517: the real-emitter runs (drive517.sh) against their no-hook bases, beside price517's analytic numbers.
Bases: c4-* vs 435's production runs (same seed); g33-* vs 467's gated33; arc-* vs 462l's arcfull. Deltas raw; overall-equiv
with 467's spans; target bootstrap (2,000) per seed; leak = call_set_unsuppressible (mean per target) and clip fraction."""
import json
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]; HERE = Path(__file__).resolve().parent
F = "base_s{s}_bps_{em}_mde0.1cap0.1min0.3coarse"
BASE = {"c4": (ROOT / "experiments/435-zs-on-coarse/work/h1", "prodem"),
        "g33": (ROOT / "experiments/467-strength-gated-programme/work/h1_gated33", "prodem"),
        "arc": (ROOT / "experiments/462-broad16-lab-signature/work/h1_arcfull", "prodem"),
        "mz": (ROOT / "experiments/435-zs-on-coarse/work/h1", "prodem"),      # block swaps are read against CLEAN4
        "mzg": (ROOT / "experiments/435-zs-on-coarse/work/h1", "prodem")}
M = {"fid": "de_wilcoxon_direction_fidelity_yield_raw", "reach": "de_wilcoxon_direction_reach_raw",
     "nmae": "de_wilcoxon_lfc_nmae", "jac": "de_wilcoxon_sig_jaccard"}
SP = {"fid": 0.30, "reach": 0.90, "nmae": -0.60, "jac": 0.37}


def per_target(d, tag):
    df = pd.read_csv(d / f"per-target-{tag}.csv")
    return df.pivot_table(index="perturbation", columns="metric", values="value")


out = {}
rng = np.random.default_rng(5174)
for d in sorted(HERE.glob("work/h1_*")):
    mode = d.name[3:]
    if mode.endswith("_zs"):
        continue                     # zs runs are compared by hand (their pre-pass tag repeats the plain run)
    block = mode.split("-")[0]
    rows = []
    for s in (0, 1, 2):
        tag = F.format(s=s, em="hookem")
        if not (d / f"metrics-{tag}.json").exists():
            continue
        bd, bem = BASE[block]
        btag = F.format(s=s, em=bem)
        a, b = json.loads((d / f"metrics-{tag}.json").read_text()), json.loads((bd / f"metrics-{btag}.json").read_text())
        dr = {m: a["raw"][m] - b["raw"][m] for m in ("pds", "mse", *M)}
        pa, pb = per_target(d, tag), per_target(bd, btag)
        D = {k: (pa[M[k]] - pb[M[k]]).reindex(pb.index).values for k in M}
        ov = lambda ix: sum(np.nanmean(D[k][ix]) / SP[k] for k in SP) / 6     # noqa: E731
        n = len(pb)
        boots = np.array([ov(rng.integers(0, n, n)) for _ in range(2000)])
        da = json.loads((d / f"diag-{tag}.json").read_text()); db = json.loads((bd / f"diag-{btag}.json").read_text())
        rows.append({"seed": s, **{k: round(v, 6) for k, v in dr.items()}, "ov": round(ov(np.arange(n)), 6),
                     "ci95": [round(float(np.quantile(boots, .025)), 5), round(float(np.quantile(boots, .975)), 5)],
                     "unsupp": [round(float(np.mean([x.get("call_set_unsuppressible", 0) for x in db])), 1),
                                round(float(np.mean([x.get("call_set_unsuppressible", 0) for x in da])), 1)],
                     "clip": [round(float(np.mean([x.get("clip_frac", 0) for x in db])), 4),
                              round(float(np.mean([x.get("clip_frac", 0) for x in da])), 4)]})
    if rows:
        out[mode] = rows
        for r in rows:
            print(f"{mode:22s} s{r['seed']} pds {r['pds']:+.1e} mse {r['mse']:+.1e} | nmae {r['nmae']:+.5f} fid {r['fid']:+.5f} "
                  f"reach {r['reach']:+.5f} jac {r['jac']:+.5f} | ov {r['ov']:+.5f} {r['ci95']} | unsupp {r['unsupp']} clip {r['clip']}")
(HERE / "work/real517.json").write_text(json.dumps(out, indent=1))
