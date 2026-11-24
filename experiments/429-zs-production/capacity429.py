"""429: what the zs pass moves on H1, base vs MDE carriers, from the stored stats (pre/post lfc of the same blocks) and
323's reference DE (lfc, q). Per target: genes moved (target gene excluded), Σ|Δl| over them, and among the moved genes
the reference-DE (q < 0.05) ones: count, delivered sign accuracy, median |l|/|r| before/after, and how many overshoot
(|l_after| > |r| with the right sign). Also the nmae efficiency: per-target Δnmae per unit Σ|Δl| on the scored targets.
usage: uv run python capacity429.py [SEEDS]  -> work/capacity429.json"""
import json
import sys
from pathlib import Path

import numpy as np
import polars as pl

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
W = HERE / "work/h1"
SEEDS = [int(s) for s in (sys.argv[1] if len(sys.argv) > 1 else "0,1,2").split(",")]
MDE = "mde0.1cap0.1min0.3"
ref = np.load(ROOT / "experiments/323-current-family-validation/work/reference.npz", allow_pickle=False)
R_l, R_q, keep, genes, targets = ref["lfc"], ref["q"], ref["keep"], ref["genes"], ref["targets"]
gk = np.asarray(genes)[keep].astype(str)


def per(tag):
    d = pl.read_csv(W / f"per-target-{tag}.csv").filter(pl.col("metric") == "de_wilcoxon_lfc_nmae")
    return dict(zip(d["perturbation"].to_list(), d["value"].to_list()))


out = {}
for arm, pre, post in (("base", "base_s{s}_bps", "base_s{s}_bps_dnk250zs"),
                       ("MDE", f"base_s{{s}}_bps_{MDE}", f"base_s{{s}}_bps_{MDE}_dnk250zs")):
    recs = []
    for s in SEEDS:
        a = np.load(W / f"stats-{pre.format(s=s)}.npz")
        b = np.load(W / f"stats-{post.format(s=s)}.npz")
        assert np.array_equal(a["targets"], targets) and np.array_equal(b["targets"], targets)
        na, nb = per(pre.format(s=s)), per(post.format(s=s))
        for i, t in enumerate(map(str, targets)):
            l0, l1 = a["lfc"][i], b["lfc"][i]
            mv = (l1 != l0) & (gk != t)
            de = (R_q[i] < 0.05) & (gk != t)
            m = mv & de
            r = R_l[i]
            right = np.sign(l0[m]) == np.sign(r[m])
            recs.append(dict(seed=s, target=t, n_moved=int(mv.sum()), dl=float(np.abs(l1[mv] - l0[mv]).sum()),
                             n_moved_de=int(m.sum()), dl_de=float(np.abs(l1[m] - l0[m]).sum()),
                             acc=float(right.mean()) if m.any() else np.nan,
                             ratio_before=float(np.median(np.abs(l0[m]) / np.abs(r[m]))) if m.any() else np.nan,
                             ratio_after=float(np.median(np.abs(l1[m]) / np.abs(r[m]))) if m.any() else np.nan,
                             overshoot=int((right & (np.abs(l1[m]) > np.abs(r[m]))).sum()),
                             dnmae=(nb[t] - na[t]) if t in na and t in nb else np.nan))
    R = recs
    sc = [x for x in R if np.isfinite(x["dnmae"])]
    tot_de = sum(x["n_moved_de"] for x in R)
    acc_w = sum(x["acc"] * x["n_moved_de"] for x in R if x["n_moved_de"]) / max(tot_de, 1)
    out[arm] = dict(blocks=len(R), n_moved_mean=float(np.mean([x["n_moved"] for x in R])),
                    dl_sum_mean=float(np.mean([x["dl"] for x in R])),
                    n_moved_de_mean=float(np.mean([x["n_moved_de"] for x in R])),
                    dl_de_sum_mean=float(np.mean([x["dl_de"] for x in R])),
                    de_share_of_dl=float(sum(x["dl_de"] for x in R) / max(sum(x["dl"] for x in R), 1e-12)),
                    sign_acc_moved_de=float(acc_w),
                    ratio_before_median=float(np.nanmedian([x["ratio_before"] for x in R])),
                    ratio_after_median=float(np.nanmedian([x["ratio_after"] for x in R])),
                    overshoot_total=int(sum(x["overshoot"] for x in R)),
                    scored_targets=len(sc),
                    dnmae_mean_scored=float(np.mean([x["dnmae"] for x in sc])),
                    dl_mean_scored=float(np.mean([x["dl"] for x in sc])),
                    dnmae_per_unit_dl=float(np.mean([x["dnmae"] for x in sc]) / np.mean([x["dl"] for x in sc])))
    print(arm, json.dumps(out[arm]))
(HERE / "work/capacity429.json").write_text(json.dumps(out, indent=1))
