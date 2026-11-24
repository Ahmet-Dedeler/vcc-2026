#!/usr/bin/env python3
"""385 step 3d: direct paired Kc1-vs-Kshuf contrast + fid mechanism split.

Reuses report385's resample construction (seed 20260909) so resamples are
identical. Reads work/results385.json for anchors/params. Prints:
 1. Kc1-Kshuf paired deltas per member (raw + H1-scaled) with 95% CIs.
 2. fid mechanism: mean directional precision and mean coverage per arm.
"""
import json
import sys
from pathlib import Path

import numpy as np
import polars as pl

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
BASE = ROOT / "experiments/323-current-family-validation"
WORK = HERE / "work"
for _p in (str(HERE), str(BASE), str(ROOT / "src")):
    if _p not in sys.path:
        sys.path.insert(0, _p)

import evaluate323 as E  # noqa: E402
import report323 as B  # noqa: E402
import report385 as R385  # noqa: E402


def per_target(name, src):
    table = pl.read_csv(src / f"per-target-{name}.csv")
    out = {}
    for m in E.MEMBERS:
        if m == "mse":
            continue
        rows = table.filter(pl.col("metric") == E.MEMBERS[m])
        out[m] = np.array([dict(zip(rows["perturbation"], rows["value"])).get(str(t), np.nan)
                           for t in TRUTH["targets"]], dtype=float)
    return out


TRUTH = {}
with np.load(BASE / "work/reference.npz", allow_pickle=False) as z:
    TRUTH["targets"] = z["targets"].astype(str)
    keep = z["keep"]
    genes = z["genes"].astype(str)
    uni = genes[keep]
    q = z["q"]
    TRUTH["n_real"] = np.array([((q[i] < 0.05) & (uni != t)).sum()
                                for i, t in enumerate(TRUTH["targets"])], dtype=float)

SRC = {"cap15": BASE / "work", "Kc1": WORK, "Kc2": WORK, "Kshuf": WORK}
N = len(TRUTH["targets"])
COUNTS = R385.resample_counts(N)
anchors = json.loads((BASE / "work/h1-anchors.json").read_text())
BM = list(E.MEMBERS)
base = np.array([anchors["baseline"][m] for m in BM])
rep = np.array([anchors["replicate"][m] for m in BM])

vals = {a: per_target(a, SRC[a]) for a in ("Kc1", "Kc2", "Kshuf", "cap15")}
print("== direct paired Kc1 - Kshuf (10,000 resamples) ==")
for m in ("fid", "jac", "reach", "nmae", "pds"):
    d = {"Kc1": None, "Kshuf": None}
    DD = []
    v1, v2 = vals["Kc1"][m], vals["Kshuf"][m]
    ok = np.isfinite(v1) & np.isfinite(v2)
    full = (COUNTS[:, ok] @ v1[ok]) / COUNTS[:, ok].sum(1) - \
           (COUNTS[:, ok] @ v2[ok]) / COUNTS[:, ok].sum(1)
    j = BM.index(m)
    s = full / abs(rep[j] - base[j])
    print(f"{m:>6}: raw {full[0]:+.6f} [{np.percentile(full[1:],2.5):+.6f},"
          f"{np.percentile(full[1:],97.5):+.6f}]  scaled {s[0]:+.5f} "
          f"[{np.percentile(s[1:],2.5):+.5f},{np.percentile(s[1:],97.5):+.5f}]")

# fid mechanism from delivered stats: need k (correct-sign calls) per target.
# Recompute cheaply: fid_p = k/max(n_pred,n_real) is in per-target csv; get n_pred
# from stats q and n_real from reference -> precision k/n_pred, coverage min(1,n_pred/n_real).
print("\n== fid mechanism (means over 85 targets) ==")
for a in ("cap15", "Kc1", "Kc2", "Kshuf"):
    with np.load(SRC[a] / f"stats-{a}.npz", allow_pickle=False) as z:
        qq = z["q"]
    genes_a = genes
    own_pos = np.array([np.where(uni == t)[0][0] if t in uni else -1 for t in TRUTH["targets"]])
    n_pred = np.array([(qq[i] < 0.05).sum() - (1 if own_pos[i] >= 0 and qq[i][own_pos[i]] < 0.05
                                                else 0) for i in range(N)], dtype=float)
    f = vals[a]["fid"]
    k = f * np.maximum(n_pred, TRUTH["n_real"])
    with np.errstate(invalid="ignore", divide="ignore"):
        prec = np.where(n_pred > 0, k / np.maximum(n_pred, 1), np.nan)
        cov = np.minimum(1, n_pred / np.maximum(TRUTH["n_real"], 1))
    print(f"{a:>6}: mean fid {np.nanmean(f):.4f}  mean precision {np.nanmean(prec):.4f}  "
          f"mean coverage {np.nanmean(cov):.4f}  mean n_pred {n_pred.mean():.0f} "
          f"(n_real mean {TRUTH['n_real'].mean():.0f})")
