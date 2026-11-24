"""518 — CD4 T-cell response strength (Marson GWCD4i DESeq2 stats, 3 conditions) as per-target magnitude information.

CD4 is a production INPUT on every PLOLO fold and never a truth, so its per-target strengths are leak-free predictors of the
folds' residual ideal scale u_p (feat518's response: the magnitude oracle's multiplier on 446bj2's fold bulk, from arms540's
saved `V1_0.5_3@R16` arrays). Features per condition (rest / stim8hr / stim48hr) on non-panel genes:
    nde   number of genes with adj p < 0.05            (DESeq2 on pseudobulk replicates: a well-powered strength)
    nde1  ... and |log2FC| > 0.5
    nrm   the delta row norm (delta = ln2 x log2FC)
and their means over the conditions carried. Screen: Spearman(u, f) per fold and the analytic tilt
w -> w (f/median)^beta (headroom540's first-order route), beta in [-0.5, 0.5]; targets CD4 does not carry keep w.

    uv run python experiments/518-bulk-mse-gap/cd4feat518.py
"""
from __future__ import annotations

import sys

sys.dont_write_bytecode = True
import json  # noqa: E402
from pathlib import Path  # noqa: E402

import numpy as np  # noqa: E402

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
WORK = HERE / "work"
CD4 = ROOT / "data/external/cd4t/cloud-compact/cd4_compact"
CONDS = ("cd4_rest", "cd4_stim8hr", "cd4_stim48hr")
W540 = ROOT / "experiments/540-lodo-clip-scan/work"
FOLDS = ("k562_gwps", "orion_hct116", "orion_hek293t", "kolf2")
CANCER = FOLDS[:3]
BASE = "V1_0.5_3@R16"
PANEL = set(l.split(",")[0] for l in (ROOT / "data/challenge/pert_counts.csv").read_text().splitlines()[1:] if l.strip())


def cd4_features(targets):
    feats = {}
    per = {k: [] for k in ("nde", "nde1", "nrm")}
    for c in CONDS:
        d = CD4 / c
        perts = [l.strip() for l in (d / "perts.txt").read_text().splitlines() if l.strip()]
        genes = [l.strip() for l in (d / "genes.txt").read_text().splitlines() if l.strip()]
        nonp = np.asarray([g not in PANEL for g in genes])
        pi = {p: i for i, p in enumerate(perts)}
        delta = np.load(d / "delta.npy", mmap_mode="r")
        adj = np.load(d / "adj_p_value.npy", mmap_mode="r")
        lfc = np.load(d / "log_fc.npy", mmap_mode="r")
        nde = np.full(len(targets), np.nan)
        nde1 = np.full(len(targets), np.nan)
        nrm = np.full(len(targets), np.nan)
        for k, t in enumerate(targets):
            i = pi.get(t)
            if i is None:
                continue
            a = np.asarray(adj[i], float)[nonp]
            l2 = np.asarray(lfc[i], float)[nonp]
            nde[k] = float((a < 0.05).sum())
            nde1[k] = float(((a < 0.05) & (np.abs(l2) > 0.5)).sum())
            nrm[k] = float(np.linalg.norm(np.asarray(delta[i], float)[nonp]))
        feats[f"nde_{c}"], feats[f"nde1_{c}"], feats[f"nrm_{c}"] = nde, nde1, nrm
        per["nde"].append(nde); per["nde1"].append(nde1); per["nrm"].append(nrm)
    with np.errstate(all="ignore"):
        feats["lognde_mean"] = np.nanmean(np.log1p(np.vstack(per["nde"])), 0)
        feats["lognde1_mean"] = np.nanmean(np.log1p(np.vstack(per["nde1"])), 0)
        feats["lognrm_mean"] = np.nanmean(np.log(np.vstack(per["nrm"])), 0)
    return feats


def spearman(u, v):
    ok = np.isfinite(u) & np.isfinite(v)
    return float(np.corrcoef(np.argsort(np.argsort(u[ok])), np.argsort(np.argsort(v[ok])))[0, 1]), int(ok.sum())


def cosagg(x, a, d, s):
    return float(np.mean([(s * x[j]).sum() / np.sqrt((s * s * a[j]).sum() * d.sum()) for j in range(3)]))


def main():
    out = {}
    for f in FOLDS:
        z = np.load(W540 / f"arrays540-{f}.npz")
        x, a, d = z[f"{BASE}__x"], z[f"{BASE}__a"], z["d"]
        targets = [str(t) for t in np.load(WORK / f"dump518-{f}.npz")["targets"]]
        u = (x.sum(0) / a.sum(0)) / (x.sum() / a.sum())
        F = cd4_features(targets)
        base = cosagg(x, a, d, np.ones(len(targets)))
        rec = {}
        for k, v in F.items():
            sp, nok = spearman(u, v)
            ok = np.isfinite(v)
            tl = {}
            for be in (-0.5, -0.25, 0.25, 0.5):
                if k.startswith("log"):
                    med, sd = np.nanmedian(v[ok]), np.nanstd(v[ok]) or 1.0
                    m = np.where(ok, np.exp(be * (np.where(ok, v, med) - med) / sd), 1.0)
                else:
                    vv = np.where(ok, np.maximum(v, 0) + 1.0, np.nan)
                    med = np.nanmedian(vv)
                    m = np.where(ok, (np.where(ok, vv, med) / med) ** be, 1.0)
                tl[f"{be:+g}"] = cosagg(x, a, d, m) / base - 1
            rec[k] = dict(spearman=sp, n=nok, tilt=tl)
        out[f] = rec
    keys = list(out[FOLDS[0]])
    print(f"{'CD4 feature':22s} " + " ".join(f"{f[:8]:>9s}" for f in FOLDS) + "   tilt beta -0.5/-0.25/+0.25/+0.5: cancer mean r (kolf2)")
    for k in keys:
        sp = " ".join(f"{out[f][k]['spearman']:+9.3f}" for f in FOLDS)
        tl = " ".join(f"{100 * np.mean([out[f][k]['tilt'][b] for f in CANCER]):+.2f}%({100 * out['kolf2'][k]['tilt'][b]:+.1f})"
                      for b in ("-0.5", "-0.25", "+0.25", "+0.5"))
        print(f"{k:22s} {sp}   {tl}   n {out[FOLDS[0]][k]['n']}")
    (WORK / "cd4feat518.json").write_text(json.dumps(out, indent=1, default=float))


if __name__ == "__main__":
    main()
