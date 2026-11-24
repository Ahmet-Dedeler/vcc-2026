#!/usr/bin/env python3
"""700: is PIE's CONTEXT-SPECIFIC part informative? This decides whether a per-context PIE blend (672) can beat the shared
one. Molina & Zhang (2026) report that the context x perturbation interaction is not zero-shot predictable by any model they
tested.

Surface: Replogle-Nadig essential (K562, RPE1, Jurkat, HepG2). PIE's xdataset checkpoint never trained on it and pooled no
evidence from it. Perturbations: those in PIE's replogle_xdataset test_seen split measured in >= 2 of the 4 lines (679).
PIE predicted all 4 contexts for each (work/query_ctxspec.json, infer700.py with PIE's own Replogle dir).
Truth: the dir's delta_p (mean log-normalised expression shift, every measured gene) per measured (context, pert).
Per perturbation with >= 2 measured contexts:
  - interaction cosine: cos(PIE_c - mean_c' PIE_c', truth_c - mean_c' truth_c') over the measured contexts;
  - own-vs-other: is cos(PIE_c, truth_c) > mean over c' != c of cos(PIE_c', truth_c)? (the right context's prediction wins)
Reports the means over perturbations with a sign-flip null (context labels permuted within each perturbation).
usage: .venv/bin/python ctxspec700.py <pred_ctxspec.parquet>
"""
from __future__ import annotations

import glob
import json
import sys
from pathlib import Path

import numpy as np
import pyarrow.parquet as pq

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]


def main():
    d = glob.glob(str(ROOT / "data/external/pie/data/hf/datasets/arcinstitute/pie_replogle_nadig_essential/*/preprocessed"))[0]
    meta = json.load(open(f"{d}/meta.json"))
    dp = np.load(f"{d}/delta_p.npy", mmap_mode="r")
    ctx_ids = np.load(f"{d}/ctx_ids.npy"); pert_ids = np.load(f"{d}/pert_ids.npy")
    id2ctx = {v: k for k, v in meta["context_to_id"].items()}; id2pert = {v: k for k, v in meta["pert_to_id"].items()}
    truth = {}
    for r in range(len(ctx_ids)):
        truth[(id2ctx[int(ctx_ids[r])], id2pert[int(pert_ids[r])])] = r
    t = pq.read_table(sys.argv[1]); df = t.to_pandas()
    P = {(r["context"], r["perturbation"]): np.asarray(r["delta_p_pred"], np.float64) for _, r in df.iterrows()}
    L = {(r["context"], r["perturbation"]): np.asarray(r["lfc_pred"], np.float64) for _, r in df.iterrows()}
    ctxs = sorted({c for c, _ in P}); perts = sorted({p for _, p in P})
    rng = np.random.default_rng(700)

    def cos(a, b):
        na, nb = np.linalg.norm(a), np.linalg.norm(b)
        return float(a @ b / (na * nb)) if na > 0 and nb > 0 else np.nan

    inter, inter_null, own_win, n_used = [], [], [], 0
    for p in perts:
        cs = [c for c in ctxs if (c, p) in truth and (c, p) in P]
        if len(cs) < 2:
            continue
        n_used += 1
        T = np.stack([np.asarray(dp[truth[(c, p)]], np.float64) for c in cs])
        Q = np.stack([P[(c, p)] for c in cs])
        Tc, Qc = T - T.mean(0), Q - Q.mean(0)
        inter.append(np.nanmean([cos(Qc[i], Tc[i]) for i in range(len(cs))]))
        perm = rng.permutation(len(cs))
        inter_null.append(np.nanmean([cos(Qc[perm[i]], Tc[i]) for i in range(len(cs))]))
        for i in range(len(cs)):
            own = cos(Q[i], T[i])
            other = np.nanmean([cos(Q[j], T[i]) for j in range(len(cs)) if j != i])
            own_win.append(own > other)
    res = {"perturbations_used": n_used, "contexts": ctxs,
           "interaction_cos_mean": float(np.nanmean(inter)), "interaction_cos_null_mean": float(np.nanmean(inter_null)),
           "own_context_wins_frac": float(np.mean(own_win)), "n_own_tests": len(own_win)}
    (HERE / "work/ctxspec700.json").write_text(json.dumps(res, indent=1))
    print(json.dumps(res, indent=1))


if __name__ == "__main__":
    main()
