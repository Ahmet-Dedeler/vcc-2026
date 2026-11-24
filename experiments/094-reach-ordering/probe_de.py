"""094 — verify on the REAL challenge controls that an arm moved the ordering.

`reach` is decided by the order the scorer walks our genes in, and that order is
set by the Wilcoxon p-value against each context's real control cells. This
recomputes exactly that, per arm, so an intervention can be checked BEFORE a
submission slot is spent on it. No ground truth is needed or used: every quantity
here is a property of our own prediction and the real controls.

Usage:  python probe_de.py <pred.h5ad> <label> [n_perts_per_context]
"""
from __future__ import annotations
import json, os, sys
import numpy as np
import anndata as ad
import h5py
from scipy.stats import false_discovery_control

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                "..", "092-call-set-truncation"))
from de_kernel import ctrl_tiesums, wilcoxon_z  # noqa: E402

PRED, LABEL = sys.argv[1], sys.argv[2]
NPERT = int(sys.argv[3]) if len(sys.argv) > 3 else 12
BUNDLE = os.environ.get("BUNDLE", "data/challenge")
ALPHA, MIN_CPM_SCORER = 0.05, 5.0

def norm_log(x):
    """CPM-normalize each cell to 1e6 then log1p — the scorer's transform."""
    s = np.asarray(x.sum(1)).ravel()
    s[s == 0] = 1.0
    return np.log1p(x * (1e6 / s[:, None]))

pred = ad.read_h5ad(PRED, backed="r")
obs = pred.obs
pcol = "target_gene" if "target_gene" in obs else obs.columns[0]
print(f"[{LABEL}] pred {pred.shape}  pert column {pcol!r}")
rng = np.random.default_rng(0)
out = {}

for ctx in ("A", "B", "C"):
    real = ad.read_h5ad(f"{BUNDLE}/context_{ctx}.h5ad")
    X = real.X
    ctrl = norm_log(np.asarray(X.todense()) if hasattr(X, "todense") else np.asarray(X))
    cpm = np.expm1(ctrl).mean(0)
    keep = cpm > MIN_CPM_SCORER
    cs = np.ascontiguousarray(np.sort(ctrl[:, keep].T, axis=1))
    tie = ctrl_tiesums(cs)
    # cell_eval2 takes the ARITHMETIC mean of per-cell CPM, not expm1 of the mean log1p --
    # the latter is a geometric mean, is systematically smaller, and biases every lfc up.
    # n_pred is unaffected either way because it comes from the rank test, not the lfc.
    ref_lin = np.expm1(ctrl[:, keep]).mean(0)
    del ctrl, real, X
    print(f"[{LABEL}] ctx {ctx}: {cs.shape[1]} controls, {keep.sum()} genes past the 5-CPM gate")

    labels = obs[pcol].astype(str).values
    ctxcol = [c for c in obs.columns if obs[c].astype(str).isin(["A","B","C"]).all()]
    mask = (obs[ctxcol[0]].astype(str).values == ctx) if ctxcol else np.ones(len(obs), bool)
    tg = sorted(set(labels[mask]) - {"non-targeting", "NTC", "control"})
    pick = rng.choice(len(tg), min(NPERT, len(tg)), replace=False)
    rows = []
    with h5py.File(PRED, "r") as h:
        for t in (tg[i] for i in pick):
            sel = np.flatnonzero(mask & (labels == t))
            blk = norm_log(pred[sel].to_memory().X.toarray().astype(np.float64))
            z = wilcoxon_z(cs, tie, np.ascontiguousarray(blk[:, keep].T))
            from scipy.stats import norm as _n
            p = 2 * _n.sf(np.abs(z)); q = false_discovery_control(p, method="bh")
            lfc = np.log2((np.expm1(blk[:, keep]).mean(0) + 1e-9) / (ref_lin + 1e-9))
            sig = q < ALPHA
            r = np.argsort(np.argsort(q))          # 0 = most confident
            dec = (r * 10 // r.size).clip(0, 9)
            rows.append(dict(
                target=t, n_pred=int(sig.sum()),
                frac_up_call=float((lfc[sig] > 0).mean()) if sig.any() else float("nan"),
                up_by_dec=[float((lfc[dec == k] > 0).mean()) for k in range(10)],
                cpm_by_dec=[float(cpm[keep][dec == k].mean()) for k in range(10)],
            ))
            del blk
    A = np.array([x["up_by_dec"] for x in rows])
    C = np.array([x["cpm_by_dec"] for x in rows])
    npred = np.array([x["n_pred"] for x in rows])
    out[ctx] = dict(n=len(rows), n_pred=npred.tolist(),
                    up_by_dec=A.mean(0).tolist(), cpm_by_dec=C.mean(0).tolist(),
                    frac_up_call=float(np.nanmean([x["frac_up_call"] for x in rows])))
    print(f"[{LABEL}] ctx {ctx}: n_pred med {np.median(npred):.0f} "
          f"[{npred.min()}, {npred.max()}]   call set {out[ctx]['frac_up_call']*100:.1f}% UP")
    print(f"[{LABEL}] ctx {ctx}: frac UP by confidence decile (1 = most confident)")
    print("          " + "  ".join(f"{v:.3f}" for v in A.mean(0)))
    print(f"[{LABEL}] ctx {ctx}: mean ctrl CPM by decile")
    print("          " + "  ".join(f"{v:8.1f}" for v in C.mean(0)))
    del cs, tie

dest = f"/work/exp094/probe-{LABEL}.json" if os.path.isdir("/work/exp094") else f"probe-{LABEL}.json"
json.dump(out, open(dest, "w"), indent=2)
print(f"[{LABEL}] wrote {dest}")
