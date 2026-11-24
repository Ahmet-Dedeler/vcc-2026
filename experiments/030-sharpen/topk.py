"""030c — sharpening is soft top-k selection, not distribution matching.

The magnitude-matching explanation is dead: the predicted distribution is barely flatter than
the truth's (sd of log|effect| 1.138 vs 1.196, matching exponent **q* = 1.05**) while the
PDS optimum sits at q = 2.0-2.5. At q = 2 the prediction is far *spikier* than the truth and
PDS improves anyway, so the transform is not making the prediction look more like the target.

The remaining explanation is selection. Raising |e| to a power suppresses the long tail of
small per-gene entries, which are the least reliably estimated part of the donor effect. That
is a **soft top-k**, and it has a hard counterpart worth testing directly — not least because
`Jessica A` sits on the board with the model string `K562 top128 transfer v1` and the only raw
`mse` below 1.

This matters beyond one knob, because it resolves an apparent contradiction in the notebook:

  - experiment 004: truncating the donor **SVD** loses monotonically -- the low-variance
    *directions across perturbations* carry signal
  - here: suppressing small **per-gene entries within a perturbation** gains a lot

Both are true and they are not in conflict. A low-variance SVD direction is a coherent pattern
shared across perturbations; a small per-gene entry is one number in one prediction. The first
is signal spread thin, the second is noise.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "experiments" / "018-proxy6"))
from proxy6 import pds_cosine, topk_mask, nmae_curve, scale   # noqa: E402
import run as proxy_run                                        # noqa: E402

GAMMA, K_REAL = 0.85, 300
ALPHAS = np.array([0.0, 0.1, 0.2, 0.3, 0.5, 0.7, 0.92, 1.2, 1.6, 2.2])
J_OP = int(np.argmin(np.abs(ALPHAS - 0.92)))

folds, common = proxy_run.build_folds("1:1")
n_genes = list(folds.values())[0]["P"].shape[1]
print(f"panel {len(common)}, {n_genes} scored genes, gamma={GAMMA}\n")

configs = [("full", None, 1.0)]
configs += [(f"top{k}", k, 1.0) for k in (128, 300, 1000, 2000, 4000)]
configs += [(f"q={q:g}", None, q) for q in (1.5, 2.0, 2.5)]
configs += [(f"top1000 q=1.5", 1000, 1.5)]

print(f"{'variant':16} {'pds':>8} {'s_pds':>7} {'nmae@.92':>9} {'s_nmae':>7} {'sum':>8}")
out = {}
for name, k, q in configs:
    acc = {"pds": [], "nmae": []}
    for h, f in folds.items():
        P = f["W"].apply(f["P"], GAMMA)
        if q != 1.0:
            P = np.sign(P) * np.abs(P) ** q
        if k is not None:
            P = P * topk_mask(P, k)
        T = f["T"]
        R = topk_mask(T, K_REAL)
        acc["pds"].append(pds_cosine(P, T))
        Pn = P / (np.median(np.linalg.norm(P, axis=1)) + 1e-12)
        Tn = T / (np.median(np.linalg.norm(T, axis=1)) + 1e-12)
        acc["nmae"].append(nmae_curve(Pn, Tn, R, ALPHAS)[J_OP])
    pv, nv = float(np.mean(acc["pds"])), float(np.mean(acc["nmae"]))
    sp, sn = scale("pds", pv), scale("nmae", nv)
    out[name] = dict(pds=pv, nmae=nv, s_pds=sp, s_nmae=sn, sum=sp + sn)
    print(f"{name:16} {pv:8.4f} {sp:7.3f} {nv:9.4f} {sn:7.3f} {sp + sn:8.4f}", flush=True)

Path("experiments/030-sharpen/topk.json").write_text(json.dumps(out, indent=2))
best = max(out, key=lambda k: out[k]["pds"])
print(f"\nbest pds: {best} -> {out[best]['pds']:.4f}")
