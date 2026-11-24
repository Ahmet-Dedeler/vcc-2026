"""030d — joint (gamma, q) grid: do the two transforms overlap or compose?

Whitening reweights *directions across perturbations*; sharpening reweights *magnitudes within
a perturbation*. If they are attacking the same redundancy from two sides, the joint optimum
will sit at a lower gamma than 0.85 and the gains will not add. If they are orthogonal, both
optima hold and the gains compose.

That distinction decides whether tomorrow's candidate is "whitening plus a second knob" or
"one knob that was mis-parameterised".
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "experiments" / "018-proxy6"))
from proxy6 import pds_cosine   # noqa: E402
import run as proxy_run          # noqa: E402

GAMMAS = [0.0, 0.4, 0.65, 0.85, 1.05]
QS = [1.0, 1.25, 1.5, 2.0, 2.5, 3.0]

folds, common = proxy_run.build_folds("1:1")
print(f"panel {len(common)}\n")
print(f"{'gamma \\ q':>10} " + "".join(f"{q:>9.2f}" for q in QS))
out = {}
for g in GAMMAS:
    row = []
    for q in QS:
        vals = []
        for h, f in folds.items():
            P = f["W"].apply(f["P"], g)
            if q != 1.0:
                P = np.sign(P) * np.abs(P) ** q
            vals.append(pds_cosine(P, f["T"]))
        v = float(np.mean(vals))
        out[f"g={g:g} q={q:g}"] = v
        row.append(v)
    print(f"{g:>10.2f} " + "".join(f"{v:9.4f}" for v in row), flush=True)

best = max(out, key=out.get)
print(f"\nbest: {best} -> {out[best]:.4f}")
print(f"whitening alone (g=0.85,q=1): {out['g=0.85 q=1']:.4f}")
print(f"sharpening alone (g=0,q=2.5): {out['g=0 q=2.5']:.4f}")
print(f"neither (g=0,q=1):            {out['g=0 q=1']:.4f}")
Path("experiments/030-sharpen/joint.json").write_text(json.dumps(out, indent=2))
