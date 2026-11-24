"""Pick alpha for v5 so that emitted effect magnitude is unchanged from v4b.

Experiment 017/019 change the *direction* of the prediction (equal donor weights plus
spectral whitening). Magnitude is a separate lever that moves `nmae`, `fid` and `jac`, and
mixing the two in one submission would make the result unattributable -- which is exactly how
v3 was lost. So alpha is re-derived here to hold the emitted median effect norm at v4b's
value, leaving direction as the single change.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "scripts"))
from build_submission import ChallengeSpec          # noqa: E402
from vcc2026.donor_transfer import Donor            # noqa: E402
from vcc2026.multi_donor import MultiDonorModel     # noqa: E402
from vcc2026.whiten import Whitener                 # noqa: E402

spec = ChallengeSpec.load("data/challenge")
paths = {"k562_gwps": "data/donors/k562_gwps",
         "orion_hct116": "data/donors/orion_hct116",
         "orion_hek293t": "data/donors/orion_hek293t"}
donors = {n: Donor.load(Path(p)) for n, p in paths.items()}
cpp = {n: json.loads((Path(p) / "meta.json").read_text())["cells_per_pert"]["median"]
       for n, p in paths.items()}
model = MultiDonorModel(donors, spec.genes, gene_norm=False, eb=False, cells_per_pert=cpp)

for label, w, a in [("v4b 4:1:1", {"k562_gwps": 4, "orion_hct116": 1, "orion_hek293t": 1}, 0.92),
                    ("v5  1:1:1", {"k562_gwps": 1, "orion_hct116": 1, "orion_hek293t": 1}, 1.0)]:
    E = np.vstack([model.effect(p, a, 0.0, weights=w) for p in spec.perts])
    print(f"{label} alpha={a}: median effect norm {np.median(np.linalg.norm(E, axis=1)):.4f}")

# Whitening axis and operator, then the magnitude it produces before renormalisation.
shared = set(model.target_genes[model.tcols[model.names[0]]])
for n in model.names[1:]:
    shared &= set(model.target_genes[model.tcols[n]])
shared -= set(spec.perts)
tindex = {g: i for i, g in enumerate(model.target_genes)}
wcols = np.array(sorted(tindex[g] for g in shared), dtype=np.int64)
print(f"whitening axis: {len(wcols)} genes (donor intersection, panel targets removed)")

cache = Path("data/donors/whitener.npz")
if cache.exists():
    W = Whitener.load(cache)
else:
    rng = np.random.default_rng(0)
    rows = []
    for n in model.names:
        S = np.zeros((model.s[n].shape[0], len(model.target_genes)), np.float32)
        S[:, model.tcols[n]] = model.s[n][:, model.dcols[n]]
        if S.shape[0] > 6000:
            S = S[np.sort(rng.choice(S.shape[0], 6000, False))]
        rows.append(model._renorm[n] * S[:, wcols].astype(np.float64))
    B = np.vstack(rows)
    print(f"fitting whitener on {B.shape} ...", flush=True)
    W = Whitener.fit(B)
    W.save(cache)
    cache.with_suffix(".json").write_text(json.dumps(
        {"n_cols": int(len(wcols)), "basis": list(B.shape), "donors": model.names}))
    del B

w = {"k562_gwps": 1, "orion_hct116": 1, "orion_hek293t": 1}
E = np.vstack([model.effect(p, 1.0, 0.0, weights=w) for p in spec.perts])
before = float(np.median(np.linalg.norm(E, axis=1)))
E[:, wcols] = W.apply(E[:, wcols], 0.65, 0.0, renorm=False)
after = float(np.median(np.linalg.norm(E, axis=1)))
print(f"1:1:1 alpha=1 pre-whiten {before:.4f} -> post-whiten {after:.4f} "
      f"(renormalised back to {before:.4f})")
