"""410 — the per-target scale artifact for the REAL A/B/C panel (VCC_PERT_SCALE).

key_p  = whitened shared-energy key (key410.shared_energy_key) over the three production donors
         (k562_gwps, orion_hct116, orion_hek293t), raw bank rows, production operator whitener_union_pruned
         (sha e2f65ecf..., the file the A/B/C whitener stage reproduces), inner products over the operator axis.
w_raw  = clip(key / median(key over carried), lo, hi); targets carried by < 2 donors get 1.0.
The builder multiplies each final effect row by w_raw * c, with c the one constant that preserves the whole-matrix
effect norm over non-panel columns (energy-preserving: the arm moves the per-target DISTRIBUTION of magnitude only).

usage: uv run python experiments/410-pertscale-revisit/build_weights410.py LO HI
writes work/pert_scale_410_lo{LO}_hi{HI}.npz (+ .json sidecar with sha256 and diagnostics)
"""
from __future__ import annotations

import hashlib
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "experiments/410-pertscale-revisit"))
sys.path.insert(0, str(ROOT / "experiments/162-retrodicting-validator"))
import plolo  # noqa: E402
from key410 import shared_energy_key, weights  # noqa: E402

DONORS = ("k562_gwps", "orion_hct116", "orion_hek293t")
OP_PATH = ROOT / "data/donors/whitener_union_pruned.npz"
OP_SHA = "e2f65ecf69b7e0b0b817d02d850c7be07a78b442644a239c4b0e520cda1db123"


def main(lo: float, hi: float) -> None:
    targets = [str(t) for t in pd.read_csv(ROOT / "data/challenge/pert_counts.csv")["target_gene"]]
    assert len(targets) == 300
    wcols = np.asarray(np.load(ROOT / "experiments/162-retrodicting-validator/work/fold-kolf2.npz")["wcols"])
    op = plolo.Operator(OP_PATH, name="union_pruned")
    assert op.sha256 == OP_SHA
    # OCT-22 GUARD (DEF-RUNBOOK 8.3): the operator axis is union(donor coverage) - PANEL. A new panel on the old axis puts
    # panel genes inside wcols and every weight silently comes out on the wrong geometry. Refuse instead.
    in_axis = sorted(set(targets) & {plolo.OFFICIAL[i] for i in wcols})
    if in_axis:
        raise SystemExit(f"{len(in_axis)} panel genes sit on the operator axis (e.g. {in_axis[:5]}): this operator was fitted "
                         "for a different panel. Refit it for this panel (DEF-RUNBOOK) and point OP_PATH/wcols at the refit.")
    key, carried, nd = shared_energy_key(targets, DONORS, op, wcols, cols=None)
    w, med = weights(key, carried, lo, hi)
    tag = f"lo{lo:g}_hi{hi:g}"
    out = ROOT / f"experiments/410-pertscale-revisit/work/pert_scale_410_{tag}.npz"
    np.savez(out, targets=np.asarray(targets), w=w.astype(np.float64), key=np.nan_to_num(key, nan=0.0),
             carried=carried, n_donors=nd.astype(np.int64))
    sha = hashlib.sha256(out.read_bytes()).hexdigest()
    q = lambda x: np.round(np.quantile(x, [0, .1, .25, .5, .75, .9, 1]), 4).tolist()
    side = dict(sha256=sha, lo=lo, hi=hi, donors=DONORS, operator_sha256=OP_SHA, gamma=0.65,
                key_median_carried=med, carried=int(carried.sum()), n_donor_hist={int(k): int((nd == k).sum()) for k in np.unique(nd)},
                key_quantiles_carried=q(key[carried]), frac_key_le0=float((key[carried] <= 0).mean()),
                w_quantiles=q(w), frac_at_floor=float((w[carried] <= lo + 1e-12).mean()),
                frac_at_cap=float((w[carried] >= hi - 1e-12).mean()),
                top10=[(targets[i], round(float(key[i] / med), 2)) for i in np.argsort(-np.nan_to_num(key, nan=-9))[:10]],
                bottom10=[(targets[i], round(float(key[i] / med), 2)) for i in np.argsort(np.nan_to_num(key, nan=9))[:10]])
    out.with_suffix(".json").write_text(json.dumps(side, indent=1))
    print(json.dumps(side, indent=1))


if __name__ == "__main__":
    main(float(sys.argv[1]), float(sys.argv[2]))
