"""432 — 410's per-target bulk key (VCC_PERT_SCALE `wshare_lo0.3_hi3`) rebuilt for ANY panel.

`experiments/410-pertscale-revisit/build_weights410.py` is the production recipe, but it hard-codes three A/B/C things:
  1. the panel: `data/challenge/pert_counts.csv` (no --bundle);
  2. the operator: `OP_PATH = data/donors/whitener_union_pruned.npz`, pinned to e2f65ecf..., and its axis `wcols`
     from `experiments/162-retrodicting-validator/work/fold-kolf2.npz` (the A/B/C axis);
  3. the OUTPUT path: `experiments/410-pertscale-revisit/work/pert_scale_410_lo{LO}_hi{HI}.npz` -- the very file the
     builder's `PERT_SCALE_ARTIFACTS["wshare_lo0.3_hi3"]` points at. Run as-is on Oct 22 (after the bundle swap) it would
     OVERWRITE the A/B/C artifact in place (the builder then refuses it on sha, so fail-closed, but the A/B/C file is gone).

This script is the same arithmetic (key410.shared_energy_key + key410.weights, the same npz fields in the same order),
with the panel, the operator and the output as arguments, and the axis recomputed from the bundle exactly as
`scripts/def_pipeline.py:bundle_axis` / `build_union_whitener.py` form it: sorted official indices of
union(production donors' genes.txt) - panel.

KNOWN ANSWER (--ka): on data/challenge with the production operator it must reproduce the pinned A/B/C artifact
`pert_scale_410_lo0.3_hi3.npz` byte for byte (sha d7f76a22...), and the recomputed axis must equal fold-kolf2's wcols.

usage: uv run python key432.py --bundle <dir> --op <whitener.npz> --lo 0.3 --hi 3 --out <out.npz> [--ka]
"""
from __future__ import annotations

import argparse
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
ABC_OP_SHA = "e2f65ecf69b7e0b0b817d02d850c7be07a78b442644a239c4b0e520cda1db123"
ABC_ARTIFACT_SHA = {(0.3, 3.0): "d7f76a22f62ad6991d7220ce8a65218db38b2d71dcfe3c3781ec6a8f70409e0e"}


def bundle_axis(bundle: Path) -> np.ndarray:
    official = pd.read_csv(bundle / "gene_names.csv")["gene_name"].astype(str).tolist()
    if tuple(official) != tuple(plolo.OFFICIAL):
        raise SystemExit("bundle gene axis differs from plolo.OFFICIAL (data/challenge/gene_names.csv)")
    panel = set(pd.read_csv(bundle / "pert_counts.csv")["target_gene"].astype(str))
    tindex = {g: i for i, g in enumerate(official)}
    union: set[str] = set()
    for donor in DONORS:
        union |= {g for g in (ROOT / "data/donors" / donor / "genes.txt").read_text().split("\n") if g in tindex}
    return np.asarray(sorted(tindex[g] for g in union - panel), dtype=np.int64)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--bundle", type=Path, required=True)
    ap.add_argument("--op", type=Path, required=True)
    ap.add_argument("--lo", type=float, default=0.3)
    ap.add_argument("--hi", type=float, default=3.0)
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--ka", action="store_true")
    a = ap.parse_args()
    bundle = a.bundle.resolve()
    targets = [str(t) for t in pd.read_csv(bundle / "pert_counts.csv")["target_gene"]]
    assert len(targets) == 300
    wcols = bundle_axis(bundle)
    op = plolo.Operator(a.op)          # registered (A/B/C) or self-described by its own json (a day-of refit)
    if a.ka:
        assert op.sha256 == ABC_OP_SHA, f"KA needs the production operator, got {op.sha256}"
        kolf = np.asarray(np.load(ROOT / "experiments/162-retrodicting-validator/work/fold-kolf2.npz")["wcols"])
        assert np.array_equal(kolf.astype(np.int64), wcols), "KA: recomputed axis != fold-kolf2 wcols"
        print("KA axis: recomputed wcols == fold-kolf2 wcols", flush=True)
    in_axis = sorted(set(targets) & {plolo.OFFICIAL[i] for i in wcols})
    if in_axis:
        raise SystemExit(f"{len(in_axis)} panel genes on the operator axis: wrong operator/panel pairing")
    op._check_axis(wcols)               # refuses unless OFFICIAL[wcols] hashes to the operator's axis record
    print(f"operator {op.sha256[:16]} axis check: {op.axis_checked}", flush=True)
    key, carried, nd = shared_energy_key(targets, DONORS, op, wcols, cols=None)
    w, med = weights(key, carried, a.lo, a.hi)
    a.out.parent.mkdir(parents=True, exist_ok=True)
    np.savez(a.out, targets=np.asarray(targets), w=w.astype(np.float64), key=np.nan_to_num(key, nan=0.0),
             carried=carried, n_donors=nd.astype(np.int64))
    sha = hashlib.sha256(a.out.read_bytes()).hexdigest()
    q = lambda x: np.round(np.quantile(x, [0, .1, .25, .5, .75, .9, 1]), 4).tolist()
    side = dict(sha256=sha, lo=a.lo, hi=a.hi, donors=DONORS, bundle=str(bundle), operator=str(a.op.resolve()),
                operator_sha256=op.sha256, operator_axis=op.axis_checked, n_wcols=int(wcols.size), gamma=0.65,
                key_median_carried=med, carried=int(carried.sum()),
                n_donor_hist={int(k): int((nd == k).sum()) for k in np.unique(nd)},
                key_quantiles_carried=q(key[carried]), frac_key_le0=float((key[carried] <= 0).mean()),
                w_quantiles=q(w), frac_at_floor=float((w[carried] <= a.lo + 1e-12).mean()),
                frac_at_cap=float((w[carried] >= a.hi - 1e-12).mean()),
                top10=[(targets[i], round(float(key[i] / med), 2)) for i in np.argsort(-np.nan_to_num(key, nan=-9))[:10]],
                bottom10=[(targets[i], round(float(key[i] / med), 2)) for i in np.argsort(np.nan_to_num(key, nan=9))[:10]])
    a.out.with_suffix(".json").write_text(json.dumps(side, indent=1))
    print(json.dumps(side, indent=1))
    if a.ka:
        want = ABC_ARTIFACT_SHA[(a.lo, a.hi)]
        print(f"KA {'PASS' if sha == want else 'FAIL'}: {sha} vs pinned {want}", flush=True)
        if sha != want:
            raise SystemExit(1)


if __name__ == "__main__":
    main()
