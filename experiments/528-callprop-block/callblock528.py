#!/usr/bin/env python3
"""528 — the per-context call-propensity block key artifact for the builder (`VCC_RC_CALLSET_BLOCK=prop:<npz>:<sha256>`),
computed from each context's OWN control cells only. The October form is the same command on the day's D/E/F bundle.

  key_<c>[g] = log10(ref_lin_c[g]) + BETA * x_c[g] * (MEDIAN_X_H1 / median_c(x))
    ref_lin_c : the emitter's control CPM (mean per-cell CPM over the context's non-targeting cells; the CPM block's key)
    x_c       : excess dropout = Poisson-expected detection (the gene's pooled fraction, 10 control-library bins) minus the
                observed detection rate (517's rule; feat528.features)
    median_c  : the median of x_c over the context's 5-CPM universe; MEDIAN_X_H1 = h1-test's (0.0087745...), so on H1 the
                key is exactly the priced `lcap` key (lc + 22 x). The median rescale makes one beta transfer across lines whose
                dropout scale differs (fold528b: it moved all three public folds the right way vs the raw form).
  BETA = 22: the value maximising the capture of h1-trainval's 126 perturbations' DE genes (keys528.py; no h1-test truth).
Genes outside a context's 5-CPM universe keep finite values (the builder only reads the universe: rc_control.keep).

    PYTHONDONTWRITEBYTECODE=1 uv run python experiments/528-callprop-block/callblock528.py \
        --bundle data/challenge --contexts A,B,C --out experiments/528-callprop-block/work/callblock528_abc.npz
KA (--ka): the H1 key from h1-test's controls gives exactly the priced `lcap` blocks on all 85 H1 targets (the float key is
within 1e-15 of keys528's `lcap`). The artifact is deterministic (same bytes, same sha256, on a rerun).
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import sys
from pathlib import Path

sys.dont_write_bytecode = True
import numpy as np  # noqa: E402

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(HERE))
from feat528 import features  # noqa: E402

BETA = 22.0
# h1-test's median excess dropout over its 5-CPM universe (feat528 on its 18,400 controls), pinned so the October form needs
# no H1 file; checked against the cached features whenever they exist
MEDIAN_X_H1 = 0.00877448551593013


def median_x_h1(load_features=False):
    p = HERE / "work/feat528_H1.npz"
    F = None
    if p.exists():
        F = dict(np.load(p, allow_pickle=True))
        got = float(np.median(F["x"][F["keep"]]))
        if got != MEDIAN_X_H1:
            raise SystemExit(f"feat528_H1's median x {got!r} != the pinned {MEDIAN_X_H1!r}")
    elif load_features:
        F = features(ROOT / "data/harness/h1-test.h5ad")
    return MEDIAN_X_H1, F


def key_from_features(F, med_h1, beta=BETA):
    med = float(np.median(F["x"][F["keep"]]))
    lc = np.log10(np.maximum(F["mean_cpm"], 1e-3))
    return lc + beta * F["x"] * (med_h1 / med), med


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--bundle", default=str(ROOT / "data/challenge"))
    ap.add_argument("--contexts", default="A,B,C")
    ap.add_argument("--out", default=str(HERE / "work/callblock528_abc.npz"))
    ap.add_argument("--ka", action="store_true")
    a = ap.parse_args()
    med_h1, FH = median_x_h1(load_features=a.ka)
    if a.ka:
        z = np.load(ROOT / "experiments/323-current-family-validation/work/reference.npz", allow_pickle=False)
        keep = np.asarray(z["keep"], bool)
        k, med = key_from_features(FH, med_h1)
        lcap = np.load(HERE / "work/keys528.npz")["lcap"][0]
        d = float(np.max(np.abs(k[keep] - lcap)))
        # the float key differs by ulps (feat528's streamed mean CPM vs the reference's emit_ref_lin, 1.7e-15 relative);
        # what the emitter reads is the ORDER: the per-target blocks on H1 must be identical for all 85 targets
        sys.path.insert(0, str(ROOT / "experiments/521-strength-nreal"))
        import lib521 as L  # noqa: PLC0415
        from k_emit528 import block_from_key  # noqa: PLC0415
        S = L.load_state(0)
        same = all(np.array_equal(block_from_key(S["head"][i], S["ref_lin"], 2200, k[keep], 0),
                                  block_from_key(S["head"][i], S["ref_lin"], 2200, lcap, 0)) for i in range(85))
        print(f"KA: H1 key from h1-test controls vs keys528 lcap: max |diff| {d:.3e}; blocks identical on 85/85 targets: "
              f"{same}; median x {med:.6f} == {med_h1:.6f}")
        if d > 1e-12 or not same:
            raise SystemExit("KA FAILED")
        return
    bundle = Path(a.bundle)
    official = [r["gene_name"] for r in csv.DictReader(open(bundle / "gene_names.csv"))]
    out, meta = {"genes": np.array(official)}, {"beta": BETA, "median_x_h1": med_h1, "contexts": {}}
    for c in [s.strip() for s in a.contexts.split(",") if s.strip()]:
        path = bundle / f"context_{c}.h5ad"
        F = features(path)
        if list(map(str, F["genes"])) != official:
            raise SystemExit(f"context_{c}: gene order differs from gene_names.csv")
        k, med = key_from_features(F, med_h1)
        if not np.all(np.isfinite(k)):
            raise SystemExit(f"context_{c}: non-finite key")
        out[f"key_{c}"] = k
        e = F["keep"]
        # the swap this key makes at the default block size (no per-target head: approximate)
        cpm_top = set(np.argsort(-np.where(e, F["mean_cpm"], -np.inf), kind="stable")[:2200].tolist())
        key_top = set(np.argsort(-np.where(e, k, -np.inf), kind="stable")[:2200].tolist())
        meta["contexts"][c] = dict(n_ctrl=int(F["n_ctrl"]), n_universe=int(e.sum()), median_x=med, scale=med_h1 / med,
                                   swapped_vs_cpm_block_no_head=len(key_top - cpm_top),
                                   control_file_sha256=hashlib.sha256(path.read_bytes()).hexdigest())
        print(c, json.dumps(meta["contexts"][c]), flush=True)
    out["meta"] = np.array(json.dumps(meta))
    np.savez(a.out, **out)
    sha = hashlib.sha256(Path(a.out).read_bytes()).hexdigest()
    Path(a.out + ".json").write_text(json.dumps({**meta, "sha256": sha, "path": str(a.out)}, indent=1))
    print(f"wrote {a.out} sha256 {sha}")


if __name__ == "__main__":
    main()
