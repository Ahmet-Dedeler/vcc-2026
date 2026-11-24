#!/usr/bin/env python3
"""410 on H1: the WELDED per-target bulk magnitude through the shipped 395 construction, six members end to end.

Reuses 400's harness (323's reference + cap15 effects + the k3c5bq80 block = CLEAN4 q80 prior, K=3 classes) read-only;
writes only into experiments/410-pertscale-revisit/work/h1/. Arms:
  base395        cap15 effects, k3c5bq80 block -- must reproduce 394's metrics-k3c5bq80.json digests (known answer)
  ps_lo03_hi3    effects scaled per target by the H1 key (key410, H1-axis operator 323), clip(0.3, 3), energy-preserving
  ps_lo03_inf    the same, uncapped
  ps_shuf        ps_lo03_hi3's weights permuted across targets (seed 20260922): the redistribution without the key
  ps_med_lo03_hi3 the same weights WITHOUT energy normalisation (median target unchanged): the DE channel's own response
usage: uv run python run410h.py ARM
"""
from __future__ import annotations

import json
import os
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
WORK = HERE / "work" / "h1"
WORK.mkdir(parents=True, exist_ok=True)
os.environ.setdefault("NUMBA_CACHE_DIR", str(HERE / "work" / "numba-cache"))
sys.path.insert(0, str(ROOT / "experiments/400-uncalled-knn"))
sys.path.insert(0, str(HERE))

import numpy as np  # noqa: E402
import run400 as H  # noqa: E402  -- read-only import, never edited

H.WORK = WORK
_orig_block = H.block_for_arm
H.block_for_arm = lambda arm, *a, **k: _orig_block("k3c5bq80", *a, **k)

OP_H1 = ROOT / "experiments/323-current-family-validation/work/whitener_h1_pruned.npz"
STAGE = ROOT / "experiments/323-current-family-validation/work/stage.npz"
DONORS = ("k562_gwps", "orion_hct116", "orion_hek293t")


class H1Op:
    """Minimal operator (evals/evecs) on the 323 H1 axis; the axis is taken from 323's stage.npz and length-checked."""
    def __init__(self):
        with np.load(OP_H1) as z:
            self.evals = np.asarray(z["evals"], float)
            self.evecs = np.asarray(z["evecs"], float)
        with np.load(STAGE, allow_pickle=False) as z:
            self.cols = np.asarray(z["columns"], dtype=np.int64)
            self.genes = [str(g) for g in z["genes"]]
        assert self.evecs.shape[0] == len(self.cols) == 18157

    def _check_axis(self, wcols):
        assert np.array_equal(np.asarray(wcols), self.cols)


def h1_weights(targets, lo, hi):
    import key410  # noqa: PLC0415
    op = H1Op()
    key, carried, nd = key410.shared_energy_key(list(map(str, targets)), DONORS, op, op.cols, cols=None)
    w, med = key410.weights(key, carried, lo, hi)
    return w, key, carried, nd


def main():
    arm = sys.argv[1]
    b = H.builder389()
    eff, genes, targets = H.native_cap15(b)
    if arm != "base395":
        lo, hi = {"ps_lo03_hi3": (0.3, 3.0), "ps_lo03_inf": (0.3, 1e9), "ps_shuf": (0.3, 3.0), "ps_med_lo03_hi3": (0.3, 3.0)}[arm]
        w, key, carried, nd = h1_weights(targets, lo, hi)
        if arm == "ps_shuf":
            w = np.random.default_rng(20260922).permutation(w)
        panel = set(map(str, targets))
        nonp = np.array([str(g) not in panel for g in genes])
        e2 = np.einsum("ij,ij->i", eff[:, nonp], eff[:, nonp])
        c = 1.0 if arm == "ps_med_lo03_hi3" else float(np.sqrt(e2.sum() / (w * w * e2).sum()))
        final = w * c
        eff = eff * final[:, None]
        rec = dict(arm=arm, lo=lo, hi=hi, carried=int(carried.sum()), n_donor_hist={int(k): int((nd == k).sum()) for k in np.unique(nd)},
                   key_quantiles=np.quantile(key[carried], [0, .1, .25, .5, .75, .9, 1]).tolist(),
                   final_quantiles=np.quantile(final, [0, .1, .25, .5, .75, .9, 1]).tolist(), energy_constant=c,
                   weights=dict(zip(map(str, targets), map(float, final))))
        (WORK / f"weights-{arm}.json").write_text(json.dumps(rec, indent=1))
        print(json.dumps({k: v for k, v in rec.items() if k != "weights"}), flush=True)
    receipt = H.emit_arm(arm, eff, genes)
    if arm == "base395":
        exp = json.loads((ROOT / "experiments/394-protected-block-signs/work/metrics-k3c5bq80.json").read_text())
        ok = (receipt["emitted_count_digest_sha256"] == exp["emitted_count_digest_sha256"]
              and receipt["panel_total_digest_sha256"] == exp["panel_total_digest_sha256"])
        diffs = {m: abs(receipt["raw"][m] - exp["raw"][m]) for m in exp["raw"]}
        (WORK / "known-answer-base395.json").write_text(json.dumps(dict(pass_=ok, diffs=diffs), indent=1))
        print("KNOWN ANSWER", "PASS" if ok else "FAIL", diffs, flush=True)
        if not ok:
            raise SystemExit("base395 did not reproduce 394's k3c5bq80 digests")


if __name__ == "__main__":
    main()
