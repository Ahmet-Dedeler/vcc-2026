"""Known-answer gate for 359: end-to-end replay of 341's cross-line diagonal.

Reproduces pooled weighted sign .589398 and pooled retrieval .794614
(341 RESULTS.md) with own driver code over the same deterministic path,
and checks bit-identity against 341's saved crossline_diag predictions.
Writes nothing outside this directory (results to work/ka359.json).
"""

import json
import os
import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
PRIOR = ROOT / "experiments/334-text-conditioned-transfer"
EXP341 = ROOT / "experiments/341-known-line-calibration"
sys.path.insert(0, str(PRIOR))
sys.path.insert(0, str(EXP341))

os.environ.setdefault("NUMBA_NUM_THREADS", "4")
os.environ.setdefault("OMP_NUM_THREADS", "4")

import run334 as B  # noqa: E402
import core334 as P  # noqa: E402
import core341 as C  # noqa: E402

EXP_W_SIGN = 0.589398
EXP_RETRIEVAL = 0.794614


def main():
    t0 = time.monotonic()
    selection = json.loads((PRIOR / "selection.json").read_text())
    targets, genes = selection["targets"], selection["genes"]
    assert len(targets) == 1024 and len(genes) == 512
    assert len(set(targets) & set(genes)) == 0

    checks = B.known()
    print("historical checks:", [(r["line"], r["maximum_historical_difference"]) for r in checks])

    loaded = [B.load(bank, targets, genes) for bank in B.BANKS]
    raw = np.array([a for a, _ in loaded])
    basal = np.array([b for _, b in loaded])
    print("raw", raw.shape, "basal", basal.shape)
    # basal medians over the 512 frame (threshold inputs for the gate arms)
    for line, b in zip(B.LINES, basal):
        print(f"basal {line}: median {float(np.median(b)):.4f} "
              f"IQR [{float(np.quantile(b, .25)):.4f},{float(np.quantile(b, .75)):.4f}] "
              f"min {float(b.min()):.4f} max {float(b.max()):.4f}")
    direction, rms = P.rms_direction(raw)
    del loaded, raw

    order = P.order(targets, "334-split|")
    blocks = [np.array(order[::2]), np.array(order[1::2])]

    w_all, r_all, per_line = [], [], {}
    max_saved_diff = 0.0
    for outer, line in enumerate(B.LINES):
        donors = [i for i in range(6) if i != outer]
        baseline = direction[donors].mean(0)
        truth = direction[outer]
        full = np.full_like(baseline, np.nan)
        for block in (0, 1):
            train, test = blocks[1 - block], blocks[block]
            inner = P.source_blocks(direction, basal, rms, outer, train)
            gain = P.diagonal(inner)
            full[test] = baseline[test] * gain
        assert np.isfinite(full).all()
        w = C.sign_endpoints(full, truth)["weighted_sign"]
        r = P.endpoints(full, truth)["retrieval"]
        w_all.append(w)
        r_all.append(r)
        per_line[line] = {"weighted_sign": float(w.mean()), "retrieval": float(r.mean())}
        # bit-identity against 341's saved crossline_diag
        with np.load(EXP341 / "work" / line / "predictions.npz") as z:
            saved = np.asarray(z["crossline_diag"])
            assert list(z["targets"]) == targets
            max_saved_diff = max(max_saved_diff, float(np.abs(saved - full).max()))

    w_all = np.array(w_all)
    r_all = np.array(r_all)
    mean_w, mean_r = float(w_all.mean()), float(r_all.mean())
    print(f"pooled weighted_sign {mean_w:.12f} (expected {EXP_W_SIGN})")
    print(f"pooled retrieval     {mean_r:.12f} (expected {EXP_RETRIEVAL})")
    print(f"max |ours - 341 saved crossline_diag| = {max_saved_diff:.3e}")
    for line, d in per_line.items():
        print(f"  {line}: wsign {d['weighted_sign']:.6f} retrieval {d['retrieval']:.6f}")

    dw, dr = abs(mean_w - EXP_W_SIGN), abs(mean_r - EXP_RETRIEVAL)
    ok = dw < 5e-7 and dr < 5e-7 and max_saved_diff == 0.0
    print(f"|dw|={dw:.3e} |dr|={dr:.3e} saved_exact={max_saved_diff == 0.0} -> {'PASS' if ok else 'FAIL'}")

    (HERE / "work").mkdir(parents=True, exist_ok=True)
    (HERE / "work" / "ka359.json").write_text(json.dumps({
        "pooled_weighted_sign": mean_w, "pooled_retrieval": mean_r,
        "dw": dw, "dr": dr, "max_saved_diff": max_saved_diff,
        "per_line": per_line, "gate": "PASS" if ok else "FAIL",
        "seconds": time.monotonic() - t0,
    }, indent=2) + "\n")
    if not ok:
        raise SystemExit("KNOWN-ANSWER GATE FAILED — stop, no arms.")


if __name__ == "__main__":
    main()
