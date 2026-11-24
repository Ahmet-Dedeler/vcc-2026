"""193 Phase 0, known-answer row 2 — reproduce 162's a15 -> a23 mass-fix `pds` price.

169 §2 records this as `+0.004014 scaled` and 162 is the source.  Unlike gate 1 (which only
re-reads the stored `arms2.json`), this recomputes BOTH arms end to end through this
experiment's own code path, so it also exercises `mass_consistent_effect` -- the one part of
the renderer a pure operator swap would otherwise never touch.
"""
from __future__ import annotations

import json
import os
import sys
from pathlib import Path

os.environ.setdefault("OPENBLAS_NUM_THREADS", "8")
os.environ.setdefault("OMP_NUM_THREADS", "8")

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(ROOT / "experiments/162-retrodicting-validator"))
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(HERE))

import plolo                                                        # noqa: E402
from plolo import build_fold, deliver, endpoints, solve_fraction    # noqa: E402
from run_arms import CONFIGS                                        # noqa: E402
from screen193 import FOLD_ORDER, SPAN_PDS, cd4_cached, frames_for, sha256_file  # noqa: E402


def score(fold, cfg, evals, evecs, wcols, cd4, frames, truth):
    pre = np.asarray(fold["pre"], dtype=np.float64) * cfg["alpha"]
    sub = pre[:, wcols]
    before = float(np.median(np.linalg.norm(sub, axis=1))) or 1.0
    w = (evals + cfg["lam"]) ** (-cfg["gamma"] / 2.0)
    t = ((sub @ evecs) * w) @ evecs.T
    after = float(np.median(np.linalg.norm(t, axis=1))) or 1.0
    whitened = pre.copy()
    whitened[:, wcols] = t * (before / after)
    frac, _ = solve_fraction(whitened, cd4, cfg["cd4_f"])
    final, _ = plolo.E062.add_post_whiten_residual(whitened, cd4, frac, wcols)
    final[:, fold["panel_positions"]] = 0.0
    out = {}
    for fname, b, gate in frames:
        mat = deliver(final, b, gate, cfg["mass"])
        out[fname] = endpoints(mat[:, fold["scored"]], truth)["pds"]
    return float(np.mean(list(out.values()))), out


def main() -> None:
    rows = {}
    for held in FOLD_ORDER:
        fold = build_fold(held)
        truth = np.asarray(fold["truth"], dtype=np.float64)[:, fold["scored"]]
        wcols = np.asarray(fold["wcols"])
        cd4 = cd4_cached(fold["targets"])
        frames = frames_for(fold)
        op = Path(fold["op_path"])
        with np.load(op) as z:
            evals = np.asarray(z["evals"], float)
            evecs = np.asarray(z["evecs"], float)
        a15, _ = score(fold, CONFIGS["a15"], evals, evecs, wcols, cd4, frames, truth)
        a23, _ = score(fold, CONFIGS["a23"], evals, evecs, wcols, cd4, frames, truth)
        rows[held] = dict(op=str(op.relative_to(ROOT)), op_sha256=sha256_file(op),
                          a15=a15, a23=a23, delta=a23 - a15)
        print(f"{held:<14} op {sha256_file(op)[:8]}  a15 {a15:.6f}  a23 {a23:.6f}  "
              f"d {a23-a15:+.6f}", flush=True)
        del evals, evecs

    d = float(np.mean([r["delta"] for r in rows.values()]))
    scaled = d / SPAN_PDS
    print(f"\n4-fold unweighted mean d_pds = {d:+.6f}   scaled = {scaled:+.6f}")
    print(f"162/169 expected scaled       = +0.004014   |diff| = {abs(scaled-0.004014):.2e}")
    ok = abs(scaled - 0.004014) < 5e-6
    (HERE / "work/gate2-massfix.json").write_text(json.dumps(
        dict(folds=rows, mean_dpds=d, scaled=scaled, expected=0.004014, pass_=bool(ok)),
        indent=1))
    print("GATE 2 PASSED" if ok else "GATE 2 FAILED")
    raise SystemExit(0 if ok else 1)


if __name__ == "__main__":
    main()
