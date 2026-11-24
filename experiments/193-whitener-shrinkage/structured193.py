"""193 — the one estimator family in the mandate that is NOT rotationally invariant.

Every map in `maps193.py` keeps production's eigenVECTORS and changes only the eigenVALUES, so
none of them can test the mandate's fourth bullet: *shrink toward a STRUCTURED target (diagonal,
block, or a factor model) rather than `I`*.

  * a factor-model target is `clip_rank` and IS rotationally invariant -- already screened;
  * a DIAGONAL target is not.  `C_rho = (1-rho) C + rho diag(C)` has different eigenvectors from
    `C`, so it needs a fresh eigendecomposition and a fresh operator.

`diag(C)` is the per-gene variance of the donor effect vectors, i.e. the "genes vary by
different amounts but are uncorrelated" model.  It is the natural structured target here because
the correlation structure it discards is exactly the shared-programme structure that causes the
hubness whitening exists to fix -- so this arm asks, directly, how much of the operator's value
is in the correlations rather than in the per-gene scales.

KNOWN-ANSWER ROW (rule 10): `rho = 0` reconstructs `C` from the frozen operator, re-solves its
eigensystem from scratch, and must reproduce production's PLOLO row.  That is a real test --
`V diag(ev) V^T` followed by a second `eigh` recovers the same subspaces only if both the
reconstruction and the solver are sound, and it is the only way to know that a `rho > 0` number
is a property of the target rather than of the pipeline.

Run on `kolf2`: the only PLOLO fold where the production operator is leakage-free.
Peak RSS is dominated by two `p x p` float64 arrays (2.6 GiB each) plus LAPACK workspace.
"""
from __future__ import annotations

import argparse
import gc
import json
import os
import resource
import sys
from pathlib import Path

os.environ.setdefault("OPENBLAS_NUM_THREADS", "8")
os.environ.setdefault("OMP_NUM_THREADS", "8")

import numpy as np
from scipy.linalg import eigh as sp_eigh

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(ROOT / "experiments/162-retrodicting-validator"))
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(HERE))

import plolo                                                        # noqa: E402
from plolo import build_fold, deliver, endpoints, solve_fraction    # noqa: E402
from run_arms import CONFIGS                                        # noqa: E402
from screen193 import SPAN_PDS, cd4_cached, frames_for, sha256_file  # noqa: E402

A15 = CONFIGS["a15"]
FLOOR_REL = 1e-12          # `Whitener.fit`'s floor, reapplied so the arms stay comparable


def peak_gib() -> float:
    p = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    return p / 2**30 if sys.platform == "darwin" else p / 2**20


def rebuild(op_path: Path, rho: float, block: int = 2048):
    """C_rho = (1-rho) V diag(ev) V^T + rho diag(...), then its eigensystem."""
    with np.load(op_path) as z:
        ev = np.asarray(z["evals"], float)
        V = np.asarray(z["evecs"], float)
    p = len(ev)
    C = np.empty((p, p), dtype=np.float64)
    Ve = V * ev[None, :]
    for lo in range(0, p, block):
        hi = min(p, lo + block)
        C[lo:hi] = Ve[lo:hi] @ V.T
    del Ve, V
    gc.collect()
    C = 0.5 * (C + C.T)                       # symmetrise the reconstruction
    if rho > 0.0:
        d = np.diag(C).copy()
        C *= (1.0 - rho)
        C[np.diag_indices(p)] += rho * d
    print(f"   reconstructed C ({p}x{p}) rho={rho}, peakRSS {peak_gib():.1f}G; eigh...",
          flush=True)
    w, U = sp_eigh(C, overwrite_a=True, driver="evr")
    del C
    gc.collect()
    w = np.maximum(w, w.max() * FLOOR_REL)
    print(f"   eigh done, maxev {w.max():.4f} minev {w.min():.3e}, peakRSS {peak_gib():.1f}G",
          flush=True)
    return w, U


def score(fold, evals, evecs, wcols, cd4, frames, truth):
    pre = np.asarray(fold["pre"], dtype=np.float64) * A15["alpha"]
    sub = pre[:, wcols]
    before = float(np.median(np.linalg.norm(sub, axis=1))) or 1.0
    w = (evals + A15["lam"]) ** (-A15["gamma"] / 2.0)
    t = ((sub @ evecs) * w) @ evecs.T
    after = float(np.median(np.linalg.norm(t, axis=1))) or 1.0
    whitened = pre.copy()
    whitened[:, wcols] = t * (before / after)
    del t, sub
    frac, _ = solve_fraction(whitened, cd4, A15["cd4_f"])
    final, _ = plolo.E062.add_post_whiten_residual(whitened, cd4, frac, wcols)
    final[:, fold["panel_positions"]] = 0.0
    del whitened
    pds = []
    for fname, b, gate in frames:
        mat = deliver(final, b, gate, A15["mass"])
        pds.append(endpoints(mat[:, fold["scored"]], truth)["pds"])
        del mat
    return float(np.mean(pds)), pds


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--rhos", nargs="*", type=float, default=[0.0, 0.05, 0.25, 0.75])
    ap.add_argument("--fold", default="kolf2")
    args = ap.parse_args()

    fold = build_fold(args.fold)
    truth = np.asarray(fold["truth"], dtype=np.float64)[:, fold["scored"]]
    wcols = np.asarray(fold["wcols"])
    cd4 = cd4_cached(fold["targets"])
    frames = frames_for(fold)
    op = Path(fold["op_path"])
    print(f"[{args.fold}] operator {op.relative_to(ROOT)} sha {sha256_file(op)[:8]}", flush=True)

    out_path = HERE / f"work/structured-{args.fold}.json"
    out = json.loads(out_path.read_text()) if out_path.exists() else {}
    for rho in args.rhos:
        key = f"rho={rho:g}"
        if key in out:
            continue
        w, U = rebuild(op, rho)
        m, per = score(fold, w, U, wcols, cd4, frames, truth)
        out[key] = dict(rho=rho, chal_pds=m, per_frame=per, max_ev=float(w.max()),
                        op=str(op.relative_to(ROOT)), op_sha256=sha256_file(op),
                        peak_rss_gib=peak_gib())
        del w, U
        gc.collect()
        print(f"   {key:12s} chal_pds {m:.6f}", flush=True)
        out_path.write_text(json.dumps(out, indent=1))

    base = out.get("rho=0", {}).get("chal_pds")
    print("\nrho     chal_pds     d vs rho=0    scaled")
    for key, r in sorted(out.items(), key=lambda kv: kv[1]["rho"]):
        d = r["chal_pds"] - base if base else float("nan")
        print(f"{r['rho']:<7g} {r['chal_pds']:.6f}   {d:+.6f}   {d/SPAN_PDS:+.6f}")
    if base is not None:
        # gate: rho = 0 must reproduce the frozen operator's own PLOLO row
        import json as _j
        p1 = _j.loads((HERE / "work/phase1.json").read_text())
        ref = [r for r in p1 if r["fold"] == args.fold and r["map"] == "prod(lam=0)"]
        if ref:
            d = abs(base - ref[0]["chal_pds"])
            print(f"\nGATE 4 (rho=0 vs frozen operator): |d_pds| = {d:.3e} "
                  f"({'PASS' if d < 1e-5 else 'FAIL'})")
    print(f"wrote {out_path}")


if __name__ == "__main__":
    main()
