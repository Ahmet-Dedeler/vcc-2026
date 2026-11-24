"""Estimate the cost of rebuilding the whitening operator for the D/E/F panel.

The operator's axis is `union(donor coverage on the official axis) - panel`
(`scripts/build_union_whitener.py:269-273`), so it is PANEL-DEPENDENT and a D/E/F build
needs a new one. The full fit cannot be run here under the 10 GB RSS cap alongside the
sibling experiments, so its two dominant kernels are timed at small n on this machine and
extrapolated by their known cubic scaling. An estimate inside ~2x is what the Oct-22 plan
needs; anything sharper would require running the real thing.

    uv run python experiments/199-def-runbook/whitenercost199.py
"""

from __future__ import annotations

import json
import time
from pathlib import Path

import numpy as np

OUT = Path(__file__).resolve().parent / "work"
N_TARGET = 17942          # measured by `build_union_whitener.py --dry-run` on the D/E/F bundle
BASIS_ROWS = 29394        # ditto


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(0)
    rows = []
    for n in (1500, 2200, 3000):
        A = rng.standard_normal((n, n))
        C = A @ A.T / n
        t0 = time.perf_counter()
        np.linalg.eigh(C)
        t_eigh = time.perf_counter() - t0
        m = max(BASIS_ROWS * n // N_TARGET, n)
        B = rng.standard_normal((m, n))
        t0 = time.perf_counter()
        B.T @ B
        t_gram = time.perf_counter() - t0
        rows.append({"n": n, "gram_rows": m, "eigh_s": round(t_eigh, 3),
                     "gram_s": round(t_gram, 3)})
        print(f"n={n:>5}  eigh {t_eigh:7.3f}s   gram ({m}x{n}) {t_gram:7.3f}s", flush=True)
        del A, C, B

    # both kernels are cubic in n at fixed aspect ratio; fit the constant on the largest point
    last = rows[-1]
    scale = (N_TARGET / last["n"]) ** 3
    est_eigh = last["eigh_s"] * scale
    est_gram = last["gram_s"] * scale
    # peak RSS: covariance + eigenvectors + LAPACK workspace, all n^2 float64
    est_rss = 3 * N_TARGET ** 2 * 8 / 2 ** 30
    out = {
        "measurements": rows,
        "extrapolated_to": {"n": N_TARGET, "basis_rows": BASIS_ROWS},
        "cubic_scale_factor": round(scale, 1),
        "estimated_eigh_seconds": round(est_eigh, 1),
        "estimated_gram_seconds": round(est_gram, 1),
        "estimated_fit_minutes": round((est_eigh + est_gram) / 60, 1),
        "estimated_peak_rss_gib": round(est_rss, 1),
        "basis_scratch_gib": round(BASIS_ROWS * N_TARGET * 8 / 2 ** 30, 2),
        "output_npz_bytes": 1_287_805_498,
        "caveat": ("extrapolated from n<=3000 by n^3 at fixed aspect ratio; the basis "
                   "assembly (streaming three donor memmaps) and the memmap-backed gram "
                   "are I/O bound and are NOT in this estimate"),
    }
    (OUT / "whitenercost.json").write_text(json.dumps(out, indent=2) + "\n")
    print(f"\nextrapolated to n={N_TARGET}: eigh ~{est_eigh / 60:.1f} min, "
          f"gram ~{est_gram / 60:.1f} min, fit ~{out['estimated_fit_minutes']} min, "
          f"peak RSS ~{est_rss:.1f} GiB, scratch {out['basis_scratch_gib']} GiB")


if __name__ == "__main__":
    main()
