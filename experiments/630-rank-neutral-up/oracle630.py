#!/usr/bin/env python3
"""630 — the capacity CEILING: what the UP pass would price at if every selected gene received its full target.

The pass delivers ~25–35 % of what it asks (frac_got; A/B/C 0.28–0.33). Is more capacity worth engineering? Here each
selected gene's delivered lfc is set to its target t = max(l, K · 0.24 / sqrt(ref_lin · lib_med / 1e6)) exactly
(selection as src_rank_up.apply: q >= 0.05, finite, l > 0, not the target's own column), everything else (sums, p, q, jk)
is the zs base's. That is an upper bound for a rank-neutral pass: calls, fid, jac, reach and the bulk cannot move, only
nmae. Writes work/h1/stats-rnup_s<seed>_oracle_k<K>.npz for price630.py.
  uv run python experiments/630-rank-neutral-up/oracle630.py 050 075 100 150
Analysis only; ~4 GB (the H1 reference).
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import rnup630 as U  # noqa: E402


def main() -> int:
    ks = [int(a) for a in sys.argv[1:]] or [100]
    r, c, _ = U.E.load_reference()
    keep = np.asarray(c.keep, dtype=bool)
    kpos = np.cumsum(keep) - 1
    gpos = {str(g): i for i, g in enumerate(r["genes"])}
    lib_med = float(np.median(r["libraries"]))
    mde = 0.24 / np.sqrt(np.maximum(np.asarray(c.ref_lin, float) * lib_med / 1e6, 1e-12))
    for seed in (0, 1, 2):
        f = U.WORK / f"stats-{U.TAG.format(seed=seed)}_dnk250zs.npz"
        if not f.exists():
            continue
        z = np.load(f, allow_pickle=False)
        q, l = z["q"], z["lfc"]
        sel = (q >= 0.05) & np.isfinite(l) & (l > 0)
        for i, t in enumerate(z["targets"]):
            full = gpos.get(str(t), -1)
            if full >= 0 and keep[full]:
                sel[i, kpos[full]] = False
        for K in ks:
            tgt = np.maximum(l, (K / 100.0) * mde[None, :])
            lo = np.where(sel, tgt, l)
            out = {k: z[k] for k in z.files}
            out["lfc"] = lo
            np.savez(U.WORK / f"stats-rnup_s{seed}_oracle_k{K:03d}.npz", **out)
            d = (lo - l)[sel]
            print(f"seed {seed} K {K/100:.2f}: selected {sel.sum(1).mean():.0f}/target, mean Δlfc on selected {d.mean():.4f}, "
                  f"genes raised {(d > 0).sum(1).mean() if d.ndim > 1 else (lo > l).sum(1).mean():.0f}/target", flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
