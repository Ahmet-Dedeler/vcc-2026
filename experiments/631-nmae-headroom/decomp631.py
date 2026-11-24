#!/usr/bin/env python3
"""631 — where is our nmae loss, and how much of it can a MAGNITUDE lever reach? (H1, real emitter's delivered lfc)

The board gap to teams #2-#5 is largest on nmae (scaled 0.168 vs 0.247-0.288, ~0.06 raw) and fid. 630 lifts uncalled
UP genes (-0.009 ... -0.015 raw expected on the board). This decomposes the rest.

Delivered lfc P: 630's zs-base stats (446bj2z's H1 construction, seeds 0-2) and the v3 630 arm on top of it. Truth and
pools: 521's PowerSurface (lib521: rl = the reference lfc; pools(S, f) = the gate at regime f; the own target removed).
Regions (design call set, lib521): head (400 by |key|), CPM block (2,200 by control expression), uncalled (the rest);
each split by the DELIVERED direction (P > 0 UP, P < 0 DOWN, P == 0 none).
Per region x direction, mean over targets of each target's nmae numerator share (sum |P - T| / sum |T| over the pool):
  share   the region's part of the current nmae;
  oracle  Δ if the region's magnitudes became |T| at our delivered signs (the 480 oracle, restricted to the region);
  m*      Δ at the best GLOBAL magnitude factor m (P -> m P on that region only), m chosen on odd targets and scored on
          even targets and vice versa (cross-fitted, so it is a realistic single-knob lever); m grid 0.25 ... 6;
  sign    Δ if the region's signs were perfect at today's magnitudes (what a sign source could buy there).
Regimes f 1 and f 0.35. Raw nmae units; negative = better. Analysis only, < 3 GB.
  uv run python experiments/631-nmae-headroom/decomp631.py
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(ROOT / "experiments/521-strength-nreal"))
import lib521 as L  # noqa: E402

W630 = ROOT / "experiments/630-rank-neutral-up/work/h1"
TAG = "base_s{seed}_bps_prodem_mde0.1cap0.1min0.3coarse"
MS = [0.25, 0.5, 0.75, 1.0, 1.25, 1.5, 2.0, 3.0, 4.0, 6.0]
REG = (1.0, 0.35)


def nmae_rows(P, T, pool):
    """Per-target nmae (sum |P - T| / sum |T| over the pool) and the numerator matrix."""
    num = np.where(pool, np.abs(P - T), 0.0)
    den = np.where(pool, np.abs(T), 0.0).sum(1)
    ok = pool.any(1) & (den > 0)
    return np.where(ok, num.sum(1) / np.where(den > 0, den, 1), np.nan), den, ok


def main() -> int:
    out = {}
    lines = []

    def say(s=""):
        print(s, flush=True)
        lines.append(s)

    for arm in ("zs_base", "v3_630"):
        acc = {}
        for seed in (0, 1, 2):
            S = L.load_state(seed)
            T = np.nan_to_num(S["rl"])
            f0 = W630 / f"stats-{TAG.format(seed=seed)}_dnk250zs.npz"
            f1 = W630 / f"stats-rnup_s{seed}_b_k100_q20r3p10.npz"
            P = np.nan_to_num(np.asarray(np.load(f0 if arm == "zs_base" else f1)["lfc"], float))
            assert P.shape == T.shape, (P.shape, T.shape)
            regions = {"head": S["head"], "block": S["block"], "uncalled": ~S["C"]}
            for f in REG:
                pool = L.pools(S, f)
                base_rows, den, ok = nmae_rows(P, T, pool)
                even = (np.arange(T.shape[0]) % 2) == 0
                for rname, R in regions.items():
                    for dname, D in (("UP", P > 0), ("DOWN", P < 0), ("zero", P == 0)):
                        M = R & D & pool
                        share = np.where(ok, np.where(M, np.abs(P - T), 0).sum(1) / np.where(den > 0, den, 1), np.nan)
                        Po = np.where(M, np.sign(P) * np.abs(T), P)
                        orc = nmae_rows(Po, T, pool)[0] - base_rows
                        Ps = np.where(M, np.sign(T) * np.abs(P), P)
                        sgn = nmae_rows(Ps, T, pool)[0] - base_rows
                        dm = {m: nmae_rows(np.where(M, m * P, P), T, pool)[0] - base_rows for m in MS}
                        # cross-fit: choose m on one half, score on the other
                        m_odd = min(MS, key=lambda m: np.nanmean(dm[m][~even]))
                        m_even = min(MS, key=lambda m: np.nanmean(dm[m][even]))
                        cf = np.where(even, dm[m_odd], dm[m_even])
                        key = (f, rname, dname)
                        a = acc.setdefault(key, {"share": [], "oracle": [], "mstar": [], "sign": [], "m": [], "cells": []})
                        a["share"].append(np.nanmean(share)); a["oracle"].append(np.nanmean(orc))
                        a["mstar"].append(np.nanmean(cf)); a["sign"].append(np.nanmean(sgn))
                        a["m"].append((m_odd, m_even)); a["cells"].append(float(M.sum(1).mean()))
        out[arm] = {}
        say(f"\n==== {arm} (H1, 3 seeds; raw nmae; negative = better)")
        for f in REG:
            say(f"  f {f:g}:  region   dir    cells/target  share    oracle(|T|)  m* (cross-fit)   m chosen        perfect sign")
            for (ff, rname, dname), a in acc.items():
                if ff != f:
                    continue
                rec = {k: float(np.mean(v)) for k, v in a.items() if k != "m"}
                rec["m"] = a["m"]
                out[arm][f"{f:g}|{rname}|{dname}"] = rec
                say(f"          {rname:9s} {dname:5s} {rec['cells']:9.0f}   {rec['share']:+.4f}  {rec['oracle']:+.4f}     "
                    f"{rec['mstar']:+.4f}          {str(sorted(set(sum(map(list, a['m']), []))))[:16]:16s} {rec['sign']:+.4f}")
    (HERE / "work/decomp631.json").write_text(json.dumps(out, indent=1, default=float))
    (HERE / "work/decomp631.out").write_text("\n".join(lines) + "\n")
    return 0


if __name__ == "__main__":
    sys.exit(main())
