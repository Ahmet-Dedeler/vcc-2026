#!/usr/bin/env python3
"""624b — EXACT check of the stage-1 screen: fold522.Oct(panel, fold).score(w) (q 1.75, render, ballast KMAX 16), on the
9 October-like folds, for the percentile cap, 522's strength rider (pct x exp(0.25 z) above the floor) and that rider x
the ridge tilt (theta 0.1, LOPO ridge as in screen624). Reports cos_agg and the relative change vs the percentile cap,
next to stage 1's approximation. Each fold re-runs 522's known answers on load (abort unless |d| < 1e-12).
~1-2 GB per fold (fold522 docstring); run under gate602.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(ROOT / "experiments/522-october-bulk-key"))
import screen624 as S6  # noqa: E402
import fold522 as FD  # noqa: E402


def main() -> None:
    S = {(p, f): S6.load(p, f) for p in S6.PANELS for f in S6.FOLDS}
    out = []
    for hp in S6.PANELS:
        train = [k for k in S if k[0] != hp]
        X = np.concatenate([S[j]["X"][S[j]["car"]] for j in train]); y = np.concatenate([S[j]["y"][S[j]["car"]] for j in train])
        m = S6.ridge_fit(X, y, 100.0)
        for f in S6.FOLDS:
            T = S[(hp, f)]; car = T["car"]
            pred = S6.ridge_pred(m, T["X"])
            lk = T["X"][:, 0]; ok = car & np.isfinite(lk)
            A_ = np.vstack([np.ones(ok.sum()), lk[ok]]).T
            coef, *_ = np.linalg.lstsq(A_, pred[ok], rcond=None)
            r_ = np.zeros_like(pred); r_[ok] = pred[ok] - A_ @ coef
            z = np.zeros_like(pred); z[ok] = np.clip(r_[ok] / (r_[ok].std() + S6.EPS), -2, 2)
            a522 = np.load(S6.W522 / f"arms522-{hp}-{f}.npz")
            above = car & (T["kv"] / T["med"] > 0.5)
            wP = T["w0"]; wS = wP * np.exp(0.25 * np.where(above, a522["z"].astype(float), 0.0)); wST = wS * np.exp(0.1 * z)
            O = FD.Oct(hp, f, check=True)
            assert np.allclose(np.asarray(O.kv, float), T["kv"])
            sP, _ = O.score(wP); sS, _ = O.score(wS); sST, _ = O.score(wST)
            row = dict(panel=hp, fold=f, cos_pct=sP["cos_agg"], cos_S=sS["cos_agg"], cos_ST=sST["cos_agg"],
                       S_vs_pct=100 * (sS["cos_agg"] / sP["cos_agg"] - 1), ST_vs_pct=100 * (sST["cos_agg"] / sP["cos_agg"] - 1),
                       ST_vs_S=100 * (sST["cos_agg"] / sS["cos_agg"] - 1),
                       approx_S_vs_pct=100 * (S6.agg_cos(T, wS) / S6.agg_cos(T, wP) - 1),
                       approx_ST_vs_pct=100 * (S6.agg_cos(T, wST) / S6.agg_cos(T, wP) - 1),
                       pds_pct=sP["pds"], pds_S=sS["pds"], pds_ST=sST["pds"])
            out.append(row)
            print(f"{hp:4s} {f:14s} exact: S {row['S_vs_pct']:+.2f}%  S+tilt {row['ST_vs_pct']:+.2f}% (incr {row['ST_vs_S']:+.2f}%) | "
                  f"approx: S {row['approx_S_vs_pct']:+.2f}%  S+tilt {row['approx_ST_vs_pct']:+.2f}% | pds {sP['pds']:.4f} {sS['pds']:.4f} {sST['pds']:.4f}", flush=True)
    for k in ("S_vs_pct", "ST_vs_pct", "ST_vs_S"):
        v = [r[k] for r in out]
        print(f"{k:10s} mean {np.mean(v):+.2f}%  folds > 0: {sum(x > 0 for x in v)}/9")
    (S6.WORK / "exact624.json").write_text(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()
