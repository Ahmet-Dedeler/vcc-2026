#!/usr/bin/env python3
"""624c — 522's adoption bar for the ridge tilt: the ARC-LIKE re-draws (essential-like targets removed, re-drawn to A/B/C's
strength mix, 20 draws per fold; 522's subpanel522 logic and seeds), EXACT fold522 scores with the mask, on the 9
October-like folds. Arms: the percentile cap (re-derived on the sub-panel), 522's rider pct x exp(0.25 z_sub), and that rider
x exp(0.1 z_res) with z_res = the LOPO ridge prediction residualised on log key and standardised over the KEPT carried
targets. Reads: per fold, mean relative cos_agg vs the percentile cap and the fraction of draws > 0; 522's bar asks
>= cap on 8/9 folds. ~1-2 GB; under gate602.
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
import subpanel522 as SP  # noqa: E402


def main() -> None:
    st = np.load(SP.WORK / "strength522.npz")
    S = {(p, f): S6.load(p, f) for p in S6.PANELS for f in S6.FOLDS}
    out = []
    for hp in S6.PANELS:
        train = [k for k in S if k[0] != hp]
        X = np.concatenate([S[j]["X"][S[j]["car"]] for j in train]); y = np.concatenate([S[j]["y"][S[j]["car"]] for j in train])
        m = S6.ridge_fit(X, y, 100.0)
        for f in S6.FOLDS:
            T = S[(hp, f)]
            O = FD.load(hp, f, check=True)
            pred = S6.ridge_pred(m, T["X"])
            val = SP.src_strength(st, hp, O.targets, O.sources)
            pt = [str(t) for t in st[f"{hp}__targets"]]
            ii = np.asarray([pt.index(t) for t in O.targets])
            essl = (np.asarray(st[f"{hp}__dep_mean"], float)[ii] < -0.5) | np.asarray(st[f"{hp}__k562_essential"])[ii]
            ref = SP.src_strength(st, "abc", [str(t) for t in st["abc__targets"]], O.sources)
            ms, info = SP.draws(val, ref, SP.NDRAW, 5220, True, ~essl)
            rS, rST = [], []
            for msk in ms:
                z = SP.zsub(st, O, msk)
                kv, car = O.kv, O.carried
                hpi = FD.pct_hi(kv, car, FD.P_PROD, msk)
                wp = FD.wclip(kv, car, 0.5, hpi, msk)
                wS = wp * np.exp(0.25 * z)
                lk = T["X"][:, 0]
                ok = car & msk & np.isfinite(lk)
                A_ = np.vstack([np.ones(ok.sum()), lk[ok]]).T
                coef, *_ = np.linalg.lstsq(A_, pred[ok], rcond=None)
                zr = np.zeros_like(pred); r_ = pred[ok] - A_ @ coef
                zr[ok] = np.clip(r_ / (r_.std() + S6.EPS), -2, 2)
                wST = wS * np.exp(0.1 * zr)
                sP, _ = O.score(wp, mask=msk); sS, _ = O.score(wS, mask=msk); sT, _ = O.score(wST, mask=msk)
                rS.append(sS["cos_agg"] / sP["cos_agg"] - 1); rST.append(sT["cos_agg"] / sP["cos_agg"] - 1)
            rS, rST = np.array(rS), np.array(rST)
            row = dict(panel=hp, fold=f, n_sub=info["n_sub"], S_mean=float(100 * rS.mean()), ST_mean=float(100 * rST.mean()),
                       incr_mean=float(100 * (rST - rS).mean()), S_frac_pos=float((rS > 0).mean()),
                       ST_frac_pos=float((rST > 0).mean()), incr_frac_pos=float((rST > rS).mean()))
            out.append(row)
            print(f"{hp:4s} {f:14s} n_sub {info['n_sub']:3d} | S {row['S_mean']:+.2f}% ({row['S_frac_pos']:.2f})  "
                  f"S+tilt {row['ST_mean']:+.2f}% ({row['ST_frac_pos']:.2f})  incr {row['incr_mean']:+.2f}% ({row['incr_frac_pos']:.2f})", flush=True)
    for k in ("S_mean", "ST_mean", "incr_mean"):
        v = [r[k] for r in out]
        print(f"{k:10s} mean {np.mean(v):+.2f}%  folds >= 0: {sum(x >= 0 for x in v)}/9")
    (S6.WORK / "arc624.json").write_text(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()
