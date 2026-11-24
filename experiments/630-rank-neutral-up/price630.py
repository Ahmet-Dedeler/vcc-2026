#!/usr/bin/env python3
"""630 — price the rank-neutral UP pass at 521's regimes (f 1 / 0.5 / 0.35 / cmp), exactly as 537's price537.py does
(its functions imported read-only: 528's rows_at / boot, 521's member functions, regime pools and board columns).

Comparison = (630 arm stats, base = 630's own zs stats = 446bj2z's H1 construction, KA'd to 435), per seed; per regime the
four DE members' Δ, the overall-equivalent value (521's SGN), the target-bootstrap CI, and the board columns
(fid/jac at 0.72x; reach x1 | x0 | x-0.76). pds/mse pinned check: per-gene sums equal.
usage: uv run python experiments/630-rank-neutral-up/price630.py <variant> [<variant> ...]
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(ROOT / "experiments/528-callprop-block"))
sys.path.insert(0, str(ROOT / "experiments/521-strength-nreal"))
sys.dont_write_bytecode = True
import lib521 as L  # noqa: E402
import price528 as P5  # noqa: E402
from global521 import board_cols  # noqa: E402
from policy521 import SGN, means, per_target_contrib, top3_share  # noqa: E402
from table521 import regime_pool  # noqa: E402

W = HERE / "work" / "h1"
TAG = "base_s{seed}_bps_prodem_mde0.1cap0.1min0.3coarse"
REG = ("1.0", "0.5", "0.35", "cmp")
MEMBERS = ("fid", "jac", "reach", "nmae")


def load(path):
    z = np.load(path, allow_pickle=True)
    return dict(p=np.asarray(z["p"], float), q=np.asarray(z["q"], float), lfc=np.asarray(z["lfc"], float),
                sums=np.asarray(z["sums"]))


def main():
    args = sys.argv[1:]
    base_dir, name = W, ""
    if "--base" in args:                     # the stacked construction's own zs stats directory (e.g. 612's h1_idkgd)
        i = args.index("--base"); base_dir = Path(args[i + 1]); del args[i:i + 2]
    if "--name" in args:
        i = args.index("--name"); name = args[i + 1] + "_"; del args[i:i + 2]
    global TAG
    if "--tag" in args:                      # another construction's block / stats tag (e.g. 654's x125 zr dumps)
        i = args.index("--tag"); TAG = args[i + 1]; del args[i:i + 2]
    variants = args
    seeds = [s for s in (0, 1, 2) if (base_dir / f"stats-{TAG.format(seed=s)}_dnk250zs.npz").exists()]
    states = {s: L.load_state(s) for s in seeds}
    pools = {(s, r): regime_pool(states[s], r)[0] for s in seeds for r in REG}
    out = {}
    for v in variants:
        sv = [s for s in seeds if (W / f"stats-rnup_{name}s{s}_{v}.npz").exists()]
        if not sv:
            print(f"{v}: no runs yet")
            continue
        rec = {"seeds": sv}
        for reg in REG:
            per, RA, RB = [], [], []
            for s in sv:
                A = load(W / f"stats-rnup_{name}s{s}_{v}.npz")
                B = load(base_dir / f"stats-{TAG.format(seed=s)}_dnk250zs.npz")
                assert np.array_equal(A["sums"], B["sums"]), "per-gene sums moved (pds/mse not pinned)"
                Ra = P5.rows_at(states[s], A, pools[(s, reg)])
                Rb = P5.rows_at(states[s], B, pools[(s, reg)])
                dm = means(Ra) - means(Rb)
                per.append(dict(value=float(SGN @ dm), d=dict(zip(MEMBERS, dm.tolist())), board=list(board_cols(dm))))
                RA.append(Ra)
                RB.append(Rb)
            Ra, Rb = np.nanmean(np.stack(RA), 0), np.nanmean(np.stack(RB), 0)
            dm = np.array([np.mean([p["d"][m] for p in per]) for m in MEMBERS])
            lo, hi, p0 = P5.boot(Ra, Rb)
            contrib = per_target_contrib(Ra, Rb)
            rec[reg] = dict(per_seed=per, value_mean=float(np.mean([p["value"] for p in per])), ci=[lo, hi], p_le0=p0,
                            d_mean=dict(zip(MEMBERS, dm.tolist())), board=list(board_cols(dm)),
                            top3_share=top3_share(contrib), n_gain=int((contrib > 0).sum()))
            r = rec[reg]
            print(f"{v:10s} f {reg:4s} n{len(sv)} ov {r['value_mean']:+.5f} [{lo:+.5f}, {hi:+.5f}] (seeds "
                  f"{' '.join(f'{p['value']:+.5f}' for p in per)}) | fid {dm[0]:+.6f} jac {dm[1]:+.6f} reach {dm[2]:+.6f} "
                  f"nmae {dm[3]:+.5f} | board [{r['board'][0]:+.5f} | {r['board'][1]:+.5f} | {r['board'][2]:+.5f}] "
                  f"top3 {r['top3_share']:.2f} gain {r['n_gain']}", flush=True)
        out[v] = rec
    # merge into the existing file (one key per variant), so pricing a new variant never drops earlier ones
    jp = W / f"price630{('_' + name.rstrip('_')) if name else ''}.json"
    old = json.loads(jp.read_text()) if jp.exists() else {}
    old.update(out)
    jp.write_text(json.dumps(old, indent=1, default=float))


if __name__ == "__main__":
    sys.exit(main())
