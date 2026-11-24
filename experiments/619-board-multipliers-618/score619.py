#!/usr/bin/env python3
"""619: re-derive 611's H1 -> board multipliers on 618's validated h1-trainval surface (coordinator 10-02 ~23:40Z).

For every board DE pair with a 618 run: the raw member deltas at 521's four regimes (f 1 / 0.5 / 0.35 / cmp) with the
SAME functions 611 used through 537's price537 (528's rows_at: fid / jac / reach / nmae per target on 521's truth pools,
realized order; means over targets; 3 seeds averaged). Two changes only, both stated:
  * the state is 618's (lib521.REF / EFF / ST435 -> ref_tv, the tv faithful key, 618's g2mde runs = 435's exact flags),
    as 618/price618.py patches them;
  * no "per-gene sums equal" assert: price630 has it for rank-neutral arms, but the coarse-MDE pair moves the sums by design,
    and the DE members do not read them.
Then board / 618 per member and regime, beside 611's board / h1-test.
cmp keeps 521's H1-calibrated constants (as price630 does); f 1 / 0.5 / 0.35 are the surface's own.
    uv run python score619.py   ->  work/score619.json, work/score619.md
"""
import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
W618 = ROOT / "experiments/618-h1-trainval-surface/work"
sys.dont_write_bytecode = True
sys.path.insert(0, str(ROOT / "experiments/521-strength-nreal"))
sys.path.insert(0, str(ROOT / "experiments/528-callprop-block"))
import lib521 as L  # noqa: E402

L.REF = W618 / "ref_tv/reference.npz"
L.EFF = W618 / "tv/effects.npz"
L.ST435 = W618 / "g2mde/stats-base_s{seed}_bps_prodem_mde0.1cap0.1min0.3coarse.npz"
L.load_state.cache_clear()
import price528 as P5  # noqa: E402
from policy521 import means  # noqa: E402
from table521 import regime_pool  # noqa: E402

REG = ("1.0", "0.5", "0.35", "cmp")
MEM = ("fid", "jac", "reach", "nmae")
ZR = "base_s{s}_bps_zr-5_prodem_mde0.1cap0.1min0.3coarse"
PAIRS = {   # name: (arm stats, base stats, board raw delta per member, 611's h1-test value per regime)
    "zs (446bj2z - 446bj2)": (W618 / "g1" / f"stats-{ZR}_dnk250zs.npz", W618 / "g1" / f"stats-{ZR}.npz",
                              dict(nmae=-0.00218, fid=0.0, jac=0.0, reach=0.0),
                              dict(nmae={"1.0": -0.00183, "0.5": -0.00182, "0.35": -0.00185, "cmp": -0.00318})),
    "coarse MDE (433a - 417e)": (W618 / "g2mde" / "stats-base_s{s}_bps_prodem_mde0.1cap0.1min0.3coarse.npz",
                                 W618 / "g2nomde" / "stats-base_s{s}_bps_prodem.npz",
                                 dict(nmae=-0.01679, reach=+0.00240),
                                 dict(nmae={"1.0": -0.00966, "0.5": -0.01101, "0.35": -0.01059, "cmp": -0.01263},
                                      reach={"1.0": -0.00549, "0.5": -0.00127, "0.35": +0.00207, "cmp": -0.00055})),
    "cp call set (446bj2zcp - 446bj2z), WORN on 618": (W618 / "g3cp" / f"stats-{ZR}_dnk250zs.npz",
                                                      W618 / "g1" / f"stats-{ZR}_dnk250zs.npz",
                                                      dict(fid=-0.00358, jac=+0.00399),
                                                      dict(fid={"1.0": -0.00126, "0.5": -0.00176, "0.35": -0.00183, "cmp": -0.00183},
                                                           jac={"1.0": +0.00910, "0.5": +0.00887, "0.35": +0.00829, "cmp": +0.00918})),
}


def load(path):
    z = np.load(path, allow_pickle=True)
    return dict(p=np.asarray(z["p"], float), q=np.asarray(z["q"], float), lfc=np.asarray(z["lfc"], float))


def main():
    out, md = {}, ["# 619 — board / 618 multipliers beside 611's board / h1-test (raw Δ, mean of 3 seeds)\n",
                   "| pair | member | board Δ | 618 f 1 | f 0.5 | f 0.35 | cmp | 611 h1-test board/H1 at f 0.35 |",
                   "|---|---|---|---|---|---|---|---|"]
    seeds = (0, 1, 2)
    states = {s: L.load_state(s) for s in seeds}
    pools = {(s, r): regime_pool(states[s], r)[0] for s in seeds for r in REG}
    for name, (arm, base, board, h1) in PAIRS.items():
        rec = {}
        for reg in REG:
            dms = []
            for s in seeds:
                A, B = load(str(arm).format(s=s)), load(str(base).format(s=s))
                dms.append(means(P5.rows_at(states[s], A, pools[(s, reg)])) - means(P5.rows_at(states[s], B, pools[(s, reg)])))
            dm = np.mean(np.stack(dms), 0)
            rec[reg] = dict(zip(MEM, map(float, dm)), per_seed=[dict(zip(MEM, map(float, d))) for d in dms])
        out[name] = dict(regimes=rec, board=board)
        for m, bv in board.items():
            if bv == 0.0:
                continue
            cells = []
            for reg in REG:
                v = rec[reg][m]
                ratio = bv / v if abs(v) > 1e-9 else float("nan")
                cells.append(f"{v:+.5f} ({ratio:.2f}×)")
            h = h1.get(m, {}).get("0.35")
            h1r = f"{bv / h:.2f}×" if h else "—"
            md.append(f"| {name} | {m} | {bv:+.5f} | " + " | ".join(cells) + f" | {h1r} |")
            out[name].setdefault("board_over_618", {})[m] = {reg: (bv / rec[reg][m] if abs(rec[reg][m]) > 1e-9 else None) for reg in REG}
            out[name].setdefault("board_over_h1test", {})[m] = {reg: (bv / hv if hv else None) for reg, hv in h1.get(m, {}).items()}
    (HERE / "work/score619.json").write_text(json.dumps(out, indent=1))
    (HERE / "work/score619.md").write_text("\n".join(md) + "\n")
    print("\n".join(md))


if __name__ == "__main__":
    main()
