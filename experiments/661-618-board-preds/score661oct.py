#!/usr/bin/env python3
"""661 (4): 625's single-octile prior-only flips priced on 618 (PRED-661-625.md registers the reading).

For octile k = 0..7: Δ = (tv_oct<k>) − (tv_octnone) per seed, members fid / jac / reach / nmae at 521's regimes
(f 1 / 0.5 / 0.35 / cmp) with 528's rows_at on 618's state (lib521's REF / EFF / ST435 -> ref_tv, the tv key, g2mde),
plus realised calls per target (q < 0.05, own gene excluded).
625's §6 rule (Amendment 3) mapped onto 618, with nmae × 1.0–1.2 for the board (619):
  flip k iff x_k^nmae < 0  AND  |x_k^nmae| >= max(0.25 × median_k |x_k^nmae|, 2 σ_k)  AND  x_k^fid >= −0.002 raw,
  where σ_k = 618's per-seed sd of Δnmae × sqrt(126 / 900) (board pooling of 3 × 300 targets, as 625 does with
  √(85/900) for H1). The regime read is f 0.35 (619), with f 1 beside it.
Measurability (coordinator): does any octile move nmae by more than 2 per-seed sd on 618?
-> work/score661oct.json"""
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
from policy521 import SGN, means  # noqa: E402
from table521 import regime_pool  # noqa: E402

REG = ("1.0", "0.5", "0.35", "cmp")
MEM = ("fid", "jac", "reach", "nmae")
TAG = "stats-base_s{s}_bps_zr-5_prodem_mde0.1cap0.1min0.3coarse_dnk250zs.npz"
SEEDS = (0, 1, 2)


def load(path):
    z = np.load(path, allow_pickle=True)
    return dict(p=np.asarray(z["p"], float), q=np.asarray(z["q"], float), lfc=np.asarray(z["lfc"], float),
                sums=np.asarray(z["sums"]))


def main():
    states = {s: L.load_state(s) for s in SEEDS}
    pools = {(s, r): regime_pool(states[s], r)[0] for s in SEEDS for r in REG}
    base = {s: load(HERE / "work/tv_octnone" / TAG.format(s=s)) for s in SEEDS}
    own = states[0]["own"]
    calls0 = {s: float(((base[s]["q"] < 0.05) & ~own).sum(1).mean()) for s in SEEDS}
    Rb = {(s, r): P5.rows_at(states[s], base[s], pools[(s, r)]) for s in SEEDS for r in REG}
    out = {}
    for k in range(8):
        d = HERE / f"work/tv_oct{k}"
        seeds = [s for s in SEEDS if (d / TAG.format(s=s)).exists()]
        if not seeds:
            continue
        rec = {"seeds": seeds}
        A = {s: load(d / TAG.format(s=s)) for s in seeds}
        rec["sums_equal_base"] = all(np.array_equal(A[s]["sums"], base[s]["sums"]) for s in seeds)
        rec["d_calls_per_target"] = float(np.mean([((A[s]["q"] < 0.05) & ~own).sum(1).mean() - calls0[s] for s in seeds]))
        for r in REG:
            dms = np.stack([means(P5.rows_at(states[s], A[s], pools[(s, r)])) - means(Rb[(s, r)]) for s in seeds])
            rec[r] = dict(zip(MEM, map(float, dms.mean(0))), sd=dict(zip(MEM, map(float, dms.std(0, ddof=1) if len(seeds) > 1 else [np.nan] * 4))),
                          overall=float(np.mean(dms @ SGN)), per_seed_nmae=[float(x) for x in dms[:, 3]])
        out[str(k)] = rec
    # 625's §6 on 618 at f 0.35 (and f 1 beside it)
    res = {}
    for r in ("0.35", "1.0"):
        xs = {k: out[k][r]["nmae"] for k in out}
        med = float(np.median(np.abs(list(xs.values())))) if xs else float("nan")
        rows = {}
        for k, x in xs.items():
            sd = out[k][r]["sd"]["nmae"]
            sig = sd * np.sqrt(126 / 900) if np.isfinite(sd) else float("nan")
            bar = max(0.25 * med, 2 * sig) if np.isfinite(sig) else 0.25 * med
            passes = (x < 0) and (abs(x) >= bar) and (out[k][r]["fid"] >= -0.002)
            rows[k] = dict(x_nmae=x, sd_seed=sd, sigma_board=sig, bar=bar, measurable=bool(abs(x) > 2 * sd) if np.isfinite(sd) else None,
                           fid=out[k][r]["fid"], passes_625_rule=bool(passes))
        res[r] = dict(median_abs=med, rows=rows, any_measurable=any(v["measurable"] for v in rows.values() if v["measurable"] is not None),
                      passing=[k for k, v in rows.items() if v["passes_625_rule"]])
    out["rule625"] = res
    out["calls_per_target_base"] = calls0
    (HERE / "work/score661oct.json").write_text(json.dumps(out, indent=1))
    for k in sorted(k for k in out if k.isdigit()):
        o = out[k]
        print(f"oct {k}: seeds {o['seeds']} sums==base {o['sums_equal_base']} dcalls {o['d_calls_per_target']:+.1f} | "
              + " | ".join(f"f {r} nmae {o[r]['nmae']:+.5f}±{o[r]['sd']['nmae']:.5f} fid {o[r]['fid']:+.5f} reach {o[r]['reach']:+.5f}" for r in ("1.0", "0.35")))
    for r, v in res.items():
        print(f"625 rule on 618 at f {r}: passing {v['passing']}; any octile measurable (|x| > 2 seed sd): {v['any_measurable']}")


if __name__ == "__main__":
    main()
