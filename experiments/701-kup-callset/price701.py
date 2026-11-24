#!/usr/bin/env python3
"""701 — price K-UP (CPM 2,800 / 3,000 / 3,400) on 537's real production-emitter H1 runs, beside the K pair (K 2,000 /
1,600) whose board reads exist, at 521's regimes plus f 0.87 / 0.75 / 0.61 (the board's disputed size regime).

Read-only reuse of 537's pricing code (price537.load / prescribed, price528.rows_at / boot, lib521's member functions).
Stats come from experiments/537-arm-pipeline/work/h1_<mode>/ (run537.py writes there; 701's runs add the cpm-n2800 /
n3000 / n3400 dirs).

Known answer (asserted before any K-UP number is printed): "P2 K 2,000 alone" and "P3 K 1,600 alone" at f 1 / 0.5 / 0.35 /
cmp equal 537's recorded work/price537.json values (value_mean and fid_jac_overall) to 1e-12.

Board calibration: for each regime, the least-squares multiplier m that maps H1's Δ(fid + jac) of the K pair onto the board's
reads (K20 −0.00086, K16 −0.00234; 523's BOARD-RESULT files) and its residual; the K-UP prediction at that regime is m × H1.
Writes work/price701.json; prints one line per (comparison, regime).

usage: uv run python experiments/701-kup-callset/price701.py
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

sys.dont_write_bytecode = True
import numpy as np  # noqa: E402

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(ROOT / "experiments/537-arm-pipeline"))
sys.path.insert(0, str(ROOT / "experiments/528-callprop-block"))
sys.path.insert(0, str(ROOT / "experiments/521-strength-nreal"))
import lib521 as L  # noqa: E402
import price528 as P5  # noqa: E402
import price537 as P7  # noqa: E402
from global521 import board_cols  # noqa: E402
from policy521 import SGN, means  # noqa: E402
from table521 import regime_pool  # noqa: E402

REG = ("1.0", "0.87", "0.75", "0.61", "0.5", "0.35", "cmp")
MEMBERS = ("fid", "jac", "reach", "nmae")
BOARD = {"cpm-n1600": -0.00086, "cpm-n1200": -0.00234}     # Δ(fid + jac) overall vs 446bj2z (K20, K16 board reads)
COMPARISONS = [
    ("P2 K 2,000 alone (446bj2zk20's construction)", "cpm-n1600", "base"),
    ("P3 K 1,600 alone (446bj2zk16's construction)", "cpm-n1200", "base"),
    ("701 K 3,200 (CPM 2,800)", "cpm-n2800", "base"),
    ("701 K 3,400 (CPM 3,000) = the arm", "cpm-n3000", "base"),
    ("701 K 3,800 (CPM 3,400)", "cpm-n3400", "base"),
]


def main():
    seeds_all = (0, 1, 2)
    states = {s: L.load_state(s) for s in seeds_all}
    pools = {(s, r): regime_pool(states[s], r)[0] for s in seeds_all for r in REG}
    keys = dict(np.load(ROOT / "experiments/528-callprop-block/work/keys528.npz", allow_pickle=False))
    rec537 = json.loads((ROOT / "experiments/537-arm-pipeline/work/price537.json").read_text())["comparisons"]
    out = {"comparisons": {}, "KA": {}, "calibration": {}}
    for label, arm, b0 in COMPARISONS:
        seeds = [s for s in seeds_all if P7.load(arm, s) is not None and P7.load(b0, s) is not None]
        if not seeds:
            print(f"{label:52s} -- no runs yet", flush=True)
            continue
        rec = {"arm": arm, "base": b0, "seeds": seeds}
        rec["pinned"] = {s: dict(sums_equal=bool(np.array_equal(P7.load(arm, s)["sums"], P7.load(b0, s)["sums"])),
                                 dpds=P7.load(arm, s)["raw"]["pds"] - P7.load(b0, s)["raw"]["pds"],
                                 dmse=P7.load(arm, s)["raw"]["mse"] - P7.load(b0, s)["raw"]["mse"]) for s in seeds}
        for reg in REG:
            per, RA, RB = [], [], []
            for s in seeds:
                Ra = P5.rows_at(states[s], P7.load(arm, s), pools[(s, reg)])
                Rb = P5.rows_at(states[s], P7.load(b0, s), pools[(s, reg)])
                dm = means(Ra) - means(Rb)
                per.append(dict(value=float(SGN @ dm), d=dict(zip(MEMBERS, dm.tolist())),
                                fid_jac=float(SGN[0] * dm[0] + SGN[1] * dm[1])))
                RA.append(Ra)
                RB.append(Rb)
            Ra, Rb = np.nanmean(np.stack(RA), 0), np.nanmean(np.stack(RB), 0)
            dm = np.array([np.mean([p["d"][m] for p in per]) for m in MEMBERS])
            lo, hi, p0 = P5.boot(Ra, Rb)
            r = dict(per_seed=per, value_mean=float(np.mean([p["value"] for p in per])), ci=[lo, hi], p_le0=p0,
                     d_mean=dict(zip(MEMBERS, dm.tolist())), fid_jac_overall=float(SGN[0] * dm[0] + SGN[1] * dm[1]),
                     board=list(board_cols(dm)))
            rec[reg] = r
            print(f"{label[:52]:52s} f {reg:4s} n{len(seeds)} ov {r['value_mean']:+.5f} [{lo:+.5f}, {hi:+.5f}] | fid {dm[0]:+.5f} "
                  f"jac {dm[1]:+.5f} reach {dm[2]:+.5f} nmae {dm[3]:+.5f} | fid+jac {r['fid_jac_overall']:+.5f} "
                  f"(seeds {' '.join(f'{p['fid_jac']:+.5f}' for p in per)})", flush=True)
        s = seeds[0]
        A, B = P7.load(arm, s), P7.load(b0, s)
        own = states[s]["own"]
        C_a, C_b = P7.prescribed(arm, states[s], keys) & ~own, P7.prescribed(b0, states[s], keys) & ~own
        D_a, D_b = (A["q"] < L.ALPHA) & ~own, (B["q"] < L.ALPHA) & ~own
        rec["deliverability"] = dict(prescribed=[float(C_a.sum(1).mean()), float(C_b.sum(1).mean())],
                                     realized=[float(D_a.sum(1).mean()), float(D_b.sum(1).mean())],
                                     missed=[float((C_a & ~D_a).sum(1).mean()), float((C_b & ~D_b).sum(1).mean())],
                                     unsuppressible=[float((D_a & ~C_a).sum(1).mean()), float((D_b & ~C_b).sum(1).mean())])
        dl = rec["deliverability"]
        print(f"   deliverability s{s} (arm | base): prescribed {dl['prescribed'][0]:.0f} | {dl['prescribed'][1]:.0f}; realized "
              f"{dl['realized'][0]:.0f} | {dl['realized'][1]:.0f}; missed {dl['missed'][0]:.1f} | {dl['missed'][1]:.1f}; "
              f"sums equal {all(v['sums_equal'] for v in rec['pinned'].values())}", flush=True)
        out["comparisons"][label] = rec
        # known answer: the K pair at 537's regimes must equal 537's record
        if label in rec537:
            ok = all(abs(rec[reg]["value_mean"] - rec537[label][reg]["value_mean"]) < 1e-12
                     and abs(rec[reg]["fid_jac_overall"] - rec537[label][reg]["fid_jac_overall"]) < 1e-12
                     for reg in ("1.0", "0.5", "0.35", "cmp"))
            out["KA"][label] = ok
            print(f"   KA vs 537's price537.json: {'PASS' if ok else 'FAIL'}", flush=True)
            assert ok, label
    # board calibration of the size regime: m * H1 fid+jac of (K20, K16) ~ board (K20, K16)
    pair = [l for l in out["comparisons"] if out["comparisons"][l]["arm"] in BOARD]
    if len(pair) == 2:
        for reg in REG:
            h = np.array([out["comparisons"][l][reg]["fid_jac_overall"] for l in pair])
            b = np.array([BOARD[out["comparisons"][l]["arm"]] for l in pair])
            m = float(h @ b / (h @ h)) if h @ h > 0 else float("nan")
            resid = float(np.sqrt(np.mean((m * h - b) ** 2)))
            preds = {out["comparisons"][l]["arm"]: m * out["comparisons"][l][reg]["fid_jac_overall"]
                     for l in out["comparisons"] if out["comparisons"][l]["arm"].startswith("cpm-n") and
                     out["comparisons"][l]["arm"] not in BOARD}
            out["calibration"][reg] = dict(m=m, rms_resid=resid, kup_pred=preds)
            print(f"CAL f {reg:4s}: m {m:+.3f} rms resid {resid:.5f} (board K20 -0.00086, K16 -0.00234) -> K-UP board pred "
                  + " ".join(f"{k} {v:+.5f}" for k, v in preds.items()), flush=True)
    (HERE / "work/price701.json").write_text(json.dumps(out, indent=1, default=float))


if __name__ == "__main__":
    main()
