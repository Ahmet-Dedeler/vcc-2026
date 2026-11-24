#!/usr/bin/env python3
"""662 — 611's two DELIVERY board pairs (zs: 435 coarse vs coarse+zs, seeds 0-2; coarse MDE: 426 no-MDE vs 435 coarse, seed 0)
re-scored on 662's fine PowerSurface grid, with desk611b.py's exact scoring (521's member functions, each run in its own
realized order). Question: does the regime that fits the call-set pairs (cp, K20) still retrodict the delivery pairs?
Known answer: at f 1 / 0.5 / 0.35 / cmp the deltas equal 610's work/desk611b.json to 1e-12 (asserted).
Writes work/delivery662.json.   PYTHONDONTWRITEBYTECODE=1 uv run python experiments/662-callset-regime/delivery662.py"""
import json
import sys
from pathlib import Path

sys.dont_write_bytecode = True
import numpy as np  # noqa: E402

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(ROOT / "experiments/521-strength-nreal"))
import lib521 as L  # noqa: E402
import table521 as TB  # noqa: E402
from regime662 import REGS  # noqa: E402

S = L.load_state(0)
own, rl = S["own"], S["rl"]
adj = L.adjud(S)
H435 = str(ROOT / "experiments/435-zs-on-coarse/work/h1/stats-base_s{s}_bps_prodem_mde0.1cap0.1min0.3coarse{z}.npz")
H426 = str(ROOT / "experiments/426-mde-magnitude/work/h1/stats-base_s0_bps_prodem.npz")


def load(p):
    z = np.load(p, allow_pickle=True)
    assert list(z["targets"].astype(str)) == list(S["targets"])
    return dict(q=np.asarray(z["q"], float), p=np.asarray(z["p"], float), lfc=np.asarray(z["lfc"], float), sums=z["sums"])


def score(run, pool):
    call = (run["q"] < 0.05) & ~own
    fr = L.fid_rows(call, run["lfc"], pool, adj, rl)
    jr = L.jac_rows(call, pool)
    nr = L.nmae_rows(run["lfc"], pool, rl)
    St = dict(S)
    St["oq"], St["op"], St["ol"] = run["q"], run["p"], run["lfc"]
    rr = np.array([L.reach_row(L.realized_order_pos(St, i), pool[i], adj[i], run["lfc"][i], rl[i]) for i in range(len(S["targets"]))])
    return dict(fid=np.nanmean(fr), jac=np.mean(jr), reach=np.nanmean(rr), nmae=np.nanmean(nr))


def main():
    ref = json.load(open(ROOT / "experiments/610-board-forensics/work/desk611b.json"))
    pairs = {"zs (446bj2z-446bj2)": [(H435.format(s=s, z="_dnk250zs"), H435.format(s=s, z="")) for s in (0, 1, 2)],
             "coarse MDE (433a-417e)": [(H435.format(s=0, z=""), H426)]}
    out, kerr = {}, 0.0
    for name, prs in pairs.items():
        runs = [(load(a), load(b)) for a, b in prs]
        assert all(np.array_equal(a["sums"], b["sums"]) for a, b in runs)
        out[name] = {}
        for reg in REGS:
            pool = TB.regime_pool(S, reg)[0]
            ds = []
            for a, b in runs:
                sa, sb = score(a, pool), score(b, pool)
                ds.append({m: float(sa[m] - sb[m]) for m in sa})
            out[name][reg] = {m: float(np.mean([x[m] for x in ds])) for m in ds[0]}
            if reg in ref[name]:
                kerr = max(kerr, max(abs(ref[name][reg][m] - out[name][reg][m]) for m in out[name][reg]))
    assert kerr < 1e-12, kerr
    out["KA_max_err"] = kerr
    (HERE / "work/delivery662.json").write_text(json.dumps(out, indent=1))
    for name in pairs:
        print(name)
        for reg in REGS:
            d = out[name][reg]
            print(f"  f {reg:4s} " + " ".join(f"{m} {d[m]:+.5f}" for m in d))
    print(f"KA vs desk611b.json: max |err| {kerr:.1e}")


if __name__ == "__main__":
    main()
