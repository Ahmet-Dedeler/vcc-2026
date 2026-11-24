#!/usr/bin/env python3
"""662 — 496a's own H1 pair (496 τ 0.05 hook on the PLAIN coarse base vs 435's production run, seeds 0-2: exactly 496a's board
construction, 496a vs 446bj2) on 662's grid. It is the denominator of T = 496a's board Δnmae / H1 Δnmae at the regime an
instrument reads, which scales CPH's nmae (537's prereg used f 1's −0.00146).
Known answer: at f 1 each run's members equal its own metrics json (raw) to 1e-12, and the f 1 Δnmae rounds to 496's −0.00146.
Writes work/t496_662.json.   PYTHONDONTWRITEBYTECODE=1 uv run python experiments/662-callset-regime/t496_662.py"""
import json
import sys
from pathlib import Path

sys.dont_write_bytecode = True
import numpy as np  # noqa: E402

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(ROOT / "experiments/528-callprop-block"))
sys.path.insert(0, str(ROOT / "experiments/521-strength-nreal"))
import lib521 as L  # noqa: E402
import price528 as P5  # noqa: E402
from policy521 import means  # noqa: E402
from table521 import regime_pool  # noqa: E402
from regime662 import REGS, MEMBERS  # noqa: E402

H496 = ROOT / "experiments/496-uncalled-sign-hook/work/h1/{k}-base_s{s}_bps_hookem_mde0.1cap0.1min0.3coarseas0.05.{e}"
H435 = ROOT / "experiments/435-zs-on-coarse/work/h1/{k}-base_s{s}_bps_prodem_mde0.1cap0.1min0.3coarse.{e}"


def load(tpl, s):
    z = np.load(str(tpl).format(k="stats", s=s, e="npz"), allow_pickle=True)
    raw = json.load(open(str(tpl).format(k="metrics", s=s, e="json")))["raw"]
    return dict(p=np.asarray(z["p"], float), q=np.asarray(z["q"], float), lfc=np.asarray(z["lfc"], float), sums=np.asarray(z["sums"]), raw=raw)


def main():
    seeds = (0, 1, 2)
    states = {s: L.load_state(s) for s in seeds}
    runs = {s: (load(H496, s), load(H435, s)) for s in seeds}
    assert all(np.array_equal(a["sums"], b["sums"]) for a, b in runs.values()), "496a is DE-only: bulk sums must be equal"
    out, kerr = {}, 0.0
    for reg in REGS:
        dms = []
        for s in seeds:
            pool = regime_pool(states[s], reg)[0]
            A, B = runs[s]
            ma, mb = means(P5.rows_at(states[s], A, pool)), means(P5.rows_at(states[s], B, pool))
            if reg == "1.0":
                kerr = max(kerr, max(abs(ma[i] - A["raw"][m]) for i, m in enumerate(MEMBERS)),
                           max(abs(mb[i] - B["raw"][m]) for i, m in enumerate(MEMBERS)))
            dms.append(ma - mb)
        d = np.stack(dms)
        out[reg] = dict(d_mean=dict(zip(MEMBERS, map(float, d.mean(0)))), sd=dict(zip(MEMBERS, map(float, d.std(0, ddof=1)))))
    assert kerr < 1e-12, kerr
    assert round(out["1.0"]["d_mean"]["nmae"], 5) == -0.00146, out["1.0"]["d_mean"]["nmae"]
    out["KA_max_err"] = kerr
    (HERE / "work/t496_662.json").write_text(json.dumps(out, indent=1))
    for reg in REGS:
        d = out[reg]["d_mean"]
        print(f"f {reg:4s} " + " ".join(f"{m} {d[m]:+.5f}" for m in MEMBERS))
    print(f"KA: members == metrics json at f 1 (max |err| {kerr:.1e}); f 1 nmae {out['1.0']['d_mean']['nmae']:+.6f} (496: -0.00146)")


if __name__ == "__main__":
    main()
