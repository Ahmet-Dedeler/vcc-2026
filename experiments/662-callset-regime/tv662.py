#!/usr/bin/env python3
"""662 — 619's three board pairs on 618's h1-trainval surface (zs, coarse MDE, cp [WORN: 528's beta was chosen on trainval]),
re-scored on 662's fine PowerSurface grid with 619's exact code path (score619.py: lib521 state patched to ref_tv / the tv key /
618's g2mde runs; 528's rows_at; 3 seeds). The second panel for 662's regime question.
Known answer: at f 1 / 0.5 / 0.35 / cmp the member deltas equal 619's work/score619.json to 1e-12 (asserted).
Writes work/tv662.json.   PYTHONDONTWRITEBYTECODE=1 uv run python experiments/662-callset-regime/tv662.py"""
import json
import sys
from pathlib import Path

sys.dont_write_bytecode = True
import numpy as np  # noqa: E402

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(ROOT / "experiments/619-board-multipliers-618"))
import score619 as S619  # noqa: E402  (patches lib521 to 618's state on import)
from regime662 import REGS  # noqa: E402

L, P5, means, regime_pool = S619.L, S619.P5, S619.means, S619.regime_pool
assert str(L.REF).endswith("ref_tv/reference.npz")


def main():
    ref = json.load(open(ROOT / "experiments/619-board-multipliers-618/work/score619.json"))
    seeds = (0, 1, 2)
    states = {s: L.load_state(s) for s in seeds}
    out, kerr = {}, 0.0
    for name, (arm, base, board, _h1) in S619.PAIRS.items():
        runs = {s: (S619.load(str(arm).format(s=s)), S619.load(str(base).format(s=s))) for s in seeds}
        rec = {}
        for reg in REGS:
            dms = []
            for s in seeds:
                pool = regime_pool(states[s], reg)[0]
                A, B = runs[s]
                dms.append(means(P5.rows_at(states[s], A, pool)) - means(P5.rows_at(states[s], B, pool)))
            dm = np.mean(np.stack(dms), 0)
            rec[reg] = dict(zip(S619.MEM, map(float, dm)), sd=dict(zip(S619.MEM, map(float, np.std(np.stack(dms), 0, ddof=1)))))
            if reg in ref[name]["regimes"]:
                kerr = max(kerr, max(abs(ref[name]["regimes"][reg][m] - rec[reg][m]) for m in S619.MEM))
        out[name] = dict(regimes=rec, board=board)
        print(name)
        for reg in REGS:
            print(f"  f {reg:4s} " + " ".join(f"{m} {rec[reg][m]:+.5f}" for m in S619.MEM)
                  + " | board/618 " + " ".join(f"{m} {board[m] / rec[reg][m]:.2f}" for m in board if board[m] and abs(rec[reg][m]) > 1e-9),
                  flush=True)
    assert kerr < 1e-12, kerr
    out["KA_max_err"] = kerr
    (HERE / "work/tv662.json").write_text(json.dumps(out, indent=1))
    print(f"KA vs score619.json: max |err| {kerr:.1e}")


if __name__ == "__main__":
    main()
