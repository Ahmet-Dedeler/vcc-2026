#!/usr/bin/env python3
"""622n re-centring check (10-04, after 657's calib657 found a global UP offset in the Broad-16 KO data): is identity's
low-I -> UP tail on uncalled gate genes still above each line's OWN uncalled base rate, line by line?

622n reported "bottom 10 % of I -> UP 0.781, above 0.5 on 12/12 lines". The KO gates lean UP (uncalled base rate of UP
0.733), so "> 0.5" is not the right null. Here the null is each line's own uncalled UP rate (1 - base_down), i.e. the
global offset removed per line. Read from tails622.json (no data reload).
usage: uv run python experiments/622-identity-track/recentre622.py"""
import json
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent


def main():
    d = json.loads((HERE / "work/tails622.json").read_text())
    lines = [k for k in d if not k.startswith("_")]
    rows, out = [], []
    say = lambda x="": (print(x), out.append(x))  # noqa: E731
    say("line | uncalled base UP | bottom-10% I UP (n) | increment over base | bottom-5% increment | top-10% DOWN - base DOWN")
    for k in lines:
        r = d[k]
        base_up = 1 - r["base_down"]
        inc10, inc5 = r["bottom10_up"] - base_up, r["bottom5_up"] - base_up
        top = r["top10_down"] - r["base_down"]
        rows.append(dict(line=k, base_up=base_up, bottom10_up=r["bottom10_up"], n10=r["bottom10_n"], inc10=inc10, inc5=inc5,
                         top10_down_minus_base=top))
        say(f"{k:9s} | {base_up:.3f} | {r['bottom10_up']:.3f} ({r['bottom10_n']}) | {inc10:+.3f} | {inc5:+.3f} | {top:+.3f}")
    inc = np.array([x["inc10"] for x in rows])
    w = np.array([x["n10"] for x in rows], float)
    s = dict(lines_above_base=int((inc > 0).sum()), n_lines=len(rows), mean_increment=float(inc.mean()),
             median_increment=float(np.median(inc)), pair_weighted_increment=float((inc * w).sum() / w.sum()),
             lines_above_base_bottom5=int((np.array([x["inc5"] for x in rows]) > 0).sum()),
             lines_top10_down_above_base=int((np.array([x["top10_down_minus_base"] for x in rows]) > 0).sum()))
    say("")
    say(f"bottom-10% I -> UP above the line's own base: {s['lines_above_base']}/{s['n_lines']} lines; increment mean "
        f"{s['mean_increment']:+.3f}, median {s['median_increment']:+.3f}, pair-weighted {s['pair_weighted_increment']:+.3f} "
        f"(622n's '> 0.5' read: 12/12, 0.781)")
    say(f"bottom-5%: {s['lines_above_base_bottom5']}/{s['n_lines']} above base; top-10% DOWN above base DOWN: "
        f"{s['lines_top10_down_above_base']}/{s['n_lines']}")
    (HERE / "work/recentre622.json").write_text(json.dumps(dict(rows=rows, summary=s), indent=1))
    (HERE / "work/recentre622.out").write_text("\n".join(out) + "\n")


if __name__ == "__main__":
    main()
