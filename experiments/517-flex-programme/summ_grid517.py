"""517: summarise grid517.json -- means over the 3 seeds per (r, w, gated, q, tau); pooled seed range."""
import json
from collections import defaultdict
from pathlib import Path
import numpy as np
W = Path(__file__).resolve().parent / "work"
d = json.loads((W / "grid517.json").read_text())
agg = defaultdict(list)
for r in d["rows"]:
    agg[(r["r"], r["w"], r["gated"], r["q"], r["tau"])].append(r)
def show(filt, title):
    print(f"\n== {title}")
    print(f"{'r':>4} {'w':>5} {'gate':>4} {'q':>4} {'tau':>5} | {'fid':>8} {'reach':>8} {'nmae':>8} | {'ov':>8} ({'min':>8}..{'max':>8}) | {'ne_ov':>8} | {'w/o4':>8}")
    for k, rows in sorted(agg.items(), key=lambda kv: (kv[0][0], -kv[0][1], not kv[0][2], kv[0][3], -1 if kv[0][4] is None else kv[0][4])):
        if not filt(*k):
            continue
        m = lambda f: np.mean([x[f] for x in rows])
        ovs = [x["ov"] for x in rows]
        print(f"{k[0]:4.1f} {k[1]:5.2f} {'g33' if k[2] else 'all':>4} {k[3]:4.1f} {str(k[4]):>5} | {m('fid'):+8.5f} {m('reach'):+8.5f} {m('nmae'):+8.5f} | "
              f"{np.mean(ovs):+8.5f} ({min(ovs):+8.5f}..{max(ovs):+8.5f}) | {m('ne_ov'):+8.5f} | {m('without_four'):+8.5f}")
show(lambda r, w, g, q, tau: w == 1.0 and tau is None, "protect_q (block only), Arc programme at attenuation r (w 1)")
show(lambda r, w, g, q, tau: q == 0.8 and tau is None, "dose w (block only, q 0.8)")
show(lambda r, w, g, q, tau: q == 0.8 and tau is not None, "block + hook (q 0.8), gated hook follows the block gate")
