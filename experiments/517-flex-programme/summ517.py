"""517: summarise price517 JSONs (mean over seeds; seed range of overall-equiv; seed-0 CI)."""
import json, sys
from pathlib import Path
import numpy as np
W = Path(__file__).resolve().parent / "work"
for f in sys.argv[1:]:
    d = json.loads((W / f).read_text())
    print(f"== {f}: base {d['base']} (strong {d.get('n_strong')})")
    seeds = sorted(d["seeds"])
    names = list(d["seeds"][seeds[0]])
    for n in names:
        rows = [d["seeds"][s][n] for s in seeds]
        m = lambda k: np.mean([r[k] for r in rows])
        ovs = [r["overall_equiv"] for r in rows]
        acc = rows[0].get("comp_gate_acc", {})
        print(f"{n:24s} nmae {m('nmae'):+.5f} fid {m('fid'):+.5f} reach {m('reach'):+.5f} jac {m('jac'):+.5f} | ov {np.mean(ovs):+.5f} "
              f"({min(ovs):+.5f}..{max(ovs):+.5f}) s0CI [{rows[0]['ci95'][0]:+.4f},{rows[0]['ci95'][1]:+.4f}] | ne ov {m('ne_overall_equiv'):+.5f} "
              f"ne nmae {m('ne_nmae'):+.5f} | flips {m('flips_median') if 'flips_median' in rows[0] else float('nan'):.0f} "
              f"fd {m('fd_median') if 'fd_median' in rows[0] else float('nan'):.3f} | acc {acc.get('acc', float('nan')):.4f} vs {acc.get('acc_delivered', float('nan')):.4f} n {acc.get('n', 0)}")
