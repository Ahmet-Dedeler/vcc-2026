"""517: summarise rank517.json (means over 3 seeds; seed-0 target-bootstrap CI)."""
import json
from pathlib import Path
import numpy as np
W = Path(__file__).resolve().parent / "work"
d = json.loads((W / "rank517.json").read_text())
print(f"CLEAN4 tau .05 hook coverage: {d['n_c4_hook_genes']} gated genes")
print(f"{'source':7s} {'cov':>5} {'down':>5} | {'fidlike':>7} {'g-blk':>6} {'g-cmp':>6} | " + " | ".join(f"{k:>40s}" for k in ("block_all", "block_g33", "hook_matched", "hook_all")))
for name, rec in d["sources"].items():
    a = rec["acc"]
    cells = []
    for k in ("block_all", "block_g33", "hook_matched", "hook_all"):
        rows = rec[k]
        m = lambda f: np.mean([x[f] for x in rows])
        cells.append(f"ov {m('overall_equiv'):+.4f} [{rows[0]['ci95'][0]:+.4f},{rows[0]['ci95'][1]:+.4f}] n{m('nmae'):+.4f} f{m('fid'):+.4f} r{m('reach'):+.4f}")
    print(f"{name:7s} {rec['coverage']:5.3f} {rec['frac_down']:5.3f} | {a['fidlike']['acc']:7.4f} {a['gate_block']['acc']:6.4f} {a['gate_comp']['acc']:6.4f} | " + " | ".join(cells))
