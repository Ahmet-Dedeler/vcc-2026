"""Which control-expression vector does the pilot's spearman(r_HE, b_ctrl) = +0.5215 use?"""
import os, sys
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from common189 import bvec, centre, gate_mask, load_table, pair_agreement, spearman, submatrix

LINES = ["orion_hct116", "orion_hek293t", "k562_gwps"]
tabs = {t: load_table(t) for t in LINES}
genes = sorted(set.intersection(*[set(t["genes"]) for t in tabs.values()]))
targets = sorted(set.intersection(*[set(t["perts"]) for t in tabs.values()]))
CEN, GAT, BV = {}, {}, {}
for t in LINES:
    A = centre(submatrix(tabs[t], targets, genes).astype(np.float32))
    CEN[t], GAT[t], BV[t] = A, gate_mask(A), bvec(tabs[t], genes)
kHE, nHE = pair_agreement(CEN["orion_hct116"], CEN["orion_hek293t"], GAT["orion_hct116"], GAT["orion_hek293t"])
kHK, nHK = pair_agreement(CEN["orion_hct116"], CEN["k562_gwps"], GAT["orion_hct116"], GAT["k562_gwps"])
kEK, nEK = pair_agreement(CEN["orion_hek293t"], CEN["k562_gwps"], GAT["orion_hek293t"], GAT["k562_gwps"])
sup = (nHE > 50) & (nHK > 50) & (nEK > 50)
rHE = kHE[sup] / nHE[sup]
c = {}
for t in LINES:
    c[f"b/{t}"] = BV[t][sup]
    c[f"cpm/{t}"] = np.expm1(BV[t])[sup]
c["b/mean3"] = np.mean([BV[t] for t in LINES], axis=0)[sup]
c["b/meanHE"] = np.mean([BV[t] for t in ["orion_hct116","orion_hek293t"]], axis=0)[sup]
c["cpm/mean3"] = np.mean([np.expm1(BV[t]) for t in LINES], axis=0)[sup]
c["cpm/meanHE"] = np.mean([np.expm1(BV[t]) for t in ["orion_hct116","orion_hek293t"]], axis=0)[sup]
c["cpm/min3"] = np.min([np.expm1(BV[t]) for t in LINES], axis=0)[sup]
c["b/median3"] = np.median([BV[t] for t in LINES], axis=0)[sup]
for k, v in sorted(c.items()):
    s = spearman(rHE, v)
    print(f"{k:<22}{s:+.4f}{'   <== pilot +0.5215' if abs(s-0.5215)<0.004 else ''}")
