"""Which definition of 'mean |effect|' and 'across-target sd' does the pilot use?

The pilot table quotes spearman(r_HE, mean|effect|) = +0.1125 and (across-target sd) = +0.1578.
My first reading -- mean over the three lines of the per-gene statistic of the CENTRED block --
gives +0.0268 / +0.0748.  Everything else in the table reproduces, so the difference is a
convention in these two scalars only.  Enumerate the plausible conventions and report which
lands on the pilot's numbers.  (Rank statistics, so monotone transforms are irrelevant.)
"""
from __future__ import annotations

import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from common189 import (bvec, centre, gate_mask, load_table, pair_agreement, spearman,  # noqa
                       submatrix)

LINES = ["orion_hct116", "orion_hek293t", "k562_gwps"]
tabs = {t: load_table(t) for t in LINES}
genes = sorted(set.intersection(*[set(t["genes"]) for t in tabs.values()]))
targets = sorted(set.intersection(*[set(t["perts"]) for t in tabs.values()]))

RAW, CEN, GAT, BV = {}, {}, {}, {}
for t in LINES:
    A = submatrix(tabs[t], targets, genes).astype(np.float32)
    RAW[t] = A
    CEN[t] = centre(A)
    GAT[t] = gate_mask(CEN[t])
    BV[t] = bvec(tabs[t], genes)

kHE, nHE = pair_agreement(CEN["orion_hct116"], CEN["orion_hek293t"],
                          GAT["orion_hct116"], GAT["orion_hek293t"])
kHK, nHK = pair_agreement(CEN["orion_hct116"], CEN["k562_gwps"],
                          GAT["orion_hct116"], GAT["k562_gwps"])
kEK, nEK = pair_agreement(CEN["orion_hek293t"], CEN["k562_gwps"],
                          GAT["orion_hek293t"], GAT["k562_gwps"])
sup = (nHE > 50) & (nHK > 50) & (nEK > 50)
rHE = kHE[sup] / nHE[sup]

PAIR = ["orion_hct116", "orion_hek293t"]
cands = {}
for tag, M in (("cen", CEN), ("raw", RAW)):
    for who, ls in (("3line", LINES), ("HE", PAIR)):
        cands[f"mag/{tag}/{who}/mean-of-lines"] = np.mean(
            [np.abs(M[t]).mean(axis=0) for t in ls], axis=0)[sup]
        cands[f"mag/{tag}/{who}/pooled"] = np.mean(
            [np.abs(M[t]) for t in ls], axis=0).mean(axis=0)[sup]
        cands[f"sd/{tag}/{who}/mean-of-lines"] = np.mean(
            [M[t].std(axis=0) for t in ls], axis=0)[sup]
        cands[f"sd/{tag}/{who}/rms-of-lines"] = np.sqrt(np.mean(
            [M[t].std(axis=0) ** 2 for t in ls], axis=0))[sup]
for t in LINES:
    cands[f"mag/cen/{t}"] = np.abs(CEN[t]).mean(axis=0)[sup]
    cands[f"sd/cen/{t}"] = CEN[t].std(axis=0)[sup]
    cands[f"mag/raw/{t}"] = np.abs(RAW[t]).mean(axis=0)[sup]
    cands[f"sd/raw/{t}"] = RAW[t].std(axis=0)[sup]
# gated-only variants: statistic over the entries that actually enter r_g
for t in LINES:
    a = np.abs(CEN[t])
    m = GAT[t]
    cands[f"mag/cen-gated/{t}"] = (np.where(m, a, 0).sum(axis=0) /
                                   np.maximum(m.sum(axis=0), 1))[sup]
# the uncentred column mean itself (the shared programme that centring removes)
cands["mag/shared-mean/3line"] = np.mean(
    [np.abs(RAW[t].mean(axis=0)) for t in LINES], axis=0)[sup]

print(f"{'candidate':<40}{'spearman(r_HE, .)':>20}")
for name, v in sorted(cands.items()):
    s = spearman(rHE, v)
    flag = ""
    if abs(s - 0.1125) < 0.004:
        flag = "  <== pilot's mean|effect| +0.1125"
    if abs(s - 0.1578) < 0.004:
        flag = "  <== pilot's across-target sd +0.1578"
    print(f"{name:<40}{s:>+20.4f}{flag}")
