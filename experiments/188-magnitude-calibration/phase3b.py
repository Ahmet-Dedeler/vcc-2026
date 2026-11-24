"""188 Phase 3b — WHY the lfc-space `pds` inverts on h1-test but not on PLOLO.

Phase 3 shows the sign flip is caused by the SPACE, not the decoy pool: on the same h1-test
table, same 85 decoys, same 10,694 genes, `q1.15L` reads -0.02339 in lfc space and +0.00504 in
the scorer's `h` space.  But PLOLO's lfc-space reading does NOT invert (+0.02250 scaled), so
"lfc space inverts it" is not by itself the mechanism.

The candidate mechanism is TRUTH DEPTH.  h1-test's reference log2FC is a ~500-cell arithmetic
mean; PLOLO's truth is a full donor pseudobulk.  Low-expression genes carry most of the lfc
norm (h1-test: mean |lfc_ref| 0.80 in the bottom expression octile vs 0.19 in the top) and are
exactly where a 500-cell mean is noisiest.  `sharpen_lfc` concentrates magnitude there.  So on
h1-test the lfc statistic should punish it and recover as the low-expression genes are dropped.

Test: recompute both spaces on h1-test with an expression floor swept.  No new emission.
"""
from __future__ import annotations
import os, sys, numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, "..", ".."))
E169 = os.path.join(REPO, "experiments/169-sharpening-six-member")
sys.path.insert(0, HERE); sys.path.insert(0, os.path.join(REPO, "experiments/153-de-bench"))
from debench import Surface                     # noqa: E402
from budget188 import PDS_SPAN, pds_of          # noqa: E402
from phase3 import lfc_to_h                     # noqa: E402

E = np.load(os.path.join(E169, "work/h1test-effect.npz"), allow_pickle=True)
pos = {str(g): i for i, g in enumerate(E["genes"])}
S = {t: Surface(os.path.join(E169, f"de169.{t}.npz"))
     for t in ("q1.00_mf", "q1.15_mf", "q1.15L_mf", "q1.25_mf", "q1.25L_mf")}
base = S["q1.00_mf"]
cpm = np.array([E["cpm"][pos[g]] for g in base.genes], dtype=np.float64)
sidx = set(int(i) for i in base.self_idx if i >= 0)
RT_l = np.nan_to_num(base.real_lfc)
RT_h = lfc_to_h(RT_l, cpm)
print(f"{'CPM floor':>10}{'n genes':>9}  " + "".join(
    f"{t.replace('_mf',''):>22}" for t in S if t != "q1.00_mf"))
print(f"{'':>19}  " + "".join(f"{'d lfc':>11}{'d h':>11}" for _ in range(len(S) - 1)))
for floor in (5.0, 10.0, 20.0, 50.0, 100.0, 200.0):
    g = np.array([i for i in range(len(base.genes))
                  if i not in sidx and cpm[i] > floor])
    bl = pds_of(np.nan_to_num(base.pred_lfc)[:, g], RT_l[:, g])
    bh = pds_of(lfc_to_h(np.nan_to_num(base.pred_lfc), cpm)[:, g], RT_h[:, g])
    line = f"{floor:>10.0f}{g.size:>9}  "
    for t in S:
        if t == "q1.00_mf":
            continue
        P = np.nan_to_num(S[t].pred_lfc)
        line += (f"{(pds_of(P[:, g], RT_l[:, g]) - bl)/PDS_SPAN:>+11.5f}"
                 f"{(pds_of(lfc_to_h(P, cpm)[:, g], RT_h[:, g]) - bh)/PDS_SPAN:>+11.5f}")
    print(line, flush=True)
print("\nshare of the reference's total |lfc| and |h| norm carried by genes BELOW each floor")
for floor in (5.0, 10.0, 20.0, 50.0, 100.0):
    m = cpm <= floor
    print(f"  CPM <= {floor:>5.0f}: {100*m.mean():5.1f}% of genes, "
          f"{100*np.abs(RT_l[:, m]).sum()/np.abs(RT_l).sum():5.1f}% of |lfc_ref|, "
          f"{100*np.abs(RT_h[:, m]).sum()/np.abs(RT_h).sum():5.1f}% of |h_ref|")
