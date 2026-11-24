"""Separate the SURFACE effect from the MODEL effect, which were perfectly confounded.

orion-thinned's 132 targets are a 100% subset of the challenge's 300, on a bit-identical 18,533
gene axis. So the SHIPPED a8/Wildflower artifact -- the same .vcc the board scored -- can be
scored on orion too. That gives three of the four cells of a 2x2:

                        orion-thinned surface        board surface
    a8 (production)          [measured here]         1.3627, cos 0.0740
    031-somatic              1.9405, cos 0.0325      (not submitted)

    a8-on-orion  vs  a8-on-board     = pure SURFACE effect (same artifact)
    031-on-orion vs  a8-on-orion     = pure MODEL effect   (same surface)

a8 predicts per context, so all three of its context blocks are scored against the same orion
reference and the spread is reported: applying context A's profile to HCT116 is a cross-context
use, but the ARTIFACT is identical, which is the criterion the comparison needs.
"""
import sys
import numpy as np
ROOT = "~/Code/virtual-cell-challenge-2026"
sys.path.insert(0, f"{ROOT}/experiments/096-mse-floor")
from decompose import terms
from msekit import CTRL, group_stats, score

rp, rs, rc, rj, rg = group_stats(f"{ROOT}/data/harness/orion-hct116-thinned.h5ad")
z = np.load(f"{ROOT}/experiments/096-mse-floor/work/a8-stats.npz", allow_pickle=True)
assert np.array_equal(rg, z["genes"].astype(str)), "gene axis differs"
perts = z["perts"].astype(str)
real_t = {str(p) for p in rp if str(p) != CTRL}

print(f"{'arm':<34}{'mse':>9}{'A':>10}{'X':>11}{'sumD':>10}{'claim':>9}{'rho':>7}"
      f"{'deconf cos':>12}")

# the 031-somatic model on orion, restricted to the same 132 rows (it already covers exactly them)
pp, ps, pc, pj, pg = group_stats(
    f"{ROOT}/experiments/095-clean-harness/out/runs/g0.65_q1.0_a1.08_spread0.25/pred.h5ad")
m, _, _ = score(pp, ps, pc, pj, rp, rs, rc, rj, rg)
t = terms(pp, ps, pc, pj, rp, rs, rc, rj, rg)
print(f"{'031-somatic  ON ORION':<34}{m:>9.4f}{t['A']:>10.5f}{t['X']:>11.6f}{t['sumD']:>10.5f}"
      f"{t['claim']:>9.5f}{t['rho']:>7.3f}"
      f"{t['X']/np.sqrt(t['A']*t['sumD']):>12.5f}")

for ctx in "ABC":
    sel = np.flatnonzero(np.char.startswith(perts, f"{ctx}|"))
    tgt = np.asarray([perts[i].split("|", 1)[1] for i in sel])
    ok = np.array([t_ in real_t for t_ in tgt])
    sel, tgt = sel[ok], tgt[ok]
    m, _, _ = score(tgt, z["sums"][sel], z["counts"][sel], z["jk"][sel], rp, rs, rc, rj, rg)
    t = terms(tgt, z["sums"][sel], z["counts"][sel], z["jk"][sel], rp, rs, rc, rj, rg)
    print(f"{f'a8 (SHIPPED) ctx {ctx} ON ORION  n={sel.size}':<34}{m:>9.4f}{t['A']:>10.5f}"
          f"{t['X']:>11.6f}{t['sumD']:>10.5f}{t['claim']:>9.5f}{t['rho']:>7.3f}"
          f"{t['X']/np.sqrt(t['A']*t['sumD']):>12.5f}")

bs = np.load(f"{ROOT}/experiments/096-mse-floor/work/board_solve.npz")
A, X, sD = float(bs["a8_A"]), float(bs["X"]), float(bs["sumD"])
print(f"{'a8 (SHIPPED) ON THE BOARD':<34}{1.3627:>9.4f}{A:>10.5f}{X:>11.6f}{sD:>10.5f}"
      f"{float(bs['a8_claim']):>9.5f}{1.0:>7.3f}{X/np.sqrt(A*sD):>12.5f}")
