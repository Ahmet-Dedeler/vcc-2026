"""208 — cheap consistency check before the heavy emission run.

Does `recon196.load_h1()` (which is what the Hamilton rider needs, because it reconstructs
the pre-rounding real-valued `x`) agree, target for target and count for count, with
`harness185.load_totals()` (which is what `ka202_emit.py` drove `rc_callset` from)?  If it
does not, the composed arm is not sitting on top of the incumbent and the whole experiment
is measuring the wrong pair.
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

ROOT = Path("~/Code/virtual-cell-challenge-2026")
sys.path.insert(0, str(ROOT / "experiments/196-emitted-depth"))
sys.path.insert(0, str(ROOT / "experiments/185-rank-controlled-emission"))
sys.path.insert(0, str(ROOT / "src"))

import harness185 as H  # noqa: E402
import recon196  # noqa: E402

D = recon196.load_h1()
TOT = H.load_totals()
tmask = TOT["targets"] != "non-targeting"
t185 = list(TOT["targets"][tmask])
tot185 = TOT["totals"][tmask]
t196 = list(D["targets"])

print(f"196 targets {len(t196)}   185 targets {len(t185)}   same order: "
      f"{t196 == t185}")
print(f"genes identical: {np.array_equal(D['genes'], TOT['genes'])}")
print(f"totals identical: {np.array_equal(D['sums'].astype(np.int64), tot185.astype(np.int64))}")
print(f"sums dtype {D['sums'].dtype}  shape {D['sums'].shape}")

C = H.load_controls()
keep, cpm = C["keep"], C["cpm"]
print(f"\nscorer gate |keep| = {int(keep.sum())}   builder gate |D.gate| = "
      f"{int(D['gate'].sum())}   overlap {int((keep & D['gate']).sum())}")

R = recon196.Recon(D["sums"][0], D["r"], D["gate"], n=400, seed=1000)
print(f"\nRecon[0]: L = {R.L:.3f}  emit_L = {R.emit_L(1.0)}  nL = {R.nL:,.1f}  "
      f"off-gate floor==T {R.off_ok:.4f}  m=1 identity "
      f"{np.array_equal(R.totals(1.0), D['sums'][0].astype(np.int64))}")

de = np.load(ROOT / "experiments/202-callset-rider/work/de202.rc_callset.npz",
             allow_pickle=True)
print(f"\nde202.rc_callset.npz: perts {de['perts'].shape}  genes {de['genes'].shape}  "
      f"same target order as 196: {list(map(str, de['perts'])) == t196}")
print(f"  median n_pred {np.median((de['pred_padj'] < 0.05).sum(1)):.0f}")
print(f"  genes == full[keep]: {np.array_equal(de['genes'], C['genes'][keep])}")
