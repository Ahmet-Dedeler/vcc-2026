"""Reproduce 153 §8's production row before building anything."""
import numpy as np, sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from bench import load_frame, score, HDR, ROW

s, inc, gdst, pok = load_frame()
print(f"frame: {pok.sum()}/{s.perts.size} perturbations x {gdst.size} genes "
      f"(surface carries {s.genes.size} gated genes)")
nr = np.array([int(s.real_sig[p].sum()) for p in range(len(s.perts))])
print(f"reference DE pool: median {np.median(nr):.0f}, mean {nr.mean():.1f}, "
      f"total {nr.sum()}")
print()
print(HDR)
r = score(s, inc)
print(ROW.format(name="incumbent (104 baseline, g=0.65)", **r))
print()
print("153 §8 target:                        0.5366       0.6646   0.09040   0.7953")
