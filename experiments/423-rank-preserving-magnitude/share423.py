"""Share of the nmae denominator (sum |ref lfc| over the reference DE gate) sitting on genes a variant moved."""
import sys, numpy as np
from pathlib import Path
HERE = Path(__file__).resolve().parent
import run423 as RN
z = np.load(RN.ROOT / "experiments/323-current-family-validation/work/reference.npz")
genes_k = z["genes"][z["keep"]]; targets = z["targets"]; RQ = z["q"]; RL = z["lfc"]
tag = sys.argv[2] if len(sys.argv) > 2 else "s0"
B = np.load(HERE / f"work/h1/stats-base_{tag}.npz")["lfc"]
for v in sys.argv[1].split(","):
    V = np.load(HERE / f"work/h1/stats-{v}_{tag}.npz")["lfc"]
    sh, nm, ns, gain_per = [], [], [], []
    for i in range(len(targets)):
        S = (RQ[i] < .05) & np.isfinite(RL[i]) & (RL[i] != 0) & (genes_k != targets[i])
        if S.sum() < 10: continue
        m = S & (np.nan_to_num(V[i]) != np.nan_to_num(B[i]))
        sh.append(np.abs(RL[i][m]).sum() / np.abs(RL[i][S]).sum()); nm.append(m.sum()); ns.append(S.sum())
    print(f"{v}: moved ref-DE genes/target median {np.median(nm):.0f} (of median |S| {np.median(ns):.0f}); share of sum|r| mean {np.mean(sh):.3f}, median {np.median(sh):.3f}")
