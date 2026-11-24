"""513f: are our CALLED genes' delivered magnitudes too small? Analytic nmae on H1 when the delivered lfc of called genes
(our q < 0.05) is multiplied by s, uncalled genes untouched (signs, calls and order unchanged; nmae only).
Also split by the gene's control-CPM class (block = top 2,200) and by pool band. In-sample on H1: an upper bound.
    uv run python experiments/513-reach-walk-anatomy/scale513.py
"""
import numpy as np

from walk513 import ALPHA, load

d = load(); G = d["genes"]; T = d["targets"]
ctrl = d["ctrl"]; rank = np.empty(len(ctrl), int); rank[np.argsort(-ctrl)] = np.arange(len(ctrl)); block = rank < 2200
gate = np.zeros_like(d["rl"], bool)
for k, t in enumerate(T):
    gate[k] = (d["rq"][k] < ALPHA) & np.isfinite(d["rl"][k]) & (G != t)
scored = gate.sum(1) >= 10
called = d["oq"] < ALPHA
ol = np.nan_to_num(d["ol"], nan=0.0, posinf=0.0, neginf=0.0)
rl = np.where(gate, d["rl"], 0.0)


def nmae(pred):
    num = np.where(gate, np.abs(pred - rl), 0).sum(1); den = np.abs(rl).sum(1)
    return float(np.mean((num / np.where(den > 0, den, np.nan))[scored]))


print(f"base nmae {nmae(ol):.5f} (known answer 0.94270)")
for sel_name, sel in [("called", called), ("called&block", called & block[None, :]), ("called&rest", called & ~block[None, :])]:
    row = []
    for s in (0.5, 1.5, 2.0, 3.0, 5.0):
        row.append(f"s{s}: {nmae(np.where(sel, ol * s, ol)) - nmae(ol):+.5f}")
    print(f"{sel_name:13s} " + "  ".join(row))
# how big are our called-gene magnitudes relative to truth on gate genes, where the sign is right?
m = gate & called & (np.sign(ol) == np.sign(rl)) & (ol != 0)
ratio = np.abs(ol[m]) / np.abs(rl[m])
print(f"called gate genes with the right sign: n {m.sum()}, |ours|/|truth| median {np.median(ratio):.3f} "
      f"(IQR {np.percentile(ratio,25):.3f}-{np.percentile(ratio,75):.3f})")
acc = (np.sign(ol[gate & called & (ol != 0)]) == np.sign(rl[gate & called & (ol != 0)])).mean()
print(f"sign accuracy of called gate genes {acc:.3f}; of uncalled gate genes "
      f"{(np.sign(ol[gate & ~called & (ol != 0)]) == np.sign(rl[gate & ~called & (ol != 0)])).mean():.3f}")
