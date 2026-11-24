#!/usr/bin/env python3
"""390e -- feasibility of CLEAN4 consensus signs at |E0| magnitude on the called CPM block (H1, analytic).
Target A_t = ref_lin * 2^(s*|E0|). Required move from delivered A0 = ref_lin*2^E0 is s|E0| - E0 (0 or -2E0).
Channel range with 3 classes at h1 control quartiles: A in [bulk*L1/Lbar, bulk*L3/Lbar] -> log2 half-widths.
Also the compositional shift absorbed by all other genes and the head/other calls it would flip."""
import numpy as np, json
from pathlib import Path
REPO = Path("~/Code/virtual-cell-challenge-2026"); W = REPO / "experiments/323-current-family-validation/work"
OUT = REPO / "experiments/390-distortion-sign-correction/work"
s = np.load(W / "stats-cap15.npz"); r = np.load(W / "reference.npz"); e = np.load(W / "effects.npz")
cp = np.load(REPO / "experiments/371-consensus-block-prior-lolo/work/consensus_prior.npz", allow_pickle=True)
keep = r["keep"].astype(bool); genes = r["genes"].astype(str)[keep]; targets = s["targets"].astype(str)
E0 = np.nan_to_num(s["lfc"].astype(float)); sig = np.isfinite(s["q"]) & (s["q"] < .05)
eg = {g: i for i, g in enumerate(e["genes"].astype(str))}; col = np.array([eg.get(g, -1) for g in genes]); ok = col >= 0
E1 = e["cap15"].astype(float)[:, np.where(ok, col, 0)]; E1[:, ~ok] = 0
ref_lin = r["ref_lin"].astype(float); sums = s["sums"].astype(float); A0 = (sums / sums.sum(1, keepdims=True) * 1e6)[:, keep]
names = list(cp["bank_names"].astype(str)); PB = cp["P_banks"].astype(float); cg = {g: j for j, g in enumerate(cp["genes"].astype(str))}
cmap = np.array([cg.get(g, -1) for g in genes]); cov = cmap >= 0
idx = [names.index(n) for n in ("k562_gwps", "k562_essential", "nadig_jurkat", "nadig_hepg2")]
P = np.full(len(genes), np.nan); P[cov] = PB[idx][:, cmap[cov]].mean(0)
psign = np.where(np.isfinite(P), np.where(P > 0.5, -1.0, 1.0), -1.0)
L = np.array([16078.75, 20190.5, 25330.5]); Lbar = L.mean()
lo, hi = np.log2(L[0] / Lbar), np.log2(L[2] / Lbar)          # achievable log2 move around bulk (optimistic caps)
d = np.log2(A0.sum(0) * 0 + 1)  # placeholder to keep shapes
nT = len(targets); flipfrac, clip, shift, headflip, othflip, movq = [], [], [], [], [], []
for i in range(nT):
    head = np.zeros(len(genes), bool); head[np.argsort(-np.abs(E1[i]), kind="stable")[:400]] = True
    block = np.zeros(len(genes), bool); block[np.argsort(-np.where(head, -np.inf, ref_lin), kind="stable")[:2200]] = True
    cb = block & sig[i]
    tgt = psign * np.abs(E0[i])                                # target LFC vs ref_lin
    move = np.where(cb, tgt - E0[i], 0.0)                      # log2 move of the arithmetic mean from its welded value
    flipfrac.append(float((move[cb] != 0).mean())); movq.append(np.abs(move[cb][move[cb] != 0]))
    clip.append(float(((move < lo) | (move > hi))[cb].mean()))
    At = A0[i] * 2 ** move
    R = (At[cb] - A0[i, cb]).sum(); O = 1e6 - A0[i, cb].sum(); sh = np.log2((O - R) / O); shift.append(sh)
    oth = sig[i] & ~cb; new = E0[i] + sh; fl = (np.sign(new) != np.sign(E0[i])) & oth
    headflip.append(float((fl & head).sum() / max((oth & head).sum(), 1))); othflip.append(float((fl & ~head).sum() / max((oth & ~head).sum(), 1)))
allmov = np.concatenate(movq)
q = lambda v: [round(float(x), 5) for x in np.quantile(v, [.1, .5, .9])]
out = {"channel_log2_move_bounds": [round(lo, 4), round(hi, 4)], "block_calls_flipped_frac_q10_50_90": q(flipfrac),
       "abs_log2_move_of_flipped_q50_90_99": [round(float(x), 4) for x in np.quantile(allmov, [.5, .9, .99])],
       "block_calls_clipped_frac_q10_50_90": q(clip), "other_gene_log2_shift_q10_50_90": q(shift),
       "head_calls_flipped_by_shift_q10_50_90": q(headflip), "other_calls_flipped_by_shift_q10_50_90": q(othflip)}
print(json.dumps(out, indent=1)); (OUT / "feasible390.json").write_text(json.dumps(out, indent=1))
