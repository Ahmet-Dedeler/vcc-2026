"""Q9 — price the whitening operator against `reach`, the member it was never measured on.

AGENTS.md 4.2: whitening LOWERS the transfer cosine 32% while raising `pds`.  CURRENT_STATE
prices the gamma sweep as "each +0.01 of kappa costs 0.018-0.052 raw pds" and, in the same
paragraph, "reach was unreadable again."  So production's WHITEN_GAMMA = 0.65 was chosen on a
surface where the member with the largest gap to the leader could not be measured.  We are
98th percentile on `pds` and 27th on `reach` among competitive teams, and whitening is the
only lever in the pipeline that changes per-gene SIGNS at scale -- which 153 §6 and the peer's
154 both identify as the binding constraint.

Powers add in the eigenbasis, so applying W(delta) on top of the production gamma=0.65
artifact gives effective gamma 0.65+delta exactly.  The gene axis is reconstructed and
VERIFIED: sorted(canonical index of (donor union & canonical) - panel) == 17,942 columns,
matching `whitener_union.json`'s n_cols.

SCOPE: reach, head purity and sign accuracy only -- all scale-invariant, so the `renorm`/alpha
question cannot confound them.  `pds`, `mse` and `nmae` are NOT priced here and gamma moves
all three.  This measures the missing half of a trade, not the trade.
"""
import numpy as np, csv
from debench import Surface, overall_delta

canon = [r["gene_name"] for r in csv.DictReader(open("../../data/challenge/gene_names.csv"))]
wcols = np.load("/private/tmp/claude-501/-Users-ahmet-Code-virtual-cell-challenge-2026/"
                "aa534728-8641-48d7-a125-17e75d9c9399/scratchpad/wcols.npy")
wgenes = [canon[i] for i in wcols]
print(f"whitener axis {len(wcols)} genes (verified against whitener_union.json n_cols=17942)")

z104 = np.load("../104-context-heldout-credit-adapter/work/hct116/frozen.npz", allow_pickle=True)
g104 = {str(g): i for i, g in enumerate(z104["genes"])}
base = z104["baseline"][:, [g104[g] for g in wgenes]].astype(np.float32)   # (132, 17942)
t104 = [str(t) for t in z104["targets"]]

W = np.load("../095-clean-harness/out/whitener_union.npz")
evals, evecs = W["evals"], W["evecs"]          # float64 (G,), float32 (G,G)
print(f"operator {evecs.shape}, evals in [{evals.min():.3e}, {evals.max():.3e}]")

s = Surface("../095-clean-harness/de_thinned.npz")
gpos = {g: i for i, g in enumerate(s.genes)}
tpos = {t: i for i, t in enumerate(t104)}
gsel = np.array([i for i, g in enumerate(wgenes) if g in gpos])
gdst = np.array([gpos[wgenes[i]] for i in gsel])
prow = np.array([tpos.get(p, -1) for p in s.perts]); pok = prow >= 0
keep = np.zeros_like(s.real_sig); keep[np.ix_(pok, gdst)] = True
s.real_sig &= keep; s.adj &= keep
print(f"scored on {pok.sum()} perturbations x {gsel.size} genes\n")

Q = base @ evecs                                # project once, reuse for every delta
print(f"{'eff. gamma':>11} {'Δ':>7} {'sign acc':>9} {'head purity':>12} {'reach':>9} {'Δreach':>9} {'Δoverall':>10}")
b_r = None
for d in (-0.65, -0.45, -0.30, -0.15, 0.0, 0.15, 0.35):
    P = base if d == 0.0 else ((Q * (evals ** (-d / 2.0))).astype(np.float32) @ evecs.T)
    F = np.zeros_like(s.pred_lfc); F[np.ix_(pok, gdst)] = P[np.ix_(prow[pok], gsel)]
    s.pred_lfc = F
    acc = s.sign_accuracy()
    r = s.reach(key=[np.abs(F[p]) for p in range(len(s.perts))])
    num = den = 0
    for p in range(len(s.perts)):
        i = np.flatnonzero(s.real_sig[p] & s.adj[p])
        if i.size < 10: continue
        m = i[np.abs(F[p, i]) >= np.percentile(np.abs(F[p, i]), 90)]
        num += int((np.sign(F[p, m]) == np.sign(s.real_lfc[p, m])).sum()); den += m.size
    if b_r is None and d == -0.65: pass
    if d == 0.0: b_r = r
    print(f"{0.65+d:11.2f} {d:+7.2f} {acc:9.4f} {num/max(den,1):12.4f} {r:9.5f} "
          f"{'' if b_r is None else f'{r-b_r:+9.5f}'} {'' if b_r is None else f'{overall_delta(d_reach=r-b_r):+10.5f}'}")
