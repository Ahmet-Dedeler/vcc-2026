"""Q11 — the capped-weight profile IS the ridge parameter, and production has always left it 0.

whiten.py: W = (C + lam*I)^(-gamma/2), and `apply(..., lam=0.0)` is the default everywhere.
Capping the weight at c is exactly flooring the eigenvalues, and flooring additively IS the
ridge: w = (ev + lam)^(-gamma/2) has max weight lam^(-gamma/2).

So the q10 result is not an exotic operator -- it is a one-parameter change on a code path
that already exists, has always been available, and has never been swept.  Both members priced
on the same surface and the same object.
"""
import numpy as np, csv
from debench import Surface, overall_delta

canon = [r["gene_name"] for r in csv.DictReader(open("../../data/challenge/gene_names.csv"))]
wcols = np.load("/private/tmp/claude-501/-Users-ahmet-Code-virtual-cell-challenge-2026/"
                "aa534728-8641-48d7-a125-17e75d9c9399/scratchpad/wcols.npy")
wgenes = [canon[i] for i in wcols]
z = np.load("../104-context-heldout-credit-adapter/work/hct116/frozen.npz", allow_pickle=True)
g104 = {str(g): i for i, g in enumerate(z["genes"])}
base = z["baseline"][:, [g104[g] for g in wgenes]].astype(np.float32)
t104 = [str(t) for t in z["targets"]]
W = np.load("../095-clean-harness/out/whitener_union.npz"); ev, V = W["evals"], W["evecs"]

s = Surface("../095-clean-harness/de_thinned.npz")
gpos = {g: i for i, g in enumerate(s.genes)}; tpos = {t: i for i, t in enumerate(t104)}
gsel = np.array([i for i, g in enumerate(wgenes) if g in gpos])
gdst = np.array([gpos[wgenes[i]] for i in gsel])
prow = np.array([tpos.get(p, -1) for p in s.perts]); pok = prow >= 0
keep = np.zeros_like(s.real_sig); keep[np.ix_(pok, gdst)] = True
s.real_sig &= keep; s.adj &= keep
U_q = (base @ V) * (ev ** (0.65 / 2.0))          # un-whitened, in the eigenbasis
med0 = None

def score(gamma, lam):
    global med0
    w = (ev + lam) ** (-gamma / 2.0)
    P = ((U_q * w).astype(np.float32) @ V.T)
    med = float(np.median(np.linalg.norm(P, axis=1)))
    if med0 is None: med0 = med
    P = P * (med0 / med)                          # renorm, so alpha keeps its meaning (whiten.py)
    F = np.zeros_like(s.pred_lfc); F[np.ix_(pok, gdst)] = P[np.ix_(prow[pok], gsel)]
    s.pred_lfc = F
    r = s.reach(key=[np.abs(F[p]) for p in range(len(s.perts))])
    num = den = 0
    for p in range(len(s.perts)):
        i = np.flatnonzero(s.real_sig[p] & s.adj[p])
        if i.size < 10: continue
        m = i[np.abs(F[p, i]) >= np.percentile(np.abs(F[p, i]), 90)]
        num += int((np.sign(F[p, m]) == np.sign(s.real_lfc[p, m])).sum()); den += m.size
    return s.pds(), r, s.sign_accuracy(), num / max(den, 1), s.nmae()

print(f"{'gamma':>6} {'lam':>8} {'pds_proxy':>10} {'reach':>9} {'Δreach':>9} {'sign':>7} {'head':>7} {'nmae':>8}")
b = None
for gamma, lam in [(0.65, 0.0), (0.65, 0.001), (0.65, 0.003), (0.65, 0.005), (0.65, 0.007),
                   (0.65, 0.010), (0.65, 0.020), (0.65, 0.040), (0.65, 0.080),
                   (0.80, 0.007), (0.80, 0.020), (1.00, 0.020), (1.00, 0.040), (0.50, 0.007)]:
    p_, r_, a_, h_, n_ = score(gamma, lam)
    if b is None: b = (p_, r_)
    tag = "  <- PRODUCTION" if (gamma, lam) == (0.65, 0.0) else \
          f"  Δpds{p_-b[0]:+.4f}  ->reach {overall_delta(d_reach=r_-b[1]):+.5f} overall"
    print(f"{gamma:6.2f} {lam:8.3f} {p_:10.5f} {r_:9.5f} {r_-b[1]:+9.5f} {a_:7.4f} {h_:7.4f} {n_:8.5f}{tag}")
