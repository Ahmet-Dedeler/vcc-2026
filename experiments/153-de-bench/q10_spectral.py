"""Q10 — is the whitening trade a FRONTIER, or is production sitting inside it?

Mechanism (AGENTS.md 4.2): hubness comes from a few dominant shared programmes -- growth
arrest, ribosome, stress.  Suppressing THOSE is what buys pds.  But uniform whitening also
AMPLIFIES ~17,000 near-null eigendirections (evals down to 8e-12), and 4.2 says exactly what
that costs: "whitening amplifies low-variance directions, so |predicted effect| becomes a
WORSE per-gene confidence signal afterwards."  153 §8 measured that cost for the first time:
0.0432 raw reach and 0.024 of sign accuracy.

Nobody has ever separated the two.  Only the global scalar gamma has been swept, and it moves
suppression and amplification together.  This screens weight profiles that decouple them.

Both members are priced on the SAME surface and the SAME predicted object, which 088 could not
do.  `pds` here is an lfc-space PROXY (reads 0.7195 where cell-eval2 reads 0.7953); it is used
only to RANK profiles, and its calibration is checked against 088's known gamma=0 collapse.
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
W = np.load("../095-clean-harness/out/whitener_union.npz")
ev, V = W["evals"], W["evecs"]
order = np.argsort(-ev)                      # descending: index 0 = dominant programme
print(f"evals descending: top5 {ev[order[:5]].round(3)}  median {np.median(ev):.3e}  min {ev.min():.3e}")

s = Surface("../095-clean-harness/de_thinned.npz")
gpos = {g: i for i, g in enumerate(s.genes)}
tpos = {t: i for i, t in enumerate(t104)}
gsel = np.array([i for i, g in enumerate(wgenes) if g in gpos])
gdst = np.array([gpos[wgenes[i]] for i in gsel])
prow = np.array([tpos.get(p, -1) for p in s.perts]); pok = prow >= 0
keep = np.zeros_like(s.real_sig); keep[np.ix_(pok, gdst)] = True
s.real_sig &= keep; s.adj &= keep

Qb = base @ V                                 # production (gamma=0.65) in the eigenbasis
U_q = Qb * (ev ** (0.65 / 2.0))               # un-whitened, still in the eigenbasis

def score(w):
    P = ((U_q * w).astype(np.float32) @ V.T)
    F = np.zeros_like(s.pred_lfc); F[np.ix_(pok, gdst)] = P[np.ix_(prow[pok], gsel)]
    s.pred_lfc = F
    r = s.reach(key=[np.abs(F[p]) for p in range(len(s.perts))])
    num = den = 0
    for p in range(len(s.perts)):
        i = np.flatnonzero(s.real_sig[p] & s.adj[p])
        if i.size < 10: continue
        m = i[np.abs(F[p, i]) >= np.percentile(np.abs(F[p, i]), 90)]
        num += int((np.sign(F[p, m]) == np.sign(s.real_lfc[p, m])).sum()); den += m.size
    return s.pds(), r, s.sign_accuracy(), num / max(den, 1)

rows = []
rows.append(("uniform g=0.65  PRODUCTION", ev ** (-0.65 / 2)))
rows.append(("uniform g=0.00  (088: pds collapses)", np.ones_like(ev)))
rows.append(("uniform g=0.45", ev ** (-0.45 / 2)))
for k in (25, 100, 400, 1600):
    w = np.ones_like(ev); w[order[:k]] = ev[order[:k]] ** (-0.65 / 2)
    rows.append((f"suppress top-{k} only, g=0.65", w))
for k in (100, 400):
    w = np.ones_like(ev); w[order[:k]] = ev[order[:k]] ** (-1.0 / 2)
    rows.append((f"suppress top-{k} only, g=1.00", w))
for c in (2.0, 5.0, 20.0):
    rows.append((f"uniform g=0.65, weight cap {c}", np.minimum(ev ** (-0.65 / 2), c)))

print(f"\n{'profile':38} {'pds_proxy':>10} {'reach':>9} {'sign acc':>9} {'head pur':>9}")
b = None
for name, w in rows:
    p_, r_, a_, h_ = score(np.asarray(w, dtype=np.float64))
    if b is None: b = (p_, r_)
    print(f"{name:38} {p_:10.5f} {r_:9.5f} {a_:9.4f} {h_:9.4f}"
          + ("   <- baseline" if b == (p_, r_) else f"   Δpds{p_-b[0]:+.4f} Δreach{r_-b[1]:+.4f}"))
