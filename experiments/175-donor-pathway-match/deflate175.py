"""175b -- per-DONOR rank-k deflation before the consensus: does removing each donor's own
dominant direction (the toxicity x blob axis) improve sign accuracy on H1?  Fit the PCs on ALL of
the donor's perturbations (not the 85 h1 targets), so the 'blob' is the donor's own, not h1's."""
import numpy as np, sys
sys.path.insert(0, "experiments/153-de-bench"); sys.path.insert(0, "experiments/172-sign-confidence-key")
from debench import Surface
from donors172 import DonorHandle
S = Surface("experiments/164-callset-artefact/de164.a15.npz")
pool = S.real_sig & S.adj; sT = np.sign(S.real_lfc)
genes = np.asarray([str(g) for g in S.genes]); targets = [str(p) for p in S.perts]
def acc(key, m): return np.mean(key[m] == sT[m])
D = {}; PCS = {}
rng = np.random.default_rng(0)
for n in ("k562_gwps", "orion_hct116", "orion_hek293t"):
    h = DonorHandle(n); Sd, carried = h.rows(targets, genes); D[n] = (Sd.astype(np.float64), carried)
    # the donor's own dominant directions, from a 4000-perturbation subsample of ALL its rows, on the shared gene axis
    tcols, dcols = h.project(genes)
    sub = np.sort(rng.choice(h.n_perts, min(4000, h.n_perts), replace=False))
    M = np.zeros((sub.size, len(genes)))
    for r, i in enumerate(sub): M[r, tcols] = (np.asarray(h.delta[i], dtype=np.float64) - h.m)[dcols]
    U, s_, Vt = np.linalg.svd(M, full_matrices=False)
    PCS[n] = (Vt, s_**2 / np.sum(s_**2))
    print(f"{n:<15} variance share of top PCs: " + " ".join(f"{v:.3f}" for v in PCS[n][1][:6]) + f"   (top-1 {PCS[n][1][0]:.1%})")
allc = np.all([c for _, c in D.values()], axis=0)
m3 = pool & allc[:, None] & np.all([np.sign(Sd) != 0 for Sd, _ in D.values()], axis=0)
print(f"\nSIGN ACCURACY on the all-three pool (n={m3.sum()}) after deflating each donor's own top-k PCs:")
print(f"{'k':>3} {'K562':>8} {'HCT116':>8} {'HEK293T':>9} {'equal mean':>11} {'(2,1,1)':>8}")
for k in (0, 1, 2, 3, 5, 10, 20):
    defl = {}
    for n, (Sd, _) in D.items():
        Vt = PCS[n][0][:k]
        defl[n] = Sd - (Sd @ Vt.T) @ Vt if k > 0 else Sd
    cons = sum(defl.values()) / 3; w211 = (2*defl["k562_gwps"] + defl["orion_hct116"] + defl["orion_hek293t"]) / 4
    print(f"{k:>3} {acc(np.sign(defl['k562_gwps']), m3):>8.4f} {acc(np.sign(defl['orion_hct116']), m3):>8.4f} {acc(np.sign(defl['orion_hek293t']), m3):>9.4f} {acc(np.sign(cons), m3):>11.4f} {acc(np.sign(w211), m3):>8.4f}")
# does the top PC of each donor look like a shared 'toxicity' axis?  cosine between donors' top PCs
import itertools
print("\ncosine between donors' top-1 PCs (is the blob direction shared?):")
for a, b in itertools.combinations(PCS, 2):
    print(f"  {a} x {b}: {abs(float(PCS[a][0][0] @ PCS[b][0][0])):.3f}")
# how much of the 85 h1 targets' effect energy sits on each donor's own top PC?
print("\nshare of the h1-target rows' energy on the donor's own top-1 / top-3 PCs:")
for n, (Sd, _) in D.items():
    Vt = PCS[n][0]; e = (Sd**2).sum()
    print(f"  {n:<15} top-1 {((Sd @ Vt[:1].T)**2).sum()/e:.3f}   top-3 {((Sd @ Vt[:3].T)**2).sum()/e:.3f}")
