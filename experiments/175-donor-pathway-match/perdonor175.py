"""175 -- per-donor sign accuracy on the h1-test reference DE pool, overall and by target.
Hypothesis: donors differ in TP53 status (K562 null, HEK293T large-T inactivated, HCT116 WT);
H1 hESC is TP53-WT, so the p53-competent donor should transfer better on p53-mediated responses."""
import numpy as np, sys
sys.path.insert(0, "experiments/153-de-bench"); sys.path.insert(0, "experiments/172-sign-confidence-key")
from debench import Surface
from donors172 import DonorHandle
S = Surface("experiments/164-callset-artefact/de164.a15.npz")
pool = S.real_sig & S.adj; sT = np.sign(S.real_lfc)
genes = np.asarray([str(g) for g in S.genes]); targets = [str(p) for p in S.perts]
D = {}
for n in ("k562_gwps", "orion_hct116", "orion_hek293t"):
    h = DonorHandle(n); Sd, carried = h.rows(targets, genes)
    D[n] = (Sd.astype(np.float64), carried)
    print(f"{n:<15} carries {carried.sum()}/85 targets, {len(h.genes)} genes, {h.n_perts} perts")
def acc(key, m): 
    v = (key[m] == sT[m]) & (key[m] != 0); nz = (key[m] != 0)
    return v.sum() / max(nz.sum(), 1), int(nz.sum())
print("\nSIGN ACCURACY on the h1-test reference DE pool (where the donor has a nonzero entry):")
for n, (Sd, c) in D.items():
    a, k = acc(np.sign(Sd), pool & c[:, None]); print(f"  {n:<15} {a:.4f}  (n={k})")
cons = sum(Sd for Sd, _ in D.values()) / 3
a, k = acc(np.sign(cons), pool); print(f"  {'equal mean':<15} {a:.4f}  (n={k})")
# pairwise: all three carried
allc = np.all([c for _, c in D.values()], axis=0)
m3 = pool & allc[:, None] & np.all([np.sign(Sd) != 0 for Sd, _ in D.values()], axis=0)
print(f"\nrestricted to pairs where ALL THREE donors have a nonzero entry (n={m3.sum()}):")
for n, (Sd, _) in D.items(): print(f"  {n:<15} {np.mean(np.sign(Sd)[m3] == sT[m3]):.4f}")
print(f"  {'equal mean':<15} {np.mean(np.sign(cons)[m3] == sT[m3]):.4f}")
maj = np.sign(sum(np.sign(Sd) for Sd, _ in D.values()))
print(f"  {'majority vote':<15} {np.mean(maj[m3] == sT[m3]):.4f}")
un = np.all([np.sign(Sd) == np.sign(D['k562_gwps'][0]) for Sd, _ in D.values()], axis=0)
print(f"  {'unanimous only':<15} {np.mean(np.sign(cons)[m3 & un] == sT[m3 & un]):.4f}  (n={int((m3&un).sum())}, {100*(m3&un).sum()/m3.sum():.0f}% of pairs)")
# p53 axis: a canonical p53 target set -- is HCT116 better THERE specifically?
P53 = {"CDKN1A","MDM2","BAX","BBC3","PMAIP1","TP53I3","GADD45A","SESN1","SESN2","FAS","TNFRSF10B","RRM2B","DDB2","XPC","ZMAT3","TIGAR","PLK3","BTG2","FDXR","AEN","TRIAP1","CCNG1","PHLDA3","RPS27L","POLH","APAF1","PERP","SFN","TP53INP1","EDA2R","ACER2","SPATA18","KITLG","DRAM1","ISCU","PIDD1","TNFRSF10D","CD82","GDF15","TOB1"}
pm = np.isin(genes, list(P53)); print(f"\np53 target genes on the axis: {pm.sum()}; pool pairs on them: {(pool & pm[None,:]).sum()}")
print("sign accuracy on p53-TARGET pool pairs (all three carried):")
for n, (Sd, _) in D.items(): print(f"  {n:<15} {np.mean(np.sign(Sd)[m3 & pm[None,:]] == sT[m3 & pm[None,:]]):.4f}")
print(f"  {'equal mean':<15} {np.mean(np.sign(cons)[m3 & pm[None,:]] == sT[m3 & pm[None,:]]):.4f}")
# per-target: which donor wins, and by how much
print("\nper-target sign accuracy (pool, all-three-carried), targets sorted by (HCT116 - mean of others):")
rows = []
for p in range(85):
    m = m3[p]
    if m.sum() < 15: continue
    a = {n: np.mean(np.sign(Sd)[p, m] == sT[p, m]) for n, (Sd, _) in D.items()}
    rows.append((a["orion_hct116"] - (a["k562_gwps"] + a["orion_hek293t"]) / 2, targets[p], int(m.sum()), a, np.mean(np.sign(cons)[p, m] == sT[p, m])))
rows.sort(reverse=True)
print(f"  {'target':<10}{'n':>5}{'K562':>8}{'HCT116':>8}{'HEK293T':>9}{'mean':>7}")
for row in rows[:8] + [None] + rows[-6:]:
    if row is None: print("  ..."); continue
    d, t, n, a, c = row
    print(f"  {t:<10}{n:>5}{a['k562_gwps']:>8.3f}{a['orion_hct116']:>8.3f}{a['orion_hek293t']:>9.3f}{c:>7.3f}")
import statistics
print(f"\n  across {len(rows)} targets: HCT116 wins {sum(1 for r in rows if r[0]>0)}, loses {sum(1 for r in rows if r[0]<0)}; median advantage {statistics.median(r[0] for r in rows):+.4f}")
# oracle per-target donor selection: pick the best donor per target (upper bound on any routing)
best = np.mean([max(r[3].values()) for r in rows]); meanc = np.mean([r[4] for r in rows])
print(f"  per-target ORACLE best-single-donor mean acc {best:.4f}  vs equal-mean {meanc:.4f}  (routing ceiling)")

# weighted means: K562 heavier?
print("\nWEIGHTED consensus sign accuracy on the all-three pool (w_K562, w_HCT, w_HEK):")
K, Hc, He = D["k562_gwps"][0], D["orion_hct116"][0], D["orion_hek293t"][0]
for w in [(1,1,1),(2,1,1),(3,1,1),(1,0,0),(1,1,0),(1,0,1),(0,1,1),(2,1,0),(4,1,1)]:
    c = (w[0]*K + w[1]*Hc + w[2]*He) / sum(w)
    print(f"  {str(w):<12} {np.mean(np.sign(c)[m3] == sT[m3]):.4f}")
# is K562's edge uniform across targets, or concentrated?
adv = [(r[3]["k562_gwps"] - r[4], r[1], r[2]) for r in rows]
print(f"\nK562-alone minus equal-mean, per target: wins {sum(1 for a,_,_ in adv if a>0)}/{len(adv)}, median {statistics.median(a for a,_,_ in adv):+.4f}")
# per-gene magnitude: are the Orion lines noisier (smaller |effect| relative to spread)?
for n,(Sd,c) in D.items():
    v = Sd[allc][:, :]; print(f"  {n:<15} median |s| on pool {np.median(np.abs(Sd[pool & c[:,None]])):.4f}   frac of pool entries with |s|<0.01: {np.mean(np.abs(Sd[pool & c[:,None]])<0.01):.3f}")
