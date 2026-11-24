"""Q33 — we HAVE 0.705 sign accuracy in hand and ship 0.574.  Where is it lost?

Q32: HEK293T's measured effect predicts HCT116's sign at 0.7051 (0.7974 on the donor's strongest
quintile).  Our shipped prediction on the same pool: 0.5736.

The critical fact is that `orion_hek293t` **is one of the two donors the LODO prediction uses**
for HCT116.  So the pipeline is handed 0.705 and emits 0.574.  This is not a data ceiling and
not a paradigm limit -- it is 0.13 of sign accuracy destroyed in processing.

Decompose: donor table raw -> blended -> whitened -> shipped.
"""
import numpy as np
from debench import Surface

hct = Surface("../095-clean-harness/de_thinned.npz")
DD = "../../data/donors"

def donor_signacc(name):
    genes = [l.strip() for l in open(f"{DD}/{name}/genes.txt") if l.strip()]
    perts = [l.strip() for l in open(f"{DD}/{name}/perts.txt") if l.strip()]
    delta = np.load(f"{DD}/{name}/delta.npy", mmap_mode="r")
    gi = {g: i for i, g in enumerate(genes)}; ti = {t: i for i, t in enumerate(perts)}
    gcol = np.array([gi.get(g, -1) for g in hct.genes])
    prow = np.array([ti.get(t, -1) for t in hct.perts])
    num = den = 0; cov = 0
    for p in range(hct.perts.size):
        if prow[p] < 0: continue
        cov += 1
        pool = np.flatnonzero(hct.real_sig[p] & hct.adj[p] & (gcol >= 0))
        if pool.size < 10: continue
        d = np.asarray(delta[prow[p], gcol[pool]])
        a = np.sign(hct.real_lfc[p, pool]); b = np.sign(d)
        num += int(((a == b) & (b != 0)).sum()); den += int((b != 0).sum())
    return num / max(den, 1), den, cov

print(f"{'source':34} {'sign acc':>9} {'n pairs':>10} {'targets covered':>16}")
for name in ("orion_hek293t", "k562_gwps"):
    acc, n, cov = donor_signacc(name)
    print(f"{'donor table  ' + name:34} {acc:9.4f} {n:10,} {cov:16d}")

# the blended, un-whitened, and shipped objects, on the same pool
z = np.load("../104-context-heldout-credit-adapter/work/hct116/frozen.npz", allow_pickle=True)
g104 = {str(g): i for i, g in enumerate(z["genes"])}; t104 = {str(t): i for i, t in enumerate(z["targets"])}
gcol = np.array([g104.get(g, -1) for g in hct.genes]); prow = np.array([t104.get(t, -1) for t in hct.perts])
W = np.load("../062-cd4-post-whiten/cloud-output/whiteners/held_orion_hct116.npz")
ev, V = W["evals"], W["evecs"]
import csv
canon = [r["gene_name"] for r in csv.DictReader(open("../../data/challenge/gene_names.csv"))]
wc = np.load("/private/tmp/claude-501/-Users-ahmet-Code-virtual-cell-challenge-2026/"
             "aa534728-8641-48d7-a125-17e75d9c9399/scratchpad/wcols.npy")
wg = [canon[i] for i in wc]
base = z["baseline"][:, [g104[g] for g in wg]].astype(np.float32)
unw = ((base @ V) * (ev ** (0.65 / 2.0))).astype(np.float32) @ V.T     # remove the operator
wpos = {g: i for i, g in enumerate(wg)}

def acc_of(mat, colmap):
    num = den = 0
    for p in range(hct.perts.size):
        if prow[p] < 0: continue
        pool = np.flatnonzero(hct.real_sig[p] & hct.adj[p])
        pool = np.array([g for g in pool if hct.genes[g] in colmap])
        if pool.size < 10: continue
        b = np.sign(mat[prow[p], [colmap[hct.genes[g]] for g in pool]])
        a = np.sign(hct.real_lfc[p, pool])
        num += int(((a == b) & (b != 0)).sum()); den += int((b != 0).sum())
    return num / max(den, 1), den

a1, n1 = acc_of(unw, wpos)
print(f"{'104 baseline, operator REMOVED':34} {a1:9.4f} {n1:10,}")
a2, n2 = acc_of(base, wpos)
print(f"{'104 baseline, as shipped (whitened)':34} {a2:9.4f} {n2:10,}")
print(f"\nQ32 truth-to-truth ceiling (HEK293T measured -> HCT116): 0.7051")
