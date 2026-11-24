"""Is our directional accuracy CONCENTRATED in a few perturbations, and is that concentration
visible WITHOUT the reference?

The intel session's open tension: we hold the board's best pds_cosine (0.8318) while being behind
by >=28% on aggregate cosine. Those are different statistics -- a normalised RANK of the correct
match versus bulk directional agreement -- so there is no contradiction, but the combination
suggests our accuracy may be thick in some perturbations and thin in others. If so, a
per-perturbation confidence weighting has real headroom, and it is the same object as Miller's
alpha. Tested here.

Candidate confidence scores, all buildable at submission time (no reference):
  self_norm   ||dh_p||                        our own predicted effect size
  uniq        min_q!=p (1 - cos(dh_p, dh_q))  how distinct our prediction for p is from the rest
  donor       mean signed donor agreement     external, k562_gwps x orion_hek293t
"""
import sys
import numpy as np
ROOT = "~/Code/virtual-cell-challenge-2026"
sys.path.insert(0, f"{ROOT}/experiments/096-mse-floor")
from decompose import terms
from msekit import group_stats

def load(real, pred):
    rp, rs, rc, rj, rg = group_stats(real); pp, ps, pc, pj, pg = group_stats(pred)
    t = terms(pp, ps, pc, pj, rp, rs, rc, rj, rg)
    dh, d, w = t["dh"], t["d"], t["w"]
    for k, c in t["cols"].items(): dh[k, c] = 0.0; d[k, c] = 0.0
    return t, dh, d, w, rg, [str(x) for x in pp[t["keep"]]]

def donor_key(genes, targets):
    DON = {}
    for name in ("k562_gwps", "orion_hek293t"):
        root = f"{ROOT}/data/donors/{name}"
        dg = np.array(open(f"{root}/genes.txt").read().split())
        dp = np.array(open(f"{root}/perts.txt").read().split())
        gi = {g: i for i, g in enumerate(dg)}
        DON[name] = (np.load(f"{root}/delta.npy", mmap_mode="r"),
                     np.array([gi.get(g, -1) for g in genes]), {p: i for i, p in enumerate(dp)})
    out = np.zeros((len(targets), genes.size))
    for r, t in enumerate(targets):
        vs = []
        for n in DON:
            delta, cols, pi = DON[n]
            if t not in pi: vs = []; break
            row = np.asarray(delta[pi[t]], float); v = np.zeros(genes.size); ok = cols >= 0
            v[ok] = row[cols[ok]]; vs.append(v)
        if len(vs) == 2:
            a, b = vs
            out[r] = np.sqrt(np.abs(a)*np.abs(b)) * ((np.sign(a)==np.sign(b)) & (a!=0) & (b!=0))
    return out

def W_of(dh, d, w):
    n = float((np.einsum("ij,ij->i", dh, d)*w).sum()); q = float((np.einsum("ij,ij->i", dh, dh)*w).sum())
    return n*n/q if q > 0 else 0.0

for real, pred, name in [
    (f"{ROOT}/data/harness/h1-test.h5ad",
     f"{ROOT}/experiments/095-clean-harness/out-h1test/runs/g0.65_q1.0_a1.08_spread0.25/pred.h5ad", "h1-test"),
    (f"{ROOT}/data/harness/orion-hct116-thinned.h5ad",
     f"{ROOT}/experiments/095-clean-harness/out/runs/g0.65_q1.0_a1.08_spread0.25/pred.h5ad", "orion-thinned")]:
    t, dh, d, w, genes, targets = load(real, pred)
    P = dh.shape[0]
    nh = np.linalg.norm(dh, axis=1); nd = np.linalg.norm(d, axis=1)
    cos_p = np.einsum("ij,ij->i", dh, d) / np.maximum(nh*nd, 1e-12)
    agg = t["X"]/np.sqrt(t["A"]*t["sumD"])
    print(f"\n===== {name}: {P} perturbations.  aggregate cos {agg:.4f} =====")
    q = np.percentile(cos_p, [5,25,50,75,95])
    print(f"  per-perturbation cosine: mean {cos_p.mean():+.4f}  sd {cos_p.std():.4f}  "
          f"p5/25/50/75/95 {q[0]:+.3f}/{q[1]:+.3f}/{q[2]:+.3f}/{q[3]:+.3f}/{q[4]:+.3f}")
    print(f"  fraction with cos <= 0: {100*(cos_p<=0).mean():.1f} %   "
          f"top-quartile share of total X: "
          f"{100*np.sort(np.einsum('ij,ij->i',dh,d)*w)[-P//4:].sum()/ (np.einsum('ij,ij->i',dh,d)*w).sum():.1f} %")
    # buildable confidence scores
    U = dh/np.maximum(nh[:,None],1e-12)
    C = U@U.T; np.fill_diagonal(C, -np.inf)
    uniq = 1.0 - C.max(1)
    dk = donor_key(genes, targets); donor = dk.sum(1)
    print(f"  {'confidence score':<22}{'spearman r with cos_p':>24}{'W mult if used as weight':>26}")
    from scipy.stats import spearmanr
    base = W_of(dh, d, w)
    for lbl, sc in (("||dh_p||  (self)", nh), ("uniqueness (self)", uniq),
                    ("donor agreement", donor), ("ORACLE cos_p (leak)", cos_p)):
        r = spearmanr(sc, cos_p).statistic
        z = (sc - sc.min())/max(sc.max()-sc.min(), 1e-12)
        v = W_of(z[:,None]*dh, d, w)
        print(f"  {lbl:<22}{r:>24.3f}{np.sqrt(v/base):>25.3f}x")
