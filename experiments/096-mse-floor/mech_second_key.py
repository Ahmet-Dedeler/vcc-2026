"""Is it DONOR agreement that works, or EXTERNAL agreement? Test with a second, independent key.

The donor key (k562_gwps x orion_hek293t) is the only confidence score in this experiment that has
ever transferred. Two readings of that, and they imply different next steps:
  (a) donor screens specifically carry the signal  -> tune the donor key harder
  (b) any external prior does                      -> add more keys, they may be complementary

A second key that is external but NOT a perturbation screen separates them. STRING network
embeddings are already pinned in this repo (data/external/slim/gene_string_embeddings.v0.3.h5,
21,688 symbols x 64 features) by experiment 087, which preregistered this exact ablation and was
never run.

Two buildable per-perturbation scores from it, neither touching any reference:
  kg_density   mean cosine similarity of p's target gene to its k nearest STRING neighbours
               -- "is this gene in a well-characterised module"
  kg_consensus cosine between our predicted effect for p and the mean predicted effect of p's
               STRING neighbours -- "does our prediction for p look like our predictions for
               functionally similar genes"
"""
import sys
import numpy as np
import h5py
from scipy.stats import spearmanr
ROOT = "~/Code/virtual-cell-challenge-2026"
sys.path.insert(0, f"{ROOT}/experiments/096-mse-floor")
from decompose import terms
from msekit import group_stats

K = 15

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

with h5py.File(f"{ROOT}/data/external/slim/gene_string_embeddings.v0.3.h5", "r") as h:
    EMB = {k: np.asarray(h[k][:], dtype=np.float64) for k in h.keys()}
print(f"STRING embeddings: {len(EMB)} symbols x {len(next(iter(EMB.values())))} features")

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
    cos_p = np.einsum("ij,ij->i", dh, d)/np.maximum(nh*nd, 1e-12)
    have = np.array([t_ in EMB for t_ in targets])
    E = np.zeros((P, 64))
    for i, t_ in enumerate(targets):
        if have[i]: E[i] = EMB[t_]
    En = E/np.maximum(np.linalg.norm(E, axis=1, keepdims=True), 1e-12)
    S = En@En.T; np.fill_diagonal(S, -np.inf)
    kg_density = np.sort(S, axis=1)[:, -K:].mean(1)
    U = dh/np.maximum(nh[:, None], 1e-12)
    kg_cons = np.zeros(P)
    for i in range(P):
        nb = np.argsort(-S[i])[:K]
        nb = nb[have[nb]]
        if nb.size and have[i]:
            m = dh[nb].mean(0); mn = np.linalg.norm(m)
            if mn > 0: kg_cons[i] = float(U[i]@(m/mn))
    donor = donor_key(genes, targets).sum(1)
    print(f"\n===== {name}: {P} perturbations, {int(have.sum())} with STRING embeddings =====")
    base = W_of(dh, d, w)
    print(f"  {'score':<34}{'spearman r':>12}{'W mult':>10}")
    scores = {"donor agreement (screens)": donor,
              "kg_density (STRING, external)": kg_density,
              "kg_consensus (STRING, external)": kg_cons}
    for lbl, sc in scores.items():
        z = (sc-sc.min())/max(sc.max()-sc.min(), 1e-12)
        print(f"  {lbl:<34}{spearmanr(sc, cos_p).statistic:>12.3f}"
              f"{np.sqrt(W_of(z[:, None]*dh, d, w)/base):>9.3f}x")
    # additivity: does the STRING key add to the donor key, or duplicate it?
    zd = (donor-donor.min())/max(donor.max()-donor.min(), 1e-12)
    for lbl, sc in (("kg_density", kg_density), ("kg_consensus", kg_cons)):
        zk = (sc-sc.min())/max(sc.max()-sc.min(), 1e-12)
        comb = 0.5*(zd+zk)
        print(f"  {'donor + '+lbl+' (mean)':<34}{spearmanr(comb, cos_p).statistic:>12.3f}"
              f"{np.sqrt(W_of(comb[:, None]*dh, d, w)/base):>9.3f}x"
              f"   [corr with donor key: {spearmanr(zk, zd).statistic:+.3f}]")
