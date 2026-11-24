"""Miller et al. 2025 (bioRxiv 2025.10.20.683304) per-gene shrinkage, in the (star) framework.

    mu_ID = alpha (*) mu_TD + (1 - alpha) (*) mu_all,   alpha_pg = 1 - p_DEG(p, g)

In effect coordinates (everything relative to b_ctrl) that is
    dh_p' = alpha (*) dh_p + (1 - alpha) (*) dh_bar
i.e. it shrinks each perturbation's DEVIATION FROM THE PANEL MEAN, keeping the common mode --
which is not the same as shrinking toward the control, and is exactly why it interacts with B.

Evaluated the way the metric actually pays for it: re-optimise the global scale onto the rho
cliff AFTER shrinking, so A', B' and X' all move together.
    s_c = sqrt(claim/B'),  MSE = 1 + (s_c^2 (A'-B') + Kc - 2 s_c X')/sumD
Reported as the multiplier on sqrt(W), W = X'^2/A', which is what the master formula scales by.

alpha must be BUILDABLE (rule 4). Two p_DEG proxies:
  self   -- per-perturbation rank of |dh_pg|, our own model's confidence, no reference
  donor  -- signed agreement between k562_gwps and orion_hek293t, genuinely external
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

def W_of(dh, d, w):
    n = float((np.einsum("ij,ij->i", dh, d) * w).sum())
    q = float((np.einsum("ij,ij->i", dh, dh) * w).sum())
    return (n * n / q if q > 0 else 0.0), q, n

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
            out[r] = np.sqrt(np.abs(a) * np.abs(b)) * ((np.sign(a) == np.sign(b)) & (a != 0) & (b != 0))
    return out

for real, pred, name in [
    (f"{ROOT}/data/harness/h1-test.h5ad",
     f"{ROOT}/experiments/095-clean-harness/out-h1test/runs/g0.65_q1.0_a1.08_spread0.25/pred.h5ad", "h1-test"),
    (f"{ROOT}/data/harness/orion-hct116-thinned.h5ad",
     f"{ROOT}/experiments/095-clean-harness/out/runs/g0.65_q1.0_a1.08_spread0.25/pred.h5ad", "orion-thinned")]:
    t, dh, d, w, genes, targets = load(real, pred)
    base, _, _ = W_of(dh, d, w)
    dbar = dh.mean(0)
    self_rank = np.argsort(np.argsort(-np.abs(dh), axis=1), axis=1)     # 0 = strongest
    dk = donor_key(genes, targets)
    donor_rank = np.argsort(np.argsort(-dk, axis=1), axis=1)
    print(f"\n===== {name}: {dh.shape[0]} perts, {dh.shape[1]} genes. W_global = {base:.6g} =====")
    print(f"  donor key covers {int((dk.sum(1) > 0).sum())}/{len(targets)} perturbations")
    print(f"  {'family':<44}{'W':>13}{'sqrt(W/W_glob)':>17}")
    for key, rank in (("self |dh| rank", self_rank), ("donor agreement rank", donor_rank)):
        for K in (500, 2000, 8000):
            a = (rank < K).astype(np.float64)                 # hard alpha: 1 on top-K, 0 else
            v, _, _ = W_of(a * dh + (1 - a) * dbar[None, :], d, w)
            print(f"  {f'Miller hard, {key}, K={K}':<44}{v:>13.6g}{np.sqrt(v/base):>17.3f}")
        for g in (0.5, 1.0):                                   # soft alpha: exp decay in rank
            a = np.exp(-g * rank / 2000.0)
            v, _, _ = W_of(a * dh + (1 - a) * dbar[None, :], d, w)
            print(f"  {f'Miller soft, {key}, gamma={g}':<44}{v:>13.6g}{np.sqrt(v/base):>17.3f}")
    # oracle ceiling for the same family shape
    a = (np.sign(dh) == np.sign(d)).astype(np.float64)
    v, _, _ = W_of(a * dh + (1 - a) * dbar[None, :], d, w)
    print(f"  {'ORACLE alpha = 1[sign(dh)==sign(d)]  (leak)':<44}{v:>13.6g}{np.sqrt(v/base):>17.3f}")
