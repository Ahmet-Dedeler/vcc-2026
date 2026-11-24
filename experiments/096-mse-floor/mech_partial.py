"""Is kg_density predicting our ACCURACY, or just TRUE EFFECT MAGNITUDE wearing a graph costume?

Dense STRING neighbourhoods belong to well-studied hub genes, which are disproportionately
essential, and essential genes have large effects (TRADE: ~45 genes moved typically, >500 for
essentials). Raju reports per-perturbation directional coherence vs effect magnitude at Spearman
0.977. So a graph-density score could be a magnitude proxy rather than an accuracy signal -- and
those are different objects for a confidence weight.

Test: partial Spearman of each key with the true per-perturbation cosine, controlling for the
TRUE effect magnitude ||d_p||.
"""
import sys
import numpy as np
import h5py
from scipy.stats import spearmanr
ROOT = "~/Code/virtual-cell-challenge-2026"
sys.path.insert(0, f"{ROOT}/experiments/096-mse-floor")
from decompose import terms
from msekit import group_stats
from mech_second_key import donor_key, EMB, K

def partial(x, y, z):
    """Spearman r_xy.z"""
    rxy = spearmanr(x, y).statistic; rxz = spearmanr(x, z).statistic; ryz = spearmanr(y, z).statistic
    return (rxy - rxz*ryz)/np.sqrt(max((1-rxz**2)*(1-ryz**2), 1e-12))

for real, pred, name in [
    (f"{ROOT}/data/harness/h1-test.h5ad",
     f"{ROOT}/experiments/095-clean-harness/out-h1test/runs/g0.65_q1.0_a1.08_spread0.25/pred.h5ad", "h1-test"),
    (f"{ROOT}/data/harness/orion-hct116-thinned.h5ad",
     f"{ROOT}/experiments/095-clean-harness/out/runs/g0.65_q1.0_a1.08_spread0.25/pred.h5ad", "orion-thinned")]:
    rp, rs, rc, rj, rg = group_stats(real); pp, ps, pc, pj, pg = group_stats(pred)
    t = terms(pp, ps, pc, pj, rp, rs, rc, rj, rg)
    dh, d = t["dh"], t["d"]
    for k, c in t["cols"].items(): dh[k, c] = 0.0; d[k, c] = 0.0
    targets = [str(x) for x in pp[t["keep"]]]
    P = dh.shape[0]
    nh = np.linalg.norm(dh, axis=1); nd = np.linalg.norm(d, axis=1)
    cos_p = np.einsum("ij,ij->i", dh, d)/np.maximum(nh*nd, 1e-12)
    E = np.zeros((P, 64)); have = np.array([x in EMB for x in targets])
    for i, x in enumerate(targets):
        if have[i]: E[i] = EMB[x]
    En = E/np.maximum(np.linalg.norm(E, axis=1, keepdims=True), 1e-12)
    S = En@En.T; np.fill_diagonal(S, -np.inf)
    kg = np.sort(S, axis=1)[:, -K:].mean(1)
    donor = donor_key(rg, targets).sum(1)
    print(f"\n===== {name}: {P} perturbations =====")
    print(f"  true effect magnitude ||d_p||   vs cos_p : spearman {spearmanr(nd, cos_p).statistic:+.3f}")
    print(f"  {'key':<26}{'r with cos_p':>14}{'r with ||d_p||':>16}{'PARTIAL r.||d||':>18}")
    for lbl, sc in (("kg_density (STRING)", kg), ("donor agreement", donor),
                    ("||dh_p|| (self, control)", nh)):
        print(f"  {lbl:<26}{spearmanr(sc, cos_p).statistic:>14.3f}"
              f"{spearmanr(sc, nd).statistic:>16.3f}{partial(sc, cos_p, nd):>18.3f}")
