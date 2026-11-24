"""Route (b): per-perturbation confidence weighting.  How much of the oracle is buildable?

12-29 % of perturbations contribute NEGATIVE X and the top quartile carries 80-85 % of it
(096 §concentration).  The family is  dh'_p = z_p * dh_p , and because the global amplitude is
re-optimised onto the rho cliff afterwards, only the SHAPE of z matters:

    W(z) = (sum_p w_p z_p x_p)^2 / (sum_p w_p z_p^2 a_p),      x_p = <dh_p,d_p>, a_p = ||dh_p||^2

sqrt(W/W_base) is the multiplier on the deconfounded cosine and therefore on the mse gain.
Three ceilings, in increasing honesty:
    ORACLE-signed   z_p propto x_p/a_p, unconstrained (may flip a row's sign)
    ORACLE-damped   z_p propto max(0, x_p)/a_p, i.e. drop the wrong-way rows, never flip
    ORACLE-gate     z_p in {0,1} on sign(x_p): the crudest "drop the wrong-way ones"
and then the buildable keys, which never touch the reference.
"""
from __future__ import annotations
import sys, json
import numpy as np
from scipy.stats import spearmanr
ROOT = "~/Code/virtual-cell-challenge-2026"
sys.path.insert(0, f"{ROOT}/experiments/096-mse-floor")
sys.path.insert(0, f"{ROOT}/src")
from decompose import terms
from msekit import group_stats
from vcc2026.donor_transfer import Donor

HARNESS = [
    ("h1-test", f"{ROOT}/data/harness/h1-test.h5ad",
     f"{ROOT}/experiments/095-clean-harness/out-h1test/runs/g0.65_q1.0_a1.08_spread0.25/pred.h5ad"),
    ("orion", f"{ROOT}/data/harness/orion-hct116-thinned.h5ad",
     f"{ROOT}/experiments/095-clean-harness/out/runs/g0.65_q1.0_a1.08_spread0.25/pred.h5ad"),
]
DONORS = ("k562_gwps", "orion_hek293t")   # exactly what the harness model carries


def load(real, pred):
    rp, rs, rc, rj, rg = group_stats(real)
    pp, ps, pc, pj, pg = group_stats(pred)
    t = terms(pp, ps, pc, pj, rp, rs, rc, rj, rg)
    for k, c in t["cols"].items():
        t["dh"][k, c] = 0.0
        t["d"][k, c] = 0.0
    t["genes"] = rg
    t["targets"] = np.array([str(x) for x in pp[t["keep"]]])
    return t


def donor_rows(genes, targets):
    """Per-donor, per-target effect vectors on the harness axis, renormalised to a common
    median norm exactly as MultiDonorModel does, plus the per-target availability mask."""
    out, have = {}, {}
    for name in DONORS:
        d = Donor.load(f"{ROOT}/data/donors/{name}")
        m, s = d.decompose()
        gi = {g: i for i, g in enumerate(d.genes)}
        cols = np.array([gi.get(str(g), -1) for g in genes])
        ok = cols >= 0
        pi = {p: i for i, p in enumerate(d.perts)}
        M = np.zeros((len(targets), genes.size))
        hv = np.zeros(len(targets), bool)
        for r, tg in enumerate(targets):
            if tg in pi:
                M[r, ok] = s[pi[tg]][cols[ok]]
                hv[r] = True
        scale = float(np.median(np.linalg.norm(s[:, cols[ok]], axis=1))) or 1.0
        out[name] = M / scale
        have[name] = hv
        del d, s, M
    return out, have


def W_of(z, x, a, w):
    num = float((w * z * x).sum())
    den = float((w * z * z * a).sum())
    return (num * num / den) if den > 0 else 0.0


def build_scores(t):
    """Every per-perturbation confidence score that can be computed WITHOUT the reference."""
    dh, genes, targets = t["dh"], t["genes"], t["targets"]
    P = dh.shape[0]
    S, meta = {}, {}
    nh = np.linalg.norm(dh, axis=1)
    S["self_norm"] = nh
    U = dh / np.maximum(nh[:, None], 1e-12)
    Cm = U @ U.T
    np.fill_diagonal(Cm, -np.inf)
    S["uniqueness"] = 1.0 - Cm.max(1)

    D, have = donor_rows(genes, targets)
    A_, B_ = D[DONORS[0]], D[DONORS[1]]
    hA, hB = have[DONORS[0]], have[DONORS[1]]
    S["n_donors"] = (hA.astype(float) + hB.astype(float))
    meta["coverage"] = {DONORS[0]: int(hA.sum()), DONORS[1]: int(hB.sum()),
                        "both": int((hA & hB).sum()), "n": P}

    # 096's key: sum of sign-gated geometric mean.  Magnitude-confounded by construction.
    gate = (np.sign(A_) == np.sign(B_)) & (A_ != 0) & (B_ != 0)
    agree_vec = np.sqrt(np.abs(A_) * np.abs(B_)) * gate
    both = hA & hB
    S["donor096"] = np.where(both, agree_vec.sum(1), 0.0)

    # NEW, scale-free variants -- the point is to strip the magnitude confound out of the key
    na, nb = np.linalg.norm(A_, axis=1), np.linalg.norm(B_, axis=1)
    S["donor_cos"] = np.where(both, np.einsum("ij,ij->i", A_, B_) / np.maximum(na * nb, 1e-12), 0.0)
    # sign concordance over the genes that BOTH donors move appreciably
    frac = np.zeros(P)
    for i in range(P):
        if not both[i]:
            continue
        k = max(1, int(0.02 * genes.size))
        idx = np.union1d(np.argpartition(-np.abs(A_[i]), k)[:k],
                         np.argpartition(-np.abs(B_[i]), k)[:k])
        frac[i] = float(np.mean(np.sign(A_[i][idx]) == np.sign(B_[i][idx])))
    S["donor_signfrac"] = frac
    # agreement mass RELATIVE to the perturbation's own predicted size: a purity, not a total
    S["donor096_norm"] = np.where(both, S["donor096"] / np.maximum(na * nb, 1e-12), 0.0)
    # how much of OUR OWN prediction the donors' agreement explains
    av = agree_vec / np.maximum(np.linalg.norm(agree_vec, axis=1, keepdims=True), 1e-12)
    S["agree_vs_pred"] = np.where(both, np.abs(np.einsum("ij,ij->i", av, U)), 0.0)
    return S, meta


def main():
    out = {}
    for name, real, pred in HARNESS:
        t = load(real, pred)
        dh, d, w = t["dh"], t["d"], t["w"]
        x = np.einsum("ij,ij->i", dh, d)
        a = np.einsum("ij,ij->i", dh, dh)
        cos_p = x / np.maximum(np.sqrt(a) * np.linalg.norm(d, axis=1), 1e-30)
        P = len(x)
        base = W_of(np.ones(P), x, a, w)
        print(f"\n########## {name}  n={P}  base W={base:.6g}  "
              f"deconfounded cos={t['X']/np.sqrt(t['A']*t['sumD']):.5f} ##########")
        print(f"  rows with x_p<=0: {100*np.mean(x<=0):.1f} %   they carry "
              f"{100*(w*x)[x<=0].sum()/(w*x).sum():+.1f} % of X (i.e. removing them alone "
              f"raises X by {-100*(w*x)[x<=0].sum()/(w*x).sum():+.1f} %)")
        oracles = {
            "ORACLE signed  z=x/a": np.where(a > 0, x / np.maximum(a, 1e-30), 0.0),
            "ORACLE damped  z=max(0,x)/a": np.where(a > 0, np.maximum(x, 0) / np.maximum(a, 1e-30), 0.0),
            "ORACLE gate    z=1[x>0]": (x > 0).astype(float),
            "ORACLE topq    z=1[top 25% by x/a]": None,
        }
        r = np.where(a > 0, x / np.maximum(a, 1e-30), 0.0)
        thr = np.quantile(r, 0.75)
        oracles["ORACLE topq    z=1[top 25% by x/a]"] = (r >= thr).astype(float)
        print(f"\n  {'ceiling':<38}{'sqrt(W/W_base)':>16}")
        for lbl, z in oracles.items():
            print(f"  {lbl:<38}{np.sqrt(W_of(z, x, a, w)/base):>15.3f}x")

        S, meta = build_scores(t)
        print(f"\n  donor coverage of the {P} targets: {json.dumps(meta['coverage'])}")
        print(f"\n  {'buildable score':<22}{'rho w/ cos_p':>13}{'partial rho':>13}"
              f"{'best gamma':>11}{'sqrt(W/Wb) in-sample':>21}")
        # partial out the TRUE effect size, as 096 §partial does
        nd = np.linalg.norm(d, axis=1)
        def partial(u, v, ctrl):
            from scipy.stats import rankdata
            U_, V_, C_ = (rankdata(q).astype(float) for q in (u, v, ctrl))
            U_ -= U_.mean(); V_ -= V_.mean(); C_ -= C_.mean()
            U_ -= C_ * (U_ @ C_) / (C_ @ C_); V_ -= C_ * (V_ @ C_) / (C_ @ C_)
            return float(U_ @ V_ / np.sqrt((U_ @ U_) * (V_ @ V_)))
        res = {}
        GAM = [0.0, 0.25, 0.5, 0.75, 1.0, 1.5, 2.0, 3.0]
        for lbl, sc in S.items():
            rk = (np.argsort(np.argsort(sc)).astype(float) + 0.5) / len(sc)   # rank in (0,1)
            mults = [np.sqrt(W_of(rk ** g, x, a, w) / base) for g in GAM]
            j = int(np.argmax(mults))
            res[lbl] = dict(gam=GAM, mult=[float(m) for m in mults],
                            rho=float(spearmanr(sc, cos_p).statistic),
                            prho=partial(sc, cos_p, nd))
            print(f"  {lbl:<22}{res[lbl]['rho']:>+13.3f}{res[lbl]['prho']:>+13.3f}"
                  f"{GAM[j]:>11.2f}{mults[j]:>20.3f}x")
        out[name] = dict(res=res, base=base, GAM=GAM, x=x.tolist(), a=a.tolist(),
                         w=w.tolist(), scores={k: v.tolist() for k, v in S.items()},
                         cos_p=cos_p.tolist(), targets=t["targets"].tolist())
    with open(f"{ROOT}/experiments/097-direction/work/perpert.json", "w") as f:
        json.dump(out, f)

    # ---- cross-surface: the gamma tuned on one harness, applied to the other -------------
    print("\n\n########## CROSS-SURFACE (gamma fitted on the OTHER harness) ##########")
    print(f"  {'score':<22}{'gam(fit h1)':>12}{'-> orion':>10}"
          f"{'gam(fit orion)':>16}{'-> h1':>9}{'transferable':>14}")
    A, B = out["h1-test"], out["orion"]
    for lbl in A["res"]:
        ga = A["res"][lbl]["gam"][int(np.argmax(A["res"][lbl]["mult"]))]
        gb = B["res"][lbl]["gam"][int(np.argmax(B["res"][lbl]["mult"]))]
        mb = B["res"][lbl]["mult"][B["res"][lbl]["gam"].index(ga)]
        ma = A["res"][lbl]["mult"][A["res"][lbl]["gam"].index(gb)]
        print(f"  {lbl:<22}{ga:>12.2f}{mb:>10.3f}{gb:>16.2f}{ma:>9.3f}"
              f"{min(ma, mb):>13.3f}x")


if __name__ == "__main__":
    main()
