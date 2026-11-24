"""099 arm 4 -- does adding donors improve the AGREEMENT KEY, holding the PREDICTION fixed?

The prediction and the key are separable and both matter, so they are measured apart. Here the
predicted effect matrix `dh` is FROZEN at the harness's own legal base-2 build, and only the
confidence key changes:

  per gene         `alpha_g = rank01(unanimity key)^eta`     -- production's `donor_agreement`
  per perturbation `z_p     = max(rank01(donor_cos_p), 0.3)^gamma`

Reported per donor subset:
  graded fraction   median share of the 18,533 genes the unanimity key scores non-zero
  rho               Spearman of the per-perturbation key against the TRUE per-perturbation
                    cosine -- how good a predictor the key is
  M                 kappa(weighted dh) / kappa(dh) -- how good a WEIGHT it is; these are
                    different functionals and 097 measured them moving in opposite directions
                    when the third donor was added.
"""
from __future__ import annotations
import argparse, itertools, json, sys
import numpy as np
from scipy.stats import rankdata, spearmanr
sys.path.insert(0, "~/Code/virtual-cell-challenge-2026/experiments/099-donors")
from dkit import Ctx, ROOT, donor_matrices, unanimity_key, pairwise_key, donor_cos_key

POOL = ("k562_gwps", "orion_hek293t", "orion_hct116", "kolf2",
        "cd4_rest", "cd4_stim8hr", "cd4_stim48hr", "perturbmatch")


def subsets(base, extras):
    out = [("D2 " + "+".join(b[:6] for b in base), base)]
    for e in extras:
        out.append((f"D3 +{e}", base + (e,)))
    cd4 = tuple(c for c in ("cd4_rest", "cd4_stim8hr", "cd4_stim48hr") if c in extras)
    if len(cd4) == 3:
        out.append(("D5 +cd4x3", base + cd4))
        if "kolf2" in extras:
            out.append(("D6 +kolf2+cd4x3", base + ("kolf2",) + cd4))
    if "kolf2" in extras and "cd4_rest" in extras:
        out.append(("D4 +kolf2+cd4_rest", base + ("kolf2", "cd4_rest")))
    if "kolf2" in extras and "orion_hct116" in extras:
        out.append(("D4 +hct116+kolf2", base + ("orion_hct116", "kolf2")))
    if "orion_hct116" in extras and len(cd4) == 3:
        out.append(("D6 +hct116+cd4x3", base + ("orion_hct116",) + cd4))
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("harness")
    ap.add_argument("--gamma-build", type=float, default=0.65,
                    help="whitening gamma of the FROZEN prediction the key is applied to")
    a = ap.parse_args()
    C = Ctx(a.harness)
    pool = tuple(p for p in POOL if p not in C.banned)
    base = tuple(d for d in ("k562_gwps", "orion_hek293t", "orion_hct116")
                 if d not in C.banned)[:2]
    m = C.set_donors(pool)
    cov = C.coverage()
    print(f"# {a.harness}: prediction FROZEN at base2={base}, whitening gamma={a.gamma_build}")
    print("# coverage: " + ", ".join(f"{n} {cov[n]}/{C.P}" for n in pool))

    w0 = {n: (1.0 if n in base else 0.0) for n in pool}
    dh = C.dh_exact(C.effects(lambda _p: dict(w0), gamma=a.gamma_build))
    k0 = C.kappa(dh)
    n0 = C.nmae(dh)
    truecos = C.per_pert_cos(dh)
    M, H = donor_matrices(m, C.perts)
    P, G = C.P, C.G

    def alpha_of(agr, eta):
        al = np.empty_like(agr)
        for i in range(P):
            al[i] = (rankdata(agr[i], method="average") / G) ** eta
            if np.unique(agr[i]).size < 2:
                al[i] = 1.0
        return al

    print(f"\n# baseline: kappa={k0:.5f}  nmae={n0:.5f}")
    print(f"\n  {'donor set':<22}{'nd_med':>7}{'graded':>8}{'rho_cos':>9}"
          f"{'M_gene':>8}{'M_pert':>8}{'M_both':>8}{'nmae_both':>10}")
    out = {}
    for lbl, dons in subsets(base, tuple(e for e in pool if e not in base)):
        agr = unanimity_key(M, H, dons, P, G)
        dc, nd = donor_cos_key(M, H, dons, P)
        frac = float(np.median([(agr[i] > 0).mean() for i in range(P)]))
        rho = float(spearmanr(dc, truecos).statistic)
        al = alpha_of(agr, 1.0)
        z = np.maximum((rankdata(dc, method="average") - 0.5) / P, 0.3) ** 1.0
        m_gene = C.kappa(al * dh) / k0
        m_pert = C.kappa(z[:, None] * dh) / k0
        dh_both = z[:, None] * al * dh
        m_both = C.kappa(dh_both) / k0
        nm = C.nmae(dh_both)
        out[lbl] = dict(donors=list(dons), graded=frac, rho=rho, m_gene=float(m_gene),
                        m_pert=float(m_pert), m_both=float(m_both), nmae=float(nm),
                        nd_med=float(np.median(nd)))
        print(f"  {lbl:<22}{np.median(nd):>7.1f}{frac:>8.3f}{rho:>9.3f}"
              f"{m_gene:>8.3f}{m_pert:>8.3f}{m_both:>8.3f}{nm:>10.5f}")
    p = ROOT / f"experiments/099-donors/work/keys_{a.harness}.json"
    json.dump(dict(harness=a.harness, base=list(base), kappa0=k0, nmae0=n0,
                   coverage=cov, sets=out), open(p, "w"), indent=1)
    print(f"\n# nmae baseline {n0:.5f}; lower is better.  -> {p}")


if __name__ == "__main__":
    main()
