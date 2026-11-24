"""099 arm 4b -- the key, with the two checks arm 4 was missing.

(1) PAIRWISE as well as UNANIMITY. Unanimity is an INTERSECTION over donor pairs, so its graded
    fraction collapses as donors are added (0.206 at two donors -> 0.015 at six on orion).
    Pairwise is a UNION and does the opposite. 097 measured unanimity > pairwise at THREE
    donors; whether that survives five or six is a different question.

(2) A PERMUTATION NULL on the added donor. Adding any fourth vector shrinks the graded
    fraction, and 097 measured that grading fewer genes helps by itself. Shuffling the added
    donor's perturbation map keeps its coverage, magnitudes and graded-fraction effect and
    destroys only the target correspondence, so the null isolates the information.
"""
from __future__ import annotations
import argparse, json, sys
import numpy as np
from scipy.stats import rankdata, spearmanr
sys.path.insert(0, "~/Code/virtual-cell-challenge-2026/experiments/099-donors")
from dkit import Ctx, ROOT, donor_matrices, unanimity_key, pairwise_key, donor_cos_key

POOL = ("k562_gwps", "orion_hek293t", "orion_hct116", "kolf2",
        "cd4_rest", "cd4_stim8hr", "cd4_stim48hr", "perturbmatch")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("harness")
    ap.add_argument("--n-null", type=int, default=40)
    ap.add_argument("--gamma-build", type=float, default=0.65)
    ap.add_argument("--seed", type=int, default=0)
    a = ap.parse_args()
    C = Ctx(a.harness)
    pool = tuple(p for p in POOL if p not in C.banned)
    base = tuple(d for d in ("k562_gwps", "orion_hek293t", "orion_hct116")
                 if d not in C.banned)[:2]
    m = C.set_donors(pool)
    cov = C.coverage()
    w0 = {n: (1.0 if n in base else 0.0) for n in pool}
    dh = C.dh_exact(C.effects(lambda _p: dict(w0), gamma=a.gamma_build))
    k0 = C.kappa(dh)
    truecos = C.per_pert_cos(dh)
    M, H = donor_matrices(m, C.perts)
    P, G = C.P, C.G
    print(f"# {a.harness}: prediction FROZEN at base2={base}; kappa0={k0:.5f}")
    print("# coverage: " + ", ".join(f"{n} {cov[n]}/{C.P}" for n in pool))

    def alpha_of(agr):
        al = np.empty_like(agr)
        for i in range(P):
            al[i] = rankdata(agr[i], method="average") / G
            if np.unique(agr[i]).size < 2:
                al[i] = 1.0
        return al

    def score(dons, keyfn):
        agr = keyfn(M, H, dons, P, G)
        dc, _ = donor_cos_key(M, H, dons, P)
        z = np.maximum((rankdata(dc, method="average") - 0.5) / P, 0.3)
        frac = float(np.median([(agr[i] > 0).mean() for i in range(P)]))
        rho = float(spearmanr(dc, truecos).statistic)
        mb = C.kappa(z[:, None] * alpha_of(agr) * dh) / k0
        return frac, rho, float(mb)

    extras = [e for e in pool if e not in base]
    sets = [("D2", base)] + [(f"+{e}", base + (e,)) for e in extras]
    cd4 = tuple(c for c in ("cd4_rest", "cd4_stim8hr", "cd4_stim48hr") if c in extras)
    if len(cd4) == 3:
        sets.append(("+cd4x3", base + cd4))
        if "kolf2" in extras:
            sets.append(("+kolf2+cd4x3", base + ("kolf2",) + cd4))
    print(f"\n  {'set':<16}{'un_frac':>9}{'un_rho':>8}{'un_M':>7}"
          f"{'pw_frac':>9}{'pw_rho':>8}{'pw_M':>7}")
    res = {}
    for lbl, dons in sets:
        u = score(dons, unanimity_key)
        p_ = score(dons, pairwise_key)
        res[lbl] = dict(donors=list(dons), unanimity=u, pairwise=p_)
        print(f"  {lbl:<16}{u[0]:>9.3f}{u[1]:>8.3f}{u[2]:>7.3f}"
              f"{p_[0]:>9.3f}{p_[1]:>8.3f}{p_[2]:>7.3f}")

    print(f"\n# PERMUTATION NULL on the added donor (unanimity key, {a.n_null} draws)")
    print(f"  {'added donor':<16}{'obs M':>8}{'null p50':>10}{'null p95':>10}{'p':>7}"
          f"{'obs rho':>9}{'rho p50':>9}{'rho p':>7}")
    rng = np.random.default_rng(a.seed)
    nulls_out = {}
    for e in extras:
        if cov[e] < 0.5 * C.P:
            continue
        dons = base + (e,)
        obs_f, obs_rho, obs_M = score(dons, unanimity_key)
        orig = m.pert_index[e]
        nrow = m.s[e].shape[0]
        keys = [k for k in orig if k and orig[k] < nrow]
        vals = np.array([orig[k] for k in keys])
        nm, nr = [], []
        for _ in range(a.n_null):
            m.pert_index[e] = {k: int(v) for k, v in zip(keys, rng.permutation(vals))}
            Mp, Hp = donor_matrices(m, C.perts)
            agr = unanimity_key(Mp, Hp, dons, P, G)
            dc, _ = donor_cos_key(Mp, Hp, dons, P)
            z = np.maximum((rankdata(dc, method="average") - 0.5) / P, 0.3)
            nm.append(C.kappa(z[:, None] * alpha_of(agr) * dh) / k0)
            nr.append(float(spearmanr(dc, truecos).statistic))
        m.pert_index[e] = orig
        nm, nr = np.array(nm), np.array(nr)
        nulls_out[e] = dict(obs_M=obs_M, obs_rho=obs_rho,
                            M_p50=float(np.median(nm)), M_p95=float(np.percentile(nm, 95)),
                            M_p=float((nm >= obs_M).mean()), rho_p50=float(np.median(nr)),
                            rho_p=float((nr >= obs_rho).mean()))
        print(f"  {e:<16}{obs_M:>8.3f}{np.median(nm):>10.3f}{np.percentile(nm,95):>10.3f}"
              f"{(nm>=obs_M).mean():>7.3f}{obs_rho:>9.3f}{np.median(nr):>9.3f}"
              f"{(nr>=obs_rho).mean():>7.3f}")
    p = ROOT / f"experiments/099-donors/work/keys2_{a.harness}.json"
    json.dump(dict(harness=a.harness, base=list(base), kappa0=k0, sets=res,
                   null=nulls_out), open(p, "w"), indent=1)
    print(f"-> {p}")


if __name__ == "__main__":
    main()
