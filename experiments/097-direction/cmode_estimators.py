"""Route (c), properly posed: can any BUILDABLE vector estimate the reference's mean response?

common_mode.py established that the numerator splits orthogonally,

    num(u, t) = [ t^2*B + Sw*||u||^2 - 2t*X_d - 2<u, Sd> ] + Kc - rho*claim ,

so a constant offset u added to every predicted effect is a SEPARABLE knob: it does not touch
B (a centred sum of squares), does not touch `claim`, and therefore does not move the rho cliff
or the optimal amplitude.  Its optimum is u = dbar (the reference's own w-weighted mean effect)
and the numerator gain there is Sw*||dbar||^2.  Using a candidate direction v with a free scale
captures exactly cos^2(v, dbar) of that.

So the whole route reduces to one number per candidate: cos(v, dbar).  Candidates measured here
are all buildable at submission time -- none touches the reference.
"""
from __future__ import annotations
import sys
import numpy as np
ROOT = "~/Code/virtual-cell-challenge-2026"
sys.path.insert(0, f"{ROOT}/experiments/096-mse-floor")
sys.path.insert(0, f"{ROOT}/src")
from decompose import terms
from msekit import group_stats
from vcc2026.donor_transfer import Donor

HARNESS = [
    ("h1-test", f"{ROOT}/data/harness/h1-test.h5ad",
     f"{ROOT}/experiments/095-clean-harness/out-h1test/runs/g0.65_q1.0_a1.08_spread0.25/pred.h5ad",
     ["k562_gwps", "orion_hek293t", "orion_hct116", "rpe1", "nadig_jurkat", "k562_essential"]),
    ("orion-thinned", f"{ROOT}/data/harness/orion-hct116-thinned.h5ad",
     f"{ROOT}/experiments/095-clean-harness/out/runs/g0.65_q1.0_a1.08_spread0.25/pred.h5ad",
     ["k562_gwps", "orion_hek293t", "rpe1", "nadig_jurkat", "k562_essential", "h1_2025"]),
]
# production's own two/three donors, so the "as shipped" row is unambiguous
PROD = ("k562_gwps", "orion_hek293t")


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


def donor_mean(name, genes, restrict=None):
    """`m` = that donor's mean delta over its perturbations, on the harness gene axis.

    `restrict` limits the average to a named subset of perturbations (e.g. only the harness's
    own targets) -- which is what a submission could actually do, since it knows its 300.
    """
    d = Donor.load(f"{ROOT}/data/donors/{name}")
    gi = {g: i for i, g in enumerate(d.genes)}
    cols = np.array([gi.get(str(g), -1) for g in genes])
    ok = cols >= 0
    if restrict is None:
        m = np.asarray(d.delta.mean(axis=0), dtype=np.float64)
    else:
        pi = {p: i for i, p in enumerate(d.perts)}
        rows = [pi[p] for p in restrict if p in pi]
        if len(rows) < 5:
            return None, 0, 0
        m = np.asarray(d.delta[np.array(rows)].mean(axis=0), dtype=np.float64)
    out = np.zeros(genes.size)
    out[ok] = m[cols[ok]]
    ncov = int(ok.sum())
    npert = d.delta.shape[0] if restrict is None else len(rows)
    del d
    return out, ncov, npert


def main():
    for name, real, pred, donors in HARNESS:
        t = load(real, pred)
        dh, d, w, genes = t["dh"], t["d"], t["w"], t["genes"]
        Sw = float(w.sum())
        dbar = (w[:, None] * d).sum(0) / Sw
        dhbar = (w[:, None] * dh).sum(0) / Sw
        Sd = (w[:, None] * d).sum(0)                    # = Sw * dbar
        sD, Rc = t["sumD"], Sw * float(dbar @ dbar)
        print(f"\n########## {name}   n={t['n']}  sumD={sD:.6g} ##########")
        print(f"  oracle common-mode numerator gain Sw*||dbar||^2 = {Rc:.6g} "
              f"= {Rc/sD:.4f} raw mse   (R_common/R = {100*Rc/t['R']:.2f} %)")
        print(f"  dbar direction balance: {100*np.mean(dbar[dbar!=0] < 0):.1f} % of nonzero "
              f"genes move DOWN;  ||dbar|| = {np.linalg.norm(dbar):.4f}")
        cands = [("our own dhbar (Miller-to-own-mean)", dhbar, "-", "-")]
        for dn in donors:
            v, nc, npt = donor_mean(dn, genes)
            if v is None:
                continue
            cands.append((f"donor m_full: {dn}", v, nc, npt))
        # the production combination, exactly as MultiDonorModel builds m_full
        vs = [donor_mean(dn, genes)[0] for dn in PROD]
        cands.append((f"m_full({'+'.join(PROD)}) / n  [BETA lever]",
                      sum(vs) / len(vs), "-", "-"))
        # restricted to the harness's own targets -- buildable, we know our 300
        vr = [donor_mean(dn, genes, restrict=set(t["targets"]))[0] for dn in PROD]
        vr = [x for x in vr if x is not None]
        if vr:
            cands.append((f"m over OUR targets only ({'+'.join(PROD)})",
                          sum(vr) / len(vr), "-", "-"))
        print(f"\n  {'candidate v (all buildable)':<44}{'cos(v,dbar)':>13}"
              f"{'capture':>10}{'gain raw mse':>14}{'genes':>8}{'perts':>8}")
        for lbl, v, nc, npt in cands:
            nv = float(np.linalg.norm(v))
            if nv == 0:
                continue
            cs = float(v @ dbar) / (nv * float(np.linalg.norm(dbar)))
            print(f"  {lbl:<44}{cs:>+13.4f}{cs*cs:>10.4f}{cs*cs*Rc/sD:>14.4f}"
                  f"{str(nc):>8}{str(npt):>8}")
        np.savez(f"{ROOT}/experiments/097-direction/work/{name}-cmode.npz",
                 dbar=dbar, dhbar=dhbar, Rc=Rc, Sw=Sw, sumD=sD, genes=genes)


if __name__ == "__main__":
    main()
