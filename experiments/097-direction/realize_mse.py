"""REALISE the recommended arms: re-render, re-emit through the production emitter, re-jackknife,
re-score.  Nothing frozen, so `claim`, the #247 cap and the rho bound all move as the scorer sees.

This is the gate CURRENT_STATE item 6 demands: "a reordering that cannot be produced by an actual
spread vector is not a result".  Everything in exact.py is an analytic evaluation of the identity
on a FIXED dh; this drives `emit_counts` and `jackknife_correction` for real.

Emission policy is held IDENTICAL across arms, and set to a cap-reaching one so the operating
point matches a12 (claim/B near the challenge's 0.88 rather than the shipped 0.18).

    uv run python realize_mse.py <h1|orion> [spread] [phi]
"""
from __future__ import annotations
import os, sys, time
import numpy as np
from scipy.sparse import csr_matrix
ROOT = "~/Code/virtual-cell-challenge-2026"
sys.path.insert(0, f"{ROOT}/src")
sys.path.insert(0, f"{ROOT}/experiments/096-mse-floor")
from cell_eval2.moments import jackknife_correction            # noqa: E402
from cell_eval2.prep import bulk_lognorm_means                 # noqa: E402
from cell_eval2.metrics.delta import _across_pert_budget       # noqa: E402
from decompose import terms                                    # noqa: E402
from msekit import CTRL, TS, group_stats, score                # noqa: E402
from vcc2026.donor_transfer import emit_counts, Donor          # noqa: E402
from scipy.stats import rankdata                                # noqa: E402

CFG = {
    "h1": (f"{ROOT}/data/harness/h1-test.h5ad",
           f"{ROOT}/experiments/095-clean-harness/out-h1test/runs/g0.65_q1.0_a1.08_spread0.25/pred.h5ad"),
    "orion": (f"{ROOT}/data/harness/orion-hct116-thinned.h5ad",
              f"{ROOT}/experiments/095-clean-harness/out/runs/g0.65_q1.0_a1.08_spread0.25/pred.h5ad"),
}
DONORS = ("k562_gwps", "orion_hek293t")


def rank01(s):
    """MIDRANKS. Ordinal ranks break ties by gene index, and ~79 % of every agreement row is
    tied at exactly zero (k562_gwps covers 7,679 of 18,533 panel genes, and sign disagreements
    zero more), so an ordinal rank hands 79 % of the genes an arbitrary index-ordered ramp.
    Found in the production `agreement_shrink_vector` by the `miller` session and present in my
    own first pass; `tiefix.py` measures what it was worth."""
    return (rankdata(s, method="average") - 0.5) / len(s)


def keys(genes, targets):
    Ms, hv = {}, {}
    P, G = len(targets), genes.size
    for nmd in DONORS:
        dn = Donor.load(f"{ROOT}/data/donors/{nmd}")
        m, s = dn.decompose()
        gi = {g: i for i, g in enumerate(dn.genes)}
        cols = np.array([gi.get(str(g), -1) for g in genes]); ok = cols >= 0
        pi = {p: i for i, p in enumerate(dn.perts)}
        M = np.zeros((P, G)); h = np.zeros(P, bool)
        for r, x in enumerate(targets):
            if x in pi:
                M[r, ok] = s[pi[x]][cols[ok]]; h[r] = True
        Ms[nmd] = M / (float(np.median(np.linalg.norm(s[:, cols[ok]], axis=1))) or 1.0)
        hv[nmd] = h
        del dn, s
    A_, B_ = Ms[DONORS[0]], Ms[DONORS[1]]
    both = hv[DONORS[0]] & hv[DONORS[1]]
    na, nb = np.linalg.norm(A_, axis=1), np.linalg.norm(B_, axis=1)
    dcos = np.where(both, np.einsum("ij,ij->i", A_, B_) / np.maximum(na * nb, 1e-12), 0.0)
    ag = np.sqrt(np.abs(A_) * np.abs(B_)) * ((np.sign(A_) == np.sign(B_)) & (A_ != 0) & (B_ != 0))
    al = np.empty_like(ag)
    for i in range(P):
        al[i] = rankdata(ag[i], method="average") / G
        if np.unique(ag[i]).size < 2:      # <2 donors: no cross-donor evidence, do not damp
            al[i] = 1.0
    # SINGLE-DONOR POLICY, settled by measurement in singledonor.py: a target carried by fewer
    # than two donors keeps donor_cos = 0, so it ranks LOW and is damped -- on h1-test those 4
    # rows carry 13.1 % of the predicted amplitude A and +0.5 % of X (mean cos_p -0.0012 against
    # +0.0369 for the rest). Holding them at z = 1 instead turns a 1.24x gain into 0.96x.
    # The FLOOR keeps them emitting rather than pinning them at zero.
    return rank01(dcos), al


def main(which, spread, phi):
    REAL, PRED = CFG[which]
    rp, rs, rc, rj, rg = group_stats(REAL)
    pp, ps, pc, pj, pg = group_stats(PRED)
    assert np.array_equal(rg, pg)
    G = rg.size
    ci = int(np.flatnonzero(rp == CTRL)[0])
    bctrl = bulk_lognorm_means(rs[ci][None, :].astype(np.float64), TS)[0]
    dh0 = bulk_lognorm_means(ps.astype(np.float64), TS) - bctrl[None, :]
    S = ps.sum(1).astype(np.float64)
    NCELL = int(round(float(np.median(pc))))
    tb = terms(pp, ps, pc, pj, rp, rs, rc, rj, rg)
    Kcap = tb["K"]
    print(f"# {which}: P={pp.size} NCELL={NCELL} G={G}  sumD={tb['sumD']:.6f}  "
          f"K (#247 cap on claim)={Kcap:.6f}  shipped claim={tb['claim']:.6f}", flush=True)

    # keys are defined on the SCORED rows only; the control row (and any row the reference
    # does not carry) keeps weight 1 so the emitted control block is untouched.
    keep = tb["keep"]
    zk, alk = keys(rg, np.array([str(x) for x in pp[keep]]))
    def full(vec_or_mat, per_gene):
        M = np.ones_like(dh0)
        M[keep] = vec_or_mat if per_gene else vec_or_mat[:, None]
        return M
    FL = 0.3
    zf = np.maximum(zk, FL) ** 2
    ARMS = {
        "base": np.ones_like(dh0),
        "pert^2 fl.3": full(zf, False),
        "gene^1": full(alk, True),
        "STACK p2f3g1": full(zf[:, None] * alk, True),
    }
    print(f"# mechanism check: {keep.size} scored rows; z=max(rank01,0.3)^2 in "
          f"[{np.maximum(zk,0.3).min()**2:.4f}, {np.maximum(zk,0.3).max()**2:.4f}], "
          f"mean {np.mean(np.maximum(zk,0.3)**2):.4f} -> alpha compensation "
          f"{1/np.sqrt(np.mean(np.maximum(zk,0.3)**2)):.3f}x; per-gene alpha mean "
          f"{alk.mean():.3f}; control row weight 1.0", flush=True)

    def run(dh, s, tag):
        t0 = time.time()
        sums = np.zeros_like(ps); jks = np.zeros_like(pj)
        for i in range(pp.size):
            prof = np.expm1(bctrl + s * dh[i]); np.clip(prof, 0.0, None, out=prof)
            tot = prof.sum()
            if tot <= 0:
                continue
            prof = prof / tot
            libs = np.full(NCELL, S[i] / NCELL)
            v, col, ip = emit_counts(prof, libs, 0.0, np.random.default_rng(7919 + i),
                                     mode="quasi", spread=spread, overdispersion=phi)
            Xm = csr_matrix((v.astype(np.float64), col.astype(np.int64), ip.astype(np.int64)),
                            shape=(NCELL, G))
            sums[i] = np.asarray(Xm.sum(0)).ravel().astype(np.int64)
            jks[i] = jackknife_correction(Xm, np.zeros(NCELL, np.intp), 1, TS, chunk=64)[0]
        cnt = np.full(pp.size, float(NCELL))
        m, _, _ = score(pp, sums, cnt, jks, rp, rs, rc, rj, rg)
        t = terms(pp, sums, cnt, jks, rp, rs, rc, rj, rg)
        print(f"  {tag:<14} s={s:<6.3f} mse={m:.5f}  A={t['A']:.6f} B={t['budget']:.6f} "
              f"claim={t['claim']:.6f} (cap {Kcap:.4f}) rho={t['rho']:.3f} X={t['X']:.6f} "
              f"eps={(t['A']-t['budget'])/t['budget']:+.5f}  [{time.time()-t0:.0f}s]", flush=True)
        return m, t

    best = {}
    for tag, Wm in ARMS.items():
        dh = Wm * dh0
        B1 = float(_across_pert_budget(dh[tb["keep"]], tb["cols"], tb["w"]))
        # one probe at s=1 fixes claim (flat in s, 096 measured 0.01 %), then the cliff
        _, t1 = run(dh, 1.0, tag + " probe")
        sc = float(np.sqrt(t1["claim"] / B1))
        grid = sorted({round(sc * f, 4) for f in (0.92, 1.06)})
        rows = [run(dh, s, tag) for s in grid]
        j = int(np.argmin([r[0] for r in rows]))
        best[tag] = (grid[j], rows[j][0], rows[j][1])
        print(f"  ==> {tag:<14} cliff s_c={sc:.4f}  realised min mse={rows[j][0]:.5f} "
              f"at s={grid[j]:.4f}", flush=True)
    print(f"\n# {which} SUMMARY (realised, same emission policy spread={spread} phi={phi})")
    b = best["base"][1]
    for tag, (s, m, t) in best.items():
        kap = t["X"] / np.sqrt(t["budget"] * t["sumD"])
        print(f"  {tag:<14} mse={m:.5f}  d vs base {m-b:+.5f}   kappa={kap:.5f}  "
              f"eps={(t['A']-t['budget'])/t['budget']:+.5f}  claim/B={t['claim']/t['budget']:.3f}")


if __name__ == "__main__":
    main(sys.argv[1], float(sys.argv[2]) if len(sys.argv) > 2 else 0.5,
         float(sys.argv[3]) if len(sys.argv) > 3 else 4.0)
