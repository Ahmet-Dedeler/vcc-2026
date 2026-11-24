"""185 Part E — `mse` through 096's `claim` route, on the REAL challenge panel.

Same machinery and the same known-answer gates as `claimx171.py` / `claim177.py`:
  flat 0.25/phi4 (a8) -> claim 0.16652 (+-3%),   a15 as shipped -> 0.24581 (+-2%).

The scored numerator is 1 + (s^2 A - 2 s X + Kc - min(s^2 B, claim)) / sumD and a15 sits at
claim 0.24581 against the knee s^2 B = 0.23268 (171 §5): `mse` is EXACTLY flat for any arm
whose claim stays >= the knee and falls linearly below it.  `claim` reads the delete-one
jackknife C_hat of OUR emitted cells, so a new emission law must be priced here even though
its pseudobulk is pinned.

Arms priced: a15 (KA), flat (KA), rc_hist_a, rc_hist_s, rc_perm, rc_k2000, rc_k4000.  The
rank-controlled arms are built exactly as on h1-test: a15's own |z| histogram on the SAME
target against the context's REAL controls, re-keyed by |effect| (the delivered log-profile
effect `dh`, the challenge-side analogue of |H|), realised by the exact-marginal emitter with
the sub-gate genes on a15's tail law.
"""
from __future__ import annotations

import os
import sys
import time

import anndata as ad
import numpy as np
from scipy.sparse import csr_matrix

ROOT = "~/Code/virtual-cell-challenge-2026"
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, f"{ROOT}/src")
sys.path.insert(0, f"{ROOT}/experiments/096-mse-floor")
sys.path.insert(0, f"{ROOT}/experiments/092-call-set-truncation")
from cell_eval2.moments import jackknife_correction  # noqa: E402
import mech_alpha_chal as M  # noqa: E402
from de_kernel import ctrl_tiesums, wilcoxon_z  # noqa: E402
from eme185 import EME  # noqa: E402

bs = np.load(f"{ROOT}/experiments/096-mse-floor/work/board_solve.npz")
sD, A, B = float(bs["sumD"]), float(bs["v11_A"]), float(bs["v11_B"])
Kc, X, K = float(bs["v11_Kc"]), float(bs["X"]), float(bs["K"])
D = M.load(f"{ROOT}/experiments/096-mse-floor/work/v11-stats.npz")
cpn = dict(np.load(f"{ROOT}/experiments/096-mse-floor/work/cp-proxy-ntc.npz").items())
NSUB = 16
S15 = 0.982 / 1.0796
KNEE = S15 * S15 * B
PHI, MIN_CPM_PB, MIN_CPM_DE, n = 4.0, 10.0, 5.0, 400
POLICIES = ["KNOWN flat 0.25/phi4 (a8)", "a15 AS SHIPPED", "rc_hist_a", "rc_hist_s", "rc_perm",
            "rc_k2000", "rc_k4000"]


def a15_spread(e, bctrl, ktop=8000):
    """claimx171's exact construction (production's `emission_spread_vector` shape)."""
    m = e.size
    sv = np.full(m, 1.0)
    expr = np.expm1(bctrl)
    cpm = 1e6 * expr / max(float(expr.sum()), 1e-12)
    gated = (cpm > MIN_CPM_PB) & (bctrl > 0)
    sc = np.where(gated, np.abs(np.nan_to_num(e)), -np.inf)
    k = min(int(ktop), int(gated.sum()))
    if k:
        sv[np.lexsort((np.arange(m), -sc))[:k]] = 0.25
    return sv


def raw_mse(s, claim):
    return 1 + (s * s * A - 2 * s * X + Kc - min(s * s * B, claim)) / sD


def score(raw):
    return min(max((0.987862879 - raw) / 0.950162566, 0.0), 1.0)


def to_block(v, col, ip, G):
    Bm = np.zeros((n, G), np.int64)
    for r in range(n):
        sl = slice(ip[r], ip[r + 1])
        Bm[r, col[sl]] = v[sl]
    return Bm


def load_ctx_controls(c, G):
    """Sorted log1p-CPM control block over the scorer's DE universe (mean CPM > 5)."""
    t0 = time.time()
    Ad = ad.read_h5ad(f"{ROOT}/data/challenge/context_{c}.h5ad")
    Xc = Ad.X.tocsr()
    assert Xc.shape[1] == G
    libs = np.asarray(Xc.sum(1)).ravel()
    libs[libs == 0] = 1.0
    cpm = np.asarray(Xc.multiply(1e6 / libs[:, None]).mean(0)).ravel()
    keep = cpm > MIN_CPM_DE
    Gk = int(keep.sum())
    m = Xc.shape[0]
    cs = np.empty((Gk, m))
    for i in range(0, m, 2000):
        blk = np.asarray(Xc[i:i + 2000].todense(), dtype=np.float64)
        s = blk.sum(1)
        s[s == 0] = 1.0
        cs[:, i:i + blk.shape[0]] = np.log1p(blk[:, keep] * (1e6 / s[:, None])).T
    cs.sort(axis=1)
    cs = np.ascontiguousarray(cs)
    tie = ctrl_tiesums(cs)
    print(f"  ctx {c}: {m} controls, {Gk}/{G} genes in the DE universe [{time.time()-t0:.0f}s]",
          flush=True)
    return cs, tie, keep


def main():
    print(f"knee s^2 B = {KNEE:.5f}   cap K = {K:.5f}   a15 raw mse "
          f"{raw_mse(S15, 0.24581):.5f}  score {score(raw_mse(S15, 0.24581)):.5f}")
    out = {p: [] for p in POLICIES}          # per context (C_hat mean, C_p, claim part)
    npred = {p: [] for p in POLICIES}
    for c in M.CTX:
        d = D[c]
        G = d["G"]
        cs, tie, keep = load_ctx_controls(c, G)
        eme = EME(cs, tie, 20000, n)
        gidx = np.flatnonzero(keep)
        ref_lin = None
        idx = np.linspace(0, len(d["tgt"]) - 1, NSUB).astype(int)
        vals = {p: [] for p in POLICIES}
        for j in range(idx.size):
            e = d["dh"][idx[j]]
            prof = np.expm1(d["bctrl"] + e)
            np.clip(prof, 0, None, out=prof)
            prof /= prof.sum()
            L = d["S"][idx[j]] / M.NCELL
            Li = int(L)
            totals = np.floor(n * Li * prof).astype(np.int64)   # what the quasi law delivers
            S = int(totals.sum())
            T = S // n
            Tc = np.full(n, T, np.int64)
            Tc[: S % n] += 1
            sv15 = a15_spread(e, d["bctrl"])
            blocks = {}
            for name in POLICIES:
                rng = np.random.default_rng(3000 + j)
                if name.startswith("KNOWN") or name.startswith("a15"):
                    sv = np.full(e.size, 0.25) if name.startswith("KNOWN") else sv15
                    v, col, ip = M.emit_counts(prof, np.full(n, L), 0.0, rng, mode="quasi",
                                               spread=sv, overdispersion=PHI)
                    Bm = to_block(v, col, ip, G)
                    if name.startswith("a15"):
                        # template |z| against the REAL context controls
                        s_ = Bm.sum(1).astype(np.float64)
                        s_[s_ == 0] = 1.0
                        nb = np.log1p(Bm[:, keep] * (1e6 / s_[:, None]))
                        z15 = wilcoxon_z(cs, tie, np.ascontiguousarray(nb.T))
                        blocks["z15"] = z15
                else:
                    z15 = blocks["z15"]
                    Gk = gidx.size
                    zabs_sorted = np.sort(np.abs(z15))[::-1]
                    key = np.abs(np.nan_to_num(e))[keep]
                    jit = np.random.default_rng(1000 + j).random(Gk) * 1e-12
                    if name == "rc_perm":
                        ro = np.random.default_rng(5000 + j).permutation(Gk)
                    else:
                        ro = np.argsort(-(key + jit), kind="stable")
                    pos = np.empty(Gk, np.int64)
                    pos[ro] = np.arange(Gk)
                    zt = zabs_sorted[pos]
                    if name.startswith("rc_k"):
                        zt = np.where(pos < int(name[4:]), zt, 0.0)
                    sgn = np.sign(np.nan_to_num(e))[keep]
                    sgn[sgn == 0] = 1.0
                    dsg = eme.design(totals[keep], sgn * zt,
                                     mode="signed" if name == "rc_hist_s" else "abs", T=T)
                    prof_sub = np.where(keep, 0.0, (totals + 0.5) / (n * Li))
                    v2, c2, p2 = M.emit_counts(prof_sub, np.full(n, float(Li)), 0.0, rng,
                                               mode="quasi", spread=1.0, overdispersion=PHI)
                    B_sub = to_block(v2, c2, p2, G)
                    assert B_sub[:, keep].sum() == 0
                    Bm, dg = eme.place(dsg, gidx, B_sub, Tc, rng)
                    assert dg["unresolved"] == 0
                assert np.array_equal(Bm.sum(0), totals), name
                Xm = csr_matrix(Bm.astype(np.float64))
                vals[name].append(jackknife_correction(Xm, np.zeros(n, np.intp), 1, M.TS,
                                                       chunk=64)[0])
                if name != "KNOWN flat 0.25/phi4 (a8)":
                    s_ = Bm.sum(1).astype(np.float64)
                    nb = np.log1p(Bm[:, keep] * (1e6 / s_[:, None]))
                    zz = wilcoxon_z(cs, tie, np.ascontiguousarray(nb.T))
                    from scipy.stats import false_discovery_control, norm
                    q = false_discovery_control(2 * norm.sf(np.abs(zz)), method="bh")
                    npred[name].append(int((q < 0.05).sum()))
            print(f"    ctx {c} target {j+1}/{NSUB} done", flush=True)
        cp = float(cpn[c].mean())
        for name in POLICIES:
            arr = np.asarray(vals[name])
            out[name].append((float(arr.mean()), cp, float(np.minimum(arr, cp).mean())))
        del cs, eme

    print(f"\n{'policy':<28}{'C_hat':>9}{'/C_p':>7}{'claim':>9}{'branch':>11}"
          f"{'raw mse':>10}{'score_mse':>11}{'d OVERALL':>11}{'n_pred':>8}")
    claims = {}
    for name in POLICIES:
        rows = out[name]
        ch = float(np.mean([r[0] for r in rows]))
        cp = float(np.mean([r[1] for r in rows]))
        claim = float(np.mean([300 * r[2] / 18532 for r in rows]))
        claims[name] = claim
        raw, sc = raw_mse(S15, claim), score(raw_mse(S15, claim))
        npd = f"{np.median(npred[name]):.0f}" if npred[name] else "-"
        print(f"{name:<28}{ch:>9.3f}{ch/cp:>7.3f}{claim:>9.5f}"
              f"{'SATURATED' if claim >= KNEE else 'linear':>11}{raw:>10.5f}{sc:>11.5f}"
              f"{(sc - score(raw_mse(S15, 0.24581)))/6:>+11.6f}{npd:>8}", flush=True)
    print("\n--- known-answer gates ---")
    bad = []
    for nm, want, tol in (("KNOWN flat 0.25/phi4 (a8)", 0.16652, 0.03),
                          ("a15 AS SHIPPED", 0.24581, 0.02)):
        got = claims[nm]
        ok = abs(got - want) / want < tol
        print(f"  {nm:<28} claim {got:.5f}  expect {want:.5f}  "
              f"({100*(got-want)/want:+.2f}%)  {'PASS' if ok else 'FAIL'}")
        if not ok:
            bad.append(nm)
    np.savez(os.path.join(HERE, "work", "claim185.npz"), policies=np.array(POLICIES),
             claims=np.array([claims[p] for p in POLICIES]))
    if bad:
        sys.exit(f"ABORT: {bad}")


if __name__ == "__main__":
    main()
