"""171x — `mse` price of the expression-compensating spread schedule, challenge panel.

Same machinery and the same known-answer gates as `claim171.py`.  The prediction to falsify:
the schedule only ever RAISES `s` above 0.25, and `msebranch171.py` measured a15 sitting
+0.01313 ABOVE the knee `s^2 B` = 0.23268, so `claim` should land at or near the reference cap
K = 0.24653 and raw `mse` should be EXACTLY a15's.  If it does not, the arm is dead and this
is where it dies.
"""
from __future__ import annotations

import sys

import numpy as np
from scipy.sparse import csr_matrix

ROOT = "~/Code/virtual-cell-challenge-2026"
sys.path.insert(0, f"{ROOT}/src")
sys.path.insert(0, f"{ROOT}/experiments/096-mse-floor")
from cell_eval2.moments import jackknife_correction  # noqa: E402
import mech_alpha_chal as M  # noqa: E402

bs = np.load(f"{ROOT}/experiments/096-mse-floor/work/board_solve.npz")
sD, A, B = float(bs["sumD"]), float(bs["v11_A"]), float(bs["v11_B"])
Kc, X, K = float(bs["v11_Kc"]), float(bs["X"]), float(bs["K"])
D = M.load(f"{ROOT}/experiments/096-mse-floor/work/v11-stats.npz")
cpn = dict(np.load(f"{ROOT}/experiments/096-mse-floor/work/cp-proxy-ntc.npz").items())
NSUB = 16
S15 = 0.982 / 1.0796
KNEE = S15 * S15 * B
PHI, S_LO, S_HI, A_CTRL = 4.0, 0.25, 1.0, 2.0
MIN_CPM = 10.0


def a15_spread(e, bctrl, ktop=8000):
    n = e.size
    sv = np.full(n, 1.0)
    expr = np.expm1(bctrl)
    cpm = 1e6 * expr / max(float(expr.sum()), 1e-12)
    gated = (cpm > MIN_CPM) & (bctrl > 0)
    sc = np.where(gated, np.abs(np.nan_to_num(e)), -np.inf)
    k = min(int(ktop), int(gated.sum()))
    if k:
        sv[np.lexsort((np.arange(n), -sc))[:k]] = 0.25
    return sv


def xcomp_spread(m, gated, q, rng=None, permute=False):
    """Identical to `realizex171.xcomp_spread` -- one definition, two surfaces."""
    sv = np.full(m.size, S_HI)
    g = np.flatnonzero(gated)
    mg = np.maximum(m[g], 1e-12)
    c = float(np.quantile(np.sqrt(mg) / np.sqrt(A_CTRL + S_LO), q))
    c2 = c * c
    disc = c2 * c2 - 4.0 * PHI * c2 * mg * (c2 * A_CTRL - mg)
    s = np.full(mg.size, S_LO)
    ok = disc >= 0
    s[ok] = (-c2 + np.sqrt(disc[ok])) / (2.0 * PHI * c2 * mg[ok])
    s = np.clip(s, S_LO, S_HI)
    if permute:
        s = s[rng.permutation(s.size)]
    sv[g] = s
    return sv


def rank_spread(m, gated):
    sv = np.full(m.size, S_HI)
    g = np.flatnonzero(gated)
    o = np.lexsort((g, m[g]))
    r = np.empty(g.size)
    r[o] = np.arange(g.size) / max(g.size - 1, 1)
    sv[g] = S_LO + (S_HI - S_LO) * r
    return sv


POLICIES = ["KNOWN flat 0.25/phi4 (a8)", "a15 AS SHIPPED", "xcomp10", "xcomp25",
            "xcomp10_perm", "sgrad_expr"]


def raw_mse(s, claim):
    return 1 + (s * s * A - 2 * s * X + Kc - min(s * s * B, claim)) / sD


def score(raw):
    return min(max((0.987862879 - raw) / 0.950162566, 0.0), 1.0)


def main():
    print(f"knee s^2 B = {KNEE:.5f}   cap K = {K:.5f}   a15 raw mse "
          f"{raw_mse(S15, 0.24581):.5f}  score {score(raw_mse(S15, 0.24581)):.5f}")
    print(f"\n{'policy':<28}{'C_hat':>9}{'/C_p':>7}{'claim':>9}{'branch':>11}"
          f"{'raw mse':>10}{'score_mse':>11}{'d OVERALL':>11}{'med s':>8}")
    out = {}
    for name in POLICIES:
        rows, meds = [], []
        for c in M.CTX:
            d = D[c]
            idx = np.linspace(0, len(d["tgt"]) - 1, NSUB).astype(int)
            vals = []
            for j in range(idx.size):
                e = d["dh"][idx[j]]
                prof = np.expm1(d["bctrl"] + e)
                np.clip(prof, 0, None, out=prof)
                prof /= prof.sum()
                L = d["S"][idx[j]] / M.NCELL
                m_cell = np.floor(M.NCELL * L * prof) / M.NCELL
                expr = np.expm1(d["bctrl"])
                cc = 1e6 * expr / expr.sum()
                gated = (cc > MIN_CPM) & (d["bctrl"] > 0)
                if name.startswith("KNOWN"):
                    sv = np.full(e.size, 0.25)
                elif name.startswith("a15"):
                    sv = a15_spread(e, d["bctrl"])
                elif name == "xcomp10":
                    sv = xcomp_spread(m_cell, gated, 0.10)
                elif name == "xcomp25":
                    sv = xcomp_spread(m_cell, gated, 0.25)
                elif name == "xcomp10_perm":
                    sv = xcomp_spread(m_cell, gated, 0.10,
                                      np.random.default_rng(11 + j), True)
                else:
                    sv = rank_spread(m_cell, gated)
                meds.append(float(np.median(sv[gated])))
                v, col, ip = M.emit_counts(
                    prof, np.full(M.NCELL, L), 0.0, np.random.default_rng(3000 + j),
                    mode="quasi", spread=sv, overdispersion=PHI)
                Xm = csr_matrix((v.astype(np.float64), col.astype(np.int64),
                                 ip.astype(np.int64)), shape=(M.NCELL, d["G"]))
                vals.append(jackknife_correction(Xm, np.zeros(M.NCELL, np.intp), 1,
                                                 M.TS, chunk=64)[0])
            cp = float(cpn[c].mean())
            rows.append((float(np.mean(vals)), cp, float(np.mean(np.minimum(vals, cp)))))
        ch = float(np.mean([r[0] for r in rows]))
        cp = float(np.mean([r[1] for r in rows]))
        claim = float(np.mean([300 * r[2] / 18532 for r in rows]))
        raw, sc = raw_mse(S15, claim), score(raw_mse(S15, claim))
        out[name] = claim
        print(f"{name:<28}{ch:>9.3f}{ch/cp:>7.3f}{claim:>9.5f}"
              f"{'SATURATED' if claim >= KNEE else 'linear':>11}{raw:>10.5f}{sc:>11.5f}"
              f"{(sc - score(raw_mse(S15, 0.24581)))/6:>+11.6f}{np.mean(meds):>8.3f}",
              flush=True)

    print("\n--- known-answer gates ---")
    bad = []
    for nm, want, tol in (("KNOWN flat 0.25/phi4 (a8)", 0.16652, 0.03),
                          ("a15 AS SHIPPED", 0.24581, 0.02)):
        got = out[nm]
        ok = abs(got - want) / want < tol
        print(f"  {nm:<28} claim {got:.5f}  expect {want:.5f}  "
              f"({100*(got-want)/want:+.2f}%)  {'PASS' if ok else 'FAIL'}")
        if not ok:
            bad.append(nm)
    if bad:
        sys.exit(f"ABORT: {bad}")


if __name__ == "__main__":
    main()
