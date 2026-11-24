"""171b — the `keep_top` ladder: can the dispersion live ENTIRELY outside the DE universe?

The scorer's DE table admits only genes whose CONTROL mean exceeds 5 CPM (BRIEF section 3).
On the challenge that is ~9,900 / 9,600 / 10,100 genes of 18,533.  Production's gate is
PSEUDOBULK CPM > 10, which 164 section 2 proved is a strict subset (8,976 genes; `prod-gated
& not DE` is 0 in all three contexts).  So there are ~8,400-9,500 genes per context that:

  * NO DE member can see -- they are outside the 5-CPM control filter, so they never enter
    `n_pred`, `n_real`, the reach pool, the jaccard union or the nmae gate; but
  * the `mse` jackknife DOES see -- `jackknife_correction` runs over the full 18,533-gene
    axis at `bulk_target_sum`.

`emission_spread_vector` initialises the WHOLE vector to `tail_spread` and only then promotes
the top `keep_top` GATED genes to `keep_spread` (164 section 5 calls this the trap).  So those
sub-gate genes already sit at `tail_spread = 1.0, phi = 4` in the shipped build.

That means `--keep-top` has a distinguished value -- anything at or above the gated count --
at which EVERY DE-visible gene sits at `keep_spread` and ALL of the tail-block dispersion
lives outside the DE universe.  If claim still clears the knee there, the whole reach/mse
conflict 164 documented dissolves: a8's flat-0.25 DE behaviour (the best reach law we have
measured, 0.20291) with a15's dispersion credit.

This ladder measures exactly how claim splits between the two blocks.

Same machinery, same known-answer gates as `claim171.py`.
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


def spread_phi(e, bctrl, min_cpm, tail, ktop, kspread, phi_keep, phi_tail):
    """Exactly `emission_spread_vector` + `emit_counts`'s per-gene phi (164's function)."""
    n = e.size
    if tail is None:
        return np.full(n, float(kspread)), np.full(n, float(phi_keep))
    sv = np.full(n, float(tail))                       # trap: initialised to tail_spread
    pv = np.full(n, float(phi_tail))
    expr = np.expm1(bctrl)
    cpm = 1e6 * expr / max(float(expr.sum()), 1e-12)
    gated = (cpm > min_cpm) & (bctrl > 0)
    score_ = np.where(gated, np.abs(np.nan_to_num(e)), -np.inf)
    k = min(int(ktop), int(gated.sum()))
    if k > 0:
        order = np.lexsort((np.arange(n), -score_))[:k]
        sv[order] = float(kspread)
        pv[order] = float(phi_keep)
    return sv, pv


bs = np.load(f"{ROOT}/experiments/096-mse-floor/work/board_solve.npz")
sD, A, B = float(bs["sumD"]), float(bs["v11_A"]), float(bs["v11_B"])
Kc, X, K = float(bs["v11_Kc"]), float(bs["X"]), float(bs["K"])
D = M.load(f"{ROOT}/experiments/096-mse-floor/work/v11-stats.npz")
cpn = dict(np.load(f"{ROOT}/experiments/096-mse-floor/work/cp-proxy-ntc.npz").items())
NSUB = 16
S15 = 0.982 / 1.0796
KNEE = S15 * S15 * B

POLICIES = [
    ("KNOWN flat 0.25/phi4 (a8)", None, 0, 0.25),
    ("a15  keep_top 8000", 1.0, 8000, 0.25),
    ("keep_top 8150", 1.0, 8150, 0.25),
    ("keep_top 8300", 1.0, 8300, 0.25),
    ("keep_top 8400", 1.0, 8400, 0.25),
]


def raw_mse(s, claim):
    return 1 + (s * s * A - 2 * s * X + Kc - min(s * s * B, claim)) / sD


def score(raw):
    return min(max((0.987862879 - raw) / 0.950162566, 0.0), 1.0)


def main():
    print(f"knee s^2 B = {KNEE:.5f}   reference cap K = {K:.5f}   a15 raw mse "
          f"{raw_mse(S15, 0.24581):.5f} score {score(raw_mse(S15, 0.24581)):.5f}")
    print(f"\n{'policy':<28}{'gated':>7}{'C_hat':>9}{'/C_p':>7}{'claim':>9}"
          f"{'branch':>11}{'raw mse':>10}{'score_mse':>11}{'d OVERALL':>11}")
    base = None
    for name, tail, ktop, kspread in POLICIES:
        rows, ng = [], []
        for c in M.CTX:
            d = D[c]
            idx = np.linspace(0, len(d["tgt"]) - 1, NSUB).astype(int)
            vals = []
            for j in range(idx.size):
                e = d["dh"][idx[j]]
                prof = np.expm1(d["bctrl"] + e)
                np.clip(prof, 0, None, out=prof)
                prof /= prof.sum()
                sv, pv = spread_phi(e, d["bctrl"], 10.0, tail, ktop, kspread, 4.0, 4.0)
                if j == 0:
                    expr = np.expm1(d["bctrl"])
                    cc = 1e6 * expr / expr.sum()
                    ng.append(int(((cc > 10.0) & (d["bctrl"] > 0)).sum()))
                libs = np.full(M.NCELL, d["S"][idx[j]] / M.NCELL)
                v, col, ip = M.emit_counts(
                    prof, libs, 0.0, np.random.default_rng(3000 + j), mode="quasi",
                    spread=sv, overdispersion=pv)
                Xm = csr_matrix((v.astype(np.float64), col.astype(np.int64),
                                 ip.astype(np.int64)), shape=(M.NCELL, d["G"]))
                vals.append(jackknife_correction(Xm, np.zeros(M.NCELL, np.intp), 1,
                                                 M.TS, chunk=64)[0])
            cp = float(cpn[c].mean())
            rows.append((float(np.mean(vals)), cp, float(np.mean(np.minimum(vals, cp)))))
        ch = float(np.mean([r[0] for r in rows]))
        cp = float(np.mean([r[1] for r in rows]))
        claim = float(np.mean([300 * r[2] / 18532 for r in rows]))
        raw = raw_mse(S15, claim)
        sc = score(raw)
        if base is None and name.startswith("a15"):
            pass
        print(f"{name:<28}{int(np.mean(ng)):>7}{ch:>9.3f}{ch/cp:>7.3f}{claim:>9.5f}"
              f"{'SATURATED' if claim >= KNEE else 'linear':>11}{raw:>10.5f}{sc:>11.5f}"
              f"{(sc - score(raw_mse(S15, 0.24581)))/6:>+11.6f}", flush=True)

    print("\nGATE: 'KNOWN flat 0.25/phi4 (a8)' must read claim 0.16652 (164 section 5) and")
    print("      'a15  keep_top 8000' must read 0.24581.")


if __name__ == "__main__":
    main()
