"""Is the `C_p` proxy right? Measure it where real perturbation groups ARE visible.

The `#247` cap is `min(C_hat_p, C_p)` with `C_p` the REFERENCE's own delete-1 jackknife. On the
challenge that is withheld, so 096 substitutes a proxy: the jackknife of 400 CONTROL cells, one
ntc construct at a time. `mech_challenge.py` calls it "a mild UNDER-estimate ... conservative in
the direction that matters", and that sentence is an ARGUMENT, not a measurement -- it was never
tested, because until a12 no arm's claim was anywhere near the cap.

a12 is the first arm whose claim IS the cap, so the proxy is now load-bearing and testable.
`data/harness/h1-test.h5ad` has exactly the challenge's geometry -- 85 perturbations x 400 cells
plus 18,400 control cells -- so the same substitution can be made there and checked against the
real perturbation groups.

Reports, per surface: mean and spread of the real per-perturbation `C_p`, the same for 400-cell
control draws, the bias ratio, the across-perturbation CV, and the per-cell depth of each side
(RESEARCH-METHOD rule 4: orion's controls are ~38% deeper than its perturbation cells, and depth
moves a sampling correction directly).

Then it re-prices a12's `claim` under the measured bias and CV.

    uv run python mech_cp_bias.py
"""
from __future__ import annotations

import sys

import numpy as np

ROOT = "~/Code/virtual-cell-challenge-2026"
sys.path.insert(0, f"{ROOT}/src")
sys.path.insert(0, f"{ROOT}/experiments/096-mse-floor")
from cell_eval2.moments import jackknife_correction  # noqa: E402

TS = 50_000.0
CTX = ("A", "B", "C")
W = f"{ROOT}/experiments/096-mse-floor/work"


def surface(path, n_draw=46, seed=0):
    import anndata as ad
    from scipy.sparse import csr_matrix, issparse
    print(f"\n===== {path.split('/')[-1]} =====", flush=True)
    A = ad.read_h5ad(path)
    lab = np.asarray(A.obs["target_gene"].astype(str).values)
    X = A.X.tocsr() if issparse(A.X) else csr_matrix(A.X)
    perts, codes = np.unique(lab, return_inverse=True)
    P = perts.size
    jk = jackknife_correction(X, codes, P, TS, chunk=64)
    ci = int(np.flatnonzero(perts == "non-targeting")[0])
    keep = np.setdiff1d(np.arange(P), [ci])
    n_pert = int((codes == keep[0]).sum())
    Cp = jk[keep]
    # 400-cell (or n_pert-cell) draws from the control pool -- the proxy, on this surface
    crows = np.flatnonzero(codes == ci)
    rng = np.random.default_rng(seed)
    Cc = []
    for _ in range(n_draw):
        r = np.sort(rng.choice(crows, n_pert, replace=False))
        sub = csr_matrix(X[r])
        Cc.append(jackknife_correction(sub, np.zeros(n_pert, np.intp), 1, TS, chunk=64)[0])
    Cc = np.asarray(Cc)
    dep = np.asarray(X.sum(axis=1)).ravel()
    dp, dc = dep[codes != ci], dep[codes == ci]
    print(f"  {len(Cp)} perturbation groups of {n_pert} cells; {crows.size} control cells")
    print(f"  depth  median: perturbation {np.median(dp):9.1f}   control {np.median(dc):9.1f}"
          f"   ratio ctrl/pert {np.median(dc)/np.median(dp):.4f}")
    print(f"  C_p  (REAL perturbation groups): mean {Cp.mean():8.4f}  sd {Cp.std():7.4f}  "
          f"CV {Cp.std()/Cp.mean():.4f}  p5 {np.percentile(Cp,5):8.4f}  "
          f"p95 {np.percentile(Cp,95):8.4f}")
    print(f"  C    ({n_pert}-cell CONTROL draws, the proxy): mean {Cc.mean():8.4f}  "
          f"sd {Cc.std():7.4f}  CV {Cc.std()/Cc.mean():.4f}")
    print(f"  >>> bias  mean(C_p real) / mean(C proxy) = {Cp.mean()/Cc.mean():.4f}")
    print(f"      the proxy is {'an UNDER' if Cp.mean()>Cc.mean() else 'an OVER'}-estimate "
          f"by {abs(100*(Cp.mean()/Cc.mean()-1)):.2f}%")
    print(f"      across-perturbation CV of the real C_p: {Cp.std()/Cp.mean():.4f}  "
          f"(the proxy replaces this spread with a constant)")
    return dict(Cp=Cp, Cc=Cc, ratio=float(Cp.mean() / Cc.mean()),
                cv=float(Cp.std() / Cp.mean()), n=n_pert)


def reprice(bias, cv, arm="a12", seed=0, n_mc=200):
    """`claim` for a CAPPED arm when `C_p` is (bias x proxy) with across-perturbation CV `cv`.

    Two separable effects, reported apart:
      * LEVEL   -- the proxy's mean is off by `bias`;
      * SPREAD  -- `min(C_hat, C_p)` is concave, so replacing a spread `C_p` by its mean
                   OVER-states the claim by `E[(C_p - C_hat)+]`, always.

    `X` and `sumD` come from the LS solve on a8+a10+a11 -- the three arms whose `claim` is
    their own measured `C_hat_p` and whose cap therefore never binds, so those constants carry
    no `C_p` risk and the sensitivity below is not circular. #348's budget is taken PER
    CONTEXT, as the scorer does.
    """
    import sys as _s
    _s.path.insert(0, f"{ROOT}/experiments/096-mse-floor")
    from mech_a12_audit import ALPHA, ALPHA_PROD, cp_proxies, measure, panel
    X, sD = 0.019240, 0.229141
    board = {"a12": 0.9013492}
    cps = cp_proxies()
    rng = np.random.default_rng(seed)
    print(f"\n===== {arm}'s claim re-priced (per-context branch; X={X:.6f} sumD={sD:.6f}) =====")
    print(f"  {'bias':>7}{'cv':>7}{'claim':>10}{'vs proxy':>10}{'rho_min':>9}{'raw mse':>10}"
          f"{'vs board':>10}")
    out = []
    for b, c in ((1.0, 0.0), (bias, 0.0), (1.0, cv), (bias, cv)):
        rows = measure(f"{W}/{arm}-stats.npz", {x: cps["ntc"][x][0] * b for x in CTX})
        if c > 0:
            new = []
            for r in rows:
                vals = [float(np.minimum(r["Chat"],
                                         r["Cp"] * (1 + c * rng.standard_normal(r["n"]))).sum()
                              / (r["G"] - 1.0)) for _ in range(n_mc)]
                r = dict(r); r["claim"] = float(np.mean(vals)); new.append(r)
            rows = new
        m = panel(rows)
        s = ALPHA[arm] / ALPHA_PROD
        raw = 1.0 + (m["A"] + m["Kc"] - m["rho_claim"] - 2 * s * X) / sD
        out.append((b, c, m["claim"], raw))
        print(f"  {b:>7.4f}{c:>7.4f}{m['claim']:>10.6f}{m['claim']-0.245917:>10.6f}"
              f"{m['rho_min']:>9.4f}{raw:>10.6f}{raw-board[arm]:>10.6f}")
    return out


def main():
    h1 = surface(f"{ROOT}/data/harness/h1-test.h5ad")
    orion = surface(f"{ROOT}/data/harness/orion-hct116-thinned.h5ad")
    print("\n  RULE 1: h1-test has the challenge's exact geometry (400 cells/pert, 18,400")
    print("  controls); orion is 150 cells at 12k UMI and carries the known depth confound, so")
    print("  it is reported as a second reading, not as an equal one.")
    reprice(h1["ratio"], h1["cv"])
    reprice(orion["ratio"], orion["cv"])
    np.savez(f"{W}/cp_bias.npz", h1_Cp=h1["Cp"], h1_Cc=h1["Cc"],
             orion_Cp=orion["Cp"], orion_Cc=orion["Cc"],
             h1_ratio=h1["ratio"], h1_cv=h1["cv"],
             orion_ratio=orion["ratio"], orion_cv=orion["cv"])


if __name__ == "__main__":
    main()
